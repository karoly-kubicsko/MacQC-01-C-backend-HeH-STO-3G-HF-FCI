#include "macqc/gto_derivatives.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <map>
#include <stdexcept>

#include "macqc/job.h"
#include "macqc/parallel.h"
#include "macqc/rhf.h"

namespace macqc {

namespace {

inline int ncart(int l) { return (l + 1) * (l + 2) / 2; }
inline int cart_index(int /*x*/, int y, int z) { return ((y + z) * (y + z + 1)) / 2 + z; }

std::vector<std::array<int, 3>> cartesian_components(int l) {
    std::vector<std::array<int, 3>> c;
    for (int x = l; x >= 0; --x)
        for (int y = l - x; y >= 0; --y) c.push_back({x, y, l - x - y});
    return c;
}

// Component shell with angular momentum l + dl and coefficients d_i (2 a_i)^p.
Shell component_shell(const Shell& s, int dl, int power) {
    Shell c = s;
    c.l = s.l + dl;
    for (std::size_t i = 0; i < c.coefficients.size(); ++i)
        c.coefficients[i] = s.coefficients[i] * std::pow(2.0 * s.exponents[i], power);
    return c;
}

// Accumulates rows: row r of the component transform gets coefficient v at column col.
struct RowBuilder {
    std::vector<std::map<int, double>> rows;
    explicit RowBuilder(int n) : rows(static_cast<std::size_t>(n)) {}
    void add(int r, int col, double v) {
        if (v != 0.0) rows[static_cast<std::size_t>(r)][col] += v;
    }
    ShellTransform build(int l) const {
        ShellTransform t;
        t.l = l;
        t.rows.resize(rows.size());
        for (std::size_t r = 0; r < rows.size(); ++r)
            for (const auto& e : rows[r])
                if (e.second != 0.0) t.rows[r].emplace_back(e.first, e.second);
        return t;
    }
};

// ---------------------------------------------------------------------------------------
// one-electron blocks
// ---------------------------------------------------------------------------------------

enum class OneKind { Overlap, Kinetic, Nuclear, Dipole };

// Transformed block [TA rows][TB rows] (Dipole: [3][TA rows][TB rows]) of two (component) shells.
void one_block(OneKind kind, const Shell& A, const ShellTransform& TA, const Shell& B, const ShellTransform& TB,
               const std::vector<Atom>* nuclei, const std::array<double, 3>* origin, std::vector<double>& cart,
               std::vector<double>& tmp, std::vector<double>& out) {
    const int nca = ncart(A.l), ra = TA.size(), rb = TB.size();
    int nblock = 1;
    switch (kind) {
        case OneKind::Overlap: overlap_cartesian(A, B, cart); break;
        case OneKind::Kinetic: kinetic_cartesian(A, B, cart); break;
        case OneKind::Nuclear: nuclear_cartesian(A, B, *nuclei, cart); break;
        case OneKind::Dipole:
            dipole_cartesian(A, B, *origin, cart);
            nblock = 3;
            break;
    }
    tmp.resize(static_cast<std::size_t>(nblock) * nca * rb);
    out.resize(static_cast<std::size_t>(nblock) * ra * rb);
    apply_transform(TB, nblock * nca, 1, cart.data(), tmp.data());
    for (int k = 0; k < nblock; ++k)
        apply_transform(TA, 1, rb, tmp.data() + static_cast<std::size_t>(k) * nca * rb,
                        out.data() + static_cast<std::size_t>(k) * ra * rb);
}

// Sum over the component pairs of two shells (each a list of components, or the shell itself).
struct ShellView {
    const Shell* shell = nullptr;
    const std::vector<DerivativeComponent>* comps = nullptr;  // nullptr: the plain shell
};

void one_block_sum(OneKind kind, const ShellView& a, const ShellView& b, const std::vector<Atom>* nuclei,
                   const std::array<double, 3>* origin, std::vector<double>& acc, std::vector<double>& cart,
                   std::vector<double>& tmp, std::vector<double>& out) {
    acc.clear();
    const auto run = [&](const Shell& A, const ShellTransform& TA, const Shell& B, const ShellTransform& TB) {
        one_block(kind, A, TA, B, TB, nuclei, origin, cart, tmp, out);
        if (acc.empty())
            acc = out;
        else
            for (std::size_t i = 0; i < acc.size(); ++i) acc[i] += out[i];
    };
    const ShellTransform& sa = spherical_transform(a.shell->l);
    const ShellTransform& sb = spherical_transform(b.shell->l);
    if (!a.comps && !b.comps) {
        run(*a.shell, sa, *b.shell, sb);
    } else if (a.comps && !b.comps) {
        for (const auto& c : *a.comps) run(c.shell, c.transform, *b.shell, sb);
    } else if (!a.comps && b.comps) {
        for (const auto& c : *b.comps) run(*a.shell, sa, c.shell, c.transform);
    } else {
        for (const auto& ca : *a.comps)
            for (const auto& cb : *b.comps) run(ca.shell, ca.transform, cb.shell, cb.transform);
    }
}

struct ShellComponents {
    std::vector<std::vector<DerivativeComponent>> d1, d2;
    ShellComponents(const BasisSet& basis, int order) {
        d1.resize(basis.shells.size());
        if (order >= 2) d2.resize(basis.shells.size());
        for (std::size_t s = 0; s < basis.shells.size(); ++s) {
            d1[s] = shell_derivative_components(basis.shells[s], 1);
            if (order >= 2) d2[s] = shell_derivative_components(basis.shells[s], 2);
        }
    }
};

std::vector<Atom> single_nucleus(const Atom& a) { return {a}; }

// Adds to the per-thread result through `sink(i, j, kind-specific data)`; loops over all
// ordered shell pairs (i, j) in parallel over i.
template <class Body>
void for_ordered_pairs(const BasisSet& basis, int threads, Body body) {
    const std::size_t ns = basis.shells.size();
    std::atomic<std::size_t> next{0};
    run_parallel(resolve_threads(threads), [&](int t) {
        for (std::size_t i; (i = next++) < ns;) {
            check_cancelled();
            for (std::size_t j = 0; j < ns; ++j) body(t, i, j);
        }
    });
}

}  // namespace

// ---------------------------------------------------------------------------------------
// derivative components
// ---------------------------------------------------------------------------------------

std::vector<DerivativeComponent> shell_derivative_components(const Shell& s, int order) {
    const int l = s.l, ns = 2 * l + 1;
    const std::vector<double>& M = solid_harmonic_matrix(l);
    const auto comps = cartesian_components(l);
    const int nc = static_cast<int>(comps.size());
    std::vector<DerivativeComponent> out;
    if (order == 1) {
        RowBuilder up(3 * ns), down(3 * ns);
        for (int k = 0; k < 3; ++k)
            for (int m = 0; m < ns; ++m)
                for (int c = 0; c < nc; ++c) {
                    const double v = M[static_cast<std::size_t>(m) * nc + c];
                    if (v == 0.0) continue;
                    auto e = comps[static_cast<std::size_t>(c)];
                    auto eu = e;
                    ++eu[k];
                    up.add(k * ns + m, cart_index(eu[0], eu[1], eu[2]), v);
                    if (e[k] > 0) {
                        auto ed = e;
                        --ed[k];
                        down.add(k * ns + m, cart_index(ed[0], ed[1], ed[2]), -e[k] * v);
                    }
                }
        out.push_back({component_shell(s, +1, 1), up.build(l + 1)});
        if (l >= 1) out.push_back({component_shell(s, -1, 0), down.build(l - 1)});
        return out;
    }
    if (order != 2) throw std::invalid_argument("shell_derivative_components: order must be 1 or 2");
    RowBuilder up2(6 * ns), mid(6 * ns), down2(6 * ns);
    for (int k = 0; k < 3; ++k)
        for (int kp = k; kp < 3; ++kp) {
            const int pr = direction_pair(k, kp);
            for (int m = 0; m < ns; ++m)
                for (int c = 0; c < nc; ++c) {
                    const double v = M[static_cast<std::size_t>(m) * nc + c];
                    if (v == 0.0) continue;
                    const auto e = comps[static_cast<std::size_t>(c)];
                    const int row = pr * ns + m;
                    auto idx = [](const std::array<int, 3>& x) { return cart_index(x[0], x[1], x[2]); };
                    if (k == kp) {
                        auto e2 = e;
                        e2[k] += 2;
                        up2.add(row, idx(e2), v);
                        mid.add(row, idx(e), -(2.0 * e[k] + 1.0) * v);
                        if (e[k] >= 2) {
                            auto em = e;
                            em[k] -= 2;
                            down2.add(row, idx(em), e[k] * (e[k] - 1.0) * v);
                        }
                    } else {
                        auto e2 = e;
                        ++e2[k];
                        ++e2[kp];
                        up2.add(row, idx(e2), v);
                        if (e[kp] > 0) {
                            auto x = e;
                            ++x[k];
                            --x[kp];
                            mid.add(row, idx(x), -e[kp] * v);
                        }
                        if (e[k] > 0) {
                            auto x = e;
                            --x[k];
                            ++x[kp];
                            mid.add(row, idx(x), -e[k] * v);
                        }
                        if (e[k] > 0 && e[kp] > 0) {
                            auto x = e;
                            --x[k];
                            --x[kp];
                            down2.add(row, idx(x), static_cast<double>(e[k]) * e[kp] * v);
                        }
                    }
                }
        }
    out.push_back({component_shell(s, +2, 2), up2.build(l + 2)});
    out.push_back({component_shell(s, 0, 1), mid.build(l)});
    if (l >= 2) out.push_back({component_shell(s, -2, 0), down2.build(l - 2)});
    return out;
}

// ---------------------------------------------------------------------------------------
// one-electron derivative matrices and contractions
// ---------------------------------------------------------------------------------------

namespace {

// First derivatives with respect to the centre of the first shell: block [3][sa][sb]
// (Dipole: [3 directions][3 components][sa][sb]).
void first_a(OneKind kind, const BasisSet& basis, const ShellComponents& sc, std::size_t i, std::size_t j,
             const std::vector<Atom>* nuclei, const std::array<double, 3>* origin, std::vector<double>& acc,
             std::vector<double>& cart, std::vector<double>& tmp, std::vector<double>& out) {
    one_block_sum(kind, ShellView{&basis.shells[i], &sc.d1[i]}, ShellView{&basis.shells[j], nullptr}, nuclei, origin,
                  acc, cart, tmp, out);
}

}  // namespace

std::vector<Matrix> overlap_derivative_matrices(const BasisSet& basis, int n_atoms, int threads) {
    const ShellComponents sc(basis, 1);
    const std::size_t N = static_cast<std::size_t>(basis.n_functions);
    std::vector<Matrix> out(static_cast<std::size_t>(3 * n_atoms), Matrix(N, N));
    // Each ordered pair (i, j) writes only the elements (mu in i, nu in j) and (nu, mu) of the
    // atom of i; the (j, i) iteration supplies the derivative with respect to the atom of j.
    // Rows of different i never overlap within an iteration over i, but the transposed
    // elements do, so accumulate per thread.
    const int nt = resolve_threads(threads);
    std::vector<std::vector<Matrix>> parts(static_cast<std::size_t>(nt), out);
    for_ordered_pairs(basis, nt, [&](int t, std::size_t i, std::size_t j) {
        thread_local std::vector<double> acc, cart, tmp, buf;
        first_a(OneKind::Overlap, basis, sc, i, j, nullptr, nullptr, acc, cart, tmp, buf);
        const Shell& A = basis.shells[i];
        const Shell& B = basis.shells[j];
        const int sa = A.size(), sb = B.size();
        for (int k = 0; k < 3; ++k) {
            Matrix& M = parts[static_cast<std::size_t>(t)][static_cast<std::size_t>(3 * A.atom + k)];
            for (int a = 0; a < sa; ++a)
                for (int b = 0; b < sb; ++b) {
                    const double v = acc[(static_cast<std::size_t>(k) * sa + a) * sb + b];
                    M(A.first_function + a, B.first_function + b) += v;
                    M(B.first_function + b, A.first_function + a) += v;
                }
        }
    });
    for (auto& p : parts)
        for (std::size_t x = 0; x < out.size(); ++x)
            for (std::size_t r = 0; r < N; ++r)
                for (std::size_t c = 0; c < N; ++c) out[x](r, c) += p[x](r, c);
    return out;
}

std::vector<Matrix> core_hamiltonian_derivative_matrices(const BasisSet& basis, const Molecule& mol, int threads) {
    const int n_atoms = static_cast<int>(mol.atoms.size());
    const ShellComponents sc(basis, 1);
    const std::size_t N = static_cast<std::size_t>(basis.n_functions);
    std::vector<Matrix> out(static_cast<std::size_t>(3 * n_atoms), Matrix(N, N));
    const int nt = resolve_threads(threads);
    std::vector<std::vector<Matrix>> parts(static_cast<std::size_t>(nt), out);
    for_ordered_pairs(basis, nt, [&](int t, std::size_t i, std::size_t j) {
        thread_local std::vector<double> acc, cart, tmp, buf;
        const Shell& A = basis.shells[i];
        const Shell& B = basis.shells[j];
        const int sa = A.size(), sb = B.size();
        auto& P = parts[static_cast<std::size_t>(t)];
        auto add = [&](int atom, double sign) {
            for (int k = 0; k < 3; ++k) {
                Matrix& M = P[static_cast<std::size_t>(3 * atom + k)];
                for (int a = 0; a < sa; ++a)
                    for (int b = 0; b < sb; ++b) {
                        const double v = sign * acc[(static_cast<std::size_t>(k) * sa + a) * sb + b];
                        M(A.first_function + a, B.first_function + b) += v;
                        M(B.first_function + b, A.first_function + a) += v;
                    }
            }
        };
        first_a(OneKind::Kinetic, basis, sc, i, j, nullptr, nullptr, acc, cart, tmp, buf);
        add(A.atom, 1.0);
        for (int c = 0; c < n_atoms; ++c) {
            const std::vector<Atom> nuc = single_nucleus(mol.atoms[static_cast<std::size_t>(c)]);
            first_a(OneKind::Nuclear, basis, sc, i, j, &nuc, nullptr, acc, cart, tmp, buf);
            add(A.atom, 1.0);
            add(c, -1.0);  // operator: d/dC = -(d/dA + d/dB); the (j, i) iteration supplies -d/dB
        }
    });
    for (auto& p : parts)
        for (std::size_t x = 0; x < out.size(); ++x)
            for (std::size_t r = 0; r < N; ++r)
                for (std::size_t c = 0; c < N; ++c) out[x](r, c) += p[x](r, c);
    return out;
}

std::vector<double> overlap_gradient(const BasisSet& basis, int n_atoms, const Matrix& D, int threads) {
    const ShellComponents sc(basis, 1);
    const int nt = resolve_threads(threads);
    std::vector<std::vector<double>> parts(static_cast<std::size_t>(nt), std::vector<double>(3 * n_atoms, 0.0));
    for_ordered_pairs(basis, nt, [&](int t, std::size_t i, std::size_t j) {
        thread_local std::vector<double> acc, cart, tmp, buf;
        first_a(OneKind::Overlap, basis, sc, i, j, nullptr, nullptr, acc, cart, tmp, buf);
        const Shell& A = basis.shells[i];
        const Shell& B = basis.shells[j];
        const int sa = A.size(), sb = B.size();
        for (int k = 0; k < 3; ++k) {
            double s = 0.0;
            for (int a = 0; a < sa; ++a)
                for (int b = 0; b < sb; ++b)
                    s += D(A.first_function + a, B.first_function + b) *
                         acc[(static_cast<std::size_t>(k) * sa + a) * sb + b];
            parts[static_cast<std::size_t>(t)][static_cast<std::size_t>(3 * A.atom + k)] += 2.0 * s;
        }
    });
    std::vector<double> g(3 * n_atoms, 0.0);
    for (const auto& p : parts)
        for (std::size_t x = 0; x < g.size(); ++x) g[x] += p[x];
    return g;
}

std::vector<double> core_hamiltonian_gradient(const BasisSet& basis, const Molecule& mol, const Matrix& D,
                                              int threads) {
    const int n_atoms = static_cast<int>(mol.atoms.size());
    const ShellComponents sc(basis, 1);
    const int nt = resolve_threads(threads);
    std::vector<std::vector<double>> parts(static_cast<std::size_t>(nt), std::vector<double>(3 * n_atoms, 0.0));
    for_ordered_pairs(basis, nt, [&](int t, std::size_t i, std::size_t j) {
        thread_local std::vector<double> acc, cart, tmp, buf;
        const Shell& A = basis.shells[i];
        const Shell& B = basis.shells[j];
        const int sa = A.size(), sb = B.size();
        auto& g = parts[static_cast<std::size_t>(t)];
        auto contract = [&](int k) {
            double s = 0.0;
            for (int a = 0; a < sa; ++a)
                for (int b = 0; b < sb; ++b)
                    s += D(A.first_function + a, B.first_function + b) *
                         acc[(static_cast<std::size_t>(k) * sa + a) * sb + b];
            return s;
        };
        first_a(OneKind::Kinetic, basis, sc, i, j, nullptr, nullptr, acc, cart, tmp, buf);
        for (int k = 0; k < 3; ++k) g[static_cast<std::size_t>(3 * A.atom + k)] += 2.0 * contract(k);
        for (int c = 0; c < n_atoms; ++c) {
            const std::vector<Atom> nuc = single_nucleus(mol.atoms[static_cast<std::size_t>(c)]);
            first_a(OneKind::Nuclear, basis, sc, i, j, &nuc, nullptr, acc, cart, tmp, buf);
            for (int k = 0; k < 3; ++k) {
                const double s = 2.0 * contract(k);
                g[static_cast<std::size_t>(3 * A.atom + k)] += s;
                g[static_cast<std::size_t>(3 * c + k)] -= s;
            }
        }
    });
    std::vector<double> g(3 * n_atoms, 0.0);
    for (const auto& p : parts)
        for (std::size_t x = 0; x < g.size(); ++x) g[x] += p[x];
    return g;
}

namespace {

// Second derivatives of the element block (mu in i, nu in j) of a one-electron operator with
// respect to the centres A (of i) and B (of j), contracted with D: blocks
//   AA[k][k'], AB[k][k'] (d/dA_k d/dB_k'), BB[k][k'].
struct PairHessian {
    double AA[3][3] = {}, AB[3][3] = {}, BB[3][3] = {};
};

PairHessian pair_hessian(OneKind kind, const BasisSet& basis, const ShellComponents& sc, std::size_t i,
                         std::size_t j, const std::vector<Atom>* nuclei, const Matrix& D, std::vector<double>& acc,
                         std::vector<double>& cart, std::vector<double>& tmp, std::vector<double>& out) {
    PairHessian h;
    const Shell& A = basis.shells[i];
    const Shell& B = basis.shells[j];
    const int sa = A.size(), sb = B.size();
    auto dm = [&](int a, int b) { return D(A.first_function + a, B.first_function + b); };
    // AA: [6][sa][sb]
    one_block_sum(kind, ShellView{&A, &sc.d2[i]}, ShellView{&B, nullptr}, nuclei, nullptr, acc, cart, tmp, out);
    for (int k = 0; k < 3; ++k)
        for (int kp = 0; kp < 3; ++kp) {
            const int pr = direction_pair(k, kp);
            double s = 0.0;
            for (int a = 0; a < sa; ++a)
                for (int b = 0; b < sb; ++b) s += dm(a, b) * acc[(static_cast<std::size_t>(pr) * sa + a) * sb + b];
            h.AA[k][kp] = s;
        }
    // BB: [sa][6][sb]
    one_block_sum(kind, ShellView{&A, nullptr}, ShellView{&B, &sc.d2[j]}, nuclei, nullptr, acc, cart, tmp, out);
    for (int k = 0; k < 3; ++k)
        for (int kp = 0; kp < 3; ++kp) {
            const int pr = direction_pair(k, kp);
            double s = 0.0;
            for (int a = 0; a < sa; ++a)
                for (int b = 0; b < sb; ++b) s += dm(a, b) * acc[(static_cast<std::size_t>(a) * 6 + pr) * sb + b];
            h.BB[k][kp] = s;
        }
    // AB: [3][sa][3][sb]
    one_block_sum(kind, ShellView{&A, &sc.d1[i]}, ShellView{&B, &sc.d1[j]}, nuclei, nullptr, acc, cart, tmp, out);
    for (int k = 0; k < 3; ++k)
        for (int kp = 0; kp < 3; ++kp) {
            double s = 0.0;
            for (int a = 0; a < sa; ++a)
                for (int b = 0; b < sb; ++b)
                    s += dm(a, b) * acc[((static_cast<std::size_t>(k) * sa + a) * 3 + kp) * sb + b];
            h.AB[k][kp] = s;
        }
    return h;
}

// Adds a 3x3 block (rows atom p, columns atom q) with factor f.
void add_block(Matrix& H, int p, int q, const double b[3][3], double f, bool transpose = false) {
    for (int k = 0; k < 3; ++k)
        for (int kp = 0; kp < 3; ++kp) H(3 * p + k, 3 * q + kp) += f * (transpose ? b[kp][k] : b[k][kp]);
}

Matrix one_electron_hessian(OneKind kind, const BasisSet& basis, const Molecule* mol, int n_atoms, const Matrix& D,
                            int threads) {
    const ShellComponents sc(basis, 2);
    const int nt = resolve_threads(threads);
    const std::size_t n3 = static_cast<std::size_t>(3 * n_atoms);
    std::vector<Matrix> parts(static_cast<std::size_t>(nt), Matrix(n3, n3));
    for_ordered_pairs(basis, nt, [&](int t, std::size_t i, std::size_t j) {
        thread_local std::vector<double> acc, cart, tmp, buf;
        Matrix& H = parts[static_cast<std::size_t>(t)];
        const int a = basis.shells[i].atom, b = basis.shells[j].atom;
        auto distribute = [&](const PairHessian& h, int c) {
            // centres A (atom a), B (atom b) and, for the nuclear attraction, C (atom c >= 0);
            // d/dC = -(d/dA + d/dB).
            add_block(H, a, a, h.AA, 1.0);
            add_block(H, b, b, h.BB, 1.0);
            add_block(H, a, b, h.AB, 1.0);
            add_block(H, b, a, h.AB, 1.0, true);
            if (c < 0) return;
            double AC[3][3], BC[3][3], CC[3][3];
            for (int k = 0; k < 3; ++k)
                for (int kp = 0; kp < 3; ++kp) {
                    AC[k][kp] = -(h.AA[k][kp] + h.AB[k][kp]);
                    BC[k][kp] = -(h.AB[kp][k] + h.BB[k][kp]);
                    CC[k][kp] = h.AA[k][kp] + h.AB[k][kp] + h.AB[kp][k] + h.BB[k][kp];
                }
            add_block(H, a, c, AC, 1.0);
            add_block(H, c, a, AC, 1.0, true);
            add_block(H, b, c, BC, 1.0);
            add_block(H, c, b, BC, 1.0, true);
            add_block(H, c, c, CC, 1.0);
        };
        if (kind == OneKind::Nuclear) {
            for (int c = 0; c < n_atoms; ++c) {
                const std::vector<Atom> nuc = single_nucleus(mol->atoms[static_cast<std::size_t>(c)]);
                distribute(pair_hessian(kind, basis, sc, i, j, &nuc, D, acc, cart, tmp, buf), c);
            }
        } else {
            distribute(pair_hessian(kind, basis, sc, i, j, nullptr, D, acc, cart, tmp, buf), -1);
        }
    });
    Matrix H(n3, n3);
    for (const auto& p : parts)
        for (std::size_t r = 0; r < n3; ++r)
            for (std::size_t c = 0; c < n3; ++c) H(r, c) += p(r, c);
    return H;
}

}  // namespace

Matrix overlap_hessian(const BasisSet& basis, int n_atoms, const Matrix& D, int threads) {
    return one_electron_hessian(OneKind::Overlap, basis, nullptr, n_atoms, D, threads);
}

Matrix core_hamiltonian_hessian(const BasisSet& basis, const Molecule& mol, const Matrix& D, int threads) {
    const int n = static_cast<int>(mol.atoms.size());
    Matrix H = one_electron_hessian(OneKind::Kinetic, basis, nullptr, n, D, threads);
    const Matrix V = one_electron_hessian(OneKind::Nuclear, basis, &mol, n, D, threads);
    for (std::size_t r = 0; r < H.rows(); ++r)
        for (std::size_t c = 0; c < H.cols(); ++c) H(r, c) += V(r, c);
    return H;
}

std::vector<std::vector<Matrix>> dipole_derivative_matrices(const BasisSet& basis, int n_atoms,
                                                            const std::array<double, 3>& origin, int threads) {
    const ShellComponents sc(basis, 1);
    const std::size_t N = static_cast<std::size_t>(basis.n_functions);
    std::vector<std::vector<Matrix>> out(static_cast<std::size_t>(3 * n_atoms), std::vector<Matrix>(3, Matrix(N, N)));
    const int nt = resolve_threads(threads);
    std::vector<std::vector<std::vector<Matrix>>> parts(static_cast<std::size_t>(nt), out);
    for_ordered_pairs(basis, nt, [&](int t, std::size_t i, std::size_t j) {
        thread_local std::vector<double> acc, cart, tmp, buf;
        const Shell& A = basis.shells[i];
        const Shell& B = basis.shells[j];
        const int sa = A.size(), sb = B.size();
        // Dipole blocks of the derivative components: [3 moments][3 dirs][sa][sb]
        // (one_block puts the moment index outermost).
        first_a(OneKind::Dipole, basis, sc, i, j, nullptr, &origin, acc, cart, tmp, buf);
        for (int k = 0; k < 3; ++k)
            for (int m = 0; m < 3; ++m) {
                Matrix& M = parts[static_cast<std::size_t>(t)][static_cast<std::size_t>(3 * A.atom + k)]
                                 [static_cast<std::size_t>(m)];
                for (int a = 0; a < sa; ++a)
                    for (int b = 0; b < sb; ++b) {
                        const double v = acc[((static_cast<std::size_t>(m) * 3 + k) * sa + a) * sb + b];
                        M(A.first_function + a, B.first_function + b) += v;
                        M(B.first_function + b, A.first_function + a) += v;
                    }
            }
    });
    for (auto& p : parts)
        for (std::size_t x = 0; x < out.size(); ++x)
            for (std::size_t m = 0; m < 3; ++m)
                for (std::size_t r = 0; r < N; ++r)
                    for (std::size_t c = 0; c < N; ++c) out[x][m](r, c) += p[x][m](r, c);
    return out;
}

// ---------------------------------------------------------------------------------------
// two-electron derivatives
// ---------------------------------------------------------------------------------------

namespace {

constexpr double kDerivativeScreen = 1.0e-15;  // Schwarz product below which a quartet is skipped

struct UniquePair {
    int i = 0, j = 0;  // shells, i >= j
    std::vector<ShellPair> p00;                      // the pair itself (one element)
    std::vector<ShellPair> p10, p01, p20, p02, p11;  // component pairs (order: A comps x B comps)
};

struct DerivativePairs {
    ShellComponents comps;
    std::vector<UniquePair> pairs;  // index ShellPairList::index(i, j)
    std::vector<double> schwarz;

    DerivativePairs(const BasisSet& basis, int order, int threads) : comps(basis, order) {
        const int ns = static_cast<int>(basis.shells.size());
        pairs.resize(static_cast<std::size_t>(ns) * (ns + 1) / 2);
        schwarz.assign(pairs.size(), 0.0);
        for (int i = 0; i < ns; ++i)
            for (int j = 0; j <= i; ++j) {
                UniquePair& u = pairs[ShellPairList::index(i, j)];
                u.i = i;
                u.j = j;
                const Shell& A = basis.shells[static_cast<std::size_t>(i)];
                const Shell& B = basis.shells[static_cast<std::size_t>(j)];
                u.p00.push_back(make_shell_pair(A, nullptr, B, nullptr, true, i, j));
                if (order >= 0) {
                    for (const auto& c : comps.d1[static_cast<std::size_t>(i)])
                        u.p10.push_back(make_shell_pair(c.shell, &c.transform, B, nullptr, true, i, j));
                    for (const auto& c : comps.d1[static_cast<std::size_t>(j)])
                        u.p01.push_back(make_shell_pair(A, nullptr, c.shell, &c.transform, true, i, j));
                }
                if (order >= 2) {
                    for (const auto& c : comps.d2[static_cast<std::size_t>(i)])
                        u.p20.push_back(make_shell_pair(c.shell, &c.transform, B, nullptr, true, i, j));
                    for (const auto& c : comps.d2[static_cast<std::size_t>(j)])
                        u.p02.push_back(make_shell_pair(A, nullptr, c.shell, &c.transform, true, i, j));
                    for (const auto& ca : comps.d1[static_cast<std::size_t>(i)])
                        for (const auto& cb : comps.d1[static_cast<std::size_t>(j)])
                            u.p11.push_back(make_shell_pair(ca.shell, &ca.transform, cb.shell, &cb.transform, true, i, j));
                }
            }
        // Schwarz bounds of the undifferentiated pairs.
        std::atomic<std::size_t> next{0};
        run_parallel(resolve_threads(threads), [&](int) {
            EriEngine engine;
            for (std::size_t p; (p = next++) < pairs.size();) {
                const ShellPair& sp = pairs[p].p00[0];
                if (sp.prims.empty()) continue;
                const double* v = engine.compute(sp, sp);
                const std::size_t n = static_cast<std::size_t>(sp.rows_a()) * sp.rows_b();
                double mx = 0.0;
                for (std::size_t ij = 0; ij < n; ++ij) mx = std::max(mx, std::fabs(v[ij * n + ij]));
                schwarz[p] = std::sqrt(mx);
            }
        });
    }
};

// Gamma block [a][b][c][d] for the shells of a quartet.
void gamma_block(const BasisSet& basis, const TwoParticleDensity& g, int si, int sj, int sk, int sl,
                 std::vector<double>& out) {
    const Shell& A = basis.shells[static_cast<std::size_t>(si)];
    const Shell& B = basis.shells[static_cast<std::size_t>(sj)];
    const Shell& C = basis.shells[static_cast<std::size_t>(sk)];
    const Shell& D = basis.shells[static_cast<std::size_t>(sl)];
    const int na = A.size(), nb = B.size(), nc = C.size(), nd = D.size();
    out.assign(static_cast<std::size_t>(na) * nb * nc * nd, 0.0);
    std::size_t x = 0;
    for (int a = 0; a < na; ++a) {
        const std::size_t m = static_cast<std::size_t>(A.first_function + a);
        for (int b = 0; b < nb; ++b) {
            const std::size_t n = static_cast<std::size_t>(B.first_function + b);
            for (int c = 0; c < nc; ++c) {
                const std::size_t l = static_cast<std::size_t>(C.first_function + c);
                for (int d = 0; d < nd; ++d, ++x) {
                    const std::size_t s = static_cast<std::size_t>(D.first_function + d);
                    double v = 0.0;
                    for (const SeparableDensity& t : g.separable) {
                        if (t.coulomb != 0.0) v += 0.5 * t.coulomb * (t.A(m, n) * t.B(l, s) + t.B(m, n) * t.A(l, s));
                        if (t.exchange != 0.0)
                            v += 0.25 * t.exchange *
                                 (t.A(m, l) * t.B(n, s) + t.A(m, s) * t.B(n, l) + t.B(m, l) * t.A(n, s) +
                                  t.B(m, s) * t.A(n, l));
                    }
                    if (g.packed) {
                        const std::size_t mn = FockBuilder::pair_index(m, n), ls = FockBuilder::pair_index(l, s);
                        v += (*g.packed)[FockBuilder::pair_index(mn, ls)];
                    }
                    out[x] = v;
                }
            }
        }
    }
}

// Contracts I (layout dims with one direction index of size nd inserted before dims[pos])
// with Gamma (layout dims): r[k] += f * sum Gamma I.
void contract1(const double* I, const std::vector<double>& G, const int dims[4], int pos, int nd, double* r) {
    std::size_t pre = 1, post = 1;
    for (int q = 0; q < pos; ++q) pre *= static_cast<std::size_t>(dims[q]);
    for (int q = pos; q < 4; ++q) post *= static_cast<std::size_t>(dims[q]);
    for (int k = 0; k < nd; ++k) {
        double s = 0.0;
        for (std::size_t p = 0; p < pre; ++p) {
            const double* g = G.data() + p * post;
            const double* v = I + (p * nd + static_cast<std::size_t>(k)) * post;
            for (std::size_t q = 0; q < post; ++q) s += g[q] * v[q];
        }
        r[k] += s;
    }
}

// Two direction indices (3 each) inserted before dims[pos1] and dims[pos2] (pos1 < pos2).
void contract2(const double* I, const std::vector<double>& G, const int dims[4], int pos1, int pos2,
               double r[3][3]) {
    std::size_t n0 = 1, n1 = 1, n2 = 1;
    for (int q = 0; q < pos1; ++q) n0 *= static_cast<std::size_t>(dims[q]);
    for (int q = pos1; q < pos2; ++q) n1 *= static_cast<std::size_t>(dims[q]);
    for (int q = pos2; q < 4; ++q) n2 *= static_cast<std::size_t>(dims[q]);
    for (int k = 0; k < 3; ++k)
        for (int kp = 0; kp < 3; ++kp) {
            double s = 0.0;
            for (std::size_t i0 = 0; i0 < n0; ++i0)
                for (std::size_t i1 = 0; i1 < n1; ++i1) {
                    const double* g = G.data() + (i0 * n1 + i1) * n2;
                    const double* v = I + ((((i0 * 3 + static_cast<std::size_t>(k)) * n1 + i1) * 3 +
                                            static_cast<std::size_t>(kp)) *
                                           n2);
                    for (std::size_t i2 = 0; i2 < n2; ++i2) s += g[i2] * v[i2];
                }
            r[k][kp] += s;
        }
}

// Sum of the engine results over component pairs.
template <class Fn>
void sum_components(EriEngine& engine, const std::vector<ShellPair>& bras, const std::vector<ShellPair>& kets,
                    std::vector<double>& acc, std::size_t size, Fn&& use) {
    acc.assign(size, 0.0);
    for (const ShellPair& b : bras) {
        if (b.prims.empty()) continue;
        for (const ShellPair& k : kets) {
            if (k.prims.empty()) continue;
            const double* v = engine.compute(b, k);
            for (std::size_t x = 0; x < size; ++x) acc[x] += v[x];
        }
    }
    use(acc.data());
}

template <class Body>
void for_unique_quartets(const DerivativePairs& dp, int threads, Body body) {
    const std::size_t np = dp.pairs.size();
    std::atomic<std::size_t> next{0};
    run_parallel(threads, [&](int t) {
        EriEngine engine;
        for (std::size_t k; (k = next++) < np;) {
            check_cancelled();
            const std::size_t P = np - 1 - k;  // heaviest rows first
            if (dp.pairs[P].p00[0].prims.empty()) continue;
            for (std::size_t Q = 0; Q <= P; ++Q) {
                if (dp.pairs[Q].p00[0].prims.empty()) continue;
                if (dp.schwarz[P] * dp.schwarz[Q] < kDerivativeScreen) continue;
                body(t, engine, dp.pairs[P], dp.pairs[Q], P == Q);
            }
        }
    });
}

}  // namespace

double eri_contract(const BasisSet& basis, const TwoParticleDensity& gamma, int threads) {
    const int nt = resolve_threads(threads);
    const DerivativePairs dp(basis, -1, nt);
    std::vector<double> parts(static_cast<std::size_t>(nt), 0.0);
    for_unique_quartets(dp, nt, [&](int t, EriEngine& engine, const UniquePair& bra, const UniquePair& ket, bool same) {
        thread_local std::vector<double> G;
        gamma_block(basis, gamma, bra.i, bra.j, ket.i, ket.j, G);
        const double deg = (bra.i == bra.j ? 1.0 : 2.0) * (ket.i == ket.j ? 1.0 : 2.0) * (same ? 1.0 : 2.0);
        const double* v = engine.compute(bra.p00[0], ket.p00[0]);
        double s = 0.0;
        for (std::size_t x = 0; x < G.size(); ++x) s += G[x] * v[x];
        parts[static_cast<std::size_t>(t)] += deg * s;
    });
    double e = 0.0;
    for (double p : parts) e += p;
    return e;
}

std::vector<double> eri_gradient(const BasisSet& basis, int n_atoms, const TwoParticleDensity& gamma, int threads) {
    const int nt = resolve_threads(threads);
    const DerivativePairs dp(basis, 1, nt);
    std::vector<std::vector<double>> parts(static_cast<std::size_t>(nt), std::vector<double>(3 * n_atoms, 0.0));
    for_unique_quartets(dp, nt, [&](int t, EriEngine& engine, const UniquePair& bra, const UniquePair& ket, bool same) {
        thread_local std::vector<double> G, acc;
        gamma_block(basis, gamma, bra.i, bra.j, ket.i, ket.j, G);
        const Shell& A = basis.shells[static_cast<std::size_t>(bra.i)];
        const Shell& B = basis.shells[static_cast<std::size_t>(bra.j)];
        const Shell& C = basis.shells[static_cast<std::size_t>(ket.i)];
        const Shell& D = basis.shells[static_cast<std::size_t>(ket.j)];
        const int dims[4] = {A.size(), B.size(), C.size(), D.size()};
        const std::size_t size = 3 * G.size();
        double gA[3] = {0, 0, 0}, gB[3] = {0, 0, 0}, gC[3] = {0, 0, 0};
        const std::vector<ShellPair>& k00 = ket.p00;
        const std::vector<ShellPair>& b00 = bra.p00;
        sum_components(engine, bra.p10, k00, acc, size, [&](const double* I) { contract1(I, G, dims, 0, 3, gA); });
        sum_components(engine, bra.p01, k00, acc, size, [&](const double* I) { contract1(I, G, dims, 1, 3, gB); });
        sum_components(engine, b00, ket.p10, acc, size, [&](const double* I) { contract1(I, G, dims, 2, 3, gC); });
        const double deg = (bra.i == bra.j ? 1.0 : 2.0) * (ket.i == ket.j ? 1.0 : 2.0) * (same ? 1.0 : 2.0);
        auto& g = parts[static_cast<std::size_t>(t)];
        for (int k = 0; k < 3; ++k) {
            g[static_cast<std::size_t>(3 * A.atom + k)] += deg * gA[k];
            g[static_cast<std::size_t>(3 * B.atom + k)] += deg * gB[k];
            g[static_cast<std::size_t>(3 * C.atom + k)] += deg * gC[k];
            g[static_cast<std::size_t>(3 * D.atom + k)] -= deg * (gA[k] + gB[k] + gC[k]);
        }
    });
    std::vector<double> g(3 * n_atoms, 0.0);
    for (const auto& p : parts)
        for (std::size_t x = 0; x < g.size(); ++x) g[x] += p[x];
    return g;
}

Matrix eri_hessian(const BasisSet& basis, int n_atoms, const TwoParticleDensity& gamma, int threads) {
    const int nt = resolve_threads(threads);
    const DerivativePairs dp(basis, 2, nt);
    const std::size_t n3 = static_cast<std::size_t>(3 * n_atoms);
    std::vector<Matrix> parts(static_cast<std::size_t>(nt), Matrix(n3, n3));
    for_unique_quartets(dp, nt, [&](int t, EriEngine& engine, const UniquePair& bra, const UniquePair& ket, bool same) {
        thread_local std::vector<double> G, acc;
        gamma_block(basis, gamma, bra.i, bra.j, ket.i, ket.j, G);
        const Shell& A = basis.shells[static_cast<std::size_t>(bra.i)];
        const Shell& B = basis.shells[static_cast<std::size_t>(bra.j)];
        const Shell& C = basis.shells[static_cast<std::size_t>(ket.i)];
        const Shell& D = basis.shells[static_cast<std::size_t>(ket.j)];
        const int dims[4] = {A.size(), B.size(), C.size(), D.size()};
        const std::vector<ShellPair>& k00 = ket.p00;
        const std::vector<ShellPair>& b00 = bra.p00;
        // Blocks H[X][Y][k][k'] over the centres X, Y = A, B, C (D by translational invariance).
        double H[4][4][3][3] = {};
        double r6[6];
        auto from_pairs = [&](double out[3][3]) {
            for (int k = 0; k < 3; ++k)
                for (int kp = 0; kp < 3; ++kp) out[k][kp] = r6[direction_pair(k, kp)];
        };
        const std::size_t s6 = 6 * G.size(), s9 = 9 * G.size();
        std::fill(r6, r6 + 6, 0.0);
        sum_components(engine, bra.p20, k00, acc, s6, [&](const double* I) { contract1(I, G, dims, 0, 6, r6); });
        from_pairs(H[0][0]);
        std::fill(r6, r6 + 6, 0.0);
        sum_components(engine, bra.p02, k00, acc, s6, [&](const double* I) { contract1(I, G, dims, 1, 6, r6); });
        from_pairs(H[1][1]);
        std::fill(r6, r6 + 6, 0.0);
        sum_components(engine, b00, ket.p20, acc, s6, [&](const double* I) { contract1(I, G, dims, 2, 6, r6); });
        from_pairs(H[2][2]);
        sum_components(engine, bra.p11, k00, acc, s9, [&](const double* I) { contract2(I, G, dims, 0, 1, H[0][1]); });
        sum_components(engine, bra.p10, ket.p10, acc, s9,
                       [&](const double* I) { contract2(I, G, dims, 0, 2, H[0][2]); });
        sum_components(engine, bra.p01, ket.p10, acc, s9,
                       [&](const double* I) { contract2(I, G, dims, 1, 2, H[1][2]); });
        for (int x = 0; x < 3; ++x)
            for (int y = 0; y < x; ++y)
                for (int k = 0; k < 3; ++k)
                    for (int kp = 0; kp < 3; ++kp) H[x][y][k][kp] = H[y][x][kp][k];
        for (int x = 0; x < 3; ++x)
            for (int k = 0; k < 3; ++k)
                for (int kp = 0; kp < 3; ++kp) {
                    double s = 0.0;
                    for (int y = 0; y < 3; ++y) s += H[x][y][k][kp];
                    H[x][3][k][kp] = -s;
                }
        for (int y = 0; y < 3; ++y)
            for (int k = 0; k < 3; ++k)
                for (int kp = 0; kp < 3; ++kp) H[3][y][k][kp] = H[y][3][kp][k];
        for (int k = 0; k < 3; ++k)
            for (int kp = 0; kp < 3; ++kp) {
                double s = 0.0;
                for (int x = 0; x < 3; ++x)
                    for (int y = 0; y < 3; ++y) s += H[x][y][k][kp];
                H[3][3][k][kp] = s;
            }
        const double deg = (bra.i == bra.j ? 1.0 : 2.0) * (ket.i == ket.j ? 1.0 : 2.0) * (same ? 1.0 : 2.0);
        const int atom[4] = {A.atom, B.atom, C.atom, D.atom};
        Matrix& M = parts[static_cast<std::size_t>(t)];
        for (int x = 0; x < 4; ++x)
            for (int y = 0; y < 4; ++y)
                for (int k = 0; k < 3; ++k)
                    for (int kp = 0; kp < 3; ++kp)
                        M(static_cast<std::size_t>(3 * atom[x] + k), static_cast<std::size_t>(3 * atom[y] + kp)) +=
                            deg * H[x][y][k][kp];
    });
    Matrix H(n3, n3);
    for (const auto& p : parts)
        for (std::size_t r = 0; r < n3; ++r)
            for (std::size_t c = 0; c < n3; ++c) H(r, c) += p(r, c);
    return H;
}

std::vector<Matrix> eri_derivative_fock(const BasisSet& basis, int n_atoms, const Matrix& Dm, int threads) {
    const int nt = resolve_threads(threads);
    const DerivativePairs dp(basis, 1, nt);
    const std::size_t N = static_cast<std::size_t>(basis.n_functions);
    const std::size_t nx = static_cast<std::size_t>(3 * n_atoms);
    std::vector<std::vector<Matrix>> parts(static_cast<std::size_t>(nt), std::vector<Matrix>(nx, Matrix(N, N)));
    for_unique_quartets(dp, nt, [&](int t, EriEngine& engine, const UniquePair& bra, const UniquePair& ket, bool same) {
        thread_local std::vector<double> acc, dX[4];
        const Shell& A = basis.shells[static_cast<std::size_t>(bra.i)];
        const Shell& B = basis.shells[static_cast<std::size_t>(bra.j)];
        const Shell& C = basis.shells[static_cast<std::size_t>(ket.i)];
        const Shell& D = basis.shells[static_cast<std::size_t>(ket.j)];
        const int na = A.size(), nb = B.size(), nc = C.size(), nd = D.size();
        const std::size_t nq = static_cast<std::size_t>(na) * nb * nc * nd;
        const std::vector<ShellPair>& k00 = ket.p00;
        const std::vector<ShellPair>& b00 = bra.p00;
        // dX[k][abcd] for X = A, B, C, D.
        const int dims[4] = {na, nb, nc, nd};
        auto to_k_outer = [&](const double* I, int pos, std::vector<double>& outv) {
            std::size_t pre = 1, post = 1;
            for (int q = 0; q < pos; ++q) pre *= static_cast<std::size_t>(dims[q]);
            for (int q = pos; q < 4; ++q) post *= static_cast<std::size_t>(dims[q]);
            outv.assign(3 * nq, 0.0);
            for (int k = 0; k < 3; ++k)
                for (std::size_t p = 0; p < pre; ++p)
                    for (std::size_t q = 0; q < post; ++q)
                        outv[static_cast<std::size_t>(k) * nq + p * post + q] =
                            I[(p * 3 + static_cast<std::size_t>(k)) * post + q];
        };
        sum_components(engine, bra.p10, k00, acc, 3 * nq, [&](const double* I) { to_k_outer(I, 0, dX[0]); });
        sum_components(engine, bra.p01, k00, acc, 3 * nq, [&](const double* I) { to_k_outer(I, 1, dX[1]); });
        sum_components(engine, b00, ket.p10, acc, 3 * nq, [&](const double* I) { to_k_outer(I, 2, dX[2]); });
        dX[3].assign(3 * nq, 0.0);
        for (std::size_t x = 0; x < 3 * nq; ++x) dX[3][x] = -(dX[0][x] + dX[1][x] + dX[2][x]);
        const double deg = (bra.i == bra.j ? 1.0 : 2.0) * (ket.i == ket.j ? 1.0 : 2.0) * (same ? 1.0 : 2.0);
        const int atom[4] = {A.atom, B.atom, C.atom, D.atom};
        auto& F = parts[static_cast<std::size_t>(t)];
        for (int X = 0; X < 4; ++X) {
            bool first = true;
            for (int Y = 0; Y < X; ++Y)
                if (atom[Y] == atom[X]) first = false;
            if (!first) continue;  // atom already handled
            for (int k = 0; k < 3; ++k) {
                // derivative with respect to coordinate k of atom[X]: sum of the centres on it
                thread_local std::vector<double> I;
                I.assign(nq, 0.0);
                for (int Y = X; Y < 4; ++Y)
                    if (atom[Y] == atom[X])
                        for (std::size_t q = 0; q < nq; ++q) I[q] += dX[Y][static_cast<std::size_t>(k) * nq + q];
                Matrix& G = F[static_cast<std::size_t>(3 * atom[X] + k)];
                const double* v = I.data();
                for (int a = 0; a < na; ++a) {
                    const std::size_t p = static_cast<std::size_t>(A.first_function + a);
                    for (int b = 0; b < nb; ++b) {
                        const std::size_t q = static_cast<std::size_t>(B.first_function + b);
                        const double dpq = Dm(p, q);
                        double gpq = 0.0;
                        for (int c = 0; c < nc; ++c) {
                            const std::size_t r = static_cast<std::size_t>(C.first_function + c);
                            for (int d = 0; d < nd; ++d) {
                                const double x = deg * *v++;
                                const std::size_t s = static_cast<std::size_t>(D.first_function + d);
                                gpq += Dm(r, s) * x;
                                G(r, s) += dpq * x;
                                const double xk = 0.25 * x;
                                G(p, r) -= Dm(q, s) * xk;
                                G(q, s) -= Dm(p, r) * xk;
                                G(p, s) -= Dm(q, r) * xk;
                                G(q, r) -= Dm(p, s) * xk;
                            }
                        }
                        G(p, q) += gpq;
                    }
                }
            }
        }
    });
    std::vector<Matrix> out(nx, Matrix(N, N));
    for (std::size_t x = 0; x < nx; ++x) {
        Matrix G(N, N);
        for (const auto& p : parts)
            for (std::size_t r = 0; r < N; ++r)
                for (std::size_t c = 0; c < N; ++c) G(r, c) += p[x](r, c);
        for (std::size_t r = 0; r < N; ++r)
            for (std::size_t c = 0; c < N; ++c) out[x](r, c) = 0.25 * (G(r, c) + G(c, r));
    }
    return out;
}

// ---------------------------------------------------------------------------------------
// geometry helpers
// ---------------------------------------------------------------------------------------

BasisSet basis_at_geometry(const BasisSet& basis, const Molecule& mol) {
    BasisSet b = basis;
    for (Shell& s : b.shells) s.center = mol.atoms.at(static_cast<std::size_t>(s.atom)).position;
    return b;
}

std::vector<double> nuclear_repulsion_gradient(const Molecule& mol) {
    const std::size_t n = mol.atoms.size();
    std::vector<double> g(3 * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            double d[3], r2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                d[k] = mol.atoms[i].position[k] - mol.atoms[j].position[k];
                r2 += d[k] * d[k];
            }
            const double r = std::sqrt(r2), f = mol.atoms[i].Z * mol.atoms[j].Z / (r2 * r);
            for (int k = 0; k < 3; ++k) g[3 * i + k] -= f * d[k];
        }
    return g;
}

Matrix nuclear_repulsion_hessian(const Molecule& mol) {
    const std::size_t n = mol.atoms.size();
    Matrix H(3 * n, 3 * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            double d[3], r2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                d[k] = mol.atoms[i].position[k] - mol.atoms[j].position[k];
                r2 += d[k] * d[k];
            }
            const double r = std::sqrt(r2), zz = mol.atoms[i].Z * mol.atoms[j].Z;
            // d^2/dRi_k dRi_l of Z Z / r = Z Z (3 d_k d_l / r^5 - delta_kl / r^3); the (i, j) block is minus that.
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l) {
                    const double v = zz * (3.0 * d[k] * d[l] / (r2 * r2 * r) - (k == l ? 1.0 / (r2 * r) : 0.0));
                    H(3 * i + k, 3 * i + l) += v;
                    H(3 * i + k, 3 * j + l) -= v;
                }
        }
    return H;
}

}  // namespace macqc
