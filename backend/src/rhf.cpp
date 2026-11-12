#include "macqc/rhf.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// Runs body(thread_index) on `threads` threads and rethrows the first exception.
template <class Body>
void parallel(int threads, Body body) {
    if (threads <= 1) {
        body(0);
        return;
    }
    std::vector<std::thread> pool;
    std::exception_ptr error;
    std::mutex mutex;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t] {
            try {
                body(t);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex);
                if (!error) error = std::current_exception();
            }
        });
    for (auto& th : pool) th.join();
    if (error) std::rethrow_exception(error);
}

// C = A B (i-k-j order).
Matrix mm(const Matrix& a, const Matrix& b) {
    const std::size_t n = a.rows(), m = b.cols(), kk = a.cols();
    Matrix c(n, m);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < kk; ++k) {
            const double aik = a(i, k);
            if (aik == 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) c(i, j) += aik * b(k, j);
        }
    return c;
}

// C = A^T B.
Matrix mtm(const Matrix& a, const Matrix& b) {
    const std::size_t n = a.cols(), m = b.cols(), kk = a.rows();
    Matrix c(n, m);
    for (std::size_t k = 0; k < kk; ++k)
        for (std::size_t i = 0; i < n; ++i) {
            const double aki = a(k, i);
            if (aki == 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) c(i, j) += aki * b(k, j);
        }
    return c;
}

double dot(const Matrix& a, const Matrix& b) {  // sum_ij a_ij b_ij
    double s = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) s += a(i, j) * b(i, j);
    return s;
}

Matrix add(const Matrix& a, const Matrix& b, double fb = 1.0) {
    Matrix c(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) c(i, j) = a(i, j) + fb * b(i, j);
    return c;
}

// Solves A x = b (Gaussian elimination, partial pivoting); false if singular.
bool solve_linear(std::vector<double> a, std::vector<double> b, std::size_t n, std::vector<double>& x) {
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (std::fabs(a[r * n + col]) > std::fabs(a[piv * n + col])) piv = r;
        if (std::fabs(a[piv * n + col]) < 1.0e-14) return false;
        if (piv != col) {
            for (std::size_t c = 0; c < n; ++c) std::swap(a[piv * n + c], a[col * n + c]);
            std::swap(b[piv], b[col]);
        }
        for (std::size_t r = col + 1; r < n; ++r) {
            const double f = a[r * n + col] / a[col * n + col];
            if (f == 0.0) continue;
            for (std::size_t c = col; c < n; ++c) a[r * n + c] -= f * a[col * n + c];
            b[r] -= f * b[col];
        }
    }
    x.assign(n, 0.0);
    for (std::size_t r = n; r-- > 0;) {
        double s = b[r];
        for (std::size_t c = r + 1; c < n; ++c) s -= a[r * n + c] * x[c];
        x[r] = s / a[r * n + r];
    }
    return true;
}

// Density P = 2 C_occ C_occ^T.
Matrix density(const Matrix& C, int nocc) {
    const std::size_t n = C.rows();
    Matrix P(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double s = 0.0;
            for (int k = 0; k < nocc; ++k) s += C(i, k) * C(j, k);
            P(i, j) = 2.0 * s;
            P(j, i) = 2.0 * s;
        }
    return P;
}

// Sums the per-thread contributions, symmetrises and halves (density = total P).
Matrix reduce_g(const std::vector<Matrix>& parts, std::size_t n) {
    Matrix G(n, n);
    for (const Matrix& g : parts)
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) G(i, j) += g(i, j);
    Matrix out(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) out(i, j) = 0.25 * (G(i, j) + G(j, i));
    return out;
}

constexpr double kStoreScreen = 1.0e-14;   // Schwarz bound below which a quartet is zero
constexpr double kDirectScreen = 1.0e-13;  // Schwarz x density bound in direct builds

}  // namespace

std::size_t unique_eri_count(std::size_t n) {
    const std::size_t m = n * (n + 1) / 2;
    return m * (m + 1) / 2;
}

FockBuilder::FockBuilder(const BasisSet& basis, int threads, double memory_limit_mb)
    : basis_(basis), threads_(resolve_threads(threads)) {
    const auto t0 = Clock::now();
    pairs_ = make_shell_pairs(basis);
    const std::size_t npairs = pairs_.pairs.size();

    // Schwarz bounds, in parallel.
    {
        std::atomic<std::size_t> next{0};
        parallel(threads_, [&](int) {
            EriEngine engine(basis_);
            for (std::size_t p; (p = next++) < npairs;) {
                const ShellPair& sp = pairs_.pairs[p];
                double mx = 0.0;
                if (!sp.prims.empty()) {
                    const double* v = engine.compute(sp, sp);
                    const std::size_t nab =
                        static_cast<std::size_t>(basis_.shells[sp.a].size()) * basis_.shells[sp.b].size();
                    for (std::size_t ij = 0; ij < nab; ++ij) mx = std::max(mx, std::fabs(v[ij * nab + ij]));
                }
                pairs_.schwarz[p] = std::sqrt(mx);
            }
        });
        pairs_.max_schwarz = 0.0;
        for (double q : pairs_.schwarz) pairs_.max_schwarz = std::max(pairs_.max_schwarz, q);
    }

    const std::size_t n = static_cast<std::size_t>(basis.n_functions);
    unique_ = unique_eri_count(n);
    const double mb = static_cast<double>(unique_) * 8.0 / (1024.0 * 1024.0);
    stored_ = memory_limit_mb > 0.0 && mb <= memory_limit_mb;

    std::atomic<std::size_t> quartets{0};
    if (stored_) {
        eri_.assign(unique_, 0.0);
        std::atomic<std::size_t> next{0};
        parallel(threads_, [&](int) {
            EriEngine engine(basis_);
            std::size_t local = 0;
            for (std::size_t k; (k = next++) < npairs;) {
                const std::size_t P = npairs - 1 - k;  // heaviest rows first
                const ShellPair& bra = pairs_.pairs[P];
                if (bra.prims.empty()) continue;
                const Shell& A = basis_.shells[bra.a];
                const Shell& B = basis_.shells[bra.b];
                for (std::size_t Q = 0; Q <= P; ++Q) {
                    const ShellPair& ket = pairs_.pairs[Q];
                    if (ket.prims.empty() || pairs_.schwarz[P] * pairs_.schwarz[Q] < kStoreScreen) continue;
                    ++local;
                    const Shell& C = basis_.shells[ket.a];
                    const Shell& D = basis_.shells[ket.b];
                    const double* v = engine.compute(bra, ket);
                    const int na = A.size(), nb = B.size(), nc = C.size(), nd = D.size();
                    for (int a = 0; a < na; ++a)
                        for (int b = 0; b < nb; ++b) {
                            const std::size_t ab = pair_index(A.first_function + a, B.first_function + b);
                            for (int c = 0; c < nc; ++c)
                                for (int d = 0; d < nd; ++d) {
                                    const std::size_t cd = pair_index(C.first_function + c, D.first_function + d);
                                    eri_[pair_index(ab, cd)] = *v++;
                                }
                        }
                }
            }
            quartets += local;
        });
    } else {
        std::size_t count = 0;
        for (std::size_t P = 0; P < npairs; ++P)
            for (std::size_t Q = 0; Q <= P; ++Q)
                if (pairs_.schwarz[P] * pairs_.schwarz[Q] >= kStoreScreen) ++count;
        quartets = count;
    }
    quartets_ = quartets;
    setup_seconds_ = seconds_since(t0);
}

Matrix FockBuilder::build(const Matrix& P) const { return stored_ ? build_stored(P) : build_direct(P); }

Matrix FockBuilder::build_stored(const Matrix& D) const {
    const std::size_t n = static_cast<std::size_t>(basis_.n_functions);
    std::vector<Matrix> parts(static_cast<std::size_t>(threads_), Matrix(n, n));
    std::atomic<std::size_t> next{0};
    parallel(threads_, [&](int t) {
        Matrix& G = parts[static_cast<std::size_t>(t)];
        for (std::size_t k; (k = next++) < n;) {
            const std::size_t p = n - 1 - k;
            for (std::size_t q = 0; q <= p; ++q) {
                const std::size_t pq = p * (p + 1) / 2 + q;
                const double* v = &eri_[pq * (pq + 1) / 2];
                const double dpq = D(p, q);
                double gpq = 0.0;
                for (std::size_t r = 0; r <= p; ++r) {
                    const std::size_t smax = (r == p) ? q : r;
                    for (std::size_t s = 0; s <= smax; ++s) {
                        double x = *v++;
                        if (x == 0.0) continue;
                        const std::size_t rs = r * (r + 1) / 2 + s;
                        if (p != q) x *= 2.0;
                        if (r != s) x *= 2.0;
                        if (pq != rs) x *= 2.0;
                        gpq += D(r, s) * x;
                        G(r, s) += dpq * x;
                        const double xk = 0.25 * x;
                        G(p, r) -= D(q, s) * xk;
                        G(q, s) -= D(p, r) * xk;
                        G(p, s) -= D(q, r) * xk;
                        G(q, r) -= D(p, s) * xk;
                    }
                }
                G(p, q) += gpq;
            }
        }
    });
    return reduce_g(parts, n);
}

Matrix FockBuilder::build_direct(const Matrix& D) const {
    const std::size_t n = static_cast<std::size_t>(basis_.n_functions);
    const std::size_t ns = basis_.shells.size();
    // Largest |D| per shell block, for screening.
    std::vector<double> dmax(ns * ns, 0.0);
    for (std::size_t a = 0; a < ns; ++a)
        for (std::size_t b = 0; b < ns; ++b) {
            const Shell& A = basis_.shells[a];
            const Shell& B = basis_.shells[b];
            double m = 0.0;
            for (int i = 0; i < A.size(); ++i)
                for (int j = 0; j < B.size(); ++j)
                    m = std::max(m, std::fabs(D(A.first_function + i, B.first_function + j)));
            dmax[a * ns + b] = m;
        }
    const std::size_t npairs = pairs_.pairs.size();
    std::vector<Matrix> parts(static_cast<std::size_t>(threads_), Matrix(n, n));
    std::atomic<std::size_t> next{0};
    parallel(threads_, [&](int t) {
        Matrix& G = parts[static_cast<std::size_t>(t)];
        EriEngine engine(basis_);
        for (std::size_t k; (k = next++) < npairs;) {
            const std::size_t P = npairs - 1 - k;
            const ShellPair& bra = pairs_.pairs[P];
            if (bra.prims.empty()) continue;
            const std::size_t sa = static_cast<std::size_t>(bra.a), sb = static_cast<std::size_t>(bra.b);
            const Shell& A = basis_.shells[sa];
            const Shell& B = basis_.shells[sb];
            for (std::size_t Q = 0; Q <= P; ++Q) {
                const ShellPair& ket = pairs_.pairs[Q];
                if (ket.prims.empty()) continue;
                const double bound = pairs_.schwarz[P] * pairs_.schwarz[Q];
                if (bound < kStoreScreen) continue;
                const std::size_t sc = static_cast<std::size_t>(ket.a), sd = static_cast<std::size_t>(ket.b);
                const double dm = std::max({dmax[sa * ns + sb], dmax[sc * ns + sd], dmax[sa * ns + sc],
                                            dmax[sb * ns + sd], dmax[sa * ns + sd], dmax[sb * ns + sc]});
                if (bound * dm < kDirectScreen) continue;
                const Shell& C = basis_.shells[sc];
                const Shell& Dd = basis_.shells[sd];
                const double* v = engine.compute(bra, ket);
                const double deg = (sa == sb ? 1.0 : 2.0) * (sc == sd ? 1.0 : 2.0) * (P == Q ? 1.0 : 2.0);
                const int na = A.size(), nb = B.size(), nc = C.size(), nd = Dd.size();
                for (int a = 0; a < na; ++a) {
                    const std::size_t p = static_cast<std::size_t>(A.first_function + a);
                    for (int b = 0; b < nb; ++b) {
                        const std::size_t q = static_cast<std::size_t>(B.first_function + b);
                        const double dpq = D(p, q);
                        double gpq = 0.0;
                        for (int c = 0; c < nc; ++c) {
                            const std::size_t r = static_cast<std::size_t>(C.first_function + c);
                            for (int d = 0; d < nd; ++d) {
                                const double x = deg * *v++;
                                const std::size_t s = static_cast<std::size_t>(Dd.first_function + d);
                                gpq += D(r, s) * x;
                                G(r, s) += dpq * x;
                                const double xk = 0.25 * x;
                                G(p, r) -= D(q, s) * xk;
                                G(q, s) -= D(p, r) * xk;
                                G(p, s) -= D(q, r) * xk;
                                G(q, r) -= D(p, s) * xk;
                            }
                        }
                        G(p, q) += gpq;
                    }
                }
            }
        }
    });
    return reduce_g(parts, n);
}

RHFResult run_rhf(const Molecule& mol, const BasisSet& basis, const RHFOptions& opt, const RHFProgress& progress) {
    const int nel = mol.electrons();
    if (nel <= 0 || nel % 2 != 0)
        throw std::runtime_error("restricted Hartree-Fock needs an even, positive number of electrons");
    const int nocc = nel / 2;
    const std::size_t n = static_cast<std::size_t>(basis.n_functions);
    if (n < static_cast<std::size_t>(nocc)) throw std::runtime_error("the basis set has fewer functions than occupied orbitals");

    RHFResult r;
    r.n_basis = basis.n_functions;
    r.n_occupied = nocc;
    r.nuclear_repulsion = mol.nuclear_repulsion();
    const auto t0 = Clock::now();

    const Matrix S = overlap_matrix(basis);
    const Matrix T = kinetic_matrix(basis);
    const Matrix V = nuclear_attraction_matrix(basis, mol);
    const Matrix H = add(T, V);
    r.S = S;
    r.H = H;

    // Canonical orthogonalisation X = U s^-1/2 (eigenvalues below lindep_tol removed).
    const EigenResult se = symmetric_eigen(S);
    r.smallest_overlap_eigenvalue = se.values.empty() ? 0.0 : se.values.front();
    std::vector<std::size_t> keep;
    for (std::size_t k = 0; k < n; ++k)
        if (se.values[k] >= opt.lindep_tol) keep.push_back(k);
    const std::size_t nmo = keep.size();
    r.removed_functions = static_cast<int>(n - nmo);
    r.n_mo = static_cast<int>(nmo);
    if (nmo < static_cast<std::size_t>(nocc)) throw std::runtime_error("too many linearly dependent basis functions");
    Matrix X(n, nmo);
    for (std::size_t k = 0; k < nmo; ++k) {
        const double f = 1.0 / std::sqrt(se.values[keep[k]]);
        for (std::size_t i = 0; i < n; ++i) X(i, k) = se.vectors(i, keep[k]) * f;
    }

    FockBuilder fock(basis, opt.threads, opt.memory_limit_mb);
    r.eri_stored = fock.stored();
    r.unique_integrals = fock.unique_integrals();
    r.threads = fock.threads();
    r.integral_seconds = seconds_since(t0);
    const auto t1 = Clock::now();

    auto diagonalize = [&](const Matrix& F, Matrix& C, std::vector<double>& eps) {
        const Matrix Fp = mtm(X, mm(F, X));
        const EigenResult e = symmetric_eigen(Fp);
        C = mm(X, e.vectors);
        eps = e.values;
    };

    // Core-Hamiltonian guess.
    Matrix C;
    std::vector<double> eps;
    diagonalize(H, C, eps);
    Matrix P = density(C, nocc);

    std::vector<Matrix> diis_f, diis_e;
    Matrix G(n, n), P_built(n, n);
    double e_prev = 0.0;
    Matrix F;
    int since_full = 0;
    for (int it = 1; it <= opt.max_iter; ++it) {
        // Fock matrix; integral-direct builds are incremental in the density change.
        if (fock.stored() || since_full >= 8 || it == 1) {
            G = fock.build(P);
            since_full = 0;
        } else {
            G = add(G, fock.build(add(P, P_built, -1.0)));
            ++since_full;
        }
        P_built = P;
        F = add(H, G);
        const double e_el = 0.5 * (dot(P, H) + dot(P, F));
        const double e_tot = e_el + r.nuclear_repulsion;

        // Orbital gradient X^T (F P S - S P F) X.
        const Matrix FPS = mm(F, mm(P, S));
        Matrix err = mtm(X, mm(add(FPS, FPS.transposed(), -1.0), X));
        double gmax = 0.0;
        for (std::size_t i = 0; i < nmo; ++i)
            for (std::size_t j = 0; j < nmo; ++j) gmax = std::max(gmax, std::fabs(err(i, j)));

        RHFIteration rec;
        rec.iteration = it;
        rec.energy = e_tot;
        rec.delta_e = it == 1 ? 0.0 : e_tot - e_prev;
        rec.gradient = gmax;
        r.electronic_energy = e_el;
        r.total_energy = e_tot;
        r.iterations = it;

        const bool done = it > 1 && std::fabs(e_tot - e_prev) < opt.energy_tol && gmax < opt.gradient_tol;
        e_prev = e_tot;
        if (done) {
            r.converged = true;
            r.history.push_back(rec);
            if (progress) progress(rec);
            break;
        }

        // DIIS extrapolation.
        Matrix Fx = F;
        diis_f.push_back(F);
        diis_e.push_back(err);
        if (static_cast<int>(diis_f.size()) > opt.diis_size) {
            diis_f.erase(diis_f.begin());
            diis_e.erase(diis_e.begin());
        }
        while (diis_f.size() >= 2) {
            const std::size_t m = diis_f.size(), dim = m + 1;
            std::vector<double> B(dim * dim, 0.0), rhs(dim, 0.0), c;
            double scale = 0.0;
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j <= i; ++j) {
                    const double v = dot(diis_e[i], diis_e[j]);
                    B[i * dim + j] = B[j * dim + i] = v;
                    if (i == j) scale = std::max(scale, v);
                }
            if (scale > 0.0)
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t j = 0; j < m; ++j) B[i * dim + j] /= scale;
            for (std::size_t i = 0; i < m; ++i) B[i * dim + m] = B[m * dim + i] = -1.0;
            rhs[m] = -1.0;
            if (solve_linear(B, rhs, dim, c)) {
                Fx = Matrix(n, n);
                for (std::size_t i = 0; i < m; ++i) Fx = add(Fx, diis_f[i], c[i]);
                rec.diis = true;
                break;
            }
            diis_f.erase(diis_f.begin());
            diis_e.erase(diis_e.begin());
        }
        r.history.push_back(rec);
        if (progress) progress(rec);

        diagonalize(Fx, C, eps);
        P = density(C, nocc);
    }

    // Final orbitals from the (non-extrapolated) Fock matrix of the last density.
    if (r.converged) diagonalize(F, C, eps);
    r.C = C;
    r.orbital_energies = eps;
    r.P = P_built;
    r.kinetic_energy = dot(P_built, T);
    r.nuclear_attraction_energy = dot(P_built, V);
    r.two_electron_energy = 0.5 * dot(P_built, G);
    const double potential = r.total_energy - r.kinetic_energy;
    r.virial_ratio = r.kinetic_energy != 0.0 ? -potential / r.kinetic_energy : 0.0;
    r.scf_seconds = seconds_since(t1);
    return r;
}

}  // namespace macqc
