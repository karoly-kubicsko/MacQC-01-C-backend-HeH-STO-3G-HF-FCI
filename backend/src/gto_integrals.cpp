#include "macqc/gto_integrals.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <thread>

namespace macqc {

namespace {

constexpr double kPiExact = 3.14159265358979323846;
constexpr int kMaxLL = 2 * kMaxAngularMomentum;      // highest l of a shell pair, + headroom
constexpr int kMaxM = 4 * kMaxAngularMomentum + 2;   // highest Boys order supported

// ---------------------------------------------------------------- Boys function

constexpr double kBoysStep = 0.05;
constexpr double kBoysTMax = 45.0;
constexpr int kBoysGrid = 901;               // T = 0, 0.05, ..., 45
constexpr int kBoysCols = kMaxM + 8;          // orders 0 .. kMaxM + 7 (7 Taylor terms)

struct BoysTable {
    std::vector<double> v;
    BoysTable() : v(static_cast<std::size_t>(kBoysGrid) * kBoysCols) {
        for (int k = 0; k < kBoysGrid; ++k) {
            const double T = k * kBoysStep;
            double* f = &v[static_cast<std::size_t>(k) * kBoysCols];
            const int top = kBoysCols - 1;
            // F_top(T) = exp(-T) sum_i (2T)^i / ((2top+1)(2top+3)...(2top+2i+1)): positive terms.
            double term = 1.0 / (2 * top + 1), sum = term;
            for (int i = 1; i < 4000; ++i) {
                term *= 2.0 * T / (2 * top + 2 * i + 1);
                sum += term;
                if (term < 1.0e-18 * sum) break;
            }
            const double e = std::exp(-T);
            f[top] = e * sum;
            for (int m = top - 1; m >= 0; --m) f[m] = (2.0 * T * f[m + 1] + e) / (2 * m + 1);
        }
    }
};

const BoysTable& boys_table() {
    static const BoysTable table;
    return table;
}

// ---------------------------------------------------------------- Cartesian tables

inline int ncart(int l) { return (l + 1) * (l + 2) / 2; }
inline int cart_index(int /*x*/, int y, int z) { return ((y + z) * (y + z + 1)) / 2 + z; }

struct CartTables {
    std::vector<std::array<int, 3>> comp[kMaxLL + 2];
    std::vector<std::array<int, 3>> down[kMaxLL + 2];  // index in l-1 of c - 1_i, or -1
    std::vector<std::array<int, 3>> up[kMaxLL + 2];    // index in l+1 of c + 1_i
    std::vector<int> dir[kMaxLL + 2];                  // recursion direction (l >= 1)
    // For each l and direction i: the components c with c_i > 0, as (c, index of c - 1_i, c_i).
    std::vector<std::array<int, 3>> lowered[kMaxLL + 2][3];
    CartTables() {
        for (int l = 0; l <= kMaxLL + 1; ++l) {
            for (int x = l; x >= 0; --x)
                for (int y = l - x; y >= 0; --y) comp[l].push_back({x, y, l - x - y});
            for (const auto& c : comp[l]) {
                std::array<int, 3> d{-1, -1, -1}, u{};
                for (int i = 0; i < 3; ++i) {
                    auto cc = c;
                    if (cc[i] > 0) {
                        --cc[i];
                        d[i] = cart_index(cc[0], cc[1], cc[2]);
                    }
                    auto cu = c;
                    ++cu[i];
                    u[i] = cart_index(cu[0], cu[1], cu[2]);
                }
                down[l].push_back(d);
                up[l].push_back(u);
                // Direction with the smallest positive exponent: the parent then has the
                // fewest terms in the recurrence.
                int best = -1;
                for (int i = 0; i < 3; ++i)
                    if (c[i] > 0 && (best < 0 || c[i] < c[best])) best = i;
                dir[l].push_back(best);
            }
            for (int i = 0; i < 3; ++i)
                for (std::size_t c = 0; c < comp[l].size(); ++c)
                    if (comp[l][c][i] > 0)
                        lowered[l][i].push_back({static_cast<int>(c), down[l][c][i], comp[l][c][i]});
        }
    }
};

const CartTables& cart_tables() {
    static const CartTables t;
    return t;
}

// ---------------------------------------------------------------- solid harmonics

double factorial(int n) {
    double r = 1.0;
    for (int k = 2; k <= n; ++k) r *= k;
    return r;
}

double binomial(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    double r = 1.0;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

// Helgaker, Jorgensen, Olsen, Molecular Electronic-Structure Theory, Eqs. 6.4.47-6.4.50.
std::vector<double> make_solid_harmonics(int l) {
    const int nc = ncart(l);
    std::vector<double> M(static_cast<std::size_t>(2 * l + 1) * nc, 0.0);
    for (int m = -l; m <= l; ++m) {
        const int am = std::abs(m);
        const double norm = 1.0 / (std::pow(2.0, am) * factorial(l)) *
                            std::sqrt(2.0 * factorial(l + am) * factorial(l - am) / (m == 0 ? 2.0 : 1.0));
        const int w0 = m < 0 ? 1 : 0;  // w = 2v
        for (int t = 0; t <= (l - am) / 2; ++t)
            for (int u = 0; u <= t; ++u)
                for (int w = w0; w <= am; w += 2) {
                    const int sign_exp = t + (w - w0) / 2;
                    const double c = ((sign_exp % 2) ? -1.0 : 1.0) * std::pow(0.25, t) * binomial(l, t) *
                                     binomial(l - t, am + t) * binomial(t, u) * binomial(am, w);
                    const int x = 2 * t + am - 2 * u - w, y = 2 * u + w, z = l - 2 * t - am;
                    if (x < 0) continue;
                    M[static_cast<std::size_t>(m + l) * nc + cart_index(x, y, z)] += norm * c;
                }
    }
    return M;
}

struct SolidTables {
    std::vector<double> dense[kMaxAngularMomentum + 1];
    // Sparse rows: for each m, list of (cartesian index, coefficient).
    std::vector<std::vector<std::pair<int, double>>> rows[kMaxAngularMomentum + 1];
    SolidTables() {
        for (int l = 0; l <= kMaxAngularMomentum; ++l) {
            dense[l] = make_solid_harmonics(l);
            const int nc = ncart(l);
            rows[l].resize(2 * l + 1);
            for (int m = 0; m < 2 * l + 1; ++m)
                for (int c = 0; c < nc; ++c) {
                    const double v = dense[l][static_cast<std::size_t>(m) * nc + c];
                    if (v != 0.0) rows[l][m].emplace_back(c, v);
                }
        }
    }
};

const SolidTables& solid_tables() {
    static const SolidTables t;
    return t;
}

struct TransformTables {
    ShellTransform sph[kMaxAngularMomentum + 1];
    TransformTables() {
        const SolidTables& st = solid_tables();
        for (int l = 0; l <= kMaxAngularMomentum; ++l) {
            sph[l].l = l;
            sph[l].rows = st.rows[l];
        }
    }
};

const TransformTables& transform_tables() {
    static const TransformTables t;
    return t;
}

void transform_rows(const std::vector<std::vector<std::pair<int, double>>>& rows, int nc, int pre, int post,
                    const double* in, double* out) {
    const int ns = static_cast<int>(rows.size());
    for (int p = 0; p < pre; ++p) {
        const double* src = in + static_cast<std::size_t>(p) * nc * post;
        double* dst = out + static_cast<std::size_t>(p) * ns * post;
        for (int m = 0; m < ns; ++m) {
            double* d = dst + static_cast<std::size_t>(m) * post;
            const auto& row = rows[m];
            if (row.empty()) {
                for (int q = 0; q < post; ++q) d[q] = 0.0;
                continue;
            }
            {
                const double c0 = row[0].second;
                const double* s0 = src + static_cast<std::size_t>(row[0].first) * post;
                for (int q = 0; q < post; ++q) d[q] = c0 * s0[q];
            }
            for (std::size_t r = 1; r < row.size(); ++r) {
                const double c = row[r].second;
                const double* s = src + static_cast<std::size_t>(row[r].first) * post;
                for (int q = 0; q < post; ++q) d[q] += c * s[q];
            }
        }
    }
}

// in [pre][ncart(l)][post] -> out [pre][2l+1][post]
void to_spherical(int l, int pre, int post, const double* in, double* out) {
    transform_rows(solid_tables().rows[l], ncart(l), pre, post, in, out);
}

// The transform of a pair member: t, or the solid harmonics of l.
inline const ShellTransform& pair_transform(const ShellTransform* t, int l) {
    return t ? *t : transform_tables().sph[l];
}

// Horizontal recurrence (a|b+1_i) = (a+1_i|b) + AB_i (a|b) on blocks with nk trailing values.
// in[e - la] points to the block [c][k] of e = la .. la + lb; out receives [a][b][k].
void hrr(int la, int lb, const double* AB, int nk, const double* const* in, std::vector<double>& work,
         double* out) {
    const CartTables& T = cart_tables();
    if (lb == 0) {
        std::copy(in[0], in[0] + static_cast<std::size_t>(ncart(la)) * nk, out);
        return;
    }
    std::size_t total = 0;
    for (int j = 1; j < lb; ++j)
        for (int e = la; e <= la + lb - j; ++e) total += static_cast<std::size_t>(ncart(e)) * ncart(j) * nk;
    if (work.size() < total) work.resize(total);
    const double* prev[kMaxLL + 2];
    double* cur[kMaxLL + 2];
    for (int e = la; e <= la + lb; ++e) prev[e - la] = in[e - la];
    std::size_t pos = 0;
    for (int j = 1; j <= lb; ++j) {
        const int nj = ncart(j), nj1 = ncart(j - 1);
        for (int e = la; e <= la + lb - j; ++e) {
            if (j == lb) {
                cur[e - la] = out;
            } else {
                cur[e - la] = work.data() + pos;
                pos += static_cast<std::size_t>(ncart(e)) * nj * nk;
            }
            const int ne = ncart(e);
            const double* p_up = prev[e + 1 - la];
            const double* p_same = prev[e - la];
            double* dst = cur[e - la];
            for (int c = 0; c < ne; ++c)
                for (int d = 0; d < nj; ++d) {
                    const int i = T.dir[j][d];
                    const int pd = T.down[j][d][i];
                    const int cu = T.up[e][c][i];
                    const double ab = AB[i];
                    const double* s1 = p_up + (static_cast<std::size_t>(cu) * nj1 + pd) * nk;
                    const double* s2 = p_same + (static_cast<std::size_t>(c) * nj1 + pd) * nk;
                    double* t = dst + (static_cast<std::size_t>(c) * nj + d) * nk;
                    for (int k = 0; k < nk; ++k) t[k] = s1[k] + ab * s2[k];
                }
        }
        for (int e = la; e <= la + lb - j; ++e) prev[e - la] = cur[e - la];
    }
}

// ---------------------------------------------------------------- one-electron integrals

// 1D overlap S(i, j) (without the exp/pi prefactor), i <= imax, j <= jmax.
void overlap_1d(double p, double PA, double PB, int imax, int jmax, double* S) {
    const int nj = jmax + 1;
    const double o2p = 0.5 / p;
    S[0] = 1.0;
    for (int i = 0; i < imax; ++i) S[(i + 1) * nj] = PA * S[i * nj] + (i > 0 ? i * o2p * S[(i - 1) * nj] : 0.0);
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i <= imax; ++i) {
            double v = PB * S[i * nj + j];
            if (i > 0) v += i * o2p * S[(i - 1) * nj + j];
            if (j > 0) v += j * o2p * S[i * nj + j - 1];
            S[i * nj + j + 1] = v;
        }
}

enum class OneElectron { Overlap, Kinetic };

// Cartesian block [ca][cb] of S or T for the shell pair (a, b).
void overlap_kinetic_block(const Shell& A, const Shell& B, OneElectron kind, std::vector<double>& cart) {
    const CartTables& T = cart_tables();
    const int la = A.l, lb = B.l, na = ncart(la), nb = ncart(lb);
    cart.assign(static_cast<std::size_t>(na) * nb, 0.0);
    const double AB2 = (A.center[0] - B.center[0]) * (A.center[0] - B.center[0]) +
                       (A.center[1] - B.center[1]) * (A.center[1] - B.center[1]) +
                       (A.center[2] - B.center[2]) * (A.center[2] - B.center[2]);
    const int jmax = lb + 2, nj = jmax + 1;
    std::vector<double> S[3], K[3];
    for (int d = 0; d < 3; ++d) {
        S[d].resize(static_cast<std::size_t>(la + 1) * nj);
        K[d].resize(static_cast<std::size_t>(la + 1) * (lb + 1));
    }
    for (std::size_t pa = 0; pa < A.exponents.size(); ++pa)
        for (std::size_t pb = 0; pb < B.exponents.size(); ++pb) {
            const double a = A.exponents[pa], b = B.exponents[pb], p = a + b;
            const double pref = A.coefficients[pa] * B.coefficients[pb] * std::pow(kPiExact / p, 1.5) *
                                std::exp(-a * b / p * AB2);
            if (pref == 0.0) continue;
            for (int d = 0; d < 3; ++d) {
                const double P = (a * A.center[d] + b * B.center[d]) / p;
                overlap_1d(p, P - A.center[d], P - B.center[d], la, jmax, S[d].data());
                if (kind == OneElectron::Kinetic)
                    for (int i = 0; i <= la; ++i)
                        for (int j = 0; j <= lb; ++j) {
                            double v = -2.0 * b * b * S[d][i * nj + j + 2] + b * (2 * j + 1) * S[d][i * nj + j];
                            if (j >= 2) v -= 0.5 * j * (j - 1) * S[d][i * nj + j - 2];
                            K[d][i * (lb + 1) + j] = v;
                        }
            }
            for (int ca = 0; ca < na; ++ca) {
                const auto& ea = T.comp[la][ca];
                for (int cb = 0; cb < nb; ++cb) {
                    const auto& eb = T.comp[lb][cb];
                    const double sx = S[0][ea[0] * nj + eb[0]], sy = S[1][ea[1] * nj + eb[1]],
                                 sz = S[2][ea[2] * nj + eb[2]];
                    double v;
                    if (kind == OneElectron::Overlap) {
                        v = sx * sy * sz;
                    } else {
                        const double kx = K[0][ea[0] * (lb + 1) + eb[0]], ky = K[1][ea[1] * (lb + 1) + eb[1]],
                                     kz = K[2][ea[2] * (lb + 1) + eb[2]];
                        v = kx * sy * sz + sx * ky * sz + sx * sy * kz;
                    }
                    cart[static_cast<std::size_t>(ca) * nb + cb] += pref * v;
                }
            }
        }
}

// Cartesian block [ca][cb] of the nuclear attraction for the shell pair (a, b).
void nuclear_block(const Shell& A, const Shell& B, const std::vector<Atom>& nuclei, std::vector<double>& cart,
                   std::vector<double>& work) {
    const CartTables& T = cart_tables();
    const int la = A.l, lb = B.l, L = la + lb;
    double AB[3];
    for (int d = 0; d < 3; ++d) AB[d] = A.center[d] - B.center[d];
    const double AB2 = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
    // VRR storage: level e has ncart(e) x (L - e + 1) values.
    std::size_t off[kMaxLL + 2], pos = 0;
    for (int e = 0; e <= L; ++e) {
        off[e] = pos;
        pos += static_cast<std::size_t>(ncart(e)) * (L - e + 1);
    }
    std::vector<double> v(pos);
    std::size_t aoff[kMaxLL + 2], apos = 0;
    for (int e = la; e <= L; ++e) {
        aoff[e - la] = apos;
        apos += ncart(e);
    }
    std::vector<double> acc(apos, 0.0);
    double F[kMaxM + 2];
    for (std::size_t pa = 0; pa < A.exponents.size(); ++pa)
        for (std::size_t pb = 0; pb < B.exponents.size(); ++pb) {
            const double a = A.exponents[pa], b = B.exponents[pb], p = a + b;
            const double K = A.coefficients[pa] * B.coefficients[pb] * std::exp(-a * b / p * AB2);
            if (K == 0.0) continue;
            double P[3], PA[3];
            for (int d = 0; d < 3; ++d) {
                P[d] = (a * A.center[d] + b * B.center[d]) / p;
                PA[d] = P[d] - A.center[d];
            }
            const double o2p = 0.5 / p;
            for (const Atom& atom : nuclei) {
                double PC[3];
                for (int d = 0; d < 3; ++d) PC[d] = P[d] - atom.position[d];
                const double Tb = p * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
                boys_function(L, Tb, F);
                const double pref = -atom.Z * 2.0 * kPiExact / p * K;
                for (int m = 0; m <= L; ++m) v[m] = pref * F[m];
                for (int e = 1; e <= L; ++e) {
                    const int nm = L - e + 1, nm1 = nm + 1, nm2 = nm + 2;
                    for (int c = 0; c < ncart(e); ++c) {
                        const int i = T.dir[e][c];
                        const int pc = T.down[e][c][i];
                        const int ai = T.comp[e][c][i] - 1;
                        const double* x1 = &v[off[e - 1] + static_cast<std::size_t>(pc) * nm1];
                        double* t = &v[off[e] + static_cast<std::size_t>(c) * nm];
                        for (int m = 0; m < nm; ++m) t[m] = PA[i] * x1[m] - PC[i] * x1[m + 1];
                        if (ai > 0) {
                            const int pp = T.down[e - 1][pc][i];
                            const double* x2 = &v[off[e - 2] + static_cast<std::size_t>(pp) * nm2];
                            const double f = ai * o2p;
                            for (int m = 0; m < nm; ++m) t[m] += f * (x2[m] - x2[m + 1]);
                        }
                    }
                }
                for (int e = la; e <= L; ++e) {
                    const int nm = L - e + 1;
                    for (int c = 0; c < ncart(e); ++c) acc[aoff[e - la] + c] += v[off[e] + static_cast<std::size_t>(c) * nm];
                }
            }
        }
    const double* in[kMaxLL + 2];
    for (int e = la; e <= L; ++e) in[e - la] = acc.data() + aoff[e - la];
    cart.assign(static_cast<std::size_t>(ncart(la)) * ncart(lb), 0.0);
    hrr(la, lb, AB, 1, in, work, cart.data());
}

enum class MatrixKind { Overlap, Kinetic, Nuclear };

Matrix one_electron_matrix(const BasisSet& basis, const Molecule* mol, MatrixKind kind) {
    const int n = basis.n_functions;
    Matrix M(n, n);
    std::vector<double> cart, tmp, sph, work;
    for (std::size_t i = 0; i < basis.shells.size(); ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            const Shell& A = basis.shells[i];
            const Shell& B = basis.shells[j];
            if (kind == MatrixKind::Nuclear)
                nuclear_block(A, B, mol->atoms, cart, work);
            else
                overlap_kinetic_block(A, B, kind == MatrixKind::Overlap ? OneElectron::Overlap : OneElectron::Kinetic,
                                      cart);
            const int na = ncart(A.l), sa = A.size(), sb = B.size();
            tmp.resize(static_cast<std::size_t>(na) * sb);
            sph.resize(static_cast<std::size_t>(sa) * sb);
            to_spherical(B.l, na, 1, cart.data(), tmp.data());
            to_spherical(A.l, 1, sb, tmp.data(), sph.data());
            for (int a = 0; a < sa; ++a)
                for (int b = 0; b < sb; ++b) {
                    const double v = sph[static_cast<std::size_t>(a) * sb + b];
                    M(A.first_function + a, B.first_function + b) = v;
                    M(B.first_function + b, A.first_function + a) = v;
                }
        }
    return M;
}

}  // namespace

// ---------------------------------------------------------------- public functions

void boys_function(int mmax, double T, double* F) {
    if (mmax < 0 || mmax > kMaxM) throw std::invalid_argument("boys_function: order out of range");
    if (T < kBoysTMax) {
        const BoysTable& tab = boys_table();
        const int k = static_cast<int>(T / kBoysStep + 0.5);
        const double d = k * kBoysStep - T;  // F_m(T) = sum_j F_{m+j}(T_k) d^j / j!
        const double* f = &tab.v[static_cast<std::size_t>(k) * kBoysCols + mmax];
        double r = f[6];
        r = f[5] + d * r / 6.0;
        r = f[4] + d * r / 5.0;
        r = f[3] + d * r / 4.0;
        r = f[2] + d * r / 3.0;
        r = f[1] + d * r / 2.0;
        r = f[0] + d * r;
        F[mmax] = r;
        if (mmax > 0) {
            const double e = std::exp(-T);
            for (int m = mmax - 1; m >= 0; --m) F[m] = (2.0 * T * F[m + 1] + e) / (2 * m + 1);
        }
    } else {
        F[0] = 0.5 * std::sqrt(kPiExact / T);
        if (mmax > 0) {
            const double e = std::exp(-T), o2t = 0.5 / T;
            for (int m = 0; m < mmax; ++m) F[m + 1] = ((2 * m + 1) * F[m] - e) * o2t;
        }
    }
}

const std::vector<double>& solid_harmonic_matrix(int l) {
    if (l < 0 || l > kMaxAngularMomentum) throw std::invalid_argument("solid_harmonic_matrix: l out of range");
    return solid_tables().dense[l];
}

Matrix overlap_matrix(const BasisSet& basis) { return one_electron_matrix(basis, nullptr, MatrixKind::Overlap); }

Matrix kinetic_matrix(const BasisSet& basis) { return one_electron_matrix(basis, nullptr, MatrixKind::Kinetic); }

Matrix nuclear_attraction_matrix(const BasisSet& basis, const Molecule& mol) {
    return one_electron_matrix(basis, &mol, MatrixKind::Nuclear);
}

ShellPairList make_shell_pairs(const BasisSet& basis) {
    ShellPairList list;
    const int ns = static_cast<int>(basis.shells.size());
    list.pairs.resize(static_cast<std::size_t>(ns) * (ns + 1) / 2);
    list.schwarz.assign(list.pairs.size(), 0.0);
    for (int i = 0; i < ns; ++i)
        for (int j = 0; j <= i; ++j) {
            ShellPair& sp = list.pairs[ShellPairList::index(i, j)];
            int a = i, b = j;
            if (basis.shells[j].l > basis.shells[i].l) std::swap(a, b);
            sp.a = a;
            sp.b = b;
            const Shell& A = basis.shells[a];
            const Shell& B = basis.shells[b];
            sp.la = A.l;
            sp.lb = B.l;
            double AB2 = 0.0;
            for (int d = 0; d < 3; ++d) {
                sp.AB[d] = A.center[d] - B.center[d];
                AB2 += sp.AB[d] * sp.AB[d];
            }
            for (std::size_t pa = 0; pa < A.exponents.size(); ++pa)
                for (std::size_t pb = 0; pb < B.exponents.size(); ++pb) {
                    const double x = A.exponents[pa], y = B.exponents[pb], z = x + y;
                    PrimitivePair pp;
                    pp.zeta = z;
                    pp.K = A.coefficients[pa] * B.coefficients[pb] * std::exp(-x * y / z * AB2);
                    // Negligible primitive pairs (|overlap-like magnitude| < 1e-17) are dropped.
                    if (std::fabs(pp.K) * std::pow(kPiExact / z, 1.5) < 1.0e-17) continue;
                    for (int d = 0; d < 3; ++d) {
                        pp.P[d] = (x * A.center[d] + y * B.center[d]) / z;
                        pp.PA[d] = pp.P[d] - A.center[d];
                    }
                    sp.prims.push_back(pp);
                }
        }
    return list;
}

ShellPair make_shell_pair(const Shell& X, const ShellTransform* tx, const Shell& Y, const ShellTransform* ty,
                          bool keep_order, int index_x, int index_y) {
    ShellPair sp;
    const bool swap = Y.l > X.l;
    const Shell& A = swap ? Y : X;
    const Shell& B = swap ? X : Y;
    sp.a = swap ? index_y : index_x;
    sp.b = swap ? index_x : index_y;
    sp.la = A.l;
    sp.lb = B.l;
    sp.ta = swap ? ty : tx;
    sp.tb = swap ? tx : ty;
    sp.unswap = swap && keep_order;
    double AB2 = 0.0;
    for (int d = 0; d < 3; ++d) {
        sp.AB[d] = A.center[d] - B.center[d];
        AB2 += sp.AB[d] * sp.AB[d];
    }
    for (std::size_t pa = 0; pa < A.exponents.size(); ++pa)
        for (std::size_t pb = 0; pb < B.exponents.size(); ++pb) {
            const double x = A.exponents[pa], y = B.exponents[pb], z = x + y;
            PrimitivePair pp;
            pp.zeta = z;
            pp.K = A.coefficients[pa] * B.coefficients[pb] * std::exp(-x * y / z * AB2);
            if (std::fabs(pp.K) * std::pow(kPiExact / z, 1.5) < 1.0e-17) continue;
            for (int d = 0; d < 3; ++d) {
                pp.P[d] = (x * A.center[d] + y * B.center[d]) / z;
                pp.PA[d] = pp.P[d] - A.center[d];
            }
            sp.prims.push_back(pp);
        }
    return sp;
}

const ShellTransform& spherical_transform(int l) {
    if (l < 0 || l > kMaxAngularMomentum) throw std::invalid_argument("spherical_transform: l out of range");
    return transform_tables().sph[l];
}

void apply_transform(const ShellTransform& T, int pre, int post, const double* in, double* out) {
    transform_rows(T.rows, ncart(T.l), pre, post, in, out);
}

void overlap_cartesian(const Shell& A, const Shell& B, std::vector<double>& out) {
    overlap_kinetic_block(A, B, OneElectron::Overlap, out);
}

void kinetic_cartesian(const Shell& A, const Shell& B, std::vector<double>& out) {
    overlap_kinetic_block(A, B, OneElectron::Kinetic, out);
}

void nuclear_cartesian(const Shell& A, const Shell& B, const std::vector<Atom>& nuclei, std::vector<double>& out) {
    std::vector<double> work;
    nuclear_block(A, B, nuclei, out, work);
}

void dipole_cartesian(const Shell& A, const Shell& B, const std::array<double, 3>& origin, std::vector<double>& out) {
    // (x - O_x) = (x - A_x) + (A_x - O_x): <a|x - O_x|b> = <a + 1_x|b> + (A_x - O_x) <a|b>.
    const CartTables& T = cart_tables();
    const int na = ncart(A.l), nb = ncart(B.l);
    std::vector<double> s0, s1;
    overlap_kinetic_block(A, B, OneElectron::Overlap, s0);
    Shell up = A;
    up.l = A.l + 1;
    overlap_kinetic_block(up, B, OneElectron::Overlap, s1);
    out.assign(static_cast<std::size_t>(3) * na * nb, 0.0);
    for (int k = 0; k < 3; ++k)
        for (int ca = 0; ca < na; ++ca) {
            const int cu = T.up[A.l][ca][k];
            for (int cb = 0; cb < nb; ++cb)
                out[(static_cast<std::size_t>(k) * na + ca) * nb + cb] =
                    s1[static_cast<std::size_t>(cu) * nb + cb] +
                    (A.center[k] - origin[k]) * s0[static_cast<std::size_t>(ca) * nb + cb];
        }
}

std::vector<Matrix> dipole_matrices(const BasisSet& basis, const std::array<double, 3>& origin) {
    const int n = basis.n_functions;
    std::vector<Matrix> M(3, Matrix(n, n));
    std::vector<double> cart, tmp, sph;
    for (std::size_t i = 0; i < basis.shells.size(); ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            const Shell& A = basis.shells[i];
            const Shell& B = basis.shells[j];
            dipole_cartesian(A, B, origin, cart);
            const int na = ncart(A.l), sa = A.size(), sb = B.size();
            tmp.resize(static_cast<std::size_t>(3) * na * sb);
            sph.resize(static_cast<std::size_t>(3) * sa * sb);
            to_spherical(B.l, 3 * na, 1, cart.data(), tmp.data());
            for (int k = 0; k < 3; ++k)
                to_spherical(A.l, 1, sb, tmp.data() + static_cast<std::size_t>(k) * na * sb,
                             sph.data() + static_cast<std::size_t>(k) * sa * sb);
            for (int k = 0; k < 3; ++k)
                for (int a = 0; a < sa; ++a)
                    for (int b = 0; b < sb; ++b) {
                        const double v = sph[(static_cast<std::size_t>(k) * sa + a) * sb + b];
                        M[k](A.first_function + a, B.first_function + b) = v;
                        M[k](B.first_function + b, A.first_function + a) = v;
                    }
        }
    return M;
}

const double* EriEngine::compute(const ShellPair& bra, const ShellPair& ket) {
    const CartTables& T = cart_tables();
    const int la = bra.la, lb = bra.lb, lc = ket.la, ld = ket.lb;
    const int Lab = la + lb, Lcd = lc + ld, L = Lab + Lcd;

    // ---- VRR layout: blocks (E, F) of ncart(F) x nm x ncart(E) values, i.e. [d][m][c]
    // with the bra component innermost (contiguous loops over c).
    std::size_t off[kMaxLL + 1][kMaxLL + 1];
    int nmv[kMaxLL + 1][kMaxLL + 1];
    int emin[kMaxLL + 1];
    std::size_t pos = 0;
    for (int F = 0; F <= Lcd; ++F) {
        emin[F] = F == 0 ? 0 : std::max(0, la - (Lcd - F));
        for (int E = emin[F]; E <= Lab; ++E) {
            const int nm = F == 0 ? L - E + 1 : Lcd - F + 1;
            off[E][F] = pos;
            nmv[E][F] = nm;
            pos += static_cast<std::size_t>(ncart(E)) * ncart(F) * nm;
        }
    }
    if (vrr_.size() < pos) vrr_.resize(pos);
    double* V = vrr_.data();

    // ---- contracted accumulator: for E = la..Lab blocks [c][k], k over (F = lc..Lcd, d).
    int kof[kMaxLL + 2];
    int NK = 0;
    for (int F = lc; F <= Lcd; ++F) {
        kof[F] = NK;
        NK += ncart(F);
    }
    std::size_t accoff[kMaxLL + 2], apos = 0;
    for (int E = la; E <= Lab; ++E) {
        accoff[E - la] = apos;
        apos += static_cast<std::size_t>(ncart(E)) * NK;
    }
    acc_.assign(apos, 0.0);
    double* acc = acc_.data();

    double Fm[kMaxM + 2];
    const double two_pi_52 = 2.0 * std::pow(kPiExact, 2.5);
    for (const PrimitivePair& p : bra.prims)
        for (const PrimitivePair& q : ket.prims) {
            const double zeta = p.zeta, eta = q.zeta, ze = zeta + eta;
            const double rho = zeta * eta / ze;
            const double PQ[3] = {p.P[0] - q.P[0], p.P[1] - q.P[1], p.P[2] - q.P[2]};
            const double Tb = rho * (PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2]);
            const double pref = two_pi_52 / (zeta * eta * std::sqrt(ze)) * p.K * q.K;
            boys_function(L, Tb, Fm);
            if (L == 0) {
                acc[0] += pref * Fm[0];
                continue;
            }
            for (int m = 0; m <= L; ++m) V[m] = pref * Fm[m];
            double WP[3], WQ[3];
            for (int d = 0; d < 3; ++d) {
                const double W = (zeta * p.P[d] + eta * q.P[d]) / ze;
                WP[d] = W - p.P[d];
                WQ[d] = W - q.P[d];
            }
            const double o2z = 0.5 / zeta, roz = rho / zeta, o2e = 0.5 / eta, roe = rho / eta, o2ze = 0.5 / ze;

            // Step 1: [E0|00]^(m), E = 1..Lab; blocks [m][c].
            for (int E = 1; E <= Lab; ++E) {
                const int nm = nmv[E][0];
                const int ne = ncart(E), ne1 = ncart(E - 1), ne2 = E >= 2 ? ncart(E - 2) : 0;
                double* tgt = V + off[E][0];
                const double* s1 = V + off[E - 1][0];
                const double* s2 = E >= 2 ? V + off[E - 2][0] : nullptr;
                for (int c = 0; c < ne; ++c) {
                    const int i = T.dir[E][c];
                    const int pc = T.down[E][c][i];
                    const int ai = T.comp[E][c][i] - 1;
                    const double pa = p.PA[i], wp = WP[i];
                    if (ai > 0) {
                        const int pp = T.down[E - 1][pc][i];
                        const double f = ai * o2z;
                        for (int m = 0; m < nm; ++m)
                            tgt[m * ne + c] = pa * s1[m * ne1 + pc] + wp * s1[(m + 1) * ne1 + pc] +
                                              f * (s2[m * ne2 + pp] - roz * s2[(m + 1) * ne2 + pp]);
                    } else {
                        for (int m = 0; m < nm; ++m) tgt[m * ne + c] = pa * s1[m * ne1 + pc] + wp * s1[(m + 1) * ne1 + pc];
                    }
                }
            }
            // Step 2: [E0|F0]^(m), F = 1..Lcd; blocks [d][m][c].
            for (int F = 1; F <= Lcd; ++F) {
                const int nf = ncart(F);
                const int nm = Lcd - F + 1;
                for (int E = emin[F]; E <= Lab; ++E) {
                    const int ne = ncart(E);
                    double* tgt = V + off[E][F];
                    const double* s1 = V + off[E][F - 1];
                    const int nm1 = nmv[E][F - 1];
                    const double* s2 = F >= 2 ? V + off[E][F - 2] : nullptr;
                    const int nm2 = F >= 2 ? nmv[E][F - 2] : 0;
                    const double* s3 = E >= 1 ? V + off[E - 1][F - 1] : nullptr;
                    const int nm3 = E >= 1 ? nmv[E - 1][F - 1] : 0;
                    const int ne1 = E >= 1 ? ncart(E - 1) : 0;
                    for (int d = 0; d < nf; ++d) {
                        const int i = T.dir[F][d];
                        const int qd = T.down[F][d][i];
                        const int fi = T.comp[F][d][i] - 1;
                        const double qc = q.PA[i], wq = WQ[i];
                        const int qq = fi > 0 ? T.down[F - 1][qd][i] : 0;
                        const double f2 = fi * o2e;
                        for (int m = 0; m < nm; ++m) {
                            double* t = tgt + (static_cast<std::size_t>(d) * nm + m) * ne;
                            const double* x1 = s1 + (static_cast<std::size_t>(qd) * nm1 + m) * ne;
                            const double* x1p = x1 + ne;
                            if (fi > 0) {
                                const double* x2 = s2 + (static_cast<std::size_t>(qq) * nm2 + m) * ne;
                                const double* x2p = x2 + ne;
                                for (int c = 0; c < ne; ++c)
                                    t[c] = qc * x1[c] + wq * x1p[c] + f2 * (x2[c] - roe * x2p[c]);
                            } else {
                                for (int c = 0; c < ne; ++c) t[c] = qc * x1[c] + wq * x1p[c];
                            }
                            if (E >= 1) {
                                const double* x3p = s3 + (static_cast<std::size_t>(qd) * nm3 + m + 1) * ne1;
                                for (const auto& lw : T.lowered[E][i]) t[lw[0]] += lw[2] * o2ze * x3p[lw[1]];
                            }
                        }
                    }
                }
            }
            // Accumulate the m = 0 values needed by the horizontal recurrences.
            for (int E = la; E <= Lab; ++E) {
                const int ne = ncart(E);
                double* a = acc + accoff[E - la];
                for (int F = lc; F <= Lcd; ++F) {
                    const int nf = ncart(F), nm = nmv[E][F];
                    const double* src = V + off[E][F];
                    for (int d = 0; d < nf; ++d) {
                        const double* sd = src + static_cast<std::size_t>(d) * nm * ne;
                        double* ad = a + kof[F] + d;
                        for (int c = 0; c < ne; ++c) ad[static_cast<std::size_t>(c) * NK] += sd[c];
                    }
                }
            }
        }

    // ---- bra HRR: (la lb | k), then spherical on b and a.
    const double* in[kMaxLL + 2];
    for (int E = la; E <= Lab; ++E) in[E - la] = acc + accoff[E - la];
    const ShellTransform& Ta = pair_transform(bra.ta, la);
    const ShellTransform& Tb = pair_transform(bra.tb, lb);
    const ShellTransform& Tc = pair_transform(ket.ta, lc);
    const ShellTransform& Td = pair_transform(ket.tb, ld);
    const int na = ncart(la), nb = ncart(lb), nsa = Ta.size(), nsb = Tb.size();
    const std::size_t bra_cart = static_cast<std::size_t>(na) * nb * NK;
    if (hrr_.size() < bra_cart) hrr_.resize(bra_cart);
    hrr(la, lb, bra.AB, NK, in, tmp_, hrr_.data());
    // [a][b][k] -> [a][sb][k] -> [sa][sb][k]
    const std::size_t s1size = static_cast<std::size_t>(na) * nsb * NK;
    const int nbra = nsa * nsb;
    const std::size_t s2size = static_cast<std::size_t>(nbra) * NK;
    std::vector<double>& buf = out_;
    if (buf.size() < s1size + s2size) buf.resize(s1size + s2size);
    transform_rows(Tb.rows, nb, na, NK, hrr_.data(), buf.data());
    transform_rows(Ta.rows, na, 1, nsb * NK, buf.data(), buf.data() + s1size);
    // transpose to [k][bra]
    const double* sb = buf.data() + s1size;
    if (hrr_.size() < s2size) hrr_.resize(s2size);
    double* kt = hrr_.data();
    for (int b = 0; b < nbra; ++b)
        for (int k = 0; k < NK; ++k) kt[static_cast<std::size_t>(k) * nbra + b] = sb[static_cast<std::size_t>(b) * NK + k];

    // ---- ket HRR on [k][bra] blocks.
    const double* kin[kMaxLL + 2];
    for (int F = lc; F <= Lcd; ++F) kin[F - lc] = kt + static_cast<std::size_t>(kof[F]) * nbra;
    const int nc = ncart(lc), nd = ncart(ld), nsc = Tc.size(), nsd = Td.size();
    const std::size_t ket_cart = static_cast<std::size_t>(nc) * nd * nbra;
    const std::size_t k1 = static_cast<std::size_t>(nc) * nsd * nbra;
    const std::size_t k2 = static_cast<std::size_t>(nsc) * nsd * nbra;
    if (acc_.size() < ket_cart + k1 + k2) acc_.resize(ket_cart + k1 + k2);
    double* kc = acc_.data();
    hrr(lc, ld, ket.AB, nbra, kin, tmp_, kc);
    transform_rows(Td.rows, nd, nc, nbra, kc, kc + ket_cart);
    transform_rows(Tc.rows, nc, 1, nsd * nbra, kc + ket_cart, kc + ket_cart + k1);
    // [ket][bra] -> [bra][ket] (and [b][a] / [d][c] for pairs that are to be unswapped)
    const double* ks = kc + ket_cart + k1;
    const int nket = nsc * nsd;
    if (buf.size() < static_cast<std::size_t>(nbra) * nket) buf.resize(static_cast<std::size_t>(nbra) * nket);
    double* out = buf.data();
    if (!bra.unswap && !ket.unswap) {
        for (int k = 0; k < nket; ++k)
            for (int b = 0; b < nbra; ++b)
                out[static_cast<std::size_t>(b) * nket + k] = ks[static_cast<std::size_t>(k) * nbra + b];
        return out;
    }
    for (int c = 0; c < nsc; ++c)
        for (int d = 0; d < nsd; ++d) {
            const int k = c * nsd + d;
            const int kd = ket.unswap ? d * nsc + c : k;
            for (int a = 0; a < nsa; ++a)
                for (int b = 0; b < nsb; ++b) {
                    const int bi = a * nsb + b;
                    const int bd = bra.unswap ? b * nsa + a : bi;
                    out[static_cast<std::size_t>(bd) * nket + kd] = ks[static_cast<std::size_t>(k) * nbra + bi];
                }
        }
    return out;
}

void compute_schwarz(const BasisSet& basis, ShellPairList& list) {
    EriEngine engine(basis);
    list.max_schwarz = 0.0;
    for (std::size_t p = 0; p < list.pairs.size(); ++p) {
        const ShellPair& sp = list.pairs[p];
        const int na = basis.shells[sp.a].size(), nb = basis.shells[sp.b].size();
        double mx = 0.0;
        if (!sp.prims.empty()) {
            const double* v = engine.compute(sp, sp);
            const std::size_t n = static_cast<std::size_t>(na) * nb;
            for (std::size_t ij = 0; ij < n; ++ij) mx = std::max(mx, std::fabs(v[ij * n + ij]));
        }
        list.schwarz[p] = std::sqrt(mx);
        list.max_schwarz = std::max(list.max_schwarz, list.schwarz[p]);
    }
}

int resolve_threads(int requested) {
    if (requested > 0) return std::min(requested, 64);
    const unsigned hw = std::thread::hardware_concurrency();
    return static_cast<int>(std::max(1u, std::min(hw, 64u)));
}

}  // namespace macqc
