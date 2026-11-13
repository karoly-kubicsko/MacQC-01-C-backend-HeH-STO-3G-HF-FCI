// Closed-shell (spin-adapted) CCSD and the (T) correction.
//
// The amplitude equations are the closed-shell CCSD equations in the form of
// Hirata/Stanton-Gauss used by PySCF's rccsd (intermediates F, L, W), written out with
// matrix products. The only (vv|vv) term, the particle-particle ladder
// sum_cd (ac|bd) tau_ij^cd, is delegated to a LadderOperator; the T1 parts of the
// dressed W_abcd are evaluated with (ia|bc). The (T) correction uses
//   E(T) = sum_ijk sum_abc (4 W_abc + W_bca + W_cab)(V_abc - V_cba) / (3 D_ijk^abc),
// W = P[(ia|fb) t_kj^cf - (ia|jm) t_mk^bc], V = W + (ia|jb) t_k^c + t_ij^ab f_ck + ...,
// looping over i >= j >= k and evaluating the energy for each distinct permutation of
// (i,j,k) from the same W and V (they are symmetric under simultaneous permutations).
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "macqc/cc.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/linalg.h"
#include "macqc/parallel.h"
#include "macqc/rhf.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// Row-major copy of a Matrix (rows x cols).
std::vector<double> to_vec(const Matrix& m) {
    std::vector<double> v(m.rows() * m.cols());
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) v[i * m.cols() + j] = m(i, j);
    return v;
}

// Columns [c0, c0 + nc) of an N x n coefficient matrix as a row-major N x nc array.
std::vector<double> columns(const Matrix& C, std::size_t c0, std::size_t nc) {
    std::vector<double> v(C.rows() * nc);
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < nc; ++j) v[i * nc + j] = C(i, c0 + j);
    return v;
}

// ---------------------------------------------------------------------------------------
// MO ladder (tests, small cases): W[(a,b)][(c,d)] = (ac|bd)
// ---------------------------------------------------------------------------------------

class MOLadder : public LadderOperator {
public:
    MOLadder(const MOIntegrals& mo, int nocc) : o_(static_cast<std::size_t>(nocc)), v_(static_cast<std::size_t>(mo.n - nocc)) {
        const std::size_t v = v_, o = o_;
        w_.resize(v * v * v * v);
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t b = 0; b < v; ++b)
                for (std::size_t c = 0; c < v; ++c)
                    for (std::size_t d = 0; d < v; ++d) w_[((a * v + b) * v + c) * v + d] = mo(o + a, o + c, o + b, o + d);
    }
    void apply(const std::vector<double>& tau, std::vector<double>& z) override {
        const std::size_t vv = v_ * v_;
        z.assign(o_ * o_ * vv, 0.0);
        gemm(false, true, o_ * o_, vv, vv, 1.0, tau.data(), vv, w_.data(), vv, 0.0, z.data(), vv);
    }
    std::string description() const override { return "MO (vv|vv) integrals"; }

private:
    std::size_t o_, v_;
    std::vector<double> w_;
};

// ---------------------------------------------------------------------------------------
// AO-driven ladder
// ---------------------------------------------------------------------------------------

constexpr double kLadderScreen = 1.0e-13;  // Schwarz bound x largest |tau| (AO) per shell block

class AOLadder : public LadderOperator {
public:
    AOLadder(const BasisSet& basis, const Matrix& C, int nocc, int threads, double eri_memory_mb)
        : basis_(basis), threads_(resolve_threads(threads)) {
        N_ = static_cast<std::size_t>(basis.n_functions);
        o_ = static_cast<std::size_t>(nocc);
        v_ = C.cols() - o_;
        cv_ = columns(C, o_, v_);
        pairs_ = make_shell_pairs(basis);
        compute_schwarz(basis, pairs_);
        const double unique = static_cast<double>(unique_eri_count(N_));
        if (eri_memory_mb > 0.0 && unique * 8.0 / (1024.0 * 1024.0) <= eri_memory_mb) {
            fock_ = std::make_unique<FockBuilder>(basis, threads_, eri_memory_mb);
            if (!fock_->stored()) fock_.reset();
        }
    }

    std::string description() const override {
        return fock_ ? "AO integrals in memory" : "AO integrals recomputed in every iteration (integral-direct)";
    }

    void apply(const std::vector<double>& tau, std::vector<double>& z) override {
        const std::size_t o = o_, v = v_, N = N_;
        // Symmetric and antisymmetric parts of tau_ij for i >= j: densities d = 0..nd-1.
        std::vector<std::pair<std::size_t, std::size_t>> sym, anti;
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j <= i; ++j) {
                sym.emplace_back(i, j);
                if (i != j) anti.emplace_back(i, j);
            }
        const std::size_t ns = sym.size(), na = anti.size(), nd = ns + na;
        // MO densities, back-transformed: D_AO = Cv M Cv^T, stored [mu nu][d].
        std::vector<double> dao(N * N * nd, 0.0);
        {
            std::vector<double> m(v * v), t(v * N), d(N * N);
            for (std::size_t k = 0; k < nd; ++k) {
                const bool is_sym = k < ns;
                const auto ij = is_sym ? sym[k] : anti[k - ns];
                const double* tij = tau.data() + (ij.first * o + ij.second) * v * v;
                for (std::size_t c = 0; c < v; ++c)
                    for (std::size_t e = 0; e < v; ++e)
                        m[c * v + e] = 0.5 * (tij[c * v + e] + (is_sym ? tij[e * v + c] : -tij[e * v + c]));
                gemm(false, true, v, N, v, 1.0, m.data(), v, cv_.data(), v, 0.0, t.data(), N);   // M Cv^T
                gemm(false, false, N, N, v, 1.0, cv_.data(), v, t.data(), N, 0.0, d.data(), N);  // Cv (M Cv^T)
                for (std::size_t x = 0; x < N * N; ++x) dao[x * nd + k] = d[x];
            }
        }
        // K-type contraction Kp[mu nu][d] over the unique AO integrals.
        std::vector<double> kp = fock_ ? contract_stored(dao, nd) : contract_direct(dao, nd);
        // Z_AO = Kp +- Kp^T, back to the MO basis, Z_ij = Z_S + Z_A, Z_ji = (Z_S - Z_A)^T.
        z.assign(o * o * v * v, 0.0);
        std::vector<double> zao(N * N), t(N * v), zs(v * v);
        std::vector<std::vector<double>> zsym(ns, std::vector<double>(v * v)), zanti(na, std::vector<double>(v * v));
        for (std::size_t k = 0; k < nd; ++k) {
            const double sign = k < ns ? 1.0 : -1.0;
            for (std::size_t mu = 0; mu < N; ++mu)
                for (std::size_t nu = 0; nu < N; ++nu)
                    zao[mu * N + nu] = kp[(mu * N + nu) * nd + k] + sign * kp[(nu * N + mu) * nd + k];
            gemm(false, false, N, v, N, 1.0, zao.data(), N, cv_.data(), v, 0.0, t.data(), v);  // Z Cv
            gemm(true, false, v, v, N, 1.0, cv_.data(), v, t.data(), v, 0.0, zs.data(), v);    // Cv^T Z Cv
            if (k < ns)
                zsym[k] = zs;
            else
                zanti[k - ns] = zs;
        }
        for (std::size_t k = 0; k < ns; ++k) {
            const std::size_t i = sym[k].first, j = sym[k].second;
            const double* s = zsym[k].data();
            const double* a = nullptr;
            if (i != j) {
                const std::size_t ka = static_cast<std::size_t>(
                    std::find(anti.begin(), anti.end(), sym[k]) - anti.begin());
                a = zanti[ka].data();
            }
            double* zij = z.data() + (i * o + j) * v * v;
            double* zji = z.data() + (j * o + i) * v * v;
            for (std::size_t x = 0; x < v; ++x)
                for (std::size_t y = 0; y < v; ++y) {
                    const double sv = s[x * v + y], av = a ? a[x * v + y] : 0.0;
                    zij[x * v + y] = sv + av;
                    if (i != j) zji[x * v + y] = sv - av;  // Z_ji^ab = Z_ij^ba
                }
        }
    }

private:
    // Kp[p r][d] += x D[q s][d], Kp[q r][d] += x D[p s][d], Kp[p s][d] += x D[q r][d],
    // Kp[q s][d] += x D[p r][d] for every unique (pq|rs), x = deg (pq|rs) / 8.
    static inline void update(double* kp, const double* d, std::size_t N, std::size_t nd, std::size_t p, std::size_t q,
                              std::size_t r, std::size_t s, double x) {
        double* k1 = kp + (p * N + r) * nd;
        double* k2 = kp + (q * N + r) * nd;
        double* k3 = kp + (p * N + s) * nd;
        double* k4 = kp + (q * N + s) * nd;
        const double* d1 = d + (q * N + s) * nd;
        const double* d2 = d + (p * N + s) * nd;
        const double* d3 = d + (q * N + r) * nd;
        const double* d4 = d + (p * N + r) * nd;
        for (std::size_t k = 0; k < nd; ++k) {
            k1[k] += x * d1[k];
            k2[k] += x * d2[k];
            k3[k] += x * d3[k];
            k4[k] += x * d4[k];
        }
    }

    std::vector<double> reduce(std::vector<std::vector<double>>& parts) const {
        std::vector<double> out = std::move(parts[0]);
        for (std::size_t t = 1; t < parts.size(); ++t)
            for (std::size_t x = 0; x < out.size(); ++x) out[x] += parts[t][x];
        return out;
    }

    std::vector<double> contract_stored(const std::vector<double>& dao, std::size_t nd) const {
        const std::size_t N = N_;
        const std::vector<double>& eri = fock_->packed();
        std::vector<std::vector<double>> parts(static_cast<std::size_t>(threads_), std::vector<double>(N * N * nd, 0.0));
        std::atomic<std::size_t> next{0};
        run_parallel(threads_, [&](int t) {
            double* kp = parts[static_cast<std::size_t>(t)].data();
            for (std::size_t k; (k = next++) < N;) {
                check_cancelled();
                const std::size_t p = N - 1 - k;
                for (std::size_t q = 0; q <= p; ++q) {
                    const std::size_t pq = p * (p + 1) / 2 + q;
                    const double* g = &eri[pq * (pq + 1) / 2];
                    for (std::size_t r = 0; r <= p; ++r) {
                        const std::size_t smax = (r == p) ? q : r;
                        for (std::size_t s = 0; s <= smax; ++s) {
                            double x = *g++;
                            if (x == 0.0) continue;
                            const std::size_t rs = r * (r + 1) / 2 + s;
                            double deg = 1.0;
                            if (p != q) deg *= 2.0;
                            if (r != s) deg *= 2.0;
                            if (pq != rs) deg *= 2.0;
                            update(kp, dao.data(), N, nd, p, q, r, s, x * deg / 8.0);
                        }
                    }
                }
            }
        });
        return reduce(parts);
    }

    std::vector<double> contract_direct(const std::vector<double>& dao, std::size_t nd) const {
        const std::size_t N = N_;
        const std::size_t nsh = basis_.shells.size();
        std::vector<double> dmax(nsh * nsh, 0.0);
        for (std::size_t a = 0; a < nsh; ++a)
            for (std::size_t b = 0; b < nsh; ++b) {
                const Shell& A = basis_.shells[a];
                const Shell& B = basis_.shells[b];
                double m = 0.0;
                for (int i = 0; i < A.size(); ++i)
                    for (int j = 0; j < B.size(); ++j) {
                        const double* d = dao.data() + (static_cast<std::size_t>(A.first_function + i) * N +
                                                        static_cast<std::size_t>(B.first_function + j)) * nd;
                        for (std::size_t k = 0; k < nd; ++k) m = std::max(m, std::fabs(d[k]));
                    }
                dmax[a * nsh + b] = m;
            }
        const std::size_t npairs = pairs_.pairs.size();
        std::vector<std::vector<double>> parts(static_cast<std::size_t>(threads_), std::vector<double>(N * N * nd, 0.0));
        std::atomic<std::size_t> next{0};
        run_parallel(threads_, [&](int t) {
            double* kp = parts[static_cast<std::size_t>(t)].data();
            EriEngine engine(basis_);
            for (std::size_t k; (k = next++) < npairs;) {
                check_cancelled();
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
                    const std::size_t sc = static_cast<std::size_t>(ket.a), sd = static_cast<std::size_t>(ket.b);
                    const double dm = std::max({dmax[sb * nsh + sd], dmax[sa * nsh + sc], dmax[sb * nsh + sc],
                                                dmax[sa * nsh + sd]});
                    if (bound * dm < kLadderScreen) continue;
                    const Shell& Cs = basis_.shells[sc];
                    const Shell& Ds = basis_.shells[sd];
                    const double* g = engine.compute(bra, ket);
                    const double deg = (sa == sb ? 1.0 : 2.0) * (sc == sd ? 1.0 : 2.0) * (P == Q ? 1.0 : 2.0);
                    for (int a = 0; a < A.size(); ++a) {
                        const std::size_t p = static_cast<std::size_t>(A.first_function + a);
                        for (int b = 0; b < B.size(); ++b) {
                            const std::size_t q = static_cast<std::size_t>(B.first_function + b);
                            for (int c = 0; c < Cs.size(); ++c) {
                                const std::size_t r = static_cast<std::size_t>(Cs.first_function + c);
                                for (int d = 0; d < Ds.size(); ++d) {
                                    const double x = *g++;
                                    if (x == 0.0) continue;
                                    const std::size_t s = static_cast<std::size_t>(Ds.first_function + d);
                                    update(kp, dao.data(), N, nd, p, q, r, s, x * deg / 8.0);
                                }
                            }
                        }
                    }
                }
            }
        });
        return reduce(parts);
    }

    const BasisSet& basis_;
    int threads_ = 1;
    std::size_t N_ = 0, o_ = 0, v_ = 0;
    std::vector<double> cv_;  // N x v virtual coefficients
    ShellPairList pairs_;
    std::unique_ptr<FockBuilder> fock_;  // stored integrals, if they fit
};

// ---------------------------------------------------------------------------------------
// Integral blocks
// ---------------------------------------------------------------------------------------

void allocate_blocks(CCIntegrals& I) {
    const std::size_t o = static_cast<std::size_t>(I.o), v = static_cast<std::size_t>(I.v);
    I.oooo.assign(o * o * o * o, 0.0);
    I.ovoo.assign(o * v * o * o, 0.0);
    I.ovov.assign(o * v * o * v, 0.0);
    I.oovv.assign(o * o * v * v, 0.0);
    I.ovvv.assign(o * v * v * v, 0.0);
}

// Fills the blocks from a function eri(p, q, r, s) = (pq|rs) over MO indices.
template <class Eri>
void fill_blocks(CCIntegrals& I, Eri eri) {
    const std::size_t o = static_cast<std::size_t>(I.o), v = static_cast<std::size_t>(I.v);
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j)
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t l = 0; l < o; ++l) I.oooo[((i * o + j) * o + k) * o + l] = eri(i, j, k, l);
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t k = 0; k < o; ++k) I.ovoo[((i * v + a) * o + j) * o + k] = eri(i, o + a, j, k);
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t b = 0; b < v; ++b) I.ovov[((i * v + a) * o + j) * v + b] = eri(i, o + a, j, o + b);
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j)
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t b = 0; b < v; ++b) I.oovv[((i * o + j) * v + a) * v + b] = eri(i, j, o + a, o + b);
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t b = 0; b < v; ++b)
                for (std::size_t c = 0; c <= b; ++c) {
                    const double x = eri(i, o + a, o + b, o + c);
                    I.ovvv[((i * v + a) * v + b) * v + c] = x;
                    I.ovvv[((i * v + a) * v + c) * v + b] = x;
                }
}

// ---------------------------------------------------------------------------------------
// CCSD residual (PySCF rccsd form)
// ---------------------------------------------------------------------------------------

struct Amplitudes {
    std::vector<double> t1, t2;  // [i][a], [i][j][a][b]
};

class ClosedShellCCSD {
public:
    ClosedShellCCSD(const CCIntegrals& ints, LadderOperator& ladder, int threads)
        : I(ints), L(ladder), o(static_cast<std::size_t>(ints.o)), v(static_cast<std::size_t>(ints.v)), threads_(threads) {
        eo.resize(o);
        ev.resize(v);
        for (std::size_t i = 0; i < o; ++i) eo[i] = I.fock(i, i);
        for (std::size_t a = 0; a < v; ++a) ev[a] = I.fock(o + a, o + a);
    }

    // Index helpers
    std::size_t T1(std::size_t i, std::size_t a) const { return i * v + a; }
    std::size_t T2(std::size_t i, std::size_t j, std::size_t a, std::size_t b) const { return ((i * o + j) * v + a) * v + b; }
    double ovov(std::size_t i, std::size_t a, std::size_t j, std::size_t b) const { return I.ovov[((i * v + a) * o + j) * v + b]; }
    double oovv(std::size_t i, std::size_t j, std::size_t a, std::size_t b) const { return I.oovv[((i * o + j) * v + a) * v + b]; }
    double ovoo(std::size_t i, std::size_t a, std::size_t j, std::size_t k) const { return I.ovoo[((i * v + a) * o + j) * o + k]; }
    double oooo(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const { return I.oooo[((i * o + j) * o + k) * o + l]; }
    double ovvv(std::size_t i, std::size_t a, std::size_t b, std::size_t c) const { return I.ovvv[((i * v + a) * v + b) * v + c]; }
    double fov(std::size_t i, std::size_t a) const { return I.fock(i, o + a); }

    double energy(const Amplitudes& t) const {
        double e = 0.0;
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t a = 0; a < v; ++a) e += 2.0 * fov(i, a) * t.t1[T1(i, a)];
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b) {
                        const double tau = t.t2[T2(i, j, a, b)] + t.t1[T1(i, a)] * t.t1[T1(j, b)];
                        e += tau * (2.0 * ovov(i, a, j, b) - ovov(i, b, j, a));
                    }
        return e;
    }

    // New amplitudes (Jacobi step) from the current ones.
    Amplitudes update(const Amplitudes& t) {
        const std::vector<double>& t1 = t.t1;
        const std::vector<double>& t2 = t.t2;
        const std::size_t ov = o * v, vv = v * v;

        // tau = t2 + t1 t1 (in ij ab order) and tau2 = t2 + t1(i,c) t1(l,d) analogue
        std::vector<double> tau(o * o * vv);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b) tau[T2(i, j, a, b)] = t2[T2(i, j, a, b)] + t1[T1(i, a)] * t1[T1(j, b)];

        // L(k,c,l,d) = 2 (kc|ld) - (kd|lc), stored [(k,l,d)][c]... used in F intermediates.
        // Foo[k][i] = foo + sum_lcd L(kcld) tau(i,l,c,d)
        std::vector<double> Foo(o * o), Fvv(vv), Fov(ov);
        for (std::size_t k = 0; k < o; ++k)
            for (std::size_t i = 0; i < o; ++i) {
                double s = I.fock(k, i);
                for (std::size_t l = 0; l < o; ++l)
                    for (std::size_t c = 0; c < v; ++c)
                        for (std::size_t d = 0; d < v; ++d)
                            s += (2.0 * ovov(k, c, l, d) - ovov(k, d, l, c)) * tau[T2(i, l, c, d)];
                Foo[k * o + i] = s;
            }
        // Fvv[a][c] = fvv - sum_kld L(kcld) tau(k,l,a,d): one matrix product over (k,l,d).
        {
            std::vector<double> A(v * o * o * v), B(o * o * v * v);
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t l = 0; l < o; ++l)
                    for (std::size_t d = 0; d < v; ++d) {
                        const std::size_t kld = (k * o + l) * v + d;
                        for (std::size_t a = 0; a < v; ++a) A[a * (o * o * v) + kld] = tau[T2(k, l, a, d)];
                        for (std::size_t c = 0; c < v; ++c) B[kld * v + c] = 2.0 * ovov(k, c, l, d) - ovov(k, d, l, c);
                    }
            gemm(false, false, v, v, o * o * v, -1.0, A.data(), o * o * v, B.data(), v, 0.0, Fvv.data(), v, threads_);
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t c = 0; c < v; ++c) Fvv[a * v + c] += I.fock(o + a, o + c);
        }
        for (std::size_t k = 0; k < o; ++k)
            for (std::size_t c = 0; c < v; ++c) {
                double s = fov(k, c);
                for (std::size_t l = 0; l < o; ++l)
                    for (std::size_t d = 0; d < v; ++d) s += (2.0 * ovov(k, c, l, d) - ovov(k, d, l, c)) * t1[T1(l, d)];
                Fov[k * v + c] = s;
            }
        // Loo, Lvv
        std::vector<double> Loo(Foo), Lvv(Fvv);
        for (std::size_t k = 0; k < o; ++k)
            for (std::size_t i = 0; i < o; ++i) {
                double s = 0.0;
                for (std::size_t c = 0; c < v; ++c) s += fov(k, c) * t1[T1(i, c)];
                for (std::size_t l = 0; l < o; ++l)
                    for (std::size_t c = 0; c < v; ++c) s += (2.0 * ovoo(l, c, k, i) - ovoo(k, c, l, i)) * t1[T1(l, c)];
                Loo[k * o + i] += s;
            }
        // G[a][c] = sum_kd (2 (kd|ac) - (kc|ad)) t1(k,d)   (also used by the T1 equation)
        std::vector<double> G(vv, 0.0);
        for (std::size_t k = 0; k < o; ++k)
            for (std::size_t d = 0; d < v; ++d) {
                const double tkd = t1[T1(k, d)];
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t c = 0; c < v; ++c) G[a * v + c] += 2.0 * tkd * ovvv(k, d, a, c);
            }
        for (std::size_t k = 0; k < o; ++k)
            for (std::size_t c = 0; c < v; ++c)
                for (std::size_t a = 0; a < v; ++a) {
                    double s = 0.0;
                    for (std::size_t d = 0; d < v; ++d) s += ovvv(k, c, a, d) * t1[T1(k, d)];
                    G[a * v + c] -= s;
                }
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t c = 0; c < v; ++c) {
                double s = G[a * v + c];
                for (std::size_t k = 0; k < o; ++k) s -= fov(k, c) * t1[T1(k, a)];
                Lvv[a * v + c] += s;
            }

        // ---- T1 equation ----
        std::vector<double> r1(ov, 0.0);
        std::vector<double> FooS(Foo), FvvS(Fvv);
        for (std::size_t i = 0; i < o; ++i) FooS[i * o + i] -= eo[i];
        for (std::size_t a = 0; a < v; ++a) FvvS[a * v + a] -= ev[a];
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t a = 0; a < v; ++a) {
                double s = fov(i, a);
                for (std::size_t k = 0; k < o; ++k)
                    for (std::size_t c = 0; c < v; ++c) {
                        s += -2.0 * fov(k, c) * t1[T1(k, a)] * t1[T1(i, c)];
                        s += 2.0 * Fov[k * v + c] * t2[T2(k, i, c, a)] - Fov[k * v + c] * t2[T2(i, k, c, a)];
                        s += Fov[k * v + c] * t1[T1(i, c)] * t1[T1(k, a)];
                        s += (2.0 * ovov(k, c, i, a) - oovv(k, i, a, c)) * t1[T1(k, c)];
                    }
                for (std::size_t c = 0; c < v; ++c) s += FvvS[a * v + c] * t1[T1(i, c)] + G[a * v + c] * t1[T1(i, c)];
                for (std::size_t k = 0; k < o; ++k) s -= FooS[k * o + i] * t1[T1(k, a)];
                for (std::size_t k = 0; k < o; ++k)
                    for (std::size_t l = 0; l < o; ++l)
                        for (std::size_t c = 0; c < v; ++c) {
                            const double x = -2.0 * ovoo(l, c, k, i) + ovoo(k, c, l, i);
                            s += x * (t2[T2(k, l, a, c)] + t1[T1(l, c)] * t1[T1(k, a)]);
                        }
                r1[T1(i, a)] = s;
            }
        // + sum_kcd (2 (kd|ac) - (kc|ad)) t2(i,k,c,d): for each (k,d) a [a][c] x [c][i] product
        {
            std::vector<double> U(v * o), Y(v * o), M(vv);
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t d = 0; d < v; ++d) {
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t c = 0; c < v; ++c) M[a * v + c] = 2.0 * ovvv(k, d, a, c) - ovvv(k, c, a, d);
                    for (std::size_t c = 0; c < v; ++c)
                        for (std::size_t i = 0; i < o; ++i) U[c * o + i] = t2[T2(i, k, c, d)];
                    gemm(false, false, v, o, v, 1.0, M.data(), v, U.data(), o, 0.0, Y.data(), o);
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t i = 0; i < o; ++i) r1[T1(i, a)] += Y[a * o + i];
                }
        }

        // ---- T2 equation ----
        std::vector<double> r2(o * o * vv, 0.0);
        std::vector<double> tmp(o * o * vv, 0.0);
        auto add_sym = [&](double f) {  // r2 += f (tmp + tmp^T(1,0,3,2))
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t j = 0; j < o; ++j)
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t b = 0; b < v; ++b)
                            r2[T2(i, j, a, b)] += f * (tmp[T2(i, j, a, b)] + tmp[T2(j, i, b, a)]);
        };
        // tmp_ijab = sum_c (ia|cb) t1(j,c) - sum_k t1(k,a) sum_c (ki|bc) t1(j,c)
        std::fill(tmp.begin(), tmp.end(), 0.0);
        {
            std::vector<double> X(vv * o);  // [(a,b)][j]
            std::vector<double> t1T(v * o);
            for (std::size_t c = 0; c < v; ++c)
                for (std::size_t j = 0; j < o; ++j) t1T[c * o + j] = t1[T1(j, c)];
            for (std::size_t i = 0; i < o; ++i) {
                gemm(false, false, vv, o, v, 1.0, I.ovvv.data() + i * v * vv, v, t1T.data(), o, 0.0, X.data(), o, threads_);
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b)
                        for (std::size_t j = 0; j < o; ++j) tmp[T2(i, j, a, b)] += X[(a * v + b) * o + j];
            }
            std::vector<double> Q(o * o * v * o);  // Q[k][i][b][j] = sum_c (ki|bc) t1(j,c)
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t i = 0; i < o; ++i)
                    for (std::size_t b = 0; b < v; ++b)
                        for (std::size_t j = 0; j < o; ++j) {
                            double s = 0.0;
                            for (std::size_t c = 0; c < v; ++c) s += oovv(k, i, b, c) * t1[T1(j, c)];
                            Q[((k * o + i) * v + b) * o + j] = s;
                        }
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t j = 0; j < o; ++j)
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t b = 0; b < v; ++b) {
                            double s = 0.0;
                            for (std::size_t k = 0; k < o; ++k) s += t1[T1(k, a)] * Q[((k * o + i) * v + b) * o + j];
                            tmp[T2(i, j, a, b)] -= s;
                        }
        }
        add_sym(1.0);
        // tmp_ijab = sum_k [sum_c (kc|ia) t1(j,c) + (ia|jk)] t1(k,b);  r2 -= tmp + tmp^T
        {
            std::vector<double> P(v * o * o * o);  // P[a][k][i][j]
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t k = 0; k < o; ++k)
                    for (std::size_t i = 0; i < o; ++i)
                        for (std::size_t j = 0; j < o; ++j) {
                            double s = ovoo(i, a, j, k);
                            for (std::size_t c = 0; c < v; ++c) s += ovov(k, c, i, a) * t1[T1(j, c)];
                            P[((a * o + k) * o + i) * o + j] = s;
                        }
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t j = 0; j < o; ++j)
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t b = 0; b < v; ++b) {
                            double s = 0.0;
                            for (std::size_t k = 0; k < o; ++k) s += P[((a * o + k) * o + i) * o + j] * t1[T1(k, b)];
                            tmp[T2(i, j, a, b)] = s;
                        }
        }
        add_sym(-1.0);
        // + (ia|jb)
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b) r2[T2(i, j, a, b)] += ovov(i, a, j, b);

        // Woooo and its product with tau
        {
            std::vector<double> W(o * o * o * o);  // [k][l][i][j]
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t l = 0; l < o; ++l)
                    for (std::size_t i = 0; i < o; ++i)
                        for (std::size_t j = 0; j < o; ++j) {
                            double s = oooo(k, i, l, j);
                            for (std::size_t c = 0; c < v; ++c)
                                s += ovoo(l, c, k, i) * t1[T1(j, c)] + ovoo(k, c, l, j) * t1[T1(i, c)];
                            for (std::size_t c = 0; c < v; ++c)
                                for (std::size_t d = 0; d < v; ++d) s += ovov(k, c, l, d) * tau[T2(i, j, c, d)];
                            W[((k * o + l) * o + i) * o + j] = s;
                        }
            // r2[ij][ab] += sum_kl W[kl][ij] tau[kl][ab]
            gemm(true, false, o * o, vv, o * o, 1.0, W.data(), o * o, tau.data(), vv, 1.0, r2.data(), vv, threads_);
        }
        // Ladder with the dressed W_abcd: Z(tau) - sum_k t1(k,b) X[k][a][ij] - sum_k t1(k,a) X[k][b][ji]
        {
            std::vector<double> z;
            L.apply(tau, z);
            for (std::size_t x = 0; x < r2.size(); ++x) r2[x] += z[x];
            // X[k][a][(i,j)] = sum_cd (kd|ac) tau(i,j,c,d)
            std::vector<double> X(o * v * o * o, 0.0), Td(v * o * o), Y(v * o * o);
            for (std::size_t d = 0; d < v; ++d) {
                for (std::size_t c = 0; c < v; ++c)
                    for (std::size_t i = 0; i < o; ++i)
                        for (std::size_t j = 0; j < o; ++j) Td[(c * o + i) * o + j] = tau[T2(i, j, c, d)];
                for (std::size_t k = 0; k < o; ++k) {
                    gemm(false, false, v, o * o, v, 1.0, I.ovvv.data() + (k * v + d) * vv, v, Td.data(), o * o, 0.0,
                         Y.data(), o * o);
                    double* xk = X.data() + k * v * o * o;
                    for (std::size_t x = 0; x < v * o * o; ++x) xk[x] += Y[x];
                }
            }
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t j = 0; j < o; ++j)
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t b = 0; b < v; ++b) {
                            double s = 0.0;
                            for (std::size_t k = 0; k < o; ++k)
                                s += t1[T1(k, b)] * X[(k * v + a) * o * o + i * o + j] +
                                     t1[T1(k, a)] * X[(k * v + b) * o * o + j * o + i];
                            r2[T2(i, j, a, b)] -= s;
                        }
        }
        // Lvv and Loo terms (diagonal Fock removed)
        for (std::size_t a = 0; a < v; ++a) Lvv[a * v + a] -= ev[a];
        for (std::size_t i = 0; i < o; ++i) Loo[i * o + i] -= eo[i];
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                gemm(false, false, v, v, v, 1.0, Lvv.data(), v, t2.data() + T2(i, j, 0, 0), v, 0.0,
                     tmp.data() + T2(i, j, 0, 0), v);
        add_sym(1.0);
        std::fill(tmp.begin(), tmp.end(), 0.0);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t k = 0; k < o; ++k) {
                    const double x = Loo[k * o + i];
                    const double* src = t2.data() + T2(k, j, 0, 0);
                    double* dst = tmp.data() + T2(i, j, 0, 0);
                    for (std::size_t ab = 0; ab < vv; ++ab) dst[ab] += x * src[ab];
                }
        add_sym(-1.0);

        // Wvoov[(i,a)][(k,c)] and Wvovo[(i,a)][(k,c)] (stored with the (i,a) row index)
        std::vector<double> Wvoov(ov * ov, 0.0), Wvovo(ov * ov, 0.0);
        {
            // ovvv x t1: Wvoov += sum_d (kc|ad) t1(i,d);  Wvovo += sum_d (kd|ac) t1(i,d)
            std::vector<double> Y(v * o), t1T(v * o), Z(o * vv);
            for (std::size_t d = 0; d < v; ++d)
                for (std::size_t i = 0; i < o; ++i) t1T[d * o + i] = t1[T1(i, d)];
            for (std::size_t k = 0; k < o; ++k) {
                for (std::size_t c = 0; c < v; ++c) {
                    gemm(false, false, v, o, v, 1.0, I.ovvv.data() + (k * v + c) * vv, v, t1T.data(), o, 0.0, Y.data(), o);
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t i = 0; i < o; ++i) Wvoov[(i * v + a) * ov + k * v + c] += Y[a * o + i];
                }
                // Z[i][(a,c)] = sum_d t1(i,d) (kd|ac)
                gemm(false, false, o, vv, v, 1.0, t1.data(), v, I.ovvv.data() + k * v * vv, vv, 0.0, Z.data(), vv);
                for (std::size_t i = 0; i < o; ++i)
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t c = 0; c < v; ++c) Wvovo[(i * v + a) * ov + k * v + c] += Z[i * vv + a * v + c];
            }
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t k = 0; k < o; ++k)
                        for (std::size_t c = 0; c < v; ++c) {
                            double s1 = ovov(k, c, i, a), s2 = oovv(k, i, a, c);
                            for (std::size_t l = 0; l < o; ++l) {
                                s1 -= ovoo(k, c, l, i) * t1[T1(l, a)];
                                s2 -= ovoo(l, c, k, i) * t1[T1(l, a)];
                            }
                            Wvoov[(i * v + a) * ov + k * v + c] += s1;
                            Wvovo[(i * v + a) * ov + k * v + c] += s2;
                        }
            // ovov x t2 parts: products over (l,d)
            std::vector<double> M1(ov * ov), N1(ov * ov), M2(ov * ov), N2(ov * ov), M3(ov * ov);
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t l = 0; l < o; ++l)
                        for (std::size_t d = 0; d < v; ++d) {
                            const std::size_t row = (i * v + a) * ov + l * v + d;
                            M1[row] = t2[T2(i, l, a, d)] - 0.5 * t2[T2(i, l, d, a)] - t1[T1(i, d)] * t1[T1(l, a)];
                            M2[row] = t2[T2(i, l, a, d)];
                            M3[row] = 0.5 * t2[T2(i, l, d, a)] + t1[T1(i, d)] * t1[T1(l, a)];
                        }
            for (std::size_t l = 0; l < o; ++l)
                for (std::size_t d = 0; d < v; ++d)
                    for (std::size_t k = 0; k < o; ++k)
                        for (std::size_t c = 0; c < v; ++c) {
                            N1[(l * v + d) * ov + k * v + c] = ovov(l, d, k, c);
                            N2[(l * v + d) * ov + k * v + c] = ovov(l, c, k, d);
                        }
            gemm(false, false, ov, ov, ov, 1.0, M1.data(), ov, N1.data(), ov, 1.0, Wvoov.data(), ov, threads_);
            gemm(false, false, ov, ov, ov, -0.5, M2.data(), ov, N2.data(), ov, 1.0, Wvoov.data(), ov, threads_);
            gemm(false, false, ov, ov, ov, -1.0, M3.data(), ov, N2.data(), ov, 1.0, Wvovo.data(), ov, threads_);
        }
        // T2 products with Wvoov and Wvovo; rows (i,a), columns (j,b), contraction (k,c)
        {
            std::vector<double> Ta(ov * ov), Tb(ov * ov), R(ov * ov), Wc(ov * ov);
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t c = 0; c < v; ++c)
                    for (std::size_t j = 0; j < o; ++j)
                        for (std::size_t b = 0; b < v; ++b) {
                            Ta[(k * v + c) * ov + j * v + b] = t2[T2(k, j, c, b)];
                            Tb[(k * v + c) * ov + j * v + b] = t2[T2(k, j, b, c)];
                        }
            for (std::size_t x = 0; x < ov * ov; ++x) Wc[x] = 2.0 * Wvoov[x] - Wvovo[x];
            // tmp = 2 Wvoov t2(kjcb) - Wvovo t2(kjcb) - Wvoov t2(kjbc)
            gemm(false, false, ov, ov, ov, 1.0, Wc.data(), ov, Ta.data(), ov, 0.0, R.data(), ov, threads_);
            gemm(false, false, ov, ov, ov, -1.0, Wvoov.data(), ov, Tb.data(), ov, 1.0, R.data(), ov, threads_);
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t j = 0; j < o; ++j)
                        for (std::size_t b = 0; b < v; ++b) tmp[T2(i, j, a, b)] = R[(i * v + a) * ov + j * v + b];
            add_sym(1.0);
            // tmp_ijab = sum_kc Wvovo[(i,b)][(k,c)] t2(k,j,a,c);  r2 -= tmp + tmp^T
            gemm(false, false, ov, ov, ov, 1.0, Wvovo.data(), ov, Tb.data(), ov, 0.0, R.data(), ov, threads_);
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t b = 0; b < v; ++b)
                    for (std::size_t j = 0; j < o; ++j)
                        for (std::size_t a = 0; a < v; ++a) tmp[T2(i, j, a, b)] = R[(i * v + b) * ov + j * v + a];
            add_sym(-1.0);
        }

        Amplitudes out;
        out.t1.resize(ov);
        out.t2.resize(o * o * vv);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t a = 0; a < v; ++a) out.t1[T1(i, a)] = r1[T1(i, a)] / (eo[i] - ev[a]);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b)
                        out.t2[T2(i, j, a, b)] = r2[T2(i, j, a, b)] / (eo[i] + eo[j] - ev[a] - ev[b]);
        return out;
    }

    // (T) correction.
    double triples(const Amplitudes& t) {
        const std::size_t vv = v * v, vvv = vv * v;
        const std::vector<double>& t1 = t.t1;
        const std::vector<double>& t2 = t.t2;
        std::vector<double> W(vvv), V(vvv), w(vvv), A(v * o);
        const std::size_t total = o * (o + 1) * (o + 2) / 6;
        std::size_t done = 0;
        double et = 0.0;
        // w_xyz[a'][b'][c'] = sum_f (x a'|b' f) t2(y,z,f,c') - sum_m (x a'|y m) t2(m,z,b',c')
        auto compute_w = [&](std::size_t x, std::size_t y, std::size_t z) {
            gemm(false, false, vv, v, v, 1.0, I.ovvv.data() + x * v * vv, v, t2.data() + T2(y, z, 0, 0), v, 0.0, w.data(),
                 v, threads_);
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t m = 0; m < o; ++m) A[a * o + m] = ovoo(x, a, y, m);
            // w[a][(b,c)] -= sum_m A[a][m] t2[m][z][(b,c)]: the rows of t2 for fixed z are not
            // contiguous over m, so gather them.
            std::vector<double> Bm(o * vv);
            for (std::size_t m = 0; m < o; ++m) std::copy_n(t2.data() + T2(m, z, 0, 0), vv, Bm.data() + m * vv);
            gemm(false, false, v, vv, o, -1.0, A.data(), o, Bm.data(), vv, 1.0, w.data(), vv, threads_);
        };
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j <= i; ++j)
                for (std::size_t k = 0; k <= j; ++k) {
                    check_cancelled();
                    std::fill(W.begin(), W.end(), 0.0);
                    const std::size_t ijk[3] = {i, j, k};
                    static const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
                    for (const auto& p : perms) {
                        compute_w(ijk[p[0]], ijk[p[1]], ijk[p[2]]);
                        // W[a][b][c] += w[(a,b,c) permuted]: w index position d holds abc[p[d]]
                        run_parallel(v > 32 ? threads_ : 1, [&](int th) {
                            const int nt = v > 32 ? threads_ : 1;
                            std::size_t idx[3];
                            for (std::size_t a = static_cast<std::size_t>(th); a < v; a += static_cast<std::size_t>(nt))
                                for (std::size_t b = 0; b < v; ++b)
                                    for (std::size_t c = 0; c < v; ++c) {
                                        const std::size_t abc[3] = {a, b, c};
                                        idx[0] = abc[p[0]];
                                        idx[1] = abc[p[1]];
                                        idx[2] = abc[p[2]];
                                        W[(a * v + b) * v + c] += w[(idx[0] * v + idx[1]) * v + idx[2]];
                                    }
                        });
                    }
                    // V = W + (ia|jb) t1(k,c) + (ia|kc) t1(j,b) + (jb|kc) t1(i,a)
                    //       + t2(i,j,a,b) f(c,k) + t2(i,k,a,c) f(b,j) + t2(j,k,b,c) f(a,i)
                    for (std::size_t a = 0; a < v; ++a)
                        for (std::size_t b = 0; b < v; ++b)
                            for (std::size_t c = 0; c < v; ++c) {
                                const std::size_t x = (a * v + b) * v + c;
                                V[x] = W[x] + ovov(i, a, j, b) * t1[T1(k, c)] + ovov(i, a, k, c) * t1[T1(j, b)] +
                                       ovov(j, b, k, c) * t1[T1(i, a)] + t2[T2(i, j, a, b)] * I.fock(o + c, k) +
                                       t2[T2(i, k, a, c)] * I.fock(o + b, j) + t2[T2(j, k, b, c)] * I.fock(o + a, i);
                            }
                    // Energy for each distinct permutation of (i,j,k): W_perm[x] = W[y], y[p[d]] = x[d].
                    const double eijk = eo[i] + eo[j] + eo[k];
                    std::vector<std::array<std::size_t, 3>> seen;
                    for (const auto& p : perms) {
                        const std::array<std::size_t, 3> trip{ijk[p[0]], ijk[p[1]], ijk[p[2]]};
                        if (std::find(seen.begin(), seen.end(), trip) != seen.end()) continue;
                        seen.push_back(trip);
                        std::vector<double> part(static_cast<std::size_t>(threads_), 0.0);
                        auto at = [&](const std::vector<double>& M, std::size_t x0, std::size_t x1, std::size_t x2) {
                            std::size_t y[3];
                            const std::size_t xs[3] = {x0, x1, x2};
                            y[p[0]] = xs[0];
                            y[p[1]] = xs[1];
                            y[p[2]] = xs[2];
                            return M[(y[0] * v + y[1]) * v + y[2]];
                        };
                        run_parallel(threads_, [&](int th) {
                            double s = 0.0;
                            for (std::size_t a = static_cast<std::size_t>(th); a < v; a += static_cast<std::size_t>(threads_))
                                for (std::size_t b = 0; b < v; ++b)
                                    for (std::size_t c = 0; c < v; ++c) {
                                        const double x = 4.0 * at(W, a, b, c) + at(W, b, c, a) + at(W, c, a, b);
                                        const double y = at(V, a, b, c) - at(V, c, b, a);
                                        s += x * y / (eijk - ev[a] - ev[b] - ev[c]);
                                    }
                            part[static_cast<std::size_t>(th)] = s;
                        });
                        for (double s : part) et += s / 3.0;
                    }
                    ++done;
                    report_progress(JobStage::Triples, static_cast<int>(100.0 * done / total), 0.0, 0.0);
                }
        return et;
    }

private:
    const CCIntegrals& I;
    LadderOperator& L;
    std::size_t o, v;
    int threads_;
    std::vector<double> eo, ev;
};

// Pulay DIIS on concatenated vectors.
class DIIS {
public:
    explicit DIIS(int size) : size_(std::max(2, size)) {}
    std::vector<double> extrapolate(const std::vector<double>& x, const std::vector<double>& err) {
        xs_.push_back(x);
        es_.push_back(err);
        if (static_cast<int>(xs_.size()) > size_) {
            xs_.erase(xs_.begin());
            es_.erase(es_.begin());
        }
        const std::size_t m = xs_.size();
        if (m < 2) return x;
        const std::size_t N = m + 1;
        std::vector<double> B(N * N, -1.0), rhs(N, 0.0);
        B[N * N - 1] = 0.0;
        rhs[m] = -1.0;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j <= i; ++j) {
                double s = 0.0;
                for (std::size_t k = 0; k < x.size(); ++k) s += es_[i][k] * es_[j][k];
                B[i * N + j] = B[j * N + i] = s;
            }
        for (std::size_t col = 0; col < N; ++col) {
            std::size_t piv = col;
            for (std::size_t r = col + 1; r < N; ++r)
                if (std::fabs(B[r * N + col]) > std::fabs(B[piv * N + col])) piv = r;
            if (std::fabs(B[piv * N + col]) < 1e-300) return x;
            if (piv != col) {
                for (std::size_t c = 0; c < N; ++c) std::swap(B[col * N + c], B[piv * N + c]);
                std::swap(rhs[col], rhs[piv]);
            }
            for (std::size_t r = col + 1; r < N; ++r) {
                const double f = B[r * N + col] / B[col * N + col];
                for (std::size_t c = col; c < N; ++c) B[r * N + c] -= f * B[col * N + c];
                rhs[r] -= f * rhs[col];
            }
        }
        std::vector<double> c(N, 0.0);
        for (std::size_t r = N; r-- > 0;) {
            double s = rhs[r];
            for (std::size_t k = r + 1; k < N; ++k) s -= B[r * N + k] * c[k];
            c[r] = s / B[r * N + r];
        }
        std::vector<double> out(x.size(), 0.0);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t k = 0; k < x.size(); ++k) out[k] += c[i] * xs_[i][k];
        return out;
    }

private:
    int size_;
    std::vector<std::vector<double>> xs_, es_;
};

}  // namespace

// ---------------------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------------------

CCIntegrals cc_integrals_from_mo(const MOIntegrals& mo, int nocc) {
    CCIntegrals I;
    I.o = nocc;
    I.v = mo.n - nocc;
    const std::size_t n = static_cast<std::size_t>(mo.n);
    I.fock = Matrix(n, n);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) {
            double f = mo.h(p, q);
            for (std::size_t i = 0; i < static_cast<std::size_t>(nocc); ++i) f += 2.0 * mo(p, q, i, i) - mo(p, i, i, q);
            I.fock(p, q) = f;
        }
    allocate_blocks(I);
    fill_blocks(I, [&](std::size_t p, std::size_t q, std::size_t r, std::size_t s) { return mo(p, q, r, s); });
    return I;
}

std::unique_ptr<LadderOperator> make_mo_ladder(const MOIntegrals& mo, int nocc) {
    return std::make_unique<MOLadder>(mo, nocc);
}

std::unique_ptr<LadderOperator> make_ao_ladder(const BasisSet& basis, const Matrix& C, int nocc, int threads,
                                               double eri_memory_mb) {
    return std::make_unique<AOLadder>(basis, C, nocc, threads, eri_memory_mb);
}

CCIntegrals cc_integrals_from_ao(const BasisSet& basis, const Matrix& h_ao, const Matrix& Cm, int nocc, int threads_req) {
    const std::size_t N = static_cast<std::size_t>(basis.n_functions);
    const std::size_t n = Cm.cols();
    const std::size_t o = static_cast<std::size_t>(nocc);
    const int threads = resolve_threads(threads_req);
    report_progress(JobStage::Transformation);
    const std::vector<double> C = columns(Cm, 0, n);
    const std::vector<double> Co = columns(Cm, 0, o);
    const std::size_t npair_ao = N * (N + 1) / 2;
    const std::size_t nrow = o * n;  // ket (i, s)
    std::vector<double> X(nrow * npair_ao, 0.0);

    // Step A: X[(i,s)][mu nu] = sum_{lambda sigma} C_lambda,i C_sigma,s (mu nu|lambda sigma)
    const ShellPairList pairs = make_shell_pairs(basis);
    const std::size_t npairs = pairs.pairs.size();
    {
        std::atomic<std::size_t> next{0};
        run_parallel(threads, [&](int) {
            EriEngine engine(basis);
            std::vector<double> G, t1, out;
            for (std::size_t k; (k = next++) < npairs;) {
                check_cancelled();
                const ShellPair& bra = pairs.pairs[npairs - 1 - k];
                if (bra.prims.empty()) continue;
                const Shell& A = basis.shells[static_cast<std::size_t>(bra.a)];
                const Shell& B = basis.shells[static_cast<std::size_t>(bra.b)];
                const std::size_t na = static_cast<std::size_t>(A.size()), nb = static_cast<std::size_t>(B.size()),
                                  nab = na * nb;
                G.assign(nab * N * N, 0.0);
                for (const ShellPair& ket : pairs.pairs) {
                    if (ket.prims.empty()) continue;
                    const Shell& Cs = basis.shells[static_cast<std::size_t>(ket.a)];
                    const Shell& Ds = basis.shells[static_cast<std::size_t>(ket.b)];
                    const double* val = engine.compute(bra, ket);
                    for (std::size_t ab = 0; ab < nab; ++ab) {
                        double* g = G.data() + ab * N * N;
                        for (int c = 0; c < Cs.size(); ++c) {
                            const std::size_t l = static_cast<std::size_t>(Cs.first_function + c);
                            for (int d = 0; d < Ds.size(); ++d) {
                                const std::size_t s = static_cast<std::size_t>(Ds.first_function + d);
                                const double x = *val++;
                                g[l * N + s] = x;
                                g[s * N + l] = x;
                            }
                        }
                    }
                }
                // t1[ab][lambda][s] = G C;  out[ab][i][s] = Co^T t1[ab]
                t1.resize(nab * N * n);
                gemm(false, false, nab * N, n, N, 1.0, G.data(), N, C.data(), n, 0.0, t1.data(), n);
                out.resize(nab * o * n);
                for (std::size_t ab = 0; ab < nab; ++ab)
                    gemm(true, false, o, n, N, 1.0, Co.data(), o, t1.data() + ab * N * n, n, 0.0, out.data() + ab * o * n, n);
                for (std::size_t a = 0; a < na; ++a)
                    for (std::size_t b = 0; b < nb; ++b) {
                        const std::size_t mu = static_cast<std::size_t>(A.first_function) + a;
                        const std::size_t nu = static_cast<std::size_t>(B.first_function) + b;
                        if (bra.a == bra.b && mu < nu) continue;
                        const std::size_t col = MOIntegrals::pair(mu, nu);
                        const double* src = out.data() + (a * nb + b) * o * n;
                        for (std::size_t r = 0; r < nrow; ++r) X[r * npair_ao + col] = src[r];
                    }
            }
        });
    }
    // Step B: each row (i,s): (pq|is) = sum_{mu nu} C_mu,p C_nu,q X[(i,s)][mu nu], p >= q, in place.
    {
        std::atomic<std::size_t> next{0};
        run_parallel(threads, [&](int) {
            std::vector<double> Y(N * N), t1(N * n), out(n * n);
            for (std::size_t row; (row = next++) < nrow;) {
                check_cancelled();
                double* x = X.data() + row * npair_ao;
                for (std::size_t mu = 0; mu < N; ++mu)
                    for (std::size_t nu = 0; nu <= mu; ++nu) {
                        const double val = x[mu * (mu + 1) / 2 + nu];
                        Y[mu * N + nu] = val;
                        Y[nu * N + mu] = val;
                    }
                gemm(false, false, N, n, N, 1.0, Y.data(), N, C.data(), n, 0.0, t1.data(), n);   // Y C
                gemm(true, false, n, n, N, 1.0, C.data(), n, t1.data(), n, 0.0, out.data(), n);  // C^T Y C
                for (std::size_t p = 0; p < n; ++p)
                    for (std::size_t q = 0; q <= p; ++q) x[p * (p + 1) / 2 + q] = out[p * n + q];
            }
        });
    }
    // (pq|is) with i occupied: X[(i * n + s)][pair(p,q)]
    auto get = [&](std::size_t p, std::size_t q, std::size_t i, std::size_t s) {
        return X[(i * n + s) * npair_ao + MOIntegrals::pair(p, q)];
    };
    // (pq|rs) with at least one of r, s occupied (or p, q by symmetry).
    auto eri = [&](std::size_t p, std::size_t q, std::size_t r, std::size_t s) {
        if (r < o) return get(p, q, r, s);
        if (s < o) return get(p, q, s, r);
        if (p < o) return get(r, s, p, q);
        if (q < o) return get(r, s, q, p);
        throw std::logic_error("cc_integrals_from_ao: (vv|vv) integrals are not transformed");
    };
    CCIntegrals I;
    I.o = nocc;
    I.v = static_cast<int>(n - o);
    // Fock matrix: f_pq = h_pq + sum_i 2 (pq|ii) - (pi|iq)
    {
        const std::vector<double> h = to_vec(h_ao);
        std::vector<double> t(N * n), hm(n * n);
        gemm(false, false, N, n, N, 1.0, h.data(), N, C.data(), n, 0.0, t.data(), n);
        gemm(true, false, n, n, N, 1.0, C.data(), n, t.data(), n, 0.0, hm.data(), n);
        I.fock = Matrix(n, n);
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) {
                double f = hm[p * n + q];
                for (std::size_t i = 0; i < o; ++i) f += 2.0 * get(p, q, i, i) - get(p, i, i, q);
                I.fock(p, q) = f;
            }
    }
    allocate_blocks(I);
    fill_blocks(I, eri);
    return I;
}

CCResult run_ccsd_closed_shell(const CCIntegrals& ints, LadderOperator& ladder, double reference_energy, bool triples,
                               const CCOptions& opt, const CCProgress& progress) {
    const auto t0 = Clock::now();
    const int threads = resolve_threads(opt.threads);
    ClosedShellCCSD cc(ints, ladder, threads);
    const std::size_t o = static_cast<std::size_t>(ints.o), v = static_cast<std::size_t>(ints.v);
    CCResult res;
    res.method = triples ? "CCSD(T)" : "CCSD";
    res.engine = "closed-shell, spin-adapted";
    res.ladder = ladder.description();
    res.occupied = ints.o;
    res.virtuals = ints.v;
    res.amplitude_tol = opt.amplitude_tol;
    res.energy_tol = opt.energy_tol;
    res.reference_energy = reference_energy;
    res.amplitude_counts = {0.0, static_cast<double>(o * v), static_cast<double>(o * o * v * v)};
    report_progress(JobStage::CoupledCluster);

    Amplitudes t;
    t.t1.assign(o * v, 0.0);
    t.t2.assign(o * o * v * v, 0.0);
    DIIS diis(opt.diis_size);
    double e_old = 0.0;
    for (int it = 1; it <= opt.max_iter; ++it) {
        check_cancelled();
        const double e = cc.energy(t);
        Amplitudes next = cc.update(t);
        double dmax = 0.0;
        std::vector<double> x, err;
        x.reserve(next.t1.size() + next.t2.size());
        x.insert(x.end(), next.t1.begin(), next.t1.end());
        x.insert(x.end(), next.t2.begin(), next.t2.end());
        err.resize(x.size());
        for (std::size_t k = 0; k < next.t1.size(); ++k) err[k] = next.t1[k] - t.t1[k];
        for (std::size_t k = 0; k < next.t2.size(); ++k) err[next.t1.size() + k] = next.t2[k] - t.t2[k];
        for (double d : err) dmax = std::max(dmax, std::fabs(d));
        CCIteration rec;
        rec.iteration = it;
        rec.correlation = e;
        rec.energy = reference_energy + e;
        rec.delta_e = e - e_old;
        rec.residual = dmax;
        res.history.push_back(rec);
        res.iterations = it;
        res.correlation_energy = e;
        report_progress(JobStage::CoupledCluster, it, rec.energy, dmax);
        if (progress) progress(rec);
        if (it > 1 && dmax < opt.amplitude_tol && std::fabs(rec.delta_e) < opt.energy_tol) {
            res.converged = true;
            break;
        }
        e_old = e;
        const std::vector<double> y = diis.extrapolate(x, err);
        std::copy(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(o * v), t.t1.begin());
        std::copy(y.begin() + static_cast<std::ptrdiff_t>(o * v), y.end(), t.t2.begin());
    }
    res.cc_energy = reference_energy + res.correlation_energy;
    res.total_energy = res.cc_energy;
    res.cc_seconds = seconds_since(t0);

    // Amplitude analysis (spatial amplitudes; t2 is the alpha-beta amplitude).
    double s1 = 0.0, s2 = 0.0;
    for (double x : t.t1) s1 += x * x;
    for (double x : t.t2) s2 += x * x;
    res.t1_diagnostic = std::sqrt(s1 / (2.0 * static_cast<double>(o)));
    res.rank_norms = {0.0, std::sqrt(s1), std::sqrt(s2)};
    struct Cand {
        double value;
        int rank;
        std::size_t i, j, a, b;
    };
    std::vector<Cand> cands;
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t a = 0; a < v; ++a)
            if (std::fabs(t.t1[i * v + a]) > 1e-6) cands.push_back({t.t1[i * v + a], 1, i, 0, a, 0});
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j)
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t b = 0; b < v; ++b) {
                    if (i < j || (i == j && a < b)) continue;  // t_ij^ab = t_ji^ba: one of each pair
                    const double x = t.t2[((i * o + j) * v + a) * v + b];
                    if (std::fabs(x) > 1e-6) cands.push_back({x, 2, i, j, a, b});
                }
    std::sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) { return std::fabs(x.value) > std::fabs(y.value); });
    for (std::size_t k = 0; k < cands.size() && k < 12; ++k) {
        const Cand& c = cands[k];
        CCAmplitude a;
        a.value = c.value;
        a.rank = c.rank;
        a.occupied.push_back(static_cast<int>(c.i) + 1);
        a.virtuals.push_back(static_cast<int>(o + c.a) + 1);
        a.spins = "a";
        if (c.rank == 2) {
            a.occupied.push_back(static_cast<int>(c.j) + 1);
            a.virtuals.push_back(static_cast<int>(o + c.b) + 1);
            a.spins = "ab";
        }
        res.leading.push_back(a);
    }

    if (triples && res.converged) {
        const auto t1s = Clock::now();
        report_progress(JobStage::Triples, 0);
        res.triples_correction = cc.triples(t);
        res.has_triples = true;
        res.total_energy = res.cc_energy + res.triples_correction;
        res.triples_seconds = seconds_since(t1s);
    }
    return res;
}

}  // namespace macqc
