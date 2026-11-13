// Closed-shell CCSD / CCSD(T) Lambda equations and response densities. The tensor
// expressions follow PySCF's ccsd_lambda.py, ccsd_rdm.py, ccsd_t_lambda_slow.py and
// ccsd_t_rdm_slow.py (same intermediates and index conventions), written with einsum on
// in-memory tensors; the (vv|vv) contractions use the ladder operator, and the vvvv and
// ovvv blocks of the two-particle density are accumulated slice by slice into the packed
// symmetric array, so no v^4 array is ever stored. The result is verified against finite
// differences of the energy and against PySCF's analytic gradients (test_gradients.cpp).
#include "macqc/cc_gradient.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <atomic>
#include <memory>
#include <stdexcept>

#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/parallel.h"
#include "macqc/tensor.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

using D = std::vector<std::size_t>;

Tensor from(const std::vector<double>& v, D dims) {
    Tensor t(dims);
    if (t.size() != v.size()) throw std::invalid_argument("cc_gradient: tensor size mismatch");
    std::copy(v.begin(), v.end(), t.data());
    return t;
}

// Pulay DIIS on concatenated vectors.
class Diis {
public:
    explicit Diis(int size) : size_(std::max(2, size)) {}
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

struct Ints {
    std::size_t o = 0, v = 0;
    Tensor foo, fov, fvo, fvv, eo, ev;
    Tensor ovoo, ovov, oovv, ovvv, oooo, voov, vvov;
};

Ints make_ints(const CCIntegrals& I) {
    Ints e;
    const std::size_t o = static_cast<std::size_t>(I.o), v = static_cast<std::size_t>(I.v);
    e.o = o;
    e.v = v;
    e.foo = Tensor({o, o});
    e.fov = Tensor({o, v});
    e.fvo = Tensor({v, o});
    e.fvv = Tensor({v, v});
    e.eo = Tensor({o});
    e.ev = Tensor({v});
    for (std::size_t p = 0; p < o; ++p) {
        e.eo[p] = I.fock(p, p);
        for (std::size_t q = 0; q < o; ++q) e.foo(p, q) = I.fock(p, q);
        for (std::size_t a = 0; a < v; ++a) {
            e.fov(p, a) = I.fock(p, o + a);
            e.fvo(a, p) = I.fock(o + a, p);
        }
    }
    for (std::size_t a = 0; a < v; ++a) {
        e.ev[a] = I.fock(o + a, o + a);
        for (std::size_t b = 0; b < v; ++b) e.fvv(a, b) = I.fock(o + a, o + b);
    }
    e.ovoo = from(I.ovoo, {o, v, o, o});
    e.ovov = from(I.ovov, {o, v, o, v});
    e.oovv = from(I.oovv, {o, o, v, v});
    e.ovvv = from(I.ovvv, {o, v, v, v});
    e.oooo = from(I.oooo, {o, o, o, o});
    e.voov = e.ovov.transpose({1, 0, 2, 3});  // (ai|jb)
    e.vvov = e.ovvv.transpose({2, 3, 0, 1});
    return e;
}

Tensor outer(const Tensor& a, const Tensor& b) { return einsum("ia,jb->ijab", a, b); }

// ---------------------------------------------------------------------------------------
// Lambda intermediates (PySCF ccsd_lambda.make_intermediates)
// ---------------------------------------------------------------------------------------

struct Imds {
    Tensor w1, w2, w3, w4, woooo, wvooo, wVOov, wvOOv, wvvov;
};

Imds make_intermediates(const Ints& e, const Tensor& t1, const Tensor& t2) {
    Imds m;
    const Tensor& fov = e.fov;
    m.w1 = e.fvv - einsum("ja,jb->ba", fov, t1);
    m.w2 = e.foo + einsum("ib,jb->ij", fov, t1);
    m.w3 = einsum("kc,jkbc->bj", fov, t2, 2.0) + fov.transpose({1, 0});
    m.w3 -= einsum("kc,kjbc->bj", fov, t2);
    m.w3 += einsum3("kc,kb,jc->bj", fov, t1, t1);
    m.w4 = fov;
    const Tensor tau = t2 + outer(t1, t1);

    m.w1 += einsum("jcba,jc->ba", e.ovvv, t1, 2.0);
    m.w1 -= einsum("jabc,jc->ba", e.ovvv, t1);
    {
        const Tensor theta = t2 * 2.0 - t2.transpose({1, 0, 2, 3});
        m.w3 += einsum("jkcd,kdcb->bj", theta, e.ovvv);
    }
    m.wVOov = einsum("jbcd,kd->bjkc", e.ovvv, t1);
    m.wvOOv = einsum("cbjd,kd->cjkb", e.vvov, t1, -1.0);
    {
        const Tensor g2vovv = e.vvov.transpose({0, 2, 1, 3}) * 2.0 - e.vvov.transpose({0, 2, 3, 1});
        m.wvooo = einsum("cibd,jkbd->ckij", g2vovv, tau);
    }
    m.wvvov = einsum("jabd,jkcd->abkc", e.ovvv, t2, -1.5);
    m.wvvov += e.vvov.transpose({0, 3, 2, 1}) * 2.0;
    m.wvvov -= e.vvov;
    {
        const Tensor g2vvov = e.vvov * 2.0 - e.ovvv.transpose({1, 2, 0, 3});
        const Tensor theta = t2 * 2.0 - t2.transpose({0, 1, 3, 2});
        const Tensor vackb = einsum("acjd,kjbd->ackb", g2vvov, theta);
        m.wvvov += vackb.transpose({0, 3, 2, 1});
        m.wvvov.axpy(-0.5, vackb);
    }
    m.w2 += einsum("kbij,kb->ij", e.ovoo, t1, 2.0);
    m.w2 -= einsum("ibkj,kb->ij", e.ovoo, t1);
    {
        const Tensor theta = t2.transpose({1, 0, 2, 3}) * 2.0 - t2;
        m.w3 -= einsum("lckj,klcb->bj", e.ovoo, theta);
    }
    {
        const Tensor tmp = einsum("lc,jcik->ijkl", t1, e.ovoo);
        m.woooo = tmp + tmp.transpose({1, 0, 3, 2});
    }
    m.wvOOv += einsum("lbjk,lc->bjkc", e.ovoo, t1);
    m.wVOov -= einsum("jbkl,lc->bjkc", e.ovoo, t1);
    m.wvooo += e.ovoo.transpose({1, 3, 2, 0}) * 2.0;
    m.wvooo -= e.ovoo.transpose({1, 0, 2, 3});
    m.wvooo -= einsum("klbc,iblj->ckij", t2, e.ovoo, 1.5);
    {
        const Tensor g2ovoo = e.ovoo * 2.0 - e.ovoo.transpose({2, 1, 0, 3});
        const Tensor theta = t2 * 2.0 - t2.transpose({1, 0, 2, 3});
        const Tensor vcjik = einsum("jlcb,lbki->cjki", theta, g2ovoo);
        m.wvooo += vcjik.transpose({0, 3, 2, 1});
        m.wvooo.axpy(-0.5, vcjik);
    }
    const Tensor& voov = e.voov;
    m.woooo += einsum("cijd,klcd->ijkl", voov, tau);
    const Tensor g2voov = voov * 2.0 - voov.transpose({0, 2, 1, 3});
    const Tensor tmpw4 = einsum("ckld,ld->kc", g2voov, t1);
    m.w1 -= einsum("ckja,kjcb->ba", g2voov, t2);
    m.w1 -= einsum("ja,jb->ba", tmpw4, t1);
    m.w2 += einsum("jkbc,bikc->ij", t2, g2voov);
    m.w2 += einsum("ib,jb->ij", tmpw4, t1);
    m.w3 += einsum3("kb,kc,jc->bj", t1, tmpw4, t1);
    m.w4 += tmpw4;
    m.wvOOv += einsum3("bljd,kd,lc->bjkc", voov, t1, t1);
    m.wVOov -= einsum3("bjld,kd,lc->bjkc", voov, t1, t1);
    Tensor VOov = einsum("bjld,klcd->bjkc", g2voov, t2);
    VOov -= einsum("bjld,kldc->bjkc", voov, t2);
    VOov += voov;
    Tensor vOOv = einsum("bljd,kldc->bjkc", voov, t2);
    vOOv -= e.oovv.transpose({2, 1, 0, 3});
    m.wVOov += VOov;
    m.wvOOv += vOOv;
    const Tensor ov1 = vOOv * 2.0 + VOov;
    const Tensor ov2 = VOov * 2.0 + vOOv;
    m.wvooo -= einsum("jb,bikc->ckij", t1, ov1);
    m.wvooo += einsum("kb,bijc->ckij", t1, ov2);
    m.w3 += einsum("ckjb,kc->bj", ov2, t1);
    m.wvvov += einsum("ajkc,jb->abkc", ov1, t1);
    m.wvvov -= einsum("ajkb,jc->abkc", ov2, t1);
    {
        const Tensor g2ovoo = e.ovoo * 2.0 - e.ovoo.transpose({2, 1, 0, 3});
        m.wvvov += einsum("laki,klbc->abic", g2ovoo, tau);
    }
    m.woooo += e.oooo.transpose({0, 2, 1, 3});
    m.w3 += einsum("bc,jc->bj", m.w1, t1);
    m.w3 -= einsum("kj,kb->bj", m.w2, t1);
    return m;
}

// ---------------------------------------------------------------------------------------
// (T): w and v of PySCF's ccsd_t_*_slow on full o^3 v^3 arrays
// ---------------------------------------------------------------------------------------

Tensor p6(const Tensor& t) {
    const Tensor t1 = t + t.transpose({0, 2, 1, 3, 5, 4});
    return t1 + t1.transpose({1, 0, 2, 4, 3, 5}) + t1.transpose({1, 2, 0, 4, 5, 3});
}

Tensor r6(const Tensor& w) {
    Tensor r = w * 4.0;
    r += w.transpose({0, 1, 2, 4, 5, 3});
    r += w.transpose({0, 1, 2, 5, 3, 4});
    r.axpy(-2.0, w.transpose({0, 1, 2, 5, 4, 3}));
    r.axpy(-2.0, w.transpose({0, 1, 2, 3, 5, 4}));
    r.axpy(-2.0, w.transpose({0, 1, 2, 4, 3, 5}));
    return r;
}

Tensor as_r6(const Tensor& m) {
    Tensor r = m * 2.0;
    r -= m.transpose({0, 1, 2, 5, 4, 3});
    r -= m.transpose({0, 1, 2, 3, 5, 4});
    return r;
}

struct TriplesWV {
    Tensor w, v;  // p6-symmetrised, divided by the denominators
};

TriplesWV triples_wv(const Ints& e, const Tensor& t1, const Tensor& t2) {
    const std::size_t o = e.o, v = e.v;
    TriplesWV r;
    Tensor w = einsum("iafb,kjcf->ijkabc", e.ovvv, t2);
    w -= einsum("iajm,mkbc->ijkabc", e.ovoo, t2);
    Tensor vv = einsum("iajb,kc->ijkabc", e.ovov, t1);
    vv += einsum("ck,ijab->ijkabc", e.fvo, t2);
    std::size_t x = 0;
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j)
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b)
                        for (std::size_t c = 0; c < v; ++c, ++x) {
                            const double d3 = e.eo[i] + e.eo[j] + e.eo[k] - e.ev[a] - e.ev[b] - e.ev[c];
                            w[x] /= d3;
                            vv[x] /= d3;
                        }
    r.w = p6(w);
    r.v = p6(vv);
    return r;
}

// All (T) contributions depend on t1 and t2 only: the Lambda source terms and the density
// additions of the CCSD(T) Lagrangian.
struct TriplesTerms {
    Tensor l1_t, l2_t;                             // Lambda sources (divided by the denominators)
    Tensor doo, dvv, dvo, dovov, dooov, dovvv;     // added to the CCSD density intermediates
};

void finish_lambda_sources(const Ints& e, TriplesTerms& T, const Tensor& J1, const Tensor& J2) {
    const std::size_t o = e.o, v = e.v;
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t a = 0; a < v; ++a) T.l1_t(i, a) /= e.eo[i] - e.ev[a];
    T.l2_t = J1 + J1.transpose({1, 0, 3, 2}) + J2;
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j)
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t b = 0; b < v; ++b)
                    T.l2_t(i, j, a, b) /= e.eo[i] - e.ev[a] + e.eo[j] - e.ev[b];
}

// Reference implementation on full o^3 v^3 arrays (PySCF ccsd_t_lambda_slow / ccsd_t_rdm_slow).
TriplesTerms triples_terms_full(const Ints& e, const Tensor& t1, const Tensor& t2) {
    TriplesTerms T;
    const TriplesWV wv = triples_wv(e, t1, t2);
    const Tensor rw = r6(wv.w);
    T.l1_t = einsum("jbkc,ijkabc->ia", e.ovov, rw, 0.5);
    const Tensor mm = as_r6(wv.w * 2.0 + wv.v * 0.5);
    Tensor J1 = einsum("kfbe,ijkaef->ijab", e.ovvv, mm);
    J1 -= einsum("ncmj,imnabc->ijab", e.ovoo, mm);
    const Tensor J2 = einsum("kc,ijkabc->ijab", e.fov, as_r6(wv.w));
    finish_lambda_sources(e, T, J1, J2);
    const Tensor wvh = wv.w + wv.v * 0.5;
    T.doo = einsum("iklabc,jklabc->ij", wvh, rw, -0.5);
    T.dvv = einsum("ijkacd,ijkbcd->ab", wvh, rw, 0.5);
    T.dvo = einsum("ijab,ijkabc->ck", t2, rw, 0.5);
    const Tensor rwv = r6(wv.w * 2.0 + wv.v * 0.5);
    T.dovov = einsum("kc,ijkabc->iajb", t1, rw, 0.5);
    T.dooov = einsum("mkbc,ijkabc->jmia", t2, rwv, -1.0);
    T.dovvv = einsum("kjcf,ijkabc->iafb", t2, rwv);
    return T;
}

// The same terms with v^3 blocks per ordered triple of occupied orbitals: the (k, l) pairs
// are distributed over the threads; for each, the blocks of all i are built from the raw
// w and v of the six permutations of (i, k, l).
TriplesTerms triples_terms_blocked(const Ints& e, const Tensor& t1, const Tensor& t2, int threads) {
    const std::size_t o = e.o, v = e.v, v2 = v * v, v3 = v2 * v;
    // ovvv_abf[p][a][b][f] = (pa|fb);  qt[k][b][e][f] = (kf|be)
    Tensor ovvv_abf({o, v, v, v}), qt({o, v, v, v});
    for (std::size_t p = 0; p < o; ++p)
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t b = 0; b < v; ++b)
                for (std::size_t f = 0; f < v; ++f) {
                    ovvv_abf[((p * v + a) * v + b) * v + f] = e.ovvv[((p * v + a) * v + f) * v + b];
                    qt[((p * v + a) * v + b) * v + f] = e.ovvv[((p * v + f) * v + a) * v + b];
                }
    const double* ovoo = e.ovoo.data();
    const double* T2 = t2.data();
    auto raw = [&](std::size_t p, std::size_t q, std::size_t r, double* W, double* V) {
        // W[a][b][c] = sum_f (pa|fb) t2[r][q][c][f] - sum_m (pa|qm) t2[m][r][b][c]
        gemm(false, true, v2, v, v, 1.0, ovvv_abf.data() + p * v3, v, T2 + (r * o + q) * v2, v, 0.0, W, v);
        gemm(false, false, v, v2, o, -1.0, ovoo + (p * v * o + q) * o, o * o, T2 + r * v2, o * v2, 1.0, W, v2);
        const double eijk = e.eo[p] + e.eo[q] + e.eo[r];
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t b = 0; b < v; ++b) {
                const double ovab = e.ovov[((p * v + a) * o + q) * v + b];
                const double tab = T2[((p * o + q) * v + a) * v + b];
                for (std::size_t c = 0; c < v; ++c) {
                    const std::size_t x = (a * v + b) * v + c;
                    const double d3 = eijk - e.ev[a] - e.ev[b] - e.ev[c];
                    W[x] /= d3;
                    V[x] = (ovab * t1(r, c) + e.fvo(c, r) * tab) / d3;
                }
            }
    };
    struct Acc {
        Tensor l1t, J1, J2, goo, gvv, dvo, dovov, dooov, dovvv;
    };
    const int nt = std::max(1, std::min(threads, static_cast<int>(o * o)));
    std::vector<Acc> acc(static_cast<std::size_t>(nt));
    for (Acc& A : acc) {
        A.l1t = Tensor({o, v});
        A.J1 = Tensor({o, o, v, v});
        A.J2 = Tensor({o, o, v, v});
        A.goo = Tensor({o, o});
        A.gvv = Tensor({v, v});
        A.dvo = Tensor({v, o});
        A.dovov = Tensor({o, v, o, v});
        A.dooov = Tensor({o, o, o, v});
        A.dovvv = Tensor({o, v, v, v});
    }
    std::atomic<std::size_t> next{0};
    run_parallel(nt, [&](int th) {
        Acc& A = acc[static_cast<std::size_t>(th)];
        std::vector<double> rw6(6 * v3), rv6(6 * v3), W(v3), V(v3), RW(v3), WVh(v3), RWV(v3), M(v3), AR(v3), tmp(v3);
        std::vector<double> WVall(o * v3), RWall(o * v3), tmp2(v2 * o);
        for (std::size_t pair; (pair = next++) < o * o;) {
            check_cancelled();
            const std::size_t J = pair / o, K = pair % o;
            for (std::size_t I = 0; I < o; ++I) {
                // raw blocks of the six permutations, p6 symmetrisation
                const std::size_t trip[6][3] = {{I, J, K}, {I, K, J}, {J, I, K}, {J, K, I}, {K, I, J}, {K, J, I}};
                for (int s = 0; s < 6; ++s)
                    raw(trip[s][0], trip[s][1], trip[s][2], rw6.data() + s * v3, rv6.data() + s * v3);
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b)
                        for (std::size_t c = 0; c < v; ++c) {
                            const std::size_t abc = (a * v + b) * v + c, acb = (a * v + c) * v + b,
                                              bac = (b * v + a) * v + c, bca = (b * v + c) * v + a,
                                              cab = (c * v + a) * v + b, cba = (c * v + b) * v + a;
                            W[abc] = rw6[abc] + rw6[v3 + acb] + rw6[2 * v3 + bac] + rw6[3 * v3 + bca] + rw6[4 * v3 + cab] +
                                     rw6[5 * v3 + cba];
                            V[abc] = rv6[abc] + rv6[v3 + acb] + rv6[2 * v3 + bac] + rv6[3 * v3 + bca] + rv6[4 * v3 + cab] +
                                     rv6[5 * v3 + cba];
                        }
                // r6, as_r6 on (a, b, c)
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b)
                        for (std::size_t c = 0; c < v; ++c) {
                            const std::size_t abc = (a * v + b) * v + c, acb = (a * v + c) * v + b,
                                              bac = (b * v + a) * v + c, bca = (b * v + c) * v + a,
                                              cab = (c * v + a) * v + b, cba = (c * v + b) * v + a;
                            auto w2v = [&](std::size_t x) { return 2.0 * W[x] + 0.5 * V[x]; };
                            RW[abc] = 4.0 * W[abc] + W[cab] + W[bca] - 2.0 * (W[cba] + W[acb] + W[bac]);
                            RWV[abc] = 4.0 * w2v(abc) + w2v(cab) + w2v(bca) - 2.0 * (w2v(cba) + w2v(acb) + w2v(bac));
                            M[abc] = 2.0 * w2v(abc) - w2v(cba) - w2v(acb);
                            AR[abc] = 2.0 * W[abc] - W[cba] - W[acb];
                            WVh[abc] = W[abc] + 0.5 * V[abc];
                        }
                // (a) l1t[I][a] += 1/2 sum_bc (Jb|Kc) RW[a][bc]
                for (std::size_t a = 0; a < v; ++a) {
                    double sum = 0.0;
                    for (std::size_t b = 0; b < v; ++b)
                        for (std::size_t c = 0; c < v; ++c) sum += e.ovov[((J * v + b) * o + K) * v + c] * RW[(a * v + b) * v + c];
                    A.l1t(I, a) += 0.5 * sum;
                }
                // (b1) J1[I][J][a][b] += sum_ef M[a][ef] (Kf|be)
                gemm(false, true, v, v, v2, 1.0, M.data(), v2, qt.data() + K * v3, v2, 1.0, A.J1.data() + (I * o + J) * v2, v);
                // (b2) J1[I][j'][a][b] -= sum_c M[ab][c] (Kc|J j')
                gemm(false, false, v2, o, v, 1.0, M.data(), v, ovoo + (K * v * o + J) * o, o * o, 0.0, tmp2.data(), o);
                for (std::size_t ab = 0; ab < v2; ++ab)
                    for (std::size_t jp = 0; jp < o; ++jp) A.J1[(I * o + jp) * v2 + ab] -= tmp2[ab * o + jp];
                // (b3) J2[I][J][a][b] += sum_c AR[ab][c] f(K, c)
                for (std::size_t ab = 0; ab < v2; ++ab) {
                    double sum = 0.0;
                    for (std::size_t c = 0; c < v; ++c) sum += AR[ab * v + c] * e.fov(K, c);
                    A.J2[(I * o + J) * v2 + ab] += sum;
                }
                // (d) gvv[a][b] += sum_cd WVh[a][cd] RW[b][cd]
                gemm(false, true, v, v, v2, 1.0, WVh.data(), v2, RW.data(), v2, 1.0, A.gvv.data(), v);
                // (e) dvo[c][K] += 1/2 sum_ab t2[I][J][ab] RW[ab][c]; (f) dovov[I][a][J][b] += 1/2 sum_c RW[ab][c] t1[K][c]
                for (std::size_t ab = 0; ab < v2; ++ab) {
                    const double t = T2[(I * o + J) * v2 + ab];
                    double sum = 0.0;
                    for (std::size_t c = 0; c < v; ++c) {
                        A.dvo(c, K) += 0.5 * t * RW[ab * v + c];
                        sum += RW[ab * v + c] * t1(K, c);
                    }
                    const std::size_t a = ab / v, b = ab % v;
                    A.dovov[((I * v + a) * o + J) * v + b] += 0.5 * sum;
                }
                // (g) dooov[J][m][I][a] -= sum_bc RWV[a][bc] t2[m][K][bc]
                gemm(false, true, v, o, v2, 1.0, RWV.data(), v2, T2 + K * v2, o * v2, 0.0, tmp2.data(), o);
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t m = 0; m < o; ++m) A.dooov[((J * o + m) * o + I) * v + a] -= tmp2[a * o + m];
                // (h) dovvv[I][a][f][b] += sum_c RWV[ab][c] t2[K][J][c][f]
                gemm(false, false, v2, v, v, 1.0, RWV.data(), v, T2 + (K * o + J) * v2, v, 0.0, tmp.data(), v);
                for (std::size_t a = 0; a < v; ++a)
                    for (std::size_t b = 0; b < v; ++b)
                        for (std::size_t f = 0; f < v; ++f) A.dovvv[((I * v + a) * v + f) * v + b] += tmp[(a * v + b) * v + f];
                std::copy(WVh.begin(), WVh.end(), WVall.begin() + static_cast<std::ptrdiff_t>(I * v3));
                std::copy(RW.begin(), RW.end(), RWall.begin() + static_cast<std::ptrdiff_t>(I * v3));
            }
            // (c) goo[i][j] += sum_abc WVh[i J K][abc] RW[j J K][abc]
            gemm(false, true, o, o, v3, 1.0, WVall.data(), v3, RWall.data(), v3, 1.0, A.goo.data(), o);
        }
    });
    TriplesTerms T;
    Tensor J1({o, o, v, v}), J2({o, o, v, v}), goo({o, o}), gvv({v, v});
    T.l1_t = Tensor({o, v});
    T.dvo = Tensor({v, o});
    T.dovov = Tensor({o, v, o, v});
    T.dooov = Tensor({o, o, o, v});
    T.dovvv = Tensor({o, v, v, v});
    for (const Acc& A : acc) {
        T.l1_t += A.l1t;
        J1 += A.J1;
        J2 += A.J2;
        goo += A.goo;
        gvv += A.gvv;
        T.dvo += A.dvo;
        T.dovov += A.dovov;
        T.dooov += A.dooov;
        T.dovvv += A.dovvv;
    }
    finish_lambda_sources(e, T, J1, J2);
    T.doo = goo * -0.5;
    T.dvv = gvv * 0.5;
    return T;
}

// ---------------------------------------------------------------------------------------
// Lambda update (PySCF ccsd_lambda.update_lambda)
// ---------------------------------------------------------------------------------------

void update_lambda(const Ints& e, LadderOperator& ladder, const Tensor& t1, const Tensor& t2, const Tensor& l1,
                   const Tensor& l2, const Imds& m, const TriplesTerms* triples, Tensor& l1new, Tensor& l2new) {
    const std::size_t o = e.o, v = e.v;
    const Tensor theta = t2 * 2.0 - t2.transpose({0, 1, 3, 2});
    const Tensor mba = einsum("klca,klcb->ba", l2, theta);
    const Tensor mij = einsum("ikcd,jkcd->ij", l2, theta);
    const Tensor mba1 = einsum("jc,jb->bc", l1, t1) + mba;
    const Tensor mij1 = einsum("kb,jb->kj", l1, t1) + mij;
    Tensor mia1 = t1 + einsum("kc,jkbc->jb", l1, t2, 2.0);
    mia1 -= einsum("kc,jkcb->jb", l1, t2);
    mia1 -= einsum3("ia,ja,jb->ib", t1, l1, t1);
    mia1 -= einsum("bd,jd->jb", mba, t1);
    mia1 -= einsum("lj,lb->jb", mij, t1);

    {
        std::vector<double> z;
        ladder.apply(l2.vec(), z);
        l2new = from(z, {o, o, v, v});
    }
    l1new = einsum("ijab,jb->ia", l2new, t1, 2.0);
    l1new -= einsum("jiab,jb->ia", l2new, t1);
    l2new *= 0.5;

    Tensor w1 = m.w1, w2 = m.w2;
    for (std::size_t a = 0; a < v; ++a) w1(a, a) -= e.ev[a];
    for (std::size_t i = 0; i < o; ++i) w2(i, i) -= e.eo[i];

    l1new += e.fov;
    l1new += einsum("ib,ba->ia", l1, w1);
    l1new -= einsum("ja,ij->ia", l1, w2);
    l1new -= einsum("ik,ka->ia", mij, m.w4);
    l1new -= einsum("ca,ic->ia", mba, m.w4);
    l1new += einsum("ijab,bj->ia", l2, m.w3, 2.0);
    l1new -= einsum("ijba,bj->ia", l2, m.w3);

    l2new += einsum("ia,jb->ijab", l1, m.w4);
    l2new += einsum("jibc,ca->jiba", l2, w1);
    l2new -= einsum("jk,kiba->jiba", w2, l2);

    l1new -= einsum("iajk,kj->ia", e.ovoo, mij1, 2.0);
    l1new += einsum("jaik,kj->ia", e.ovoo, mij1);
    l2new -= einsum("jbki,ka->jiba", e.ovoo, l1);

    const Tensor tau = t2 + outer(t1, t1);
    const Tensor l2tau = einsum("ijcd,klcd->ijkl", l2, tau);
    const Tensor l2t1 = einsum("jidc,kc->ijkd", l2, t1);

    l1new -= einsum("jb,jiab->ia", l1, e.oovv);
    l1new += einsum("iabc,bc->ia", e.ovvv, mba1, 2.0);
    l1new -= einsum("ibca,bc->ia", e.ovvv, mba1);
    l2new += einsum("jbac,ic->jiba", e.ovvv, l1);
    {
        const Tensor m4 = einsum("ijkd,kadb->ijab", l2t1, e.ovvv);
        l2new -= m4;
        l1new -= einsum("ijab,jb->ia", m4, t1, 2.0);
        l1new -= einsum("ijab,ia->jb", m4, t1, 2.0);
        l1new += einsum("jiab,jb->ia", m4, t1);
        l1new += einsum("jiab,ia->jb", m4, t1);
    }
    const Tensor& voov = e.voov;
    l1new += einsum("jb,aijb->ia", l1, voov, 2.0);
    l2new += voov.transpose({1, 2, 0, 3}) * 0.5;
    l2new -= einsum("bjic,ca->jiba", voov, mba1);
    l2new -= einsum("bjka,ik->jiba", voov, mij1);
    l1new += einsum("aijb,jb->ia", voov, mia1, 2.0);
    l1new -= einsum("bija,jb->ia", voov, mia1);
    {
        const Tensor m4 = einsum("ijkl,aklb->ijab", l2tau, voov);
        l2new.axpy(0.5, m4);
        l1new += einsum("ijab,jb->ia", m4, t1, 2.0);
        l1new -= einsum("ijba,jb->ia", m4, t1);
    }
    l1new -= einsum("ckij,jkca->ia", m.wvooo, l2);
    l1new += einsum("abkc,kibc->ia", m.wvvov, l2);
    {
        const Tensor tmp_voov = m.wVOov * 2.0 + m.wvOOv;
        const Tensor tmp = l2.transpose({0, 2, 1, 3}) - l2.transpose({0, 3, 1, 2}) * 0.5;
        l2new += einsum("iakc,bjkc->jiba", tmp, tmp_voov);
    }
    {
        const Tensor tmp = einsum("jkca,bikc->jiba", l2, m.wvOOv);
        l2new += tmp;
        l2new.axpy(0.5, tmp.transpose({1, 0, 2, 3}));
    }
    {
        const Tensor m3 = einsum("ijkl,klab->ijab", m.woooo, l2);
        l2new.axpy(0.5, m3);
        l1new += einsum("ijab,jb->ia", m3, t1, 2.0);
        l1new -= einsum("ijba,jb->ia", m3, t1);
    }
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t a = 0; a < v; ++a) l1new(i, a) /= e.eo[i] - e.ev[a];
    Tensor sym = l2new + l2new.transpose({1, 0, 3, 2});
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j)
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t b = 0; b < v; ++b)
                    sym(i, j, a, b) /= e.eo[i] - e.ev[a] + e.eo[j] - e.ev[b];
    l2new = sym;
    if (triples) {
        l1new += triples->l1_t;
        l2new += triples->l2_t;
    }
}

// ---------------------------------------------------------------------------------------
// densities (PySCF ccsd_rdm._gamma1/_gamma2 intermediates, _make_rdm1/_make_rdm2)
// ---------------------------------------------------------------------------------------

// Adds f * block[x...] at the MO positions (offsets per index) to the packed density.
void add_block(Packed4& G, const Tensor& B, std::size_t o0, std::size_t o1, std::size_t o2, std::size_t o3, double f) {
    const std::size_t d0 = B.dim(0), d1 = B.dim(1), d2 = B.dim(2), d3 = B.dim(3);
    std::size_t x = 0;
    for (std::size_t p = 0; p < d0; ++p)
        for (std::size_t q = 0; q < d1; ++q)
            for (std::size_t r = 0; r < d2; ++r)
                for (std::size_t s = 0; s < d3; ++s, ++x)
                    if (B[x] != 0.0) G.add_raw(o0 + p, o1 + q, o2 + r, o3 + s, f * B[x]);
}

void build_densities(const Ints& e, const Tensor& t1, const Tensor& t2, const Tensor& l1, const Tensor& l2,
                     const TriplesTerms* triples, CorrelatedDensity& out) {
    const std::size_t o = e.o, v = e.v, n = o + v;
    // ---- gamma1
    Tensor doo = einsum("ja,ia->ij", t1, l1, -1.0);
    Tensor dvv = einsum("ia,ib->ab", t1, l1);
    const Tensor xtv = einsum("ie,me->im", t1, l1);
    Tensor dvo = t1.transpose({1, 0}) - einsum("im,ma->ai", xtv, t1);
    {
        const Tensor theta = t2 * 2.0 - t2.transpose({0, 1, 3, 2});
        doo -= einsum("jkab,ikab->ij", theta, l2);
        dvv += einsum("jica,jicb->ab", theta, l2);
        const Tensor xt1 = einsum("mnef,inef->mi", l2, theta);
        const Tensor xt2 = einsum("mnaf,mnef->ea", l2, theta);
        dvo += einsum("imae,me->ai", theta, l1);
        dvo -= einsum("mi,ma->ai", xt1, t1);
        dvo -= einsum("ie,ae->ai", t1, xt2);
    }
    const Tensor& dov = l1;

    // ---- gamma2 (single block)
    const Tensor pvOOv = einsum("ikca,jkcb->aijb", l2, t2);
    Tensor moo = einsum("dljd->jl", pvOOv, 2.0);
    Tensor mvv = einsum("blld->db", pvOOv, 2.0);
    Tensor gooov = einsum("kc,cija->jkia", t1, pvOOv);
    Tensor pvoOV = einsum("ikca,jkbc->aijb", l2, t2, -1.0);
    {
        const Tensor theta = t2 * 2.0 - t2.transpose({0, 1, 3, 2});
        pvoOV += einsum("ikac,jkbc->aijb", l2, theta);
    }
    moo += einsum("dljd->jl", pvoOV);
    mvv += einsum("blld->db", pvoOV);
    gooov -= einsum("jc,cika->jkia", t1, pvoOV);
    const Tensor mia = einsum("kc,ikac->ia", l1, t2, 2.0) - einsum("kc,ikca->ia", l1, t2);
    const Tensor mab = einsum("kc,kb->cb", l1, t1);
    const Tensor mij = einsum("kc,jc->jk", l1, t1) + moo * 0.5;
    Tensor tau = outer(t1, t1) + t2;
    const Tensor goooo = einsum("ijab,klab->ijkl", tau, l2, 0.5);
    const Tensor doooo = goooo.transpose({0, 2, 1, 3}) * 2.0 - goooo.transpose({0, 3, 1, 2});
    gooov += einsum("ji,ka->jkia", moo, t1, -0.5);
    gooov += einsum("la,jkil->jkia", t1, goooo, 2.0);
    gooov -= einsum("ib,jkba->jkia", l1, tau);
    gooov -= einsum("jkba,ib->jkia", l2, t1);
    Tensor dooov = gooov.transpose({0, 2, 1, 3}) * 2.0 - gooov.transpose({1, 2, 0, 3});

    Tensor goovv = einsum("ia,jb->ijab", mia, t1);
    {
        Tensor tmp = einsum("ijkl,klab->ijab", goooo, tau);
        tmp -= einsum("jk,ikab->ijab", mij, tau);
        tmp -= einsum("cb,ijac->ijab", mab, t2);
        tmp -= einsum("bd,ijad->ijab", mvv * 0.5, tau);
        tmp.axpy(0.5, tau);
        tmp.axpy(0.5, l2);
        goovv += tmp;
    }
    Tensor gOvvO = einsum3("kiac,jc,kb->iabj", l2, t1, t1);
    gOvvO += einsum("aijb->iabj", pvOOv);
    Tensor govVO = einsum("ia,jb->iabj", l1, t1);
    govVO -= einsum3("ikac,jc,kb->iabj", l2, t1, t1);
    govVO += einsum("aijb->iabj", pvoOV);
    Tensor dovvo = govVO * 2.0 + gOvvO;
    Tensor doovv = (gOvvO * -2.0 - govVO).transpose({3, 0, 1, 2});
    {
        Tensor tau2 = tau - t2 * 0.5;
        goovv += einsum("dlib,jlda->ijab", pvOOv, tau2);
        goovv -= einsum("dlia,jldb->ijab", pvoOV, tau2);
        const Tensor tmp = pvoOV + pvOOv * 0.5;
        goovv += einsum("dlia,jlbd->ijab", tmp, t2);
    }
    Tensor dovov = goovv.transpose({0, 2, 1, 3}) * 2.0 - goovv.transpose({1, 2, 0, 3});

    // ovvv and vvvv: gvvvv[a] slices; gvovv accumulated in full (o v^3)
    Tensor gvovv({v, o, v, v});
    {
        const Tensor jabc = einsum("ijab,ic->jabc", l2, t1);  // [j][a][b][c]
        // t2 as a (ij) x (cd) matrix, l2 slice for a fixed a as (ij) x b
        Tensor l2a({o * o, v}), gv({v, v, v}), dv({v, v, v});
        Tensor t2m = t2;  // (ij)(cd) row-major already
        for (std::size_t a = 0; a < v; ++a) {
            check_cancelled();
            for (std::size_t ij = 0; ij < o * o; ++ij)
                for (std::size_t b = 0; b < v; ++b) l2a[ij * v + b] = l2[(ij * v + a) * v + b];
            // gvvvv[a][b][c][d] = sum_ij l2[ij a b] t2[ij c d] + sum_j jabc[j a b c] t1[j d]
            gemm(true, false, v, v * v, o * o, 1.0, l2a.data(), v, t2m.data(), v * v, 0.0, gv.data(), v * v);
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t b = 0; b < v; ++b)
                    for (std::size_t c = 0; c < v; ++c) {
                        const double x = jabc[((j * v + a) * v + b) * v + c];
                        if (x == 0.0) continue;
                        for (std::size_t d = 0; d < v; ++d) gv[(b * v + c) * v + d] += x * t1(j, d);
                    }
            // dvvvv[a][x][y][z] = g[a][y][x][z] - 1/2 g[a][y][z][x]; dm2_vvvv = 2 (d + d^T(1,0,3,2)),
            // i.e. 4 d in the symmetric storage (both terms lie in the same orbit).
            for (std::size_t x = 0; x < v; ++x)
                for (std::size_t y = 0; y < v; ++y)
                    for (std::size_t z = 0; z < v; ++z) {
                        const double dval = gv[(y * v + x) * v + z] - 0.5 * gv[(y * v + z) * v + x];
                        dv[(x * v + y) * v + z] = dval;
                    }
            for (std::size_t x = 0; x < v; ++x)
                for (std::size_t y = 0; y < v; ++y)
                    for (std::size_t z = 0; z < v; ++z) {
                        const double dval = dv[(x * v + y) * v + z];
                        if (dval != 0.0) out.G.add_raw(o + a, o + x, o + y, o + z, 0.5 * 4.0 * dval);
                    }
            // gvovv[a][i][b][c] = -sum_d g[a][d][b][c] t1[i][d]
            gemm(false, false, o, v * v, v, -1.0, t1.data(), v, gv.data(), v * v, 0.0,
                 gvovv.data() + a * o * v * v, v * v);
        }
    }
    gvovv += einsum("akic,kb->aibc", pvoOV, t1);
    gvovv -= einsum("akib,kc->aibc", pvOOv, t1);
    gvovv += einsum("ja,jibc->aibc", l1, t2);
    gvovv += einsum3("ja,jb,ic->aibc", l1, t1, t1);
    gvovv += einsum("ba,ic->aibc", mvv * 0.5, t1);
    gvovv += einsum("ja,jibc->aibc", t1, l2);
    Tensor dovvv = gvovv.transpose({1, 3, 0, 2}) * 2.0 - gvovv.transpose({1, 2, 0, 3});

    // ---- (T) contributions
    if (triples) {
        doo += triples->doo;
        dvv += triples->dvv;
        dvo += triples->dvo;
        dovov += triples->dovov;
        dooov += triples->dooov;
        dovvv += triples->dovvv;
    }

    // ---- assemble (PySCF _make_rdm1 / _make_rdm2; G = dm2 / 2, each block with the number
    //      of its symmetric copies)
    Matrix& Dm = out.D;
    Dm = Matrix(n, n);
    Tensor dm1c({n, n});
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j) dm1c(i, j) = doo(i, j) + doo(j, i);
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t a = 0; a < v; ++a) {
            dm1c(i, o + a) = dov(i, a) + dvo(a, i);
            dm1c(o + a, i) = dm1c(i, o + a);
        }
    for (std::size_t a = 0; a < v; ++a)
        for (std::size_t b = 0; b < v; ++b) dm1c(o + a, o + b) = dvv(a, b) + dvv(b, a);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) Dm(p, q) = dm1c(p, q) + (p == q && p < o ? 2.0 : 0.0);

    Packed4& G = out.G;
    const double h = 0.5;  // G = dm2 / 2
    add_block(G, dovov, 0, o, 0, o, 4.0 * h);
    add_block(G, doovv, 0, 0, o, o, 4.0 * h);
    add_block(G, dovvo, 0, o, o, 0, 4.0 * h);
    add_block(G, doooo, 0, 0, 0, 0, 4.0 * h);
    add_block(G, dovvv, 0, o, o, o, 4.0 * h);
    add_block(G, dooov, 0, 0, 0, o, 4.0 * h);
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) {
                const double d = dm1c(p, q);
                if (d == 0.0) continue;
                G.add_raw(i, i, p, q, h * 2.0 * d);
                G.add_raw(p, q, i, i, h * 2.0 * d);
                G.add_raw(p, i, i, q, -h * d);
                G.add_raw(i, q, p, i, -h * d);  // dm2[i,:,:,i] -= dm1.T: (i, q, p, i) gets dm1[p, q]
            }
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j) {
            G.add_raw(i, i, j, j, h * 4.0);
            G.add_raw(i, j, j, i, -h * 2.0);
        }
}

}  // namespace

double compare_triples_implementations(const CCIntegrals& ints, const std::vector<double>& t1v,
                                       const std::vector<double>& t2v, int threads) {
    const Ints e = make_ints(ints);
    const Tensor t1 = from(t1v, {e.o, e.v}), t2 = from(t2v, {e.o, e.o, e.v, e.v});
    const TriplesTerms a = triples_terms_full(e, t1, t2), b = triples_terms_blocked(e, t1, t2, resolve_threads(threads));
    double w = 0.0;
    auto cmp = [&](const Tensor& x, const Tensor& y) {
        for (std::size_t i = 0; i < x.size(); ++i) w = std::max(w, std::fabs(x[i] - y[i]));
    };
    cmp(a.l1_t, b.l1_t);
    cmp(a.l2_t, b.l2_t);
    cmp(a.doo, b.doo);
    cmp(a.dvv, b.dvv);
    cmp(a.dvo, b.dvo);
    cmp(a.dovov, b.dovov);
    cmp(a.dooov, b.dooov);
    cmp(a.dovvv, b.dovvv);
    return w;
}

double ccsd_density_bytes(int oi, int vi, bool triples) {
    const double o = oi, v = vi, n = o + v;
    double b = 8.0 * (12.0 * o * o * v * v + 6.0 * o * v * v * v + 3.0 * v * v * v + 2.0 * o * o * o * v);
    b += 8.0 * static_cast<double>(Packed4::size_for(static_cast<std::size_t>(n)));
    if (triples) b += 8.0 * (std::max(1.0, std::min(4.0, o * o)) * (o * v * v * v + 2.0 * o * o * v * v) + 20.0 * v * v * v +
                             2.0 * o * v * v * v);
    return b;
}

CCDensityResult ccsd_densities(const CCIntegrals& ints, LadderOperator& ladder, const std::vector<double>& t1v,
                               const std::vector<double>& t2v, bool triples, const CCOptions& opt,
                               const LambdaProgress& progress) {
    const auto t0 = Clock::now();
    set_tensor_threads(opt.threads);
    const Ints e = make_ints(ints);
    const std::size_t o = e.o, v = e.v;
    const Tensor t1 = from(t1v, {o, v}), t2 = from(t2v, {o, o, v, v});
    CCDensityResult res;
    const Imds m = make_intermediates(e, t1, t2);
    std::unique_ptr<TriplesTerms> tt;
    if (triples) {
        const auto t3 = Clock::now();
        report_progress(JobStage::Triples, 0);
        tt = std::make_unique<TriplesTerms>(triples_terms_blocked(e, t1, t2, resolve_threads(opt.threads)));
        res.triples_seconds = seconds_since(t3);
    }
    Tensor l1 = t1, l2 = t2;
    Diis diis(opt.diis_size);
    report_progress(JobStage::Lambda, 0);
    for (int it = 1; it <= opt.max_iter; ++it) {
        check_cancelled();
        Tensor n1, n2;
        update_lambda(e, ladder, t1, t2, l1, l2, m, tt.get(), n1, n2);
        double dmax = 0.0;
        std::vector<double> x(n1.size() + n2.size()), err(x.size());
        for (std::size_t k = 0; k < n1.size(); ++k) {
            x[k] = n1[k];
            err[k] = n1[k] - l1[k];
        }
        for (std::size_t k = 0; k < n2.size(); ++k) {
            x[n1.size() + k] = n2[k];
            err[n1.size() + k] = n2[k] - l2[k];
        }
        for (double d : err) dmax = std::max(dmax, std::fabs(d));
        res.lambda_iterations = it;
        res.lambda_residual = dmax;
        res.history.push_back({it, dmax});
        report_progress(JobStage::Lambda, it, 0.0, dmax);
        if (progress) progress(res.history.back());
        if (dmax < opt.amplitude_tol) {
            l1 = n1;
            l2 = n2;
            res.lambda_converged = true;
            break;
        }
        const std::vector<double> y = diis.extrapolate(x, err);
        std::copy(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(n1.size()), l1.data());
        std::copy(y.begin() + static_cast<std::ptrdiff_t>(n1.size()), y.end(), l2.data());
    }
    res.lambda_seconds = seconds_since(t0);
    const auto t1s = Clock::now();
    res.density.G = Packed4(o + v);
    build_densities(e, t1, t2, l1, l2, tt.get(), res.density);
    res.density_seconds = seconds_since(t1s);
    return res;
}

}  // namespace macqc
