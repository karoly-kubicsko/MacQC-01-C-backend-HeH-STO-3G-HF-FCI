#include "macqc/correlated_gradient.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

Matrix mul(const Matrix& A, bool ta, const Matrix& B, bool tb, int threads) {
    const std::size_t m = ta ? A.cols() : A.rows();
    const std::size_t k = ta ? A.rows() : A.cols();
    const std::size_t n = tb ? B.rows() : B.cols();
    Matrix C(m, n);
    if (m && n && k) gemm(ta, tb, m, n, k, 1.0, A.data(), A.cols(), B.data(), B.cols(), 0.0, C.data(), n, threads);
    return C;
}

Matrix columns(const Matrix& C, std::size_t c0, std::size_t nc) {
    Matrix out(C.rows(), nc);
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < nc; ++j) out(i, j) = C(i, c0 + j);
    return out;
}

Matrix symmetrized(const Matrix& A) {
    Matrix S(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) S(i, j) = 0.5 * (A(i, j) + A(j, i));
    return S;
}

// Packed MO two-particle density -> packed AO density: Gamma_mnls = sum C_mp C_nq C_lr C_ms G_pqrs.
std::vector<double> back_transform(const Packed4& G, const Matrix& C, int threads) {
    const std::size_t N = C.rows(), n = C.cols();
    const std::size_t np = n * (n + 1) / 2, NP = N * (N + 1) / 2;
    std::vector<double> half(NP * np, 0.0);  // [mn][rs]
    {
        Matrix M(n, n);
        for (std::size_t r = 0, rs = 0; r < n; ++r)
            for (std::size_t s = 0; s <= r; ++s, ++rs) {
                check_cancelled();
                for (std::size_t p = 0; p < n; ++p)
                    for (std::size_t q = 0; q <= p; ++q) M(p, q) = M(q, p) = G.get(p, q, r, s);
                const Matrix A = mul(C, false, mul(M, false, C, true, threads), false, threads);
                for (std::size_t m = 0, mn = 0; m < N; ++m)
                    for (std::size_t k = 0; k <= m; ++k, ++mn) half[mn * np + rs] = A(m, k);
            }
    }
    std::vector<double> out(NP * (NP + 1) / 2, 0.0);
    Matrix M(n, n);
    for (std::size_t m = 0, mn = 0; m < N; ++m)
        for (std::size_t k = 0; k <= m; ++k, ++mn) {
            check_cancelled();
            const double* h = half.data() + mn * np;
            for (std::size_t r = 0, rs = 0; r < n; ++r)
                for (std::size_t s = 0; s <= r; ++s, ++rs) M(r, s) = M(s, r) = h[rs];
            const Matrix A = mul(C, false, mul(M, false, C, true, threads), false, threads);
            for (std::size_t l = 0, ls = 0; l < N && ls <= mn; ++l)
                for (std::size_t s = 0; s <= l && ls <= mn; ++s, ++ls) out[mn * (mn + 1) / 2 + ls] = A(l, s);
        }
    return out;
}

}  // namespace

double density_energy(const CorrelatedDensity& d, const MOIntegrals& mo, double nuclear_repulsion) {
    const std::size_t n = static_cast<std::size_t>(mo.n);
    double e = nuclear_repulsion;
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) e += d.D(p, q) * mo.h(p, q);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q <= p; ++q)
            for (std::size_t r = 0; r <= p; ++r)
                for (std::size_t s = 0; s <= (r == p ? q : r); ++s)
                    e += Packed4::degeneracy(p, q, r, s) * d.G.get(p, q, r, s) * mo(p, q, r, s);
    return e;
}

CorrelatedDensity rhf_density(int n_mo, int n_occ) {
    CorrelatedDensity d;
    const std::size_t n = static_cast<std::size_t>(n_mo), o = static_cast<std::size_t>(n_occ);
    d.D = Matrix(n, n);
    for (std::size_t i = 0; i < o; ++i) d.D(i, i) = 2.0;
    d.G = Packed4(n);
    // Gamma_pqrs = 1/2 D_pq D_rs - 1/4 D_ps D_rq (spin-summed, closed shell): E_2 = sum_ij 2(ii|jj) - (ij|ji)
    for (std::size_t i = 0; i < o; ++i)
        for (std::size_t j = 0; j < o; ++j) {
            d.G.add_raw(i, i, j, j, 2.0);
            d.G.add_raw(i, j, j, i, -1.0);
        }
    return d;
}

double correlated_gradient_bytes(std::size_t n, std::size_t N) {
    const double np = static_cast<double>(n) * (n + 1) / 2.0, NP = static_cast<double>(N) * (N + 1) / 2.0;
    const double mo_ints = np * np, g_mo = np * (np + 1) / 2.0, half = NP * np, g_ao = NP * (NP + 1) / 2.0;
    return 8.0 * std::max({mo_ints + g_mo, g_mo + half, half + g_ao}) + 64.0 * static_cast<double>(N) * N * 8.0;
}

CorrelatedGradientResult correlated_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                             const MOIntegrals& mo, const CorrelatedDensity& dens, int threads,
                                             double memory_mb, const CPHFCallback& progress) {
    if (rhf.removed_functions > 0)
        throw std::runtime_error("analytic correlated gradient with removed near-linear dependencies");
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const int nat = static_cast<int>(mol.atoms.size());
    const Matrix& C = rhf.C;
    const std::size_t n = C.cols(), N = C.rows(), o = static_cast<std::size_t>(rhf.n_occupied), v = n - o;
    if (static_cast<std::size_t>(mo.n) != n || dens.D.rows() != n || dens.G.n != n)
        throw std::invalid_argument("correlated_gradient: dimensions of the densities and integrals differ");
    const std::vector<double>& eps = rhf.orbital_energies;
    CorrelatedGradientResult res;
    res.density_energy = density_energy(dens, mo, rhf.nuclear_repulsion);
    const Matrix Ds = symmetrized(dens.D);

    // ---- orbital Lagrangian X_tp = 2 sum_q h_tq D_qp + 4 sum_qrs (tq|rs) G_pqrs
    report_progress(JobStage::Gradient, 1);
    Matrix X = mul(mo.h, false, Ds, false, nt);
    for (std::size_t t = 0; t < n; ++t)
        for (std::size_t p = 0; p < n; ++p) X(t, p) *= 2.0;
    {
        const std::size_t np = n * (n + 1) / 2;
        std::vector<double> A(n * np), B(n * np), Xq(n * n);
        for (std::size_t q = 0; q < n; ++q) {
            check_cancelled();
            for (std::size_t t = 0; t < n; ++t)
                for (std::size_t r = 0, rs = 0; r < n; ++r)
                    for (std::size_t s = 0; s <= r; ++s, ++rs) {
                        const double w = r == s ? 1.0 : 2.0;
                        A[t * np + rs] = mo(t, q, r, s);
                        B[t * np + rs] = w * dens.G.get(t, q, r, s);
                    }
            gemm(false, true, n, n, np, 4.0, A.data(), np, B.data(), np, 0.0, Xq.data(), n, nt);
            for (std::size_t t = 0; t < n; ++t)
                for (std::size_t p = 0; p < n; ++p) X(t, p) += Xq[t * n + p];
        }
    }
    res.lagrangian_seconds = seconds_since(t0);
    // Invariance check: X - X' vanishes in the occupied-occupied and virtual-virtual blocks.
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < p; ++q)
            if ((p < o) == (q < o)) res.orbital_lagrangian_oo_vv = std::max(res.orbital_lagrangian_oo_vv, std::fabs(X(p, q) - X(q, p)));

    // ---- Z-vector
    const auto t1 = Clock::now();
    FockBuilder fock(basis, nt, memory_mb);
    Matrix Y(v, o);
    for (std::size_t a = 0; a < v; ++a)
        for (std::size_t i = 0; i < o; ++i) Y(a, i) = X(o + a, i) - X(i, o + a);
    const CPHFSolution zs = solve_cphf(fock, C, eps, static_cast<int>(o), {Y}, 1.0e-9, 80, progress, nt);
    res.zvector_converged = zs.converged;
    res.zvector_iterations = zs.iterations;
    res.zvector_residual = zs.residual;
    const Matrix& Z = zs.U[0];
    const Matrix Co = columns(C, 0, o), Cv = columns(C, o, v);
    const Matrix DZ = symmetrized(mul(mul(Cv, false, Z, false, nt), false, Co, true, nt));  // sym(Cv Z Co')
    const Matrix GDZ = fock.build(DZ);
    const Matrix GDZoo = transform_matrix(Co, GDZ, Co, nt);
    res.zvector_seconds = seconds_since(t1);

    // ---- energy-weighted density (MO) and the AO densities
    Matrix w(n, n);
    for (std::size_t t = 0; t < n; ++t)
        for (std::size_t p = 0; p < n; ++p) w(t, p) = -0.5 * X(t, p);
    for (std::size_t a = 0; a < v; ++a)
        for (std::size_t i = 0; i < o; ++i) w(o + a, i) += 0.5 * Y(a, i) + Z(a, i) * eps[i];
    for (std::size_t k = 0; k < o; ++k)
        for (std::size_t l = 0; l < o; ++l) w(k, l) += 2.0 * GDZoo(k, l);
    const Matrix Wao = symmetrized(mul(mul(C, false, w, false, nt), false, C, true, nt));
    Matrix D1 = mul(mul(C, false, Ds, false, nt), false, C, true, nt);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) D1(i, j) -= DZ(i, j);
    D1 = symmetrized(D1);
    res.relaxed_density_ao = D1;
    res.relaxed_dipole = dipole_moment(mol, basis, D1);
    Matrix Pref = mul(Co, false, Co, true, nt);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) Pref(i, j) *= 2.0;

    // ---- AO two-particle density
    const auto t2 = Clock::now();
    report_progress(JobStage::Gradient, 2);
    const std::vector<double> gao = back_transform(dens.G, C, nt);
    res.backtransform_seconds = seconds_since(t2);
    TwoParticleDensity gamma;
    gamma.packed = &gao;
    SeparableDensity zterm;
    zterm.A = DZ;
    zterm.B = Pref;
    zterm.coulomb = -1.0;
    zterm.exchange = 0.5;
    gamma.separable.push_back(zterm);

    // ---- derivative integrals
    const auto t3 = Clock::now();
    report_progress(JobStage::Gradient, 3);
    res.nuclear = nuclear_repulsion_gradient(mol);
    res.one_electron = core_hamiltonian_gradient(basis, mol, D1, nt);
    res.two_electron = eri_gradient(basis, nat, gamma, nt);
    res.overlap = overlap_gradient(basis, nat, Wao, nt);
    res.integral_seconds = seconds_since(t3);
    res.gradient.assign(static_cast<std::size_t>(3 * nat), 0.0);
    for (std::size_t x = 0; x < res.gradient.size(); ++x)
        res.gradient[x] = res.nuclear[x] + res.one_electron[x] + res.two_electron[x] + res.overlap[x];
    res.seconds = seconds_since(t0);
    return res;
}

}  // namespace macqc
