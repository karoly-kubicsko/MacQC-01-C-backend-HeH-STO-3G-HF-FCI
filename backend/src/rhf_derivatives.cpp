#include "macqc/rhf_derivatives.h"

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

// op(A) op(B) with the blocked gemm.
Matrix mul(const Matrix& A, bool ta, const Matrix& B, bool tb, int threads = 1) {
    const std::size_t m = ta ? A.cols() : A.rows();
    const std::size_t k = ta ? A.rows() : A.cols();
    const std::size_t n = tb ? B.rows() : B.cols();
    Matrix C(m, n);
    if (m == 0 || n == 0 || k == 0) return C;
    gemm(ta, tb, m, n, k, 1.0, A.data(), A.cols(), B.data(), B.cols(), 0.0, C.data(), n, threads);
    return C;
}

Matrix columns(const Matrix& C, std::size_t c0, std::size_t nc) {
    Matrix out(C.rows(), nc);
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < nc; ++j) out(i, j) = C(i, c0 + j);
    return out;
}

Matrix block(const Matrix& M, std::size_t r0, std::size_t nr, std::size_t c0, std::size_t nc) {
    Matrix out(nr, nc);
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) out(i, j) = M(r0 + i, c0 + j);
    return out;
}

void axpy(Matrix& y, double a, const Matrix& x) {
    for (std::size_t i = 0; i < y.rows(); ++i)
        for (std::size_t j = 0; j < y.cols(); ++j) y(i, j) += a * x(i, j);
}

double dot(const Matrix& a, const Matrix& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) s += a(i, j) * b(i, j);
    return s;
}

// Symmetric part of x y^T + y x^T scaled: returns f (X Y^T + Y X^T).
Matrix sym_outer(const Matrix& X, const Matrix& Y, double f, int threads) {
    Matrix A = mul(X, false, Y, true, threads);
    Matrix out(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) out(i, j) = f * (A(i, j) + A(j, i));
    return out;
}

// Dense solve by Gaussian elimination with partial pivoting; false if singular.
bool solve_dense(std::vector<double> a, std::vector<double> b, std::size_t n, std::vector<double>& x) {
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (std::fabs(a[r * n + col]) > std::fabs(a[piv * n + col])) piv = r;
        if (std::fabs(a[piv * n + col]) < 1.0e-300) return false;
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

Matrix occupied_density(const Matrix& C, int nocc, int threads) {
    const Matrix Co = columns(C, 0, static_cast<std::size_t>(nocc));
    Matrix P = mul(Co, false, Co, true, threads);
    for (std::size_t i = 0; i < P.rows(); ++i)
        for (std::size_t j = 0; j < P.cols(); ++j) P(i, j) *= 2.0;
    return P;
}

Matrix energy_weighted_density(const Matrix& C, const std::vector<double>& eps, int nocc, int threads) {
    Matrix Co = columns(C, 0, static_cast<std::size_t>(nocc));
    Matrix Ce = Co;
    for (std::size_t i = 0; i < Ce.rows(); ++i)
        for (int k = 0; k < nocc; ++k) Ce(i, static_cast<std::size_t>(k)) *= 2.0 * eps[static_cast<std::size_t>(k)];
    Matrix W = mul(Ce, false, Co, true, threads);
    for (std::size_t i = 0; i < W.rows(); ++i)
        for (std::size_t j = 0; j < i; ++j) W(i, j) = W(j, i) = 0.5 * (W(i, j) + W(j, i));
    return W;
}

TwoParticleDensity rhf_gamma(const Matrix& P) {
    TwoParticleDensity g;
    SeparableDensity t;
    t.A = P;
    t.B = P;
    t.coulomb = 0.5;
    t.exchange = -0.25;
    g.separable.push_back(t);
    return g;
}

}  // namespace

Matrix transform_matrix(const Matrix& left, const Matrix& M, const Matrix& right, int threads) {
    return mul(mul(left, true, M, false, threads), false, right, false, threads);
}

RHFGradientResult rhf_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf, int threads,
                               double lindep_tol) {
    const auto t0 = Clock::now();
    const int nat = static_cast<int>(mol.atoms.size());
    const int nt = resolve_threads(threads);
    const int nocc = rhf.n_occupied;
    const Matrix P = occupied_density(rhf.C, nocc, nt);
    const Matrix W = energy_weighted_density(rhf.C, rhf.orbital_energies, nocc, nt);

    RHFGradientResult r;
    r.nuclear = nuclear_repulsion_gradient(mol);
    report_progress(JobStage::Gradient, 1);
    r.one_electron = core_hamiltonian_gradient(basis, mol, P, nt);
    report_progress(JobStage::Gradient, 2);
    r.two_electron = eri_gradient(basis, nat, rhf_gamma(P), nt);
    r.overlap = overlap_gradient(basis, nat, W, nt);
    for (double& v : r.overlap) v = -v;
    r.lindep.assign(static_cast<std::size_t>(3 * nat), 0.0);

    if (rhf.removed_functions > 0) {
        // Rotation of the kept orbital space with the geometry (see the header).
        const EigenResult se = symmetric_eigen(rhf.S);
        const std::size_t N = rhf.S.rows();
        std::vector<std::size_t> kept, removed;
        for (std::size_t k = 0; k < N; ++k) (se.values[k] >= lindep_tol ? kept : removed).push_back(k);
        FockBuilder fock(basis, nt, 0.0);
        Matrix F = fock.build(P);
        axpy(F, 1.0, rhf.H);
        const Matrix FPS = mul(mul(F, false, P, false, nt), false, rhf.S, false, nt);
        const Matrix E = se.vectors;  // columns
        const Matrix X = mul(E, true, mul(FPS, false, E, false, nt), false, nt);  // e_l' F P S e_k
        Matrix Y(N, N);
        for (std::size_t l : removed)
            for (std::size_t k : kept) {
                const double c = 2.0 * X(l, k) / (se.values[k] * (se.values[k] - se.values[l]));
                for (std::size_t m = 0; m < N; ++m)
                    for (std::size_t n = 0; n < N; ++n) Y(m, n) += c * E(m, l) * E(n, k);
            }
        Matrix Ys(N, N);
        for (std::size_t m = 0; m < N; ++m)
            for (std::size_t n = 0; n < N; ++n) Ys(m, n) = 0.5 * (Y(m, n) + Y(n, m));
        r.lindep = overlap_gradient(basis, nat, Ys, nt);
    }
    r.gradient.assign(static_cast<std::size_t>(3 * nat), 0.0);
    for (std::size_t x = 0; x < r.gradient.size(); ++x)
        r.gradient[x] = r.nuclear[x] + r.one_electron[x] + r.two_electron[x] + r.overlap[x] + r.lindep[x];
    r.seconds = seconds_since(t0);
    return r;
}

std::array<double, 3> dipole_moment(const Molecule& mol, const BasisSet& basis, const Matrix& P,
                                    const std::array<double, 3>& origin) {
    const std::vector<Matrix> r = dipole_matrices(basis, origin);
    std::array<double, 3> mu{};
    for (int k = 0; k < 3; ++k) {
        mu[static_cast<std::size_t>(k)] = -dot(P, r[static_cast<std::size_t>(k)]);
        for (const Atom& a : mol.atoms)
            mu[static_cast<std::size_t>(k)] += a.Z * (a.position[static_cast<std::size_t>(k)] - origin[static_cast<std::size_t>(k)]);
    }
    return mu;
}

CPHFSolution solve_cphf(const FockBuilder& fock, const Matrix& C, const std::vector<double>& eps, int nocc,
                        const std::vector<Matrix>& rhs, double tol, int max_iter, const CPHFCallback& progress,
                        int threads) {
    const int nt = resolve_threads(threads);
    const std::size_t o = static_cast<std::size_t>(nocc);
    const std::size_t nmo = C.cols(), v = nmo - o;
    const Matrix Co = columns(C, 0, o), Cv = columns(C, o, v);
    Matrix Dg(v, o);
    for (std::size_t a = 0; a < v; ++a)
        for (std::size_t i = 0; i < o; ++i) Dg(a, i) = eps[o + a] - eps[i];
    auto apply = [&](const Matrix& U) {  // A U
        const Matrix CvU = mul(Cv, false, U, false, nt);  // N x o
        const Matrix D = sym_outer(CvU, Co, 2.0, nt);   // 2 (Cv U Co' + Co U' Cv')
        const Matrix G = fock.build(D);
        Matrix AU = transform_matrix(Cv, G, Co, nt);
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t i = 0; i < o; ++i) AU(a, i) += Dg(a, i) * U(a, i);
        return AU;
    };
    auto precondition = [&](const Matrix& R) {
        Matrix Z(v, o);
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t i = 0; i < o; ++i) Z(a, i) = R(a, i) / Dg(a, i);
        return Z;
    };

    const std::size_t nx = rhs.size();
    CPHFSolution sol;
    sol.U.assign(nx, Matrix(v, o));
    std::vector<std::vector<Matrix>> V(nx), AV(nx);
    std::vector<Matrix> next(nx);
    std::vector<bool> done(nx, false);
    for (std::size_t x = 0; x < nx; ++x) {
        const double bn = std::sqrt(dot(rhs[x], rhs[x]));
        if (bn < 1.0e-14) {
            done[x] = true;
            continue;
        }
        next[x] = precondition(rhs[x]);
    }
    double worst = 0.0;
    for (int it = 1; it <= max_iter; ++it) {
        check_cancelled();
        bool any = false;
        for (std::size_t x = 0; x < nx; ++x) {
            if (done[x]) continue;
            // Orthonormalise the new vector against the subspace (twice for stability).
            Matrix w = next[x];
            for (int pass = 0; pass < 2; ++pass)
                for (const Matrix& q : V[x]) axpy(w, -dot(q, w), q);
            const double nw = std::sqrt(dot(w, w));
            if (nw < 1.0e-14) {
                done[x] = true;
                continue;
            }
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t i = 0; i < o; ++i) w(a, i) /= nw;
            V[x].push_back(w);
            AV[x].push_back(apply(w));
            any = true;
        }
        if (!any) break;
        worst = 0.0;
        for (std::size_t x = 0; x < nx; ++x) {
            if (V[x].empty()) continue;
            const std::size_t m = V[x].size();
            std::vector<double> M(m * m), b(m), c;
            for (std::size_t i = 0; i < m; ++i) {
                b[i] = dot(V[x][i], rhs[x]);
                for (std::size_t j = 0; j < m; ++j) M[i * m + j] = 0.5 * (dot(V[x][i], AV[x][j]) + dot(V[x][j], AV[x][i]));
            }
            if (!solve_dense(M, b, m, c)) throw std::runtime_error("CPHF: singular subspace matrix");
            Matrix U(v, o), R = rhs[x];
            for (std::size_t i = 0; i < m; ++i) {
                axpy(U, c[i], V[x][i]);
                axpy(R, -c[i], AV[x][i]);
            }
            sol.U[x] = U;
            const double rn = std::sqrt(dot(R, R));
            worst = std::max(worst, rn);
            if (rn < tol) {
                done[x] = true;
            } else {
                next[x] = precondition(R);
            }
        }
        sol.iterations = it;
        sol.residual = worst;
        report_progress(JobStage::CPHF, it, 0.0, worst);
        if (progress) progress(CPHFProgress{it, worst});
        if (std::all_of(done.begin(), done.end(), [](bool d) { return d; })) break;
    }
    sol.converged = worst < tol;
    return sol;
}

RHFHessianResult rhf_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf, int threads,
                             double memory_mb, double cphf_tol, const CPHFCallback& progress) {
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const int nat = static_cast<int>(mol.atoms.size());
    const std::size_t n3 = static_cast<std::size_t>(3 * nat);
    const int nocc = rhf.n_occupied;
    const std::size_t o = static_cast<std::size_t>(nocc);
    const Matrix& C = rhf.C;
    const std::size_t nmo = C.cols(), v = nmo - o;
    const std::vector<double>& eps = rhf.orbital_energies;
    const Matrix P = occupied_density(C, nocc, nt);
    const Matrix W = energy_weighted_density(C, eps, nocc, nt);
    const Matrix Co = columns(C, 0, o), Cv = columns(C, o, v);

    RHFHessianResult res;
    // ---- derivative integrals
    report_progress(JobStage::Hessian, 1);
    const std::vector<Matrix> Sx = overlap_derivative_matrices(basis, nat, nt);
    std::vector<Matrix> Fx = core_hamiltonian_derivative_matrices(basis, mol, nt);
    {
        const std::vector<Matrix> Gx = eri_derivative_fock(basis, nat, P, nt);
        for (std::size_t x = 0; x < n3; ++x) axpy(Fx[x], 1.0, Gx[x]);
    }
    report_progress(JobStage::Hessian, 2);
    Matrix H = core_hamiltonian_hessian(basis, mol, P, nt);
    axpy(H, 1.0, eri_hessian(basis, nat, rhf_gamma(P), nt));
    axpy(H, -1.0, overlap_hessian(basis, nat, W, nt));
    axpy(H, 1.0, nuclear_repulsion_hessian(mol));
    res.integral_seconds = seconds_since(t0);

    // ---- MO quantities and CPHF right-hand sides
    const auto t1 = Clock::now();
    FockBuilder fock(basis, nt, memory_mb);
    std::vector<Matrix> Smo(n3), Fmo(n3), GS(n3), rhs(n3);
    for (std::size_t x = 0; x < n3; ++x) {
        check_cancelled();
        Smo[x] = transform_matrix(C, Sx[x], C, nt);
        Fmo[x] = transform_matrix(C, Fx[x], C, nt);
        const Matrix Soo = block(Smo[x], 0, o, 0, o);
        const Matrix DS = mul(mul(Co, false, Soo, false, nt), false, Co, true, nt);  // Co S_oo Co'
        GS[x] = fock.build(DS);
        const Matrix GSvo = transform_matrix(Cv, GS[x], Co, nt);
        Matrix B(v, o);
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t i = 0; i < o; ++i)
                B(a, i) = -Fmo[x](o + a, i) + eps[i] * Smo[x](o + a, i) + 2.0 * GSvo(a, i);
        rhs[x] = B;
    }
    report_progress(JobStage::CPHF, 0);
    const CPHFSolution sol = solve_cphf(fock, C, eps, nocc, rhs, cphf_tol, 60, progress, nt);
    res.cphf_converged = sol.converged;
    res.cphf_iterations = sol.iterations;
    res.cphf_residual = sol.residual;

    // ---- full response U^y (nmo x o) and the derivative of the occupied Fock block
    std::vector<Matrix> U(n3), Foo(n3);
    for (std::size_t y = 0; y < n3; ++y) {
        check_cancelled();
        Matrix Uy(nmo, o);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j) Uy(i, j) = -0.5 * Smo[y](i, j);
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t i = 0; i < o; ++i) Uy(o + a, i) = sol.U[y](a, i);
        U[y] = Uy;
        // P^y = 2 (Cv U Co' + Co U' Cv') - 2 Co S_oo Co'
        const Matrix CvU = mul(Cv, false, sol.U[y], false, nt);
        Matrix Py = sym_outer(CvU, Co, 2.0, nt);
        const Matrix Soo = block(Smo[y], 0, o, 0, o);
        axpy(Py, -2.0, mul(mul(Co, false, Soo, false, nt), false, Co, true, nt));
        const Matrix Gy = fock.build(Py);
        const Matrix Goo = transform_matrix(Co, Gy, Co, nt);
        Matrix F(o, o);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                F(i, j) = -0.5 * Smo[y](i, j) * (eps[i] + eps[j]) + Fmo[y](i, j) + Goo(i, j);
        Foo[y] = F;
    }
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t y = 0; y < n3; ++y) {
            double s = 0.0;
            for (std::size_t p = 0; p < nmo; ++p)
                for (std::size_t i = 0; i < o; ++i) s += U[y](p, i) * (Fmo[x](p, i) - eps[i] * Smo[x](p, i));
            double t = 0.0;
            for (std::size_t i = 0; i < o; ++i)
                for (std::size_t j = 0; j < o; ++j) t += Foo[y](i, j) * Smo[x](i, j);
            H(x, y) += 4.0 * s - 2.0 * t;
        }
    // Symmetrise (the two halves agree to the CPHF convergence).
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t y = 0; y < x; ++y) H(x, y) = H(y, x) = 0.5 * (H(x, y) + H(y, x));
    res.hessian = H;

    // ---- dipole moment and its derivatives (atomic polar tensors)
    const std::array<double, 3> origin{0.0, 0.0, 0.0};
    const std::vector<Matrix> r = dipole_matrices(basis, origin);
    res.dipole = dipole_moment(mol, basis, P, origin);
    const std::vector<std::vector<Matrix>> rx = dipole_derivative_matrices(basis, nat, origin, nt);
    std::vector<Matrix> rmo(3);
    for (int m = 0; m < 3; ++m) rmo[static_cast<std::size_t>(m)] = transform_matrix(C, r[static_cast<std::size_t>(m)], C, nt);
    res.dipole_derivatives.assign(n3, {0.0, 0.0, 0.0});
    for (std::size_t x = 0; x < n3; ++x)
        for (int m = 0; m < 3; ++m) {
            const Matrix& R = rmo[static_cast<std::size_t>(m)];
            double resp = 0.0;  // sum P^x r^m = 4 sum_ai U_ai r_ai - 2 sum_kl S_kl r_kl
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t i = 0; i < o; ++i) resp += 4.0 * U[x](o + a, i) * R(o + a, i);
            for (std::size_t k = 0; k < o; ++k)
                for (std::size_t l = 0; l < o; ++l) resp -= 2.0 * Smo[x](k, l) * R(k, l);
            double val = -resp - dot(P, rx[x][static_cast<std::size_t>(m)]);
            if (static_cast<int>(x % 3) == m) val += mol.atoms[x / 3].Z;
            res.dipole_derivatives[x][static_cast<std::size_t>(m)] = val;
        }
    res.cphf_seconds = seconds_since(t1);
    res.seconds = seconds_since(t0);
    return res;
}

}  // namespace macqc
