#include "macqc/hessian.h"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "macqc/constants.h"
#include "macqc/gto_derivatives.h"
#include "macqc/internal_coordinates.h"
#include "macqc/job.h"
#include "macqc/rhf_derivatives.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::string format(const char* fmt, ...) {
    char buf[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    return buf;
}

struct Direction {
    std::string name;
    std::vector<double> d;  // in (R1, R2, angle) or (R)
    bool antisymmetric = false;  // C2v B2 coordinate
};

// Displacement directions: symmetry-adapted for C2v H2O.
std::vector<Direction> directions(const Molecule& mol, const SimpleInternals& si) {
    std::vector<Direction> out;
    if (si.values.size() == 1) {
        out.push_back({si.names[0], {1.0}, false});
        return out;
    }
    const bool c2v = std::fabs(si.values[0] - si.values[1]) < 1.0e-6 && mol.atoms[1].Z == mol.atoms[2].Z;
    if (c2v) {
        const double s = 1.0 / std::sqrt(2.0);
        out.push_back({"S1 (A1) = [R(O1-H2) + R(O1-H3)]/sqrt(2), symmetric stretch", {s, s, 0.0}, false});
        out.push_back({"S2 (A1) = A(H2-O1-H3), bend", {0.0, 0.0, 1.0}, false});
        out.push_back({"S3 (B2) = [R(O1-H2) - R(O1-H3)]/sqrt(2), antisymmetric stretch", {s, -s, 0.0}, true});
    } else {
        for (std::size_t k = 0; k < 3; ++k) {
            std::vector<double> d(3, 0.0);
            d[k] = 1.0;
            out.push_back({si.names[k], d, false});
        }
    }
    return out;
}

bool is_c2v(const std::vector<Direction>& dirs) { return dirs.size() == 3 && dirs[2].antisymmetric; }

// sum_i c_i d2q_i/dx dx' of the simple internal coordinates (central differences of B).
Matrix internal_second_derivatives(const Molecule& mol, const std::vector<double>& c) {
    const std::size_t n3 = 3 * mol.atoms.size();
    Matrix K(n3, n3);
    const double h = 1.0e-5;
    for (std::size_t j = 0; j < n3; ++j) {
        Molecule mp = mol, mm = mol;
        mp.atoms[j / 3].position[j % 3] += h;
        mm.atoms[j / 3].position[j % 3] -= h;
        const SimpleInternals p = simple_internals(mp), m = simple_internals(mm);
        for (std::size_t k = 0; k < n3; ++k) {
            double s = 0.0;
            for (std::size_t i = 0; i < c.size(); ++i) s += c[i] * (p.B[i][k] - m.B[i][k]) / (2.0 * h);
            K(j, k) = s;
        }
    }
    for (std::size_t j = 0; j < n3; ++j)
        for (std::size_t k = 0; k < j; ++k) K(j, k) = K(k, j) = 0.5 * (K(j, k) + K(k, j));
    return K;
}

// Internal gradient g_q = (B B')^-1 B g_x of the (non-redundant) simple internals.
std::vector<double> simple_internal_gradient(const Molecule& mol, const std::vector<double>& gx) {
    const SimpleInternals si = simple_internals(mol);
    const std::size_t nq = si.values.size(), n3 = gx.size();
    Matrix B(nq, n3), G(nq, nq);
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t x = 0; x < n3; ++x) B(i, x) = si.B[i][x];
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t j = 0; j < nq; ++j)
            for (std::size_t x = 0; x < n3; ++x) G(i, j) += B(i, x) * B(j, x);
    return internal_gradient(B, generalized_inverse(G), gx);
}

// Cartesian Hessian from the Hessian in the displacement directions and the internal gradient.
Matrix cartesian_hessian(const Molecule& mol, const SimpleInternals& si, const std::vector<Direction>& dirs,
                         const Matrix& HS, const std::vector<double>& gq) {
    const std::size_t nq = si.values.size(), ns = dirs.size(), n3 = 3 * mol.atoms.size();
    Matrix Hq(nq, nq);
    for (std::size_t a = 0; a < nq; ++a)
        for (std::size_t b = 0; b < nq; ++b)
            for (std::size_t k = 0; k < ns; ++k)
                for (std::size_t l = 0; l < ns; ++l) Hq(a, b) += dirs[k].d[a] * HS(k, l) * dirs[l].d[b];
    Matrix H = internal_second_derivatives(mol, gq);
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t y = 0; y < n3; ++y)
            for (std::size_t a = 0; a < nq; ++a)
                for (std::size_t b = 0; b < nq; ++b) H(x, y) += si.B[a][x] * Hq(a, b) * si.B[b][y];
    return H;
}

// Atomic polar tensor B' D dmu/dS from dipole derivatives along the directions.
std::vector<std::array<double, 3>> apt_from_directions(const SimpleInternals& si, const std::vector<Direction>& dirs,
                                                       const std::vector<std::array<double, 3>>& dmu) {
    const std::size_t n3 = si.B[0].size();
    std::vector<std::array<double, 3>> P(n3, {0.0, 0.0, 0.0});
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t k = 0; k < dirs.size(); ++k) {
            double bd = 0.0;
            for (std::size_t a = 0; a < si.values.size(); ++a) bd += si.B[a][x] * dirs[k].d[a];
            for (std::size_t c = 0; c < 3; ++c) P[x][c] += bd * dmu[k][c];
        }
    return P;
}

Molecule displaced(const Molecule& mol, const SimpleInternals& si, const Direction& d, double h) {
    std::vector<double> q = si.values;
    for (std::size_t a = 0; a < q.size(); ++a) q[a] += h * d.d[a];
    return align_to(mol, geometry_from_internals(mol, q));
}

}  // namespace

Molecule align_to(const Molecule& ref, const Molecule& mv) {
    const std::size_t n = ref.atoms.size();
    std::vector<double> w(n);
    double W = 0.0;
    std::array<double, 3> ca{}, cb{};
    for (std::size_t i = 0; i < n; ++i) {
        w[i] = isotope_mass(ref.atoms[i].Z);
        W += w[i];
        for (std::size_t k = 0; k < 3; ++k) {
            ca[k] += w[i] * ref.atoms[i].position[k];
            cb[k] += w[i] * mv.atoms[i].position[k];
        }
    }
    for (std::size_t k = 0; k < 3; ++k) {
        ca[k] /= W;
        cb[k] /= W;
    }
    double S[3][3] = {};
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b)
                S[a][b] += w[i] * (mv.atoms[i].position[a] - cb[a]) * (ref.atoms[i].position[b] - ca[b]);
    Matrix N(4, 4);
    N(0, 0) = S[0][0] + S[1][1] + S[2][2];
    N(0, 1) = N(1, 0) = S[1][2] - S[2][1];
    N(0, 2) = N(2, 0) = S[2][0] - S[0][2];
    N(0, 3) = N(3, 0) = S[0][1] - S[1][0];
    N(1, 1) = S[0][0] - S[1][1] - S[2][2];
    N(1, 2) = N(2, 1) = S[0][1] + S[1][0];
    N(1, 3) = N(3, 1) = S[2][0] + S[0][2];
    N(2, 2) = -S[0][0] + S[1][1] - S[2][2];
    N(2, 3) = N(3, 2) = S[1][2] + S[2][1];
    N(3, 3) = -S[0][0] - S[1][1] + S[2][2];
    const EigenResult e = symmetric_eigen(N);
    const double q0 = e.vectors(0, 3), q1 = e.vectors(1, 3), q2 = e.vectors(2, 3), q3 = e.vectors(3, 3);
    const double R[3][3] = {{q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)},
                            {2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)},
                            {2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3}};
    Molecule out = mv;
    for (std::size_t i = 0; i < n; ++i) {
        double p[3];
        for (std::size_t k = 0; k < 3; ++k) p[k] = mv.atoms[i].position[k] - cb[k];
        for (std::size_t a = 0; a < 3; ++a)
            out.atoms[i].position[a] = ca[a] + R[a][0] * p[0] + R[a][1] * p[1] + R[a][2] * p[2];
    }
    return out;
}

NumericalHessianCost numerical_hessian_cost(const Molecule& mol, bool analytic_gradients, bool dipoles) {
    NumericalHessianCost c;
    const SimpleInternals si = simple_internals(mol);
    const auto dirs = directions(mol, si);
    for (const auto& d : dirs) c.coordinates.push_back(d.name);
    const std::size_t n = dirs.size();
    if (analytic_gradients) {
        c.gradients = static_cast<int>(2 * n);
        return c;
    }
    if (n == 1) {
        c.energies = 2;
        const bool polar = mol.atoms[0].Z != mol.atoms[1].Z;
        if (dipoles && polar) c.field_energies = 4;
    } else if (is_c2v(dirs)) {
        c.energies = 4 + 4 + 1;  // S1 +-, S2 +-, (S1, S2) mixed 4, S3 + (S3 - by symmetry)
        if (dipoles) c.field_energies = 4 + 4 + 2;
    } else {
        c.energies = 6 + 12;
        if (dipoles) c.field_energies = 6 * 2 * 2;
    }
    return c;
}

HessianResult run_hessian(const Molecule& mol, const BasisSet& basis, const std::string& method,
                          const ElectronicOptions& eo, const HessianOptions& ho, const ElectronicCallbacks& cb,
                          const Matrix* guess, const DisplacementProgress& progress) {
    const auto t0 = Clock::now();
    HessianResult r;
    const std::size_t n3 = 3 * mol.atoms.size();
    r.reference = ho.reference ? *ho.reference : run_gradient(mol, basis, method, eo, guess, cb);
    const GradientResult& ref = r.reference;
    if (!ref.point.converged) throw std::runtime_error("the calculation at the reference geometry did not converge");
    r.dipole = ref.dipole;
    r.has_dipole = ref.has_dipole;
    const Matrix C0 = ref.point.rhf.C;

    // ---- analytic RHF
    if (method == "RHF" && ref.point.rhf.removed_functions == 0 && !ho.force_numerical && !ho.from_energies) {
        report_progress(JobStage::Hessian, 1);
        const RHFHessianResult h = rhf_hessian(mol, basis, ref.point.rhf, eo.threads, eo.memory_mb, 1.0e-9, cb.cphf);
        r.hessian = h.hessian;
        r.analytic = true;
        r.description = "analytic (second-derivative integrals, coupled-perturbed Hartree-Fock for all 3N coordinates)";
        r.has_dipole_derivatives = true;
        r.dipole_derivatives = h.dipole_derivatives;
        r.dipole_description = "analytic (CPHF response of the SCF density)";
        r.cphf_iterations = h.cphf_iterations;
        r.cphf_converged = h.cphf_converged;
        r.seconds = seconds_since(t0);
        return r;
    }
    if (method == "RHF" && ref.point.rhf.removed_functions > 0)
        r.numerical_reason = format("%d near-linear dependencies were removed; the analytic RHF Hessian does not contain "
                                    "the second-order rotation of the kept orbital space",
                                    ref.point.rhf.removed_functions);
    else if (ho.force_numerical || ho.from_energies)
        r.numerical_reason = "requested";
    else
        r.numerical_reason = "no analytic Hessian for " + method + " in this version";

    const SimpleInternals si = simple_internals(mol);
    const bool linear3 = si.values.size() == 3 && std::fabs(std::sin(si.values[2])) < 1.0e-6;
    const double h = ho.step;
    ElectronicOptions eod = eo;

    // ---- finite differences of analytic gradients
    if (ref.analytic && !ho.from_energies) {
        if (linear3) {
            // Cartesian displacements (the bend of a linear triatomic is degenerate)
            std::vector<std::vector<double>> gp(n3), gm(n3);
            std::vector<std::array<double, 3>> dp(n3), dm(n3);
            int idx = 0;
            for (std::size_t x = 0; x < n3; ++x) {
                for (int s : {1, -1}) {
                    Molecule m = mol;
                    m.atoms[x / 3].position[x % 3] += s * h;
                    ++idx;
                    if (progress) progress(idx, static_cast<int>(2 * n3), format("Cartesian coordinate %zu %s", x + 1, s > 0 ? "+" : "-"));
                    report_progress(JobStage::NumericalHessian, idx, ref.point.energy, 2.0 * n3);
                    const GradientResult g = run_gradient(m, basis_at_geometry(basis, m), method, eod, &C0);
                    ++r.gradient_evaluations;
                    (s > 0 ? gp[x] : gm[x]) = g.gradient;
                    (s > 0 ? dp[x] : dm[x]) = g.dipole;
                }
            }
            r.hessian = Matrix(n3, n3);
            for (std::size_t x = 0; x < n3; ++x)
                for (std::size_t y = 0; y < n3; ++y) r.hessian(x, y) = (gp[x][y] - gm[x][y]) / (2.0 * h);
            for (std::size_t x = 0; x < n3; ++x)
                for (std::size_t y = 0; y < x; ++y) r.hessian(x, y) = r.hessian(y, x) = 0.5 * (r.hessian(x, y) + r.hessian(y, x));
            r.dipole_derivatives.assign(n3, {0.0, 0.0, 0.0});
            for (std::size_t x = 0; x < n3; ++x)
                for (std::size_t c = 0; c < 3; ++c) r.dipole_derivatives[x][c] = (dp[x][c] - dm[x][c]) / (2.0 * h);
            r.has_dipole_derivatives = ref.has_dipole;
            r.description = format("numerical (central differences of analytic gradients, %d Cartesian displacements of "
                                   "%.4f bohr; the molecule is linear)",
                                   r.gradient_evaluations, h);
            r.dipole_description = "central differences of the relaxed dipole moments of the gradients";
            r.seconds = seconds_since(t0);
            return r;
        }
        const auto dirs = directions(mol, si);
        const std::size_t ns = dirs.size(), nq = si.values.size();
        for (const auto& d : dirs) r.coordinates.push_back(d.name);
        std::vector<std::vector<double>> dg(ns, std::vector<double>(nq, 0.0));
        std::vector<std::array<double, 3>> dmu(ns, {0.0, 0.0, 0.0});
        int idx = 0;
        const int total = static_cast<int>(2 * ns);
        for (std::size_t k = 0; k < ns; ++k) {
            std::vector<double> gq[2];
            std::array<double, 3> mu[2];
            for (int s = 0; s < 2; ++s) {
                const Molecule m = displaced(mol, si, dirs[k], s == 0 ? h : -h);
                ++idx;
                if (progress) progress(idx, total, dirs[k].name + (s == 0 ? " +" : " -"));
                report_progress(JobStage::NumericalHessian, idx, ref.point.energy, total);
                const GradientResult g = run_gradient(m, basis_at_geometry(basis, m), method, eod, &C0);
                if (!g.point.converged) throw std::runtime_error("a displaced calculation of the numerical Hessian did not converge");
                ++r.gradient_evaluations;
                gq[s] = simple_internal_gradient(m, g.gradient);
                mu[s] = g.dipole;
            }
            for (std::size_t a = 0; a < nq; ++a) dg[k][a] = (gq[0][a] - gq[1][a]) / (2.0 * h);
            for (std::size_t c = 0; c < 3; ++c) dmu[k][c] = (mu[0][c] - mu[1][c]) / (2.0 * h);
        }
        Matrix HS(ns, ns);
        for (std::size_t k = 0; k < ns; ++k)
            for (std::size_t l = 0; l < ns; ++l) {
                double s = 0.0;
                for (std::size_t a = 0; a < nq; ++a) s += dirs[l].d[a] * dg[k][a];
                HS(l, k) = s;
            }
        for (std::size_t k = 0; k < ns; ++k)
            for (std::size_t l = 0; l < k; ++l) HS(k, l) = HS(l, k) = 0.5 * (HS(k, l) + HS(l, k));
        const std::vector<double> gq0 = simple_internal_gradient(mol, ref.gradient);
        r.hessian = cartesian_hessian(mol, si, dirs, HS, gq0);
        r.has_dipole_derivatives = ref.has_dipole && ho.dipoles;
        if (r.has_dipole_derivatives) {
            r.dipole_derivatives = apt_from_directions(si, dirs, dmu);
            r.dipole_description = method == "RHF" ? "central differences of the SCF dipole moments of the gradients"
                                                   : "central differences of the relaxed dipole moments of the gradients";
        }
        r.description = format("numerical (central differences of analytic gradients along %zu %sinternal coordinate%s, "
                               "%d gradients, step %.4f bohr / rad)",
                               ns, is_c2v(dirs) ? "symmetry-adapted " : "", ns == 1 ? "" : "s", r.gradient_evaluations, h);
        r.seconds = seconds_since(t0);
        return r;
    }

    // ---- finite differences of energies (CCSDT, CCSDTQ)
    if (linear3)
        throw std::runtime_error("a numerical Hessian from energies at a linear H2O geometry is not available; "
                                 "optimise first or use a method with analytic gradients");
    const auto dirs = directions(mol, si);
    const std::size_t ns = dirs.size(), nq = si.values.size();
    for (const auto& d : dirs) r.coordinates.push_back(d.name);
    const bool c2v = is_c2v(dirs);
    const double E0 = ref.point.energy;
    const NumericalHessianCost cost = numerical_hessian_cost(mol, false, ho.dipoles);
    const int total = cost.energies + cost.field_energies;
    int idx = 0, field_energies = 0;
    auto energy = [&](const Molecule& m, const std::string& what, const std::array<double, 3>& field) {
        ++idx;
        if (progress) progress(idx, total, what);
        report_progress(JobStage::NumericalHessian, idx, E0, total);
        ElectronicOptions e2 = eod;
        e2.field = field;
        const SinglePoint p = run_single_point(m, basis_at_geometry(basis, m), method, e2, &C0);
        if (!p.converged) throw std::runtime_error("a displaced calculation of the numerical Hessian did not converge");
        ++r.energy_evaluations;
        if (field[0] != 0.0 || field[1] != 0.0 || field[2] != 0.0) ++field_energies;
        return p.energy;
    };
    // molecular axes for the finite-field dipoles: e_z along the bond (diatomic) or the C2 axis,
    // e_x in the plane
    std::array<double, 3> ez{}, ex{};
    {
        double b1[3], b2[3] = {0, 0, 0};
        for (std::size_t k = 0; k < 3; ++k) b1[k] = mol.atoms[1].position[k] - mol.atoms[0].position[k];
        const double n1 = std::sqrt(b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2]);
        for (double& v : b1) v /= n1;
        if (mol.atoms.size() == 2) {
            ez = {b1[0], b1[1], b1[2]};
        } else {
            for (std::size_t k = 0; k < 3; ++k) b2[k] = mol.atoms[2].position[k] - mol.atoms[0].position[k];
            const double n2 = std::sqrt(b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]);
            for (double& v : b2) v /= n2;
            double s[3] = {b1[0] + b2[0], b1[1] + b2[1], b1[2] + b2[2]}, d[3] = {b1[0] - b2[0], b1[1] - b2[1], b1[2] - b2[2]};
            const double ns_ = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]),
                         nd_ = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            ez = {s[0] / ns_, s[1] / ns_, s[2] / ns_};
            ex = {d[0] / nd_, d[1] / nd_, d[2] / nd_};
        }
    }
    const double F = ho.field;
    auto dipole_along = [&](const Molecule& m, const std::array<double, 3>& e, const std::string& what) {
        const std::array<double, 3> fp{F * e[0], F * e[1], F * e[2]}, fm{-F * e[0], -F * e[1], -F * e[2]};
        const double ep = energy(m, what + ", field +", fp), em = energy(m, what + ", field -", fm);
        return -(ep - em) / (2.0 * F);
    };
    const bool polar = mol.atoms.size() == 3 || mol.atoms[0].Z != mol.atoms[1].Z;
    const bool want_mu = ho.dipoles && polar;
    Matrix HS(ns, ns);
    std::vector<double> gS(ns, 0.0);
    std::vector<std::array<double, 3>> dmu(ns, {0.0, 0.0, 0.0});
    std::vector<double> Ep(ns), Em(ns);
    for (std::size_t k = 0; k < ns; ++k) {
        const Molecule mp = displaced(mol, si, dirs[k], h);
        Ep[k] = energy(mp, dirs[k].name + " +", {0.0, 0.0, 0.0});
        Molecule mm;
        if (c2v && dirs[k].antisymmetric) {
            Em[k] = Ep[k];  // E(-S3) = E(+S3) by symmetry
        } else {
            mm = displaced(mol, si, dirs[k], -h);
            Em[k] = energy(mm, dirs[k].name + " -", {0.0, 0.0, 0.0});
        }
        HS(k, k) = (Ep[k] - 2.0 * E0 + Em[k]) / (h * h);
        gS[k] = (Ep[k] - Em[k]) / (2.0 * h);
        if (want_mu) {
            if (mol.atoms.size() == 2) {
                const double up = dipole_along(mp, ez, dirs[k].name + " +"), um = dipole_along(mm, ez, dirs[k].name + " -");
                for (std::size_t c = 0; c < 3; ++c) dmu[k][c] = (up - um) / (2.0 * h) * ez[c];
            } else if (c2v && dirs[k].antisymmetric) {
                const double up = dipole_along(mp, ex, dirs[k].name + " +");  // mu_x(-S3) = -mu_x(+S3)
                for (std::size_t c = 0; c < 3; ++c) dmu[k][c] = up / h * ex[c];
            } else if (c2v) {
                const double up = dipole_along(mp, ez, dirs[k].name + " +"), um = dipole_along(mm, ez, dirs[k].name + " -");
                for (std::size_t c = 0; c < 3; ++c) dmu[k][c] = (up - um) / (2.0 * h) * ez[c];
            } else {
                for (const auto& e : {ez, ex}) {
                    const double up = dipole_along(mp, e, dirs[k].name + " +"), um = dipole_along(mm, e, dirs[k].name + " -");
                    for (std::size_t c = 0; c < 3; ++c) dmu[k][c] += (up - um) / (2.0 * h) * e[c];
                }
            }
        }
    }
    for (std::size_t k = 0; k < ns; ++k)
        for (std::size_t l = k + 1; l < ns; ++l) {
            if (c2v && (dirs[k].antisymmetric || dirs[l].antisymmetric)) continue;  // A1 x B2 vanish
            std::vector<double> q = si.values;
            auto at = [&](double sk, double sl, const char* tag) {
                std::vector<double> qq = q;
                for (std::size_t a = 0; a < nq; ++a) qq[a] += sk * h * dirs[k].d[a] + sl * h * dirs[l].d[a];
                return energy(align_to(mol, geometry_from_internals(mol, qq)),
                              dirs[k].name.substr(0, 2) + "," + dirs[l].name.substr(0, 2) + " " + tag, {0.0, 0.0, 0.0});
            };
            const double epp = at(1, 1, "++"), epm = at(1, -1, "+-"), emp = at(-1, 1, "-+"), emm = at(-1, -1, "--");
            HS(k, l) = HS(l, k) = (epp - epm - emp + emm) / (4.0 * h * h);
        }
    std::vector<double> gq(nq, 0.0);
    for (std::size_t k = 0; k < ns; ++k)
        for (std::size_t a = 0; a < nq; ++a) gq[a] += gS[k] * dirs[k].d[a];
    r.hessian = cartesian_hessian(mol, si, dirs, HS, gq);
    r.has_dipole_derivatives = ho.dipoles;
    if (ho.dipoles) {
        r.dipole_derivatives = apt_from_directions(si, dirs, dmu);
        r.dipole_description = polar ? format("central differences of finite-field dipole moments (field %.0e a.u.)", F)
                                     : "zero by symmetry (homonuclear diatomic)";
    }
    r.description = format("numerical (central differences of %s energies along %zu %sinternal coordinate%s, %d energies%s, "
                           "step %.4f bohr / rad)",
                           method.c_str(), ns, c2v ? "symmetry-adapted " : "", ns == 1 ? "" : "s",
                           r.energy_evaluations - field_energies,
                           field_energies ? format(" + %d finite-field energies for the dipoles", field_energies).c_str() : "",
                           h);
    r.seconds = seconds_since(t0);
    return r;
}

}  // namespace macqc
