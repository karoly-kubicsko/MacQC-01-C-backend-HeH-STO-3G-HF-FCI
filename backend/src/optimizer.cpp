#include "macqc/optimizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/job.h"
#include "macqc/rhf_derivatives.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }
constexpr double kPi = 3.14159265358979323846;

std::string format(const char* fmt, ...) {
    char buf[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    return buf;
}

int period_row(int Z) { return Z <= 2 ? 0 : (Z <= 10 ? 1 : 2); }

double distance(const Molecule& m, int a, int b) {
    double s = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
        const double d = m.atoms[static_cast<std::size_t>(a)].position[k] - m.atoms[static_cast<std::size_t>(b)].position[k];
        s += d * d;
    }
    return std::sqrt(s);
}

Matrix mat_mul(const Matrix& A, const Matrix& B) {
    Matrix C(A.rows(), B.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t k = 0; k < A.cols(); ++k) {
            const double a = A(i, k);
            if (a == 0.0) continue;
            for (std::size_t j = 0; j < B.cols(); ++j) C(i, j) += a * B(k, j);
        }
    return C;
}

double norm(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s);
}

// Symmetric difference |R(O-H1) - R(O-H2)| of a triatomic with the heavy atom first.
double c2v_deviation(const Molecule& m) {
    if (m.atoms.size() != 3) return 0.0;
    return std::fabs(distance(m, 0, 1) - distance(m, 0, 2));
}

}  // namespace

ConvergenceCriteria convergence_criteria(bool tight) {
    ConvergenceCriteria c;
    if (tight) {
        c.max_force = 1.5e-5;
        c.rms_force = 1.0e-5;
        c.max_displacement = 6.0e-5;
        c.rms_displacement = 4.0e-5;
        c.energy = 1.0e-8;
    }
    return c;
}

std::vector<double> lindh_model_hessian(const Molecule& m, const InternalCoordinates& ic) {
    // Lindh, Bernhardsson, Karlström & Malmqvist, Chem. Phys. Lett. 241, 423 (1995).
    static const double alpha[3][3] = {{1.0000, 0.3949, 0.3949}, {0.3949, 0.2800, 0.2800}, {0.3949, 0.2800, 0.2800}};
    static const double rref[3][3] = {{1.35, 2.10, 2.53}, {2.10, 2.87, 3.40}, {2.53, 3.40, 3.40}};
    auto rho = [&](int a, int b) {
        const int ra = period_row(m.atoms[static_cast<std::size_t>(a)].Z), rb = period_row(m.atoms[static_cast<std::size_t>(b)].Z);
        const double r = distance(m, a, b);
        return std::exp(alpha[ra][rb] * (rref[ra][rb] * rref[ra][rb] - r * r));
    };
    std::vector<double> h;
    for (const auto& c : ic.coordinates()) {
        double k = 0.0;
        if (c.type == InternalType::Stretch) k = 0.45 * rho(c.atoms[0], c.atoms[1]);
        else k = 0.15 * rho(c.atoms[0], c.atoms[1]) * rho(c.atoms[1], c.atoms[2]);
        h.push_back(std::max(k, 0.005));
    }
    return h;
}

std::vector<double> schlegel_model_hessian(const Molecule& m, const InternalCoordinates& ic) {
    // Schlegel, Theor. Chim. Acta 66, 333 (1984): F_r = 1.734 / (r - B)^3; bends 0.160 (H) / 0.250.
    static const double Btab[3][3] = {{-0.244, 0.352, 0.660}, {0.352, 1.085, 1.522}, {0.660, 1.522, 2.068}};
    std::vector<double> h;
    for (const auto& c : ic.coordinates()) {
        if (c.type == InternalType::Stretch) {
            const int ra = period_row(m.atoms[static_cast<std::size_t>(c.atoms[0])].Z),
                      rb = period_row(m.atoms[static_cast<std::size_t>(c.atoms[1])].Z);
            const double d = distance(m, c.atoms[0], c.atoms[1]) - Btab[ra][rb];
            h.push_back(std::max(0.005, 1.734 / (d * d * d)));
        } else {
            const bool hydrogen = m.atoms[static_cast<std::size_t>(c.atoms[0])].Z == 1 ||
                                  m.atoms[static_cast<std::size_t>(c.atoms[2])].Z == 1;
            h.push_back(hydrogen ? 0.160 : 0.250);
        }
    }
    return h;
}

std::vector<GeometricParameter> geometric_parameters(const Molecule& mol) {
    std::vector<GeometricParameter> out;
    const InternalCoordinates ic(mol);
    const std::vector<double> q = ic.values(mol);
    for (std::size_t i = 0; i < ic.size(); ++i) {
        const InternalCoordinate& c = ic.coordinates()[i];
        if (c.type == InternalType::Stretch) {
            out.push_back({c.definition(mol), q[i], false});
        } else if (c.type == InternalType::Bend) {
            out.push_back({c.definition(mol), q[i] * 180.0 / kPi, true});
        } else if (c.component == 0) {
            InternalCoordinate b = c;
            b.type = InternalType::Bend;
            const std::string def = "A(" + mol.atoms[static_cast<std::size_t>(c.atoms[0])].symbol + std::to_string(c.atoms[0] + 1) +
                                    "-" + mol.atoms[static_cast<std::size_t>(c.atoms[1])].symbol + std::to_string(c.atoms[1] + 1) +
                                    "-" + mol.atoms[static_cast<std::size_t>(c.atoms[2])].symbol + std::to_string(c.atoms[2] + 1) + ")";
            // angle from the bond vectors
            double u[3], v[3], uu = 0, vv = 0, uv = 0;
            for (std::size_t k = 0; k < 3; ++k) {
                u[k] = mol.atoms[static_cast<std::size_t>(c.atoms[0])].position[k] - mol.atoms[static_cast<std::size_t>(c.atoms[1])].position[k];
                v[k] = mol.atoms[static_cast<std::size_t>(c.atoms[2])].position[k] - mol.atoms[static_cast<std::size_t>(c.atoms[1])].position[k];
                uu += u[k] * u[k];
                vv += v[k] * v[k];
                uv += u[k] * v[k];
            }
            out.push_back({def, std::acos(std::max(-1.0, std::min(1.0, uv / std::sqrt(uu * vv)))) * 180.0 / kPi, true});
        }
    }
    return out;
}

OptimizationResult optimize_geometry(const Molecule& start, const BasisSet& basis, const std::string& method,
                                     const ElectronicOptions& eo, const OptimizationOptions& oo,
                                     const ElectronicCallbacks& cb, const OptimizationCallback& on_step) {
    const auto t0 = Clock::now();
    OptimizationResult res;
    res.criteria = convergence_criteria(oo.tight);
    res.initial = start;
    res.initial_parameters = geometric_parameters(start);
    res.symmetric_start = start.atoms.size() == 3 && c2v_deviation(start) < 1.0e-8;
    const ConvergenceCriteria& cc = res.criteria;

    Molecule x = start;
    InternalCoordinates ic(x);
    Matrix H;
    bool need_initial = true;
    std::vector<double> q_prev, g_prev;
    double e_prev = 0.0, predicted_prev = 0.0, step_prev = 0.0;
    double trust = oo.trust_radius;
    Matrix orbitals;
    const Matrix* guess = nullptr;
    int removed0 = -1;
    std::vector<double> negative_direction;  // of the initial analytic Hessian, if it had a negative eigenvalue
    double negative_curvature = 0.0;

    for (int step = 1; step <= oo.max_steps; ++step) {
        check_cancelled();
        const auto ts = Clock::now();
        OptimizationStep rec;
        rec.step = step;
        rec.geometry = x;
        const BasisSet bx = basis_at_geometry(basis, x);
        GradientResult g = run_gradient(x, bx, method, eo, guess, cb);
        rec.warnings = g.point.warnings;
        if (!g.point.converged) {
            res.message = format("The optimisation stopped at step %d: the %s calculation did not converge.", step,
                                 method.c_str());
            res.final_point = g;
            res.final_geometry = x;
            res.final_energy = g.point.energy;
            res.steps = step;
            break;
        }
        if (removed0 < 0) {
            removed0 = g.point.rhf.removed_functions;
        } else if (g.point.rhf.removed_functions != removed0) {
            res.removed_functions_changed = true;
            rec.warnings.push_back(format(
                "The number of removed near-linear dependencies changed from %d to %d: the potential energy surface is "
                "not continuous here.",
                removed0, g.point.rhf.removed_functions));
            removed0 = g.point.rhf.removed_functions;
        }
        orbitals = g.point.rhf.C;
        guess = &orbitals;
        rec.energy = g.point.energy;
        rec.delta_e = step == 1 ? 0.0 : rec.energy - e_prev;
        rec.gradient_description = g.description;
        rec.cartesian_gradient = g.gradient;

        // ---- internal coordinates, gradient
        if (ic.needs_rebuild(x)) {
            ic = InternalCoordinates(x);
            need_initial = true;
            q_prev.clear();
            rec.warnings.push_back("The internal coordinates were rebuilt (a bend crossed the linear threshold).");
        }
        const std::vector<double> q = ic.values(x);
        const std::size_t nq = q.size(), n3 = 3 * x.atoms.size();
        const Matrix B = ic.wilson_b(x);
        Matrix G(nq, nq);
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nq; ++j) {
                double s = 0.0;
                for (std::size_t k = 0; k < n3; ++k) s += B(i, k) * B(j, k);
                G(i, j) = s;
            }
        int rank = 0;
        const Matrix Gi = generalized_inverse(G, 1.0e-8, &rank);
        const Matrix P = mat_mul(G, Gi);
        const std::vector<double> gq0 = internal_gradient(B, Gi, g.gradient);
        std::vector<double> gq(nq, 0.0);
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nq; ++j) gq[i] += P(i, j) * gq0[j];
        for (const auto& c : ic.coordinates()) {
            rec.names.push_back(c.name());
            rec.definitions.push_back(c.definition(x));
        }
        rec.q = q;
        rec.gq = gq;

        // ---- Hessian: initial or BFGS
        if (need_initial || q_prev.size() != nq) {
            std::string kind = oo.initial_hessian;
            if (kind == "auto")
                kind = (method == "RHF" && bx.n_functions <= oo.analytic_hessian_max_functions) ? "analytic" : "lindh";
            if (kind == "analytic" && method != "RHF") kind = "lindh";
            H = Matrix(nq, nq);
            if (kind == "analytic") {
                const RHFHessianResult hx = rhf_hessian(x, bx, g.point.rhf, eo.threads, eo.memory_mb, 1e-8);
                Matrix Hc = hx.hessian;
                const Matrix K = ic.b_derivative_contracted(x, gq0);
                for (std::size_t a = 0; a < n3; ++a)
                    for (std::size_t b = 0; b < n3; ++b) Hc(a, b) -= K(a, b);
                // H_q = G^- B (H_x - K) B' G^-
                const Matrix GB = mat_mul(Gi, B);
                const Matrix T = mat_mul(GB, Hc);
                for (std::size_t i = 0; i < nq; ++i)
                    for (std::size_t j = 0; j < nq; ++j) {
                        double s = 0.0;
                        for (std::size_t k = 0; k < n3; ++k) s += T(i, k) * GB(j, k);
                        H(i, j) = s;
                    }
                rec.hessian_update = "initial Hessian: analytic RHF Hessian in internal coordinates";
                res.initial_hessian = "analytic RHF Hessian";
                // A minimisation needs a positive definite Hessian: negative eigenvalues (in the
                // non-redundant space) are replaced by their magnitude; the most negative direction
                // is kept, to leave a saddle point along it if the gradient there vanishes.
                const Matrix PHP = mat_mul(mat_mul(P, H), P);
                const EigenResult eh = symmetric_eigen(PHP);
                Matrix Hn(nq, nq);
                int negative = 0;
                for (std::size_t k = 0; k < nq; ++k) {
                    double lam = eh.values[k];
                    double proj = 0.0;
                    for (std::size_t i = 0; i < nq; ++i)
                        for (std::size_t j = 0; j < nq; ++j) proj += eh.vectors(i, k) * P(i, j) * eh.vectors(j, k);
                    if (proj < 0.5) continue;  // redundant direction
                    if (lam < 0.0) {
                        ++negative;
                        if (negative_direction.empty()) {
                            negative_direction.assign(nq, 0.0);
                            for (std::size_t i = 0; i < nq; ++i) negative_direction[i] = eh.vectors(i, k);
                            negative_curvature = lam;
                        }
                        lam = std::max(-lam, 0.005);
                    }
                    lam = std::max(lam, 0.005);
                    for (std::size_t i = 0; i < nq; ++i)
                        for (std::size_t j = 0; j < nq; ++j) Hn(i, j) += lam * eh.vectors(i, k) * eh.vectors(j, k);
                }
                H = Hn;
                if (negative)
                    rec.hessian_update += format("; %d negative eigenvalue%s (lowest %.4f) replaced by %s magnitude for the "
                                                 "minimisation",
                                                 negative, negative == 1 ? "" : "s", negative_curvature,
                                                 negative == 1 ? "its" : "their");
            } else {
                const std::vector<double> d = kind == "schlegel" ? schlegel_model_hessian(x, ic) : lindh_model_hessian(x, ic);
                for (std::size_t i = 0; i < nq; ++i) H(i, i) = d[i];
                rec.hessian_update = std::string("initial Hessian: ") +
                                     (kind == "schlegel" ? "Schlegel model Hessian" : "Lindh model Hessian");
                res.initial_hessian = kind == "schlegel" ? "Schlegel model Hessian" : "Lindh model Hessian";
            }
            need_initial = false;
        } else {
            const std::vector<double> s = ic.difference(q, q_prev);
            std::vector<double> y(nq);
            for (std::size_t i = 0; i < nq; ++i) y[i] = gq[i] - g_prev[i];
            double sy = 0.0;
            for (std::size_t i = 0; i < nq; ++i) sy += s[i] * y[i];
            std::vector<double> Hs(nq, 0.0);
            for (std::size_t i = 0; i < nq; ++i)
                for (std::size_t j = 0; j < nq; ++j) Hs[i] += H(i, j) * s[j];
            double sHs = 0.0;
            for (std::size_t i = 0; i < nq; ++i) sHs += s[i] * Hs[i];
            if (sy > 1.0e-12 * norm(s) * norm(y) && sy > 0.0 && sHs > 1.0e-14) {
                for (std::size_t i = 0; i < nq; ++i)
                    for (std::size_t j = 0; j < nq; ++j) H(i, j) += y[i] * y[j] / sy - Hs[i] * Hs[j] / sHs;
                rec.hessian_update = format("BFGS update (s.y = %.3e)", sy);
            } else {
                rec.hessian_update = format("BFGS update skipped (curvature condition s.y = %.3e not positive)", sy);
            }
            // trust radius from the quality of the previous step
            if (std::fabs(predicted_prev) > 1.0e-14) {
                rec.ratio = rec.delta_e / predicted_prev;
                rec.has_ratio = true;
                if (rec.ratio < 0.25) trust = std::max(oo.min_trust, 0.25 * step_prev);
                else if (rec.ratio > 0.75 && step_prev > 0.8 * trust) trust = std::min(oo.max_trust, 2.0 * trust);
            }
            if (rec.delta_e > 1.0e-6)
                rec.warnings.push_back(format("The energy rose by %.2e E_h; the trust radius is reduced.", rec.delta_e));
        }

        // ---- projected Hessian and RFO step
        Matrix Ht(nq, nq);
        {
            const Matrix PH = mat_mul(mat_mul(P, H), P);
            for (std::size_t i = 0; i < nq; ++i)
                for (std::size_t j = 0; j < nq; ++j) Ht(i, j) = PH(i, j) + 1000.0 * ((i == j ? 1.0 : 0.0) - P(i, j));
        }
        Matrix A(nq + 1, nq + 1);
        for (std::size_t i = 0; i < nq; ++i) {
            for (std::size_t j = 0; j < nq; ++j) A(i, j) = Ht(i, j);
            A(i, nq) = A(nq, i) = gq[i];
        }
        const EigenResult ea = symmetric_eigen(A);
        std::vector<double> dq(nq, 0.0);
        const double last = ea.vectors(nq, 0);
        rec.rfo_lambda = ea.values[0];
        if (std::fabs(last) > 1.0e-8) {
            for (std::size_t i = 0; i < nq; ++i) dq[i] = ea.vectors(i, 0) / last;
        } else {
            for (std::size_t i = 0; i < nq; ++i) dq[i] = -gq[i];
            rec.warnings.push_back("The RFO step was degenerate; a steepest-descent step is used.");
        }
        {
            std::vector<double> pdq(nq, 0.0);
            for (std::size_t i = 0; i < nq; ++i)
                for (std::size_t j = 0; j < nq; ++j) pdq[i] += P(i, j) * dq[j];
            dq = pdq;
        }
        if (step == 1 && !negative_direction.empty()) {
            double gv = 0.0;
            for (std::size_t i = 0; i < nq; ++i) gv += gq[i] * negative_direction[i];
            if (std::fabs(gv) < 1.0e-5) {
                // Saddle point along the negative-curvature direction with no gradient there
                // (e.g. a linear H2O started at 180 degrees): leave it along that direction.
                for (std::size_t i = 0; i < nq; ++i) dq[i] += 0.1 * negative_direction[i];
                rec.warnings.push_back(format(
                    "The starting geometry is a saddle point along a direction with negative curvature (%.4f E_h/unit^2) "
                    "and no gradient; the first step follows that direction.",
                    negative_curvature));
            }
        }
        double sn = norm(dq);
        if (sn > trust) {
            for (double& d : dq) d *= trust / sn;
            rec.step_scaled = true;
            sn = trust;
        }
        rec.trust_radius = trust;
        rec.step_norm = sn;
        double pred = 0.0;
        for (std::size_t i = 0; i < nq; ++i) {
            pred += gq[i] * dq[i];
            for (std::size_t j = 0; j < nq; ++j) pred += 0.5 * dq[i] * Ht(i, j) * dq[j];
        }
        rec.predicted_change = pred;
        rec.dq = dq;

        // ---- convergence
        double mf = 0.0, rf = 0.0, md = 0.0, rd = 0.0;
        for (std::size_t i = 0; i < nq; ++i) {
            mf = std::max(mf, std::fabs(gq[i]));
            rf += gq[i] * gq[i];
            md = std::max(md, std::fabs(dq[i]));
            rd += dq[i] * dq[i];
        }
        const double nd = static_cast<double>(std::max<std::size_t>(1, nq));
        rf = std::sqrt(rf / nd);
        rd = std::sqrt(rd / nd);
        rec.convergence = {{"Maximum force", mf, cc.max_force, mf < cc.max_force, true},
                           {"RMS force", rf, cc.rms_force, rf < cc.rms_force, true},
                           {"Maximum displacement", md, cc.max_displacement, md < cc.max_displacement, true},
                           {"RMS displacement", rd, cc.rms_displacement, rd < cc.rms_displacement, true},
                           {"Energy change", std::fabs(rec.delta_e), cc.energy,
                            step > 1 && std::fabs(rec.delta_e) < cc.energy, step > 1}};
        const bool four = rec.convergence[0].met && rec.convergence[1].met && rec.convergence[2].met &&
                          rec.convergence[3].met;
        const bool energy_ok = step == 1 || rec.convergence[4].met;
        const bool tiny_forces = mf < 0.01 * cc.max_force && rf < 0.01 * cc.rms_force;
        rec.converged = (four && energy_ok) || tiny_forces;
        report_progress(JobStage::Optimization, step, rec.energy, mf);

        if (rec.converged) {
            rec.dq.assign(nq, 0.0);
            rec.seconds = seconds_since(ts);
            res.history.push_back(rec);
            if (on_step) on_step(res.history.back());
            res.converged = true;
            res.steps = step;
            res.final_point = std::move(g);
            res.final_geometry = x;
            res.final_energy = rec.energy;
            res.message = format("Stationary point found after %d step%s.", step, step == 1 ? "" : "s");
            break;
        }
        // ---- take the step
        const BackTransform bt = back_transform(ic, x, dq);
        rec.backtransform_iterations = bt.iterations;
        rec.backtransform_residual = bt.residual;
        if (!bt.converged)
            rec.warnings.push_back(format("The back-transformation to Cartesian coordinates stopped after %d iterations "
                                          "(residual %.1e).",
                                          bt.iterations, bt.residual));
        rec.seconds = seconds_since(ts);
        res.history.push_back(rec);
        if (on_step) on_step(res.history.back());
        q_prev = q;
        g_prev = gq;
        e_prev = rec.energy;
        predicted_prev = pred;
        step_prev = sn;
        res.final_point = std::move(g);
        res.final_geometry = x;
        res.final_energy = rec.energy;
        res.steps = step;
        if (step == oo.max_steps) {
            res.message = format("The optimisation did not converge within %d steps.", oo.max_steps);
            break;
        }
        x = bt.molecule;
    }
    res.final_parameters = geometric_parameters(res.final_geometry);
    res.initial_energy = res.history.empty() ? 0.0 : res.history.front().energy;
    res.symmetry_deviation = c2v_deviation(res.final_geometry);
    res.symmetric_end = res.final_geometry.atoms.size() == 3 && res.symmetry_deviation < 1.0e-6;
    res.seconds = seconds_since(t0);
    return res;
}

}  // namespace macqc
