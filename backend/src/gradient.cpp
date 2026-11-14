#include "macqc/gradient.h"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/mo_integrals.h"
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

RHFOptions rhf_options(const ElectronicOptions& o, const Matrix* guess) {
    RHFOptions r;
    r.max_iter = o.scf_max_iter;
    r.energy_tol = o.scf_energy_tol;
    r.gradient_tol = o.scf_gradient_tol;
    r.lindep_tol = o.lindep_tol;
    r.threads = o.threads;
    r.memory_limit_mb = o.memory_mb;
    r.guess_orbitals = guess;
    r.field = o.field;
    return r;
}

CCOptions cc_options(const ElectronicOptions& o) {
    CCOptions c;
    c.threads = o.threads;
    c.memory_limit_mb = o.fci_memory_mb;
    c.eri_memory_mb = o.memory_mb;
    c.max_iter = o.cc_max_iter;
    c.energy_tol = o.cc_energy_tol;
    c.amplitude_tol = o.cc_amplitude_tol;
    return c;
}

FCIOptions fci_options(const ElectronicOptions& o, bool keep) {
    FCIOptions f;
    f.threads = o.threads;
    f.memory_limit_mb = o.fci_memory_mb;
    f.residual_tol = o.fci_residual_tol;
    f.energy_tol = 1.0e-12;
    f.max_iter = 200;
    f.keep_vector = keep;
    return f;
}

double norm3(const double* a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

void check_rhf(SinglePoint& sp) {
    if (!sp.rhf.converged)
        sp.warnings.push_back(format("The SCF did not converge in %d iterations.", sp.rhf.iterations));
}

// Correlated energy after RHF (MO integrals given when available).
void correlated_energy(SinglePoint& sp, const BasisSet& basis, const ElectronicOptions& o,
                       const ElectronicCallbacks& cb, const MOIntegrals* mo_in) {
    const std::string& m = sp.method;
    const int nocc = sp.rhf.n_occupied;
    if (m == "RHF") {
        sp.energy = sp.rhf.total_energy;
        sp.converged = sp.rhf.converged;
        return;
    }
    if (!sp.rhf.converged) {
        sp.energy = sp.rhf.total_energy;
        sp.converged = false;
        return;
    }
    if (m == "FCI") {
        MOIntegrals local;
        if (!mo_in) local = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
        const MOIntegrals& mo = mo_in ? *mo_in : local;
        sp.fci = run_fci_determinants(mo, nocc, nocc, sp.rhf.nuclear_repulsion, sp.rhf.total_energy,
                                      fci_options(o, false), cb.fci);
        sp.energy = sp.fci->total_energy;
        sp.converged = sp.fci->converged;
        if (!sp.fci->converged)
            sp.warnings.push_back(format("The Full CI (Davidson) did not converge in %d iterations.", sp.fci->iterations));
        return;
    }
    CCMethod ccm = CCMethod::CCSD;
    if (m == "CCSD(T)") ccm = CCMethod::CCSD_T;
    else if (m == "CCSDT") ccm = CCMethod::CCSDT;
    else if (m == "CCSDTQ") ccm = CCMethod::CCSDTQ;
    else if (m != "CCSD") throw std::invalid_argument("unknown method " + m);
    const CCOptions co = cc_options(o);
    if (cc_max_rank(ccm) == 2) {
        const CCIntegrals ints = mo_in ? cc_integrals_from_mo(*mo_in, nocc)
                                       : cc_integrals_from_ao(basis, sp.rhf.H, sp.rhf.C, nocc, o.threads);
        auto ladder = make_ao_ladder(basis, sp.rhf.C, nocc, o.threads, o.memory_mb);
        sp.cc = run_ccsd_closed_shell(ints, *ladder, sp.rhf.total_energy, ccm == CCMethod::CCSD_T, co, cb.cc);
    } else {
        MOIntegrals local;
        if (!mo_in) local = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
        const MOIntegrals& mo = mo_in ? *mo_in : local;
        sp.cc = run_cc_spin_orbital(mo, nocc, cc_max_rank(ccm), sp.rhf.nuclear_repulsion, co, cb.cc);
        sp.cc->reference_energy = sp.rhf.total_energy;
        sp.cc->cc_energy = sp.cc->reference_energy + sp.cc->correlation_energy;
        sp.cc->total_energy = sp.cc->cc_energy;
    }
    sp.energy = sp.cc->total_energy;
    sp.converged = sp.cc->converged;
    if (!sp.cc->converged)
        sp.warnings.push_back(format("The %s iterations did not converge in %d iterations.", sp.cc->method.c_str(),
                                     sp.cc->iterations));
}

}  // namespace

bool analytic_gradient_method(const std::string& m) {
    return m == "RHF" || m == "CCSD" || m == "CCSD(T)" || m == "FCI";
}

double analytic_gradient_bytes(const std::string& m, int nocc, int nmo, int nbasis) {
    if (m == "RHF") return 0.0;
    const double base = correlated_gradient_bytes(static_cast<std::size_t>(nmo), static_cast<std::size_t>(nbasis));
    if (m == "CCSD" || m == "CCSD(T)") return base + ccsd_density_bytes(nocc, nmo - nocc, m == "CCSD(T)");
    if (m == "FCI") return base + fci_density_bytes(nmo, nocc);
    return 0.0;
}

SinglePoint run_single_point(const Molecule& mol, const BasisSet& basis, const std::string& method,
                             const ElectronicOptions& o, const Matrix* guess, const ElectronicCallbacks& cb) {
    const auto t0 = Clock::now();
    SinglePoint sp;
    sp.molecule = mol;
    sp.method = method;
    sp.rhf = run_rhf(mol, basis, rhf_options(o, guess), cb.scf);
    check_rhf(sp);
    correlated_energy(sp, basis, o, cb, nullptr);
    sp.seconds = seconds_since(t0);
    return sp;
}

SimpleInternals simple_internals(const Molecule& mol) {
    SimpleInternals s;
    const std::size_t n = mol.atoms.size();
    auto stretch = [&](std::size_t a, std::size_t b, const std::string& name) {
        double d[3];
        for (int k = 0; k < 3; ++k) d[k] = mol.atoms[b].position[static_cast<std::size_t>(k)] - mol.atoms[a].position[static_cast<std::size_t>(k)];
        const double r = norm3(d);
        std::vector<double> row(3 * n, 0.0);
        for (int k = 0; k < 3; ++k) {
            row[3 * a + static_cast<std::size_t>(k)] = -d[k] / r;
            row[3 * b + static_cast<std::size_t>(k)] = d[k] / r;
        }
        s.names.push_back(name);
        s.values.push_back(r);
        s.B.push_back(row);
    };
    if (n == 2) {
        stretch(0, 1, "R(" + mol.atoms[0].symbol + "1-" + mol.atoms[1].symbol + "2)");
        return s;
    }
    if (n != 3) throw std::invalid_argument("simple_internals: two or three atoms");
    stretch(0, 1, "R(O1-H2)");
    stretch(0, 2, "R(O1-H3)");
    double u[3], v[3];
    for (int k = 0; k < 3; ++k) {
        u[k] = mol.atoms[1].position[static_cast<std::size_t>(k)] - mol.atoms[0].position[static_cast<std::size_t>(k)];
        v[k] = mol.atoms[2].position[static_cast<std::size_t>(k)] - mol.atoms[0].position[static_cast<std::size_t>(k)];
    }
    const double ru = norm3(u), rv = norm3(v);
    for (int k = 0; k < 3; ++k) {
        u[k] /= ru;
        v[k] /= rv;
    }
    const double c = std::max(-1.0, std::min(1.0, u[0] * v[0] + u[1] * v[1] + u[2] * v[2]));
    const double th = std::acos(c), sn = std::sin(th);
    std::vector<double> row(9, 0.0);
    if (sn > 1.0e-8)
        for (int k = 0; k < 3; ++k) {
            const double da = (c * u[k] - v[k]) / (ru * sn), db = (c * v[k] - u[k]) / (rv * sn);
            row[3 + static_cast<std::size_t>(k)] = da;
            row[6 + static_cast<std::size_t>(k)] = db;
            row[static_cast<std::size_t>(k)] = -(da + db);
        }
    s.names.push_back("A(H2-O1-H3)");
    s.values.push_back(th);
    s.B.push_back(row);
    return s;
}

Molecule geometry_from_internals(const Molecule& ref, const std::vector<double>& q) {
    Molecule m = ref;
    const std::size_t n = ref.atoms.size();
    double e1[3];
    for (int k = 0; k < 3; ++k) e1[k] = ref.atoms[1].position[static_cast<std::size_t>(k)] - ref.atoms[0].position[static_cast<std::size_t>(k)];
    const double r = norm3(e1);
    for (double& x : e1) x /= r;
    const auto& O = ref.atoms[0].position;
    if (n == 2) {
        for (int k = 0; k < 3; ++k) m.atoms[1].position[static_cast<std::size_t>(k)] = O[static_cast<std::size_t>(k)] + q[0] * e1[k];
        return m;
    }
    double w[3];
    for (int k = 0; k < 3; ++k) w[k] = ref.atoms[2].position[static_cast<std::size_t>(k)] - O[static_cast<std::size_t>(k)];
    // in-plane unit vector perpendicular to e1, towards atom 3
    const double proj = w[0] * e1[0] + w[1] * e1[1] + w[2] * e1[2];
    double p[3];
    for (int k = 0; k < 3; ++k) p[k] = w[k] - proj * e1[k];
    double pn = norm3(p);
    if (pn < 1.0e-10) {  // linear: any perpendicular direction
        const double a[3] = {std::fabs(e1[0]) < 0.9 ? 1.0 : 0.0, std::fabs(e1[0]) < 0.9 ? 0.0 : 1.0, 0.0};
        const double d = a[0] * e1[0] + a[1] * e1[1] + a[2] * e1[2];
        for (int k = 0; k < 3; ++k) p[k] = a[k] - d * e1[k];
        pn = norm3(p);
    }
    for (double& x : p) x /= pn;
    for (int k = 0; k < 3; ++k) {
        m.atoms[1].position[static_cast<std::size_t>(k)] = O[static_cast<std::size_t>(k)] + q[0] * e1[k];
        m.atoms[2].position[static_cast<std::size_t>(k)] =
            O[static_cast<std::size_t>(k)] + q[1] * (std::cos(q[2]) * e1[k] + std::sin(q[2]) * p[k]);
    }
    return m;
}

GradientResult run_gradient(const Molecule& mol, const BasisSet& basis, const std::string& method,
                            const ElectronicOptions& o, const Matrix* guess, const ElectronicCallbacks& cb,
                            bool force_numerical) {
    const auto t0 = Clock::now();
    GradientResult g;
    const std::size_t n3 = 3 * mol.atoms.size();
    g.gradient.assign(n3, 0.0);
    SinglePoint& sp = g.point;
    sp.molecule = mol;
    sp.method = method;
    sp.rhf = run_rhf(mol, basis, rhf_options(o, guess), cb.scf);
    check_rhf(sp);
    const int nocc = sp.rhf.n_occupied, nmo = sp.rhf.n_mo, nbf = basis.n_functions;

    // ---- analytic or numerical?
    std::string reason;
    if (force_numerical) reason = "requested";
    else if (!analytic_gradient_method(method))
        reason = "no analytic gradient for " + method + " in this version (the general-order engine has no "
                 "Lambda equations)";
    else if (method != "RHF" && sp.rhf.removed_functions > 0)
        reason = format("%d near-linear dependencies were removed; the analytic %s gradient does not contain the "
                        "rotation of the kept orbital space",
                        sp.rhf.removed_functions, method.c_str());
    else if (method != "RHF") {
        const double need = analytic_gradient_bytes(method, nocc, nmo, nbf) / (1024.0 * 1024.0);
        if (need > o.fci_memory_mb)
            reason = format("the analytic gradient would need about %.1f GB (limit %.1f GB)", need / 1024.0,
                            o.fci_memory_mb / 1024.0);
    }

    if (reason.empty() && sp.rhf.converged) {
        g.analytic = true;
        if (method == "RHF") {
            sp.energy = sp.rhf.total_energy;
            sp.converged = true;
            g.rhf_gradient = rhf_gradient(mol, basis, sp.rhf, o.threads, o.lindep_tol);
            g.gradient = g.rhf_gradient->gradient;
            g.description = "analytic (derivative integrals, energy-weighted density" +
                            std::string(sp.rhf.removed_functions ? ", rotation of the kept orbital space" : "") + ")";
            const std::size_t N = sp.rhf.C.rows();
            Matrix P(N, N);
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < nocc; ++k) s += sp.rhf.C(i, static_cast<std::size_t>(k)) * sp.rhf.C(j, static_cast<std::size_t>(k));
                    P(i, j) = 2.0 * s;
                }
            g.dipole = dipole_moment(mol, basis, P);
            g.has_dipole = true;
        } else {
            report_progress(JobStage::Transformation);
            const MOIntegrals mo = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
            CorrelatedDensity dens;
            if (method == "FCI") {
                sp.fci = run_fci_determinants(mo, nocc, nocc, sp.rhf.nuclear_repulsion, sp.rhf.total_energy,
                                              fci_options(o, true), cb.fci);
                sp.energy = sp.fci->total_energy;
                sp.converged = sp.fci->converged;
                if (cb.note) cb.note("Full CI densities from the CI vector");
                dens = fci_densities(nmo, nocc, nocc, sp.fci->vector, o.threads);
                sp.fci->vector.clear();
                sp.fci->vector.shrink_to_fit();
                g.description = "analytic (Full CI one- and two-particle densities, orbital response by Z-vector)";
            } else {
                const bool triples = method == "CCSD(T)";
                const CCIntegrals ints = cc_integrals_from_mo(mo, nocc);
                auto ladder = make_ao_ladder(basis, sp.rhf.C, nocc, o.threads, o.memory_mb);
                const CCOptions co = cc_options(o);
                sp.cc = run_ccsd_closed_shell(ints, *ladder, sp.rhf.total_energy, triples, co, cb.cc);
                sp.energy = sp.cc->total_energy;
                sp.converged = sp.cc->converged;
                if (sp.cc->converged) {
                    CCDensityResult d = ccsd_densities(ints, *ladder, sp.cc->t1, sp.cc->t2, triples, co, cb.lambda);
                    g.lambda_converged = d.lambda_converged;
                    g.lambda_iterations = d.lambda_iterations;
                    g.lambda_residual = d.lambda_residual;
                    if (!d.lambda_converged)
                        sp.warnings.push_back(format("The %s Lambda equations did not converge in %d iterations.",
                                                     method.c_str(), d.lambda_iterations));
                    dens = std::move(d.density);
                }
                g.description = std::string("analytic (") + method +
                                " Lambda equations, relaxed one- and two-particle densities, orbital response by "
                                "Z-vector)";
            }
            if (sp.converged) {
                g.correlated = correlated_gradient(mol, basis, sp.rhf, mo, dens, o.threads, o.memory_mb, cb.cphf);
                g.gradient = g.correlated->gradient;
                g.dipole = g.correlated->relaxed_dipole;
                g.has_dipole = true;
                if (!g.correlated->zvector_converged)
                    sp.warnings.push_back("The Z-vector equations did not converge.");
            }
        }
        g.seconds = seconds_since(t0);
        return g;
    }

    // ---- numerical: central differences of the energy along internal coordinates
    correlated_energy(sp, basis, o, cb, nullptr);
    g.numerical_reason = reason;
    if (!sp.converged) {
        g.seconds = seconds_since(t0);
        return g;
    }
    const SimpleInternals si = simple_internals(mol);
    const double h = o.numerical_step;
    const bool linear = mol.atoms.size() == 3 && std::fabs(std::sin(si.values[2])) < 1.0e-6;
    auto energy_at = [&](const Molecule& m) {
        const SinglePoint p = run_single_point(m, basis_at_geometry(basis, m), method, o, &sp.rhf.C, ElectronicCallbacks());
        ++g.energy_evaluations;
        if (!p.converged) throw std::runtime_error("a displaced calculation of the numerical gradient did not converge");
        return p.energy;
    };
    if (linear) {
        // Cartesian central differences (the bend of a linear triatomic is degenerate).
        for (std::size_t x = 0; x < n3; ++x) {
            Molecule mp = mol, mm = mol;
            mp.atoms[x / 3].position[x % 3] += h;
            mm.atoms[x / 3].position[x % 3] -= h;
            const double ep = energy_at(mp), em = energy_at(mm);
            g.gradient[x] = (ep - em) / (2.0 * h);
            if (cb.note) cb.note(format("  numerical gradient: Cartesian coordinate %zu of %zu", x + 1, n3));
        }
        g.description = format("numerical (central differences of %s energies, %d Cartesian displacements of %.4f bohr)",
                               method.c_str(), g.energy_evaluations, h);
        g.seconds = seconds_since(t0);
        return g;
    }
    const std::size_t nq = si.values.size();
    std::vector<double> gq(nq, 0.0);
    const bool symmetric = nq == 3 && std::fabs(si.values[0] - si.values[1]) < 1.0e-9;
    struct Disp {
        std::string name;
        std::vector<double> dir;  // displacement of (q1, q2, ...) per unit step
    };
    std::vector<Disp> disps;
    if (nq == 1) {
        disps.push_back({si.names[0], {1.0}});
    } else if (symmetric) {
        const double s = 1.0 / std::sqrt(2.0);
        disps.push_back({"S1 = [R(O1-H2) + R(O1-H3)]/sqrt(2) (symmetric stretch)", {s, s, 0.0}});
        disps.push_back({"S2 = A(H2-O1-H3) (bend)", {0.0, 0.0, 1.0}});
    } else {
        for (std::size_t k = 0; k < nq; ++k) {
            std::vector<double> d(nq, 0.0);
            d[k] = 1.0;
            disps.push_back({si.names[k], d});
        }
    }
    for (std::size_t k = 0; k < disps.size(); ++k) {
        std::vector<double> qp = si.values, qm = si.values;
        for (std::size_t j = 0; j < nq; ++j) {
            qp[j] += h * disps[k].dir[j];
            qm[j] -= h * disps[k].dir[j];
        }
        const double ep = energy_at(geometry_from_internals(mol, qp));
        const double em = energy_at(geometry_from_internals(mol, qm));
        const double dEdS = (ep - em) / (2.0 * h);
        // dE/dq_j = sum_k dE/dS_k dS_k/dq_j with orthonormal directions: dS_k/dq = dir_k
        for (std::size_t j = 0; j < nq; ++j) gq[j] += dEdS * disps[k].dir[j];
        g.internal_coordinates.push_back(disps[k].name);
        if (cb.note) cb.note(format("  numerical gradient: %s, dE/dq = %+.8f", disps[k].name.c_str(), dEdS));
    }
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t j = 0; j < nq; ++j) g.gradient[x] += si.B[j][x] * gq[j];
    g.description = format("numerical (central differences of %s energies along %zu %sinternal coordinate%s, %d "
                           "energies, step %.4f bohr / rad)",
                           method.c_str(), disps.size(), symmetric ? "symmetry-adapted " : "",
                           disps.size() == 1 ? "" : "s", g.energy_evaluations, h);
    g.seconds = seconds_since(t0);
    return g;
}

}  // namespace macqc
