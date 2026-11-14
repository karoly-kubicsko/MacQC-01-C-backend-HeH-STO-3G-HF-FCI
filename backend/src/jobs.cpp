#include "macqc/jobs.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>

#include "macqc/cc.h"
#include "macqc/constants.h"
#include "macqc/gto_derivatives.h"
#include "macqc/job.h"
#include "macqc/output.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::string human_time(double s) {
    if (s < 1.0) return "less than a second";
    if (s < 90.0) return strprintf("about %.0f s", s);
    if (s < 5400.0) return strprintf("about %.0f min", s / 60.0);
    if (s < 2.0 * 86400.0) return strprintf("about %.1f h", s / 3600.0);
    return strprintf("about %.1f days", s / 86400.0);
}

std::string method_label(const std::string& m) { return m == "RHF" ? "HF" : (m == "FCI" ? "Full CI" : m); }

ElectronicOptions electronic_options(const MolecularInput& in) {
    ElectronicOptions o;
    o.threads = in.threads;
    o.memory_mb = in.memory_mb;
    o.fci_memory_mb = in.fci_memory_mb;
    o.scf_max_iter = in.max_iter;
    o.scf_energy_tol = std::min(in.energy_tol, 1.0e-10);
    o.scf_gradient_tol = std::min(in.gradient_tol, 1.0e-9);
    return o;
}

// Cost of one gradient (seconds); sets *how. Calibrated against measured times (4 threads,
// x86-64): an RHF gradient costs about 4 RHF energies (3.2-4.0 measured for H2O aug-cc-pVDZ to
// aug-cc-pVQZ); correlated analytic gradients 1.1-2 times the sum of their parts, so that sum
// is scaled by 1.3 plus 0.2 s of set-up; numerical gradients (n + 1) energies, each with a fixed
// set-up time (about 0.25 s for the general-order CCSDT/CCSDTQ engine).
double gradient_seconds(const MolecularResult& p, bool* analytic, std::string* how, int* energies) {
    const std::string& m = p.input.method;
    const CostEstimate& e = p.estimate;
    const int nbf = p.basis.n_functions, nocc = p.molecule.electrons() / 2;
    const bool c2v_or_diatomic = true;  // the sidebar geometries are C2v
    const int numerical_energies = p.molecule.atoms.size() == 2 ? 2 : (c2v_or_diatomic ? 4 : 6);
    *energies = 0;
    if (analytic_gradient_method(m)) {
        const double need = analytic_gradient_bytes(m, nocc, nbf, nbf) / (1024.0 * 1024.0);
        if (need <= p.input.fci_memory_mb) {
            *analytic = true;
            if (m == "RHF") {
                *how = "analytic gradient (derivative integrals)";
                return 4.0 * e.rhf_seconds;
            }
            const double n = nbf;
            const double n5 = 6.0 * std::pow(n, 5) / (4.0e9 * std::max(1, e.threads));
            if (m == "FCI") {
                *how = "analytic gradient (Full CI densities, Z-vector)";
                return 1.3 * (e.total_seconds * 1.3 + 2.5 * e.rhf_seconds + n5) + 0.2;
            }
            *how = "analytic gradient (" + m + " Lambda equations, relaxed densities, Z-vector)";
            double t = e.total_seconds + 2.5 * e.rhf_seconds + n5;
            if (e.cc) t += e.cc->iterations * e.cc->iteration_seconds + 6.0 * e.cc->triples_seconds;
            return 1.3 * t + 0.2;
        }
        *how = strprintf("numerical gradient (%d %s energies; the analytic gradient would need %.1f GB)",
                         numerical_energies, m.c_str(), need / 1024.0);
    } else {
        *how = strprintf("numerical gradient (%d %s energies per gradient, central differences)", numerical_energies,
                         m.c_str());
    }
    *analytic = false;
    *energies = numerical_energies;
    const double setup = (m == "CCSDT" || m == "CCSDTQ") ? 0.25 : 0.05;
    return (numerical_energies + 1) * (e.total_seconds + setup);
}

}  // namespace

std::string canonical_job_type(const std::string& type) {
    std::string t;
    for (char c : type)
        if (!std::isspace(static_cast<unsigned char>(c)) && c != '_' && c != '-')
            t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (t == "energy" || t == "sp" || t == "singlepoint" || t == "1") return "energy";
    if (t == "opt" || t == "optimization" || t == "optimisation" || t == "geometryoptimization" ||
        t == "geometryoptimisation" || t == "2")
        return "opt";
    if (t == "freq" || t == "frequencies" || t == "frequency" || t == "3") return "freq";
    if (t == "opt+freq" || t == "optfreq" || t == "4") return "opt+freq";
    if (t == "opt+freq+sp" || t == "optfreqsp" || t == "5") return "opt+freq+sp";
    return "";
}

std::string job_title(const std::string& t) {
    if (t == "opt") return "Geometry optimisation";
    if (t == "freq") return "Frequencies";
    if (t == "opt+freq") return "Optimisation + frequencies";
    if (t == "opt+freq+sp") return "Optimisation + frequencies + single-point energy";
    return "Energy";
}

bool job_has_optimization(const std::string& t) { return t == "opt" || t == "opt+freq" || t == "opt+freq+sp"; }
bool job_has_frequencies(const std::string& t) { return t == "freq" || t == "opt+freq" || t == "opt+freq+sp"; }
bool job_has_single_point(const std::string& t) { return t == "opt+freq+sp"; }

JobEstimate estimate_job(const MolecularResult& p, const JobSettings& s, const MolecularResult* sp) {
    JobEstimate est;
    const std::string& m = p.input.method;
    const std::string& t = s.type;
    const CostEstimate& e = p.estimate;
    if (!p.refusal.empty()) est.refusal = p.refusal;
    const int nbf = p.basis.n_functions, nocc = p.molecule.electrons() / 2;
    if (t == "energy") {
        JobStepEstimate st{"Energy", m, p.input.basis, "single point", true, 1, "energies", e.total_seconds, 0.0};
        est.steps.push_back(st);
        est.total_seconds = e.total_seconds;
        return est;
    }
    bool analytic = false;
    std::string how;
    int grad_energies = 0;
    const double tg = gradient_seconds(p, &analytic, &how, &grad_energies);
    const double grad_mb = analytic && m != "RHF" ? analytic_gradient_bytes(m, nocc, nbf, nbf) / (1024.0 * 1024.0) : 0.0;
    est.memory_mb = grad_mb;
    if (e.fci) est.memory_mb = std::max(est.memory_mb, e.fci->memory_bytes / (1024.0 * 1024.0));
    if (e.cc) est.memory_mb = std::max(est.memory_mb, e.cc->memory_bytes / (1024.0 * 1024.0));
    const int typical_steps = 6;
    if (job_has_optimization(t)) {
        JobStepEstimate st;
        st.name = "Geometry optimisation";
        st.method = m;
        st.basis = p.input.basis;
        st.analytic = analytic;
        st.how = how + strprintf("; typically %d steps (at most %d)", typical_steps, s.max_steps);
        st.evaluations = typical_steps;
        st.evaluation_kind = "steps";
        st.seconds = typical_steps * tg;
        if (m == "RHF" && nbf <= 100) st.seconds += 6.0 * tg;  // analytic initial Hessian
        st.memory_mb = grad_mb;
        est.steps.push_back(st);
    }
    if (job_has_frequencies(t)) {
        JobStepEstimate st;
        st.name = "Hessian";
        st.method = m;
        st.basis = p.input.basis;
        const bool diatomic = p.molecule.atoms.size() == 2;
        if (m == "RHF") {
            st.analytic = true;
            st.how = "analytic Hessian (second-derivative integrals, CPHF), analytic dipole derivatives";
            st.evaluations = 1;
            st.evaluation_kind = "Hessians";
            st.seconds = 7.0 * tg;
        } else if (analytic) {
            st.analytic = false;
            st.evaluations = diatomic ? 2 : 6;
            st.evaluation_kind = "gradients";
            st.how = strprintf("numerical Hessian: central differences of %d analytic gradients along %s; IR intensities "
                               "from the relaxed dipole moments",
                               st.evaluations, diatomic ? "the bond" : "the symmetry-adapted coordinates S1, S2, S3");
            st.seconds = (st.evaluations + 1) * tg;
        } else {
            st.analytic = false;
            const NumericalHessianCost c = numerical_hessian_cost(p.molecule, false, true);
            st.evaluations = c.energies + c.field_energies;
            st.evaluation_kind = "energies";
            st.how = c.field_energies > 0
                         ? strprintf("numerical Hessian: %d %s energies (second differences) + %d finite-field energies for "
                                     "the dipole derivatives",
                                     c.energies, m.c_str(), c.field_energies)
                         : strprintf("numerical Hessian: %d %s energies (second differences); no IR intensities needed "
                                     "(no dipole moment change by symmetry)",
                                     c.energies, m.c_str());
            st.seconds = (st.evaluations + 1) * (e.total_seconds + ((m == "CCSDT" || m == "CCSDTQ") ? 0.25 : 0.05));
        }
        est.steps.push_back(st);
    }
    if (job_has_single_point(t) && sp) {
        JobStepEstimate st;
        st.name = "Single-point energy";
        st.method = sp->input.method;
        st.basis = sp->input.basis;
        st.how = "single point at the optimised geometry";
        st.evaluations = 1;
        st.evaluation_kind = "energies";
        st.seconds = sp->estimate.total_seconds;
        if (sp->estimate.fci) st.memory_mb = sp->estimate.fci->memory_bytes / (1024.0 * 1024.0);
        if (sp->estimate.cc) st.memory_mb = std::max(st.memory_mb, sp->estimate.cc->memory_bytes / (1024.0 * 1024.0));
        est.memory_mb = std::max(est.memory_mb, st.memory_mb);
        est.steps.push_back(st);
        if (!sp->refusal.empty() && est.refusal.empty())
            est.refusal = "The single-point step is not possible: " + sp->refusal;
    }
    for (const auto& st : est.steps) est.total_seconds += st.seconds;
    // refusals
    if (est.refusal.empty() && m == "CCSDTQ" && p.molecule.name == "H2O") {
        const NumericalHessianCost c = numerical_hessian_cost(p.molecule, false, true);
        est.refusal = strprintf(
            "CCSDTQ %s for H2O is not possible: CCSDTQ has no analytic gradient, so %s would need %s CCSDTQ energies "
            "of %s each (%.1f GB of memory). Use CCSD(T) or CCSDT for the geometry and the frequencies, and CCSDTQ in "
            "the single-point step of job 5.",
            job_has_frequencies(t) ? "frequencies" : "optimisation",
            job_has_frequencies(t) ? "the Hessian" : "every optimisation step",
            job_has_frequencies(t) ? std::to_string(c.energies + c.field_energies).c_str() : "4",
            human_time(e.total_seconds).c_str(), e.cc ? e.cc->memory_bytes / 1073741824.0 : 0.0);
    }
    if (est.refusal.empty() && est.total_seconds > kMaxCCSeconds)
        est.refusal = strprintf("The %s would take %s (limit: three days). Use a smaller basis set or a cheaper method.",
                                job_title(t).c_str(), human_time(est.total_seconds).c_str());
    if (est.refusal.empty() && est.total_seconds > kLongCalculationSeconds)
        est.warnings.push_back(strprintf("Long calculation: %s with %s/%s, estimated %s on %d threads.", job_title(t).c_str(),
                                         method_label(m).c_str(), p.input.basis.c_str(), human_time(est.total_seconds).c_str(),
                                         e.threads));
    for (const auto& st : est.steps)
        est.notes.push_back(st.name + ": " + st.how + (st.seconds > 0 ? " — " + human_time(st.seconds) : ""));
    return est;
}

namespace {

// Text shared by both job paths.
void write_header(OutputWriter& w, const JobSettings& s, const std::string& echo, int threads) {
    const std::string rule(88, '=');
    w.line(rule);
    w.line("                                ΨLab " + std::string(kProgramVersion));
    w.line("        restricted Hartree-Fock, coupled cluster and Full CI for H2, LiH and H2O");
    w.line(rule);
    w.line("  Program      ΨLab " + std::string(kProgramVersion) + " (C++17 backend; integrals, SCF, CC, CI, derivatives)");
    w.line("  Started      " + local_time_string());
    w.line("  Job          " + job_title(s.type) + " (" + s.type + ")");
    w.line(strprintf("  Threads      %d", threads));
    w.section("Input");
    w.write(echo);
    if (!echo.empty() && echo.back() != '\n') w.line();
}

void write_molecule(OutputWriter& w, const MolecularResult& p, const std::string& what) {
    const Molecule& m = p.molecule;
    w.section("Molecule");
    w.line(strprintf("    %s, %d electrons, charge %d, closed-shell singlet", m.name.c_str(), m.electrons(), m.charge));
    w.line("    " + what);
    w.write(geometry_table(m));
    w.write(parameter_lines(m));
    w.line(strprintf("    Nuclear repulsion energy %20.10f E_h", m.nuclear_repulsion()));
    w.section("Basis set");
    w.write(basis_block(p.basis, p.notes));
    const CostEstimate& e = p.estimate;
    w.line(strprintf("    two-electron integrals: %zu unique, %s", e.unique_eris,
                     e.eri_stored ? strprintf("kept in memory (%.1f MB)", e.eri_memory_mb).c_str()
                                  : strprintf("integral-direct (storing them would need %.1f GB)", e.eri_memory_mb / 1024.0).c_str()));
}

void write_estimate(OutputWriter& w, const MolecularResult& p, const JobEstimate& est) {
    w.section("Size, memory and time (estimates)");
    const CostEstimate& e = p.estimate;
    auto count = [](double n) { return n < 1.0e9 ? strprintf("%.0f", n) : strprintf("%.3g", n); };
    auto bytes = [](double b) {
        const double mb = b / 1048576.0;
        if (mb < 1024.0) return strprintf("%.0f MB", mb);
        if (mb < 1048576.0) return strprintf("%.1f GB", mb / 1024.0);
        return strprintf("%.1f TB", mb / 1048576.0);
    };
    if (e.fci)
        w.line(strprintf("    Full CI: %s determinants (%d electrons in %d orbitals), about %s", count(e.fci->determinants).c_str(),
                         e.fci->alpha + e.fci->beta, e.fci->orbitals, bytes(e.fci->memory_bytes).c_str()));
    if (e.cc)
        w.line(strprintf("    %s: %d occupied, %d virtual orbitals, %s amplitudes, about %s", e.cc->method.c_str(),
                         e.cc->occupied, e.cc->virtuals, count(e.cc->amplitudes).c_str(), bytes(e.cc->memory_bytes).c_str()));
    for (const auto& n : est.notes) w.line("    " + n);
    if (est.memory_mb > 0) w.line(strprintf("    memory of the derivative steps: about %.0f MB", est.memory_mb));
    w.line("    total: " + human_time(est.total_seconds) + strprintf(" on %d threads", e.threads));
    const bool job_level = !(est.steps.size() == 1 && est.steps[0].name == "Energy");
    for (const auto& x : p.warnings)  // the job-level estimate replaces the single-energy "long calculation" warning
        if (!(job_level && x.rfind("Long calculation", 0) == 0)) w.line("    Warning: " + x);
    for (const auto& x : est.warnings) w.line("    Warning: " + x);
    if (!est.refusal.empty()) w.line("    NOT POSSIBLE: " + est.refusal);
}

void write_correlated(OutputWriter& w, const SinglePoint& sp) {
    if (sp.cc) w.write(cc_report_section(*sp.cc, 1));
    if (sp.fci) w.write(fci_report_section(*sp.fci, sp.rhf.total_energy, 1));
}

void write_termination(OutputWriter& w, JobResult& r, double seconds) {
    w.line();
    const std::string rule(88, '=');
    w.line(rule);
    if (r.termination == "normal")
        w.line("  Normal termination of ΨLab " + std::string(kProgramVersion) + " at " + local_time_string() +
               "; total time " + duration_string(seconds) + ".");
    else if (r.termination == "cancelled")
        w.line("  Calculation cancelled at " + local_time_string() + " after " + duration_string(seconds) + ".");
    else
        w.line("  Error termination of ΨLab " + std::string(kProgramVersion) + " at " + local_time_string() + ": " + r.message);
    w.line(rule);
}

struct Streams {
    OutputWriter* w = nullptr;
    std::string context;  // progress-text prefix, e.g. "Optimisation step 3 · "
    std::string method;
    bool scf_started = false, cc_started = false, fci_started = false, lambda_started = false, cphf_started = false;
    double scf_energy = 0.0;
    int scf_iterations = 0;
    void reset() { scf_started = cc_started = fci_started = lambda_started = cphf_started = false; }
    void scf_done_line() {
        if (scf_started) w->line(strprintf("    SCF converged: E(RHF) = %.10f E_h after %d iterations", scf_energy, scf_iterations));
    }
};

ElectronicCallbacks stream_callbacks(Streams& s) {
    ElectronicCallbacks cb;
    cb.scf = [&s](const RHFIteration& it) {
        if (!s.scf_started) {
            s.w->line("\n    SCF (restricted Hartree-Fock, DIIS)");
            s.w->write(scf_header());
            s.scf_started = true;
        }
        s.w->write(scf_iteration_line(it));
        s.scf_energy = it.energy;
        s.scf_iterations = it.iteration;
        report_progress_text(s.context + strprintf("SCF iteration %d · E = %.10f E_h", it.iteration, it.energy));
    };
    cb.cc = [&s](const CCIteration& it) {
        if (!s.cc_started) {
            s.scf_done_line();
            s.w->line();
            s.w->write(cc_header(s.method == "CCSD(T)" ? "CCSD" : s.method));
            s.cc_started = true;
        }
        s.w->write(cc_iteration_line(it));
        report_progress_text(s.context + strprintf("%s iteration %d · E = %.10f E_h · max |dt| %.1e",
                                                   s.method == "CCSD(T)" ? "CCSD" : s.method.c_str(), it.iteration,
                                                   it.energy, it.residual));
    };
    cb.fci = [&s](const FCIIteration& it) {
        if (!s.fci_started) {
            s.scf_done_line();
            s.w->line();
            s.w->write(fci_header());
            s.fci_started = true;
        }
        s.w->write(fci_iteration_line(it));
        report_progress_text(s.context + strprintf("Davidson iteration %d · E = %.10f E_h · residual %.1e", it.iteration,
                                                   it.energy, it.residual));
    };
    cb.lambda = [&s](const LambdaIteration& it) {
        if (!s.lambda_started) {
            s.w->line("\n    " + s.method + " Lambda equations (Jacobi updates with DIIS)\n    iter   max |d lambda|");
            s.lambda_started = true;
        }
        s.w->line(strprintf("    %4d   %12.3e", it.iteration, it.residual));
        report_progress_text(s.context + strprintf("Lambda iteration %d · max change %.1e", it.iteration, it.residual));
    };
    cb.cphf = [&s](const CPHFProgress& it) {
        if (!s.cphf_started) {
            s.w->line("\n    Orbital response (coupled-perturbed HF / Z-vector, Krylov subspace)\n    iter   largest residual");
            s.cphf_started = true;
        }
        s.w->line(strprintf("    %4d   %12.3e", it.iteration, it.residual));
        report_progress_text(s.context + strprintf("Orbital response iteration %d · residual %.1e", it.iteration, it.residual));
    };
    cb.note = [&s](const std::string& text) {
        s.w->line("    " + text);
        report_progress_text(s.context + text);
    };
    return cb;
}

void write_point_details(OutputWriter& w, const SinglePoint& sp, const BasisSet& b, bool full_orbitals) {
    w.line();
    w.write(scf_result_block(sp.rhf, sp.molecule, b, full_orbitals));
    write_correlated(w, sp);
    for (const auto& x : sp.warnings) w.line("    Warning: " + x);
}

std::string final_summary(const JobResult& r, const MolecularResult& final) {
    std::string s;
    const Molecule& m = r.final_geometry;
    s += strprintf("    Job: %s\n", job_title(r.settings.type).c_str());
    s += "    Final geometry:\n" + parameter_lines(m, "      ");
    s += strprintf("    Final energy (%s): %.10f E_h\n", r.final_energy_label.c_str(), r.final_energy);
    if (r.optimization)
        s += strprintf("    Optimisation: %s, %d steps\n", r.optimization->converged ? "converged" : "NOT converged",
                       r.optimization->steps);
    if (r.harmonic) {
        s += "    Harmonic wavenumbers (cm-1):";
        for (const auto& md : r.harmonic->modes)
            s += md.imaginary ? strprintf(" %.2fi", std::fabs(md.wavenumber)) : strprintf(" %.2f", md.wavenumber);
        s += "\n";
    }
    if (r.thermo) {
        const ThermoResult& t = *r.thermo;
        const double base = r.final_energy;
        s += strprintf("    Thermochemistry at %.2f K and %.0f Pa (corrections from %s):\n", t.temperature, t.pressure,
                       (method_label(r.point ? r.point->method : final.input.method) + "/" + final.input.basis).c_str());
        s += strprintf("      E + ZPE %20.10f E_h\n", base + joule_per_mol_to_hartree(t.zero_point_energy));
        s += strprintf("      U       %20.10f E_h\n", base + joule_per_mol_to_hartree(t.thermal_energy));
        s += strprintf("      H       %20.10f E_h\n", base + joule_per_mol_to_hartree(t.enthalpy));
        s += strprintf("      G       %20.10f E_h\n", base + joule_per_mol_to_hartree(t.gibbs));
        s += strprintf("      S       %20.4f J/(mol K)\n", t.total.entropy);
    }
    return s;
}

std::string summary_lines(const MolecularResult& m) {
    // The short summary of versions 3.0 and 4.0 (also at the top of psilab's output file).
    std::string s;
    if (m.rhf) {
        s += strprintf("    Hartree-Fock (RHF) total energy:  %.10f E_h (%d iterations)\n", m.rhf->total_energy, m.rhf->iterations);
        if (m.cc) s += strprintf("    %s total energy:  %.10f E_h\n", m.cc->method.c_str(), m.cc->total_energy);
        if (m.fci) s += strprintf("    Full CI total energy:  %.10f E_h\n", m.fci->total_energy);
    }
    return s;
}

}  // namespace

std::string job_estimate_text(const MolecularResult& p, const JobEstimate& est) {
    std::string s;
    OutputWriter w([&s](const std::string& t) { s += t; });
    write_estimate(w, p, est);
    return s;
}

MolecularResult final_point_result(const JobResult& job, const MolecularInput& input) {
    MolecularResult r;
    r.input = input;
    const SinglePoint* sp = job.single_point ? &*job.single_point : (job.point ? &*job.point : nullptr);
    r.input.method = sp ? sp->method : input.method;
    if (job.single_point) r.input.basis = job.single_point_basis;
    r.molecule = sp ? sp->molecule : job.final_geometry;
    if (r.molecule.name.empty()) r.molecule.name = canonical_molecule_name(input.molecule);
    // geometry parameters of the final point (R of the first bond, the H-O-H angle)
    const SimpleInternals si = simple_internals(r.molecule);
    r.input.R = si.values[0];
    if (si.values.size() == 3) r.input.angle = si.values[2] * 180.0 / 3.14159265358979323846;
    r.basis = build_basis(r.molecule, r.input.basis, input.basis_directory);
    for (const auto& a : r.basis.atoms)
        if (!a.note.empty() && std::find(r.notes.begin(), r.notes.end(), a.note) == r.notes.end()) r.notes.push_back(a.note);
    if (sp) {
        r.rhf = sp->rhf;
        r.cc = sp->cc;
        r.fci = sp->fci;
        r.warnings = sp->warnings;
    }
    r.report = molecular_report(r);
    return r;
}

JobResult run_energy_job(MolecularResult& p, const JobSettings& s, const OutputSink& sink, const std::string& echo) {
    const auto t0 = Clock::now();
    JobResult r;
    r.settings = s;
    OutputWriter w(sink);
    write_header(w, s, echo, p.estimate.threads);
    write_molecule(w, p, "Geometry");
    r.estimate = estimate_job(p, s, nullptr);
    write_estimate(w, p, r.estimate);
    r.initial_geometry = r.final_geometry = p.molecule;
    if (!p.refusal.empty()) {
        r.termination = "refused";
        r.message = p.refusal;
        write_termination(w, r, seconds_since(t0));
        r.output = w.text();
        return r;
    }
    Streams st;
    st.w = &w;
    st.method = p.input.method;
    w.section("Energy: " + method_label(p.input.method) + "/" + p.input.basis);
    const ElectronicCallbacks cb = stream_callbacks(st);
    try {
        run_molecular(p, cb.scf, cb.fci, cb.cc);
        SinglePoint sp;
        sp.molecule = p.molecule;
        sp.method = p.input.method;
        sp.rhf = *p.rhf;
        sp.cc = p.cc;
        sp.fci = p.fci;
        sp.energy = p.cc ? p.cc->total_energy : (p.fci ? p.fci->total_energy : p.rhf->total_energy);
        sp.converged = p.rhf->converged && (!p.cc || p.cc->converged) && (!p.fci || p.fci->converged);
        sp.warnings = p.warnings;
        w.subsection("Results");
        write_point_details(w, sp, p.basis, true);
        r.point = sp;
        r.final_energy = sp.energy;
        r.final_energy_label = method_label(sp.method) + "/" + p.input.basis;
        w.section("Summary");
        w.write(summary_lines(p));
        w.write(final_summary(r, p));
        r.completed = sp.converged;
        r.termination = sp.converged ? "normal" : "error";
        if (!sp.converged) r.message = "the calculation did not converge";
    } catch (const Cancelled&) {
        r.termination = "cancelled";
        write_termination(w, r, seconds_since(t0));
        r.output = w.text();
        throw;
    } catch (const std::exception& e) {
        r.termination = "error";
        r.message = e.what();
    }
    r.seconds = seconds_since(t0);
    write_termination(w, r, r.seconds);
    r.output = w.text();
    return r;
}

JobResult run_job(const MolecularInput& input, const JobSettings& s, const OutputSink& sink, const std::string& echo) {
    const auto t0 = Clock::now();
    JobResult r;
    r.settings = s;
    MolecularResult p = prepare_molecular(input);
    std::optional<MolecularResult> psp;
    if (job_has_single_point(s.type)) {
        MolecularInput in2 = input;
        in2.method = s.sp_method.empty() ? input.method : s.sp_method;
        in2.basis = s.sp_basis.empty() ? input.basis : s.sp_basis;
        psp = prepare_molecular(in2);
        r.single_point_method = psp->input.method;
        r.single_point_basis = psp->input.basis;
    }
    OutputWriter w(sink);
    write_header(w, s, echo, p.estimate.threads);
    write_molecule(w, p, job_has_optimization(s.type) ? "Starting geometry" : "Geometry (input)");
    r.estimate = estimate_job(p, s, psp ? &*psp : nullptr);
    write_estimate(w, p, r.estimate);
    r.initial_geometry = r.final_geometry = p.molecule;
    if (!r.estimate.refusal.empty()) {
        r.termination = "refused";
        r.message = r.estimate.refusal;
        write_termination(w, r, seconds_since(t0));
        r.output = w.text();
        return r;
    }
    const std::string method = p.input.method;
    const ElectronicOptions eo = electronic_options(p.input);
    Streams st;
    st.w = &w;
    st.method = method;
    ElectronicCallbacks cb = stream_callbacks(st);
    Matrix orbitals;
    const Matrix* guess = nullptr;
    try {
        Molecule x = p.molecule;
        // ---- optimisation
        if (job_has_optimization(s.type)) {
            w.section("Geometry optimisation: " + method_label(method) + "/" + p.input.basis);
            const ConvergenceCriteria c = convergence_criteria(s.tight);
            w.line(strprintf("    Redundant internal coordinates, RFO steps with a trust radius, BFGS updates; %s "
                             "convergence, at most %d steps",
                             s.tight ? "tight" : "normal", s.max_steps));
            w.line(strprintf("    thresholds: max force %.1e, RMS force %.1e E_h/bohr; max step %.1e, RMS step %.1e bohr; "
                             "energy change %.0e E_h",
                             c.max_force, c.rms_force, c.max_displacement, c.rms_displacement, c.energy));
            OptimizationOptions oo;
            oo.tight = s.tight;
            oo.max_steps = s.max_steps;
            oo.on_step_start = [&](int step, const Molecule& geom) {
                st.reset();
                st.context = strprintf("Optimisation step %d · ", step);
                w.subsection(strprintf("Optimisation step %d", step));
                w.write(geometry_table(geom));
                report_progress_text(st.context + "energy and gradient");
            };
            const OptimizationResult opt = optimize_geometry(x, p.basis, method, eo, oo, cb, [&](const OptimizationStep& stp) {
                w.line(stp.step == 1 ? strprintf("\n    Energy %.10f E_h", stp.energy)
                                     : strprintf("\n    Energy %.10f E_h (change %+.3e E_h)", stp.energy, stp.delta_e));
                w.write(gradient_block(stp.geometry, stp.cartesian_gradient, stp.gradient_description));
                w.write(optimization_step_block(stp, c));
                report_progress_text(strprintf("Optimisation step %d · E = %.10f E_h · max force %.2e", stp.step, stp.energy,
                                               stp.convergence[0].value));
            });
            w.write(optimization_result_block(opt));
            r.optimization = opt;
            x = opt.final_geometry;
            r.final_geometry = x;
            r.final_energy = opt.final_energy;
            r.final_energy_label = method_label(method) + "/" + p.input.basis;
            orbitals = opt.final_point.point.rhf.C;
            guess = &orbitals;
            {
                SinglePoint sp = opt.final_point.point;
                r.point = sp;
                st.reset();
                w.subsection("Electronic structure at the " + std::string(opt.converged ? "optimised" : "last") + " geometry");
                write_point_details(w, sp, basis_at_geometry(p.basis, x), true);
                if (opt.final_point.has_dipole && method != "RHF") {
                    const auto& d = opt.final_point.dipole;
                    const double nd = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                    w.line("\n    Dipole moment (relaxed " + method + " density, origin at 0)");
                    w.line(strprintf("      x %12.8f   y %12.8f   z %12.8f   |mu| %12.8f e bohr   (%.6f D)", d[0], d[1], d[2],
                                     nd, nd / codata2018::kDebyeInAtomicUnits));
                }
            }
            if (!opt.converged) {
                r.termination = "error";
                r.message = opt.message + (s.type != "opt" ? " The frequency and single-point steps were not run: they "
                                                             "need a stationary point."
                                                           : "");
                throw std::runtime_error(r.message);
            }
        }
        // ---- frequencies
        if (job_has_frequencies(s.type)) {
            const bool at_input = s.type == "freq";
            w.section("Harmonic frequencies: " + method_label(method) + "/" + p.input.basis);
            w.line(std::string("    at the ") + (at_input ? "input geometry" : "optimised geometry"));
            w.write(geometry_table(x));
            st.reset();
            st.context = "Hessian · ";
            HessianOptions ho;
            if (r.optimization && r.optimization->converged) {
                ho.reference = &r.optimization->final_point;  // energy and gradient of the last step
                w.line("    energy, gradient and orbitals at this geometry: from the last optimisation step");
            }
            int last_index = 0;
            const HessianResult h = run_hessian(x, basis_at_geometry(p.basis, x), method, eo, ho, cb, guess,
                                                [&](int i, int n, const std::string& what) {
                                                    if (i == 1) w.line("\n    Displaced calculations");
                                                    last_index = i;
                                                    w.line(strprintf("    displacement %d of %d: %s", i, n, what.c_str()));
                                                    report_progress_text(strprintf("Hessian: displacement %d of %d (%s)", i, n, what.c_str()));
                                                    st.reset();
                                                });
            (void)last_index;
            r.hessian = h;
            if (at_input) {
                r.point = h.reference.point;
                r.final_energy = h.reference.point.energy;
                r.final_energy_label = method_label(method) + "/" + p.input.basis;
                w.subsection("Electronic structure at the input geometry");
                write_point_details(w, h.reference.point, basis_at_geometry(p.basis, x), true);
                if (h.reference.has_dipole && method != "RHF") {
                    const auto& d = h.reference.dipole;
                    const double nd = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                    w.line("\n    Dipole moment (relaxed " + method + " density, origin at 0)");
                    w.line(strprintf("      x %12.8f   y %12.8f   z %12.8f   |mu| %12.8f e bohr   (%.6f D)", d[0], d[1], d[2],
                                     nd, nd / codata2018::kDebyeInAtomicUnits));
                }
                w.line();
                w.write(gradient_block(x, h.reference.gradient, h.reference.description));
                double mx = 0.0;
                for (double g : h.reference.gradient) mx = std::max(mx, std::fabs(g));
                if (mx > convergence_criteria(false).max_force) {
                    const std::string warn = strprintf(
                        "The input geometry is not a stationary point (largest gradient component %.2e E_h/bohr): the "
                        "frequencies and the thermochemistry refer to a non-stationary geometry. Optimise first "
                        "(job opt+freq).",
                        mx);
                    w.line("    Warning: " + warn);
                    r.warnings.push_back(warn);
                }
            }
            HarmonicAnalysis a = harmonic_analysis(x, h.hessian, h.has_dipole_derivatives ? &h.dipole_derivatives : nullptr);
            w.line();
            w.write(frequency_block(x, h, a));
            report_progress_text("Thermochemistry");
            ThermoResult t = thermochemistry(x, a, r.final_energy, s.temperature, s.pressure, symmetry_number(x));
            w.section("Thermochemistry (ideal gas, rigid rotor, harmonic oscillators)");
            w.write(thermochemistry_block(t));
            for (const auto& x2 : t.warnings) r.warnings.push_back(x2);
            r.harmonic = a;
            r.thermo = t;
        }
        // ---- single point
        if (job_has_single_point(s.type)) {
            const std::string m2 = psp->input.method, b2 = psp->input.basis;
            w.section("Single-point energy at the optimised geometry: " + method_label(m2) + "/" + b2);
            w.write(geometry_table(x));
            const BasisSet bs = build_basis(x, b2, p.input.basis_directory);
            w.write(basis_block(bs, psp->notes));
            st.reset();
            st.method = m2;
            st.context = "Single point · ";
            ElectronicCallbacks cb2 = stream_callbacks(st);
            const SinglePoint sp = run_single_point(x, bs, m2, eo, b2 == p.input.basis ? guess : nullptr, cb2);
            write_point_details(w, sp, bs, true);
            r.single_point = sp;
            r.final_energy = sp.energy;
            r.final_energy_label = method_label(m2) + "/" + b2;
            if (!sp.converged) throw std::runtime_error("the single-point calculation did not converge");
            if (r.thermo) {
                // the final energies: single point + thermal corrections of the frequency step
                ThermoResult& t = *r.thermo;
                t.electronic_energy = sp.energy;
                t.total_zpe = sp.energy + joule_per_mol_to_hartree(t.zero_point_energy);
                t.total_u = sp.energy + joule_per_mol_to_hartree(t.thermal_energy);
                t.total_h = sp.energy + joule_per_mol_to_hartree(t.enthalpy);
                t.total_g = sp.energy + joule_per_mol_to_hartree(t.gibbs);
                t.total_a = sp.energy + joule_per_mol_to_hartree(t.helmholtz);
                w.line();
                w.line("    Final energies: " + r.final_energy_label + " single point + thermal corrections of " +
                       method_label(method) + "/" + p.input.basis + strprintf(" (%.2f K, %.0f Pa)", t.temperature, t.pressure));
                w.line(strprintf("      E(%s)%*s %20.10f E_h", r.final_energy_label.c_str(),
                                 std::max(0, 24 - static_cast<int>(r.final_energy_label.size())), "", sp.energy));
                w.line(strprintf("      E + ZPE                      %20.10f E_h", t.total_zpe));
                w.line(strprintf("      U = E + thermal energy       %20.10f E_h", t.total_u));
                w.line(strprintf("      H = E + thermal enthalpy     %20.10f E_h", t.total_h));
                w.line(strprintf("      A = E + thermal Helmholtz    %20.10f E_h", t.total_a));
                w.line(strprintf("      G = E + thermal Gibbs        %20.10f E_h", t.total_g));
            }
        }
        r.completed = true;
        r.termination = "normal";
    } catch (const Cancelled&) {
        r.termination = "cancelled";
        r.seconds = seconds_since(t0);
        write_termination(w, r, r.seconds);
        r.output = w.text();
        throw;
    } catch (const std::exception& e) {
        r.termination = "error";
        if (r.message.empty()) r.message = e.what();
    }
    // ---- summary
    w.section("Summary");
    const MolecularResult fp = final_point_result(r, p.input);
    w.write(summary_lines(fp));
    w.write(final_summary(r, p));
    for (const auto& x : r.warnings) w.line("    Warning: " + x);
    r.seconds = seconds_since(t0);
    write_termination(w, r, r.seconds);
    r.output = w.text();
    return r;
}

}  // namespace macqc
