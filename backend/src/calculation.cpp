#include "macqc/calculation.h"

#include <cmath>
#include <exception>
#include <stdexcept>

#include "macqc/output.h"

namespace macqc {

bool is_heh_input(const Input& in) { return canonical_molecule_name(in.molecule) == "HeH+"; }

namespace {

MolecularInput molecular_input(const Input& in) {
    MolecularInput m;
    m.molecule = in.molecule;
    m.R = in.R;
    m.angle = in.angle;
    m.basis = in.basis;
    m.method = in.method;
    m.basis_directory = in.basisDirectory;
    m.threads = in.threads;
    m.memory_mb = in.memoryMB;
    m.fci_memory_mb = in.fciMemoryMB;
    m.max_iter = in.scfMaxIter;
    m.energy_tol = in.energyConvergence;
    m.gradient_tol = in.gradientConvergence;
    m.print_level = in.printLevel;
    return m;
}

JobSettings job_settings(const Input& in) {
    JobSettings s;
    s.type = canonical_job_type(in.job);
    s.sp_method = in.spMethod.empty() ? "" : canonical_method_name(in.spMethod);
    s.sp_basis = in.spBasis.empty() ? "" : canonical_basis_name(in.spBasis);
    s.temperature = in.temperature;
    s.pressure = in.pressure;
    s.tight = in.optConvergence == "tight";
    s.max_steps = in.maxSteps;
    return s;
}

std::string job_input_error(const Input& in) {
    if (canonical_job_type(in.job).empty()) return "job must be energy, opt, freq, opt+freq or opt+freq+sp";
    if (!in.spMethod.empty() && canonical_method_name(in.spMethod).empty())
        return "spMethod must be HF (RHF), CCSD, CCSD(T), CCSDT, CCSDTQ or FCI";
    if (!in.spBasis.empty() && canonical_basis_name(in.spBasis).empty()) return "unknown basis set '" + in.spBasis + "' for spBasis";
    if (!(std::isfinite(in.temperature) && in.temperature > 0.0 && in.temperature <= 1.0e5))
        return "temperature must be between 0 and 100000 K";
    if (!(std::isfinite(in.pressure) && in.pressure > 0.0 && in.pressure <= 1.0e12)) return "pressure must be positive (Pa)";
    if (in.optConvergence != "normal" && in.optConvergence != "tight") return "optConvergence must be normal or tight";
    if (in.maxSteps < 1 || in.maxSteps > 1000) return "maxSteps must be between 1 and 1000";
    return "";
}

Output run_molecular_calculation(const Input& in0, const OutputSink& sink) {
    Output out;
    out.input = in0;
    Input in = in0;
    try {
        const std::string je = job_input_error(in);
        if (!je.empty()) throw std::runtime_error(je);
        in.job = canonical_job_type(in.job);
        out.input.job = in.job;
        out.molecular = prepare_molecular(molecular_input(in));
        out.input.R = out.molecular->input.R;
        out.input.molecule = out.molecular->input.molecule;
        out.input.basis = out.molecular->input.basis;
        out.input.method = out.molecular->input.method;
        const JobSettings js = job_settings(in);
        if (in.estimateOnly) {
            out.report = molecular_report(*out.molecular);
            std::optional<MolecularResult> sp;
            if (job_has_single_point(js.type)) {
                MolecularInput m2 = molecular_input(in);
                if (!js.sp_method.empty()) m2.method = js.sp_method;
                if (!js.sp_basis.empty()) m2.basis = js.sp_basis;
                sp = prepare_molecular(m2);
            }
            JobResult jr;
            jr.settings = js;
            jr.estimate = estimate_job(*out.molecular, js, sp ? &*sp : nullptr);
            jr.termination = "estimate";
            out.job = jr;
            return out;
        }
        const std::string echo = in.echo.empty() ? input_echo(out.input) : in.echo;
        if (js.type != "energy") {
            const MolecularResult prepared = *out.molecular;
            JobResult jr = run_job(molecular_input(in), js, sink, echo);
            out.molecular = final_point_result(jr, molecular_input(in));
            out.molecular->estimate = prepared.estimate;
            out.molecular->refusal = jr.estimate.refusal;
            for (const auto& w : jr.warnings) out.molecular->warnings.push_back(w);
            out.report = out.molecular->report;
            if (jr.termination != "normal") out.error = jr.message.empty() ? "the job did not finish" : jr.message;
            out.job = std::move(jr);
            return out;
        }
        JobResult jr = run_energy_job(*out.molecular, js, sink, echo);
        if (!out.molecular->refusal.empty()) {
            out.error = out.molecular->refusal;
            out.job = std::move(jr);
            return out;
        }
        out.job = std::move(jr);
        if (!out.molecular->rhf) {
            out.error = out.job->message.empty() ? "no result" : out.job->message;
            return out;
        }
        out.report = out.molecular->report;
        const RHFResult& h = *out.molecular->rhf;
        if (!h.converged)
            out.error = "SCF did not converge in " + std::to_string(h.iterations) + " iterations";
        else if (!std::isfinite(h.total_energy))
            out.error = "SCF produced non-finite energies";
        else if (out.molecular->fci && !out.molecular->fci->converged)
            out.error = "Full CI did not converge in " + std::to_string(out.molecular->fci->iterations) + " iterations";
    } catch (const std::exception& e) {
        out.error = e.what();
    } catch (...) {
        out.error = "unknown internal error";
    }
    return out;
}

}  // namespace

std::string input_echo(const Input& in) {
    std::string s;
    auto kv = [&](const char* k, const std::string& v) { s += strprintf("    %-16s = %s\n", k, v.c_str()); };
    kv("molecule", in.molecule);
    kv("job", in.job);
    kv("method", in.method == "RHF" ? "HF" : in.method);
    kv("basis", in.basis);
    kv("R", strprintf("%.10f bohr (%.8f angstrom)", in.R, in.R * kBohrInAngstrom));
    if (canonical_molecule_name(in.molecule) == "H2O") kv("angle", strprintf("%.6f degrees", in.angle));
    if (job_has_single_point(canonical_job_type(in.job))) {
        kv("sp_method", in.spMethod.empty() ? "(same as method)" : in.spMethod);
        kv("sp_basis", in.spBasis.empty() ? "(same as basis)" : in.spBasis);
    }
    if (job_has_frequencies(canonical_job_type(in.job))) {
        kv("temperature", strprintf("%.2f K", in.temperature));
        kv("pressure", strprintf("%.2f Pa", in.pressure));
    }
    if (job_has_optimization(canonical_job_type(in.job))) {
        kv("opt_convergence", in.optConvergence);
        kv("max_steps", std::to_string(in.maxSteps));
    }
    kv("threads", in.threads == 0 ? "all cores" : std::to_string(in.threads));
    kv("memory", strprintf("%.0f MB (two-electron integrals in memory up to this size)", in.memoryMB));
    kv("fci_memory", strprintf("%.0f MB (limit of Full CI, coupled cluster and correlated gradients)", in.fciMemoryMB));
    return s;
}

std::string validate(const Input& in) {
    if (canonical_molecule_name(in.molecule).empty())
        return "molecule must be HeH+, H2, LiH or H2O";
    if (!is_heh_input(in)) {
        const std::string e = validate_molecular(molecular_input(in));
        return e.empty() ? job_input_error(in) : e;
    }
    if (canonical_job_type(in.job) != "energy") return "HeH+ (STO-nG) supports only the energy job";
    auto positive = [](double v) { return std::isfinite(v) && v > 0.0; };
    if (!positive(in.R)) return "R must be a positive, finite number (bohr)";
    if (!positive(in.zeta1)) return "zeta1 must be a positive, finite number";
    if (!positive(in.zeta2)) return "zeta2 must be a positive, finite number";
    if (!positive(in.ZA)) return "ZA must be a positive, finite number";
    if (!positive(in.ZB)) return "ZB must be a positive, finite number";
    if (in.nGaussians < 1 || in.nGaussians > 3) return "nGaussians must be 1, 2 or 3";
    if (in.maxIter < 1) return "maxIter must be at least 1";
    if (!positive(in.convergence)) return "convergence must be a positive, finite number";
    if (in.printLevel < 0 || in.printLevel > 2) return "printLevel must be 0, 1 or 2";
    return "";
}

Output run_calculation(const Input& in, const OutputSink& sink) {
    if (!canonical_molecule_name(in.molecule).empty() && !is_heh_input(in)) return run_molecular_calculation(in, sink);
    Output out;
    out.input = in;
    try {
        out.error = validate(in);
        if (!out.error.empty()) return out;

        DiatomicBasis basis;
        basis.n_gaussians = in.nGaussians;
        basis.r = in.R;
        basis.zeta1 = in.zeta1;
        basis.zeta2 = in.zeta2;
        basis.za = in.ZA;
        basis.zb = in.ZB;
        HFOptions opt;
        opt.exact_erf = in.exactErf;
        opt.max_iter = in.maxIter;
        opt.convergence = in.convergence;
        opt.print_level = in.printLevel;

        TwoElectronDiatomicHF hf(basis, opt);
        out.hf = hf.run();
        out.hf_report = out.hf->report;
        out.report = out.hf_report;

        const bool finite = std::isfinite(out.hf->total_energy);
        if (!out.hf->converged || !finite) {
            out.error = out.hf->converged ? "SCF produced non-finite energies"
                                          : "SCF did not converge in " + std::to_string(in.maxIter) +
                                                " iterations";
            return out;
        }

        out.fci = run_fci(out.hf->H, out.hf->eri, out.hf->C, out.hf->nuclear_repulsion,
                          out.hf->total_energy);
        out.report += fci_report(*out.fci, in.printLevel);
    } catch (const std::exception& e) {
        out.error = e.what();
    } catch (...) {
        out.error = "unknown internal error";
    }
    return out;
}

}  // namespace macqc
