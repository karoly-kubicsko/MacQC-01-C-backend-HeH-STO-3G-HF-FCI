#include "macqc/calculation.h"

#include <cmath>
#include <exception>

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
    m.max_iter = in.scfMaxIter;
    m.energy_tol = in.energyConvergence;
    m.gradient_tol = in.gradientConvergence;
    m.print_level = in.printLevel;
    return m;
}

Output run_molecular_calculation(const Input& in) {
    Output out;
    out.input = in;
    try {
        out.molecular = prepare_molecular(molecular_input(in));
        out.input.R = out.molecular->input.R;
        out.input.molecule = out.molecular->input.molecule;
        out.input.basis = out.molecular->input.basis;
        if (in.estimateOnly) {
            out.report = molecular_report(*out.molecular);
            return out;
        }
        run_molecular(*out.molecular);
        out.report = out.molecular->report;
        const RHFResult& h = *out.molecular->rhf;
        if (!h.converged)
            out.error = "SCF did not converge in " + std::to_string(h.iterations) + " iterations";
        else if (!std::isfinite(h.total_energy))
            out.error = "SCF produced non-finite energies";
    } catch (const std::exception& e) {
        out.error = e.what();
    } catch (...) {
        out.error = "unknown internal error";
    }
    return out;
}

}  // namespace

std::string validate(const Input& in) {
    if (canonical_molecule_name(in.molecule).empty())
        return "molecule must be HeH+, H2, LiH or H2O";
    if (!is_heh_input(in)) return validate_molecular(molecular_input(in));
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

Output run_calculation(const Input& in) {
    if (!canonical_molecule_name(in.molecule).empty() && !is_heh_input(in)) return run_molecular_calculation(in);
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
