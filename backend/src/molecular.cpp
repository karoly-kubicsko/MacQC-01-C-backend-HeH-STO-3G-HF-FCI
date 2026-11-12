#include "macqc/molecular.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "macqc/gto_integrals.h"

namespace macqc {

namespace {

std::string format(const char* fmt, ...) {
    char buf[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    return buf;
}

inline int nc(int l) { return (l + 1) * (l + 2) / 2; }

// Work units of one shell quartet class (la >= lb, lc >= ld as in the engine).
struct ClassCost {
    double per_primitive = 0.0;  // vertical recurrence + accumulation, per primitive quartet
    double per_quartet = 0.0;    // horizontal recurrences + spherical transformation
};

ClassCost class_cost(int la, int lb, int lc, int ld) {
    const int Lab = la + lb, Lcd = lc + ld, L = Lab + Lcd;
    double vrr = 0.0;
    for (int E = 1; E <= Lab; ++E) vrr += nc(E) * (L - E + 1);
    for (int F = 1; F <= Lcd; ++F) {
        const int emin = std::max(0, la - (Lcd - F));
        double se = 0.0;
        for (int E = emin; E <= Lab; ++E) se += nc(E);
        vrr += (Lcd - F + 1) * nc(F) * se;
    }
    double nk = 0.0;
    for (int F = lc; F <= Lcd; ++F) nk += nc(F);
    double acc = 0.0;
    for (int E = la; E <= Lab; ++E) acc += nc(E) * nk;
    double post = 0.0;
    for (int j = 1; j <= lb; ++j)
        for (int e = la; e <= la + lb - j; ++e) post += nc(e) * nc(j) * nk;
    post += 2.0 * nc(la) * nc(lb) * nk;
    const double nbra = (2 * la + 1) * (2 * lb + 1);
    for (int j = 1; j <= ld; ++j)
        for (int f = lc; f <= lc + ld - j; ++f) post += nc(f) * nc(j) * nbra;
    post += 2.0 * nc(lc) * nc(ld) * nbra;
    ClassCost c;
    c.per_primitive = L == 0 ? 8.0 : vrr + acc + 40.0;
    c.per_quartet = post + 150.0;
    return c;
}

// Seconds per work unit and per stored integral in a Fock build (one thread), calibrated
// with H2O and LiH in aug-cc-pV5Z on a 2.1 GHz x86-64 core (Apple silicon cores are faster).
constexpr double kSecondsPerUnit = 1.14e-9;
constexpr double kSecondsPerStoredIntegral = 3.5e-9;

std::string angstrom_bohr(double bohr) {
    return format("%.6f \xC3\x85 (%.6f bohr)", bohr * kBohrInAngstrom, bohr);
}

std::string human_seconds(double s) {
    if (s < 1.0) return "less than a second";
    if (s < 90.0) return format("about %.0f s", s);
    if (s < 5400.0) return format("about %.0f min", s / 60.0);
    return format("about %.1f h", s / 3600.0);
}

std::string thousands(std::size_t v) {
    std::string s = std::to_string(v), out;
    int k = 0;
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        if (k && k % 3 == 0) out += ',';
        out += *it;
        ++k;
    }
    return std::string(out.rbegin(), out.rend());
}

}  // namespace

std::string validate_molecular(const MolecularInput& in) {
    const std::string m = canonical_molecule_name(in.molecule);
    if (m.empty() || m == "HeH+") return "molecule must be H2, LiH or H2O";
    if (canonical_basis_name(in.basis).empty()) {
        std::string list;
        for (const auto& b : supported_basis_sets()) list += (list.empty() ? "" : ", ") + b;
        return "unknown basis set '" + in.basis + "' (available: " + list + ")";
    }
    std::string method;
    for (char c : in.method) method += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (method != "RHF" && method != "HF")
        return "method must be RHF (Full CI for H2, LiH and H2O follows in a later version)";
    if (!std::isfinite(in.R) || in.R < 0.0) return "R must be a positive, finite number (bohr)";
    if (in.R > 0.0 && in.R < 0.3) return "R must be at least 0.3 bohr (0.16 \xC3\x85)";
    if (in.R > 50.0) return "R must be at most 50 bohr (26.5 \xC3\x85)";
    if (m == "H2O" && !(std::isfinite(in.angle) && in.angle >= 30.0 && in.angle <= 180.0))
        return "the H-O-H angle must be between 30 and 180 degrees";
    if (in.threads < 0 || in.threads > 64) return "threads must be between 0 (automatic) and 64";
    if (!std::isfinite(in.memory_mb) || in.memory_mb < 0.0) return "memoryMB must be >= 0";
    if (in.max_iter < 1 || in.max_iter > 1000) return "scfMaxIter must be between 1 and 1000";
    if (!(in.energy_tol > 0.0) || !(in.gradient_tol > 0.0)) return "SCF convergence thresholds must be positive";
    if (in.print_level < 0 || in.print_level > 2) return "printLevel must be 0, 1 or 2";
    return "";
}

int expected_scf_iterations(const std::string& molecule) {
    // Observed with the core-Hamiltonian guess and DIIS near equilibrium (7, 10, 13), plus one.
    if (molecule == "H2") return 8;
    if (molecule == "LiH") return 11;
    return 14;
}

double estimate_rhf_seconds(const BasisSet& basis, bool eri_stored, int threads, int iterations) {
    const ShellPairList pairs = make_shell_pairs(basis);
    // Group the pairs by (la, lb, number of primitive pairs) to keep the double loop cheap.
    struct Group {
        int la, lb;
        double count = 0.0, prims = 0.0;
    };
    std::vector<Group> groups;
    for (const ShellPair& p : pairs.pairs) {
        const int la = basis.shells[p.a].l, lb = basis.shells[p.b].l;
        const double k = static_cast<double>(p.prims.size());
        if (k == 0.0) continue;
        auto it = std::find_if(groups.begin(), groups.end(),
                               [&](const Group& g) { return g.la == la && g.lb == lb; });
        if (it == groups.end()) {
            groups.push_back({la, lb, 0.0, 0.0});
            it = groups.end() - 1;
        }
        it->count += 1.0;
        it->prims += k;
    }
    // Sum over unordered pairs of pairs: ~ half of all ordered combinations.
    double units = 0.0;
    for (const Group& g : groups)
        for (const Group& h : groups) {
            const ClassCost c = class_cost(g.la, g.lb, h.la, h.lb);
            units += 0.5 * (g.prims * h.prims * c.per_primitive + g.count * h.count * c.per_quartet);
        }
    const double t = std::max(1, threads);
    const double one_pass = units * kSecondsPerUnit / t;
    if (eri_stored) {
        const double n = static_cast<double>(unique_eri_count(static_cast<std::size_t>(basis.n_functions)));
        return one_pass + iterations * n * kSecondsPerStoredIntegral / t;
    }
    return iterations * one_pass;
}

MolecularResult prepare_molecular(const MolecularInput& input) {
    const std::string problem = validate_molecular(input);
    if (!problem.empty()) throw std::runtime_error(problem);
    MolecularResult r;
    r.input = input;
    r.input.molecule = canonical_molecule_name(input.molecule);
    r.input.basis = canonical_basis_name(input.basis);
    r.input.method = "RHF";
    if (!(r.input.R > 0.0)) r.input.R = default_bond_length(r.input.molecule);
    r.molecule = make_molecule(r.input.molecule, r.input.R, r.input.angle);
    r.basis = build_basis(r.molecule, r.input.basis, r.input.basis_directory);
    for (const auto& a : r.basis.atoms)
        if (!a.note.empty() && std::find(r.notes.begin(), r.notes.end(), a.note) == r.notes.end())
            r.notes.push_back(a.note);

    CostEstimate& e = r.estimate;
    e.basis_functions = r.basis.n_functions;
    e.shells = static_cast<int>(r.basis.shells.size());
    e.max_l = r.basis.max_l;
    e.unique_eris = unique_eri_count(static_cast<std::size_t>(r.basis.n_functions));
    e.eri_memory_mb = static_cast<double>(e.unique_eris) * 8.0 / (1024.0 * 1024.0);
    e.eri_stored = r.input.memory_mb > 0.0 && e.eri_memory_mb <= r.input.memory_mb;
    e.threads = resolve_threads(r.input.threads);
    e.rhf_seconds =
        estimate_rhf_seconds(r.basis, e.eri_stored, e.threads, expected_scf_iterations(r.molecule.name));
    if (e.rhf_seconds > kLongCalculationSeconds)
        r.warnings.push_back(format(
            "Long calculation: %d basis functions (up to %c functions); %s two-electron integrals %s. "
            "Estimated time %s on %d threads.",
            e.basis_functions, angular_letter(e.max_l), thousands(e.unique_eris).c_str(),
            e.eri_stored ? format("kept in memory (%.0f MB)", e.eri_memory_mb).c_str()
                         : format("recomputed in every SCF iteration (storing them would need %.1f GB)",
                                  e.eri_memory_mb / 1024.0)
                               .c_str(),
            human_seconds(e.rhf_seconds).c_str(), e.threads));
    return r;
}

void run_molecular(MolecularResult& r, const RHFProgress& progress) {
    RHFOptions opt;
    opt.max_iter = r.input.max_iter;
    opt.energy_tol = r.input.energy_tol;
    opt.gradient_tol = r.input.gradient_tol;
    opt.threads = r.input.threads;
    opt.memory_limit_mb = r.input.memory_mb;
    r.rhf = run_rhf(r.molecule, r.basis, opt, progress);
    if (r.rhf->removed_functions > 0)
        r.warnings.push_back(format("%d near-linearly dependent combination(s) of basis functions removed "
                                    "(overlap eigenvalues below %.0e).",
                                    r.rhf->removed_functions, opt.lindep_tol));
    if (!r.rhf->converged)
        r.warnings.push_back(format("The SCF did not converge in %d iterations.", r.rhf->iterations));
    r.report = molecular_report(r);
}

std::string molecular_report(const MolecularResult& r) {
    const int level = r.input.print_level;
    if (level <= 0) return "";
    const std::string rule(78, '=');
    const Molecule& mol = r.molecule;
    std::string s;
    s += "\xCE\xA8Lab \xE2\x80\x93 restricted Hartree-Fock\n" + rule + "\n";
    s += format("Molecule        %s, %d electrons, closed shell (charge %d, singlet)\n", mol.name.c_str(),
                mol.electrons(), mol.charge);
    if (mol.name == "H2O")
        s += "Geometry        R(O-H) = " + angstrom_bohr(r.input.R) + format(", angle H-O-H = %.4f deg\n", r.input.angle);
    else if (mol.name == "LiH")
        s += "Geometry        R(Li-H) = " + angstrom_bohr(r.input.R) + "\n";
    else
        s += "Geometry        R(H-H) = " + angstrom_bohr(r.input.R) + "\n";
    s += "Cartesian coordinates (bohr)\n";
    for (const Atom& a : mol.atoms)
        s += format("    %-2s  %16.10f  %16.10f  %16.10f\n", a.symbol.c_str(), a.position[0], a.position[1],
                    a.position[2]);
    s += format("Nuclear repulsion energy  %20.10f E_h\n\n", mol.nuclear_repulsion());

    const BasisSet& b = r.basis;
    s += "Basis set       " + b.name + " (spherical harmonics; basissetexchange.org, Gaussian94 format)\n";
    std::vector<std::string> seen;
    for (const auto& a : b.atoms) {
        if (std::find(seen.begin(), seen.end(), a.symbol) != seen.end()) continue;
        seen.push_back(a.symbol);
        s += format("    %-2s  %-13s %-16s -> %-16s %4d functions\n", a.symbol.c_str(), a.basis.c_str(),
                    a.primitives.c_str(), a.contraction.c_str(), a.functions);
    }
    for (const auto& n : r.notes) s += "    Note: " + n + "\n";
    s += format("    %d basis functions, %zu shells, highest angular momentum %c (l = %d)\n", b.n_functions,
                b.shells.size(), angular_letter(b.max_l), b.max_l);
    const CostEstimate& e = r.estimate;
    s += format("    two-electron integrals: %s unique, ", thousands(e.unique_eris).c_str());
    s += e.eri_stored ? format("kept in memory (%.1f MB)\n", e.eri_memory_mb)
                      : format("integral-direct (storing them would need %.1f GB)\n", e.eri_memory_mb / 1024.0);
    for (const auto& w : r.warnings) s += "Warning: " + w + "\n";
    if (!r.rhf) return s;

    const RHFResult& h = *r.rhf;
    s += format("    overlap matrix: smallest eigenvalue %.3e; %s\n", h.smallest_overlap_eigenvalue,
                h.removed_functions ? format("%d functions removed", h.removed_functions).c_str()
                                    : "no functions removed");
    s += "\nSCF             core-Hamiltonian guess, DIIS (up to 8 Fock matrices)\n";
    s += format("    converged when |dE| < %.1e E_h and max |X'(FPS - SPF)X| < %.1e\n", r.input.energy_tol,
                r.input.gradient_tol);
    if (level >= 2) {
        s += "    iter        total energy (E_h)        dE (E_h)      gradient\n";
        for (const auto& it : h.history)
            s += format("    %4d  %24.12f  %14.3e  %12.3e%s\n", it.iteration, it.energy, it.delta_e, it.gradient,
                        it.diis ? "  DIIS" : "");
    }
    s += h.converged ? format("    SCF converged in %d iterations\n", h.iterations)
                     : format("    SCF NOT converged after %d iterations\n", h.iterations);
    s += format("    time: integrals %.2f s, SCF %.2f s, %d threads\n\n", h.integral_seconds, h.scf_seconds, h.threads);

    s += "Energies (E_h)\n";
    s += format("    total energy (RHF)        %20.10f\n", h.total_energy);
    s += format("    electronic energy         %20.10f\n", h.electronic_energy);
    s += format("    nuclear repulsion         %20.10f\n", h.nuclear_repulsion);
    s += format("    kinetic energy            %20.10f\n", h.kinetic_energy);
    s += format("    nuclear attraction        %20.10f\n", h.nuclear_attraction_energy);
    s += format("    electron repulsion        %20.10f\n", h.two_electron_energy);
    s += format("    virial ratio -V/T         %20.10f\n\n", h.virial_ratio);

    const int nocc = h.n_occupied;
    const int nshow = level >= 2 ? static_cast<int>(h.orbital_energies.size())
                                 : std::min(static_cast<int>(h.orbital_energies.size()), nocc + 5);
    s += format("Orbital energies (E_h), %d occupied (*)%s\n", nocc,
                nshow < static_cast<int>(h.orbital_energies.size()) ? ", lowest virtual orbitals" : "");
    for (int k = 0; k < nshow; ++k) {
        s += format("  %4d%s %15.8f", k + 1, k < nocc ? "*" : " ", h.orbital_energies[k]);
        s += (k % 4 == 3 || k + 1 == nshow) ? "\n" : "";
    }
    if (nocc >= 1 && nocc < static_cast<int>(h.orbital_energies.size()))
        s += format("    HOMO %.8f E_h, LUMO %.8f E_h, gap %.8f E_h\n", h.orbital_energies[nocc - 1],
                    h.orbital_energies[nocc], h.orbital_energies[nocc] - h.orbital_energies[nocc - 1]);
    return s;
}

}  // namespace macqc
