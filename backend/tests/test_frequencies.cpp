// Tests of harmonic frequencies, IR intensities and thermochemistry of ΨLab 5.0 (CTest):
//  - RHF frequencies, reduced masses and force constants from the analytic Hessian against
//    PySCF (pyscf.hessian.thermo.harmonic_analysis) at the same geometries;
//  - thermochemistry against PySCF (pyscf.hessian.thermo.thermo, older constants) and against a
//    hand evaluation of the textbook formulas with CODATA 2018 and AME 2020
//    (backend/tests/data/frequency_reference.json, reference/pyscf_frequency_reference.py);
//  - numerical Hessians: finite differences of analytic gradients against the analytic RHF
//    Hessian; for CCSD(T) and Full CI against finite differences of energies; IR intensities
//    from relaxed dipoles against finite-field dipoles; relaxed dipoles against -dE/dF;
//  - imaginary frequencies left out of the thermochemistry; Eckart alignment.
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/constants.h"
#include "macqc/frequencies.h"
#include "macqc/gradient.h"
#include "macqc/gto_derivatives.h"
#include "macqc/hessian.h"
#include "macqc/optimizer.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

static void check_max(const std::string& name, double worst, double tol) {
    ++g_checks;
    const bool ok = worst <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-72s %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
}

#define CHECK(cond)                                                     \
    do {                                                                \
        ++g_checks;                                                     \
        if (!(cond)) {                                                  \
            ++g_failures;                                               \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                               \
    } while (0)

using Clock = std::chrono::steady_clock;
static double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static Molecule molecule_from_json(const json& atoms) {
    Molecule mol;
    for (const auto& a : atoms) {
        Atom at;
        at.symbol = a[0].get<std::string>();
        at.Z = atomic_number(at.symbol);
        at.position = {a[1][0].get<double>(), a[1][1].get<double>(), a[1][2].get<double>()};
        mol.atoms.push_back(at);
    }
    mol.name = mol.atoms.size() == 3 ? "H2O" : (mol.atoms[0].Z == 3 ? "LiH" : "H2");
    return mol;
}

static ElectronicOptions tight() {
    ElectronicOptions o;
    o.scf_energy_tol = 1e-12;
    o.scf_gradient_tol = 1e-10;
    o.cc_amplitude_tol = 1e-10;
    o.cc_energy_tol = 1e-12;
    o.fci_residual_tol = 1e-9;
    return o;
}

static double kj_to_h(double x) { return joule_per_mol_to_hartree(x); }

static void test_constants_and_alignment() {
    std::printf("\n== Constants, masses, Eckart alignment\n");
    check_max("AME 2020 mass of 16O", std::fabs(isotope_mass(8) - 15.994914619257), 1e-12);
    check_max("AME 2020 mass of 1H", std::fabs(isotope_mass(1) - 1.007825031898), 1e-12);
    check_max("1 E_h = 627.5094740631 kcal/mol (CODATA 2018)", std::fabs(codata2018::kHartreeKcalPerMol - 627.5094740631), 1e-7);
    // a translated and rotated copy is aligned back
    Molecule m = make_molecule("H2O", 1.81, 104.5);
    Molecule r = m;
    const double c = std::cos(0.7), s = std::sin(0.7);
    for (auto& a : r.atoms) {
        const auto p = a.position;
        a.position = {c * p[0] - s * p[1] + 0.3, s * p[0] + c * p[1] - 1.2, p[2] + 0.5};
    }
    const Molecule b = align_to(m, r);
    double w = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 3; ++k) w = std::max(w, std::fabs(b.atoms[i].position[k] - m.atoms[i].position[k]));
    check_max("align_to undoes a rotation and translation", w, 1e-12);
}

static void test_pyscf() {
    std::printf("\n== RHF frequencies and thermochemistry against PySCF and the hand evaluation\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/frequency_reference.json"));
    std::printf("     reference: %s\n", ref["program"].get<std::string>().c_str());
    for (const auto& c : ref["cases"]) {
        const std::string label = c["molecule"].get<std::string>() + " " + c["basis"].get<std::string>() + " RHF";
        const Molecule mol = molecule_from_json(c["atoms"]);
        const BasisSet b = build_basis(mol, c["basis"].get<std::string>(), default_basis_directory());
        const auto t0 = Clock::now();
        const HessianResult h = run_hessian(mol, b, "RHF", tight());
        CHECK(h.analytic);
        const HarmonicAnalysis a = harmonic_analysis(mol, h.hessian, &h.dipole_derivatives);
        const auto f = c["frequencies"].get<std::vector<double>>();
        const auto mu = c["reduced_masses"].get<std::vector<double>>();
        const auto k = c["force_constants"].get<std::vector<double>>();
        CHECK(a.modes.size() == f.size());
        double wf = 0, wm = 0, wk = 0;
        for (std::size_t i = 0; i < f.size() && i < a.modes.size(); ++i) {
            wf = std::max(wf, std::fabs(a.modes[i].wavenumber - f[i]));
            wm = std::max(wm, std::fabs(a.modes[i].reduced_mass - mu[i]));
            wk = std::max(wk, std::fabs(a.modes[i].force_constant - k[i]));
        }
        check_max(label + ": wavenumbers / cm-1", wf, 2e-3);
        check_max(label + ": reduced masses / u", wm, 1e-6);
        check_max(label + ": force constants / mdyn A-1", wk, 1e-5);
        const ThermoResult t = thermochemistry(mol, a, h.reference.point.energy, 298.15, 101325.0, symmetry_number(mol));
        CHECK(t.symmetry_number == c["symmetry_number"].get<int>());
        const json& hd = c["hand_thermo"];
        check_max(label + ": ZPE vs hand / E_h", std::fabs(kj_to_h(t.zero_point_energy) - hd["ZPE"].get<double>()), 1e-9);
        check_max(label + ": U, H, G, A corrections vs hand / E_h",
                  std::max({std::fabs(kj_to_h(t.thermal_energy) - hd["U"].get<double>()),
                            std::fabs(kj_to_h(t.enthalpy) - hd["H"].get<double>()),
                            std::fabs(kj_to_h(t.gibbs) - hd["G"].get<double>()),
                            std::fabs(kj_to_h(t.helmholtz) - hd["A"].get<double>())}),
                  1e-9);
        check_max(label + ": S, S_trans, S_rot, S_vib, Cv vs hand / J mol-1 K-1",
                  std::max({std::fabs(t.total.entropy - hd["S"].get<double>()),
                            std::fabs(t.translation.entropy - hd["S_trans"].get<double>()),
                            std::fabs(t.rotation.entropy - hd["S_rot"].get<double>()),
                            std::fabs(t.vibration.entropy - hd["S_vib"].get<double>()),
                            std::fabs(t.total.cv - hd["Cv"].get<double>())}),
                  5e-5);
        const json& p = c["pyscf_thermo"];
        const double jk = codata2018::kHartree * codata2018::kAvogadro;  // E_h/K -> J/(mol K)
        check_max(label + ": ZPE and totals vs PySCF / E_h",
                  std::max({std::fabs(kj_to_h(t.zero_point_energy) - p["ZPE"].get<double>()),
                            std::fabs(t.total_u - p["E_tot"].get<double>()),
                            std::fabs(t.total_h - p["H_tot"].get<double>()),
                            std::fabs(t.total_g - p["G_tot"].get<double>())}),
                  1e-5);
        check_max(label + ": S and Cv vs PySCF (older constants, masses) / J mol-1 K-1",
                  std::max(std::fabs(t.total.entropy - p["S_tot"].get<double>() * jk), std::fabs(t.total.cv - p["Cv_tot"].get<double>() * jk)),
                  2e-2);
        std::printf("     %s: %.1f s; S = %.4f J/(mol K)\n", label.c_str(), since(t0), t.total.entropy);
    }
}

static void test_numerical() {
    std::printf("\n== Numerical Hessians\n");
    const ElectronicOptions o = tight();
    // RHF: finite differences of analytic gradients against the analytic Hessian
    {
        const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/frequency_reference.json"));
        const Molecule mol = molecule_from_json(ref["cases"][0]["atoms"]);
        const BasisSet b = build_basis(mol, "aug-cc-pVDZ", default_basis_directory());
        const HessianResult ha = run_hessian(mol, b, "RHF", o);
        HessianOptions hn;
        hn.force_numerical = true;
        const HessianResult hg = run_hessian(mol, b, "RHF", o, hn);
        CHECK(!hg.analytic && hg.gradient_evaluations == 6);
        const HarmonicAnalysis A = harmonic_analysis(mol, ha.hessian, &ha.dipole_derivatives);
        const HarmonicAnalysis N = harmonic_analysis(mol, hg.hessian, &hg.dipole_derivatives);
        double wf = 0, wi = 0;
        for (std::size_t i = 0; i < A.modes.size(); ++i) {
            wf = std::max(wf, std::fabs(A.modes[i].wavenumber - N.modes[i].wavenumber));
            wi = std::max(wi, std::fabs(A.modes[i].ir_intensity - N.modes[i].ir_intensity));
        }
        check_max("H2O RHF: wavenumbers, FD of gradients vs analytic / cm-1", wf, 0.2);
        check_max("H2O RHF: IR intensities, FD of SCF dipoles vs analytic / km mol-1", wi, 0.01);
        CHECK(A.modes[0].symmetry == "A1" && A.modes[1].symmetry == "A1" && A.modes[2].symmetry == "B2");
        std::printf("     analytic: %.3f %.3f %.3f cm-1, %.3f %.3f %.3f km/mol\n", A.modes[0].wavenumber,
                    A.modes[1].wavenumber, A.modes[2].wavenumber, A.modes[0].ir_intensity, A.modes[1].ir_intensity,
                    A.modes[2].ir_intensity);
    }
    // CCSD(T) and Full CI: FD of analytic gradients against FD of energies (+ finite-field dipoles)
    struct Case {
        const char* name;
        const char* method;
        Molecule mol;
    };
    std::vector<Case> cases;
    cases.push_back({"H2O CCSD(T)", "CCSD(T)", make_molecule("H2O", 0.966 / kBohrInAngstrom, 104.0)});
    cases.push_back({"LiH FCI", "FCI", make_molecule("LiH", 1.585 / kBohrInAngstrom, 0.0)});
    for (const Case& c : cases) {
        const BasisSet b = build_basis(c.mol, "aug-cc-pVDZ", default_basis_directory());
        const auto t0 = Clock::now();
        const HessianResult hg = run_hessian(c.mol, b, c.method, o);
        HessianOptions he;
        he.from_energies = true;
        const HessianResult hE = run_hessian(c.mol, b, c.method, o, he);
        const HarmonicAnalysis G = harmonic_analysis(c.mol, hg.hessian, &hg.dipole_derivatives);
        const HarmonicAnalysis E = harmonic_analysis(c.mol, hE.hessian, &hE.dipole_derivatives);
        double wf = 0, wi = 0;
        for (std::size_t i = 0; i < G.modes.size(); ++i) {
            wf = std::max(wf, std::fabs(G.modes[i].wavenumber - E.modes[i].wavenumber));
            wi = std::max(wi, std::fabs(G.modes[i].ir_intensity - E.modes[i].ir_intensity));
        }
        check_max(std::string(c.name) + ": wavenumbers, FD of gradients vs FD of energies / cm-1", wf, 0.5);
        check_max(std::string(c.name) + ": IR, relaxed dipoles vs finite-field dipoles / km mol-1", wi, 0.05);
        // relaxed dipole at the reference geometry against -dE/dF (4-point formula: the
        // hyperpolarizability term beta F^2 / 6 of the 2-point formula is 1e-6 for LiH)
        double wd = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double F = 2e-4;
            auto E = [&](double f) {
                ElectronicOptions of = o;
                of.field[static_cast<std::size_t>(k)] = f;
                return run_single_point(c.mol, b, c.method, of).energy;
            };
            const double dEdF = (8.0 * (E(F) - E(-F)) - (E(2 * F) - E(-2 * F))) / (12.0 * F);
            wd = std::max(wd, std::fabs(hg.reference.dipole[static_cast<std::size_t>(k)] + dEdF));
        }
        check_max(std::string(c.name) + ": relaxed dipole vs -dE/dF / e bohr", wd, 2e-7);
        std::printf("     %s: %d gradients vs %d energies; %.1f s; %s\n", c.name, hg.gradient_evaluations,
                    hE.energy_evaluations, since(t0), G.modes.back().symmetry.c_str());
    }
}

static void test_imaginary() {
    std::printf("\n== Saddle point: imaginary frequency left out of the thermochemistry\n");
    const Molecule mol = make_molecule("H2O", 1.75, 180.0);
    const BasisSet b = build_basis(mol, "aug-cc-pVDZ", default_basis_directory());
    const HessianResult h = run_hessian(mol, b, "RHF", tight());
    const HarmonicAnalysis a = harmonic_analysis(mol, h.hessian, &h.dipole_derivatives);
    CHECK(a.linear);
    CHECK(a.n_vib == 4);
    int imag = 0;
    for (const auto& m : a.modes) imag += m.imaginary ? 1 : 0;
    CHECK(imag == 2);  // the degenerate bend of linear H2O
    const ThermoResult t = thermochemistry(mol, a, h.reference.point.energy, 298.15, 101325.0, 2);
    CHECK(t.excluded_modes == 2);
    CHECK(t.warnings.size() == 2);
    std::printf("     linear H2O: %.1fi cm-1 (twice), %.1f, %.1f cm-1; %s\n", -a.modes[0].wavenumber,
                a.modes[2].wavenumber, a.modes[3].wavenumber, t.warnings.empty() ? "" : t.warnings[0].c_str());
}

int main() {
    try {
        const auto t0 = Clock::now();
        test_constants_and_alignment();
        test_pyscf();
        test_numerical();
        test_imaginary();
        std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t0));
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    return g_failures == 0 ? 0 : 1;
}
