// Tests of the correlated gradients of ΨLab 5.0 (CTest):
//  - CCSD, CCSD(T) and Full CI analytic gradients (Lambda equations, relaxed densities,
//    Z-vector) against 4-point central finite differences of the energy, at distorted
//    geometries (H2O without symmetry);
//  - the same gradients against PySCF (pyscf.grad.ccsd, pyscf.grad.ccsd_t, CASCI with all
//    orbitals = Full CI): backend/tests/data/correlated_gradient_reference.json from
//    reference/pyscf_correlated_gradient_reference.py;
//  - the densities reproduce the energies; the Lagrangian is invariant to occupied-occupied
//    and virtual-virtual rotations; the Full CI orbital response vanishes;
//  - numerical gradients (central differences along symmetry-adapted internal coordinates)
//    against the analytic ones, and exactness (CCSDT = Full CI for two electrons).
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
#include "macqc/gradient.h"
#include "macqc/gto_derivatives.h"
#include "macqc/mo_integrals.h"
#include "macqc/molecule.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

static void check_max(const std::string& name, double worst, double tol) {
    ++g_checks;
    const bool ok = worst <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-66s max |diff| %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
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

static json g_ref;

static bool full_mode() {
    const char* e = std::getenv("PSILAB_TEST_MAX_NBF");
    return e && std::atoi(e) > 200;
}


static Molecule geometry(const std::string& name) {
    for (const auto& c : g_ref["cases"])
        if (c["geometry"] == name) return molecule_from_json(c["atoms"]);
    throw std::runtime_error("geometry " + name + " not in the reference file");
}

// Analytic gradient against 4-point finite differences of the energy for the listed coordinates.
static void fd_case(const std::string& geom, const std::string& basis_name, const std::string& method,
                    const std::vector<int>& coords) {
    const Molecule mol = geometry(geom);
    const BasisSet b = build_basis(mol, basis_name, default_basis_directory());
    const ElectronicOptions o = tight();
    const auto t0 = Clock::now();
    const GradientResult g = run_gradient(mol, b, method, o);
    const double ta = since(t0);
    CHECK(g.analytic);
    CHECK(g.point.converged);
    if (g.correlated) {
        check_max(geom + " " + method + ": densities reproduce the energy", std::fabs(g.correlated->density_energy - g.point.energy), 1e-9);
        check_max(geom + " " + method + ": Lagrangian invariant (occ-occ, vir-vir)", g.correlated->orbital_lagrangian_oo_vv, 1e-8);
        CHECK(g.correlated->zvector_converged);
    }
    double worst = 0.0;
    const double h = 2.0e-3;
    for (int x : coords) {
        auto E = [&](double d) {
            Molecule m = mol;
            m.atoms[static_cast<std::size_t>(x / 3)].position[static_cast<std::size_t>(x % 3)] += d;
            const SinglePoint p = run_single_point(m, basis_at_geometry(b, m), method, o, &g.point.rhf.C);
            return p.energy;
        };
        const double fd = (8.0 * (E(h) - E(-h)) - (E(2 * h) - E(-2 * h))) / (12.0 * h);
        worst = std::max(worst, std::fabs(fd - g.gradient[static_cast<std::size_t>(x)]));
    }
    check_max(geom + " " + basis_name + " " + method + ": analytic vs FD of energies", worst, 1e-8);
    std::printf("     (analytic gradient %.2f s%s)\n", ta,
                g.lambda_iterations ? (", " + std::to_string(g.lambda_iterations) + " Lambda iterations").c_str() : "");
}

static void test_fd() {
    std::printf("\n== Correlated analytic gradients against finite differences of the energy\n");
    fd_case("H2O-distorted", "aug-cc-pVDZ", "CCSD", {0, 1, 2, 3, 4, 5, 6, 7, 8});
    fd_case("H2O-distorted", "aug-cc-pVDZ", "CCSD(T)", {0, 4, 8});
    fd_case("LiH-3.10", "aug-cc-pVDZ", "CCSD", {5});
    fd_case("LiH-3.10", "aug-cc-pVDZ", "CCSD(T)", {5});
    fd_case("LiH-3.10", "aug-cc-pVDZ", "FCI", {5});
    fd_case("H2-1.30", "aug-cc-pVTZ", "FCI", {2});
    if (full_mode()) {
        fd_case("H2O-distorted", "aug-cc-pVTZ", "CCSD(T)", {0, 4, 8});
        fd_case("LiH-3.10", "aug-cc-pCVDZ", "FCI", {5});
    }
}

static void test_pyscf() {
    std::printf("\n== Correlated analytic gradients against PySCF\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/correlated_gradient_reference.json"));
    std::printf("     reference: %s\n", ref["program"].get<std::string>().c_str());
    for (const auto& c : ref["cases"]) {
        const std::string label = c["geometry"].get<std::string>() + " " + c["basis"].get<std::string>() + " " +
                                  c["method"].get<std::string>();
        const Molecule mol = molecule_from_json(c["atoms"]);
        const BasisSet b = build_basis(mol, c["basis"].get<std::string>(), default_basis_directory());
        const auto t0 = Clock::now();
        const GradientResult g = run_gradient(mol, b, c["method"].get<std::string>(), tight());
        CHECK(g.analytic);
        check_max(label + ": energy", std::fabs(g.point.energy - c["energy"].get<double>()), 1e-9);
        const auto rg = c["gradient"].get<std::vector<double>>();
        double w = 0.0;
        for (std::size_t x = 0; x < rg.size(); ++x) w = std::max(w, std::fabs(g.gradient[x] - rg[x]));
        // PySCF's Z-vector solvers stop at residuals of ~1e-8 (it warns), hence 5e-8. Its CASCI
        // gradient for H2/aug-cc-pVTZ differs from its own CCSD gradient (exact for two
        // electrons, equal energies) by 2.5e-7; there ΨLab's Full CI gradient is compared with
        // PySCF's CCSD gradient at 5e-8 and with the CASCI one at 5e-7.
        const bool two_electron_fci = c["method"] == "FCI" && c["geometry"] == "H2-1.30";
        check_max(label + ": gradient" + (two_electron_fci ? " (PySCF CASCI)" : ""), w, two_electron_fci ? 5e-7 : 5e-8);
        if (two_electron_fci)
            for (const auto& o : ref["cases"])
                if (o["geometry"] == c["geometry"] && o["basis"] == c["basis"] && o["method"] == "CCSD") {
                    const auto og = o["gradient"].get<std::vector<double>>();
                    double wc = 0.0;
                    for (std::size_t x = 0; x < og.size(); ++x) wc = std::max(wc, std::fabs(g.gradient[x] - og[x]));
                    check_max(label + ": gradient (PySCF CCSD, exact for 2 electrons)", wc, 5e-8);
                }
        std::printf("     %s: %.2f s\n", label.c_str(), since(t0));
    }
}

// The blocked (T) terms against PySCF's formulation on full o^3 v^3 arrays.
static void test_triples_blocks() {
    std::printf("\n== (T) Lambda sources and density terms: blocked implementation vs full arrays\n");
    for (const char* gname : {"H2O-distorted", "LiH-3.10"}) {
        const Molecule mol = geometry(gname);
        const BasisSet b = build_basis(mol, "aug-cc-pVDZ", default_basis_directory());
        ElectronicOptions o = tight();
        const SinglePoint sp = run_single_point(mol, b, "CCSD", o);
        const MOIntegrals mo = transform_to_mo(b, sp.rhf.H, sp.rhf.C, 0);
        const CCIntegrals ints = cc_integrals_from_mo(mo, sp.rhf.n_occupied);
        check_max(std::string(gname) + " aug-cc-pVDZ: blocked (T) terms = full-array terms",
                  compare_triples_implementations(ints, sp.cc->t1, sp.cc->t2), 1e-12);
    }
}

static void test_numerical() {
    std::printf("\n== Numerical gradients (internal coordinates) against analytic gradients\n");
    ElectronicOptions o = tight();
    o.numerical_step = 0.005;
    struct Case {
        std::string label, method, basis;
        Molecule mol;
        int energies;
    };
    std::vector<Case> cases;
    cases.push_back({"H2O distorted (3 internal coordinates)", "CCSD", "aug-cc-pVDZ", geometry("H2O-distorted"), 6});
    cases.push_back({"H2O C2v (2 symmetry-adapted coordinates)", "CCSD", "aug-cc-pVDZ",
                     make_molecule("H2O", 1.85, 101.0), 4});
    cases.push_back({"LiH (1 coordinate)", "CCSD(T)", "aug-cc-pVDZ", geometry("LiH-3.10"), 2});
    cases.push_back({"H2O C2v RHF", "RHF", "aug-cc-pVDZ", make_molecule("H2O", 1.75, 110.0), 4});
    for (const Case& c : cases) {
        const BasisSet b = build_basis(c.mol, c.basis, default_basis_directory());
        const GradientResult a = run_gradient(c.mol, b, c.method, o);
        const GradientResult n = run_gradient(c.mol, b, c.method, o, nullptr, ElectronicCallbacks(), true);
        CHECK(a.analytic && !n.analytic);
        CHECK(n.energy_evaluations == c.energies);
        double w = 0.0;
        for (std::size_t x = 0; x < a.gradient.size(); ++x) w = std::max(w, std::fabs(a.gradient[x] - n.gradient[x]));
        // central differences with step 0.005: error ~ h^2/6 E''' ~ 1e-5
        check_max(c.label + " " + c.method + ": numerical vs analytic", w, 3e-5);
        std::printf("     %s\n", n.description.c_str());
    }
    // Exactness: for two electrons CCSDT (numerical) equals Full CI (analytic).
    {
        const Molecule mol = geometry("H2-1.30");
        const BasisSet b = build_basis(mol, "aug-cc-pVDZ", default_basis_directory());
        const GradientResult t = run_gradient(mol, b, "CCSDT", o);
        const GradientResult f = run_gradient(mol, b, "FCI", o);
        CHECK(!t.analytic && f.analytic);
        CHECK(!t.numerical_reason.empty());
        check_max("H2 CCSDT (numerical) = Full CI (analytic) energy", std::fabs(t.point.energy - f.point.energy), 1e-9);
        double w = 0.0;
        for (std::size_t x = 0; x < t.gradient.size(); ++x) w = std::max(w, std::fabs(t.gradient[x] - f.gradient[x]));
        check_max("H2 CCSDT (numerical) = Full CI (analytic) gradient", w, 3e-5);
        std::printf("     CCSDT: %s\n     reason: %s\n", t.description.c_str(), t.numerical_reason.c_str());
    }
    // Removed near-linear dependencies: the correlated gradient falls back to numerical.
    {
        const Molecule mol = geometry("H2O-distorted");
        const BasisSet b = build_basis(mol, "aug-cc-pVDZ", default_basis_directory());
        ElectronicOptions ol = o;
        ol.lindep_tol = 1e-2;
        const GradientResult n = run_gradient(mol, b, "CCSD", ol);
        CHECK(!n.analytic);
        CHECK(n.numerical_reason.find("near-linear") != std::string::npos);
        CHECK(n.energy_evaluations == 6);
        std::printf("     removed functions: %s\n     (%s)\n", n.description.c_str(), n.numerical_reason.c_str());
    }
    // Internal coordinates: B matrix against finite differences, geometry_from_internals round trip.
    {
        const Molecule mol = geometry("H2O-distorted");
        const SimpleInternals si = simple_internals(mol);
        double wb = 0.0, wr = 0.0;
        for (std::size_t x = 0; x < 9; ++x) {
            Molecule mp = mol, mm = mol;
            mp.atoms[x / 3].position[x % 3] += 1e-5;
            mm.atoms[x / 3].position[x % 3] -= 1e-5;
            const SimpleInternals p = simple_internals(mp), m = simple_internals(mm);
            for (std::size_t k = 0; k < 3; ++k) wb = std::max(wb, std::fabs((p.values[k] - m.values[k]) / 2e-5 - si.B[k][x]));
        }
        std::vector<double> q = si.values;
        q[0] += 0.1;
        q[1] -= 0.05;
        q[2] += 0.2;
        const SimpleInternals r = simple_internals(geometry_from_internals(mol, q));
        for (std::size_t k = 0; k < 3; ++k) wr = std::max(wr, std::fabs(r.values[k] - q[k]));
        check_max("Wilson B matrix of R, R, angle vs finite differences", wb, 1e-8);
        check_max("geometry_from_internals reproduces the internal coordinates", wr, 1e-12);
    }
}

int main() {
    try {
        const auto t0 = Clock::now();
        g_ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/gradient_reference.json"));
        test_fd();
        test_pyscf();
        test_triples_blocks();
        test_numerical();
        std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t0));
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    return g_failures == 0 ? 0 : 1;
}
