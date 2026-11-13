// Tests of the geometry optimiser of ΨLab 5.0 (CTest):
//  - internal coordinates: Wilson B matrix (bends, linear bends) against finite differences,
//    rank of G, back-transformation round trip;
//  - tight optimisations from ΨLab's default geometries against PySCF + geomeTRIC
//    (backend/tests/data/optimization_reference.json from reference/pyscf_optimization_reference.py):
//    RHF, CCSD, CCSD(T), Full CI for H2O, LiH, H2;
//  - difficult starts: asymmetric H2O (ends C2v), nearly linear and exactly linear H2O (leaves
//    the saddle point, keeps C2v);
//  - a numerical-gradient optimisation (CCSDT of H2, equal to Full CI) and the step limit.
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
#include "macqc/internal_coordinates.h"
#include "macqc/molecule.h"
#include "macqc/optimizer.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;
constexpr double kPi = 3.14159265358979323846;

static void check_max(const std::string& name, double worst, double tol) {
    ++g_checks;
    const bool ok = worst <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-70s %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
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

static Molecule water(double r1_ang, double r2_ang, double angle_deg) {
    const Molecule m = make_molecule("H2O", 1.8, 104.5);
    return geometry_from_internals(m, {r1_ang / kBohrInAngstrom, r2_ang / kBohrInAngstrom, angle_deg * kPi / 180.0});
}

struct Geo {
    double r1 = 0, r2 = 0, angle = 0;  // angstrom, degrees (diatomics: r1)
};

static Geo geo_of(const OptimizationResult& r) {
    Geo g;
    const auto& p = r.final_parameters;
    if (r.final_geometry.atoms.size() == 2) {
        g.r1 = p.at(0).value * kBohrInAngstrom;
        return g;
    }
    for (const auto& x : p) {
        if (x.angle) g.angle = x.value;
        else if (g.r1 == 0) g.r1 = x.value * kBohrInAngstrom;
        else g.r2 = x.value * kBohrInAngstrom;
    }
    return g;
}

static OptimizationResult optimize(const Molecule& start, const std::string& basis, const std::string& method,
                                   bool tight, int max_steps = 50) {
    const BasisSet b = build_basis(start, basis, default_basis_directory());
    ElectronicOptions eo;
    OptimizationOptions oo;
    oo.tight = tight;
    oo.max_steps = max_steps;
    return optimize_geometry(start, b, method, eo, oo);
}

static void test_internal_coordinates() {
    std::printf("\n== Internal coordinates\n");
    for (double angle : {104.5, 178.0}) {
        const Molecule m = water(0.97, 1.02, angle);
        const InternalCoordinates ic(m);
        const Matrix B = ic.wilson_b(m);
        double w = 0.0;
        for (std::size_t x = 0; x < 9; ++x) {
            Molecule mp = m, mm = m;
            mp.atoms[x / 3].position[x % 3] += 1e-5;
            mm.atoms[x / 3].position[x % 3] -= 1e-5;
            const auto qp = ic.values(mp), qm = ic.values(mm);
            for (std::size_t i = 0; i < ic.size(); ++i) w = std::max(w, std::fabs((qp[i] - qm[i]) / 2e-5 - B(i, x)));
        }
        Matrix G(ic.size(), ic.size());
        for (std::size_t i = 0; i < ic.size(); ++i)
            for (std::size_t j = 0; j < ic.size(); ++j)
                for (std::size_t k = 0; k < 9; ++k) G(i, j) += B(i, k) * B(j, k);
        int rank = 0;
        generalized_inverse(G, 1e-8, &rank);
        const bool linear = angle > 175.0;
        check_max(std::string("H2O at ") + (linear ? "178" : "104.5") + " deg: B matrix vs finite differences", w, 1e-7);
        CHECK(static_cast<int>(ic.size()) == (linear ? 4 : 3));
        CHECK(rank == (linear ? 4 : 3));  // 3N-5 for the (nearly) linear set, 3N-6 otherwise
        std::vector<double> dq(ic.size(), 0.0);
        dq[0] = 0.05;
        dq[1] = -0.03;
        dq[2] = linear ? 0.01 : 0.1;
        const BackTransform bt = back_transform(ic, m, dq);
        CHECK(bt.converged);
        check_max(std::string("back-transformation residual (") + (linear ? "linear bends" : "bend") + ")", bt.residual, 1e-9);
    }
}

static void test_reference() {
    std::printf("\n== Tight optimisations against PySCF + geomeTRIC\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/optimization_reference.json"));
    std::printf("     reference: %s\n", ref["program"].get<std::string>().c_str());
    for (const auto& c : ref["cases"]) {
        const std::string name = c["molecule"], basis = c["basis"], method = c["method"];
        const auto t0 = Clock::now();
        const Molecule start = make_molecule(name, default_bond_length(name), kDefaultWaterAngle);
        const OptimizationResult r = optimize(start, basis, method, true);
        const std::string label = name + " " + basis + " " + method;
        CHECK(r.converged);
        const Geo g = geo_of(r);
        if (name == "H2O") {
            check_max(label + ": R(O-H) / angstrom", std::max(std::fabs(g.r1 - c["R1"].get<double>()), std::fabs(g.r2 - c["R2"].get<double>())), 2e-5);
            check_max(label + ": angle / degree", std::fabs(g.angle - c["angle"].get<double>()), 3e-3);
            CHECK(r.symmetric_start && r.symmetric_end);
        } else {
            check_max(label + ": R / angstrom", std::fabs(g.r1 - c["R"].get<double>()), 2e-5);
        }
        check_max(label + ": energy / E_h", std::fabs(r.final_energy - c["energy"].get<double>()), 1e-8);
        std::printf("     %s: %d steps (%s), %.1f s\n", label.c_str(), r.steps, r.initial_hessian.c_str(), since(t0));
    }
}

static void test_difficult_starts() {
    std::printf("\n== Difficult starts (RHF/aug-cc-pVDZ)\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/optimization_reference.json"));
    Geo want;
    for (const auto& c : ref["cases"])
        if (c["molecule"] == "H2O" && c["basis"] == "aug-cc-pVDZ" && c["method"] == "RHF") {
            want.r1 = c["R1"].get<double>();
            want.angle = c["angle"].get<double>();
        }
    struct Start {
        const char* label;
        Molecule m;
        bool symmetric;
    };
    const Start starts[] = {{"asymmetric start (1.05, 0.90 angstrom, 120 deg)", water(1.05, 0.90, 120.0), false},
                            {"nearly linear start (179 deg)", water(0.96, 0.97, 179.0), false},
                            {"linear start (180 deg, saddle point)", water(0.96, 0.96, 180.0), true},
                            {"compressed start (0.8 angstrom, 70 deg)", water(0.8, 0.8, 70.0), true}};
    for (const Start& s : starts) {
        const OptimizationResult r = optimize(s.m, "aug-cc-pVDZ", "RHF", true);
        CHECK(r.converged);
        const Geo g = geo_of(r);
        check_max(std::string(s.label) + ": R", std::max(std::fabs(g.r1 - want.r1), std::fabs(g.r2 - want.r1)), 5e-5);
        check_max(std::string(s.label) + ": angle", std::fabs(g.angle - want.angle), 1e-2);
        if (s.symmetric) CHECK(r.symmetric_end);
        std::printf("     %d steps\n", r.steps);
    }
}

static void test_numerical_and_limits() {
    std::printf("\n== Numerical-gradient optimisation and step limit\n");
    const Molecule h2 = make_molecule("H2", default_bond_length("H2"), 0.0);
    const OptimizationResult t = optimize(h2, "aug-cc-pVDZ", "CCSDT", true);
    const OptimizationResult f = optimize(h2, "aug-cc-pVDZ", "FCI", true);
    CHECK(t.converged && f.converged);
    CHECK(!t.final_point.analytic && f.final_point.analytic);
    check_max("H2 CCSDT (numerical gradient) = Full CI (analytic) optimum, angstrom",
              std::fabs(geo_of(t).r1 - geo_of(f).r1), 2e-5);
    const OptimizationResult l = optimize(water(1.2, 1.2, 70.0), "aug-cc-pVDZ", "RHF", false, 2);
    CHECK(!l.converged);
    CHECK(l.steps == 2);
    CHECK(l.message.find("did not converge within 2 steps") != std::string::npos);
    std::printf("     step limit: %s\n", l.message.c_str());
}

int main() {
    try {
        const auto t0 = Clock::now();
        test_internal_coordinates();
        test_reference();
        test_difficult_starts();
        test_numerical_and_limits();
        std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t0));
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    return g_failures == 0 ? 0 : 1;
}
