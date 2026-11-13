// Tests of the coupled-cluster methods (CTest):
//  - the general-order spin-orbital engine (CCSD, CCSDT, CCSDTQ) against brute-force
//    coupled cluster in the determinant space (backend/tests/data/cc_reference.json, from
//    gen_cc_reference.py): random integrals with a non-diagonal Fock matrix, H2O/STO-3G
//    and LiH/6-31G (CCSDTQ = Full CI for four electrons);
//  - the closed-shell CCSD and (T) against the same data and PySCF;
//  - molecules through the driver and the JSON API.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/calculation.h"
#include "macqc/cc.h"
#include "macqc/macqc_api.h"
#include "macqc/molecular.h"
#include "macqc/mo_integrals.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

#define CHECK(cond)                                                     \
    do {                                                                \
        ++g_checks;                                                     \
        if (!(cond)) {                                                  \
            ++g_failures;                                               \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                               \
    } while (0)

static void check_near(const std::string& name, double got, double want, double tol) {
    ++g_checks;
    const bool ok = std::fabs(got - want) <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-40s got %.12f expected %.12f (diff %.1e)\n", ok ? "ok  " : "FAIL", name.c_str(), got, want,
                got - want);
}

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

struct Case {
    std::string name;
    int nocc = 0;
    double enuc = 0.0;
    MOIntegrals mo;
    json energies, pyscf;
};

static std::vector<Case> load_cases() {
    const json data = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/cc_reference.json"));
    std::vector<Case> out;
    for (const json& c : data["cases"]) {
        Case k;
        k.name = c["name"].get<std::string>();
        k.nocc = c["nocc"].get<int>();
        k.enuc = c["nuclear_repulsion"].get<double>();
        const int n = c["n"].get<int>();
        k.mo.n = n;
        k.mo.h = Matrix(static_cast<std::size_t>(n), static_cast<std::size_t>(n));
        const auto h = c["h"].get<std::vector<double>>();
        for (int p = 0; p < n; ++p)
            for (int q = 0; q < n; ++q) k.mo.h(static_cast<std::size_t>(p), static_cast<std::size_t>(q)) = h[static_cast<std::size_t>(p * n + q)];
        k.mo.npair = static_cast<std::size_t>(n) * (n + 1) / 2;
        k.mo.stride = k.mo.npair;
        k.mo.eri.assign(k.mo.npair * k.mo.npair, 0.0);
        const auto packed = c["eri_packed"].get<std::vector<double>>();
        std::size_t pos = 0;
        for (std::size_t pq = 0; pq < k.mo.npair; ++pq)
            for (std::size_t rs = 0; rs <= pq; ++rs) {
                const double v = packed.at(pos++);
                k.mo.eri[pq * k.mo.stride + rs] = v;
                k.mo.eri[rs * k.mo.stride + pq] = v;
            }
        k.energies = c["energies"];
        if (c.contains("pyscf")) k.pyscf = c["pyscf"];
        out.push_back(std::move(k));
    }
    return out;
}

static void test_spin_orbital_engine(const std::vector<Case>& cases) {
    std::printf("-- general-order spin-orbital engine vs brute force\n");
    CCOptions opt;
    opt.energy_tol = 1e-12;
    opt.amplitude_tol = 1e-9;
    for (const Case& c : cases) {
        const char* names[] = {"", "", "CCSD", "CCSDT", "CCSDTQ"};
        for (int rank = 2; rank <= 4; ++rank) {
            if (!c.energies.contains(names[rank])) continue;
            const CCResult r = run_cc_spin_orbital(c.mo, c.nocc, rank, c.enuc, opt);
            CHECK(r.converged);
            check_near(c.name + " " + names[rank] + " (spin-orbital)", r.correlation_energy,
                       c.energies[names[rank]].get<double>(), 1e-9);
        }
    }
}

static void test_closed_shell(const std::vector<Case>& cases) {
    std::printf("-- closed-shell CCSD and (T) (MO ladder)\n");
    CCOptions opt;
    opt.energy_tol = 1e-12;
    opt.amplitude_tol = 1e-9;
    for (const Case& c : cases) {
        const CCIntegrals ints = cc_integrals_from_mo(c.mo, c.nocc);
        auto ladder = make_mo_ladder(c.mo, c.nocc);
        const double eref = reference_energy_from_mo(c.mo, c.nocc, c.enuc);
        const CCResult r = run_ccsd_closed_shell(ints, *ladder, eref, true, opt);
        CHECK(r.converged);
        check_near(c.name + " CCSD (closed shell)", r.correlation_energy, c.energies["CCSD"].get<double>(), 1e-9);
        if (!c.pyscf.is_null())
            check_near(c.name + " CCSD(T) vs PySCF", r.total_energy - r.reference_energy, c.pyscf["CCSD(T)"].get<double>(),
                       1e-9);
    }
}

static MolecularResult run(const std::string& molecule, const std::string& basis, const std::string& method) {
    MolecularInput in;
    in.molecule = molecule;
    in.basis = basis;
    in.method = method;
    in.print_level = 1;
    MolecularResult r = prepare_molecular(in);
    if (!r.refusal.empty()) throw std::runtime_error(r.refusal);
    run_molecular(r);
    return r;
}

// CCSD and CCSD(T) through the driver (AO-driven ladder) against PySCF.
static void test_molecules_pyscf() {
    std::printf("-- CCSD and CCSD(T) against PySCF (cc_molecules_reference.json)\n");
    const json data = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/cc_molecules_reference.json"));
    for (const json& c : data["cases"]) {
        const std::string mol = c["molecule"], basis = c["basis"];
        const MolecularResult r = run(mol, basis, "CCSD(T)");
        CHECK(r.cc && r.cc->converged && r.cc->has_triples);
        check_near(mol + "/" + basis + " RHF", r.rhf->total_energy, c["rhf_energy"].get<double>(), 2e-9);
        check_near(mol + "/" + basis + " CCSD", r.cc->cc_energy, c["ccsd_energy"].get<double>(), 2e-9);
        check_near(mol + "/" + basis + " CCSD(T)", r.cc->total_energy, c["ccsd_t_energy"].get<double>(), 2e-9);
    }
}

// CCSDT and CCSDTQ through the driver: exact for two and four electrons.
static void test_molecules_exact() {
    std::printf("-- CCSDT/CCSDTQ equal Full CI for H2 and LiH (fci_reference.json)\n");
    const json data = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/fci_reference.json"));
    auto fci = [&](const std::string& m, const std::string& b) {
        for (const json& c : data["cases"])
            if (c["molecule"] == m && c["basis"] == b) return c["fci_energy"].get<double>();
        throw std::runtime_error("no Full CI reference for " + m + "/" + b);
    };
    for (const char* method : {"CCSD", "CCSDT", "CCSDTQ"}) {
        const MolecularResult r = run("H2", "aug-cc-pVDZ", method);
        CHECK(r.cc && r.cc->converged);
        check_near(std::string("H2/aug-cc-pVDZ ") + method + " = FCI", r.cc->total_energy, fci("H2", "aug-cc-pVDZ"), 2e-9);
    }
    const MolecularResult t = run("LiH", "aug-cc-pVDZ", "CCSDT");
    CHECK(t.cc && t.cc->converged);
    const double e_fci = fci("LiH", "aug-cc-pVDZ");
    // CCSDT lies within 0.05 mE_h of Full CI for LiH; CCSDTQ is exact.
    CHECK(std::fabs(t.cc->total_energy - e_fci) < 5e-5 && std::fabs(t.cc->total_energy - e_fci) > 1e-9);
    std::printf("     LiH/aug-cc-pVDZ CCSDT - FCI = %.3e E_h\n", t.cc->total_energy - e_fci);
    const MolecularResult q = run("LiH", "aug-cc-pVDZ", "CCSDTQ");
    CHECK(q.cc && q.cc->converged);
    check_near("LiH/aug-cc-pVDZ CCSDTQ = FCI", q.cc->total_energy, e_fci, 2e-9);
    CHECK(q.report.find("CCSDTQ energies") != std::string::npos);
}

// Estimates, refusals, notes and the JSON API.
static void test_api() {
    std::printf("-- estimates, refusals and the JSON API\n");
    {
        MolecularInput in;
        in.molecule = "H2O";
        in.basis = "aug-cc-pVTZ";
        in.method = "CCSDTQ";
        const MolecularResult r = prepare_molecular(in);
        CHECK(!r.refusal.empty() && r.refusal.find("CCSDTQ is not possible") != std::string::npos);
        CHECK(r.estimate.cc && !r.estimate.cc->feasible);
        in.method = "ccsd(t)";
        const MolecularResult s = prepare_molecular(in);
        CHECK(s.refusal.empty() && s.input.method == "CCSD(T)" && s.estimate.cc && s.estimate.cc->feasible);
        in.molecule = "H2";
        in.method = "CCSDT";
        const MolecularResult h = prepare_molecular(in);
        CHECK(!h.notes.empty() && h.notes.back().find("two electrons") != std::string::npos);
        in.method = "MP2";
        CHECK(validate_molecular(in).find("CCSDTQ") != std::string::npos);
    }
    char* out = macqc_run(R"json({"molecule":"LiH","basis":"aug-cc-pVDZ","method":"CCSD(T)"})json");
    const json j = json::parse(out);
    macqc_free(out);
    CHECK(j["error"].is_null());
    CHECK(j["coupledCluster"].is_object());
    CHECK(j["coupledCluster"]["method"] == "CCSD(T)");
    CHECK(j["coupledCluster"]["triplesCorrection"].is_number());
    CHECK(j["estimate"]["coupledCluster"]["feasible"] == true);
    CHECK(j["fullCI"].is_null());
    check_near("LiH/aug-cc-pVDZ CCSD(T) via JSON", j["coupledCluster"]["totalEnergy"].get<double>(), -8.0212069372, 2e-9);
    char* est = macqc_run(R"json({"molecule":"H2O","basis":"aug-cc-pVDZ","method":"CCSDT","estimateOnly":true})json");
    const json e = json::parse(est);
    macqc_free(est);
    CHECK(e["estimate"]["coupledCluster"]["engine"] == "spin-orbital, general order");
    CHECK(e["coupledCluster"].is_null() && e["rhf"].is_null());
}

int main() {
    try {
        const std::vector<Case> cases = load_cases();
        test_spin_orbital_engine(cases);
        test_closed_shell(cases);
        test_molecules_pyscf();
        test_molecules_exact();
        test_api();
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
