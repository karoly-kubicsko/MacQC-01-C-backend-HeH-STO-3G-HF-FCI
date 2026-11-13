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

#include "macqc/cc.h"
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

int main() {
    try {
        const std::vector<Case> cases = load_cases();
        test_spin_orbital_engine(cases);
        test_closed_shell(cases);
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
