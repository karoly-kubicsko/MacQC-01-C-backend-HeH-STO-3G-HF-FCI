// MacQC backend tests (run by CTest). Plain C++17, no test framework.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include "macqc/calculation.h"
#include "macqc/macqc_api.h"
#include "macqc/report.h"
#include "nlohmann/json.hpp"

using nlohmann::json;

static int g_checks = 0, g_failures = 0;

#define CHECK(cond)                                                                    \
    do {                                                                               \
        ++g_checks;                                                                    \
        if (!(cond)) {                                                                 \
            ++g_failures;                                                              \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);                \
        }                                                                              \
    } while (0)

static void check_near(const char* name, double got, double want, double tol) {
    ++g_checks;
    if (!(std::fabs(got - want) <= tol)) {
        ++g_failures;
        std::printf("FAIL %s: got %.12f, expected %.12f (tol %.1e)\n", name, got, want, tol);
    }
}

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        std::printf("FAIL cannot open %s\n", path.c_str());
        ++g_failures;
        return "<missing>";
    }
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static macqc::Input make_input(int n, double R, int level) {
    macqc::Input in;
    in.nGaussians = n;
    in.R = R;
    in.printLevel = level;
    return in;
}

// 7a: byte-for-byte comparison of the HF text with the reference script output.
static void compare_report(const char* file, int n, double R, int level) {
    const std::string expected = read_file(std::string(MACQC_TEST_DATA_DIR) + "/" + file);
    const macqc::Output out = macqc::run_calculation(make_input(n, R, level));
    ++g_checks;
    if (out.hf_report != expected) {
        ++g_failures;
        std::size_t k = 0;
        while (k < out.hf_report.size() && k < expected.size() && out.hf_report[k] == expected[k]) ++k;
        std::printf("FAIL %s: reports differ at byte %zu (got %zu bytes, expected %zu)\n", file, k,
                    out.hf_report.size(), expected.size());
    }
    // The full report must start with the HF part.
    CHECK(out.report.compare(0, expected.size(), expected) == 0);
    if (level >= 1) CHECK(out.report.size() > expected.size());  // FCI section appended
}

static void test_reports() {
    compare_report("hf_sto3g_level0.txt", 3, 1.4632, 0);
    compare_report("hf_sto3g_level1.txt", 3, 1.4632, 1);
    compare_report("hf_sto3g_level2.txt", 3, 1.4632, 2);
    compare_report("hf_sto2g_level2.txt", 2, 1.4632, 2);
    compare_report("hf_sto1g_level2.txt", 1, 1.4632, 2);
    compare_report("hf_sto3g_R2.5_level2.txt", 3, 2.5, 2);
}

static void test_number_formats() {
    CHECK(macqc::fortran_d(0.0, 18, 10) == "  0.0000000000D+00");
    CHECK(macqc::fortran_d(-4.227529132, 20, 12).find("-0.422752913200D+01") != std::string::npos);
    CHECK(macqc::fortran_d(1.0, 18, 10) == "  0.1000000000D+01");
    CHECK(macqc::fortran_d(9.99999999999, 18, 10) == "  0.1000000000D+02");
    CHECK(macqc::fortran_d(1.0e-12, 18, 10) == "  0.1000000000D-11");
}

// 7b: numerical reference values.
static void test_numbers() {
    const macqc::Output out = macqc::run_calculation(make_input(3, 1.4632, 2));
    CHECK(out.error.empty() && out.hf && out.fci);
    if (!out.hf || !out.fci) return;
    const auto& v = out.hf->integrals;
    const double t = 1e-6;
    check_near("S12", v.S12, 0.450770, t);
    check_near("T11", v.T11, 2.164313, t);
    check_near("T12", v.T12, 0.167013, t);
    check_near("T22", v.T22, 0.760033, t);
    check_near("V11A", v.V11A, -4.139827, t);
    check_near("V12A", v.V12A, -1.102912, t);
    check_near("V22A", v.V22A, -1.265246, t);
    check_near("V11B", v.V11B, -0.677230, t);
    check_near("V12B", v.V12B, -0.411305, t);
    check_near("V22B", v.V22B, -1.226615, t);
    check_near("V1111", v.V1111, 1.307152, t);
    check_near("V2111", v.V2111, 0.437279, t);
    check_near("V2121", v.V2121, 0.177267, t);
    check_near("V2211", v.V2211, 0.605703, t);
    check_near("V2221", v.V2221, 0.311795, t);
    check_near("V2222", v.V2222, 0.774608, t);

    check_near("E_el", out.hf->electronic_energy, -4.227529132, 1e-8);
    check_near("E_tot", out.hf->total_energy, -2.860661992, 1e-8);
    CHECK(out.hf->converged);
    CHECK(out.hf->iterations == 6);
    check_near("eps1", out.hf->orbital_energies[0], -1.597448, t);
    check_near("eps2", out.hf->orbital_energies[1], -0.061669, t);

    check_near("FCI total", out.fci->total_energy, -2.880711848, 1e-8);
    check_near("FCI corr", out.fci->correlation_energy, -0.020049857, 1e-8);
    CHECK(out.fci->states.size() == 4);
    const double want[4] = {-2.880711848, -2.073840420, -1.733990390, -0.471078740};
    if (out.fci->states.size() == 4) {
        for (int k = 0; k < 4; ++k) check_near("FCI state", out.fci->states[k].energy, want[k], t);
        // spin labels: three singlets and one triplet
        int singlets = 0, triplets = 0;
        for (const auto& s : out.fci->states) {
            singlets += s.spin == "singlet";
            triplets += s.spin == "triplet";
        }
        CHECK(singlets == 3 && triplets == 1);
        CHECK(out.fci->states[0].spin == "singlet");
        check_near("S2 ground", out.fci->states[0].S2, 0.0, 1e-10);
    }
    // Ground-state weights sum to one, HF configuration dominates.
    double sum = 0.0;
    for (const auto& w : out.fci->singlet_weights) sum += w.weight;
    check_near("weights sum", sum, 1.0, 1e-12);
    CHECK(out.fci->singlet_weights.size() == 3);
    if (out.fci->singlet_weights.size() == 3) CHECK(out.fci->singlet_weights[0].weight > 0.9);

    // The exact-erf option changes the result only slightly.
    macqc::Input in = make_input(3, 1.4632, 0);
    in.exactErf = true;
    const macqc::Output ex = macqc::run_calculation(in);
    CHECK(ex.hf && ex.hf->converged);
    if (ex.hf) {
        CHECK(std::fabs(ex.hf->total_energy - out.hf->total_energy) < 1e-4);
        CHECK(ex.hf->total_energy != out.hf->total_energy);
    }
}

static void test_other_basis_sizes_and_curve() {
    for (int n = 1; n <= 2; ++n) {
        const macqc::Output out = macqc::run_calculation(make_input(n, 1.4632, 0));
        CHECK(out.error.empty());
        CHECK(out.hf && out.hf->converged);
        CHECK(out.fci && out.fci->total_energy <= out.hf->total_energy + 1e-9);
    }
    // Potential-curve scan R = 0.8 ... 5.0 bohr.
    int points = 0;
    for (int i = 0; i <= 42; ++i) {
        const double R = 0.8 + 0.1 * i;
        for (int n = 1; n <= 3; ++n) {
            const macqc::Output out = macqc::run_calculation(make_input(n, R, 0));
            ++points;
            const bool ok = out.error.empty() && out.hf && out.hf->converged && out.fci &&
                            std::isfinite(out.fci->total_energy) &&
                            out.fci->correlation_energy <= 1e-8 && out.fci->states[0].spin == "singlet";
            if (!ok) std::printf("curve point failed: n=%d R=%.2f error=%s\n", n, R, out.error.c_str());
            CHECK(ok);
        }
    }
    CHECK(points == 43 * 3);
}

static json run_json(const char* text) {
    char* p = macqc_run(text);
    CHECK(p != nullptr);
    if (!p) return json();
    json j = json::parse(p);
    macqc_free(p);
    return j;
}

static void test_api() {
    const json ok = run_json(
        "{\"R\":1.4632,\"zeta1\":2.0925,\"zeta2\":1.24,\"ZA\":2,\"ZB\":1,\"nGaussians\":3,"
        "\"exactErf\":false,\"maxIter\":25,\"convergence\":1e-4,\"printLevel\":2}");
    CHECK(ok["error"].is_null());
    CHECK(ok["hf"]["converged"].get<bool>());
    CHECK(ok["hf"]["iterations"].get<int>() == 6);
    check_near("json E_tot", ok["hf"]["totalEnergy"].get<double>(), -2.860661992, 1e-8);
    check_near("json S12", ok["hf"]["integrals"]["S12"].get<double>(), 0.450770, 1e-6);
    CHECK(ok["hf"]["scfHistory"].size() == 6);
    CHECK(ok["hf"]["C"].size() == 2 && ok["hf"]["PS"][1].size() == 2);
    check_near("json FCI", ok["fci"]["totalEnergy"].get<double>(), -2.880711848, 1e-8);
    CHECK(ok["fci"]["states"].size() == 4);
    CHECK(ok["fci"]["determinants"].size() == 4);
    CHECK(ok["fci"]["singletWeights"].size() == 3);
    CHECK(ok["fci"]["groundStateVector"].size() == 4);
    CHECK(ok["input"]["nGaussians"].get<int>() == 3);
    CHECK(ok["report"].get<std::string>().find("CALCULATION CONVERGED") != std::string::npos);
    CHECK(ok["report"].get<std::string>().find("FULL CI") != std::string::npos);

    // Defaults when the object is empty.
    CHECK(run_json("{}")["error"].is_null());

    // Invalid input must give an error JSON, never a crash.
    const char* bad[] = {
        "{\"R\":0}", "{\"R\":-1.0}", "{\"zeta1\":0}", "{\"zeta2\":-2}", "{\"nGaussians\":0}",
        "{\"nGaussians\":4}", "{\"nGaussians\":2.5}", "{\"R\":\"abc\"}", "{\"ZA\":-1}",
        "{\"maxIter\":0}", "{\"convergence\":0}", "{\"printLevel\":7}", "not json", "", "[1,2]",
        "{\"exactErf\":3}"};
    for (const char* b : bad) {
        const json j = run_json(b);
        const bool good = j.contains("error") && j["error"].is_string() && j["hf"].is_null() &&
                          j["fci"].is_null();
        if (!good) std::printf("invalid input not rejected: %s\n", b);
        CHECK(good);
    }
    char* p = macqc_run(nullptr);  // NULL must not crash
    CHECK(p != nullptr);
    macqc_free(p);
    macqc_free(nullptr);

    // Non-convergence is reported, not thrown.
    const json nc = run_json("{\"maxIter\":2}");
    CHECK(nc["error"].is_string());
    CHECK(!nc["hf"]["converged"].get<bool>());
    CHECK(nc["fci"].is_null());
}

int main() {
    test_number_formats();
    test_reports();
    test_numbers();
    test_other_basis_sizes_and_curve();
    test_api();
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
