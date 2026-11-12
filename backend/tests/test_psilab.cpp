// Tests for the ΨLab terminal version (text input -> text output).
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "psilab_text.h"

static int g_checks = 0, g_failures = 0;

#define CHECK(cond)                                                            \
    do {                                                                       \
        ++g_checks;                                                            \
        if (!(cond)) {                                                         \
            ++g_failures;                                                      \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);        \
        }                                                                      \
    } while (0)

static bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

static bool rejects(const std::string& text, const std::string& expected_part) {
    try {
        psilab::parse_input_text(text);
    } catch (const std::exception& e) {
        if (!contains(e.what(), expected_part)) std::printf("  message was: %s\n", e.what());
        return contains(e.what(), expected_part);
    }
    return false;
}

int main() {
    // The example input is the Szabo-Ostlund default case.
    const psilab::Settings ex = psilab::parse_input_text(psilab::example_input());
    CHECK(ex.input.nGaussians == 3);
    CHECK(std::fabs(ex.input.R - 1.4632) < 1e-12);
    CHECK(std::fabs(ex.input.zeta1 - 2.0925) < 1e-12 && std::fabs(ex.input.zeta2 - 1.24) < 1e-12);
    CHECK(!ex.input.exactErf && ex.scan && std::fabs(ex.scanStep - 0.1) < 1e-12);

    // Units, case, aliases, comments.
    const psilab::Settings a = psilab::parse_input_text(
        "  BASIS = sto-2g  # comment\nr = 0.7743 Angstrom\nZeta_He: 2\nzeta-h = 1.2\nExact_Erf = YES\n");
    CHECK(a.input.nGaussians == 2);
    CHECK(std::fabs(a.input.R - 0.7743 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(a.input.zeta1 == 2.0 && a.input.zeta2 == 1.2 && a.input.exactErf && !a.scan);
    CHECK(std::fabs(psilab::parse_input_text("R = 1.0 \xC3\x85\n").input.R - 1.0 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(psilab::parse_input_text("R = 2\n").input.R == 2.0);  // default unit bohr

    // Errors carry the line number and say what is wrong.
    CHECK(rejects("basis = STO-3G\nbasis = STO-2G\n", "line 2"));
    CHECK(rejects("basis = 6-31G\n", "STO-1G, STO-2G or STO-3G"));
    CHECK(rejects("R = abc\n", "needs a number"));
    CHECK(rejects("R = 1.4 parsec\n", "unknown length unit"));
    CHECK(rejects("colour = blue\n", "unknown key"));
    CHECK(rejects("just some text\n", "expected 'key = value'"));
    CHECK(rejects("exact_erf = maybe\n", "yes or no"));
    CHECK(rejects("scan_step = 0\n", "scan_step"));
    CHECK(rejects("R = -1\n", "R must be"));
    CHECK(rejects("zeta_H = 0\n", "zeta2"));

    // Scan grid: 0.8 ... 5.0 inclusive.
    const auto grid = psilab::scan_grid(0.1);
    CHECK(grid.size() == 43 && std::fabs(grid.front() - 0.8) < 1e-12 && std::fabs(grid.back() - 5.0) < 1e-9);
    CHECK(std::fabs(psilab::scan_grid(0.25).back() - 5.0) < 1e-12);  // end point appended

    // Full run of the example: energies, report, curve and minima.
    const psilab::RunResult r = psilab::run_to_text(ex);
    CHECK(r.error.empty());
    CHECK(contains(r.text, "\xCE\xA8Lab \xE2\x80\x93 HeH+ (STO-3G)"));
    CHECK(contains(r.text, "R = 0.7743 \xC3\x85 (1.4632 bohr)"));
    CHECK(contains(r.text, "total energy:         -2.86066199 E_h"));
    CHECK(contains(r.text, "total energy:         -2.88071185 E_h"));
    CHECK(contains(r.text, "correlation energy:   -0.02004986 E_h"));
    CHECK(contains(r.text, "state 2: -2.07384042 E_h  <S^2> = 2.0000  triplet"));
    CHECK(contains(r.text, "STO-3G FOR ATOMIC NUMBERS  2.00 AND  1.00"));  // FORTRAN-style report
    CHECK(contains(r.text, "CORRELATION ENERGY ="));
    CHECK(contains(r.text, "POTENTIAL CURVE"));
    CHECK(contains(r.text, "HF minimum (parabolic fit): R = 0.7323 \xC3\x85 (1.3838 bohr), E = -2.86277416 E_h"));
    CHECK(contains(r.text, "FCI minimum (parabolic fit): R = 0.7397 \xC3\x85 (1.3978 bohr), E = -2.88223178 E_h"));
    CHECK(!contains(r.text, "failed point"));

    // Without scan: no curve section.
    psilab::Settings single = ex;
    single.scan = false;
    CHECK(!contains(psilab::run_to_text(single).text, "POTENTIAL CURVE"));

    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
