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

    // General molecules: molecule, basis set, geometry defaults, keys per molecule.
    const psilab::Settings w = psilab::parse_input_text("molecule = water\nbasis = AUG-CC-PVTZ\n");
    CHECK(w.input.molecule == "H2O" && w.input.basis == "aug-cc-pVTZ");
    CHECK(std::fabs(w.input.R - 0.9572 / psilab::kBohrToAngstrom) < 1e-12 && w.input.angle == 104.52);
    const psilab::Settings w2 = psilab::parse_input_text(
        "molecule = H2O\nR = 1.0 angstrom\nangle = 100\nmethod = rhf\nthreads = 2\nmemory = 2 GB\n");
    CHECK(std::fabs(w2.input.R - 1.0 / psilab::kBohrToAngstrom) < 1e-12 && w2.input.angle == 100.0);
    CHECK(w2.input.threads == 2 && w2.input.memoryMB == 2048.0 && w2.input.basis == "aug-cc-pVDZ");
    CHECK(std::fabs(psilab::parse_input_text("molecule = LiH\n").input.R - 1.5949 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(rejects("molecule = CH4\n", "molecule must be HeH+, H2, LiH or H2O"));
    CHECK(rejects("molecule = H2\nzeta_H = 1.2\n", "line 2: 'zeta_H' applies to HeH+ only"));
    CHECK(rejects("molecule = H2\nangle = 100\n", "applies to H2O only"));
    CHECK(rejects("angle = 100\n", "applies to H2O only"));
    CHECK(rejects("method = RHF\n", "HeH+ always runs Hartree-Fock and Full CI"));
    CHECK(rejects("molecule = H2O\nbasis = STO-3G\n", "basis must be one of aug-cc-pVDZ"));
    CHECK(rejects("molecule = H2\nmethod = CCSD\n", "method must be RHF or FCI"));
    CHECK(rejects("molecule = H2\nfci_memory = 1 PB\n", "fci_memory unit"));
    const psilab::Settings fs = psilab::parse_input_text("molecule = H2\nmethod = Full-CI\nfci_memory = 2 GB\n");
    CHECK(fs.input.method == "FCI" && fs.input.fciMemoryMB == 2048.0);
    CHECK(rejects("molecule = H2O\nangle = 10\n", "angle must be between 30 and 180"));
    CHECK(rejects("molecule = H2\nthreads = 1.5\n", "threads must be an integer"));
    CHECK(rejects("molecule = H2\nmemory = 3 TB\n", "memory unit"));

    const psilab::RunResult h2 = psilab::run_to_text(psilab::parse_input_text("molecule = H2\nbasis = aug-cc-pVDZ\n"));
    CHECK(h2.error.empty());
    CHECK(contains(h2.text, "\xCE\xA8Lab \xE2\x80\x93 H2, RHF/aug-cc-pVDZ"));
    CHECK(contains(h2.text, "total energy:         -1.1287937367 E_h"));
    CHECK(contains(h2.text, "18 basis functions (H [3s2p])"));
    CHECK(contains(h2.text, "DETAILED REPORT"));
    CHECK(contains(psilab::estimate_text(psilab::parse_input_text("molecule = H2O\nbasis = aug-cc-pCVDZ\n")),
                   "aug-cc-pCVDZ is not defined for H"));
    // Full CI: H2 runs, H2O is refused with an explanation.
    const psilab::RunResult h2fci =
        psilab::run_to_text(psilab::parse_input_text("molecule = H2\nbasis = aug-cc-pVDZ\nmethod = FCI\n"));
    CHECK(h2fci.error.empty());
    CHECK(contains(h2fci.text, "H2, FCI/aug-cc-pVDZ"));
    CHECK(contains(h2fci.text, "total energy:         -1.1646244795 E_h"));
    CHECK(contains(h2fci.text, "determinants:         324 (18 orbitals, 2 electrons)"));
    CHECK(contains(h2fci.text, "Natural orbital occupation numbers"));
    const psilab::RunResult h2ofci = psilab::run_to_text(psilab::parse_input_text("molecule = H2O\nmethod = FCI\n"));
    CHECK(h2ofci.text.empty() && contains(h2ofci.error, "impossible"));
    CHECK(contains(psilab::estimate_text(psilab::parse_input_text("molecule = LiH\nbasis = aug-cc-pVQZ\nmethod = FCI\n")),
                   "NOT POSSIBLE"));
    const psilab::RunResult lih =
        psilab::run_to_text(psilab::parse_input_text("molecule = LiH\nbasis = aug-cc-pCV5Z\n"));
    CHECK(lih.text.empty() && contains(lih.error, "not available for Li"));

    // Without scan: no curve section.
    psilab::Settings single = ex;
    single.scan = false;
    CHECK(!contains(psilab::run_to_text(single).text, "POTENTIAL CURVE"));

    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
