// ΨLab 3.0 driver for general molecules (H2, LiH, H2O): basis set, size and time
// estimate, warnings and refusals, restricted Hartree-Fock, text report.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/ci.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"

namespace macqc {

struct MolecularInput {
    std::string molecule = "H2O";
    double R = 0.0;              // bond length in bohr (H2O: O-H); <= 0 means the default
    double angle = kDefaultWaterAngle;  // H-O-H angle in degrees (H2O only)
    std::string basis = "aug-cc-pVDZ";
    std::string method = "RHF";  // "RHF" or "FCI" (RHF followed by Full CI)
    std::string basis_directory;  // "" -> default_basis_directory()
    int threads = 0;              // 0 -> all hardware threads
    double memory_mb = 1024.0;    // limit for keeping the two-electron integrals in memory
    double fci_memory_mb = 4096.0;  // hard memory limit of a Full CI calculation
    int max_iter = 100;
    double energy_tol = 1.0e-10;
    double gradient_tol = 1.0e-7;
    int print_level = 2;          // 0 silent, 1 summary, 2 full
};

// Size and cost of a calculation, known before it runs.
struct CostEstimate {
    int basis_functions = 0;
    int shells = 0;
    int max_l = 0;
    std::size_t unique_eris = 0;
    double eri_memory_mb = 0.0;
    bool eri_stored = false;
    int threads = 1;
    double rhf_seconds = 0.0;    // rough estimate of the wall time
    std::optional<FCISize> fci;  // method FCI only
    double total_seconds = 0.0;  // RHF (+ Full CI)
};

struct MolecularResult {
    MolecularInput input;
    Molecule molecule;
    BasisSet basis;
    CostEstimate estimate;
    std::vector<std::string> notes;     // e.g. which basis set is used for H
    std::vector<std::string> warnings;  // e.g. long calculation, removed functions
    std::string refusal;                // non-empty: the calculation is not possible (reason)
    std::optional<RHFResult> rhf;
    std::optional<FCIResultGeneral> fci;
    std::string report;
};

// Validates the input; returns "" or a message (unknown molecule/basis/method, bad geometry).
std::string validate_molecular(const MolecularInput& in);

// Builds molecule and basis and estimates the cost without computing anything.
// Throws std::runtime_error with an explanation if the input cannot be set up (e.g. basis
// set not defined for an element). A Full CI that is too large is not an exception here:
// result.refusal explains why it is not possible (run_molecular then refuses to run).
MolecularResult prepare_molecular(const MolecularInput& in);

// Wall-time estimate (seconds) for an RHF calculation in this basis: a cost model of the
// integral work per shell quartet, times the expected number of SCF iterations when the
// integrals are recomputed (integral-direct). Typically within a factor of 2.
double estimate_rhf_seconds(const BasisSet& basis, bool eri_stored, int threads, int iterations);

// Typical number of SCF iterations for the molecule near its equilibrium geometry.
int expected_scf_iterations(const std::string& molecule);

// Runs the calculation (after prepare_molecular). Throws std::runtime_error on errors,
// including a refused calculation (result.refusal).
void run_molecular(MolecularResult& result, const RHFProgress& progress = RHFProgress(),
                   const FCIProgress& fci_progress = FCIProgress());

// Report text (also stored in result.report by run_molecular).
std::string molecular_report(const MolecularResult& result);

// Estimated-time threshold above which a warning is issued (seconds).
constexpr double kLongCalculationSeconds = 30.0;

}  // namespace macqc
