// High-level driver: validated input -> HF -> FCI -> text report / JSON.
#pragma once

#include <optional>
#include <string>

#include "macqc/fci.h"
#include "macqc/molecular.h"
#include "macqc/scf.h"

namespace macqc {

struct Input {
    // "HeH+" (default): the Szabo-Ostlund STO-nG calculation (HF + FCI), unchanged
    // since version 1.0, which uses the keys R .. printLevel below.
    // "H2", "LiH", "H2O": general molecules (version 3.0), see the keys after printLevel.
    std::string molecule = "HeH+";
    double R = 1.4632;  // bohr; for H2O the O-H distance
    double zeta1 = 2.0925;
    double zeta2 = 1.24;
    double ZA = 2.0;
    double ZB = 1.0;
    int nGaussians = 3;
    bool exactErf = false;
    int maxIter = 25;
    double convergence = 1.0e-4;
    int printLevel = 2;

    // General molecules (molecule != "HeH+").
    std::string basis = "aug-cc-pVDZ";
    double angle = kDefaultWaterAngle;  // H-O-H angle (degrees), H2O only
    std::string method = "RHF";
    std::string basisDirectory;         // "" -> default_basis_directory()
    int threads = 0;                    // 0 -> all hardware threads
    double memoryMB = 1024.0;           // limit for keeping the two-electron integrals in memory
    int scfMaxIter = 100;
    double energyConvergence = 1.0e-10;
    double gradientConvergence = 1.0e-7;
    bool estimateOnly = false;          // only molecule, basis, size and time estimate, warnings
};

// True if the input selects the HeH+ STO-nG calculation of versions 1 and 2.
bool is_heh_input(const Input& in);

struct Output {
    Input input;
    std::string error;  // empty: no error
    std::optional<HFResult> hf;
    std::optional<FCIResult> fci;
    std::string hf_report;   // HF part only (as printed by the reference program)
    std::string report;      // hf_report + FCI section (general molecules: their report)
    std::optional<MolecularResult> molecular;  // general molecules only
};

// Returns an empty string if the input is valid, else a message.
std::string validate(const Input& in);

// Never throws; problems are reported through Output::error.
Output run_calculation(const Input& in);

// Parse input JSON text (missing keys take defaults). Throws std::runtime_error
// with a descriptive message on malformed input.
Input parse_input_json(const std::string& text);

// Serialise to the documented output JSON (never throws).
std::string output_to_json(const Output& out, int indent = -1);

}  // namespace macqc
