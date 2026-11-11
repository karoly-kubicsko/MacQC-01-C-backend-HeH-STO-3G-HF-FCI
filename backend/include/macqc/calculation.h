// High-level driver: validated input -> HF -> FCI -> text report / JSON.
#pragma once

#include <optional>
#include <string>

#include "macqc/fci.h"
#include "macqc/scf.h"

namespace macqc {

struct Input {
    double R = 1.4632;
    double zeta1 = 2.0925;
    double zeta2 = 1.24;
    double ZA = 2.0;
    double ZB = 1.0;
    int nGaussians = 3;
    bool exactErf = false;
    int maxIter = 25;
    double convergence = 1.0e-4;
    int printLevel = 2;
};

struct Output {
    Input input;
    std::string error;  // empty: no error
    std::optional<HFResult> hf;
    std::optional<FCIResult> fci;
    std::string hf_report;   // HF part only (as printed by the reference program)
    std::string report;      // hf_report + FCI section
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
