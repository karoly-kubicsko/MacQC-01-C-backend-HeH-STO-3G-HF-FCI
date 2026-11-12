// ΨLab terminal version: plain-text input file -> plain-text output file.
//
// Input: one "key = value" per line, '#' starts a comment. Keys (case-insensitive):
//   basis      STO-1G | STO-2G | STO-3G
//   R          bond length, a number optionally followed by a unit: bohr (default) or angstrom
//   zeta_He    Slater exponent of He          zeta_H     Slater exponent of H
//   exact_erf  yes | no                       scan       yes | no
//   scan_step  step of the potential-curve scan (bohr), R = 0.8 ... 5.0 bohr
// The output text matches the macOS app's File > Export Results as Text.
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "macqc/calculation.h"

namespace psilab {

// 1 bohr in angstrom (CODATA 2022).
constexpr double kBohrToAngstrom = 0.529177210544;

constexpr double kScanStart = 0.8;  // bohr
constexpr double kScanEnd = 5.0;    // bohr

struct Settings {
    macqc::Input input;     // R always in bohr
    bool scan = false;
    double scanStep = 0.1;  // bohr
};

// Parses the text input. Throws std::runtime_error("line N: ...") on any problem,
// including values the calculation would reject.
Settings parse_input_text(const std::string& text);

// A commented example input file (the Szabo-Ostlund default case).
std::string example_input();

struct CurvePoint {
    double R, hf, fci;  // bohr, E_h, E_h
};

struct CurveMinimum {
    std::string method;
    double R = 0.0, energy = 0.0;
    bool refined = false;  // parabola through the lowest point and its neighbours
};

std::vector<double> scan_grid(double step);
std::optional<CurveMinimum> curve_minimum(const std::vector<std::pair<double, double>>& samples,
                                          const std::string& method);

struct RunResult {
    std::string text;   // the complete output file (empty only if nothing could be computed)
    std::string error;  // empty on success
    std::optional<macqc::Output> output;  // the single-point calculation
};

// Runs the single-point calculation (and the scan if requested) and formats the output.
RunResult run_to_text(const Settings& settings);

}  // namespace psilab
