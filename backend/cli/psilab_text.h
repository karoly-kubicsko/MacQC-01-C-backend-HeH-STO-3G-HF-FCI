// ΨLab terminal version: plain-text input file -> plain-text output file.
//
// Input: one "key = value" per line, '#' starts a comment. Keys (case-insensitive):
//   molecule   HeH+ (default) | H2 | LiH | H2O
//   basis      HeH+: STO-1G | STO-2G | STO-3G; others: aug-cc-pVXZ, aug-cc-pCVXZ (X = D, T, Q, 5)
//   R          bond length (H2O: O-H), a number optionally followed by a unit: bohr (default)
//              or angstrom; default for H2, LiH, H2O: the experimental equilibrium value
//   angle      H-O-H angle in degrees (H2O; default 104.52)
//   method     RHF (H2, LiH, H2O)
//   threads, memory (MB or GB), basis_dir      (H2, LiH, H2O)
//   zeta_He, zeta_H, exact_erf, scan, scan_step (HeH+ only):
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

// ΨLab version (same as the macOS app).
constexpr const char* kVersion = "2.0";

// 1 bohr in angstrom (CODATA 2022).
constexpr double kBohrToAngstrom = 0.529177210544;

constexpr double kScanStart = 0.8;  // bohr
constexpr double kScanEnd = 5.0;    // bohr

struct Settings {
    macqc::Input input;     // R always in bohr
    bool scan = false;
    double scanStep = 0.1;  // bohr
    bool rGiven = false;    // R was in the input (else the molecule's default)
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

// H2, LiH, H2O: size of the calculation, time estimate and warnings, without computing
// (empty for HeH+).
std::string estimate_text(const Settings& settings);
// The warnings of the estimate (e.g. a long calculation), for the terminal.
std::vector<std::string> estimate_warnings(const Settings& settings);

}  // namespace psilab
