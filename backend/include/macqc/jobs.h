// Job types of ΨLab 5.0 for H2, LiH and H2O, with the complete output file:
//   energy        single-point energy (as in 4.0)
//   opt           geometry optimisation
//   freq          harmonic frequencies and thermochemistry at the input geometry
//   opt+freq      optimisation, then frequencies at the optimised geometry
//   opt+freq+sp   optimisation, frequencies, then a single-point energy with its own method
//                 and basis set at the optimised geometry; the final energy is reported with the
//                 thermal corrections of the frequency step
// If the optimisation does not converge, jobs 4 and 5 stop after it with a message.
#pragma once

#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "macqc/frequencies.h"
#include "macqc/gradient.h"
#include "macqc/hessian.h"
#include "macqc/molecular.h"
#include "macqc/optimizer.h"

namespace macqc {

struct JobSettings {
    std::string type = "energy";      // energy | opt | freq | opt+freq | opt+freq+sp
    std::string sp_method, sp_basis;  // job 5 (default: the method and basis set of the job)
    double temperature = 298.15;      // K
    double pressure = 101325.0;       // Pa
    bool tight = false;               // optimisation convergence
    int max_steps = 50;
};

// Canonical job type for any spelling ("opt+freq", "OPT FREQ", "optfreq", "4"), or "".
std::string canonical_job_type(const std::string& type);
// "Energy", "Geometry optimisation", "Frequencies", "Optimisation + frequencies",
// "Optimisation + frequencies + single-point energy".
std::string job_title(const std::string& type);
bool job_has_optimization(const std::string& type);
bool job_has_frequencies(const std::string& type);
bool job_has_single_point(const std::string& type);

struct JobStepEstimate {
    std::string name;          // "Geometry optimisation", "Hessian", "Single-point energy", "Energy"
    std::string method, basis;
    std::string how;           // e.g. "analytic gradients (Lambda equations, relaxed densities, Z-vector)"
    bool analytic = true;
    int evaluations = 0;       // energies or gradients (per optimisation step for "opt")
    std::string evaluation_kind;  // "gradients", "energies", "steps"
    double seconds = 0.0;
    double memory_mb = 0.0;
};

struct JobEstimate {
    std::vector<JobStepEstimate> steps;
    double total_seconds = 0.0;
    double memory_mb = 0.0;
    std::string refusal;                 // non-empty: the job is not possible
    std::vector<std::string> warnings;   // long calculation etc.
    std::vector<std::string> notes;      // what is analytic / numerical, ...
};

// Uses prepare_molecular() for the cost of single energies (sp_prepared: the method and basis
// of the single-point step of job 5, or nullptr).
JobEstimate estimate_job(const MolecularResult& prepared, const JobSettings& settings,
                         const MolecularResult* sp_prepared);

struct JobResult {
    JobSettings settings;
    JobEstimate estimate;
    bool completed = false;
    std::string termination;           // "normal", "error", "refused", "cancelled"
    std::string message;               // reason for an error or a stop
    Molecule initial_geometry, final_geometry;
    std::optional<SinglePoint> point;           // energy job / single point at the final geometry
    std::optional<OptimizationResult> optimization;
    std::optional<HessianResult> hessian;
    std::optional<HarmonicAnalysis> harmonic;
    std::optional<ThermoResult> thermo;
    std::optional<SinglePoint> single_point;    // job 5
    std::string single_point_method, single_point_basis;
    double final_energy = 0.0;
    std::string final_energy_label;             // e.g. "CCSD(T)/aug-cc-pVTZ"
    std::vector<std::string> warnings;
    double seconds = 0.0;
    std::string output;                         // the complete output file
};

using OutputSink = std::function<void(const std::string&)>;

// Runs a job (not "energy": that path is run_energy_job). `input` gives molecule, starting
// geometry, method, basis and limits. Never throws for calculation problems: they end the
// job with termination "error" and the reason in the output; Cancelled is rethrown after the
// output says so.
JobResult run_job(const MolecularInput& input, const JobSettings& settings, const OutputSink& sink,
                  const std::string& input_echo);

// The energy job with the 5.0 output file: runs run_molecular (4.0 results unchanged) and
// writes the output around it.
JobResult run_energy_job(MolecularResult& prepared, const JobSettings& settings, const OutputSink& sink,
                         const std::string& input_echo);

// The "Size, memory and time (estimates)" section of the output file.
std::string job_estimate_text(const MolecularResult& prepared, const JobEstimate& estimate);

// A MolecularResult (4.0 structure) for the final energy calculation of a job.
MolecularResult final_point_result(const JobResult& job, const MolecularInput& input);

// The program version.
constexpr const char* kProgramVersion = "5.1";

}  // namespace macqc
