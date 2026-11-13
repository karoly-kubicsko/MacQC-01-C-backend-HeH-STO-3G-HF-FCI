// Geometry optimisation in redundant internal coordinates (ΨLab 5.0).
//
// Each step: energy and Cartesian gradient (analytic where available, see gradient.h);
// internal gradient g_q = G^- B g_x; projected gradient and Hessian (P = G G^-,
// H~ = P H P + 1000 (1 - P)); rational-function (RFO) step from the lowest eigenvector of the
// augmented Hessian [[H~, g~], [g~', 0]], scaled down to the trust radius; iterative
// back-transformation to Cartesian coordinates; BFGS update of the internal Hessian with the
// actual change of the internal coordinates and gradients (skipped when the curvature
// condition dq.dg > 0 fails). The trust radius adapts to the ratio of the actual and the
// predicted energy change (Fletcher). Initial Hessian: the analytic RHF Hessian transformed
// to internal coordinates (including the gradient term) when the method is RHF and the basis
// is small, otherwise Lindh's model Hessian (Chem. Phys. Lett. 241, 423 (1995)) or
// Schlegel's (Theor. Chim. Acta 66, 333 (1984)).
//
// Convergence, as in Gaussian: maximum and RMS internal force (E_h/bohr, E_h/rad) and maximum
// and RMS internal displacement of the next step (bohr, rad), plus the energy change since
// the previous step; converged when all are met, or when the forces are 100 times below
// their thresholds.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gradient.h"
#include "macqc/internal_coordinates.h"
#include "macqc/molecule.h"

namespace macqc {

struct ConvergenceCriteria {
    double max_force = 4.5e-4, rms_force = 3.0e-4, max_displacement = 1.8e-3, rms_displacement = 1.2e-3;
    double energy = 1.0e-6;
};
// "normal" (Gaussian's default) or "tight" (Gaussian's tight: 1.5e-5, 1.0e-5, 6.0e-5, 4.0e-5; 1e-8 E_h).
ConvergenceCriteria convergence_criteria(bool tight);

struct OptimizationOptions {
    bool tight = false;
    int max_steps = 50;
    double trust_radius = 0.3;   // initial (bohr, rad)
    double max_trust = 0.5, min_trust = 1.0e-3;
    std::string initial_hessian = "auto";  // auto | analytic | lindh | schlegel
    int analytic_hessian_max_functions = 100;  // "auto": analytic RHF Hessian up to this size
    // Called before the energy and gradient of each step (for the output file).
    std::function<void(int, const Molecule&)> on_step_start;
};

struct ConvergenceItem {
    std::string name;
    double value = 0.0, threshold = 0.0;
    bool met = false;
    bool applicable = true;  // the energy change has no value at the first step
};

struct OptimizationStep {
    int step = 0;
    Molecule geometry;                    // where energy and gradient were computed
    double energy = 0.0, delta_e = 0.0;   // delta_e: change since the previous step
    std::string gradient_description;     // analytic / numerical, how
    std::vector<std::string> names, definitions;  // of the internal coordinates
    std::vector<double> q, gq, dq;        // values, (projected) gradient, step taken from here
    std::vector<double> cartesian_gradient;
    std::string hessian_update;           // how the Hessian of this step was obtained
    double trust_radius = 0.0, step_norm = 0.0, rfo_lambda = 0.0, predicted_change = 0.0;
    bool step_scaled = false;
    double ratio = 0.0;                   // actual / predicted change of the previous step
    bool has_ratio = false;
    int backtransform_iterations = 0;
    double backtransform_residual = 0.0;
    std::vector<ConvergenceItem> convergence;
    bool converged = false;
    std::vector<std::string> warnings;
    double seconds = 0.0;
};

struct GeometricParameter {
    std::string name;       // e.g. "R(O1-H2)"
    double value = 0.0;     // bohr or degrees
    bool angle = false;
};

struct OptimizationResult {
    bool converged = false;
    int steps = 0;
    std::string message;
    std::string initial_hessian;
    ConvergenceCriteria criteria;
    Molecule initial, final_geometry;
    double initial_energy = 0.0, final_energy = 0.0;
    std::vector<OptimizationStep> history;
    GradientResult final_point;           // energy, gradient, orbitals at the final geometry
    std::vector<GeometricParameter> initial_parameters, final_parameters;
    bool symmetric_start = false;         // H2O: started C2v
    bool symmetric_end = false;           // H2O: C2v kept
    double symmetry_deviation = 0.0;      // |R1 - R2| at the end (bohr)
    bool removed_functions_changed = false;
    double seconds = 0.0;
};

using OptimizationCallback = std::function<void(const OptimizationStep&)>;

// `basis` may be built at any geometry of the molecule (shells are moved to each geometry).
OptimizationResult optimize_geometry(const Molecule& start, const BasisSet& basis, const std::string& method,
                                     const ElectronicOptions& electronic, const OptimizationOptions& options,
                                     const ElectronicCallbacks& callbacks = ElectronicCallbacks(),
                                     const OptimizationCallback& on_step = OptimizationCallback());

// Bond lengths (bohr) and angles (degrees) of the internal coordinates of a geometry.
std::vector<GeometricParameter> geometric_parameters(const Molecule& mol);

// Model Hessians in internal coordinates (diagonal).
std::vector<double> lindh_model_hessian(const Molecule& mol, const InternalCoordinates& ic);
std::vector<double> schlegel_model_hessian(const Molecule& mol, const InternalCoordinates& ic);

}  // namespace macqc
