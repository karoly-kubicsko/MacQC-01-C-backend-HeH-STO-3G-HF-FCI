// Molecular Hessians for harmonic frequencies (ΨLab 5.0).
//
//  - RHF: analytic (second-derivative integrals and CPHF, rhf_derivatives.h), with the analytic
//    dipole derivatives; when near-linear dependencies were removed, numerical (below).
//  - CCSD, CCSD(T), Full CI (and RHF with removed functions): central differences of analytic
//    gradients along internal coordinates; for H2O with C2v symmetry the symmetry-adapted
//    coordinates S1 (symmetric stretch), S2 (bend) and S3 (antisymmetric stretch), otherwise
//    R(O-H1), R(O-H2) and the angle; for diatomics R. Two gradients per coordinate. The
//    Cartesian Hessian is B' H_q B + sum_i g_q,i d2q_i/dx2 (exact also away from a stationary
//    point). The displaced geometries are aligned with the reference (Eckart frame, mass-weighted
//    superposition), so the relaxed dipole moments of the gradients give the dipole derivatives.
//  - CCSDT, CCSDTQ: central differences of energies along the same coordinates (second
//    derivatives, mixed terms only where symmetry allows them); dipole moments by finite fields
//    (mu = -dE/dF) at the displaced geometries.
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "macqc/gradient.h"
#include "macqc/linalg.h"

namespace macqc {

struct HessianOptions {
    double step = 0.005;        // bohr / rad
    double field = 2.0e-4;      // a.u., finite-field dipoles (CCSDT, CCSDTQ)
    bool dipoles = true;        // IR intensities
    bool force_numerical = false;
    bool from_energies = false;  // numerical from energies even where analytic gradients exist (tests)
    // Energy and gradient already computed at this geometry (e.g. the last optimisation step),
    // used instead of computing them again; nullptr: computed here.
    const GradientResult* reference = nullptr;
};

struct HessianResult {
    Matrix hessian;             // Cartesian, E_h/bohr^2
    bool analytic = false;
    std::string description;
    std::string numerical_reason;
    int gradient_evaluations = 0, energy_evaluations = 0;
    std::vector<std::string> coordinates;  // displaced (symmetry-adapted) coordinates
    bool has_dipole_derivatives = false;
    std::vector<std::array<double, 3>> dipole_derivatives;  // [3 n_atoms][xyz], e
    std::string dipole_description;
    std::array<double, 3> dipole{};
    bool has_dipole = false;
    GradientResult reference;   // energy and gradient at the geometry
    int cphf_iterations = 0;
    bool cphf_converged = true;
    double seconds = 0.0;
};

// Called before each displaced calculation: (index, total, text).
using DisplacementProgress = std::function<void(int, int, const std::string&)>;

HessianResult run_hessian(const Molecule& mol, const BasisSet& basis, const std::string& method,
                          const ElectronicOptions& electronic, const HessianOptions& options = HessianOptions(),
                          const ElectronicCallbacks& callbacks = ElectronicCallbacks(),
                          const Matrix* guess_orbitals = nullptr,
                          const DisplacementProgress& progress = DisplacementProgress());

// Number of gradient (analytic-gradient methods) or energy (CCSDT, CCSDTQ) evaluations of a
// numerical Hessian of this molecule, without the reference point (for the estimates).
struct NumericalHessianCost {
    int gradients = 0, energies = 0, field_energies = 0;
    std::vector<std::string> coordinates;
};
NumericalHessianCost numerical_hessian_cost(const Molecule& mol, bool analytic_gradients, bool dipoles);

// `moved` translated and rotated onto `reference` (mass-weighted least squares, Horn's
// quaternion method): the Eckart frame to first order.
Molecule align_to(const Molecule& reference, const Molecule& moved);

}  // namespace macqc
