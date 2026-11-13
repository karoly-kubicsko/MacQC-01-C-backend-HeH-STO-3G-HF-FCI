// Energies and nuclear gradients at arbitrary Cartesian geometries for every method of
// ΨLab (5.0): the building block of geometry optimisations and numerical Hessians.
//
// Analytic gradients: RHF (derivative integrals, energy-weighted density, exact term for
// removed near-linear dependencies); CCSD and CCSD(T) (Lambda equations, relaxed densities,
// Z-vector); Full CI (relaxed densities). CCSDT and CCSDTQ, and the correlated methods when
// near-linear dependencies were removed or the analytic gradient does not fit into memory,
// use central differences of energies along (symmetry-adapted) internal coordinates; the
// result says which was used and why.
#pragma once

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/cc.h"
#include "macqc/cc_gradient.h"
#include "macqc/ci.h"
#include "macqc/correlated_gradient.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"

namespace macqc {

struct ElectronicOptions {
    int threads = 0;                 // 0: all cores
    double memory_mb = 1024.0;       // AO two-electron integrals kept in memory up to this size
    double fci_memory_mb = 4096.0;   // memory limit of Full CI, coupled cluster and correlated gradients
    int scf_max_iter = 100;
    double scf_energy_tol = 1.0e-10;
    double scf_gradient_tol = 1.0e-9;
    double lindep_tol = 1.0e-7;
    int cc_max_iter = 150;
    double cc_amplitude_tol = 1.0e-9;
    double cc_energy_tol = 1.0e-10;
    double fci_residual_tol = 1.0e-8;
    double numerical_step = 0.005;   // bohr (stretches) and rad (bends) of numerical gradients
};

// Callbacks for progress and for the text output (each may be empty).
struct ElectronicCallbacks {
    RHFProgress scf;
    CCProgress cc;
    FCIProgress fci;
    LambdaProgress lambda;
    CPHFCallback cphf;                               // Z-vector / CPHF iterations
    std::function<void(const std::string&)> note;    // one line of text (e.g. a numerical step)
};

// Energy of one method at one geometry.
struct SinglePoint {
    Molecule molecule;
    std::string method;           // canonical: RHF, CCSD, CCSD(T), CCSDT, CCSDTQ, FCI
    double energy = 0.0;          // total energy of the method (E_h)
    RHFResult rhf;
    std::optional<CCResult> cc;
    std::optional<FCIResultGeneral> fci;
    std::vector<std::string> warnings;
    bool converged = false;       // SCF and the correlated method
    double seconds = 0.0;
};

struct GradientResult {
    SinglePoint point;            // the calculation at the geometry itself
    std::vector<double> gradient; // 3 n_atoms (E_h/bohr)
    bool analytic = false;
    std::string description;      // how it was computed (for the output and the UI)
    std::string numerical_reason; // why numerical (empty if analytic)
    int energy_evaluations = 0;   // numerical: energies computed (besides the reference)
    std::vector<std::string> internal_coordinates;  // numerical: the coordinates displaced
    std::array<double, 3> dipole{};  // e bohr: SCF dipole (RHF) or relaxed (correlated, analytic)
    bool has_dipole = false;
    std::optional<RHFGradientResult> rhf_gradient;
    std::optional<CorrelatedGradientResult> correlated;
    bool lambda_converged = false;
    int lambda_iterations = 0;
    double lambda_residual = 0.0;
    double seconds = 0.0;
};

// True for the methods with an analytic gradient: RHF, CCSD, CCSD(T), FCI.
bool analytic_gradient_method(const std::string& canonical_method);

// Memory (bytes) of an analytic gradient of `method`, or 0 for RHF.
double analytic_gradient_bytes(const std::string& canonical_method, int n_occupied, int n_mo, int n_basis);

SinglePoint run_single_point(const Molecule& mol, const BasisSet& basis, const std::string& method,
                             const ElectronicOptions& options, const Matrix* guess_orbitals = nullptr,
                             const ElectronicCallbacks& callbacks = ElectronicCallbacks());

// `force_numerical`: central differences of energies even where an analytic gradient exists
// (tests and comparisons).
GradientResult run_gradient(const Molecule& mol, const BasisSet& basis, const std::string& method,
                            const ElectronicOptions& options, const Matrix* guess_orbitals = nullptr,
                            const ElectronicCallbacks& callbacks = ElectronicCallbacks(), bool force_numerical = false);

// Natural internal coordinates of the molecules of ΨLab: a diatomic has R; H2O (atom 0 = O)
// has R(O-H1), R(O-H2) and the angle H1-O-H2. values() in bohr and rad.
struct SimpleInternals {
    std::vector<std::string> names;
    std::vector<double> values;
    std::vector<std::vector<double>> B;  // Wilson B matrix rows (dq/dx), 3 n_atoms columns
};
SimpleInternals simple_internals(const Molecule& mol);
// The geometry with the given internal-coordinate values (same orientation of the first bond).
Molecule geometry_from_internals(const Molecule& reference, const std::vector<double>& values);

}  // namespace macqc
