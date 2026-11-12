// Closed-shell restricted Hartree-Fock for general molecules in contracted
// spherical Gaussian basis sets (H2, LiH, H2O in ΨLab 3.0).
#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gto_integrals.h"
#include "macqc/linalg.h"
#include "macqc/molecule.h"

namespace macqc {

// Two-electron part of the Fock matrix, G(P) = J(P) - K(P)/2 for a total density P.
// The unique integrals (ij|kl) are kept in memory when they fit into the memory
// limit ("stored"); otherwise they are recomputed for every Fock build ("direct",
// with Schwarz and density screening).
class FockBuilder {
public:
    FockBuilder(const BasisSet& basis, int threads, double memory_limit_mb);

    Matrix build(const Matrix& P) const;

    bool stored() const { return stored_; }
    std::size_t unique_integrals() const { return unique_; }
    std::size_t significant_shell_quartets() const { return quartets_; }
    double setup_seconds() const { return setup_seconds_; }
    int threads() const { return threads_; }

    // Packed unique (ij|kl), i >= j, k >= l, ij >= kl (stored mode only).
    const std::vector<double>& packed() const { return eri_; }
    static std::size_t pair_index(std::size_t i, std::size_t j) {
        return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
    }

private:
    Matrix build_stored(const Matrix& P) const;
    Matrix build_direct(const Matrix& P) const;

    const BasisSet& basis_;
    ShellPairList pairs_;
    int threads_ = 1;
    bool stored_ = false;
    std::size_t unique_ = 0, quartets_ = 0;
    double setup_seconds_ = 0.0;
    std::vector<double> eri_;
};

// Number of unique two-electron integrals for n basis functions, n^4/8 asymptotically.
std::size_t unique_eri_count(std::size_t n);

struct RHFOptions {
    int max_iter = 100;
    double energy_tol = 1.0e-10;    // |E_n - E_(n-1)| (E_h)
    double gradient_tol = 1.0e-7;   // max |X^T (FPS - SPF) X|
    double lindep_tol = 1.0e-7;     // overlap eigenvalues below this are removed
    int diis_size = 8;
    int threads = 0;                // 0: all hardware threads
    double memory_limit_mb = 1024;  // for storing the two-electron integrals
};

struct RHFIteration {
    int iteration = 0;
    double energy = 0.0;      // total energy (E_h)
    double delta_e = 0.0;
    double gradient = 0.0;    // max |orbital gradient|
    bool diis = false;
};

struct RHFResult {
    bool converged = false;
    int iterations = 0;
    double electronic_energy = 0.0;
    double nuclear_repulsion = 0.0;
    double total_energy = 0.0;
    double kinetic_energy = 0.0;           // Tr(P T)
    double nuclear_attraction_energy = 0.0;  // Tr(P V)
    double two_electron_energy = 0.0;      // Tr(P G)/2
    double virial_ratio = 0.0;             // -V/T
    int n_basis = 0;
    int n_mo = 0;          // after removal of linear dependencies
    int n_occupied = 0;
    double smallest_overlap_eigenvalue = 0.0;
    int removed_functions = 0;
    std::vector<double> orbital_energies;  // ascending, n_mo values
    Matrix C;              // AO x MO coefficients
    Matrix P;              // total density
    Matrix H;              // core Hamiltonian
    Matrix S;              // overlap
    std::vector<RHFIteration> history;
    bool eri_stored = false;
    std::size_t unique_integrals = 0;
    int threads = 1;
    double integral_seconds = 0.0;  // one-electron + stored ERIs (or Schwarz setup)
    double scf_seconds = 0.0;
};

// Progress callback: called after each SCF iteration (may be empty).
using RHFProgress = std::function<void(const RHFIteration&)>;

// Throws std::runtime_error for impossible input (e.g. odd electron count).
// Non-convergence is reported through RHFResult::converged.
RHFResult run_rhf(const Molecule& mol, const BasisSet& basis, const RHFOptions& options,
                  const RHFProgress& progress = RHFProgress());

}  // namespace macqc
