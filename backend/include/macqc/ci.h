// Determinant-based Full CI for closed-shell molecules (ΨLab 3.0: H2, LiH).
//
// The CI vector is a matrix C(I_alpha, I_beta) over alpha and beta occupation strings
// (M_S = 0, all electrons correlated). sigma = H C is built from
//   - the same-spin parts: C F + F C, with F the one-spin string Hamiltonian;
//   - the opposite-spin part through the (N-1)-electron resolution of the identity:
//       D_jl(K_a,K_b) = <K_a|a_j|J_a><K_b|a_l|J_b> C(J_a,J_b),
//       E_ik = sum_jl (ij|kl) D_jl   (a matrix product, cost n^6 for two electrons per spin),
//       sigma(I_a,I_b) += <I_a|a+_i|K_a><I_b|a+_k|K_b> E_ik(K_a,K_b).
// The lowest root is found with the Davidson method (diagonal preconditioner); trial
// vectors are kept symmetric, C = C^T, which removes the odd-spin (triplet) components.
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "macqc/mo_integrals.h"

namespace macqc {

// Occupation strings of N electrons in n orbitals (one spin), addressed by the
// combinatorial number system: index = sum_k binom(o_k, k + 1) for o_0 < o_1 < ...
class StringSpace {
public:
    StringSpace(int n_orbitals, int n_electrons);
    std::size_t size() const { return count_; }
    int orbitals() const { return n_; }
    int electrons() const { return N_; }
    const std::uint16_t* occupation(std::size_t I) const { return occ_.data() + I * static_cast<std::size_t>(N_); }
    std::size_t index(const std::uint16_t* sorted) const;

private:
    int n_ = 0, N_ = 0;
    std::size_t count_ = 0;
    std::vector<std::uint16_t> occ_;
    std::vector<std::vector<std::size_t>> binom_;  // binom_[i][k] = C(i, k), k <= N + 1
};

// Number of strings C(n, N) as a floating-point number (no overflow).
double binomial_count(int n, int N);

struct FCISize {
    int orbitals = 0, alpha = 0, beta = 0;
    double alpha_strings = 0.0, beta_strings = 0.0;
    double determinants = 0.0;
    double memory_bytes = 0.0;   // peak, estimated
    double seconds = 0.0;        // integral transformation + Davidson, estimated
    bool feasible = true;
    std::string reason;          // why not feasible
};

// Hard limit on the determinant space (about 160 MB per CI vector).
constexpr double kMaxFCIDeterminants = 2.0e7;

// Size, memory and time of a Full CI with n orbitals (N basis functions) and the
// given electrons; `rhf_pass_seconds` is the time of one pass over all two-electron
// integrals (from the RHF estimate). Sets feasible = false with a reason above the
// determinant limit or the memory limit.
FCISize estimate_fci(int n_orbitals, int n_basis, int n_alpha, int n_beta, int threads, double memory_limit_mb,
                     double eri_pass_seconds);

struct FCIOptions {
    int threads = 0;
    double memory_limit_mb = 4096.0;
    int max_iter = 100;
    double residual_tol = 1.0e-6;
    double energy_tol = 1.0e-10;
};

struct FCIIteration {
    int iteration = 0;
    double energy = 0.0;     // total energy (E_h)
    double residual = 0.0;   // norm of H x - E x
    int subspace = 0;
};

struct FCIDeterminant {
    double coefficient = 0.0;
    std::vector<int> alpha, beta;  // occupied MOs, 1-based
    int excitation = 0;            // relative to the RHF determinant
};

struct FCIResultGeneral {
    bool converged = false;
    int iterations = 0;
    double electronic_energy = 0.0;
    double total_energy = 0.0;
    double correlation_energy = 0.0;  // E_FCI - E_RHF
    double s2 = 0.0;                  // <S^2>
    double hf_weight = 0.0;           // |C_RHF|^2
    int orbitals = 0, alpha = 0, beta = 0;
    std::size_t alpha_strings = 0, beta_strings = 0, determinants = 0;
    std::vector<FCIDeterminant> leading;          // largest |C| first
    std::vector<double> natural_occupations;      // descending, sum = number of electrons
    std::vector<FCIIteration> history;
    double transform_seconds = 0.0, ci_seconds = 0.0;
    bool dense_integrals = false;  // (ij|kl) as an n^2 x n^2 matrix for the matrix products
};

using FCIProgress = std::function<void(const FCIIteration&)>;

// Full CI in the MO basis (electronic Hamiltonian from `mo`), lowest singlet root.
FCIResultGeneral run_fci_determinants(const MOIntegrals& mo, int n_alpha, int n_beta, double nuclear_repulsion,
                                      double rhf_energy, const FCIOptions& options,
                                      const FCIProgress& progress = FCIProgress());

}  // namespace macqc
