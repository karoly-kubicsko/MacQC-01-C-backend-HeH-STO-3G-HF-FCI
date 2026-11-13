// Analytic first and second derivatives of the closed-shell RHF energy with respect to the
// nuclear coordinates (ΨLab 5.0), and the dipole moment and its derivatives.
//
// Gradient (E_h/bohr), with W = 2 sum_i e_i C_i C_i^T the energy-weighted density:
//   dE/dx = sum P h^x + sum Gamma (mn|ls)^x - sum W S^x + dV_nn/dx (+ lindep term),
// Gamma = P P / 2 - P P / 4 (Coulomb and exchange).
// When near-linearly dependent combinations of basis functions are removed (canonical
// orthogonalisation), the orbital space rotates with the geometry, which adds the exact term
//   2 sum_{l removed, k kept} (e_l' S^x e_k)(e_l' F P S e_k) / (s_k (s_k - s_l))
// (e, s: eigenvectors and eigenvalues of S). It makes the gradient the derivative of the
// energy that is actually computed, as long as the number of removed functions does not change.
//
// Hessian: derivative integrals up to second order and the coupled-perturbed HF equations
// (CPHF) for the orbital response U^x of all 3 n_atoms coordinates (Pople, Krishnan,
// Schlegel and Binkley, Int. J. Quantum Chem. Symp. 13, 225 (1979)):
//   E^{xy} = sum P h^{xy} + sum Gamma (mn|ls)^{xy} - sum W S^{xy} + V_nn^{xy}
//          + 4 sum_{pi} U^y_pi (F^(x)_pi - e_i S^x_pi) - 2 sum_{ij} F^y_ij S^x_ij,
// F^(x) = h^x + G^x(P) the explicit Fock derivative, U^y_ij = -S^y_ij / 2 and U^y_ai from
// (e_a - e_i) U_ai + sum_bj [4(ai|bj) - (ab|ij) - (aj|ib)] U_bj
//      = e_i S_ai - F^(x)_ai + sum_kl S_kl [2(ai|kl) - (ak|il)].
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/linalg.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"

namespace macqc {

struct RHFGradientResult {
    std::vector<double> gradient;      // 3 n_atoms, total (E_h/bohr)
    std::vector<double> nuclear;       // nuclear repulsion
    std::vector<double> one_electron;  // sum P h^x (kinetic + nuclear attraction, incl. operator)
    std::vector<double> two_electron;  // sum Gamma (mn|ls)^x
    std::vector<double> overlap;       // - sum W S^x
    std::vector<double> lindep;        // orbital-space rotation term (zero if nothing removed)
    double seconds = 0.0;
};

// `lindep_tol` must be the threshold used by the SCF (RHFOptions::lindep_tol).
RHFGradientResult rhf_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf, int threads = 0,
                               double lindep_tol = 1.0e-7);

// Electronic + nuclear dipole moment (e bohr) of a density P about `origin`.
std::array<double, 3> dipole_moment(const Molecule& mol, const BasisSet& basis, const Matrix& P,
                                    const std::array<double, 3>& origin = {0.0, 0.0, 0.0});

struct CPHFProgress {
    int iteration = 0;
    double residual = 0.0;  // largest residual norm over the perturbations
};

struct RHFHessianResult {
    Matrix hessian;                                  // 3 n_atoms x 3 n_atoms (E_h/bohr^2)
    std::vector<std::array<double, 3>> dipole_derivatives;  // [3 n_atoms][xyz] (e bohr / bohr = e)
    std::array<double, 3> dipole{};                  // dipole moment (e bohr), origin 0
    bool cphf_converged = false;
    int cphf_iterations = 0;
    double cphf_residual = 0.0;
    double integral_seconds = 0.0, cphf_seconds = 0.0, seconds = 0.0;
};

using CPHFCallback = std::function<void(const CPHFProgress&)>;

// The RHF must be converged tightly (orbital gradient ~1e-8 or better). `fock` may be
// passed to reuse the stored integrals of the SCF (nullptr: a new builder).
RHFHessianResult rhf_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf, int threads = 0,
                             double memory_mb = 1024.0, double cphf_tol = 1.0e-9,
                             const CPHFCallback& progress = CPHFCallback());

// Solves the RHF orbital-response equations
//   (e_a - e_i) U_ai + sum_bj [4(ai|bj) - (ab|ij) - (aj|ib)] U_bj = B_ai
// for several right-hand sides at once (B: v x o matrices); returns U (same layout).
// Krylov subspace method with the orbital-energy-difference preconditioner.
struct CPHFSolution {
    std::vector<Matrix> U;
    bool converged = false;
    int iterations = 0;
    double residual = 0.0;
};
CPHFSolution solve_cphf(const FockBuilder& fock, const Matrix& C, const std::vector<double>& eps, int nocc,
                        const std::vector<Matrix>& rhs, double tol, int max_iter = 60,
                        const CPHFCallback& progress = CPHFCallback(), int threads = 0);

// Helpers shared with the correlated gradients: C^T M C' for row-major matrices.
Matrix transform_matrix(const Matrix& left, const Matrix& M, const Matrix& right, int threads = 0);

}  // namespace macqc
