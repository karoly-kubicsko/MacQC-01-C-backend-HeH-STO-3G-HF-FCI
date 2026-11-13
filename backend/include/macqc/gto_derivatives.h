// Derivatives of the Gaussian integrals with respect to the nuclear coordinates (ΨLab 5.0).
//
// The derivative of a contracted spherical Gaussian with respect to its centre,
//   d/dA_k [x^cx y^cy z^cz exp(-a r^2)] = 2a G(c + 1_k) - c_k G(c - 1_k),
// is a sum of "component shells": the Cartesian functions of l+1 with the contraction
// coefficients multiplied by 2a, and those of l-1 with the original coefficients, each with
// an output transform (ShellTransform) whose rows are the 3(2l+1) derivatives
// (row = k (2l+1) + m). Second derivatives use l+2 (4a^2), l (2a) and l-2 components with
// 6(2l+1) rows (kk' = xx, xy, xz, yy, yz, zz). The ordinary integral engines then evaluate the
// derivative integrals; nothing is approximated.
//
// Translational invariance: two-centre integrals depend only on A - B, three-centre nuclear
// attraction on A - C and B - C, and electron repulsion on A - D, B - D, C - D, so the
// derivatives with respect to the last centre are minus the sum of the others (for second
// derivatives the corresponding sums of blocks).
//
// Perturbations are the 3 n_atoms nuclear coordinates, index 3 atom + k (k = x, y, z).
#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gto_integrals.h"
#include "macqc/linalg.h"
#include "macqc/molecule.h"

namespace macqc {

struct DerivativeComponent {
    Shell shell;               // angular momentum l + 2, l + 1, l, l - 1 or l - 2; modified coefficients
    ShellTransform transform;  // rows: 3 (2l+1) first or 6 (2l+1) second derivatives
};

// First (order 1) or second (order 2) derivative components of a shell.
std::vector<DerivativeComponent> shell_derivative_components(const Shell& shell, int order);

// Index of the symmetric direction pair (k <= k'): xx 0, xy 1, xz 2, yy 3, yz 4, zz 5.
inline int direction_pair(int k, int kp) {
    if (k > kp) {
        const int t = k;
        k = kp;
        kp = t;
    }
    static const int idx[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
    return idx[k][kp];
}

// ---- one-electron integrals -------------------------------------------------------------

// dS/dx for every nuclear coordinate x (3 n_atoms matrices N x N).
std::vector<Matrix> overlap_derivative_matrices(const BasisSet& basis, int n_atoms, int threads = 0);
// dh/dx with h = T + V, including the derivative of the nuclear-attraction operator.
std::vector<Matrix> core_hamiltonian_derivative_matrices(const BasisSet& basis, const Molecule& mol,
                                                         int threads = 0);
// Gradient contractions sum_mn D_mn dX_mn/dx (vectors of length 3 n_atoms, D symmetric).
std::vector<double> overlap_gradient(const BasisSet& basis, int n_atoms, const Matrix& D, int threads = 0);
std::vector<double> core_hamiltonian_gradient(const BasisSet& basis, const Molecule& mol, const Matrix& D,
                                              int threads = 0);
// Second-derivative contractions sum_mn D_mn d2X_mn/dx dy (3 n_atoms x 3 n_atoms, D symmetric).
Matrix overlap_hessian(const BasisSet& basis, int n_atoms, const Matrix& D, int threads = 0);
Matrix core_hamiltonian_hessian(const BasisSet& basis, const Molecule& mol, const Matrix& D, int threads = 0);
// d<mu| r_k - O_k |nu>/dx for the dipole integrals: [3 n_atoms][3 components] N x N
// (the origin O is fixed in space).
std::vector<std::vector<Matrix>> dipole_derivative_matrices(const BasisSet& basis, int n_atoms,
                                                            const std::array<double, 3>& origin, int threads = 0);

// ---- two-electron integrals -------------------------------------------------------------

// Two-particle density in the AO basis, E_2 = sum_{mn ls} Gamma_mnls (mn|ls), given as a sum
// of separable terms
//   Gamma_mnls += coulomb/2 (A_mn B_ls + B_mn A_ls)
//               + exchange/4 (A_ml B_ns + A_ms B_nl + B_ml A_ns + B_ms A_nl)
// (A, B symmetric) plus, optionally, a non-separable part stored with 8-fold symmetry in the
// layout of FockBuilder::packed() (index pair(pair(m,n), pair(l,s)); the value multiplies
// (mn|ls) once for every one of the equivalent index orders).
// RHF: one term A = B = P with coulomb 1/2 and exchange -1/4.
struct SeparableDensity {
    Matrix A, B;
    double coulomb = 0.0, exchange = 0.0;
};

struct TwoParticleDensity {
    std::vector<SeparableDensity> separable;
    const std::vector<double>* packed = nullptr;  // optional non-separable part
};

// sum Gamma (mn|ls)^x for every nuclear coordinate (3 n_atoms).
std::vector<double> eri_gradient(const BasisSet& basis, int n_atoms, const TwoParticleDensity& gamma,
                                 int threads = 0);
// sum Gamma (mn|ls)^{xy} (3 n_atoms x 3 n_atoms).
Matrix eri_hessian(const BasisSet& basis, int n_atoms, const TwoParticleDensity& gamma, int threads = 0);
// G^x(D)_mn = sum_ls D_ls [(mn|ls)^x - (ml|ns)^x / 2] for every nuclear coordinate (D symmetric):
// the explicit derivative of the two-electron part of the Fock matrix.
std::vector<Matrix> eri_derivative_fock(const BasisSet& basis, int n_atoms, const Matrix& D, int threads = 0);

// Value of the energy-like contraction sum Gamma (mn|ls) (for tests).
double eri_contract(const BasisSet& basis, const TwoParticleDensity& gamma, int threads = 0);

// The basis set moved to the atom positions of `mol` (same shells, new centres).
BasisSet basis_at_geometry(const BasisSet& basis, const Molecule& mol);

// Nuclear repulsion gradient and Hessian (3 n_atoms, 3 n_atoms x 3 n_atoms).
std::vector<double> nuclear_repulsion_gradient(const Molecule& mol);
Matrix nuclear_repulsion_hessian(const Molecule& mol);

}  // namespace macqc
