// Analytic gradients of correlated methods from relaxed densities (ΨLab 5.0).
//
// A method supplies its response densities in the RHF molecular orbitals: the one-particle
// density D and the two-particle density G such that the energy (Lagrangian) is
//   L = sum_pq D_pq h_pq + sum_pqrs G_pqrs (pq|rs) + V_nn,
// with every dependence on the Fock matrix and the orbital energies written through h and
// (pq|rs) (f = h + 2J - K of the reference). The method must be invariant to rotations among
// the occupied and among the virtual orbitals (CCSD, CCSD(T) with the full Fock matrix,
// CCSDT, CCSDTQ, Full CI are). The gradient is then
//   dE/dx = sum D h^(x) + sum G (pq|rs)^(x) - sum W S^x + dV_nn/dx
// with the orbital response of the canonical RHF orbitals handled by one Z-vector (CPHF)
// equation:
//   X_tp = 2 sum_q h_tq D_qp + 4 sum_qrs (tq|rs) G_pqrs        (orbital Lagrangian)
//   A Z = X_ai - X_ia,  A the RHF orbital Hessian of solve_cphf,
// which adds -Z_ai F^(x)_ai to the explicit terms and Z, X to the energy-weighted density.
// The relaxed density C (D_sym) C' - C_v Z C_o' (symmetrised) is also the density whose
// dipole moment is the derivative of the energy with respect to an electric field.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/linalg.h"
#include "macqc/mo_integrals.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"
#include "macqc/rhf_derivatives.h"

namespace macqc {

// Symmetric 4-index array (8-fold permutational symmetry of real (pq|rs)) over n indices,
// stored once per class: v[pair(pair(p,q), pair(r,s))]. get() returns the value of every
// index order; add_raw(p,q,r,s,x) accumulates a component of a not yet symmetrised tensor
// (the stored value becomes the average over the equivalent orders).
struct Packed4 {
    std::size_t n = 0;
    std::vector<double> v;

    Packed4() = default;
    explicit Packed4(std::size_t n_) : n(n_), v(size_for(n_), 0.0) {}
    static std::size_t pair(std::size_t a, std::size_t b) { return a >= b ? a * (a + 1) / 2 + b : b * (b + 1) / 2 + a; }
    static std::size_t size_for(std::size_t n_) {
        const std::size_t m = n_ * (n_ + 1) / 2;
        return m * (m + 1) / 2;
    }
    static std::size_t index(std::size_t p, std::size_t q, std::size_t r, std::size_t s) {
        return pair(pair(p, q), pair(r, s));
    }
    static double degeneracy(std::size_t p, std::size_t q, std::size_t r, std::size_t s) {
        const std::size_t pq = pair(p, q), rs = pair(r, s);
        return (p == q ? 1.0 : 2.0) * (r == s ? 1.0 : 2.0) * (pq == rs ? 1.0 : 2.0);
    }
    double get(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const { return v[index(p, q, r, s)]; }
    void add_raw(std::size_t p, std::size_t q, std::size_t r, std::size_t s, double x) {
        v[index(p, q, r, s)] += x / degeneracy(p, q, r, s);
    }
    double bytes() const { return 8.0 * static_cast<double>(v.size()); }
};

struct CorrelatedDensity {
    Matrix D;   // n x n (MO), E_1 = sum D h; symmetrised by the gradient code
    Packed4 G;  // MO two-particle density, E_2 = sum_pqrs G (pq|rs)
};

// E = sum D h + sum G (pq|rs) + V_nn: checks that the densities reproduce the energy.
double density_energy(const CorrelatedDensity& d, const MOIntegrals& mo, double nuclear_repulsion);

struct CorrelatedGradientResult {
    std::vector<double> gradient;      // 3 n_atoms, total (E_h/bohr)
    std::vector<double> nuclear, one_electron, two_electron, overlap;
    Matrix relaxed_density_ao;         // N x N, relaxed one-particle density (AO)
    std::array<double, 3> relaxed_dipole{};  // e bohr, origin 0
    double orbital_lagrangian_oo_vv = 0.0;   // largest |X_pq - X_qp| in the occ-occ and vir-vir blocks
    bool zvector_converged = false;
    int zvector_iterations = 0;
    double zvector_residual = 0.0;
    double density_energy = 0.0;       // sum D h + sum G g + V_nn (equals the method's energy)
    double seconds = 0.0, lagrangian_seconds = 0.0, zvector_seconds = 0.0, backtransform_seconds = 0.0,
           integral_seconds = 0.0;
};

// `mo` must hold the MO integrals of all n = rhf.C.cols() orbitals. Throws if functions were
// removed (the rotation of the kept orbital space is not included here; the caller then
// differentiates numerically).
CorrelatedGradientResult correlated_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                             const MOIntegrals& mo, const CorrelatedDensity& density,
                                             int threads = 0, double memory_mb = 1024.0,
                                             const CPHFCallback& progress = CPHFCallback());

// Memory (bytes) of an analytic correlated gradient with n orbitals and N basis functions:
// MO integrals, MO and AO two-particle densities and the half-transformed density.
double correlated_gradient_bytes(std::size_t n_mo, std::size_t n_basis);

// The reference (RHF) densities in this form (tests: they must give the RHF gradient).
CorrelatedDensity rhf_density(int n_mo, int n_occ);

}  // namespace macqc
