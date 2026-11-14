// Redundant internal coordinates (ΨLab 5.0): bond stretches, bends and linear bends, the
// Wilson B matrix, the generalised inverse of G = B B', internal gradients and Hessians, and
// the iterative back-transformation of internal steps to Cartesian coordinates
// (Pulay & Fogarasi, J. Chem. Phys. 96, 2856 (1992); Peng, Ayala, Schlegel & Frisch,
// J. Comput. Chem. 17, 49 (1996)).
//
// Primitives: every pair of atoms closer than 1.3 times the sum of their covalent radii is a
// bond; every pair of bonds at a common atom is a bend, replaced by two linear-bend
// components (in two perpendicular planes through the A...C axis) when the angle exceeds 175
// degrees. Units: bohr for stretches, radians for bends; the linear-bend components are the
// projections (u + v) . e_k of the unit vectors from the central atom (dimensionless, close to
// the deviation from 180 degrees in radians).
#pragma once

#include <array>
#include <string>
#include <vector>

#include "macqc/linalg.h"
#include "macqc/molecule.h"

namespace macqc {

enum class InternalType { Stretch, Bend, LinearBend };

struct InternalCoordinate {
    InternalType type = InternalType::Stretch;
    std::vector<int> atoms;          // 2 (stretch) or 3 (bend: A, centre, C)
    int component = 0;               // linear bend: 0 or 1
    std::array<double, 3> axis{};    // linear bend: the fixed direction e_k
    std::string name() const;        // e.g. "R1", "A1", "L1a"
    std::string definition(const Molecule& mol) const;  // e.g. "R(O1-H2)", "A(H2-O1-H3)"
};

class InternalCoordinates {
public:
    InternalCoordinates() = default;
    // Builds the primitives from the connectivity of `mol`.
    explicit InternalCoordinates(const Molecule& mol);

    const std::vector<InternalCoordinate>& coordinates() const { return coords_; }
    std::size_t size() const { return coords_.size(); }
    // True if any bend became linear (> 175 deg) or a linear bend bent back (< 170 deg): the
    // set must be rebuilt for this geometry.
    bool needs_rebuild(const Molecule& mol) const;

    std::vector<double> values(const Molecule& mol) const;
    Matrix wilson_b(const Molecule& mol) const;  // n_q x 3 n_atoms
    // Second derivatives d2q_i/dx dx' contracted with c_i: sum_i c_i d2q_i (3N x 3N),
    // by central differences of the B matrix.
    Matrix b_derivative_contracted(const Molecule& mol, const std::vector<double>& c) const;
    // q1 - q0 with angles kept in range.
    std::vector<double> difference(const std::vector<double>& q1, const std::vector<double>& q0) const;

private:
    std::vector<InternalCoordinate> coords_;
};

// Generalised inverse of a symmetric positive semidefinite matrix (eigenvalues below
// tol * largest are treated as zero); `rank` receives the number of kept eigenvalues.
Matrix generalized_inverse(const Matrix& G, double tol = 1.0e-8, int* rank = nullptr);

// Internal-coordinate gradient g_q = G^- B g_x.
std::vector<double> internal_gradient(const Matrix& B, const Matrix& Ginv, const std::vector<double>& gx);

struct BackTransform {
    Molecule molecule;
    int iterations = 0;
    double residual = 0.0;  // largest |q(x) - q_target| at the end
    bool converged = false;
};

// Cartesian geometry whose internal coordinates are q0 + dq (iterative: x += B' G^- (q_t - q(x)),
// first iteration exact to first order). Stops when the Cartesian change is below 1e-10 bohr.
BackTransform back_transform(const InternalCoordinates& ic, const Molecule& mol, const std::vector<double>& dq,
                             int max_iter = 50);

// Covalent radius (bohr) used for the connectivity (Pyykkö-Atsumi single-bond radii).
double covalent_radius(int Z);

}  // namespace macqc
