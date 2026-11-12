// Molecules of ΨLab 3.0: H2, LiH and H2O built from bond length (and angle).
#pragma once

#include <array>
#include <string>
#include <vector>

namespace macqc {

// 1 bohr in angstrom (CODATA 2022).
constexpr double kBohrInAngstrom = 0.529177210544;

struct Atom {
    int Z = 0;
    std::string symbol;
    std::array<double, 3> position{};  // bohr
};

struct Molecule {
    std::string name;  // "H2", "LiH", "H2O"
    std::vector<Atom> atoms;
    int charge = 0;

    int electrons() const;
    double nuclear_repulsion() const;
};

// Element symbol of Z = 1..10 ("H", "He", "Li", ..., "Ne"); throws otherwise.
std::string element_symbol(int Z);
// Atomic number of an element symbol (case-insensitive), 0 if unknown.
int atomic_number(const std::string& symbol);

// Canonical molecule name ("H2", "LiH", "H2O", "HeH+") for any spelling, or "".
std::string canonical_molecule_name(const std::string& name);

// Default geometries (experimental equilibrium structures), in bohr / degrees:
// H2 r_e = 0.74144 Å; LiH r_e = 1.5949 Å; H2O r(OH) = 0.9572 Å, angle 104.52°.
double default_bond_length(const std::string& molecule);  // bohr
constexpr double kDefaultWaterAngle = 104.52;              // degrees

// Builds the molecule. H2: atoms at z = -R/2, +R/2. LiH: Li at the origin, H at z = R.
// H2O: O at the origin, the hydrogens in the xz plane at distance R, symmetric about
// the z axis, enclosing `angle_degrees`. Throws std::invalid_argument on bad input.
Molecule make_molecule(const std::string& name, double bond_length_bohr, double angle_degrees);

}  // namespace macqc
