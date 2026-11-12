#include "macqc/molecule.h"

#include <cctype>
#include <cmath>
#include <stdexcept>

namespace macqc {

namespace {

const char* const kSymbols[] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"};

std::string compact_lower(const std::string& s) {
    std::string out;
    for (char c : s)
        if (!std::isspace(static_cast<unsigned char>(c)) && c != '_' && c != '-')
            out += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

}  // namespace

std::string element_symbol(int Z) {
    if (Z < 1 || Z > 10) throw std::invalid_argument("element Z = " + std::to_string(Z) + " not supported");
    return kSymbols[Z - 1];
}

int atomic_number(const std::string& symbol) {
    const std::string s = compact_lower(symbol);
    for (int z = 1; z <= 10; ++z)
        if (s == compact_lower(kSymbols[z - 1])) return z;
    return 0;
}

std::string canonical_molecule_name(const std::string& name) {
    const std::string s = compact_lower(name);
    if (s == "h2") return "H2";
    if (s == "lih") return "LiH";
    if (s == "h2o" || s == "water") return "H2O";
    if (s == "heh+" || s == "heh") return "HeH+";
    return "";
}

double default_bond_length(const std::string& molecule) {
    const std::string m = canonical_molecule_name(molecule);
    if (m == "H2") return 0.74144 / kBohrInAngstrom;
    if (m == "LiH") return 1.5949 / kBohrInAngstrom;
    if (m == "H2O") return 0.9572 / kBohrInAngstrom;
    if (m == "HeH+") return 1.4632;
    throw std::invalid_argument("unknown molecule '" + molecule + "'");
}

int Molecule::electrons() const {
    int n = -charge;
    for (const Atom& a : atoms) n += a.Z;
    return n;
}

double Molecule::nuclear_repulsion() const {
    double e = 0.0;
    for (std::size_t i = 0; i < atoms.size(); ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const double dx = atoms[i].position[0] - atoms[j].position[0];
            const double dy = atoms[i].position[1] - atoms[j].position[1];
            const double dz = atoms[i].position[2] - atoms[j].position[2];
            e += atoms[i].Z * atoms[j].Z / std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    return e;
}

Molecule make_molecule(const std::string& name, double R, double angle_degrees) {
    const std::string m = canonical_molecule_name(name);
    if (!(std::isfinite(R) && R > 0.0))
        throw std::invalid_argument("the bond length must be a positive, finite number");
    if (R < 0.3) throw std::invalid_argument("the bond length must be at least 0.3 bohr (0.16 \xC3\x85)");
    if (R > 50.0) throw std::invalid_argument("the bond length must be at most 50 bohr (26.5 \xC3\x85)");
    Molecule mol;
    mol.name = m;
    auto atom = [](int Z, double x, double y, double z) {
        Atom a;
        a.Z = Z;
        a.symbol = element_symbol(Z);
        a.position = {x, y, z};
        return a;
    };
    if (m == "H2") {
        mol.atoms = {atom(1, 0.0, 0.0, -0.5 * R), atom(1, 0.0, 0.0, 0.5 * R)};
    } else if (m == "LiH") {
        mol.atoms = {atom(3, 0.0, 0.0, 0.0), atom(1, 0.0, 0.0, R)};
    } else if (m == "H2O") {
        if (!(std::isfinite(angle_degrees) && angle_degrees >= 30.0 && angle_degrees <= 180.0))
            throw std::invalid_argument("the H-O-H angle must be between 30 and 180 degrees");
        const double half = 0.5 * angle_degrees * 3.14159265358979323846 / 180.0;
        const double x = R * std::sin(half), z = R * std::cos(half);
        mol.atoms = {atom(8, 0.0, 0.0, 0.0), atom(1, x, 0.0, z), atom(1, -x, 0.0, z)};
    } else {
        throw std::invalid_argument("unknown molecule '" + name + "' (H2, LiH or H2O)");
    }
    return mol;
}

}  // namespace macqc
