// Gaussian basis sets: Gaussian94 (.gbs) parser, the basis-set library of ΨLab and
// contracted spherical-harmonic shells placed on the atoms of a molecule.
#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

#include "macqc/molecule.h"

namespace macqc {

// Highest angular momentum the integral engine supports (h functions, l = 5;
// the tables are sized for one more).
constexpr int kMaxAngularMomentum = 6;

// One contracted shell as written in a basis-set file (before normalisation).
struct ShellData {
    int l = 0;
    std::vector<double> exponents;
    std::vector<double> coefficients;  // for normalised primitives (file convention)
};

struct ElementBasis {
    int Z = 0;
    std::vector<ShellData> shells;
};

struct BasisFile {
    std::string header;                   // the '!' comment lines at the top
    std::map<int, ElementBasis> elements;  // by atomic number
};

// Parses Gaussian94 text (as written by basissetexchange.org). SP shells are split
// into an S and a P shell; a scale factor other than 1 multiplies the exponents
// by its square. Throws std::runtime_error("line N: ...") on malformed input.
BasisFile parse_gaussian94(const std::string& text);

// The basis sets offered by ΨLab, in display order.
const std::vector<std::string>& supported_basis_sets();

// Canonical display name ("aug-cc-pVDZ") for any spelling ("AUG-CC-PVDZ"), or "".
std::string canonical_basis_name(const std::string& name);

// File name in the basis directory, e.g. "aug-cc-pvdz.gbs".
std::string basis_file_name(const std::string& canonical_name);

// Default basis directory: $PSILAB_BASIS_DIR if set, else the compile-time default
// (MACQC_BASIS_DIR, the repository's basis/ folder), else "basis".
std::string default_basis_directory();

// A contracted shell of real solid-harmonic (spherical) Gaussians on an atom.
struct Shell {
    int l = 0;
    int atom = 0;                    // index into Molecule::atoms
    std::array<double, 3> center{};  // bohr
    std::vector<double> exponents;
    // Contraction coefficients that already include the normalisation of each
    // primitive (that of x^l exp(-a r^2)) and of the contracted function.
    std::vector<double> coefficients;
    int first_function = 0;  // index of the first spherical function (m = -l .. l)

    int size() const { return 2 * l + 1; }
    int cartesian_size() const { return (l + 1) * (l + 2) / 2; }
};

// Basis-set assignment of one atom, for reports.
struct AtomBasisInfo {
    std::string symbol;
    std::string basis;        // basis set actually used for this atom
    std::string contraction;  // e.g. "[4s3p2d]"
    std::string primitives;   // e.g. "(10s5p2d)" (distinct exponents per l)
    int functions = 0;
    std::string note;         // e.g. why another basis set was used for H
};

struct BasisSet {
    std::string name;  // requested basis set (canonical)
    std::vector<Shell> shells;
    std::vector<AtomBasisInfo> atoms;
    int n_functions = 0;
    int max_l = 0;
    int n_primitives = 0;  // sum over shells of exponents x (2l+1)
    std::string source;    // header line(s) of the file(s) read
};

// Builds the basis for the molecule from the files in `directory`. For H, the
// core-valence sets aug-cc-pCVXZ are not defined (H has no core); the matching
// aug-cc-pVXZ functions are used and a note is recorded. Throws std::runtime_error
// with an explanation if a basis set is not available for an element.
BasisSet build_basis(const Molecule& mol, const std::string& basis_name, const std::string& directory);

// Same, from already parsed files: `files` maps canonical basis names to contents.
BasisSet build_basis(const Molecule& mol, const std::string& basis_name,
                     const std::map<std::string, BasisFile>& files);

// Normalised contraction coefficients (primitive and contraction normalisation).
std::vector<double> normalized_coefficients(int l, const std::vector<double>& exponents,
                                            const std::vector<double>& coefficients);

// Letter of an angular momentum (0 -> 's', 5 -> 'h').
char angular_letter(int l);

}  // namespace macqc
