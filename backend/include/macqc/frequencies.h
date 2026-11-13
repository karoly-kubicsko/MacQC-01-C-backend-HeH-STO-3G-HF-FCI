// Harmonic vibrational analysis and ideal-gas thermochemistry (ΨLab 5.0).
//
// Normal modes: the Cartesian Hessian is mass-weighted with the masses of the most abundant
// isotopes (AME 2020); translations and rotations about the centre of mass (Eckart conditions)
// are removed by working in the orthonormal complement of their 6 (linear: 5) mass-weighted
// vectors; the eigenvalues give the harmonic wavenumbers (negative = imaginary). For each
// mode: wavenumber (cm^-1), reduced mass (u), force constant (mdyn/angstrom), IR intensity
// (km/mol, from the dipole derivatives dmu/dQ), normalised Cartesian displacements, and for
// C2v H2O and diatomics the symmetry species.
//
// Thermochemistry (CODATA 2018): ideal-gas translation at the given pressure, rigid rotor
// (linear or nonlinear, symmetry number sigma), harmonic oscillators with the zero-point
// energy, closed-shell singlet ground state (q_el = 1). Imaginary modes and modes below the
// cut-off (10 cm^-1) are left out of the vibrational sums, with a warning, as in Gaussian.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "macqc/linalg.h"
#include "macqc/molecule.h"

namespace macqc {

struct NormalMode {
    double wavenumber = 0.0;      // cm^-1; negative for an imaginary frequency
    bool imaginary = false;
    double reduced_mass = 0.0;    // u
    double force_constant = 0.0;  // mdyn/angstrom (negative for imaginary modes)
    double ir_intensity = -1.0;   // km/mol (negative: not available)
    std::array<double, 3> dmu_dq{};  // dipole derivative along the mass-weighted mode (e / u^1/2)
    std::vector<std::array<double, 3>> displacement;  // normalised Cartesian displacement per atom
    std::string symmetry;         // "A1", "B2", "SG+", "S+" or ""
};

struct HarmonicAnalysis {
    std::vector<double> masses;           // u
    std::vector<std::string> isotopes;    // e.g. "16O"
    bool linear = false;
    int n_vib = 0;
    std::vector<NormalMode> modes;        // ascending wavenumber (imaginary first)
    std::vector<double> low_frequencies;  // the translations, rotations (unprojected) and lowest modes, cm^-1
    bool has_intensities = false;
    std::array<double, 3> center_of_mass{};  // bohr
};

// `dipole_derivatives`: atomic polar tensor, [3 n_atoms][xyz] (e), or nullptr.
HarmonicAnalysis harmonic_analysis(const Molecule& mol, const Matrix& hessian,
                                   const std::vector<std::array<double, 3>>* dipole_derivatives);

struct ThermoTerm {
    double energy = 0.0;   // J/mol (thermal energy; vibration includes the zero-point energy)
    double entropy = 0.0;  // J/(mol K)
    double cv = 0.0;       // J/(mol K)
    double q = 1.0;        // partition function (vibration: from the bottom of the well)
};

struct ThermoResult {
    double temperature = 298.15;  // K
    double pressure = 101325.0;   // Pa
    int symmetry_number = 1;
    bool linear = false;
    double mass = 0.0;            // u (molecule)
    std::array<double, 3> moments{};          // principal moments of inertia, u bohr^2 (ascending)
    std::array<double, 3> rotational_constants_ghz{}, rotational_constants_cm{}, rotational_temperatures{};
    std::vector<double> vibrational_temperatures;  // K, of the modes used
    std::vector<std::string> warnings;
    int excluded_modes = 0;
    double zero_point_energy = 0.0;  // J/mol
    ThermoTerm translation, rotation, vibration, electronic, total;
    double electronic_energy = 0.0;  // E_h (the energy of the method at the geometry)
    // corrections (J/mol) and totals (E_h)
    double thermal_energy = 0.0, enthalpy = 0.0, gibbs = 0.0, helmholtz = 0.0;  // corrections, J/mol
    double total_u = 0.0, total_h = 0.0, total_g = 0.0, total_a = 0.0, total_zpe = 0.0;  // E_h
};

ThermoResult thermochemistry(const Molecule& mol, const HarmonicAnalysis& h, double electronic_energy,
                             double temperature = 298.15, double pressure = 101325.0, int symmetry_number = 1,
                             double low_frequency_cutoff = 10.0);

// Symmetry number of ΨLab's molecules: homonuclear diatomic 2, heteronuclear 1, C2v H2O 2.
int symmetry_number(const Molecule& mol);

// J/mol -> E_h, kcal/mol, kJ/mol.
double joule_per_mol_to_hartree(double x);

}  // namespace macqc
