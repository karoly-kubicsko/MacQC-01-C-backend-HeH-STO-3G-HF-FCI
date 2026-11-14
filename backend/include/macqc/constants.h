// Physical constants (CODATA 2018) and isotopic masses (AME 2020) for harmonic frequencies,
// IR intensities and thermochemistry (ΨLab 5.0).
//
// The geometry input keeps the bohr/angstrom factor of version 3.0 (kBohrInAngstrom, CODATA
// 2022) so that every earlier result stays identical; the two values differ by 7e-10
// relative, far below any printed digit of a frequency or a thermodynamic quantity.
#pragma once

namespace macqc {
namespace codata2018 {

constexpr double kSpeedOfLight = 299792458.0;          // m/s (exact)
constexpr double kPlanck = 6.62607015e-34;             // J s (exact)
constexpr double kBoltzmann = 1.380649e-23;            // J/K (exact)
constexpr double kAvogadro = 6.02214076e23;            // 1/mol (exact)
constexpr double kElementaryCharge = 1.602176634e-19;  // C (exact)
constexpr double kHartree = 4.3597447222071e-18;       // J
constexpr double kBohr = 5.29177210903e-11;            // m
constexpr double kAtomicMassUnit = 1.66053906660e-27;  // kg
constexpr double kVacuumPermittivity = 8.8541878128e-12;  // F/m
constexpr double kCalorie = 4.184;                     // J (thermochemical calorie)
constexpr double kGasConstant = kAvogadro * kBoltzmann;  // J/(mol K)
constexpr double kPi = 3.14159265358979323846;

// 1 E_h in kJ/mol and kcal/mol.
constexpr double kHartreeKJPerMol = kHartree * kAvogadro / 1000.0;
constexpr double kHartreeKcalPerMol = kHartreeKJPerMol / kCalorie;
// Debye in e bohr: 1 D = 1e-21 / c  C m.
constexpr double kDebyeInAtomicUnits = 1.0e-21 / kSpeedOfLight / (kElementaryCharge * kBohr);

}  // namespace codata2018

// Mass (u) of the most abundant isotope of element Z = 1 .. 10 (AME 2020: Wang, Huang, Kondev,
// Audi & Naimi, Chin. Phys. C 45, 030003 (2021)); ¹H, ⁴He, ⁷Li, ⁹Be, ¹¹B, ¹²C, ¹⁴N, ¹⁶O, ¹⁹F, ²⁰Ne.
double isotope_mass(int Z);
// Its label, e.g. "1H", "7Li", "16O".
const char* isotope_label(int Z);

}  // namespace macqc
