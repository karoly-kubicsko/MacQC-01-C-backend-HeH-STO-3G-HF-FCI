// Primitive and contracted integrals for the two-centre, s-type STO-nG problem.
// Ported from CALCULATE_INTEGRALS and its helpers in Szabo & Ostlund, App. B.
#pragma once

namespace macqc {

// The FORTRAN DATA constant (deliberately not M_PI) so printed output matches.
extern const double kPi;

// Error function used by the original program (Abramowitz & Stegun 7.1.26).
double derfother(double arg);

// Boys function F0. With exact_erf the std::erf is used instead of derfother.
double boys_f0(double arg, bool exact_erf);

double overlap_primitive(double a, double b, double rab2);
double kinetic_primitive(double a, double b, double rab2);
double nuclear_attraction_primitive(double a, double b, double rab2, double rcp2, double zc,
                                    bool exact_erf);
double two_electron_primitive(double a, double b, double c, double d, double rab2, double rcd2,
                              double rpq2, bool exact_erf);

// Geometry/basis description for the stage-1 problem: atom A (basis function 1,
// exponent scale zeta1, charge za) and atom B (function 2, zeta2, zb), R apart.
struct DiatomicBasis {
    int n_gaussians = 3;  // STO-nG, n = 1..3
    double r = 1.4632;
    double zeta1 = 2.0925;
    double zeta2 = 1.24;
    double za = 2.0;
    double zb = 1.0;
};

// The unique integrals named exactly as in the FORTRAN program.
struct Integrals {
    double S12 = 0, T11 = 0, T12 = 0, T22 = 0;
    double V11A = 0, V12A = 0, V22A = 0;
    double V11B = 0, V12B = 0, V22B = 0;
    double V1111 = 0, V2111 = 0, V2121 = 0, V2211 = 0, V2221 = 0, V2222 = 0;
};

Integrals calculate_integrals(const DiatomicBasis& basis, bool exact_erf);

}  // namespace macqc
