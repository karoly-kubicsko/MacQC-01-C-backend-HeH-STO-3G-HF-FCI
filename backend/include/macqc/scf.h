// Restricted Hartree-Fock for a two-electron, two-basis-function diatomic.
// Ported from the Szabo & Ostlund Appendix B program (routines HFCALC,
// COLLECT, FORMG, DIAG, SCF, MATOUT).
#pragma once

#include <string>
#include <vector>

#include "macqc/integrals.h"
#include "macqc/linalg.h"

namespace macqc {

struct HFOptions {
    bool exact_erf = false;     // use std::erf instead of derfother
    int max_iter = 25;
    double convergence = 1.0e-4;
    int print_level = 2;        // FORTRAN IOP: 0 silent, 1 summary, 2 every iteration
};

struct SCFIteration {
    int iteration = 0;
    double electronic_energy = 0.0;
    double delta_p = 0.0;
};

struct HFResult {
    bool converged = false;
    int iterations = 0;
    double electronic_energy = 0.0;
    double nuclear_repulsion = 0.0;
    double total_energy = 0.0;
    std::vector<double> orbital_energies;
    Matrix C, P, S, X, H, F, G, PS;
    Tensor4 eri;  // (ij|kl), chemists' notation, AO basis
    Integrals integrals;
    std::vector<SCFIteration> history;
    std::string report;  // FORTRAN-style printed output
};

class TwoElectronDiatomicHF {
public:
    TwoElectronDiatomicHF(const DiatomicBasis& basis, const HFOptions& options);

    void calculate_integrals();
    void collect_matrices();
    void form_g();
    static void diag_2x2(const Matrix& f, Matrix& c, Matrix& e);
    void scf();

    // calculate_integrals + collect_matrices + scf
    HFResult run();

private:
    DiatomicBasis basis_;
    HFOptions opt_;
    Integrals ints_;
    Matrix S_, X_, XT_, H_, F_, G_, C_, Fprime_, Cprime_, P_, oldP_, E_;
    Tensor4 TT_;
    std::string out_;
    bool converged_ = false;
    int iterations_ = 0;
    double electronic_energy_ = 0.0, total_energy_ = 0.0;
    Matrix PS_;
    std::vector<SCFIteration> history_;
};

}  // namespace macqc
