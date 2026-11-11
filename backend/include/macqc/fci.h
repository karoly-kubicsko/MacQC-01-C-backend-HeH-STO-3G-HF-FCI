// Full CI for two electrons (M_s = 0) in the converged RHF molecular-orbital basis.
#pragma once

#include <string>
#include <vector>

#include "macqc/linalg.h"

namespace macqc {

struct FCIState {
    double energy = 0.0;  // total energy (electronic + nuclear repulsion)
    double S2 = 0.0;      // <S^2>
    std::string spin;     // "singlet", "triplet", ...
};

struct SingletWeight {
    std::string configuration;  // e.g. "(1sigma)^2"
    double weight = 0.0;
};

struct FCIResult {
    double electronic_energy = 0.0;  // ground state, electronic only
    double total_energy = 0.0;
    double correlation_energy = 0.0;  // total_energy - E_HF(total)
    std::vector<std::string> determinants;  // labels, e.g. "|1a 1b>" (a = alpha, b = beta)
    std::vector<FCIState> states;           // ascending energy
    std::vector<std::vector<double>> vectors;  // vectors[state][determinant]
    std::vector<double> ground_state_vector;
    std::vector<SingletWeight> singlet_weights;  // of the ground state
    Matrix h_mo;         // core Hamiltonian in the MO basis
    Tensor4 eri_mo;      // (pq|rs) in the MO basis
    Matrix hamiltonian;  // CI Hamiltonian (electronic energies), determinant basis
};

// h_ao, eri_ao: core Hamiltonian and (ij|kl) in the AO basis; C: MO coefficients
// (columns, AO x MO); the lowest MO is doubly occupied in the HF reference.
FCIResult run_fci(const Matrix& h_ao, const Tensor4& eri_ao, const Matrix& C,
                  double nuclear_repulsion, double hf_total_energy);

// FORTRAN-style text section appended after the HF report (empty at print_level 0).
std::string fci_report(const FCIResult& r, int print_level);

}  // namespace macqc
