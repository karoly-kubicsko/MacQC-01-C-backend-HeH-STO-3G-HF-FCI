// One- and two-electron integrals in the molecular-orbital basis, by an
// integral-direct four-index transformation from the AO shell quartets.
#pragma once

#include <cstddef>
#include <vector>

#include "macqc/basis.h"
#include "macqc/linalg.h"

namespace macqc {

struct MOIntegrals {
    int n = 0;                 // number of molecular orbitals
    Matrix h;                  // core Hamiltonian, n x n
    std::size_t npair = 0;     // n (n + 1) / 2
    std::size_t stride = 0;    // row length of `eri`
    std::vector<double> eri;   // (pq|rs) = eri[pair(p,q) * stride + pair(r,s)]

    static std::size_t pair(std::size_t p, std::size_t q) { return p >= q ? p * (p + 1) / 2 + q : q * (q + 1) / 2 + p; }
    double operator()(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const {
        return eri[pair(p, q) * stride + pair(r, s)];
    }
};

// Memory of the transformation (bytes): the half-transformed integrals,
// n(n+1)/2 x N(N+1)/2 values for N basis functions and n orbitals.
double mo_transform_bytes(std::size_t n_basis, std::size_t n_mo);

// Transforms the AO integrals of `basis` with the MO coefficients C (N x n, columns = MOs).
// h_ao is the AO core Hamiltonian. Uses `threads` threads (0: all cores).
MOIntegrals transform_to_mo(const BasisSet& basis, const Matrix& h_ao, const Matrix& C, int threads);

}  // namespace macqc
