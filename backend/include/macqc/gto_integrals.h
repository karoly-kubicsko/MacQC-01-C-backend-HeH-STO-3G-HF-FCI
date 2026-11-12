// Integrals over contracted spherical-harmonic Gaussian shells (l <= 5):
//   overlap, kinetic energy, nuclear attraction (Obara-Saika recurrences) and
//   electron repulsion (Head-Gordon-Pople: Obara-Saika vertical recurrence per
//   primitive quartet, contraction, then horizontal recurrence and the
//   Cartesian -> real-solid-harmonic transformation).
#pragma once

#include <vector>

#include "macqc/basis.h"
#include "macqc/linalg.h"
#include "macqc/molecule.h"

namespace macqc {

// Boys function F_m(T) for m = 0 .. mmax (mmax <= 4 * kMaxAngularMomentum + 2).
// Tabulated Taylor expansion (7 terms, grid 0.05) with downward recursion for
// T < 45, asymptotic F_0 with upward recursion above. Relative error ~1e-14.
void boys_function(int mmax, double T, double* F);

// Real solid harmonics in Racah normalisation, S_lm = sum_c M(m, c) x^cx y^cy z^cz,
// rows m = -l .. l, columns the Cartesian components in the order
// (l,0,0), (l-1,1,0), (l-1,0,1), (l-2,2,0), ... (0,0,l). With primitives normalised as
// x^l exp(-a r^2) the resulting spherical functions are normalised.
const std::vector<double>& solid_harmonic_matrix(int l);

Matrix overlap_matrix(const BasisSet& basis);
Matrix kinetic_matrix(const BasisSet& basis);
Matrix nuclear_attraction_matrix(const BasisSet& basis, const Molecule& mol);

// Primitive-pair data of a shell pair (bra or ket of an ERI).
struct PrimitivePair {
    double zeta = 0.0;            // a + b
    double P[3] = {0, 0, 0};      // Gaussian product centre
    double PA[3] = {0, 0, 0};     // P - A
    double K = 0.0;               // c_a c_b exp(-a b |AB|^2 / zeta)
};

struct ShellPair {
    int a = 0, b = 0;            // shell indices, l(a) >= l(b)
    double AB[3] = {0, 0, 0};    // A - B
    std::vector<PrimitivePair> prims;
};

// All shell pairs (a >= b by angular momentum, see ShellPair) of a basis.
struct ShellPairList {
    std::vector<ShellPair> pairs;       // index p = i*(i+1)/2 + j for shells i >= j
    std::vector<double> schwarz;        // sqrt(max |(ab|ab)|) per pair
    double max_schwarz = 0.0;
    static int index(int i, int j) { return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }
};

ShellPairList make_shell_pairs(const BasisSet& basis);

// Electron-repulsion integrals of one shell quartet. One engine per thread.
class EriEngine {
public:
    explicit EriEngine(const BasisSet& basis);
    // (ab|cd) in chemists' notation for the shells of `bra` and `ket`; the result has
    // the layout [a][b][c][d] (spherical functions of shells bra.a, bra.b, ket.a, ket.b).
    const double* compute(const ShellPair& bra, const ShellPair& ket);

private:
    const BasisSet& basis_;
    std::vector<double> vrr_, acc_, hrr_, tmp_, out_;
};

// Fills the Schwarz bounds of the list (needs an engine).
void compute_schwarz(const BasisSet& basis, ShellPairList& pairs);

// Number of threads used for integral work: `requested` if > 0, else the hardware
// concurrency (at least 1, at most 64).
int resolve_threads(int requested);

}  // namespace macqc
