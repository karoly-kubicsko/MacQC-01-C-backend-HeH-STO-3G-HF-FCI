// Integrals over contracted spherical-harmonic Gaussian shells (l <= 5):
//   overlap, kinetic energy, nuclear attraction (Obara-Saika recurrences) and
//   electron repulsion (Head-Gordon-Pople: Obara-Saika vertical recurrence per
//   primitive quartet, contraction, then horizontal recurrence and the
//   Cartesian -> real-solid-harmonic transformation).
#pragma once

#include <array>
#include <utility>
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

// Linear map from the Cartesian components of a shell (angular momentum l, order as in
// solid_harmonic_matrix) to output functions: out_r = sum_c coefficient * cart_c over the
// sparse entries of row r. The default for an ordinary shell is the real solid harmonics
// (2l+1 rows); derivative integrals use other maps (e.g. the 3(2l+1) x-, y-, z-derivatives
// of a spherical shell, expressed in the Cartesian functions of l+1 or l-1).
struct ShellTransform {
    int l = 0;
    std::vector<std::vector<std::pair<int, double>>> rows;
    int size() const { return static_cast<int>(rows.size()); }
};

// The real-solid-harmonic transform of angular momentum l (2l+1 rows).
const ShellTransform& spherical_transform(int l);

// out[pre][rows][post] = T applied to in[pre][ncart(T.l)][post].
void apply_transform(const ShellTransform& T, int pre, int post, const double* in, double* out);

// Primitive-pair data of a shell pair (bra or ket of an ERI).
struct PrimitivePair {
    double zeta = 0.0;            // a + b
    double P[3] = {0, 0, 0};      // Gaussian product centre
    double PA[3] = {0, 0, 0};     // P - A
    double K = 0.0;               // c_a c_b exp(-a b |AB|^2 / zeta)
};

struct ShellPair {
    int a = 0, b = 0;            // shell indices, l(a) >= l(b)
    int la = 0, lb = 0;          // angular momenta of a and b
    double AB[3] = {0, 0, 0};    // A - B
    std::vector<PrimitivePair> prims;
    // Output transforms of a and b (nullptr: real solid harmonics of la, lb).
    const ShellTransform* ta = nullptr;
    const ShellTransform* tb = nullptr;
    // True if a and b were exchanged to have la >= lb and the integrals are to be returned
    // in the original order [b][a] (pairs from make_shell_pair with keep_order).
    bool unswap = false;
    int rows_a() const { return ta ? ta->size() : 2 * la + 1; }
    int rows_b() const { return tb ? tb->size() : 2 * lb + 1; }
};

// One shell pair of two (possibly modified) shells X and Y with output transforms tx, ty
// (nullptr: spherical). The engine needs l(a) >= l(b); if l(Y) > l(X) the shells are
// exchanged internally and, with keep_order, the integrals are still returned as [X][Y].
ShellPair make_shell_pair(const Shell& X, const ShellTransform* tx, const Shell& Y, const ShellTransform* ty,
                          bool keep_order = true, int index_x = 0, int index_y = 0);

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
    EriEngine() = default;
    explicit EriEngine(const BasisSet&) {}
    // (ab|cd) in chemists' notation for the shells of `bra` and `ket`; the result has
    // the layout [a][b][c][d] (spherical functions of shells bra.a, bra.b, ket.a, ket.b,
    // or the rows of the pairs' transforms; [b][a] for a pair with unswap).
    const double* compute(const ShellPair& bra, const ShellPair& ket);

private:
    std::vector<double> vrr_, acc_, hrr_, tmp_, out_;
};

// Cartesian one-electron blocks [ncart(A.l)][ncart(B.l)] of a shell pair (no transform):
// overlap, kinetic energy, nuclear attraction of the point charges `nuclei` (Z and position;
// -sum_C Z_C / |r - C|), and the first moments (x - O_x), (y - O_y), (z - O_z) about the
// point O ([3][ncart(A.l)][ncart(B.l)]).
void overlap_cartesian(const Shell& A, const Shell& B, std::vector<double>& out);
void kinetic_cartesian(const Shell& A, const Shell& B, std::vector<double>& out);
void nuclear_cartesian(const Shell& A, const Shell& B, const std::vector<Atom>& nuclei, std::vector<double>& out);
void dipole_cartesian(const Shell& A, const Shell& B, const std::array<double, 3>& origin, std::vector<double>& out);

// Dipole integrals <mu| r_k - O_k |nu>, k = x, y, z (three N x N matrices).
std::vector<Matrix> dipole_matrices(const BasisSet& basis, const std::array<double, 3>& origin);

// Fills the Schwarz bounds of the list (needs an engine).
void compute_schwarz(const BasisSet& basis, ShellPairList& pairs);

// Number of threads used for integral work: `requested` if > 0, else the hardware
// concurrency (at least 1, at most 64).
int resolve_threads(int requested);

}  // namespace macqc
