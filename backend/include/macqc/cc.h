// Coupled-cluster methods for closed-shell molecules (ΨLab 3.1): CCSD, CCSD(T), CCSDT,
// CCSDTQ, all electrons correlated, canonical RHF reference (the Fock matrix is used in
// full, so off-diagonal elements are allowed).
//
// Two engines:
//  * closed-shell, spin-adapted CCSD and the (T) correction (cc_closed_shell.cpp): spatial
//    orbitals, MO integrals with at most two virtual indices plus (ia|bc); the
//    particle-particle ladder sum_cd (ac|bd) tau_ij^cd is evaluated in the AO basis
//    (integrals kept in memory or recomputed), so no (vv|vv) integrals are stored;
//  * general-order spin-orbital CC (cc_general.cpp) for CCSDT and CCSDTQ (and CCSD for
//    tests): amplitudes stored once per set of distinct spin orbitals (packed,
//    antisymmetric, blocked by the number of alpha electrons), residuals from the 183
//    connected diagrams generated and verified by backend/tools/cc_diagrams.py,
//    evaluated vertex by vertex as matrix products.
#pragma once
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/linalg.h"
#include "macqc/mo_integrals.h"

namespace macqc {

enum class CCMethod { CCSD, CCSD_T, CCSDT, CCSDTQ };

// "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ".
const char* cc_method_name(CCMethod method);
// Highest excitation rank of the cluster operator: 2, 2, 3, 4.
int cc_max_rank(CCMethod method);

struct CCOptions {
    int threads = 0;                  // 0: all hardware threads
    double memory_limit_mb = 4096.0;  // amplitudes, integrals, DIIS vectors, work arrays
    double eri_memory_mb = 1024.0;    // AO integrals of the ladder kept in memory up to this size
    int max_iter = 100;
    double energy_tol = 1.0e-10;      // |E_n - E_(n-1)| (E_h)
    double amplitude_tol = 1.0e-7;    // largest |t_new - t| of the Jacobi step
    int diis_size = 8;
};

struct CCIteration {
    int iteration = 0;
    double energy = 0.0;       // total CC energy (E_h)
    double correlation = 0.0;  // correlation energy (E_h)
    double delta_e = 0.0;
    double residual = 0.0;     // largest |t_new - t|
};

struct CCAmplitude {
    double value = 0.0;
    int rank = 0;
    std::vector<int> occupied, virtuals;  // spatial MOs, 1-based, in excitation order
    std::string spins;                    // spin of each occupied (= virtual) spin orbital, 'a' or 'b'
};

struct CCResult {
    std::string method;                 // "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ"
    std::string engine;                 // which implementation
    bool converged = false;
    int iterations = 0;
    int occupied = 0, virtuals = 0;     // spatial orbitals, all correlated
    double reference_energy = 0.0;      // RHF total energy
    double correlation_energy = 0.0;    // CC correlation energy without (T)
    double cc_energy = 0.0;             // reference + correlation
    bool has_triples = false;           // CCSD(T)
    double triples_correction = 0.0;    // (T)
    double total_energy = 0.0;          // cc_energy (+ (T))
    double t1_diagnostic = 0.0;         // sqrt(sum t_ia^2 / N), spatial amplitudes
    std::vector<double> amplitude_counts;  // index = rank: number of amplitudes stored
    std::vector<double> rank_norms;     // index = rank: |T_n| over the stored amplitudes
    std::vector<CCAmplitude> leading;   // largest |t| first (one of each spin-equivalent set)
    std::vector<CCIteration> history;
    std::string ladder;                 // closed-shell engine: how the (vv|vv) ladder was evaluated
    double transform_seconds = 0.0, cc_seconds = 0.0, triples_seconds = 0.0;
};

using CCProgress = std::function<void(const CCIteration&)>;

// ---- general-order spin-orbital engine -------------------------------------------------

// Coupled cluster with all excitations up to max_rank (2, 3 or 4) on the MO integrals of
// `mo` (h and (pq|rs) over all orbitals); nocc doubly occupied orbitals. The Fock matrix
// is built from the integrals.
CCResult run_cc_spin_orbital(const MOIntegrals& mo, int nocc, int max_rank, double nuclear_repulsion,
                             const CCOptions& options, const CCProgress& progress = CCProgress());

// Size of the spin-orbital engine for nocc occupied and nvir virtual spatial orbitals.
struct SpinOrbitalCCCost {
    double amplitudes = 0.0;           // stored amplitudes (all ranks)
    double flops_per_iteration = 0.0;  // multiply-adds in the matrix products
    double scatter_per_iteration = 0.0;  // elements gathered and scattered
    double integral_bytes = 0.0;       // spin-orbital Fock and <pq||rs> blocks
    double work_bytes = 0.0;           // intermediates and work buffers
    double fixed_bytes = 0.0;          // integrals + work + amplitudes and residuals
    double vector_bytes = 0.0;         // one amplitude vector (for DIIS)
};
SpinOrbitalCCCost estimate_cc_spin_orbital(int nocc, int nvir, int max_rank);

// ---- closed-shell engine ---------------------------------------------------------------

// Closed-shell MO integrals of the CCSD and (T) equations (chemists' notation, spatial
// orbitals; occupied 0..o-1, virtual 0..v-1 within their blocks).
struct CCIntegrals {
    int o = 0, v = 0;
    Matrix fock;                 // n x n MO Fock matrix
    std::vector<double> oooo;    // (ij|kl)  [i][j][k][l]
    std::vector<double> ovoo;    // (ia|jk)  [i][a][j][k]
    std::vector<double> ovov;    // (ia|jb)  [i][a][j][b]
    std::vector<double> oovv;    // (ij|ab)  [i][j][a][b]
    std::vector<double> ovvv;    // (ia|bc)  [i][a][b][c]
};

// Particle-particle ladder Z_ij^ab = sum_cd (ac|bd) tau_ij^cd for all i, j (o^2 v^2 arrays
// [i][j][a][b]).
class LadderOperator {
public:
    virtual ~LadderOperator() = default;
    virtual void apply(const std::vector<double>& tau, std::vector<double>& z) = 0;
    virtual std::string description() const = 0;
};

// Integrals from full MO integrals (tests and small cases): all blocks and an MO ladder.
CCIntegrals cc_integrals_from_mo(const MOIntegrals& mo, int nocc);
std::unique_ptr<LadderOperator> make_mo_ladder(const MOIntegrals& mo, int nocc);

// Integrals with at least one occupied index, by a transformation that keeps one
// occupied index in the ket (memory o n N(N+1)/2), and the Fock matrix from them.
CCIntegrals cc_integrals_from_ao(const BasisSet& basis, const Matrix& h_ao, const Matrix& C, int nocc,
                                 int threads);
// AO-driven ladder: tau is back-transformed with the virtual MO coefficients, contracted
// with the AO integrals (kept in memory up to eri_memory_mb, otherwise recomputed with
// Schwarz screening) and transformed back.
std::unique_ptr<LadderOperator> make_ao_ladder(const BasisSet& basis, const Matrix& C, int nocc, int threads,
                                               double eri_memory_mb);

// Closed-shell CCSD (and the (T) correction when `triples`).
CCResult run_ccsd_closed_shell(const CCIntegrals& ints, LadderOperator& ladder, double reference_energy,
                               bool triples, const CCOptions& options, const CCProgress& progress = CCProgress());

// Energy of the reference determinant from MO integrals (electronic + nuclear repulsion).
double reference_energy_from_mo(const MOIntegrals& mo, int nocc, double nuclear_repulsion);

}  // namespace macqc
