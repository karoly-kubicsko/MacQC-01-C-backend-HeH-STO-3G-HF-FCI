// Response densities of closed-shell CCSD and CCSD(T) for analytic gradients (ΨLab 5.0).
//
// Lambda equations: the closed-shell CCSD left-hand equations in the form of PySCF's
// ccsd_lambda (Stanton-Gauss intermediates), with the particle-particle ladder through the
// same LadderOperator as the amplitude equations (AO-driven: no (vv|vv) integrals stored).
// CCSD(T): the Lambda equations of the CCSD(T) Lagrangian get the (T) source terms, and
// the densities the (T) contributions, in the formulation where the triples denominators
// come from the full Fock matrix (so the Lagrangian is invariant to occupied-occupied and
// virtual-virtual rotations and only the occupied-virtual orbital response is needed).
// The one- and two-particle densities (spin-traced, all orbitals correlated) are returned in
// the packed form of correlated_gradient (E = sum D h + sum G (pq|rs) + V_nn).
#pragma once

#include <functional>
#include <vector>

#include "macqc/cc.h"
#include "macqc/correlated_gradient.h"

namespace macqc {

struct LambdaIteration {
    int iteration = 0;
    double residual = 0.0;  // largest |lambda_new - lambda|
};

struct CCDensityResult {
    CorrelatedDensity density;
    bool lambda_converged = false;
    int lambda_iterations = 0;
    double lambda_residual = 0.0;
    std::vector<LambdaIteration> history;
    double lambda_seconds = 0.0, density_seconds = 0.0, triples_seconds = 0.0;
};

using LambdaProgress = std::function<void(const LambdaIteration&)>;

// t1, t2 from run_ccsd_closed_shell (converged tightly). `triples` adds the (T) terms.
CCDensityResult ccsd_densities(const CCIntegrals& ints, LadderOperator& ladder, const std::vector<double>& t1,
                               const std::vector<double>& t2, bool triples, const CCOptions& options,
                               const LambdaProgress& progress = LambdaProgress());

// Largest difference between the (T) Lambda sources and density terms on full o^3 v^3 arrays
// (PySCF's formulation) and the blocked implementation used (tests).
double compare_triples_implementations(const CCIntegrals& ints, const std::vector<double>& t1,
                                       const std::vector<double>& t2, int threads = 0);

// Memory (bytes) of ccsd_densities for o occupied and v virtual orbitals (with (T): the
// triples work arrays).
double ccsd_density_bytes(int o, int v, bool triples);

}  // namespace macqc
