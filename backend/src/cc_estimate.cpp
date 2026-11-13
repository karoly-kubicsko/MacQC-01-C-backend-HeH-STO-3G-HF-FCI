// Size, memory and time estimates of the coupled-cluster calculations.
//
// Rates were measured with H2O in aug-cc-pVDZ .. aug-cc-pVQZ on four 2.1 GHz x86-64
// cores (the portable matrix product reaches about 2.5e9 multiply-adds per second and
// thread); Apple silicon cores are faster, so the estimates are on the safe side.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "macqc/cc.h"

namespace macqc {

namespace {

constexpr double kMultiplyAddsPerSecond = 2.4e9;  // per thread, matrix products and ladder
constexpr double kMovesPerSecond = 1.0e8;         // per thread, gather/scatter of packed amplitudes
constexpr double kElementwisePerSecond = 1.0e9;   // per thread, (T) energy loops

std::string gigabytes(double bytes) {
    char buf[64];
    if (bytes >= 1024.0 * 1024.0 * 1024.0)
        std::snprintf(buf, sizeof buf, "%.1f GB", bytes / (1024.0 * 1024.0 * 1024.0));
    else
        std::snprintf(buf, sizeof buf, "%.0f MB", bytes / (1024.0 * 1024.0));
    return buf;
}

std::string days(double s) {
    char buf[64];
    if (s >= 2.0 * 86400.0)
        std::snprintf(buf, sizeof buf, "about %.0f days", s / 86400.0);
    else
        std::snprintf(buf, sizeof buf, "about %.0f hours", s / 3600.0);
    return buf;
}

}  // namespace

CCSize estimate_cc(CCMethod method, int nocc, int nmo, int nbasis, double unique_eris, double eri_pass_seconds,
                   int threads_in, double memory_limit_mb, double eri_memory_mb) {
    CCSize c;
    c.method = cc_method_name(method);
    const double o = nocc, v = nmo - nocc, n = nmo, N = nbasis;
    const double t = std::max(1, threads_in);
    c.occupied = nocc;
    c.virtuals = nmo - nocc;
    const int rank = cc_max_rank(method);
    // AO -> MO work buffers of one thread: (ab|..) for the largest shell pair (h x h = 121
    // functions) over all AO pairs, and its half transformation.
    const double transform_buffers = t * 8.0 * 121.0 * (N * N + N * n);
    if (rank == 2) {
        c.engine = "closed-shell, spin-adapted";
        c.amplitudes = o * v + o * o * v * v;
        c.iterations = 18;
        const double half = 8.0 * o * n * N * (N + 1) / 2.0;  // (pq|is) half-transformed
        const double ints = 8.0 * (o * v * v * v + 3.0 * o * o * v * v + o * o * o * v + o * o * o * o);
        const double amps = 8.0 * o * o * v * v * (10.0 + 2.0 * 8.0);  // work arrays + DIIS
        c.ladder_stored = eri_memory_mb > 0.0 && unique_eris * 8.0 / (1024.0 * 1024.0) <= eri_memory_mb;
        const double ladder = (c.ladder_stored ? 8.0 * unique_eris : 0.0) + 8.0 * N * N * o * o * (t + 2.0);
        const double triples = method == CCMethod::CCSD_T ? 8.0 * 3.0 * v * v * v : 0.0;
        c.memory_bytes = std::max(half + ints + transform_buffers, ints + amps + ladder + triples);
        c.transform_seconds = eri_pass_seconds + (N * N * N * n / 2.0 + o * n * (N * N * n + N * n * n)) /
                                                     (kMultiplyAddsPerSecond * t);
        const double ladder_ops = 4.0 * o * o * unique_eris + 4.0 * o * o * (N * N * v + N * v * v);
        const double other_ops = 12.0 * o * o * o * v * v * v + 10.0 * o * o * v * v * v;
        c.iteration_seconds = (ladder_ops + other_ops) / (kMultiplyAddsPerSecond * t) +
                              (c.ladder_stored ? 0.0 : eri_pass_seconds);
        if (method == CCMethod::CCSD_T) {
            const double triplets = o * (o + 1.0) * (o + 2.0) / 6.0;
            c.triples_seconds = triplets * 6.0 * (v * v * v * v + o * v * v * v) / (kMultiplyAddsPerSecond * t) +
                                triplets * 70.0 * v * v * v / (kElementwisePerSecond * t);
        }
    } else {
        c.engine = "spin-orbital, general order";
        const SpinOrbitalCCCost so = estimate_cc_spin_orbital(nocc, nmo - nocc, rank);
        c.amplitudes = so.amplitudes;
        c.iterations = 18 + 2 * (rank - 2);
        const double mo_ints = 8.0 * (n * (n + 1) / 2.0) * (N * (N + 1) / 2.0);  // transform_to_mo
        c.memory_bytes = std::max(mo_ints + transform_buffers, mo_ints + so.fixed_bytes + 2.0 * 3.0 * so.vector_bytes);
        c.transform_seconds =
            eri_pass_seconds + (N * N * N * n / 2.0 + n * n / 2.0 * (N * N * n + N * n * n)) / (kMultiplyAddsPerSecond * t);
        c.iteration_seconds = so.flops_per_iteration / (kMultiplyAddsPerSecond * t) +
                              so.scatter_per_iteration / (kMovesPerSecond * t);
    }
    c.seconds = c.transform_seconds + c.iterations * c.iteration_seconds + c.triples_seconds;
    const double limit = memory_limit_mb * 1024.0 * 1024.0;
    char buf[512];
    if (c.memory_bytes > limit) {
        c.feasible = false;
        std::snprintf(buf, sizeof buf,
                      "%s is not possible here: %.0f occupied and %.0f virtual orbitals give %.3g amplitudes, and the "
                      "calculation would need about %s of memory; ΨLab's limit is %s.",
                      c.method.c_str(), o, v, c.amplitudes, gigabytes(c.memory_bytes).c_str(), gigabytes(limit).c_str());
        c.reason = buf;
    } else if (c.seconds > kMaxCCSeconds) {
        c.feasible = false;
        std::snprintf(buf, sizeof buf,
                      "%s is not possible here: %.0f occupied and %.0f virtual orbitals give %.3g amplitudes, and the "
                      "calculation would take %s on %d threads (ΨLab's limit is three days).",
                      c.method.c_str(), o, v, c.amplitudes, days(c.seconds).c_str(), static_cast<int>(t));
        c.reason = buf;
    }
    return c;
}

}  // namespace macqc
