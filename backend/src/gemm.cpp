// Portable cache-blocked matrix multiplication (row-major), see macqc/linalg.h.
#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#include "macqc/linalg.h"

// The rest of the backend is compiled without floating-point contraction (identical
// output to the reference program); the matrix product may use fused multiply-adds.
#if defined(__clang__)
#pragma clang fp contract(fast)
#endif

namespace macqc {

namespace {

constexpr std::size_t MR = 4;
#if defined(__aarch64__) || defined(__AVX__)
constexpr std::size_t NR = 8;
#else
constexpr std::size_t NR = 4;
#endif
constexpr std::size_t KC = 256;
constexpr std::size_t MC = 96;
constexpr std::size_t NC = 2048;

struct View {
    const double* p;
    std::size_t ld;
    bool trans;
    double at(std::size_t i, std::size_t j) const { return trans ? p[j * ld + i] : p[i * ld + j]; }
};

// Packs op(A)[i0:i0+mc, p0:p0+kc] into MR-row panels: panel r holds a[p*MR + i].
void pack_a(const View& A, std::size_t i0, std::size_t mc, std::size_t p0, std::size_t kc, double* out) {
    for (std::size_t r = 0; r < mc; r += MR) {
        const std::size_t mr = std::min(MR, mc - r);
        for (std::size_t p = 0; p < kc; ++p)
            for (std::size_t i = 0; i < MR; ++i) *out++ = i < mr ? A.at(i0 + r + i, p0 + p) : 0.0;
    }
}

// Packs op(B)[p0:p0+kc, j0:j0+nc] into NR-column panels: panel c holds b[p*NR + j].
void pack_b(const View& B, std::size_t p0, std::size_t kc, std::size_t j0, std::size_t nc, double* out) {
    for (std::size_t c = 0; c < nc; c += NR) {
        const std::size_t nr = std::min(NR, nc - c);
        if (!B.trans && nr == NR) {
            for (std::size_t p = 0; p < kc; ++p) {
                const double* src = B.p + (p0 + p) * B.ld + j0 + c;
                for (std::size_t j = 0; j < NR; ++j) *out++ = src[j];
            }
        } else {
            for (std::size_t p = 0; p < kc; ++p)
                for (std::size_t j = 0; j < NR; ++j) *out++ = j < nr ? B.at(p0 + p, j0 + c + j) : 0.0;
        }
    }
}

inline void micro_kernel(std::size_t kc, const double* a, const double* b, double alpha, double* c, std::size_t ldc,
                         std::size_t mr, std::size_t nr) {
    double acc[MR][NR] = {};
    for (std::size_t p = 0; p < kc; ++p) {
        const double* bp = b + p * NR;
        const double* ap = a + p * MR;
        for (std::size_t i = 0; i < MR; ++i) {
            const double ai = ap[i];
            for (std::size_t j = 0; j < NR; ++j) acc[i][j] += ai * bp[j];
        }
    }
    for (std::size_t i = 0; i < mr; ++i)
        for (std::size_t j = 0; j < nr; ++j) c[i * ldc + j] += alpha * acc[i][j];
}

void gemm_serial(const View& A, const View& B, std::size_t m, std::size_t n, std::size_t k, double alpha, double* c,
                 std::size_t ldc) {
    std::vector<double> pa(((MC + MR - 1) / MR) * MR * KC), pb(((NC + NR - 1) / NR) * NR * KC);
    for (std::size_t j0 = 0; j0 < n; j0 += NC) {
        const std::size_t nc = std::min(NC, n - j0);
        for (std::size_t p0 = 0; p0 < k; p0 += KC) {
            const std::size_t kc = std::min(KC, k - p0);
            pack_b(B, p0, kc, j0, nc, pb.data());
            for (std::size_t i0 = 0; i0 < m; i0 += MC) {
                const std::size_t mc = std::min(MC, m - i0);
                pack_a(A, i0, mc, p0, kc, pa.data());
                for (std::size_t jr = 0; jr < nc; jr += NR)
                    for (std::size_t ir = 0; ir < mc; ir += MR)
                        micro_kernel(kc, pa.data() + ir * kc, pb.data() + jr * kc, alpha,
                                     c + (i0 + ir) * ldc + j0 + jr, ldc, std::min(MR, mc - ir), std::min(NR, nc - jr));
            }
        }
    }
}

}  // namespace

void gemm(bool trans_a, bool trans_b, std::size_t m, std::size_t n, std::size_t k, double alpha, const double* a,
          std::size_t lda, const double* b, std::size_t ldb, double beta, double* c, std::size_t ldc, int threads) {
    if (m == 0 || n == 0) return;
    if (beta != 1.0)
        for (std::size_t i = 0; i < m; ++i) {
            double* row = c + i * ldc;
            if (beta == 0.0)
                std::fill(row, row + n, 0.0);
            else
                for (std::size_t j = 0; j < n; ++j) row[j] *= beta;
        }
    if (k == 0 || alpha == 0.0) return;
    const View A{a, lda, trans_a}, B{b, ldb, trans_b};
    // Split the rows of C between threads (only when each gets a reasonable share).
    const std::size_t work = m * n * k;
    int t = std::max(1, threads);
    if (work < 4000000) t = 1;
    t = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(t), (m + MR - 1) / MR));
    if (t <= 1) {
        gemm_serial(A, B, m, n, k, alpha, c, ldc);
        return;
    }
    std::vector<std::thread> pool;
    const std::size_t chunk = ((m + t - 1) / t + MR - 1) / MR * MR;
    for (int r = 0; r < t; ++r) {
        const std::size_t i0 = static_cast<std::size_t>(r) * chunk;
        if (i0 >= m) break;
        const std::size_t mi = std::min(chunk, m - i0);
        const View Ai{trans_a ? a + i0 : a + i0 * lda, lda, trans_a};
        pool.emplace_back([=] { gemm_serial(Ai, B, mi, n, k, alpha, c + i0 * ldc, ldc); });
    }
    for (auto& th : pool) th.join();
}

}  // namespace macqc
