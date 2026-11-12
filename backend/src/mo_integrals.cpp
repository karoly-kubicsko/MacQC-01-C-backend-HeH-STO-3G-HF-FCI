#include "macqc/mo_integrals.h"

#include <algorithm>
#include <atomic>
#include <stdexcept>

#include "macqc/gto_integrals.h"
#include "macqc/parallel.h"

namespace macqc {

namespace {

// out_b = C^T in_b C for nb symmetric N x N matrices stored one after another;
// out holds nb n x n matrices.
void transform_batch(const double* in, std::size_t nb, std::size_t N, const double* C, std::size_t n,
                     std::vector<double>& t1, std::vector<double>& t2, double* out) {
    t1.resize(nb * N * n);
    t2.resize(nb * n * N);
    gemm(false, false, nb * N, n, N, 1.0, in, N, C, n, 0.0, t1.data(), n);  // [b][lambda][s]
    for (std::size_t b = 0; b < nb; ++b) {
        const double* src = t1.data() + b * N * n;
        double* dst = t2.data() + b * n * N;
        for (std::size_t l = 0; l < N; ++l)
            for (std::size_t s = 0; s < n; ++s) dst[s * N + l] = src[l * n + s];  // [b][s][lambda]
    }
    gemm(false, false, nb * n, n, N, 1.0, t2.data(), N, C, n, 0.0, out, n);  // [b][s][r]
}

}  // namespace

double mo_transform_bytes(std::size_t n_basis, std::size_t n_mo) {
    return 8.0 * static_cast<double>(n_mo * (n_mo + 1) / 2) * static_cast<double>(n_basis * (n_basis + 1) / 2);
}

MOIntegrals transform_to_mo(const BasisSet& basis, const Matrix& h_ao, const Matrix& Cm, int threads_requested) {
    const std::size_t N = static_cast<std::size_t>(basis.n_functions);
    const std::size_t n = Cm.cols();
    if (Cm.rows() != N || h_ao.rows() != N) throw std::invalid_argument("transform_to_mo: dimension mismatch");
    const int threads = resolve_threads(threads_requested);
    std::vector<double> C(N * n);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < n; ++j) C[i * n + j] = Cm(i, j);

    MOIntegrals mo;
    mo.n = static_cast<int>(n);
    mo.npair = n * (n + 1) / 2;
    const std::size_t npair_ao = N * (N + 1) / 2;
    mo.stride = npair_ao;
    mo.eri.assign(mo.npair * npair_ao, 0.0);

    // h_MO = C^T h C.
    {
        std::vector<double> hin(N * N), t1, t2, out(n * n);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) hin[i * N + j] = h_ao(i, j);
        transform_batch(hin.data(), 1, N, C.data(), n, t1, t2, out.data());
        mo.h = Matrix(n, n);
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t s = 0; s < n; ++s) mo.h(r, s) = out[s * n + r];
    }

    // Step A: X[rs][mu nu] = sum_{lambda sigma} C_lambda r C_sigma s (mu nu|lambda sigma).
    const ShellPairList pairs = make_shell_pairs(basis);
    const std::size_t npairs = pairs.pairs.size();
    {
        std::atomic<std::size_t> next{0};
        run_parallel(threads, [&](int) {
            EriEngine engine(basis);
            std::vector<double> G, t1, t2, out;
            for (std::size_t k; (k = next++) < npairs;) {
                const ShellPair& bra = pairs.pairs[npairs - 1 - k];
                if (bra.prims.empty()) continue;
                const Shell& A = basis.shells[bra.a];
                const Shell& B = basis.shells[bra.b];
                const std::size_t na = A.size(), nb = B.size(), nab = na * nb;
                G.assign(nab * N * N, 0.0);
                for (const ShellPair& ket : pairs.pairs) {
                    if (ket.prims.empty()) continue;
                    const Shell& Cs = basis.shells[ket.a];
                    const Shell& Ds = basis.shells[ket.b];
                    const double* v = engine.compute(bra, ket);
                    for (std::size_t ab = 0; ab < nab; ++ab) {
                        double* g = G.data() + ab * N * N;
                        for (int c = 0; c < Cs.size(); ++c) {
                            const std::size_t l = static_cast<std::size_t>(Cs.first_function + c);
                            for (int d = 0; d < Ds.size(); ++d) {
                                const std::size_t s = static_cast<std::size_t>(Ds.first_function + d);
                                const double x = *v++;
                                g[l * N + s] = x;
                                g[s * N + l] = x;
                            }
                        }
                    }
                }
                out.resize(nab * n * n);
                transform_batch(G.data(), nab, N, C.data(), n, t1, t2, out.data());
                for (std::size_t a = 0; a < na; ++a)
                    for (std::size_t b = 0; b < nb; ++b) {
                        const std::size_t mu = static_cast<std::size_t>(A.first_function) + a;
                        const std::size_t nu = static_cast<std::size_t>(B.first_function) + b;
                        if (bra.a == bra.b && mu < nu) continue;
                        const std::size_t col = MOIntegrals::pair(mu, nu);
                        const double* o = out.data() + (a * nb + b) * n * n;
                        for (std::size_t r = 0; r < n; ++r)
                            for (std::size_t s = 0; s <= r; ++s)
                                mo.eri[(r * (r + 1) / 2 + s) * npair_ao + col] = o[s * n + r];
                    }
            }
        });
    }

    // Step B: row rs -> (pq|rs) = sum_{mu nu} C_mu p C_nu q X[rs][mu nu], in place.
    {
        const std::size_t batch = 16;
        std::atomic<std::size_t> next{0};
        run_parallel(threads, [&](int) {
            std::vector<double> Y, t1, t2, out;
            for (std::size_t b0; (b0 = batch * next++) < mo.npair;) {
                const std::size_t nb = std::min(batch, mo.npair - b0);
                Y.assign(nb * N * N, 0.0);
                for (std::size_t b = 0; b < nb; ++b) {
                    const double* row = mo.eri.data() + (b0 + b) * npair_ao;
                    double* y = Y.data() + b * N * N;
                    for (std::size_t mu = 0; mu < N; ++mu)
                        for (std::size_t nu = 0; nu <= mu; ++nu) {
                            const double x = row[mu * (mu + 1) / 2 + nu];
                            y[mu * N + nu] = x;
                            y[nu * N + mu] = x;
                        }
                }
                out.resize(nb * n * n);
                transform_batch(Y.data(), nb, N, C.data(), n, t1, t2, out.data());
                for (std::size_t b = 0; b < nb; ++b) {
                    double* row = mo.eri.data() + (b0 + b) * npair_ao;
                    const double* o = out.data() + b * n * n;
                    for (std::size_t p = 0; p < n; ++p)
                        for (std::size_t q = 0; q <= p; ++q) row[p * (p + 1) / 2 + q] = o[q * n + p];
                }
            }
        });
    }
    return mo;
}

}  // namespace macqc
