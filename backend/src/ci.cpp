#include "macqc/ci.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

#include "macqc/gto_integrals.h"
#include "macqc/linalg.h"
#include "macqc/parallel.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

std::string human_count(double v) {
    char buf[64];
    if (v >= 1e9) {
        std::snprintf(buf, sizeof buf, "%.2e", v);
        return buf;
    }
    std::snprintf(buf, sizeof buf, "%.0f", v);
    std::string s = buf, out;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (i && (s.size() - i) % 3 == 0) out += ',';
        out += s[i];
    }
    return out;
}

std::string human_bytes(double b) {
    char buf[64];
    if (b < 1024.0 * 1024.0 * 1024.0) std::snprintf(buf, sizeof buf, "%.0f MB", b / (1024.0 * 1024.0));
    else if (b < 1024.0 * 1024.0 * 1024.0 * 1024.0) std::snprintf(buf, sizeof buf, "%.1f GB", b / (1024.0 * 1024.0 * 1024.0));
    else std::snprintf(buf, sizeof buf, "%.1f TB", b / (1024.0 * 1024.0 * 1024.0 * 1024.0));
    return buf;
}

// E_pq |J> = sign |target> (p created, q annihilated; p == q included).
struct Excitation {
    std::uint32_t target;
    std::uint16_t p, q;
    double sign;
};

// a+_orb |K> = sign |target>  (or, read backwards, a_orb |target> = sign |K>).
struct Hop {
    std::uint32_t target;
    std::uint16_t orb;
    double sign;
};

template <class T>
struct Lists {
    std::vector<std::size_t> offset;  // size + 1
    std::vector<T> items;
    const T* begin(std::size_t i) const { return items.data() + offset[i]; }
    const T* end(std::size_t i) const { return items.data() + offset[i + 1]; }
};

// Inserts p into the sorted list o[0..m) -> out, returns the insertion position.
int insert_sorted(const std::uint16_t* o, int m, int p, std::uint16_t* out) {
    int pos = 0;
    while (pos < m && o[pos] < p) ++pos;
    for (int k = 0; k < pos; ++k) out[k] = o[k];
    out[pos] = static_cast<std::uint16_t>(p);
    for (int k = pos; k < m; ++k) out[k + 1] = o[k];
    return pos;
}

Lists<Excitation> single_excitations(const StringSpace& S) {
    const int n = S.orbitals(), N = S.electrons();
    Lists<Excitation> L;
    L.offset.assign(S.size() + 1, 0);
    std::vector<std::uint16_t> R(static_cast<std::size_t>(std::max(N, 1))), I(static_cast<std::size_t>(N) + 1);
    for (std::size_t J = 0; J < S.size(); ++J) {
        const std::uint16_t* o = S.occupation(J);
        for (int t = 0; t < N; ++t) {
            const int q = o[t];
            int m = 0;
            for (int k = 0; k < N; ++k)
                if (k != t) R[static_cast<std::size_t>(m++)] = o[k];
            const double s1 = (t % 2) ? -1.0 : 1.0;
            for (int p = 0; p < n; ++p) {
                if (std::binary_search(R.begin(), R.begin() + m, static_cast<std::uint16_t>(p))) continue;
                const int pos = insert_sorted(R.data(), m, p, I.data());
                const double s2 = (pos % 2) ? -1.0 : 1.0;
                L.items.push_back({static_cast<std::uint32_t>(S.index(I.data())), static_cast<std::uint16_t>(p),
                                   static_cast<std::uint16_t>(q), s1 * s2});
            }
        }
        L.offset[J + 1] = L.items.size();
    }
    return L;
}

// For every (N-1)-string K and orbital j not in K: a+_j |K> = sign |J>, J in `full`.
Lists<Hop> creations(const StringSpace& minus, const StringSpace& full) {
    const int n = minus.orbitals(), m = minus.electrons();
    Lists<Hop> L;
    L.offset.assign(minus.size() + 1, 0);
    std::vector<std::uint16_t> J(static_cast<std::size_t>(m) + 1);
    for (std::size_t K = 0; K < minus.size(); ++K) {
        const std::uint16_t* o = minus.occupation(K);
        for (int j = 0; j < n; ++j) {
            if (std::binary_search(o, o + m, static_cast<std::uint16_t>(j))) continue;
            const int pos = insert_sorted(o, m, j, J.data());
            L.items.push_back({static_cast<std::uint32_t>(full.index(J.data())), static_cast<std::uint16_t>(j),
                               (pos % 2) ? -1.0 : 1.0});
        }
        L.offset[K + 1] = L.items.size();
    }
    return L;
}

// For every N-string I and occupied orbital k: a_k |I> = sign |K>, K in `minus`.
Lists<Hop> annihilations(const StringSpace& full, const StringSpace& minus) {
    const int N = full.electrons();
    Lists<Hop> L;
    L.offset.assign(full.size() + 1, 0);
    std::vector<std::uint16_t> K(static_cast<std::size_t>(std::max(N, 1)));
    for (std::size_t I = 0; I < full.size(); ++I) {
        const std::uint16_t* o = full.occupation(I);
        for (int t = 0; t < N; ++t) {
            int m = 0;
            for (int k = 0; k < N; ++k)
                if (k != t) K[static_cast<std::size_t>(m++)] = o[k];
            L.items.push_back({static_cast<std::uint32_t>(minus.index(K.data())), o[t], (t % 2) ? -1.0 : 1.0});
        }
        L.offset[I + 1] = L.items.size();
    }
    return L;
}

// One-spin string Hamiltonian F(I,J) = <I| sum k_pq E_pq + 1/2 sum (pq|rs) E_pq E_rs |J>,
// k_pq = h_pq - 1/2 sum_r (pr|rq).
std::vector<double> string_hamiltonian(const StringSpace& S, const Lists<Excitation>& E, const MOIntegrals& mo,
                                       int threads) {
    const std::size_t Ns = S.size();
    const int n = mo.n;
    Matrix k(static_cast<std::size_t>(n), static_cast<std::size_t>(n));
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            double v = mo.h(static_cast<std::size_t>(p), static_cast<std::size_t>(q));
            for (int r = 0; r < n; ++r) v -= 0.5 * mo(p, r, r, q);
            k(static_cast<std::size_t>(p), static_cast<std::size_t>(q)) = v;
        }
    std::vector<double> F(Ns * Ns, 0.0);
    std::atomic<std::size_t> next{0};
    run_parallel(threads, [&](int) {
        for (std::size_t J; (J = next++) < Ns;) {
            for (const Excitation* e1 = E.begin(J); e1 != E.end(J); ++e1) {
                const std::size_t K = e1->target;
                F[K * Ns + J] += k(e1->p, e1->q) * e1->sign;
                const std::size_t rs = MOIntegrals::pair(e1->p, e1->q);
                for (const Excitation* e2 = E.begin(K); e2 != E.end(K); ++e2)
                    F[static_cast<std::size_t>(e2->target) * Ns + J] +=
                        0.5 * mo.eri[MOIntegrals::pair(e2->p, e2->q) * mo.stride + rs] * e1->sign * e2->sign;
            }
        }
    });
    return F;
}

// Davidson subspace size: as many vector pairs (b, sigma) as fit into the memory left after
// the fixed arrays, between 3 and 10.
int davidson_subspace(double determinants, double fixed_bytes, double memory_limit_bytes) {
    const double pair = 2.0 * 8.0 * determinants;
    const double fit = std::floor((memory_limit_bytes - fixed_bytes) / pair);
    return static_cast<int>(std::clamp(fit, 3.0, 10.0));
}

// Arrays besides the Davidson subspace: MO integrals, string Hamiltonian(s), (ij|kl) matrix,
// D and E batches, and five work vectors (diagonal, x, sigma x, residual, correction).
double fixed_fci_bytes(double n, double N, double na_strings, double nb_strings, bool same, double ka, double kb,
                       bool dense) {
    const double X = mo_transform_bytes(static_cast<std::size_t>(N), static_cast<std::size_t>(n));
    const double F = 8.0 * (na_strings * na_strings + (same ? 0.0 : nb_strings * nb_strings));
    const double M = dense ? 8.0 * n * n * n * n : 0.0;
    const double DE = std::min(2.0 * 128.0 * 1048576.0, 2.0 * 8.0 * n * n * ka * kb);
    return X + F + M + DE + 5.0 * 8.0 * na_strings * nb_strings;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y) {
    for (std::size_t i = 0; i < x.size(); ++i) y[i] += alpha * x[i];
}

// The Hamiltonian in the determinant basis (applied, never stored).
class FCIHamiltonian {
public:
    FCIHamiltonian(const MOIntegrals& mo, int na, int nb, int threads, double memory_limit_bytes)
        : mo_(mo), n_(mo.n), threads_(threads), A_(mo.n, na), B_(mo.n, nb),
          Ka_(mo.n, na - 1), Kb_(mo.n, nb - 1) {
        same_ = (na == nb);
        exA_ = single_excitations(A_);
        if (!same_) exB_ = single_excitations(B_);
        Fa_ = string_hamiltonian(A_, exA_, mo, threads);
        if (!same_) Fb_ = string_hamiltonian(B_, exB_, mo, threads);
        creA_ = creations(Ka_, A_);
        creB_ = creations(Kb_, B_);
        annB_ = annihilations(B_, Kb_);
        // Diagonal.
        const std::size_t Na = A_.size(), Nb = B_.size();
        const std::vector<double>& Fb = same_ ? Fa_ : Fb_;
        diag_.resize(Na * Nb);
        for (std::size_t I = 0; I < Na; ++I) {
            const std::uint16_t* oa = A_.occupation(I);
            for (std::size_t J = 0; J < Nb; ++J) {
                const std::uint16_t* ob = B_.occupation(J);
                double v = Fa_[I * Na + I] + Fb[J * Nb + J];
                for (int i = 0; i < na; ++i)
                    for (int j = 0; j < nb; ++j) v += mo(oa[i], oa[i], ob[j], ob[j]);
                diag_[I * Nb + J] = v;
            }
        }
        // (ij|kl) as the matrix M[(i,k)][(j,l)] when the opposite-spin products are large.
        const double nn = static_cast<double>(n_) * n_;
        dense_ = Ka_.size() * Kb_.size() > 1 && nn * nn * 8.0 <= std::min(1024.0 * 1024.0 * 1024.0, 0.5 * memory_limit_bytes);
        if (dense_) {
            const std::size_t N2 = static_cast<std::size_t>(n_) * n_;
            M_.assign(N2 * N2, 0.0);
            std::atomic<int> next{0};
            run_parallel(threads, [&](int) {
                for (int i; (i = next++) < n_;)
                    for (int k = 0; k < n_; ++k) {
                        double* row = M_.data() + (static_cast<std::size_t>(i) * n_ + k) * N2;
                        for (int j = 0; j < n_; ++j) {
                            const double* src = mo.eri.data() + MOIntegrals::pair(i, j) * mo.stride;
                            for (int l = 0; l < n_; ++l) row[j * n_ + l] = src[MOIntegrals::pair(k, l)];
                        }
                    }
            });
        }
        // Batch of K_alpha strings so that D and E stay below ~128 MB each.
        const double per_ka = nn * static_cast<double>(Kb_.size()) * 8.0;
        batch_ = static_cast<std::size_t>(std::max(1.0, std::floor(128.0 * 1024 * 1024 / per_ka)));
        batch_ = std::min(batch_, Ka_.size());
    }

    std::size_t size() const { return A_.size() * B_.size(); }
    std::size_t alpha_strings() const { return A_.size(); }
    std::size_t beta_strings() const { return B_.size(); }
    const std::vector<double>& diagonal() const { return diag_; }
    const StringSpace& alpha() const { return A_; }
    const StringSpace& beta() const { return B_; }
    const Lists<Excitation>& alpha_excitations() const { return exA_; }
    bool dense() const { return dense_; }
    bool same_spin_strings() const { return same_; }

    // S = H C (electronic). `symmetric` may be set when C = C^T (n_alpha == n_beta).
    void apply(const std::vector<double>& C, std::vector<double>& S, bool symmetric) const {
        const std::size_t Na = A_.size(), Nb = B_.size();
        S.assign(Na * Nb, 0.0);
        // Same-spin parts: sigma_1 = C F_beta, sigma_2 = F_alpha C (= sigma_1^T for symmetric C).
        const std::vector<double>& Fb = same_ ? Fa_ : Fb_;
        gemm(false, false, Na, Nb, Nb, 1.0, C.data(), Nb, Fb.data(), Nb, 0.0, S.data(), Nb, threads_);
        if (symmetric && same_) {
            // S = T + T^T with T = C F + (half of the opposite-spin part).
            opposite_spin_half(C, S);
            for (std::size_t i = 0; i < Na; ++i) {
                S[i * Nb + i] *= 2.0;
                for (std::size_t j = 0; j < i; ++j) {
                    const double v = S[i * Nb + j] + S[j * Nb + i];
                    S[i * Nb + j] = v;
                    S[j * Nb + i] = v;
                }
            }
        } else {
            gemm(false, false, Na, Nb, Na, 1.0, Fa_.data(), Na, C.data(), Nb, 1.0, S.data(), Nb, threads_);
            opposite_spin(C, S);
        }
    }

private:
    // E = (ij|kl) D for the column pairs (K_a, K_b) with K_a <= K_b only. For C = C^T,
    // E_(k,i)(K_b,K_a) = E_(i,k)(K_a,K_b), so the opposite-spin sigma is T + T^T with T the
    // scatter of these columns (the K_a = K_b columns with weight 1/2). Adds T to S.
    void opposite_spin_half(const std::vector<double>& C, std::vector<double>& S) const {
        const std::size_t n = static_cast<std::size_t>(n_), nn = n * n;
        const std::size_t Nb = B_.size(), Nk = Ka_.size();
        const std::size_t max_cols = std::max<std::size_t>(Nk, static_cast<std::size_t>(128.0 * 1048576.0 / (8.0 * nn)));
        std::vector<double> D, E;
        std::vector<std::size_t> coff;
        for (std::size_t ka0 = 0; ka0 < Nk;) {
            // Batch of K_a strings whose columns (K_b >= K_a) fit into max_cols.
            std::size_t ka1 = ka0, cols = 0;
            coff.clear();
            while (ka1 < Nk && (ka1 == ka0 || cols + (Nk - ka1) <= max_cols)) {
                coff.push_back(cols);
                cols += Nk - ka1;
                ++ka1;
            }
            const std::size_t nk = ka1 - ka0;
            D.assign(nn * cols, 0.0);
            {
                std::atomic<std::size_t> next{0};
                run_parallel(std::min<int>(threads_, static_cast<int>(nk)), [&](int) {
                    for (std::size_t kl; (kl = next++) < nk;) {
                        const std::size_t ka = ka0 + kl;
                        for (const Hop* a = creA_.begin(ka); a != creA_.end(ka); ++a) {
                            const double* crow = C.data() + static_cast<std::size_t>(a->target) * Nb;
                            for (std::size_t kb = ka; kb < Nk; ++kb)
                                for (const Hop* b = creB_.begin(kb); b != creB_.end(kb); ++b)
                                    D[(a->orb * n + b->orb) * cols + coff[kl] + (kb - ka)] =
                                        a->sign * b->sign * crow[b->target];
                        }
                    }
                });
            }
            E.assign(nn * cols, 0.0);
            contract(D, E, cols);
            std::atomic<std::size_t> next{0};
            run_parallel(threads_, [&](int) {
                for (std::size_t Ib; (Ib = next++) < Nb;)
                    for (const Hop* b = annB_.begin(Ib); b != annB_.end(Ib); ++b) {
                        const std::size_t kb = b->target;
                        const std::size_t last = std::min<std::size_t>(ka1, kb + 1);  // K_a <= K_b
                        for (std::size_t ka = ka0; ka < last; ++ka) {
                            const double w = (ka == kb ? 0.5 : 1.0) * b->sign;
                            const std::size_t col = coff[ka - ka0] + (kb - ka);
                            for (const Hop* a = creA_.begin(ka); a != creA_.end(ka); ++a)
                                S[static_cast<std::size_t>(a->target) * Nb + Ib] +=
                                    w * a->sign * E[(a->orb * n + b->orb) * cols + col];
                        }
                    }
            });
            ka0 = ka1;
        }
    }

    // E = M D: one matrix product, or a direct contraction over the packed MO integrals.
    void contract(const std::vector<double>& D, std::vector<double>& E, std::size_t cols) const {
        const std::size_t n = static_cast<std::size_t>(n_), nn = n * n;
        if (dense_) {
            gemm(false, false, nn, cols, nn, 1.0, M_.data(), nn, D.data(), cols, 0.0, E.data(), cols, threads_);
            return;
        }
        std::atomic<std::size_t> next{0};
        run_parallel(threads_, [&](int) {
            for (std::size_t i; (i = next++) < n;)
                for (std::size_t j = 0; j < n; ++j) {
                    const double* row = mo_.eri.data() + MOIntegrals::pair(i, j) * mo_.stride;
                    for (std::size_t k = 0; k < n; ++k) {
                        double* e = E.data() + (i * n + k) * cols;
                        for (std::size_t l = 0; l < n; ++l) {
                            const double w = row[MOIntegrals::pair(k, l)];
                            const double* d = D.data() + (j * n + l) * cols;
                            for (std::size_t c = 0; c < cols; ++c) e[c] += w * d[c];
                        }
                    }
                }
        });
    }

    void opposite_spin(const std::vector<double>& C, std::vector<double>& S) const {
        const std::size_t n = static_cast<std::size_t>(n_), nn = n * n;
        const std::size_t Nb = B_.size(), Nka = Ka_.size(), Nkb = Kb_.size();
        std::vector<double> D, E;
        for (std::size_t ka0 = 0; ka0 < Nka; ka0 += batch_) {
            const std::size_t nk = std::min(batch_, Nka - ka0), cols = nk * Nkb;
            D.assign(nn * cols, 0.0);
            {
                std::atomic<std::size_t> next{0};
                run_parallel(std::min<int>(threads_, static_cast<int>(nk)), [&](int) {
                    for (std::size_t kl; (kl = next++) < nk;)
                        for (const Hop* a = creA_.begin(ka0 + kl); a != creA_.end(ka0 + kl); ++a) {
                            const double* crow = C.data() + static_cast<std::size_t>(a->target) * Nb;
                            for (std::size_t kb = 0; kb < Nkb; ++kb)
                                for (const Hop* b = creB_.begin(kb); b != creB_.end(kb); ++b)
                                    D[(a->orb * n + b->orb) * cols + kl * Nkb + kb] = a->sign * b->sign * crow[b->target];
                        }
                });
            }
            E.assign(nn * cols, 0.0);
            contract(D, E, cols);
            std::atomic<std::size_t> next{0};
            run_parallel(threads_, [&](int) {
                for (std::size_t Ib; (Ib = next++) < Nb;)
                    for (const Hop* b = annB_.begin(Ib); b != annB_.end(Ib); ++b)
                        for (std::size_t kl = 0; kl < nk; ++kl)
                            for (const Hop* a = creA_.begin(ka0 + kl); a != creA_.end(ka0 + kl); ++a)
                                S[static_cast<std::size_t>(a->target) * Nb + Ib] +=
                                    a->sign * b->sign * E[(a->orb * n + b->orb) * cols + kl * Nkb + b->target];
            });
        }
    }

    const MOIntegrals& mo_;
    int n_, threads_;
    bool same_ = true, dense_ = false;
    StringSpace A_, B_, Ka_, Kb_;
    Lists<Excitation> exA_, exB_;
    Lists<Hop> creA_, creB_, annB_;
    std::vector<double> Fa_, Fb_, diag_, M_;
    std::size_t batch_ = 1;
};

void symmetrize(std::vector<double>& v, std::size_t N) {
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const double a = 0.5 * (v[i * N + j] + v[j * N + i]);
            v[i * N + j] = a;
            v[j * N + i] = a;
        }
}

// <S^2> = |S_+ Psi|^2 for M_S = 0.
double spin_squared(const std::vector<double>& C, const StringSpace& A, const StringSpace& B) {
    const int n = A.orbitals(), na = A.electrons(), nb = B.electrons();
    if (na != nb) return std::numeric_limits<double>::quiet_NaN();
    if (nb == 0 || na == n) return 0.0;
    const StringSpace Ap(n, na + 1), Bm(n, nb - 1);
    const std::size_t Nb = B.size(), Nbm = Bm.size();
    std::vector<double> T(Ap.size() * Nbm, 0.0);
    std::vector<std::uint16_t> buf(static_cast<std::size_t>(na) + 1);
    for (int p = 0; p < n; ++p) {
        std::vector<std::pair<std::size_t, std::pair<std::size_t, double>>> la, lb;  // source, (target, sign)
        for (std::size_t J = 0; J < A.size(); ++J) {
            const std::uint16_t* o = A.occupation(J);
            if (std::binary_search(o, o + na, static_cast<std::uint16_t>(p))) continue;
            const int pos = insert_sorted(o, na, p, buf.data());
            la.push_back({J, {Ap.index(buf.data()), (pos % 2) ? -1.0 : 1.0}});
        }
        for (std::size_t J = 0; J < Nb; ++J) {
            const std::uint16_t* o = B.occupation(J);
            const auto it = std::lower_bound(o, o + nb, static_cast<std::uint16_t>(p));
            if (it == o + nb || *it != p) continue;
            const int t = static_cast<int>(it - o);
            int m = 0;
            for (int k = 0; k < nb; ++k)
                if (k != t) buf[static_cast<std::size_t>(m++)] = o[k];
            lb.push_back({J, {Bm.index(buf.data()), (t % 2) ? -1.0 : 1.0}});
        }
        for (const auto& a : la)
            for (const auto& b : lb)
                T[a.second.first * Nbm + b.second.first] += a.second.second * b.second.second * C[a.first * Nb + b.first];
    }
    double s = 0.0;
    for (double x : T) s += x * x;
    return s;
}

}  // namespace

// ------------------------------------------------------------------ StringSpace

StringSpace::StringSpace(int n, int N) : n_(n), N_(N) {
    if (n < 1 || n > 65535 || N < 0 || N > n) throw std::invalid_argument("StringSpace: invalid size");
    if (binomial_count(n, N) > 4.0e9) throw std::runtime_error("StringSpace: too many strings");
    binom_.assign(static_cast<std::size_t>(n) + 1, std::vector<std::size_t>(static_cast<std::size_t>(N) + 2, 0));
    for (int i = 0; i <= n; ++i) {
        binom_[static_cast<std::size_t>(i)][0] = 1;
        for (int k = 1; k <= N + 1; ++k)
            binom_[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                i == 0 ? 0
                       : binom_[static_cast<std::size_t>(i) - 1][static_cast<std::size_t>(k) - 1] +
                             binom_[static_cast<std::size_t>(i) - 1][static_cast<std::size_t>(k)];
    }
    count_ = binom_[static_cast<std::size_t>(n)][static_cast<std::size_t>(N)];
    occ_.resize(count_ * static_cast<std::size_t>(N));
    for (std::size_t idx = 0; idx < count_; ++idx) {
        std::size_t rem = idx;
        int hi = n - 1;
        for (int k = N - 1; k >= 0; --k) {
            while (binom_[static_cast<std::size_t>(hi)][static_cast<std::size_t>(k) + 1] > rem) --hi;
            occ_[idx * static_cast<std::size_t>(N) + static_cast<std::size_t>(k)] = static_cast<std::uint16_t>(hi);
            rem -= binom_[static_cast<std::size_t>(hi)][static_cast<std::size_t>(k) + 1];
            --hi;
        }
    }
}

std::size_t StringSpace::index(const std::uint16_t* o) const {
    std::size_t idx = 0;
    for (int k = 0; k < N_; ++k) idx += binom_[o[k]][static_cast<std::size_t>(k) + 1];
    return idx;
}

double binomial_count(int n, int N) {
    if (N < 0 || N > n) return 0.0;
    double r = 1.0;
    for (int i = 1; i <= N; ++i) r = r * (n - N + i) / i;
    return std::round(r);
}

// ------------------------------------------------------------------ estimate

FCISize estimate_fci(int n, int N, int na, int nb, int threads, double memory_limit_mb, double eri_pass_seconds) {
    FCISize s;
    s.orbitals = n;
    s.alpha = na;
    s.beta = nb;
    s.alpha_strings = binomial_count(n, na);
    s.beta_strings = binomial_count(n, nb);
    s.determinants = s.alpha_strings * s.beta_strings;
    const double t = std::max(1, threads);
    const double nd = n, Nd = N, nn = nd * nd;
    const double ka = binomial_count(n, na - 1), kb = binomial_count(n, nb - 1);
    // Memory: the transformation (half-transformed integrals and per-thread buffers), then
    // the CI arrays and a Davidson subspace of at least three vector pairs.
    const bool dense = ka * kb > 1.0 && nn * nn * 8.0 <= std::min(1024.0 * 1024.0 * 1024.0, 0.5 * memory_limit_mb * 1048576.0);
    const double transform = mo_transform_bytes(static_cast<std::size_t>(N), static_cast<std::size_t>(n)) +
                             t * 4.0 * 81.0 * Nd * Nd * 8.0;
    const double fixed = fixed_fci_bytes(nd, Nd, s.alpha_strings, s.beta_strings, na == nb, ka, kb, dense);
    const int sub = davidson_subspace(s.determinants, fixed, memory_limit_mb * 1048576.0);
    s.memory_bytes = std::max(transform, fixed + 2.0 * sub * 8.0 * s.determinants);
    // Time: matrix products at ~6 GFlop/s per thread, the integral-direct contraction at ~1 GFlop/s.
    const double gemm_rate = 6.0e9 * t, loop_rate = 1.0e9 * t;
    const double transform_flops = (0.5 * Nd * Nd + 0.5 * nd * nd) * 2.0 * (Nd * Nd * nd + Nd * nd * nd);
    // Opposite-spin product: only the columns K_a <= K_b for a symmetric vector (n_alpha == n_beta).
    const double columns = na == nb ? 0.5 * ka * (ka + 1.0) : ka * kb;
    const double sigma_seconds = 2.0 * s.alpha_strings * s.beta_strings * s.beta_strings / gemm_rate +
                                 2.0 * nn * nn * columns / (dense ? gemm_rate : loop_rate);
    s.seconds = 2.0 * eri_pass_seconds + transform_flops / gemm_rate + 14.0 * sigma_seconds;

    if (s.determinants > kMaxFCIDeterminants) {
        s.feasible = false;
        char buf[512];
        std::snprintf(buf, sizeof buf,
                      "Full CI is not possible here: %d electrons (%d alpha, %d beta) in %d orbitals give %s "
                      "determinants (%s alpha x %s beta strings); one CI vector alone would need %s. "
                      "\xCE\xA8Lab's limit is %s determinants.",
                      na + nb, na, nb, n, human_count(s.determinants).c_str(), human_count(s.alpha_strings).c_str(),
                      human_count(s.beta_strings).c_str(), human_bytes(8.0 * s.determinants).c_str(),
                      human_count(kMaxFCIDeterminants).c_str());
        s.reason = buf;
    } else if (s.memory_bytes > memory_limit_mb * 1048576.0) {
        s.feasible = false;
        char buf[512];
        std::snprintf(buf, sizeof buf,
                      "Full CI would need about %s of memory (%s determinants, %d orbitals), more than the "
                      "limit of %s.",
                      human_bytes(s.memory_bytes).c_str(), human_count(s.determinants).c_str(), n,
                      human_bytes(memory_limit_mb * 1048576.0).c_str());
        s.reason = buf;
    }
    return s;
}

// ------------------------------------------------------------------ Full CI

FCIResultGeneral run_fci_determinants(const MOIntegrals& mo, int na, int nb, double nuclear_repulsion, double rhf_energy,
                                      const FCIOptions& opt, const FCIProgress& progress) {
    const auto t0 = Clock::now();
    const int threads = resolve_threads(opt.threads);
    if (na < 1 || nb < 1 || na > mo.n || nb > mo.n) throw std::runtime_error("Full CI: invalid number of electrons");
    if (binomial_count(mo.n, na) * binomial_count(mo.n, nb) > kMaxFCIDeterminants)
        throw std::runtime_error("Full CI: determinant space above the limit");

    const FCIHamiltonian H(mo, na, nb, threads, opt.memory_limit_mb * 1048576.0);
    const std::size_t Na = H.alpha_strings(), Nb = H.beta_strings(), Nd = H.size();
    const bool symmetric = H.same_spin_strings();
    const std::vector<double>& hd = H.diagonal();

    FCIResultGeneral r;
    r.orbitals = mo.n;
    r.alpha = na;
    r.beta = nb;
    r.alpha_strings = Na;
    r.beta_strings = Nb;
    r.determinants = Nd;
    r.dense_integrals = H.dense();

    // Davidson. Subspace size from the memory that is left for vectors.
    const double ka = binomial_count(mo.n, na - 1), kb = binomial_count(mo.n, nb - 1);
    const int max_sub = davidson_subspace(
        static_cast<double>(Nd),
        fixed_fci_bytes(mo.n, static_cast<double>(mo.stride > 0 ? 0.5 * (std::sqrt(8.0 * mo.stride + 1.0) - 1.0) : mo.n),
                        static_cast<double>(Na), static_cast<double>(Nb), symmetric, ka, kb, H.dense()),
        opt.memory_limit_mb * 1048576.0);
    std::vector<std::vector<double>> bv, sv;
    std::vector<double> x(Nd, 0.0), sx(Nd, 0.0), res(Nd, 0.0);
    {
        std::size_t start = 0;  // lowest diagonal element (the RHF determinant)
        for (std::size_t i = 1; i < Nd; ++i)
            if (hd[i] < hd[start]) start = i;
        std::vector<double> b0(Nd, 0.0);
        b0[start] = 1.0;
        if (symmetric) {
            const std::size_t ia = start / Nb, ib = start % Nb;
            if (ia != ib) {
                b0[start] = std::sqrt(0.5);
                b0[ib * Nb + ia] = std::sqrt(0.5);
            }
        }
        bv.push_back(b0);
        sv.emplace_back();
        H.apply(bv.back(), sv.back(), symmetric);
    }
    std::vector<std::vector<double>> G(1, std::vector<double>(1, dot(bv[0], sv[0])));
    double theta = 0.0, theta_prev = 0.0;
    std::vector<double> y_prev;  // Ritz vector of the previous iteration in the current basis
    for (int it = 1; it <= opt.max_iter; ++it) {
        const std::size_t m = bv.size();
        Matrix Gm(m, m);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) Gm(i, j) = 0.5 * (G[i][j] + G[j][i]);
        const EigenResult eg = symmetric_eigen(Gm);
        theta = eg.values[0];
        std::fill(x.begin(), x.end(), 0.0);
        std::fill(sx.begin(), sx.end(), 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            axpy(eg.vectors(i, 0), bv[i], x);
            axpy(eg.vectors(i, 0), sv[i], sx);
        }
        double rn = 0.0;
        for (std::size_t i = 0; i < Nd; ++i) {
            res[i] = sx[i] - theta * x[i];
            rn += res[i] * res[i];
        }
        rn = std::sqrt(rn);
        FCIIteration rec{it, theta + nuclear_repulsion, rn, static_cast<int>(m)};
        r.history.push_back(rec);
        r.iterations = it;
        if (progress) progress(rec);
        if (rn < opt.residual_tol && (it > 1 && std::fabs(theta - theta_prev) < opt.energy_tol)) {
            r.converged = true;
            break;
        }
        if (rn < 1e-3 * opt.residual_tol) {  // exact within the subspace
            r.converged = true;
            break;
        }
        theta_prev = theta;
        std::vector<double> y_now(m);
        for (std::size_t i = 0; i < m; ++i) y_now[i] = eg.vectors(i, 0);
        if (static_cast<int>(m) >= max_sub) {
            // Thick restart: keep the current and the previous Ritz vector (orthonormalised).
            std::vector<double> x2(Nd, 0.0), s2(Nd, 0.0);
            for (std::size_t i = 0; i < y_prev.size(); ++i) {
                axpy(y_prev[i], bv[i], x2);
                axpy(y_prev[i], sv[i], s2);
            }
            const double c = dot(x, x2);
            axpy(-c, x, x2);
            axpy(-c, sx, s2);
            const double n2 = std::sqrt(dot(x2, x2));
            bv.assign(1, x);
            sv.assign(1, sx);
            y_now.assign(1, 1.0);
            if (n2 > 1e-3) {
                for (double& v : x2) v /= n2;
                for (double& v : s2) v /= n2;
                bv.push_back(std::move(x2));
                sv.push_back(std::move(s2));
                y_now.push_back(0.0);
            }
            G.assign(bv.size(), std::vector<double>(bv.size(), 0.0));
            for (std::size_t i = 0; i < bv.size(); ++i)
                for (std::size_t j = 0; j < bv.size(); ++j) G[i][j] = dot(bv[i], sv[j]);
        }
        y_prev = y_now;
        y_prev.push_back(0.0);  // the correction vector added below
        // Correction vector (diagonal preconditioner), singlet projection, orthogonalisation.
        std::vector<double> t(Nd);
        for (std::size_t i = 0; i < Nd; ++i) {
            double d = hd[i] - theta;
            if (std::fabs(d) < 1e-4) d = d < 0 ? -1e-4 : 1e-4;
            t[i] = -res[i] / d;
        }
        if (symmetric) symmetrize(t, Na);
        for (int pass = 0; pass < 2; ++pass)
            for (const auto& b : bv) axpy(-dot(b, t), b, t);
        const double tn = std::sqrt(dot(t, t));
        if (tn < 1e-12) {
            r.converged = rn < 10.0 * opt.residual_tol;
            break;
        }
        for (double& v : t) v /= tn;
        bv.push_back(std::move(t));
        sv.emplace_back();
        H.apply(bv.back(), sv.back(), symmetric);
        const std::size_t k = bv.size() - 1;
        for (auto& row : G) row.push_back(0.0);
        G.emplace_back(k + 1, 0.0);
        for (std::size_t i = 0; i <= k; ++i) {
            G[i][k] = dot(bv[i], sv[k]);
            G[k][i] = G[i][k];
        }
    }

    r.electronic_energy = theta;
    r.total_energy = theta + nuclear_repulsion;
    r.correlation_energy = r.total_energy - rhf_energy;

    // Normalise the final vector and analyse it.
    const double xn = std::sqrt(dot(x, x));
    for (double& v : x) v /= xn;
    r.hf_weight = x[0] * x[0];
    r.s2 = spin_squared(x, H.alpha(), H.beta());

    // Leading determinants.
    {
        std::vector<std::pair<double, std::size_t>> top;
        const std::size_t keep = 10;
        for (std::size_t i = 0; i < Nd; ++i) {
            const double a = std::fabs(x[i]);
            if (top.size() < keep || a > top.back().first) {
                top.emplace_back(a, i);
                std::sort(top.begin(), top.end(), [](const auto& u, const auto& v) {
                    return u.first > v.first || (u.first == v.first && u.second < v.second);
                });
                if (top.size() > keep) top.pop_back();
            }
        }
        for (const auto& [a, i] : top) {
            if (a < 1e-12) break;
            FCIDeterminant d;
            d.coefficient = x[i];
            const std::uint16_t* oa = H.alpha().occupation(i / Nb);
            const std::uint16_t* ob = H.beta().occupation(i % Nb);
            for (int k = 0; k < na; ++k) {
                d.alpha.push_back(oa[k] + 1);
                if (oa[k] >= na) ++d.excitation;
            }
            for (int k = 0; k < nb; ++k) {
                d.beta.push_back(ob[k] + 1);
                if (ob[k] >= nb) ++d.excitation;
            }
            r.leading.push_back(d);
        }
    }

    // One-particle density matrix and natural occupation numbers.
    {
        const std::size_t n = static_cast<std::size_t>(mo.n);
        std::vector<Matrix> parts(static_cast<std::size_t>(threads), Matrix(n, n));
        const auto& ex = H.alpha_excitations();
        std::atomic<std::size_t> next{0};
        run_parallel(threads, [&](int th) {
            Matrix& g = parts[static_cast<std::size_t>(th)];
            for (std::size_t J; (J = next++) < Na;) {
                const double* cj = x.data() + J * Nb;
                for (const Excitation* e = ex.begin(J); e != ex.end(J); ++e) {
                    const double* ci = x.data() + static_cast<std::size_t>(e->target) * Nb;
                    double s = 0.0;
                    for (std::size_t b = 0; b < Nb; ++b) s += ci[b] * cj[b];
                    g(e->p, e->q) += e->sign * s;
                }
            }
        });
        Matrix gamma(n, n);
        for (const Matrix& g : parts)
            for (std::size_t p = 0; p < n; ++p)
                for (std::size_t q = 0; q < n; ++q) gamma(p, q) += g(p, q);
        // beta part: equal to the alpha part for a symmetric vector (n_alpha == n_beta)
        Matrix total(n, n);
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) total(p, q) = 0.5 * (gamma(p, q) + gamma(q, p)) * 2.0;
        const EigenResult ne = symmetric_eigen(total);
        r.natural_occupations.assign(ne.values.rbegin(), ne.values.rend());
    }
    r.ci_seconds = seconds_since(t0);
    return r;
}

}  // namespace macqc
