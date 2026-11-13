// General-order spin-orbital coupled cluster (CCSD, CCSDT, CCSDTQ) for closed-shell
// references.
//
// Amplitudes t(A;I) are stored once per pair of sorted spin-orbital sets (A virtual,
// I occupied), blocked by the number of alpha orbitals in each set; within a block the
// alpha and the beta subsets are addressed by the combinatorial number system. Every
// tensor here has four index groups, each antisymmetric and packed the same way:
//
//   op  virtual lines into the Hamiltonian vertex that are still open (annihilators),
//   oh  occupied lines out of the vertex that are still open (creators),
//   ep  external particle lines (creators),   eh  external hole lines (annihilators).
//
// Amplitudes and residuals are tensors with empty op and oh. A diagram (from
// cc_diagrams.inc, see backend/tools/cc_diagrams.py) starts from a Hamiltonian block and
// contracts one cluster vertex T_n at a time:
//
//   Y(op', oh', ep.Pt, eh.Ht) = sum_{Cp, Ch} X(Cp.op', Ch.oh', ep, eh) T_n(Cp.Pt; Ch.Ht),
//
// where "." is the concatenation of two sets (the value at the sorted union carries the
// sign of sorting the concatenation) and the sums run over sorted sets. Each step is one
// matrix product per spin case: X is gathered into A [rows x k], T_n into B [k x cols],
// and the product is scattered into Y with the merge signs. Diagrams that share their
// first vertices share the intermediates (a trie, evaluated depth first).
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/cc.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/linalg.h"
#include "macqc/parallel.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;

struct DiagramSpec {
    int target;
    char vertex;
    int x, y;
    double coef;
    int k;
    int cl[4][3];
};

const DiagramSpec kDiagrams[] = {
#include "cc_diagrams.inc"
};

// ---------------------------------------------------------------------------------------
// Combinatorics
// ---------------------------------------------------------------------------------------

constexpr int kMaxN = 2048;
constexpr int kMaxK = 8;

const std::vector<std::uint64_t>& binomial_table() {
    static const std::vector<std::uint64_t> table = [] {
        std::vector<std::uint64_t> c(static_cast<std::size_t>(kMaxN + 1) * (kMaxK + 1), 0);
        for (int n = 0; n <= kMaxN; ++n) {
            c[static_cast<std::size_t>(n) * (kMaxK + 1)] = 1;
            for (int k = 1; k <= std::min(n, kMaxK); ++k)
                c[static_cast<std::size_t>(n) * (kMaxK + 1) + k] =
                    c[static_cast<std::size_t>(n - 1) * (kMaxK + 1) + k - 1] +
                    (k <= n - 1 ? c[static_cast<std::size_t>(n - 1) * (kMaxK + 1) + k] : 0);
        }
        return c;
    }();
    return table;
}

inline std::size_t choose(int n, int k) {
    if (k < 0 || n < 0 || k > n) return 0;
    if (n > kMaxN || k > kMaxK) throw std::runtime_error("coupled cluster: orbital space too large");
    return static_cast<std::size_t>(binomial_table()[static_cast<std::size_t>(n) * (kMaxK + 1) + k]);
}

// All g-subsets of {0..n-1} in colex order: rank(S) = sum_t C(s_t, t + 1).
struct SubsetList {
    int n = 0, g = 0;
    std::size_t count = 0;
    std::vector<std::uint16_t> e;  // e[r * g + t]
    const std::uint16_t* operator[](std::size_t r) const { return e.data() + r * static_cast<std::size_t>(g); }
};

std::mutex g_cache_mutex;

const SubsetList& subsets(int n, int g) {
    static std::map<std::pair<int, int>, SubsetList> cache;
    std::lock_guard<std::mutex> lock(g_cache_mutex);
    auto it = cache.find({n, g});
    if (it != cache.end()) return it->second;
    SubsetList s;
    s.n = n;
    s.g = g;
    s.count = choose(n, g);
    s.e.resize(s.count * static_cast<std::size_t>(g));
    std::vector<int> c(static_cast<std::size_t>(g));
    for (int t = 0; t < g; ++t) c[static_cast<std::size_t>(t)] = t;
    for (std::size_t r = 0; r < s.count; ++r) {
        for (int t = 0; t < g; ++t) s.e[r * g + t] = static_cast<std::uint16_t>(c[static_cast<std::size_t>(t)]);
        // colex successor: increase the first element that can move, reset the ones below it
        for (int t = 0; t < g; ++t) {
            const int limit = (t + 1 < g) ? c[static_cast<std::size_t>(t) + 1] : n;
            if (c[static_cast<std::size_t>(t)] + 1 < limit) {
                ++c[static_cast<std::size_t>(t)];
                for (int u = 0; u < t; ++u) c[static_cast<std::size_t>(u)] = u;
                break;
            }
        }
    }
    return cache.emplace(std::make_pair(n, g), std::move(s)).first->second;
}

// Union of an g1-subset and a g2-subset of {0..n-1}: sign * (rank of the union + 1), or 0
// if they overlap. The sign is that of sorting the concatenation S1.S2.
const std::vector<std::int32_t>& union_table(int n, int g1, int g2) {
    static std::map<std::array<int, 3>, std::vector<std::int32_t>> cache;
    {
        std::lock_guard<std::mutex> lock(g_cache_mutex);
        auto it = cache.find({n, g1, g2});
        if (it != cache.end()) return it->second;
    }
    const SubsetList& s1 = subsets(n, g1);
    const SubsetList& s2 = subsets(n, g2);
    if (choose(n, g1 + g2) >= static_cast<std::size_t>(INT32_MAX))
        throw std::runtime_error("coupled cluster: orbital space too large for the index tables");
    std::vector<std::int32_t> table(s1.count * s2.count, 0);
    std::vector<int> u(static_cast<std::size_t>(g1 + g2));
    for (std::size_t r1 = 0; r1 < s1.count; ++r1) {
        const std::uint16_t* a = s1[r1];
        for (std::size_t r2 = 0; r2 < s2.count; ++r2) {
            const std::uint16_t* b = s2[r2];
            int i = 0, j = 0, k = 0, inversions = 0;
            bool overlap = false;
            while (i < g1 || j < g2) {
                if (j == g2 || (i < g1 && a[i] < b[j])) {
                    u[static_cast<std::size_t>(k++)] = a[i++];
                } else if (i == g1 || b[j] < a[i]) {
                    inversions += g1 - i;  // b[j] moves past the remaining elements of a
                    u[static_cast<std::size_t>(k++)] = b[j++];
                } else {
                    overlap = true;
                    break;
                }
            }
            if (overlap) continue;
            std::size_t rank = 0;
            for (int t = 0; t < g1 + g2; ++t) rank += choose(u[static_cast<std::size_t>(t)], t + 1);
            const std::int32_t v = static_cast<std::int32_t>(rank + 1);
            table[r1 * s2.count + r2] = (inversions % 2) ? -v : v;
        }
    }
    std::lock_guard<std::mutex> lock(g_cache_mutex);
    return cache.emplace(std::array<int, 3>{n, g1, g2}, std::move(table)).first->second;
}

// A set of spin orbitals of one space: a alpha and b beta orbitals, canonical order
// alpha (sorted) then beta (sorted); index = rank_alpha * C(n_beta, b) + rank_beta.
inline std::size_t part_count(const int* dims, int a, int b) { return choose(dims[0], a) * choose(dims[1], b); }

// Concatenation of two parts (a1,b1).(a2,b2) of the same space: sign * (index + 1) of
// the union, 0 if they overlap. Indexed [i1 * count2 + i2].
const std::vector<std::int64_t>& split_map(const int* dims, int a1, int b1, int a2, int b2) {
    static std::map<std::array<int, 6>, std::vector<std::int64_t>> cache;
    const std::array<int, 6> key{dims[0], dims[1], a1, b1, a2, b2};
    {
        std::lock_guard<std::mutex> lock(g_cache_mutex);
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
    }
    const std::vector<std::int32_t>& ua = union_table(dims[0], a1, a2);
    const std::vector<std::int32_t>& ub = union_table(dims[1], b1, b2);
    const std::size_t ca1 = choose(dims[0], a1), cb1 = choose(dims[1], b1);
    const std::size_t ca2 = choose(dims[0], a2), cb2 = choose(dims[1], b2);
    const std::size_t cbu = choose(dims[1], b1 + b2);
    const int interleave = ((b1 * a2) % 2) ? -1 : 1;  // move the alpha part of set 2 past the beta part of set 1
    std::vector<std::int64_t> map(ca1 * cb1 * ca2 * cb2, 0);
    for (std::size_t ra1 = 0; ra1 < ca1; ++ra1)
        for (std::size_t rb1 = 0; rb1 < cb1; ++rb1) {
            const std::size_t i1 = ra1 * cb1 + rb1;
            for (std::size_t ra2 = 0; ra2 < ca2; ++ra2) {
                const std::int32_t xa = ua[ra1 * ca2 + ra2];
                if (xa == 0) continue;
                for (std::size_t rb2 = 0; rb2 < cb2; ++rb2) {
                    const std::int32_t xb = ub[rb1 * cb2 + rb2];
                    if (xb == 0) continue;
                    const std::size_t idx = static_cast<std::size_t>(std::abs(xa) - 1) * cbu +
                                            static_cast<std::size_t>(std::abs(xb) - 1);
                    const int sign = (xa < 0 ? -1 : 1) * (xb < 0 ? -1 : 1) * interleave;
                    map[i1 * (ca2 * cb2) + ra2 * cb2 + rb2] = sign * static_cast<std::int64_t>(idx + 1);
                }
            }
        }
    std::lock_guard<std::mutex> lock(g_cache_mutex);
    return cache.emplace(key, std::move(map)).first->second;
}

// ---------------------------------------------------------------------------------------
// Packed tensors with four groups: op (vir), oh (occ), ep (vir), eh (occ)
// ---------------------------------------------------------------------------------------

struct Spaces {
    int occ[2] = {0, 0};  // alpha, beta
    int vir[2] = {0, 0};
    const int* dims(int group) const { return (group == 0 || group == 2) ? vir : occ; }
};

struct GBlock {
    std::array<int, 4> alpha{};
    std::array<std::size_t, 4> count{};
    std::size_t offset = 0, size = 0;
};

struct GTensor {
    std::array<int, 4> size{};
    std::vector<GBlock> blocks;
    std::vector<int> lookup;
    std::vector<double> data;
    std::size_t total = 0;

    int find(int a0, int a1, int a2, int a3) const {
        if (a0 < 0 || a1 < 0 || a2 < 0 || a3 < 0 || a0 > size[0] || a1 > size[1] || a2 > size[2] || a3 > size[3])
            return -1;
        return lookup[static_cast<std::size_t>(((a0 * (size[1] + 1) + a1) * (size[2] + 1) + a2) * (size[3] + 1) + a3)];
    }
};

GTensor make_tensor(const Spaces& sp, std::array<int, 4> size, bool allocate = true) {
    GTensor t;
    t.size = size;
    t.lookup.assign(static_cast<std::size_t>((size[0] + 1) * (size[1] + 1) * (size[2] + 1) * (size[3] + 1)), -1);
    for (int a0 = 0; a0 <= size[0]; ++a0)
        for (int a1 = 0; a1 <= size[1]; ++a1)
            for (int a2 = 0; a2 <= size[2]; ++a2)
                for (int a3 = 0; a3 <= size[3]; ++a3) {
                    if (a1 + a2 != a0 + a3) continue;  // alpha creators = alpha annihilators
                    const std::array<int, 4> a{a0, a1, a2, a3};
                    GBlock b;
                    b.alpha = a;
                    std::size_t n = 1;
                    for (int g = 0; g < 4; ++g) {
                        b.count[static_cast<std::size_t>(g)] = part_count(sp.dims(g), a[g], size[static_cast<std::size_t>(g)] - a[g]);
                        n *= b.count[static_cast<std::size_t>(g)];
                    }
                    if (n == 0) continue;
                    b.offset = t.total;
                    b.size = n;
                    t.total += n;
                    t.lookup[static_cast<std::size_t>(((a0 * (size[1] + 1) + a1) * (size[2] + 1) + a2) * (size[3] + 1) + a3)] =
                        static_cast<int>(t.blocks.size());
                    t.blocks.push_back(b);
                }
    if (allocate) t.data.assign(t.total, 0.0);
    return t;
}

// ---------------------------------------------------------------------------------------
// One contraction step
// ---------------------------------------------------------------------------------------

struct Workspace {
    std::size_t budget = std::size_t(1) << 24;  // elements per buffer (A, B, M)
    int threads = 1;
    std::vector<double> A, B, M;
};

struct Profile {
    double gather_a = 0, gather_b = 0, gemm = 0, scatter = 0, setup = 0;
    double flops = 0;
    ~Profile() {
        if (std::getenv("PSILAB_CC_PROFILE"))
            std::fprintf(stderr, "profile: setup %.2f gatherA %.2f gatherB %.2f gemm %.2f (%.2f GF/s) scatter %.2f\n", setup,
                         gather_a, gather_b, gemm, 2e-9 * flops / std::max(gemm, 1e-9), scatter);
    }
};
Profile g_prof;
inline double now_s() { return std::chrono::duration<double>(Clock::now().time_since_epoch()).count(); }

struct StepStats {
    double flops = 0.0;   // multiply-adds
    double moves = 0.0;   // gathered + scattered elements
    double peak = 0.0;    // largest intermediate (elements)
};

struct RowCfg {
    int aop, aoh, aep, aeh;
    std::size_t nop, noh, nep, neh, count, start;
};
struct ColCfg {
    int aPt, aHt;
    std::size_t nPt, nHt, count, start;
};
struct KCfg {
    int aCp, aCh;
    std::size_t nCp, nCh, count, start;
};

inline std::pair<int, int> alpha_range(int s, const int* dims) {
    return {std::max(0, s - dims[1]), std::min(s, dims[0])};
}

// Y += w * (X contracted with T over cp particle and ch hole lines). Y == nullptr: count only.
void contract_step(const Spaces& sp, const GTensor& X, const GTensor& T, int n, int cp, int ch, GTensor* Y, double w,
                   Workspace& ws, StepStats* stats) {
    const int sop = X.size[0], soh = X.size[1], sep = X.size[2], seh = X.size[3];
    const int sop2 = sop - cp, soh2 = soh - ch, sPt = n - cp, sHt = n - ch;
    const int* V = sp.vir;
    const int* O = sp.occ;
    const std::array<int, 4> ysize{sop2, soh2, sep + sPt, seh + sHt};
    GTensor shape;
    const GTensor* Ys = Y;
    if (!Y) {
        shape = make_tensor(sp, ysize, false);
        Ys = &shape;
    } else if (Y->size != ysize) {
        throw std::logic_error("contract_step: output shape");
    }

    std::vector<RowCfg> rows;
    std::vector<ColCfg> cols;
    std::vector<KCfg> ks;
    {
        const auto rop = alpha_range(sop2, V), roh = alpha_range(soh2, O), rep = alpha_range(sep, V),
                   reh = alpha_range(seh, O);
        for (int a0 = rop.first; a0 <= rop.second; ++a0)
            for (int a1 = roh.first; a1 <= roh.second; ++a1)
                for (int a2 = rep.first; a2 <= rep.second; ++a2)
                    for (int a3 = reh.first; a3 <= reh.second; ++a3) {
                        RowCfg r{a0, a1, a2, a3, part_count(V, a0, sop2 - a0), part_count(O, a1, soh2 - a1),
                                 part_count(V, a2, sep - a2), part_count(O, a3, seh - a3), 0, 0};
                        r.count = r.nop * r.noh * r.nep * r.neh;
                        if (r.count) rows.push_back(r);
                    }
        const auto rPt = alpha_range(sPt, V), rHt = alpha_range(sHt, O);
        for (int a0 = rPt.first; a0 <= rPt.second; ++a0)
            for (int a1 = rHt.first; a1 <= rHt.second; ++a1) {
                ColCfg c{a0, a1, part_count(V, a0, sPt - a0), part_count(O, a1, sHt - a1), 0, 0};
                c.count = c.nPt * c.nHt;
                if (c.count) cols.push_back(c);
            }
        const auto rCp = alpha_range(cp, V), rCh = alpha_range(ch, O);
        for (int a0 = rCp.first; a0 <= rCp.second; ++a0)
            for (int a1 = rCh.first; a1 <= rCh.second; ++a1) {
                KCfg k{a0, a1, part_count(V, a0, cp - a0), part_count(O, a1, ch - a1), 0, 0};
                k.count = k.nCp * k.nCh;
                if (k.count) ks.push_back(k);
            }
    }
    auto xblock = [&](const RowCfg& r, const KCfg& k) { return X.find(k.aCp + r.aop, k.aCh + r.aoh, r.aep, r.aeh); };
    auto tblock = [&](const KCfg& k, const ColCfg& c) { return T.find(0, 0, k.aCp + c.aPt, k.aCh + c.aHt); };
    auto yblock = [&](const RowCfg& r, const ColCfg& c) { return Ys->find(r.aop, r.aoh, r.aep + c.aPt, r.aeh + c.aHt); };

    // Group the spin cases by delta = (alpha holes - alpha particles) of the contracted lines.
    std::map<int, std::array<std::vector<int>, 3>> groups;  // delta -> rows, ks, cols
    for (std::size_t i = 0; i < rows.size(); ++i)
        groups[rows[i].aop + rows[i].aeh - rows[i].aoh - rows[i].aep][0].push_back(static_cast<int>(i));
    for (std::size_t i = 0; i < ks.size(); ++i) groups[ks[i].aCh - ks[i].aCp][1].push_back(static_cast<int>(i));
    for (std::size_t i = 0; i < cols.size(); ++i) groups[cols[i].aPt - cols[i].aHt][2].push_back(static_cast<int>(i));

    for (auto& entry : groups) {
        std::vector<int> R, K, C;
        // Keep only spin cases that take part in at least one valid (row, k, col) triple.
        {
            const auto& r0 = entry.second[0];
            const auto& k0 = entry.second[1];
            const auto& c0 = entry.second[2];
            std::vector<char> ur(r0.size(), 0), uk(k0.size(), 0), uc(c0.size(), 0);
            for (std::size_t i = 0; i < r0.size(); ++i)
                for (std::size_t j = 0; j < k0.size(); ++j) {
                    if (xblock(rows[static_cast<std::size_t>(r0[i])], ks[static_cast<std::size_t>(k0[j])]) < 0) continue;
                    for (std::size_t l = 0; l < c0.size(); ++l) {
                        const ColCfg& c = cols[static_cast<std::size_t>(c0[l])];
                        if (tblock(ks[static_cast<std::size_t>(k0[j])], c) < 0) continue;
                        if (yblock(rows[static_cast<std::size_t>(r0[i])], c) < 0) continue;
                        ur[i] = uk[j] = uc[l] = 1;
                    }
                }
            for (std::size_t i = 0; i < r0.size(); ++i)
                if (ur[i]) R.push_back(r0[i]);
            for (std::size_t j = 0; j < k0.size(); ++j)
                if (uk[j]) K.push_back(k0[j]);
            for (std::size_t l = 0; l < c0.size(); ++l)
                if (uc[l]) C.push_back(c0[l]);
        }
        if (R.empty() || K.empty() || C.empty()) continue;
        std::size_t nrows = 0, nk = 0, ncols = 0;
        for (int i : R) {
            rows[static_cast<std::size_t>(i)].start = nrows;
            nrows += rows[static_cast<std::size_t>(i)].count;
        }
        for (int i : K) {
            ks[static_cast<std::size_t>(i)].start = nk;
            nk += ks[static_cast<std::size_t>(i)].count;
        }
        for (int i : C) {
            cols[static_cast<std::size_t>(i)].start = ncols;
            ncols += cols[static_cast<std::size_t>(i)].count;
        }
        if (stats) {
            stats->flops += static_cast<double>(nrows) * nk * ncols;
            stats->moves += static_cast<double>(nrows) * nk + static_cast<double>(nk) * ncols +
                            static_cast<double>(nrows) * ncols;
        }
        if (!Y) continue;

        const std::size_t budget = std::max<std::size_t>(ws.budget, 1024);
        const std::size_t rchunk = std::min(nrows, std::max<std::size_t>(1, budget / nk));
        const std::size_t cchunk = std::min(ncols, std::max<std::size_t>(1, budget / std::max(nk, rchunk)));
        const int threads = std::max(1, ws.threads);

        // Gather A rows [r0, r1).
        auto gather_a = [&](std::size_t r0, std::size_t r1) {
            ws.A.resize((r1 - r0) * nk);
            for (int ri : R) {
                const RowCfg& r = rows[static_cast<std::size_t>(ri)];
                const std::size_t lo = std::max(r0, r.start), hi = std::min(r1, r.start + r.count);
                if (lo >= hi) continue;
                for (int ki : K) {
                    const KCfg& k = ks[static_cast<std::size_t>(ki)];
                    const int xb = xblock(r, k);
                    if (xb < 0) {
                        for (std::size_t g = lo; g < hi; ++g)
                            std::fill_n(ws.A.data() + (g - r0) * nk + k.start, k.count, 0.0);
                        continue;
                    }
                    const GBlock& b = X.blocks[static_cast<std::size_t>(xb)];
                    const std::vector<std::int64_t>& mop = split_map(V, k.aCp, cp - k.aCp, r.aop, sop2 - r.aop);
                    const std::vector<std::int64_t>& moh = split_map(O, k.aCh, ch - k.aCh, r.aoh, soh2 - r.aoh);
                    const double* xd = X.data.data() + b.offset;
                    const std::size_t xoh = b.count[1], xep = b.count[2], xeh = b.count[3];
                    run_parallel(hi - lo > 64 ? threads : 1, [&](int t) {
                        for (std::size_t g = lo + static_cast<std::size_t>(t); g < hi; g += static_cast<std::size_t>(threads > 1 && hi - lo > 64 ? threads : 1)) {
                            std::size_t l = g - r.start;
                            const std::size_t ieh = l % r.neh;
                            l /= r.neh;
                            const std::size_t iep = l % r.nep;
                            l /= r.nep;
                            const std::size_t ioh = l % r.noh;
                            const std::size_t iop = l / r.noh;
                            double* a = ws.A.data() + (g - r0) * nk + k.start;
                            for (std::size_t c1 = 0; c1 < k.nCp; ++c1) {
                                const std::int64_t e1 = mop[c1 * r.nop + iop];
                                if (e1 == 0) {
                                    std::fill_n(a + c1 * k.nCh, k.nCh, 0.0);
                                    continue;
                                }
                                const std::size_t op = static_cast<std::size_t>(std::llabs(e1) - 1);
                                const double s1 = e1 < 0 ? -1.0 : 1.0;
                                for (std::size_t c2 = 0; c2 < k.nCh; ++c2) {
                                    const std::int64_t e2 = moh[c2 * r.noh + ioh];
                                    if (e2 == 0) {
                                        a[c1 * k.nCh + c2] = 0.0;
                                        continue;
                                    }
                                    const std::size_t oh = static_cast<std::size_t>(std::llabs(e2) - 1);
                                    const double s = e2 < 0 ? -s1 : s1;
                                    a[c1 * k.nCh + c2] = s * xd[((op * xoh + oh) * xep + iep) * xeh + ieh];
                                }
                            }
                        }
                    });
                }
            }
        };
        // Gather B columns [c0, c1).
        auto gather_b = [&](std::size_t c0, std::size_t c1) {
            const std::size_t width = c1 - c0;
            ws.B.resize(nk * width);
            for (int ki : K) {
                const KCfg& k = ks[static_cast<std::size_t>(ki)];
                for (int ci : C) {
                    const ColCfg& c = cols[static_cast<std::size_t>(ci)];
                    const std::size_t lo = std::max(c0, c.start), hi = std::min(c1, c.start + c.count);
                    if (lo >= hi) continue;
                    const int tb = tblock(k, c);
                    if (tb < 0) {
                        for (std::size_t kk = 0; kk < k.count; ++kk)
                            std::fill(ws.B.begin() + static_cast<std::ptrdiff_t>((k.start + kk) * width + (lo - c0)),
                                      ws.B.begin() + static_cast<std::ptrdiff_t>((k.start + kk) * width + (hi - c0)), 0.0);
                        continue;
                    }
                    const GBlock& b = T.blocks[static_cast<std::size_t>(tb)];
                    const std::vector<std::int64_t>& mP = split_map(V, k.aCp, cp - k.aCp, c.aPt, sPt - c.aPt);
                    const std::vector<std::int64_t>& mH = split_map(O, k.aCh, ch - k.aCh, c.aHt, sHt - c.aHt);
                    const double* td = T.data.data() + b.offset;
                    const std::size_t tH = b.count[3];
                    const int nt = k.count > 16 ? threads : 1;
                    run_parallel(nt, [&](int t) {
                        for (std::size_t kk = static_cast<std::size_t>(t); kk < k.count; kk += static_cast<std::size_t>(nt)) {
                            const std::size_t c1i = kk / k.nCh, c2i = kk % k.nCh;
                            double* out = ws.B.data() + (k.start + kk) * width;
                            for (std::size_t g = lo; g < hi; ++g) {
                                const std::size_t l = g - c.start;
                                const std::size_t iPt = l / c.nHt, iHt = l % c.nHt;
                                const std::int64_t e1 = mP[c1i * c.nPt + iPt];
                                const std::int64_t e2 = mH[c2i * c.nHt + iHt];
                                if (e1 == 0 || e2 == 0) {
                                    out[g - c0] = 0.0;
                                    continue;
                                }
                                const double s = ((e1 < 0) != (e2 < 0)) ? -1.0 : 1.0;
                                out[g - c0] = s * td[static_cast<std::size_t>(std::llabs(e1) - 1) * tH +
                                                     static_cast<std::size_t>(std::llabs(e2) - 1)];
                            }
                        }
                    });
                }
            }
        };
        // Scatter M (rows [r0,r1) x cols [c0,c1)) into Y with weight w.
        auto scatter = [&](std::size_t r0, std::size_t r1, std::size_t c0, std::size_t c1) {
            const std::size_t width = c1 - c0;
            const bool parallel_ok = ysize[2] > 0 && (r1 - r0) * width > 4096;
            const int nt = parallel_ok ? threads : 1;
            run_parallel(nt, [&](int t) {
                for (int ri : R) {
                    const RowCfg& r = rows[static_cast<std::size_t>(ri)];
                    const std::size_t rlo = std::max(r0, r.start), rhi = std::min(r1, r.start + r.count);
                    if (rlo >= rhi) continue;
                    for (int ci : C) {
                        const ColCfg& c = cols[static_cast<std::size_t>(ci)];
                        const std::size_t clo = std::max(c0, c.start), chi = std::min(c1, c.start + c.count);
                        if (clo >= chi) continue;
                        const int yb = yblock(r, c);
                        if (yb < 0) continue;
                        const GBlock& b = Y->blocks[static_cast<std::size_t>(yb)];
                        const std::vector<std::int64_t>& mEP = split_map(V, r.aep, sep - r.aep, c.aPt, sPt - c.aPt);
                        const std::vector<std::int64_t>& mEH = split_map(O, r.aeh, seh - r.aeh, c.aHt, sHt - c.aHt);
                        double* yd = Y->data.data() + b.offset;
                        const std::size_t nEP = b.count[2], nEH = b.count[3];
                        for (std::size_t g = rlo; g < rhi; ++g) {
                            std::size_t l = g - r.start;
                            const std::size_t ieh = l % r.neh;
                            l /= r.neh;
                            const std::size_t iep = l % r.nep;
                            l /= r.nep;
                            const std::size_t base = l * nEP;  // l = op' * noh' + oh'
                            const double* m = ws.M.data() + (g - r0) * width;
                            const std::int64_t* rowEP = mEP.data() + iep * c.nPt;
                            const std::int64_t* rowEH = mEH.data() + ieh * c.nHt;
                            std::size_t cc = clo - c.start;
                            std::size_t iPt = cc / c.nHt, iHt = cc % c.nHt;
                            std::size_t col = clo;
                            while (col < chi) {
                                const std::int64_t e1 = rowEP[iPt];
                                const std::size_t run = std::min(c.nHt - iHt, chi - col);
                                if (e1 != 0) {
                                    const std::size_t ep = static_cast<std::size_t>(std::llabs(e1) - 1);
                                    if (nt == 1 || static_cast<int>(ep % static_cast<std::size_t>(nt)) == t) {
                                        const double s1 = e1 < 0 ? -w : w;
                                        double* yrow = yd + (base + ep) * nEH;
                                        const double* mm = m + (col - c0);
                                        for (std::size_t h = 0; h < run; ++h) {
                                            const std::int64_t e2 = rowEH[iHt + h];
                                            if (e2 == 0) continue;
                                            const double v = e2 < 0 ? -s1 * mm[h] : s1 * mm[h];
                                            yrow[static_cast<std::size_t>(std::llabs(e2) - 1)] += v;
                                        }
                                    }
                                }
                                col += run;
                                iHt = 0;
                                ++iPt;
                            }
                        }
                    }
                }
            });
        };

        const bool b_once = (cchunk == ncols);
        double tp = now_s();
        if (b_once) gather_b(0, ncols);
        g_prof.gather_b += now_s() - tp;
        for (std::size_t r0 = 0; r0 < nrows; r0 += rchunk) {
            check_cancelled();
            const std::size_t r1 = std::min(nrows, r0 + rchunk);
            tp = now_s();
            gather_a(r0, r1);
            g_prof.gather_a += now_s() - tp;
            for (std::size_t c0 = 0; c0 < ncols; c0 += cchunk) {
                const std::size_t c1 = std::min(ncols, c0 + cchunk);
                tp = now_s();
                if (!b_once) gather_b(c0, c1);
                g_prof.gather_b += now_s() - tp;
                tp = now_s();
                ws.M.resize((r1 - r0) * (c1 - c0));
                gemm(false, false, r1 - r0, c1 - c0, nk, 1.0, ws.A.data(), nk, ws.B.data(), c1 - c0, 0.0,
                     ws.M.data(), c1 - c0, threads);
                g_prof.gemm += now_s() - tp;
                g_prof.flops += static_cast<double>(r1 - r0) * (c1 - c0) * nk;
                tp = now_s();
                scatter(r0, r1, c0, c1);
                g_prof.scatter += now_s() - tp;
            }
        }
    }
    if (stats) stats->peak = std::max(stats->peak, static_cast<double>(Ys->total));
}

// ---------------------------------------------------------------------------------------
// Diagram trie
// ---------------------------------------------------------------------------------------

struct Node {
    int n = 0, cp = 0, ch = 0;   // cluster vertex contracted at this node (n = 0: Hamiltonian)
    std::array<int, 4> size{};   // shape of the intermediate after this node
    std::vector<int> children;
    int target = -1;             // leaf: residual rank
    double coef = 0.0;
};

struct Root {
    int vertex = 0;  // 0 Fock, 1 two-electron
    int x = 0, y = 0;
    int node = 0;
};

struct Trie {
    std::vector<Node> nodes;
    std::vector<Root> roots;
};

Trie build_trie(int max_rank) {
    Trie t;
    for (const DiagramSpec& d : kDiagrams) {
        if (d.target > max_rank) continue;
        bool ok = true;
        for (int i = 0; i < d.k; ++i)
            if (d.cl[i][0] > max_rank) ok = false;
        if (!ok) continue;
        const int vertex = d.vertex == 'F' ? 0 : 1;
        const int c = vertex == 0 ? 1 : 2;
        auto it = std::find_if(t.roots.begin(), t.roots.end(),
                               [&](const Root& r) { return r.vertex == vertex && r.x == d.x && r.y == d.y; });
        if (it == t.roots.end()) {
            Root r;
            r.vertex = vertex;
            r.x = d.x;
            r.y = d.y;
            r.node = static_cast<int>(t.nodes.size());
            Node root;
            root.size = {d.y, c - d.x, d.x, c - d.y};
            t.nodes.push_back(root);
            t.roots.push_back(r);
            it = t.roots.end() - 1;
        }
        int cur = it->node;
        for (int i = 0; i < d.k; ++i) {
            const int n = d.cl[i][0], cp = d.cl[i][1], ch = d.cl[i][2];
            int next = -1;
            for (int child : t.nodes[static_cast<std::size_t>(cur)].children) {
                const Node& c2 = t.nodes[static_cast<std::size_t>(child)];
                if (c2.n == n && c2.cp == cp && c2.ch == ch) next = child;
            }
            if (next < 0) {
                Node nn;
                nn.n = n;
                nn.cp = cp;
                nn.ch = ch;
                const auto& s = t.nodes[static_cast<std::size_t>(cur)].size;
                nn.size = {s[0] - cp, s[1] - ch, s[2] + n - cp, s[3] + n - ch};
                next = static_cast<int>(t.nodes.size());
                t.nodes.push_back(nn);
                t.nodes[static_cast<std::size_t>(cur)].children.push_back(next);
            }
            cur = next;
        }
        Node& leaf = t.nodes[static_cast<std::size_t>(cur)];
        if (leaf.size[0] != 0 || leaf.size[1] != 0 || leaf.size[2] != d.target || leaf.size[3] != d.target)
            throw std::logic_error("coupled cluster: inconsistent diagram table");
        if (leaf.target >= 0) throw std::logic_error("coupled cluster: duplicate diagram");
        leaf.target = d.target;
        leaf.coef = d.coef;
    }
    return t;
}

// ---------------------------------------------------------------------------------------
// The spin-orbital problem
// ---------------------------------------------------------------------------------------

struct Problem {
    Spaces sp;
    int o = 0, v = 0, max_rank = 2;
    Matrix fock;                       // spatial MO Fock matrix
    std::vector<double> eps;           // diagonal of the Fock matrix
    // Hamiltonian blocks [vertex][x][y]
    std::array<std::array<std::array<GTensor, 3>, 3>, 2> H;
    Trie trie;
};

// Spin orbitals of one part of a block: (MO index, spin) in canonical order.
void decode_part(const int* dims, int a, int b, std::size_t idx, int mo_offset, std::vector<std::pair<int, int>>& out) {
    const std::size_t cb = choose(dims[1], b);
    const std::size_t ra = idx / cb, rb = idx % cb;
    const SubsetList& sa = subsets(dims[0], a);
    const SubsetList& sb = subsets(dims[1], b);
    for (int t = 0; t < a; ++t) out.emplace_back(mo_offset + sa[ra][t], 0);
    for (int t = 0; t < b; ++t) out.emplace_back(mo_offset + sb[rb][t], 1);
}

void fill_hamiltonian(Problem& P, const MOIntegrals& mo) {
    const Spaces& sp = P.sp;
    const int o = P.o;
    auto eri = [&](int p, int q, int r, int s) { return mo(static_cast<std::size_t>(p), static_cast<std::size_t>(q),
                                                         static_cast<std::size_t>(r), static_cast<std::size_t>(s)); };
    // <pq|rs> = (pr|qs) for equal spins
    auto phys = [&](const std::pair<int, int>& p, const std::pair<int, int>& q, const std::pair<int, int>& r,
                    const std::pair<int, int>& s) {
        if (p.second != r.second || q.second != s.second) return 0.0;
        return eri(p.first, r.first, q.first, s.first);
    };
    for (int vertex = 0; vertex < 2; ++vertex) {
        const int c = vertex == 0 ? 1 : 2;
        for (int x = 0; x <= c; ++x)
            for (int y = 0; y <= c; ++y) {
                GTensor t = make_tensor(sp, {y, c - x, x, c - y});
                std::vector<std::pair<int, int>> op, oh, ep, eh, cre, ann;
                for (const GBlock& b : t.blocks) {
                    const auto& a = b.alpha;
                    std::size_t pos = b.offset;
                    for (std::size_t i0 = 0; i0 < b.count[0]; ++i0)
                        for (std::size_t i1 = 0; i1 < b.count[1]; ++i1)
                            for (std::size_t i2 = 0; i2 < b.count[2]; ++i2)
                                for (std::size_t i3 = 0; i3 < b.count[3]; ++i3, ++pos) {
                                    op.clear();
                                    oh.clear();
                                    ep.clear();
                                    eh.clear();
                                    decode_part(sp.vir, a[0], t.size[0] - a[0], i0, o, op);
                                    decode_part(sp.occ, a[1], t.size[1] - a[1], i1, 0, oh);
                                    decode_part(sp.vir, a[2], t.size[2] - a[2], i2, o, ep);
                                    decode_part(sp.occ, a[3], t.size[3] - a[3], i3, 0, eh);
                                    cre = ep;
                                    cre.insert(cre.end(), oh.begin(), oh.end());
                                    ann = op;
                                    ann.insert(ann.end(), eh.begin(), eh.end());
                                    double val;
                                    if (vertex == 0) {
                                        const auto& p = cre[0];
                                        const auto& q = ann[0];
                                        val = (p.second == q.second) ? P.fock(static_cast<std::size_t>(p.first),
                                                                              static_cast<std::size_t>(q.first))
                                                                     : 0.0;
                                        if (p == q) val = 0.0;  // the diagonal is in the denominators
                                    } else {
                                        val = phys(cre[0], cre[1], ann[0], ann[1]) - phys(cre[0], cre[1], ann[1], ann[0]);
                                    }
                                    t.data[pos] = val;
                                }
                }
                P.H[static_cast<std::size_t>(vertex)][static_cast<std::size_t>(x)][static_cast<std::size_t>(y)] = std::move(t);
            }
    }
}

// Orbital-energy sums of all subsets of size g (one spin) of a space: sums[rank].
std::vector<double> subset_energy_sums(int n, int g, const std::vector<double>& eps, int offset) {
    const SubsetList& s = subsets(n, g);
    std::vector<double> out(s.count, 0.0);
    for (std::size_t r = 0; r < s.count; ++r)
        for (int t = 0; t < g; ++t) out[r] += eps[static_cast<std::size_t>(offset + s[r][t])];
    return out;
}

// Calls f(block, index in block, D) for every amplitude of rank n.
template <class F>
void for_each_denominator(const Problem& P, const GTensor& t, F f) {
    for (const GBlock& b : t.blocks) {
        const int n = t.size[2];
        const int aP = b.alpha[2], aH = b.alpha[3];
        const std::vector<double> pa = subset_energy_sums(P.sp.vir[0], aP, P.eps, P.o);
        const std::vector<double> pb = subset_energy_sums(P.sp.vir[1], n - aP, P.eps, P.o);
        const std::vector<double> ha = subset_energy_sums(P.sp.occ[0], aH, P.eps, 0);
        const std::vector<double> hb = subset_energy_sums(P.sp.occ[1], n - aH, P.eps, 0);
        const std::size_t nH = b.count[3];
        const std::size_t cbP = choose(P.sp.vir[1], n - aP), cbH = choose(P.sp.occ[1], n - aH);
        for (std::size_t ip = 0; ip < b.count[2]; ++ip) {
            const double dp = pa[ip / cbP] + pb[ip % cbP];
            for (std::size_t ih = 0; ih < nH; ++ih) {
                const double d = dp - ha[ih / cbH] - hb[ih % cbH];
                f(b, ip * nH + ih, d);
            }
        }
    }
}

struct Evaluator {
    const Problem& P;
    const std::vector<GTensor>& T;  // index = rank (0 unused)
    std::vector<GTensor>& R;        // index = rank (0: energy)
    Workspace& ws;
    StepStats* stats = nullptr;
    bool dry = false;

    void node(int id, const GTensor& X) {
        const Node& nd = P.trie.nodes[static_cast<std::size_t>(id)];
        for (int child : nd.children) {
            const Node& c = P.trie.nodes[static_cast<std::size_t>(child)];
            const GTensor& Tn = T[static_cast<std::size_t>(c.n)];
            if (!dry && Tn.total == 0) continue;
            if (c.children.empty()) {
                contract_step(P.sp, X, Tn, c.n, c.cp, c.ch, dry ? nullptr : &R[static_cast<std::size_t>(c.target)],
                              c.coef, ws, stats);
            } else {
                GTensor Y = make_tensor(P.sp, c.size, !dry);
                contract_step(P.sp, X, Tn, c.n, c.cp, c.ch, dry ? nullptr : &Y, 1.0, ws, stats);
                node(child, Y);
            }
        }
    }

    void run() {
        for (const Root& r : P.trie.roots) {
            const GTensor& H = P.H[static_cast<std::size_t>(r.vertex)][static_cast<std::size_t>(r.x)][static_cast<std::size_t>(r.y)];
            const Node& root = P.trie.nodes[static_cast<std::size_t>(r.node)];
            if (root.target >= 0 && !dry) {  // the Hamiltonian block itself (f_ai, <ab||ij>)
                GTensor& out = R[static_cast<std::size_t>(root.target)];
                for (std::size_t i = 0; i < out.total; ++i) out.data[i] += root.coef * H.data[i];
            }
            node(r.node, H);
        }
    }
};

Problem make_problem(int o, int v, int max_rank, bool allocate) {
    Problem P;
    P.o = o;
    P.v = v;
    P.max_rank = max_rank;
    P.sp.occ[0] = P.sp.occ[1] = o;
    P.sp.vir[0] = P.sp.vir[1] = v;
    P.trie = build_trie(max_rank);
    if (!allocate)
        for (int vertex = 0; vertex < 2; ++vertex) {
            const int c = vertex == 0 ? 1 : 2;
            for (int x = 0; x <= c; ++x)
                for (int y = 0; y <= c; ++y)
                    P.H[static_cast<std::size_t>(vertex)][static_cast<std::size_t>(x)][static_cast<std::size_t>(y)] =
                        make_tensor(P.sp, {y, c - x, x, c - y}, false);
        }
    return P;
}

}  // namespace

// ---------------------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------------------

const char* cc_method_name(CCMethod m) {
    switch (m) {
        case CCMethod::CCSD: return "CCSD";
        case CCMethod::CCSD_T: return "CCSD(T)";
        case CCMethod::CCSDT: return "CCSDT";
        case CCMethod::CCSDTQ: return "CCSDTQ";
    }
    return "CC";
}

int cc_max_rank(CCMethod m) {
    switch (m) {
        case CCMethod::CCSD:
        case CCMethod::CCSD_T: return 2;
        case CCMethod::CCSDT: return 3;
        case CCMethod::CCSDTQ: return 4;
    }
    return 2;
}

double reference_energy_from_mo(const MOIntegrals& mo, int nocc, double nuclear_repulsion) {
    double e = nuclear_repulsion;
    for (int i = 0; i < nocc; ++i) {
        const std::size_t si = static_cast<std::size_t>(i);
        e += 2.0 * mo.h(si, si);
        for (int j = 0; j < nocc; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            e += 2.0 * mo(si, si, sj, sj) - mo(si, sj, sj, si);
        }
    }
    return e;
}

SpinOrbitalCCCost estimate_cc_spin_orbital(int nocc, int nvir, int max_rank) {
    Problem P = make_problem(nocc, nvir, max_rank, false);
    std::vector<GTensor> T(static_cast<std::size_t>(max_rank) + 1), R(static_cast<std::size_t>(max_rank) + 1);
    SpinOrbitalCCCost c;
    for (int n = 1; n <= max_rank; ++n) {
        T[static_cast<std::size_t>(n)] = make_tensor(P.sp, {0, 0, n, n}, false);
        c.amplitudes += static_cast<double>(T[static_cast<std::size_t>(n)].total);
    }
    for (const auto& a : P.H)
        for (const auto& b : a)
            for (const GTensor& t : b) c.integral_bytes += 8.0 * static_cast<double>(t.total);
    Workspace ws;
    StepStats stats;
    Evaluator ev{P, T, R, ws, &stats, true};
    ev.run();
    c.flops_per_iteration = stats.flops;
    c.scatter_per_iteration = stats.moves;
    c.vector_bytes = 8.0 * c.amplitudes;
    // Work buffers (A, B, M) up to 3 x 128 MB, the largest intermediate twice (nested).
    const double buffers = std::min(3.0 * 8.0 * static_cast<double>(std::size_t(1) << 24), 3.0 * 8.0 * stats.moves);
    c.work_bytes = buffers + 2.0 * 8.0 * stats.peak;
    c.fixed_bytes = c.integral_bytes + c.work_bytes + 2.0 * c.vector_bytes;
    return c;
}

CCResult run_cc_spin_orbital(const MOIntegrals& mo, int nocc, int max_rank, double nuclear_repulsion,
                             const CCOptions& opt, const CCProgress& progress) {
    if (max_rank < 2 || max_rank > 4) throw std::invalid_argument("coupled cluster: rank must be 2, 3 or 4");
    const int n = mo.n;
    const int o = nocc, v = n - nocc;
    if (o <= 0 || v <= 0) throw std::runtime_error("coupled cluster needs occupied and virtual orbitals");
    const auto t_start = Clock::now();
    report_progress(JobStage::CoupledCluster);

    Problem P = make_problem(o, v, max_rank, true);
    P.fock = Matrix(static_cast<std::size_t>(n), static_cast<std::size_t>(n));
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            double f = mo.h(static_cast<std::size_t>(p), static_cast<std::size_t>(q));
            for (int i = 0; i < o; ++i) {
                const std::size_t si = static_cast<std::size_t>(i);
                f += 2.0 * mo(static_cast<std::size_t>(p), static_cast<std::size_t>(q), si, si) -
                     mo(static_cast<std::size_t>(p), si, si, static_cast<std::size_t>(q));
            }
            P.fock(static_cast<std::size_t>(p), static_cast<std::size_t>(q)) = f;
        }
    P.eps.resize(static_cast<std::size_t>(n));
    for (int p = 0; p < n; ++p) P.eps[static_cast<std::size_t>(p)] = P.fock(static_cast<std::size_t>(p), static_cast<std::size_t>(p));
    fill_hamiltonian(P, mo);

    const int threads = resolve_threads(opt.threads);
    Workspace ws;
    ws.threads = threads;
    std::vector<GTensor> T(static_cast<std::size_t>(max_rank) + 1), R(static_cast<std::size_t>(max_rank) + 1);
    std::size_t nvec = 0;
    for (int r = 1; r <= max_rank; ++r) {
        T[static_cast<std::size_t>(r)] = make_tensor(P.sp, {0, 0, r, r});
        R[static_cast<std::size_t>(r)] = make_tensor(P.sp, {0, 0, r, r});
        nvec += T[static_cast<std::size_t>(r)].total;
    }
    R[0] = make_tensor(P.sp, {0, 0, 0, 0});

    // Work buffers: up to 1/8 of the memory limit each, at most 2^24 elements.
    ws.budget = static_cast<std::size_t>(std::clamp(opt.memory_limit_mb * 1024.0 * 1024.0 / 8.0 / 8.0, 65536.0,
                                                    static_cast<double>(std::size_t(1) << 24)));
    // DIIS: as many vector pairs as fit into what is left of the memory limit (3 to diis_size).
    const SpinOrbitalCCCost cost = estimate_cc_spin_orbital(o, v, max_rank);
    const double left = opt.memory_limit_mb * 1024.0 * 1024.0 - cost.fixed_bytes;
    const int diis_max = std::max(2, std::min(opt.diis_size, static_cast<int>(left / (2.0 * 8.0 * static_cast<double>(nvec)))));

    CCResult res;
    res.method = max_rank == 2 ? "CCSD" : (max_rank == 3 ? "CCSDT" : "CCSDTQ");
    res.engine = "spin-orbital, general order";
    res.occupied = o;
    res.virtuals = v;
    res.reference_energy = reference_energy_from_mo(mo, o, nuclear_repulsion);
    res.amplitude_counts.assign(static_cast<std::size_t>(max_rank) + 1, 0.0);
    for (int r = 1; r <= max_rank; ++r)
        res.amplitude_counts[static_cast<std::size_t>(r)] = static_cast<double>(T[static_cast<std::size_t>(r)].total);

    std::vector<std::vector<double>> diis_t, diis_e;
    std::vector<double> tvec(nvec), evec(nvec);
    double e_old = 0.0;
    for (int it = 1; it <= opt.max_iter; ++it) {
        check_cancelled();
        for (auto& r : R) std::fill(r.data.begin(), r.data.end(), 0.0);
        Evaluator ev{P, T, R, ws, nullptr, false};
        ev.run();
        const double e_corr = R[0].data[0];
        // Jacobi step: t_new = t - R / D, with R including D t.
        double dmax = 0.0;
        std::size_t pos = 0;
        for (int r = 1; r <= max_rank; ++r) {
            GTensor& t = T[static_cast<std::size_t>(r)];
            const GTensor& rr = R[static_cast<std::size_t>(r)];
            for_each_denominator(P, t, [&](const GBlock& b, std::size_t i, double d) {
                const std::size_t k = b.offset + i;
                const double full = rr.data[k] + d * t.data[k];
                const double delta = -full / d;
                dmax = std::max(dmax, std::fabs(delta));
                evec[pos + k] = delta;
                tvec[pos + k] = t.data[k] + delta;
            });
            pos += t.total;
        }
        CCIteration rec;
        rec.iteration = it;
        rec.correlation = e_corr;
        rec.energy = res.reference_energy + e_corr;
        rec.delta_e = e_corr - e_old;
        rec.residual = dmax;
        res.history.push_back(rec);
        res.iterations = it;
        report_progress(JobStage::CoupledCluster, it, rec.energy, dmax);
        if (progress) progress(rec);
        const bool done = it > 1 && dmax < opt.amplitude_tol && std::fabs(rec.delta_e) < opt.energy_tol;
        e_old = e_corr;
        res.correlation_energy = e_corr;
        if (done) {
            res.converged = true;
            break;
        }
        // DIIS extrapolation of the new amplitudes.
        diis_t.push_back(tvec);
        diis_e.push_back(evec);
        if (static_cast<int>(diis_t.size()) > diis_max) {
            diis_t.erase(diis_t.begin());
            diis_e.erase(diis_e.begin());
        }
        const std::size_t m = diis_t.size();
        std::vector<double> next = tvec;
        if (m >= 2) {
            Matrix B(m + 1, m + 1, -1.0);
            B(m, m) = 0.0;
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j <= i; ++j) {
                    double s = 0.0;
                    for (std::size_t k = 0; k < nvec; ++k) s += diis_e[i][k] * diis_e[j][k];
                    B(i, j) = B(j, i) = s;
                }
            // Solve B c = (0,..,0,-1) by Gaussian elimination with partial pivoting.
            std::vector<double> rhs(m + 1, 0.0);
            rhs[m] = -1.0;
            bool ok = true;
            const std::size_t N = m + 1;
            for (std::size_t col = 0; col < N && ok; ++col) {
                std::size_t piv = col;
                for (std::size_t r2 = col + 1; r2 < N; ++r2)
                    if (std::fabs(B(r2, col)) > std::fabs(B(piv, col))) piv = r2;
                if (std::fabs(B(piv, col)) < 1e-300) {
                    ok = false;
                    break;
                }
                if (piv != col) {
                    for (std::size_t c2 = 0; c2 < N; ++c2) std::swap(B(col, c2), B(piv, c2));
                    std::swap(rhs[col], rhs[piv]);
                }
                for (std::size_t r2 = col + 1; r2 < N; ++r2) {
                    const double f = B(r2, col) / B(col, col);
                    for (std::size_t c2 = col; c2 < N; ++c2) B(r2, c2) -= f * B(col, c2);
                    rhs[r2] -= f * rhs[col];
                }
            }
            if (ok) {
                std::vector<double> c(N, 0.0);
                for (std::size_t r2 = N; r2-- > 0;) {
                    double s = rhs[r2];
                    for (std::size_t c2 = r2 + 1; c2 < N; ++c2) s -= B(r2, c2) * c[c2];
                    c[r2] = s / B(r2, r2);
                }
                std::fill(next.begin(), next.end(), 0.0);
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t k = 0; k < nvec; ++k) next[k] += c[i] * diis_t[i][k];
            }
        }
        pos = 0;
        for (int r = 1; r <= max_rank; ++r) {
            GTensor& t = T[static_cast<std::size_t>(r)];
            std::copy(next.begin() + static_cast<std::ptrdiff_t>(pos), next.begin() + static_cast<std::ptrdiff_t>(pos + t.total),
                      t.data.begin());
            pos += t.total;
        }
    }
    res.cc_energy = res.reference_energy + res.correlation_energy;
    res.total_energy = res.cc_energy;

    // Amplitude analysis.
    res.rank_norms.assign(static_cast<std::size_t>(max_rank) + 1, 0.0);
    struct Cand {
        double value;
        int rank;
        int block;
        std::size_t index;
    };
    std::vector<Cand> cands;
    for (int r = 1; r <= max_rank; ++r) {
        const GTensor& t = T[static_cast<std::size_t>(r)];
        double s2 = 0.0;
        for (double x : t.data) s2 += x * x;
        res.rank_norms[static_cast<std::size_t>(r)] = std::sqrt(s2);
        if (r == 1) res.t1_diagnostic = std::sqrt(s2 / (2.0 * 2.0 * o));  // spin-orbital sum = 2 x spatial
        for (std::size_t bi = 0; bi < t.blocks.size(); ++bi) {
            const GBlock& b = t.blocks[bi];
            if (2 * b.alpha[2] < r) continue;  // one of each alpha <-> beta pair
            for (std::size_t i = 0; i < b.size; ++i) {
                const double x = t.data[b.offset + i];
                if (std::fabs(x) > 1e-6) cands.push_back({x, r, static_cast<int>(bi), i});
            }
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return std::fabs(a.value) > std::fabs(b.value); });
    std::vector<std::string> seen;
    for (const Cand& c : cands) {
        if (res.leading.size() >= 12) break;
        const GTensor& t = T[static_cast<std::size_t>(c.rank)];
        const GBlock& b = t.blocks[static_cast<std::size_t>(c.block)];
        std::vector<std::pair<int, int>> parts, holes;
        decode_part(P.sp.vir, b.alpha[2], c.rank - b.alpha[2], c.index / b.count[3], o, parts);
        decode_part(P.sp.occ, b.alpha[3], c.rank - b.alpha[3], c.index % b.count[3], 0, holes);
        // Skip the spin-flipped copy of an amplitude already listed (same spatial pattern and value).
        std::string key = std::to_string(c.rank) + ":";
        std::vector<int> ks;
        for (auto& h : holes) ks.push_back(h.first);
        std::sort(ks.begin(), ks.end());
        for (int x : ks) key += std::to_string(x) + ",";
        ks.clear();
        for (auto& p : parts) ks.push_back(p.first);
        std::sort(ks.begin(), ks.end());
        for (int x : ks) key += std::to_string(x) + ",";
        char buf[32];
        std::snprintf(buf, sizeof buf, "%.8f", std::fabs(c.value));
        key += buf;
        if (std::find(seen.begin(), seen.end(), key) != seen.end()) continue;
        seen.push_back(key);
        CCAmplitude a;
        a.value = c.value;
        a.rank = c.rank;
        for (auto& h : holes) {
            a.occupied.push_back(h.first + 1);
            a.spins.push_back(h.second ? 'b' : 'a');
        }
        for (auto& p : parts) a.virtuals.push_back(p.first + 1);
        res.leading.push_back(a);
    }
    res.cc_seconds = std::chrono::duration<double>(Clock::now() - t_start).count();
    return res;
}

}  // namespace macqc
