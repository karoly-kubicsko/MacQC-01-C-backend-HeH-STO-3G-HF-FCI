// Tests of the general determinant Full CI (CTest):
//  - occupation-string addressing;
//  - the string-driven sigma algorithm against a brute-force Hamiltonian matrix built with
//    spin-orbital Slater-Condon rules (1, 2 and 5 electrons per spin, truncated orbital spaces);
//  - Full CI energies against PySCF (backend/tests/data/fci_reference.json) and against
//    exact non-relativistic energies from the literature (variational bound, monotone);
//  - size estimates, refusals and the JSON API.
// PSILAB_TEST_MAX_DET (default 500000) limits the PySCF comparison; larger cases take minutes.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/ci.h"
#include "macqc/macqc_api.h"
#include "macqc/molecular.h"
#include "macqc/mo_integrals.h"
#include "macqc/rhf.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

#define CHECK(cond)                                                             \
    do {                                                                        \
        ++g_checks;                                                             \
        if (!(cond)) {                                                          \
            ++g_failures;                                                       \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);         \
        }                                                                       \
    } while (0)

static void check_near(const std::string& name, double got, double want, double tol) {
    ++g_checks;
    if (!(std::fabs(got - want) <= tol)) {
        ++g_failures;
        std::printf("FAIL %s: got %.12f, expected %.12f (|diff| %.2e > %.1e)\n", name.c_str(), got, want,
                    std::fabs(got - want), tol);
    }
}

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static double max_det() {
    const char* e = std::getenv("PSILAB_TEST_MAX_DET");
    return e ? std::atof(e) : 5.0e5;
}

struct Prepared {
    Molecule mol;
    BasisSet basis;
    RHFResult rhf;
    MOIntegrals mo;
};

static Prepared prepare(const std::string& molecule, const std::string& basis, double R) {
    Prepared p;
    p.mol = make_molecule(molecule, R, kDefaultWaterAngle);
    p.basis = build_basis(p.mol, basis, default_basis_directory());
    p.rhf = run_rhf(p.mol, p.basis, RHFOptions());
    p.mo = transform_to_mo(p.basis, p.rhf.H, p.rhf.C, 0);
    return p;
}

// The first n_act orbitals of `mo` (all electrons kept: a truncated orbital space).
static MOIntegrals truncate(const MOIntegrals& mo, int n_act) {
    MOIntegrals t;
    t.n = n_act;
    t.npair = static_cast<std::size_t>(n_act) * (n_act + 1) / 2;
    t.stride = t.npair;
    t.h = Matrix(static_cast<std::size_t>(n_act), static_cast<std::size_t>(n_act));
    t.eri.assign(t.npair * t.npair, 0.0);
    for (int p = 0; p < n_act; ++p)
        for (int q = 0; q < n_act; ++q) {
            t.h(static_cast<std::size_t>(p), static_cast<std::size_t>(q)) =
                mo.h(static_cast<std::size_t>(p), static_cast<std::size_t>(q));
            for (int r = 0; r < n_act; ++r)
                for (int s = 0; s < n_act; ++s)
                    t.eri[MOIntegrals::pair(p, q) * t.stride + MOIntegrals::pair(r, s)] = mo(p, q, r, s);
        }
    return t;
}

// ---------------------------------------------------------------- brute force

// Spin orbitals: alpha p -> p, beta p -> n + p; a determinant is a bit mask, the creators
// in ascending spin-orbital order (alpha string before beta string, as in the CI code).
static int popcount(unsigned long long x) {
    int c = 0;
    for (; x; x &= x - 1) ++c;
    return c;
}

static bool excite(unsigned long long det, int a, int i, unsigned long long& out, int& sign) {
    if (!(det >> i & 1ULL)) return false;
    unsigned long long m = det & ~(1ULL << i);
    if (m >> a & 1ULL) return false;
    const int lo = std::min(i, a), hi = std::max(i, a);
    int between = 0;
    for (int k = lo + 1; k < hi; ++k) between += static_cast<int>(m >> k & 1ULL);
    sign = (between % 2) ? -1 : 1;
    out = m | (1ULL << a);
    return true;
}

static double brute_force_ground_state(const MOIntegrals& mo, int na, int nb) {
    const int n = mo.n, nso = 2 * n;
    auto h1 = [&](int s, int t) { return (s < n) == (t < n) ? mo.h(static_cast<std::size_t>(s % n), static_cast<std::size_t>(t % n)) : 0.0; };
    // <st|uv> physicists'
    auto g = [&](int s, int t, int u, int v) {
        if ((s < n) != (u < n) || (t < n) != (v < n)) return 0.0;
        return mo(s % n, u % n, t % n, v % n);
    };
    auto anti = [&](int s, int t, int u, int v) { return g(s, t, u, v) - g(s, t, v, u); };
    std::vector<unsigned long long> dets;
    const StringSpace A(n, na), B(n, nb);
    for (std::size_t I = 0; I < A.size(); ++I)
        for (std::size_t J = 0; J < B.size(); ++J) {
            unsigned long long d = 0;
            for (int k = 0; k < na; ++k) d |= 1ULL << A.occupation(I)[k];
            for (int k = 0; k < nb; ++k) d |= 1ULL << (n + B.occupation(J)[k]);
            dets.push_back(d);
        }
    const std::size_t nd = dets.size();
    Matrix H(nd, nd);
    for (std::size_t x = 0; x < nd; ++x)
        for (std::size_t y = 0; y < nd; ++y) {
            const unsigned long long bra = dets[x], ket = dets[y], diff = bra ^ ket;
            const int nd2 = popcount(diff);
            double v = 0.0;
            if (nd2 == 0) {
                std::vector<int> occ;
                for (int k = 0; k < nso; ++k)
                    if (ket >> k & 1ULL) occ.push_back(k);
                for (int i : occ) v += h1(i, i);
                for (int i : occ)
                    for (int j : occ) v += 0.5 * anti(i, j, i, j);
            } else if (nd2 == 2) {
                int i = 0, a = 0;
                while (!((ket & diff) >> i & 1ULL)) ++i;
                while (!((bra & diff) >> a & 1ULL)) ++a;
                unsigned long long r = 0;
                int sg = 1;
                excite(ket, a, i, r, sg);
                v = h1(a, i);
                for (int k = 0; k < nso; ++k)
                    if ((bra & ket) >> k & 1ULL) v += anti(a, k, i, k);
                v *= sg;
            } else if (nd2 == 4) {
                std::vector<int> from, to;
                for (int k = 0; k < nso; ++k) {
                    if ((ket & diff) >> k & 1ULL) from.push_back(k);
                    if ((bra & diff) >> k & 1ULL) to.push_back(k);
                }
                unsigned long long d1 = 0, d2 = 0;
                int s1 = 1, s2 = 1;
                excite(ket, to[0], from[0], d1, s1);
                excite(d1, to[1], from[1], d2, s2);
                v = s1 * s2 * anti(to[0], to[1], from[0], from[1]);
            }
            H(x, y) = v;
        }
    return symmetric_eigen(H).values[0];
}

// ---------------------------------------------------------------- tests

static void test_strings() {
    const StringSpace S(9, 4);
    CHECK(S.size() == 126);
    bool ok = true;
    for (std::size_t I = 0; I < S.size(); ++I) {
        const std::uint16_t* o = S.occupation(I);
        for (int k = 1; k < 4; ++k) ok = ok && o[k - 1] < o[k];
        ok = ok && S.index(o) == I;
    }
    CHECK(ok);
    CHECK(S.occupation(0)[0] == 0 && S.occupation(0)[3] == 3);  // lowest string first
    const StringSpace E(5, 0);
    CHECK(E.size() == 1);
    CHECK(binomial_count(41, 5) == 749398.0);
    CHECK(binomial_count(160, 1) == 160.0);
}

static void test_brute_force() {
    // H2 (1 alpha, 1 beta) in all 18 orbitals; LiH (2, 2) in 8 and H2O (5, 5) in 7 orbitals.
    struct Case {
        const char* molecule;
        const char* basis;
        int n_act;
    };
    const Case cases[] = {{"H2", "aug-cc-pVDZ", 18}, {"LiH", "aug-cc-pVDZ", 8}, {"H2O", "aug-cc-pVDZ", 7}};
    for (const Case& c : cases) {
        const Prepared p = prepare(c.molecule, c.basis, default_bond_length(c.molecule));
        const MOIntegrals t = truncate(p.mo, c.n_act);
        const int nocc = p.rhf.n_occupied;
        const double e_bf = brute_force_ground_state(t, nocc, nocc) + p.rhf.nuclear_repulsion;
        FCIOptions o;
        o.residual_tol = 1e-8;
        o.energy_tol = 1e-12;
        const FCIResultGeneral r = run_fci_determinants(t, nocc, nocc, p.rhf.nuclear_repulsion, p.rhf.total_energy, o);
        std::printf("  %-3s %d electrons in %2d orbitals (%zu determinants): sigma algorithm %.12f, dense matrix %.12f\n",
                    c.molecule, 2 * nocc, c.n_act, r.determinants, r.total_energy, e_bf);
        CHECK(r.converged);
        check_near(std::string(c.molecule) + " sigma vs dense", r.total_energy, e_bf, 1e-10);
        check_near(std::string(c.molecule) + " <S^2>", r.s2, 0.0, 1e-8);
        double sum = 0.0;
        for (double x : r.natural_occupations) sum += x;
        check_near(std::string(c.molecule) + " occupation sum", sum, 2.0 * nocc, 1e-10);
        CHECK(r.correlation_energy < 0.0);
    }
}

static void test_reference() {
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/fci_reference.json"));
    int run = 0, skipped = 0;
    for (const auto& c : ref["cases"]) {
        const std::string m = c["molecule"].get<std::string>(), b = c["basis"].get<std::string>();
        if (c["determinants"].get<double>() > max_det()) {
            ++skipped;
            continue;
        }
        const Prepared p = prepare(m, b, c["R_bohr"].get<double>());
        const int nocc = p.rhf.n_occupied;
        const FCIResultGeneral r =
            run_fci_determinants(p.mo, nocc, nocc, p.rhf.nuclear_repulsion, p.rhf.total_energy, FCIOptions());
        ++run;
        std::printf("  %-4s %-13s %10zu determinants  E_FCI = %.10f (PySCF %.10f)  %d it, %.1f s\n", m.c_str(),
                    b.c_str(), r.determinants, r.total_energy, c["fci_energy"].get<double>(), r.iterations,
                    r.ci_seconds);
        std::fflush(stdout);
        CHECK(r.converged);
        CHECK(static_cast<double>(r.determinants) == c["determinants"].get<double>());
        check_near(m + "/" + b + " E_FCI", r.total_energy, c["fci_energy"].get<double>(), 2e-9);
        check_near(m + "/" + b + " E_RHF", p.rhf.total_energy, c["rhf_energy"].get<double>(), 1e-9);
        check_near(m + "/" + b + " |C0|^2", r.hf_weight, c["hf_weight"].get<double>(), 1e-6);
        check_near(m + "/" + b + " <S^2>", r.s2, 0.0, 1e-8);
    }
    std::printf("Full CI against PySCF: %d cases run, %d larger cases skipped (PSILAB_TEST_MAX_DET)\n", run, skipped);
    CHECK(run >= 4);
}

// Exact non-relativistic (Born-Oppenheimer) energies: Full CI in a finite basis lies above
// them and approaches them as the basis grows.
static void test_exact_bounds() {
    struct Exact {
        const char* molecule;
        double R, energy;
        const char* source;
        const char* bases[4];
        bool monotone;  // a hierarchy of nested basis sets
    };
    const Exact cases[] = {
        {"H2", 1.4, -1.1744757142, "Pachucki, Phys. Rev. A 82, 032509 (2010)",
         {"aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", nullptr}, true},
        {"LiH", 3.015, -8.070548, "explicitly correlated Gaussians, Tung, Pavanello & Adamowicz 2011",
         {"aug-cc-pVDZ", "aug-cc-pCVDZ", nullptr, nullptr}, false},
    };
    for (const Exact& e : cases) {
        double previous = 0.0;
        for (int k = 0; k < 4 && e.bases[k]; ++k) {
            const Prepared p = prepare(e.molecule, e.bases[k], e.R);
            const int nocc = p.rhf.n_occupied;
            if (binomial_count(p.mo.n, nocc) * binomial_count(p.mo.n, nocc) > max_det()) continue;
            const FCIResultGeneral r =
                run_fci_determinants(p.mo, nocc, nocc, p.rhf.nuclear_repulsion, p.rhf.total_energy, FCIOptions());
            std::printf("  %-3s R = %.3f bohr %-13s E_FCI = %.10f, %.2e E_h above the exact %.10f (%s)\n", e.molecule,
                        e.R, e.bases[k], r.total_energy, r.total_energy - e.energy, e.energy, e.source);
            CHECK(r.total_energy > e.energy);
            if (k > 0 && e.monotone) CHECK(r.total_energy < previous);
            previous = r.total_energy;
        }
    }
}

static json run_json(const std::string& text) {
    char* p = macqc_run(text.c_str());
    CHECK(p != nullptr);
    if (!p) return json();
    json j = json::parse(p);
    macqc_free(p);
    return j;
}

static void test_limits_and_api() {
    // Sizes: H2O is impossible in every basis set, LiH from aug-cc-pVQZ on.
    const FCISize h2o = estimate_fci(41, 41, 5, 5, 4, 4096, 1.0);
    CHECK(!h2o.feasible && h2o.determinants > 5.6e11 && h2o.reason.find("not possible") != std::string::npos);
    CHECK(!estimate_fci(126, 126, 2, 2, 4, 4096, 1.0).feasible);  // LiH aug-cc-pVQZ
    CHECK(!estimate_fci(155, 155, 2, 2, 4, 4096, 1.0).feasible);  // LiH aug-cc-pCVQZ
    CHECK(estimate_fci(82, 82, 2, 2, 4, 4096, 1.0).feasible);     // LiH aug-cc-pCVTZ
    CHECK(estimate_fci(160, 160, 1, 1, 4, 4096, 1.0).feasible);   // H2 aug-cc-pV5Z
    CHECK(!estimate_fci(82, 82, 2, 2, 4, 512, 1.0).feasible);     // memory limit

    const json ok = run_json(R"({"molecule":"LiH","basis":"aug-cc-pVDZ","method":"FCI","printLevel":1})");
    CHECK(ok["error"].is_null() && ok["refusal"].is_null());
    CHECK(ok["fullCI"]["converged"].get<bool>());
    CHECK(ok["fullCI"]["determinants"].get<double>() == 246016.0);
    CHECK(ok["fullCI"]["leadingDeterminants"].size() >= 3);
    CHECK(ok["fullCI"]["naturalOccupations"].size() == 32);
    CHECK(ok["estimate"]["fci"]["feasible"].get<bool>());
    CHECK(ok["report"].get<std::string>().find("Full CI energies") != std::string::npos);
    check_near("api FCI = RHF + correlation", ok["fullCI"]["totalEnergy"].get<double>(),
               ok["rhf"]["totalEnergy"].get<double>() + ok["fullCI"]["correlationEnergy"].get<double>(), 1e-12);

    const json water = run_json(R"({"molecule":"H2O","method":"FCI"})");
    CHECK(water["error"].is_string() && water["refusal"].is_string());
    CHECK(water["error"].get<std::string>().find("impossible") != std::string::npos);
    CHECK(water["rhf"].is_null() && water["fullCI"].is_null());

    const json est = run_json(R"({"molecule":"LiH","basis":"aug-cc-pV5Z","method":"FCI","estimateOnly":true})");
    CHECK(est["error"].is_null());  // an estimate is not an error
    CHECK(est["refusal"].is_string() && !est["estimate"]["fci"]["feasible"].get<bool>());

    const json big = run_json(R"({"molecule":"LiH","basis":"aug-cc-pCVTZ","method":"FCI","estimateOnly":true,"threads":4})");
    CHECK(big["refusal"].is_null() && big["estimate"]["longCalculation"].get<bool>() && big["warnings"].size() == 1);
    CHECK(run_json(R"({"molecule":"H2","method":"MP2"})")["error"].is_string());
}

int main() {
    try {
        test_strings();
        test_brute_force();
        test_reference();
        test_exact_bounds();
        test_limits_and_api();
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
