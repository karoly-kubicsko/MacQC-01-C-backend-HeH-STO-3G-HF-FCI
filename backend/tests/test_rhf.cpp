// Tests of the ΨLab 3.0 backend for general molecules (CTest): basis-set files and
// parser, Boys function, solid harmonics, integrals against PySCF, RHF energies against
// PySCF and against published Hartree-Fock limits, and the JSON API.
//
// Cases with more basis functions than PSILAB_TEST_MAX_NBF (default 130) are skipped;
// PSILAB_TEST_MAX_NBF=400 runs the complete set (minutes).
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "macqc/basis.h"
#include "macqc/calculation.h"
#include "macqc/gto_integrals.h"
#include "macqc/macqc_api.h"
#include "macqc/molecular.h"
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

static int max_nbf() {
    const char* e = std::getenv("PSILAB_TEST_MAX_NBF");
    return e ? std::atoi(e) : 130;
}

static bool throws_with(const std::string& text, const std::string& part) {
    try {
        parse_gaussian94(text);
    } catch (const std::exception& e) {
        if (std::string(e.what()).find(part) == std::string::npos) std::printf("  message was: %s\n", e.what());
        return std::string(e.what()).find(part) != std::string::npos;
    }
    return false;
}

// ---------------------------------------------------------------- Boys function

static long double boys_series(int m, long double T) {
    long double term = 1.0L / (2 * m + 1), sum = term;
    for (int k = 1; k < 100000; ++k) {
        term *= 2.0L * T / (2 * m + 2 * k + 1);
        sum += term;
        if (term < 1e-22L * sum) break;
    }
    return std::exp(-T) * sum;
}

static void test_boys() {
    const double Ts[] = {0.0, 1e-9, 0.024, 0.3, 1.0, 2.5, 7.77, 17.3, 29.99, 44.97, 45.0, 45.03, 61.0, 150.0, 700.0};
    double worst = 0.0;
    for (double T : Ts) {
        double F[30];
        boys_function(24, T, F);
        for (int m = 0; m <= 24; ++m) {
            const long double ref = boys_series(m, T);
            worst = std::max(worst, static_cast<double>(std::fabs((F[m] - ref) / ref)));
        }
    }
    std::printf("Boys function: largest relative error %.2e\n", worst);
    CHECK(worst < 2e-13);
}

// ---------------------------------------------------------------- solid harmonics

static double sphere_monomial(int a, int b, int c) {  // integral of x^a y^b z^c over the unit sphere
    if (a % 2 || b % 2 || c % 2) return 0.0;
    return 2.0 * std::tgamma((a + 1) / 2.0) * std::tgamma((b + 1) / 2.0) * std::tgamma((c + 1) / 2.0) /
           std::tgamma((a + b + c + 3) / 2.0);
}

static void test_solid_harmonics() {
    const double pi = 3.14159265358979323846;
    double worst = 0.0;
    for (int l = 0; l <= 5; ++l) {
        const auto& M = solid_harmonic_matrix(l);
        std::vector<std::array<int, 3>> comps;
        for (int x = l; x >= 0; --x)
            for (int y = l - x; y >= 0; --y) comps.push_back({x, y, l - x - y});
        const int nc = static_cast<int>(comps.size()), ns = 2 * l + 1;
        for (int m1 = 0; m1 < ns; ++m1)
            for (int m2 = 0; m2 < ns; ++m2) {
                double s = 0.0;
                for (int i = 0; i < nc; ++i)
                    for (int j = 0; j < nc; ++j)
                        s += M[m1 * nc + i] * M[m2 * nc + j] *
                             sphere_monomial(comps[i][0] + comps[j][0], comps[i][1] + comps[j][1],
                                             comps[i][2] + comps[j][2]);
                const double want = m1 == m2 ? 4.0 * pi / (2 * l + 1) : 0.0;
                worst = std::max(worst, std::fabs(s - want));
            }
    }
    std::printf("solid harmonics l <= 5: orthonormality error %.2e\n", worst);
    CHECK(worst < 1e-12);
    // S_10 = z, S_11 = x, S_1-1 = y; S_20 = z^2 - (x^2 + y^2)/2.
    const auto& p = solid_harmonic_matrix(1);
    CHECK(p[0 * 3 + 1] == 1.0 && p[1 * 3 + 2] == 1.0 && p[2 * 3 + 0] == 1.0);
    const auto& d = solid_harmonic_matrix(2);
    check_near("S20 zz", d[2 * 6 + 5], 1.0, 1e-15);
    check_near("S20 xx", d[2 * 6 + 0], -0.5, 1e-15);
}

// ---------------------------------------------------------------- basis files

struct ElementExpectation {
    const char* basis;
    int Z;
    const char* contraction;  // "" = not defined in the file
    int functions;
};

static void test_basis_files() {
    const ElementExpectation table[] = {
        {"aug-cc-pVDZ", 1, "[3s2p]", 9},          {"aug-cc-pVDZ", 3, "[4s3p2d]", 23},
        {"aug-cc-pVDZ", 8, "[4s3p2d]", 23},       {"aug-cc-pVTZ", 1, "[4s3p2d]", 23},
        {"aug-cc-pVTZ", 3, "[5s4p3d2f]", 46},     {"aug-cc-pVTZ", 8, "[5s4p3d2f]", 46},
        {"aug-cc-pVQZ", 1, "[5s4p3d2f]", 46},     {"aug-cc-pVQZ", 3, "[6s5p4d3f2g]", 80},
        {"aug-cc-pVQZ", 8, "[6s5p4d3f2g]", 80},   {"aug-cc-pV5Z", 1, "[6s5p4d3f2g]", 80},
        {"aug-cc-pV5Z", 3, "[7s6p5d4f3g2h]", 127}, {"aug-cc-pV5Z", 8, "[7s6p5d4f3g2h]", 127},
        {"aug-cc-pCVDZ", 1, "", 0},               {"aug-cc-pCVDZ", 3, "[5s4p2d]", 27},
        {"aug-cc-pCVDZ", 8, "[5s4p2d]", 27},      {"aug-cc-pCVTZ", 1, "", 0},
        {"aug-cc-pCVTZ", 3, "[7s6p4d2f]", 59},    {"aug-cc-pCVTZ", 8, "[7s6p4d2f]", 59},
        {"aug-cc-pCVQZ", 1, "", 0},               {"aug-cc-pCVQZ", 3, "[9s8p6d4f2g]", 109},
        {"aug-cc-pCVQZ", 8, "[9s8p6d4f2g]", 109}, {"aug-cc-pCV5Z", 1, "", 0},
        {"aug-cc-pCV5Z", 3, "", 0},               {"aug-cc-pCV5Z", 8, "[11s10p8d6f4g2h]", 181},
    };
    CHECK(supported_basis_sets().size() == 8);
    for (const auto& e : table) {
        const std::string path = default_basis_directory() + "/" + basis_file_name(e.basis);
        BasisFile f;
        try {
            f = parse_gaussian94(read_file(path));
        } catch (const std::exception& ex) {
            std::printf("FAIL %s: %s\n", path.c_str(), ex.what());
            ++g_failures;
            continue;
        }
        CHECK(f.header.find("Basis Set Exchange") != std::string::npos);
        auto it = f.elements.find(e.Z);
        if (!*e.contraction) {
            CHECK(it == f.elements.end());
            continue;
        }
        CHECK(it != f.elements.end());
        if (it == f.elements.end()) continue;
        int count[8] = {0}, functions = 0;
        for (const auto& sh : it->second.shells) {
            ++count[sh.l];
            functions += 2 * sh.l + 1;
        }
        std::string con = "[";
        for (int l = 0; l < 8; ++l)
            if (count[l]) con += std::to_string(count[l]) + angular_letter(l);
        con += "]";
        if (con != e.contraction || functions != e.functions)
            std::printf("  %s Z=%d: %s %d functions\n", e.basis, e.Z, con.c_str(), functions);
        CHECK(con == e.contraction);
        CHECK(functions == e.functions);
    }
    CHECK(canonical_basis_name("AUG-CC-PVTZ") == "aug-cc-pVTZ");
    CHECK(canonical_basis_name(" aug-cc-pcv5z ") == "aug-cc-pCV5Z");
    CHECK(canonical_basis_name("6-31G*").empty());

    // Parser details: D exponents, SP shells, scale factor, element 'H' versus h shells.
    const BasisFile f = parse_gaussian94(
        "! test\n\nH     0\nS   2   1.00\n  1.0D+01  0.5\n  1.0  0.5\nH   1   1.00\n  2.0  1.0\n****\n"
        "Li 0\nSP   1   2.00\n  0.5  0.7  0.3\n****\n");
    CHECK(f.elements.size() == 2);
    CHECK(f.elements.at(1).shells.size() == 2 && f.elements.at(1).shells[1].l == 5);
    CHECK(f.elements.at(1).shells[0].exponents[0] == 10.0);
    CHECK(f.elements.at(3).shells.size() == 2 && f.elements.at(3).shells[1].l == 1);
    CHECK(f.elements.at(3).shells[0].exponents[0] == 2.0);  // 0.5 * 2.00^2
    CHECK(f.elements.at(3).shells[1].coefficients[0] == 0.3);
    CHECK(throws_with("O 0\nS 1 1.00\n 1.0 1.0\n", "without '****'"));
    CHECK(throws_with("O 0\nX 1 1.00\n 1.0 1.0\n****\n", "line 2: unknown shell type"));
    CHECK(throws_with("Xx 0\n****\n", "line 1: unknown or unsupported element"));
    CHECK(throws_with("O 0\nS 2 1.00\n 1.0 1.0\n****\n", "line 4"));
    CHECK(throws_with("O 0\nI 1 1.00\n 1.0 1.0\n****\n", "beyond the supported maximum"));
    CHECK(throws_with("O 0\nS 1 1.00\n -1.0 1.0\n****\n", "exponents must be positive"));
    CHECK(throws_with("", "no element"));

    // Contracted functions are normalised: unit diagonal of S in every basis set.
    for (const auto& name : supported_basis_sets()) {
        const Molecule mol = make_molecule("H2O", default_bond_length("H2O"), kDefaultWaterAngle);
        const BasisSet b = build_basis(mol, name, default_basis_directory());
        const Matrix S = overlap_matrix(b);
        double worst = 0.0;
        for (int i = 0; i < b.n_functions; ++i) worst = std::max(worst, std::fabs(S(i, i) - 1.0));
        if (worst > 1e-12) std::printf("  %s: max |S_ii - 1| = %.2e\n", name.c_str(), worst);
        CHECK(worst < 1e-12);
    }

    // Element coverage: H takes aug-cc-pVXZ in the core-valence sets; Li has no aug-cc-pCV5Z.
    const Molecule h2o = make_molecule("H2O", default_bond_length("H2O"), kDefaultWaterAngle);
    const BasisSet cv = build_basis(h2o, "aug-cc-pCVDZ", default_basis_directory());
    CHECK(cv.n_functions == 27 + 2 * 9);
    CHECK(cv.atoms[1].basis == "aug-cc-pVDZ" && !cv.atoms[1].note.empty());
    CHECK(cv.atoms[0].basis == "aug-cc-pCVDZ" && cv.atoms[0].note.empty());
    bool refused = false;
    try {
        build_basis(make_molecule("LiH", 3.0, 0.0), "aug-cc-pCV5Z", default_basis_directory());
    } catch (const std::exception& e) {
        refused = std::string(e.what()).find("not available for Li") != std::string::npos;
    }
    CHECK(refused);
    bool missing = false;
    try {
        build_basis(h2o, "aug-cc-pVDZ", "/nonexistent-folder");
    } catch (const std::exception& e) {
        missing = std::string(e.what()).find("cannot read the basis-set file") != std::string::npos;
    }
    CHECK(missing);
}

// ---------------------------------------------------------------- integrals against PySCF

static BasisSet custom_basis(const json& c, Molecule& mol) {
    for (const auto& a : c["atoms"]) {
        Atom at;
        at.symbol = a[0].get<std::string>();
        at.Z = atomic_number(at.symbol);
        at.position = {a[1][0].get<double>(), a[1][1].get<double>(), a[1][2].get<double>()};
        mol.atoms.push_back(at);
    }
    BasisSet b;
    for (const auto& s : c["shells"]) {
        Shell sh;
        sh.atom = s[0].get<int>();
        sh.l = s[1].get<int>();
        sh.center = mol.atoms[static_cast<std::size_t>(sh.atom)].position;
        std::vector<double> e, k;
        for (const auto& p : s[2]) {
            e.push_back(p[0].get<double>());
            k.push_back(p[1].get<double>());
        }
        sh.exponents = e;
        sh.coefficients = normalized_coefficients(sh.l, e, k);
        sh.first_function = b.n_functions;
        b.n_functions += sh.size();
        b.max_l = std::max(b.max_l, sh.l);
        b.shells.push_back(sh);
    }
    return b;
}

static void test_integrals() {
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/integral_reference.json"));
    for (const auto& c : ref["cases"]) {
        const std::string name = c["name"].get<std::string>();
        Molecule mol;
        BasisSet b;
        if (name == "test-basis-l5") {
            b = custom_basis(c, mol);
        } else {
            mol = make_molecule("H2O", default_bond_length("H2O"), kDefaultWaterAngle);
            b = build_basis(mol, name.substr(4), default_basis_directory());
        }
        CHECK(b.n_functions == c["nbf"].get<int>());
        const std::size_t n = static_cast<std::size_t>(b.n_functions);
        auto trace_norm = [&](const Matrix& M, const char* what) {
            double tr = 0.0, nn = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                tr += M(i, i);
                for (std::size_t j = 0; j < n; ++j) nn += M(i, j) * M(i, j);
            }
            check_near(name + " trace " + what, tr, c[std::string("trace_") + what].get<double>(), 1e-10);
            check_near(name + " norm " + what, std::sqrt(nn), c[std::string("norm_") + what].get<double>(), 1e-10);
        };
        trace_norm(overlap_matrix(b), "S");
        trace_norm(kinetic_matrix(b), "T");
        trace_norm(nuclear_attraction_matrix(b, mol), "V");

        // Full ERI tensor from the unique shell quartets.
        std::vector<double> eri(n * n * n * n, 0.0);
        const ShellPairList pairs = make_shell_pairs(b);
        EriEngine engine(b);
        for (std::size_t P = 0; P < pairs.pairs.size(); ++P)
            for (std::size_t Q = 0; Q <= P; ++Q) {
                const ShellPair& bra = pairs.pairs[P];
                const ShellPair& ket = pairs.pairs[Q];
                const double* v = engine.compute(bra, ket);
                const Shell& A = b.shells[bra.a];
                const Shell& B = b.shells[bra.b];
                const Shell& C = b.shells[ket.a];
                const Shell& D = b.shells[ket.b];
                for (int i = 0; i < A.size(); ++i)
                    for (int j = 0; j < B.size(); ++j)
                        for (int k = 0; k < C.size(); ++k)
                            for (int l = 0; l < D.size(); ++l) {
                                const double x = *v++;
                                const std::size_t p = A.first_function + i, q = B.first_function + j,
                                                  r = C.first_function + k, s = D.first_function + l;
                                const std::size_t idx[8][4] = {{p, q, r, s}, {q, p, r, s}, {p, q, s, r}, {q, p, s, r},
                                                               {r, s, p, q}, {s, r, p, q}, {r, s, q, p}, {s, r, q, p}};
                                for (const auto& t : idx) eri[((t[0] * n + t[1]) * n + t[2]) * n + t[3]] = x;
                            }
            }
        double norm = 0.0, coulomb = 0.0, exchange = 0.0;
        for (double x : eri) norm += x * x;
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) {
                coulomb += eri[((p * n + p) * n + q) * n + q];
                exchange += eri[((p * n + q) * n + p) * n + q];
            }
        check_near(name + " ERI norm", std::sqrt(norm), c["norm_ERI"].get<double>(), 1e-9);
        check_near(name + " sum (pp|qq)", coulomb, c["coulomb_sum"].get<double>(), 1e-9);
        check_near(name + " sum (pq|pq)", exchange, c["exchange_sum"].get<double>(), 1e-9);
    }
}

// ---------------------------------------------------------------- RHF energies

static RHFResult rhf(const std::string& molecule, const std::string& basis, double R, double angle,
                     double memory_mb = 1024.0) {
    const Molecule mol = make_molecule(molecule, R, angle);
    const BasisSet b = build_basis(mol, basis, default_basis_directory());
    RHFOptions opt;
    opt.memory_limit_mb = memory_mb;
    return run_rhf(mol, b, opt);
}

static void test_rhf_reference() {
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/rhf_reference.json"));
    const int limit = max_nbf();
    int run = 0, skipped = 0;
    for (const auto& c : ref["cases"]) {
        const std::string m = c["molecule"].get<std::string>(), b = c["basis"].get<std::string>();
        if (c["nbf"].get<int>() > limit) {
            ++skipped;
            continue;
        }
        const double R = c["R_bohr"].get<double>();
        const double angle = c["angle"].is_null() ? kDefaultWaterAngle : c["angle"].get<double>();
        const RHFResult r = rhf(m, b, R, angle);
        ++run;
        CHECK(r.converged);
        CHECK(r.n_basis == c["nbf"].get<int>());
        check_near(m + "/" + b + " E_RHF", r.total_energy, c["energy"].get<double>(), 1e-8);
        check_near(m + "/" + b + " E_nuc", r.nuclear_repulsion, c["nuclear_repulsion"].get<double>(), 1e-10);
        check_near(m + "/" + b + " HOMO", r.orbital_energies[r.n_occupied - 1], c["homo"].get<double>(), 1e-6);
        check_near(m + "/" + b + " LUMO", r.orbital_energies[r.n_occupied], c["lumo"].get<double>(), 1e-6);
        std::printf("  %-4s %-13s %4d functions  E = %.10f  (PySCF %.10f)  %d it, %.1f s\n", m.c_str(), b.c_str(),
                    r.n_basis, r.total_energy, c["energy"].get<double>(), r.iterations,
                    r.integral_seconds + r.scf_seconds);
        std::fflush(stdout);
    }
    std::printf("RHF against PySCF: %d cases run, %d larger cases skipped (PSILAB_TEST_MAX_NBF)\n", run, skipped);
    CHECK(run >= 10);

    // Integral-direct Fock builds give the same energy as stored integrals.
    const double R = default_bond_length("H2O");
    const RHFResult stored = rhf("H2O", "aug-cc-pVDZ", R, kDefaultWaterAngle);
    const RHFResult direct = rhf("H2O", "aug-cc-pVDZ", R, kDefaultWaterAngle, 0.0);
    CHECK(stored.eri_stored && !direct.eri_stored);
    check_near("direct vs stored", direct.total_energy, stored.total_energy, 1e-10);
    // Virial theorem near equilibrium and the energy decomposition.
    check_near("virial ratio", stored.virial_ratio, 2.0, 5e-3);
    check_near("energy components",
               stored.kinetic_energy + stored.nuclear_attraction_energy + stored.two_electron_energy +
                   stored.nuclear_repulsion,
               stored.total_energy, 1e-9);
}

// Published Hartree-Fock limits (numerical HF for diatomics): the basis-set energies
// approach them from above.
static void test_hf_limits() {
    struct Limit {
        const char* molecule;
        double R;      // bohr
        double limit;  // E_h
        const char* source;
    };
    const Limit limits[] = {
        {"H2", 1.4, -1.133629571, "numerical HF; Kolos & Roothaan 1960: -1.13363"},
        {"LiH", 3.015, -7.987352, "numerical HF, Laaksonen, Pyykko & Sundholm 1986"},
    };
    for (const Limit& L : limits) {
        double previous = 0.0;
        const char* bases[] = {"aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z"};
        const double bound[] = {1.0e-2, 8.0e-4, 2.5e-4, 4.0e-5};
        for (int k = 0; k < 4; ++k) {
            const Molecule mol = make_molecule(L.molecule, L.R, 0.0);
            const BasisSet b = build_basis(mol, bases[k], default_basis_directory());
            if (b.n_functions > max_nbf()) continue;
            const RHFResult r = run_rhf(mol, b, RHFOptions());
            const double above = r.total_energy - L.limit;
            std::printf("  %-3s R = %.3f bohr %-12s E = %.10f, %.2e E_h above the HF limit %.9f (%s)\n", L.molecule,
                        L.R, bases[k], r.total_energy, above, L.limit, L.source);
            CHECK(r.converged);
            CHECK(above > 0.0 && above < bound[k]);
            if (k > 0) CHECK(r.total_energy < previous);
            previous = r.total_energy;
        }
    }
}

// ---------------------------------------------------------------- API

static json run_json(const std::string& text) {
    char* p = macqc_run(text.c_str());
    CHECK(p != nullptr);
    if (!p) return json();
    json j = json::parse(p);
    macqc_free(p);
    return j;
}

static void test_api() {
    const json ok = run_json(R"({"molecule":"H2O","basis":"aug-cc-pVDZ","printLevel":1})");
    CHECK(ok["error"].is_null());
    CHECK(ok["hf"].is_null() && ok["fci"].is_null());
    check_near("api H2O", ok["rhf"]["totalEnergy"].get<double>(), -76.0414279605, 1e-8);
    CHECK(ok["rhf"]["converged"].get<bool>());
    CHECK(ok["rhf"]["occupiedOrbitals"].get<int>() == 5);
    CHECK(ok["basis"]["functions"].get<int>() == 41);
    CHECK(ok["molecule"]["atoms"].size() == 3);
    CHECK(ok["input"]["molecule"].get<std::string>() == "H2O");
    check_near("api default R", ok["input"]["R"].get<double>(), 0.9572 / kBohrInAngstrom, 1e-12);
    CHECK(ok["report"].get<std::string>().find("restricted Hartree-Fock") != std::string::npos);
    CHECK(!ok["estimate"]["longCalculation"].get<bool>());

    // Estimate only: no SCF, warnings for a long calculation.
    const json est = run_json(R"({"molecule":"H2O","basis":"aug-cc-pCV5Z","estimateOnly":true,"threads":4})");
    CHECK(est["error"].is_null());
    CHECK(est["rhf"].is_null());
    CHECK(est["basis"]["functions"].get<int>() == 341);
    CHECK(est["basis"]["maxL"].get<int>() == 5);
    CHECK(!est["estimate"]["integralsStored"].get<bool>());
    CHECK(est["estimate"]["longCalculation"].get<bool>());
    CHECK(est["warnings"].size() == 1);
    CHECK(est["notes"].size() == 1);  // aug-cc-pV5Z for H

    // Refusals and invalid input are errors with explanations, never crashes.
    const json lih = run_json(R"({"molecule":"LiH","basis":"aug-cc-pCV5Z"})");
    CHECK(lih["error"].is_string() && lih["error"].get<std::string>().find("not available for Li") != std::string::npos);
    CHECK(lih["rhf"].is_null());
    const char* bad[] = {R"({"molecule":"CH4"})",
                         R"({"molecule":"H2O","basis":"6-31G"})",
                         R"({"molecule":"H2O","angle":200})",
                         R"({"molecule":"H2","R":-1})",
                         R"({"molecule":"H2","R":0.1})",
                         R"({"molecule":"H2","method":"CCSD"})",
                         R"({"molecule":"H2","threads":-2})",
                         R"({"molecule":"H2","basis":7})",
                         R"({"molecule":"H2","basisDirectory":"/nonexistent"})"};
    for (const char* b : bad) {
        const json j = run_json(b);
        const bool good = j["error"].is_string();
        if (!good) std::printf("invalid input not rejected: %s\n", b);
        CHECK(good);
    }
    // HeH+ (the default) keeps its original output format.
    const json heh = run_json("{}");
    CHECK(heh["error"].is_null() && heh["hf"].is_object() && heh["fci"].is_object());
    CHECK(!heh.contains("rhf"));
    CHECK(!heh["input"].contains("molecule"));
}

// Job handles: progress while running, cancellation from another thread.
static void test_jobs() {
    // Cancel an integral-direct RHF once its SCF iterations have started.
    macqc_job* job = macqc_job_create();
    CHECK(job != nullptr);
    char* raw = nullptr;
    std::thread worker([&] {
        raw = macqc_run_job(R"({"molecule":"H2O","basis":"aug-cc-pVTZ","memoryMB":0,"printLevel":0})", job);
    });
    int stage = 0, iteration = 0;
    double energy = 0.0, residual = 0.0;
    for (int k = 0; k < 20000 && stage != 2; ++k) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        macqc_job_progress(job, &stage, &iteration, &energy, &residual);
    }
    CHECK(stage == 2);  // SCF iterations running
    macqc_job_cancel(job);
    worker.join();
    const json cancelled = json::parse(raw);
    macqc_free(raw);
    CHECK(cancelled["error"].is_string() && cancelled["error"].get<std::string>() == "calculation cancelled");
    CHECK(cancelled["rhf"].is_null());
    macqc_job_free(job);

    // A Full CI job reports its stages and ends in stage 5 (done).
    job = macqc_job_create();
    int seen = 0;
    raw = nullptr;
    std::thread fci([&] {
        raw = macqc_run_job(R"({"molecule":"LiH","basis":"aug-cc-pVDZ","method":"FCI","printLevel":0})", job);
    });
    for (int k = 0; k < 60000; ++k) {
        macqc_job_progress(job, &stage, &iteration, &energy, &residual);
        if (stage == 4 && iteration > 0) seen = 1;
        if (stage == 5) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    fci.join();
    const json done = json::parse(raw);
    macqc_free(raw);
    macqc_job_progress(job, &stage, &iteration, &energy, &residual);
    CHECK(done["error"].is_null() && stage == 5);
    CHECK(seen == 1);
    macqc_job_free(job);
    // macqc_run without a job and NULL handles stay safe.
    macqc_job_cancel(nullptr);
    macqc_job_progress(nullptr, &stage, nullptr, nullptr, nullptr);
    macqc_job_free(nullptr);
}

int main() {
    try {
        test_boys();
        test_solid_harmonics();
        test_basis_files();
        test_integrals();
        test_rhf_reference();
        test_hf_limits();
        test_api();
        test_jobs();
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
