// Tests of the analytic derivatives of ΨLab 5.0 (CTest):
//  - first and second derivatives of the overlap, kinetic, nuclear-attraction, dipole and
//    electron-repulsion integrals against finite differences of the integrals, for a test
//    basis with every angular momentum up to h (l = 5) on three different atoms;
//  - the RHF gradient against finite differences of the energy (also with removed
//    near-linear dependencies) and against PySCF;
//  - the RHF Hessian and the dipole derivatives against finite differences of the analytic
//    gradient and dipole, against PySCF, and the translational sum rules.
//
// Finite differences use the 4-point central formula
//   f'(0) = [8 (f(h) - f(-h)) - (f(2h) - f(-2h))] / (12 h)   (error O(h^4)).
// Cases with more basis functions than PSILAB_TEST_MAX_NBF (default 130) are skipped; the
// CI job "full" runs all of them.
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/rhf.h"
#include "macqc/rhf_derivatives.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

static void check_near(const std::string& name, double got, double want, double tol, bool verbose = false) {
    ++g_checks;
    const bool ok = std::fabs(got - want) <= tol;
    if (!ok) ++g_failures;
    if (!ok || verbose)
        std::printf("%s %-58s got %.12f expected %.12f (diff %.1e, tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), got,
                    want, got - want, tol);
}

static void check_max(const std::string& name, double worst, double tol) {
    ++g_checks;
    const bool ok = worst <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-58s max |diff| %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
}

#define CHECK(cond)                                                     \
    do {                                                                \
        ++g_checks;                                                     \
        if (!(cond)) {                                                  \
            ++g_failures;                                               \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                               \
    } while (0)

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static int max_nbf() {
    const char* e = std::getenv("PSILAB_TEST_MAX_NBF");
    return e ? std::atoi(e) : 130;
}

static bool full_mode() { return max_nbf() > 200; }

using Clock = std::chrono::steady_clock;
static double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

template <class F>
static double fd4(F f, double h) {
    return (8.0 * (f(h) - f(-h)) - (f(2.0 * h) - f(-2.0 * h))) / (12.0 * h);
}

static Molecule displaced(const Molecule& mol, int x, double h) {
    Molecule m = mol;
    m.atoms[static_cast<std::size_t>(x / 3)].position[static_cast<std::size_t>(x % 3)] += h;
    return m;
}

static double dotm(const Matrix& a, const Matrix& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) s += a(i, j) * b(i, j);
    return s;
}

// The test basis of integral_reference.json: O, H, Li with shells s ... h.
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

static Matrix random_symmetric(std::size_t n, std::mt19937& rng) {
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    Matrix m(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) m(i, j) = m(j, i) = u(rng);
    return m;
}

// ---------------------------------------------------------------- derivative integrals

static void test_derivative_integrals() {
    std::printf("\n== Derivative integrals (test basis s ... h on O, H, Li) against finite differences\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/integral_reference.json"));
    Molecule mol;
    BasisSet b;
    for (const auto& c : ref["cases"])
        if (c["name"] == "test-basis-l5") b = custom_basis(c, mol);
    CHECK(b.max_l == 5);
    const int nat = static_cast<int>(mol.atoms.size());
    const std::size_t N = static_cast<std::size_t>(b.n_functions);
    std::mt19937 rng(20261113);
    const Matrix D = random_symmetric(N, rng);
    TwoParticleDensity gamma;
    gamma.separable.push_back({random_symmetric(N, rng), random_symmetric(N, rng), 0.7, -0.3});
    const double h = 1.0e-3;
    auto at = [&](int x, double d) {
        const Molecule m = displaced(mol, x, d);
        return std::make_pair(m, basis_at_geometry(b, m));
    };

    // ---- first derivatives, contracted
    auto t0 = Clock::now();
    const std::vector<double> gS = overlap_gradient(b, nat, D);
    const std::vector<double> gH = core_hamiltonian_gradient(b, mol, D);
    const std::vector<double> gE = eri_gradient(b, nat, gamma);
    std::printf("     analytic first derivatives: %.2f s\n", since(t0));
    double wS = 0, wH = 0, wE = 0;
    for (int x = 0; x < 3 * nat; ++x) {
        const double fS = fd4([&](double d) { auto g = at(x, d); return dotm(D, overlap_matrix(g.second)); }, h);
        const double fH = fd4(
            [&](double d) {
                auto g = at(x, d);
                return dotm(D, kinetic_matrix(g.second)) + dotm(D, nuclear_attraction_matrix(g.second, g.first));
            },
            h);
        const double fE = fd4([&](double d) { return eri_contract(at(x, d).second, gamma); }, h);
        wS = std::max(wS, std::fabs(gS[static_cast<std::size_t>(x)] - fS) / std::max(1.0, std::fabs(fS)));
        wH = std::max(wH, std::fabs(gH[static_cast<std::size_t>(x)] - fH) / std::max(1.0, std::fabs(fH)));
        wE = std::max(wE, std::fabs(gE[static_cast<std::size_t>(x)] - fE) / std::max(1.0, std::fabs(fE)));
    }
    check_max("sum D dS/dx (relative)", wS, 1e-9);
    check_max("sum D dh/dx, h = T + V incl. operator (relative)", wH, 1e-9);
    check_max("sum Gamma d(mn|ls)/dx (relative)", wE, 1e-9);

    // ---- first derivatives, matrices
    {
        const std::vector<Matrix> Sx = overlap_derivative_matrices(b, nat);
        const std::vector<Matrix> Hx = core_hamiltonian_derivative_matrices(b, mol);
        const std::vector<Matrix> Gx = eri_derivative_fock(b, nat, D);
        const std::array<double, 3> origin{0.3, -0.2, 0.5};
        const std::vector<std::vector<Matrix>> Rx = dipole_derivative_matrices(b, nat, origin);
        double ws = 0, wh = 0, wg = 0, wr = 0;
        for (int x = 0; x < 3 * nat; ++x) {
            std::vector<Matrix> S(4), Hm(4), G(4);
            std::vector<std::vector<Matrix>> R(4);
            const double steps[4] = {h, -h, 2 * h, -2 * h};
            for (int k = 0; k < 4; ++k) {
                auto g = at(x, steps[k]);
                S[static_cast<std::size_t>(k)] = overlap_matrix(g.second);
                Hm[static_cast<std::size_t>(k)] = kinetic_matrix(g.second);
                const Matrix V = nuclear_attraction_matrix(g.second, g.first);
                for (std::size_t i = 0; i < N; ++i)
                    for (std::size_t j = 0; j < N; ++j) Hm[static_cast<std::size_t>(k)](i, j) += V(i, j);
                FockBuilder fb(g.second, 0, 1024.0);
                G[static_cast<std::size_t>(k)] = fb.build(D);
                R[static_cast<std::size_t>(k)] = dipole_matrices(g.second, origin);
            }
            auto fdm = [&](const std::vector<Matrix>& M, std::size_t i, std::size_t j) {
                return (8.0 * (M[0](i, j) - M[1](i, j)) - (M[2](i, j) - M[3](i, j))) / (12.0 * h);
            };
            const std::size_t xs = static_cast<std::size_t>(x);
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j) {
                    ws = std::max(ws, std::fabs(Sx[xs](i, j) - fdm(S, i, j)));
                    wh = std::max(wh, std::fabs(Hx[xs](i, j) - fdm(Hm, i, j)));
                    wg = std::max(wg, std::fabs(Gx[xs](i, j) - fdm(G, i, j)));
                    for (std::size_t m = 0; m < 3; ++m) {
                        std::vector<Matrix> Rm = {R[0][m], R[1][m], R[2][m], R[3][m]};
                        wr = std::max(wr, std::fabs(Rx[xs][m](i, j) - fdm(Rm, i, j)));
                    }
                }
        }
        check_max("dS/dx matrices", ws, 1e-9);
        check_max("dh/dx matrices", wh, 1e-7);
        check_max("G^x(D) = sum D [(mn|ls)^x - (ml|ns)^x / 2] matrices", wg, 1e-7);
        check_max("d<m|r - O|n>/dx matrices", wr, 1e-9);
    }

    // ---- second derivatives against finite differences of the analytic first derivatives
    t0 = Clock::now();
    const Matrix hS = overlap_hessian(b, nat, D);
    const Matrix hH = core_hamiltonian_hessian(b, mol, D);
    const Matrix hE = eri_hessian(b, nat, gamma);
    std::printf("     analytic second derivatives: %.2f s\n", since(t0));
    double ws = 0, wh = 0, we = 0, asym = 0;
    for (int x = 0; x < 3 * nat; ++x) {
        const bool eri_row = full_mode() || x % 4 == 0;  // rows 0, 4, 8 (one per atom) by default
        std::vector<std::vector<double>> s(4), hh(4), e(4);
        const double steps[4] = {h, -h, 2 * h, -2 * h};
        for (int k = 0; k < 4; ++k) {
            auto g = at(x, steps[k]);
            s[static_cast<std::size_t>(k)] = overlap_gradient(g.second, nat, D);
            hh[static_cast<std::size_t>(k)] = core_hamiltonian_gradient(g.second, g.first, D);
            if (eri_row) e[static_cast<std::size_t>(k)] = eri_gradient(g.second, nat, gamma);
        }
        for (int y = 0; y < 3 * nat; ++y) {
            const std::size_t ys = static_cast<std::size_t>(y);
            auto f = [&](const std::vector<std::vector<double>>& v) {
                return (8.0 * (v[0][ys] - v[1][ys]) - (v[2][ys] - v[3][ys])) / (12.0 * h);
            };
            const std::size_t xs = static_cast<std::size_t>(x);
            ws = std::max(ws, std::fabs(hS(xs, ys) - f(s)) / std::max(1.0, std::fabs(hS(xs, ys))));
            wh = std::max(wh, std::fabs(hH(xs, ys) - f(hh)) / std::max(1.0, std::fabs(hH(xs, ys))));
            if (eri_row) we = std::max(we, std::fabs(hE(xs, ys) - f(e)) / std::max(1.0, std::fabs(hE(xs, ys))));
            asym = std::max(asym, std::fabs(hE(xs, ys) - hE(ys, xs)));
        }
    }
    check_max("sum D d2S/dxdy (relative)", ws, 1e-9);
    // The finite differences of the nuclear-attraction gradient (values up to ~50 with tight
    // core exponents) have a rounding floor of a few 1e-9 that does not shrink with the step.
    check_max("sum D d2h/dxdy (relative)", wh, 1e-8);
    check_max("sum Gamma d2(mn|ls)/dxdy (relative)", we, 1e-8);
    check_max("ERI second derivatives symmetric", asym, 1e-12);
}

// ---------------------------------------------------------------- RHF gradient

static RHFOptions tight_options(double lindep = 1e-7) {
    RHFOptions o;
    o.energy_tol = 1e-12;
    o.gradient_tol = 1e-10;
    o.max_iter = 300;
    o.lindep_tol = lindep;
    return o;
}

static Molecule molecule_from_json(const json& atoms) {
    Molecule mol;
    for (const auto& a : atoms) {
        Atom at;
        at.symbol = a[0].get<std::string>();
        at.Z = atomic_number(at.symbol);
        at.position = {a[1][0].get<double>(), a[1][1].get<double>(), a[1][2].get<double>()};
        mol.atoms.push_back(at);
    }
    mol.name = mol.atoms.size() == 3 ? "H2O" : (mol.atoms[0].Z == 3 ? "LiH" : "H2");
    return mol;
}

static Matrix density_of(const RHFResult& r) {
    const std::size_t n = r.C.rows();
    Matrix P(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < r.n_occupied; ++k) s += r.C(i, static_cast<std::size_t>(k)) * r.C(j, static_cast<std::size_t>(k));
            P(i, j) = 2.0 * s;
        }
    return P;
}

static void gradient_fd_case(const std::string& label, const Molecule& mol, const std::string& basis_name,
                             double lindep, double tol) {
    const BasisSet b = build_basis(mol, basis_name, default_basis_directory());
    const RHFOptions o = tight_options(lindep);
    const RHFResult r = run_rhf(mol, b, o);
    CHECK(r.converged);
    const RHFGradientResult g = rhf_gradient(mol, b, r, 0, lindep);
    double worst = 0.0, lind = 0.0;
    bool removed_changed = false;
    for (std::size_t x = 0; x < g.gradient.size(); ++x) {
        const double fd = fd4(
            [&](double d) {
                const Molecule m = displaced(mol, static_cast<int>(x), d);
                const RHFResult rr = run_rhf(m, basis_at_geometry(b, m), o);
                if (rr.removed_functions != r.removed_functions) removed_changed = true;
                return rr.total_energy;
            },
            2.0e-3);
        worst = std::max(worst, std::fabs(g.gradient[x] - fd));
        lind = std::max(lind, std::fabs(g.lindep[x]));
    }
    CHECK(!removed_changed);
    check_max(label + ": RHF gradient vs FD of the energy", worst, tol);
    if (r.removed_functions > 0)
        std::printf("     %d functions removed; largest orbital-space rotation term %.2e E_h/bohr\n",
                    r.removed_functions, lind);
}

static void test_rhf_gradient_fd() {
    std::printf("\n== RHF gradient against finite differences of the energy\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/gradient_reference.json"));
    auto geometry = [&](const std::string& name) {
        for (const auto& c : ref["cases"])
            if (c["geometry"] == name) return molecule_from_json(c["atoms"]);
        throw std::runtime_error("geometry " + name + " not in gradient_reference.json");
    };
    gradient_fd_case("H2O distorted aug-cc-pVDZ", geometry("H2O-distorted"), "aug-cc-pVDZ", 1e-7, 1e-8);
    gradient_fd_case("H2O distorted aug-cc-pCVDZ", geometry("H2O-distorted"), "aug-cc-pCVDZ", 1e-7, 1e-8);
    gradient_fd_case("LiH R = 3.10 aug-cc-pCVDZ", geometry("LiH-3.10"), "aug-cc-pCVDZ", 1e-7, 1e-8);
    gradient_fd_case("H2 R = 1.30 aug-cc-pVTZ", geometry("H2-1.30"), "aug-cc-pVTZ", 1e-7, 1e-8);
    // Removed near-linear dependencies: threshold 1e-2 removes two combinations.
    gradient_fd_case("H2O distorted aug-cc-pVDZ, lindep 1e-2", geometry("H2O-distorted"), "aug-cc-pVDZ", 1e-2, 1e-8);
    if (full_mode()) gradient_fd_case("H2O distorted aug-cc-pVTZ", geometry("H2O-distorted"), "aug-cc-pVTZ", 1e-7, 1e-8);
}

// ---------------------------------------------------------------- against PySCF

static void test_pyscf() {
    std::printf("\n== RHF energy, gradient, dipole moment and Hessian against PySCF\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/gradient_reference.json"));
    std::printf("     reference: %s\n", ref["program"].get<std::string>().c_str());
    for (const auto& c : ref["cases"]) {
        const std::string label = c["geometry"].get<std::string>() + " " + c["basis"].get<std::string>() +
                                  (c["lindep_threshold"].get<double>() > 1e-6 ? " (lindep 1e-2)" : "");
        if (c["nbf"].get<int>() > max_nbf()) {
            std::printf("     skipped %s (%d functions; PSILAB_TEST_MAX_NBF=%d)\n", label.c_str(), c["nbf"].get<int>(),
                        max_nbf());
            continue;
        }
        const Molecule mol = molecule_from_json(c["atoms"]);
        const double lindep = c["lindep_threshold"].get<double>();
        const BasisSet b = build_basis(mol, c["basis"].get<std::string>(), default_basis_directory());
        CHECK(b.n_functions == c["nbf"].get<int>());
        const auto t0 = Clock::now();
        const RHFResult r = run_rhf(mol, b, tight_options(lindep));
        CHECK(r.converged);
        CHECK(r.n_mo == c["orbitals"].get<int>());
        check_near(label + ": energy", r.total_energy, c["energy"].get<double>(), 1e-9);
        const RHFGradientResult g = rhf_gradient(mol, b, r, 0, lindep);
        const auto ref_g = c["gradient"].get<std::vector<double>>();
        double w = 0.0;
        for (std::size_t x = 0; x < ref_g.size(); ++x) w = std::max(w, std::fabs(g.gradient[x] - g.lindep[x] - ref_g[x]));
        check_max(label + (lindep > 1e-6 ? ": gradient without the rotation term" : ": gradient"), w, 1e-8);
        const std::array<double, 3> mu = dipole_moment(mol, b, density_of(r));
        const auto ref_mu = c["dipole"].get<std::vector<double>>();
        double wm = 0.0;
        for (int k = 0; k < 3; ++k) wm = std::max(wm, std::fabs(mu[static_cast<std::size_t>(k)] - ref_mu[static_cast<std::size_t>(k)]));
        check_max(label + ": dipole moment", wm, 1e-8);
        double th = 0.0;
        if (c.contains("hessian")) {
            const auto t1 = Clock::now();
            const RHFHessianResult H = rhf_hessian(mol, b, r, 0, 1024.0, 1e-10);
            th = since(t1);
            CHECK(H.cphf_converged);
            double wh = 0.0;
            for (std::size_t i = 0; i < H.hessian.rows(); ++i)
                for (std::size_t j = 0; j < H.hessian.cols(); ++j)
                    wh = std::max(wh, std::fabs(H.hessian(i, j) - c["hessian"][i][j].get<double>()));
            check_max(label + ": Hessian", wh, 2e-8);
        }
        std::printf("     %s: %d functions, RHF + gradient%s %.2f s (Hessian %.2f s)\n", label.c_str(), b.n_functions,
                    c.contains("hessian") ? " + Hessian" : "", since(t0), th);
    }
}

// ---------------------------------------------------------------- RHF Hessian

static void test_rhf_hessian_fd() {
    std::printf("\n== RHF Hessian and dipole derivatives against finite differences of analytic gradients\n");
    const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/gradient_reference.json"));
    Molecule mol;
    for (const auto& c : ref["cases"])
        if (c["geometry"] == "H2O-distorted") mol = molecule_from_json(c["atoms"]);
    const BasisSet b = build_basis(mol, "aug-cc-pVDZ", default_basis_directory());
    const RHFOptions o = tight_options();
    const RHFResult r = run_rhf(mol, b, o);
    const RHFHessianResult H = rhf_hessian(mol, b, r, 0, 1024.0, 1e-10);
    std::printf("     CPHF: %d iterations, residual %.1e; integrals %.2f s, CPHF %.2f s\n", H.cphf_iterations,
                H.cphf_residual, H.integral_seconds, H.cphf_seconds);
    const std::size_t n3 = H.hessian.rows();
    double wh = 0.0, wd = 0.0;
    const double h = 2.0e-3;
    for (std::size_t x = 0; x < n3; ++x) {
        std::vector<std::vector<double>> g(4);
        std::vector<std::array<double, 3>> mu(4);
        const double steps[4] = {h, -h, 2 * h, -2 * h};
        for (int k = 0; k < 4; ++k) {
            const Molecule m = displaced(mol, static_cast<int>(x), steps[k]);
            const BasisSet bb = basis_at_geometry(b, m);
            const RHFResult rr = run_rhf(m, bb, o);
            g[static_cast<std::size_t>(k)] = rhf_gradient(m, bb, rr).gradient;
            mu[static_cast<std::size_t>(k)] = dipole_moment(m, bb, density_of(rr));
        }
        for (std::size_t y = 0; y < n3; ++y)
            wh = std::max(wh, std::fabs(H.hessian(x, y) - (8.0 * (g[0][y] - g[1][y]) - (g[2][y] - g[3][y])) / (12.0 * h)));
        for (std::size_t m = 0; m < 3; ++m)
            wd = std::max(wd, std::fabs(H.dipole_derivatives[x][m] -
                                        (8.0 * (mu[0][m] - mu[1][m]) - (mu[2][m] - mu[3][m])) / (12.0 * h)));
    }
    check_max("H2O distorted aug-cc-pVDZ: Hessian vs FD of gradients", wh, 1e-8);
    check_max("H2O distorted aug-cc-pVDZ: dipole derivatives vs FD of dipoles", wd, 1e-8);
    // Translational invariance: sum over the atoms B of H(A k, B l) = 0, sum of dmu/dR_A = 0.
    double wt = 0.0, wq = 0.0;
    const std::size_t nat = n3 / 3;
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t l = 0; l < 3; ++l) {
            double s = 0.0;
            for (std::size_t B = 0; B < nat; ++B) s += H.hessian(x, 3 * B + l);
            wt = std::max(wt, std::fabs(s));
        }
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t m = 0; m < 3; ++m) {
            double s = 0.0;
            for (std::size_t A = 0; A < nat; ++A) s += H.dipole_derivatives[3 * A + k][m];
            wq = std::max(wq, std::fabs(s));
        }
    check_max("translational sum rule of the Hessian", wt, 1e-8);
    check_max("translational sum rule of the dipole derivatives", wq, 1e-8);
}

int main() {
    try {
        const auto t0 = Clock::now();
        test_derivative_integrals();
        test_rhf_gradient_fd();
        test_pyscf();
        test_rhf_hessian_fd();
        std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t0));
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    return g_failures == 0 ? 0 : 1;
}
