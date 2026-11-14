// Tests of the job types, the output file and the C API of ΨLab 5.0 (CTest):
//  - job-type spellings, input validation (HeH+ only runs the energy job);
//  - every job type through the C API for H2, LiH and H2O with RHF/aug-cc-pVDZ, against the
//    PySCF + geomeTRIC geometries (optimization_reference.json) and PySCF frequencies and
//    thermochemistry (frequency_reference.json);
//  - opt+freq with every feasible correlated method for H2 and LiH; exact identities:
//    H2: CCSD = CCSDT = CCSDTQ = Full CI and (T) = 0 (2 electrons); LiH: CCSDTQ = Full CI
//    (4 electrons); numerical-gradient methods against analytic ones;
//  - H2O opt+freq with CCSD and CCSD(T) (aug-cc-pVDZ) against PySCF + geomeTRIC;
//  - opt+freq+sp: the single point equals an energy job at the optimised geometry, and the
//    final G, H, U are the single-point energy plus the thermal corrections of the frequencies;
//  - refusals (H2O CCSDTQ frequencies and optimisation, H2O Full CI) with the estimate's counts
//    of energies and gradients of numerical Hessians;
//  - the streamed output file: read while the job runs with macqc_job_output, it equals the
//    "output" of the result; progress text; cancellation during an optimisation.
// PSILAB_TEST_MAX_NBF >= 100 (the manual "full" CI job) adds H2O CCSDT frequencies (numerical
// Hessian from 19 energies).
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "macqc/calculation.h"
#include "macqc/jobs.h"
#include "macqc/macqc_api.h"
#include "nlohmann/json.hpp"

using nlohmann::json;

static int g_checks = 0, g_failures = 0;

static void check_max(const std::string& name, double worst, double tol) {
    ++g_checks;
    const bool ok = worst <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-78s %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
}

#define CHECK(cond)                                                     \
    do {                                                                \
        ++g_checks;                                                     \
        if (!(cond)) {                                                  \
            ++g_failures;                                               \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                               \
    } while (0)

using Clock = std::chrono::steady_clock;
static double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

static bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static json run(const json& input, macqc_job* job = nullptr) {
    char* r = macqc_run_job(input.dump().c_str(), job);
    const json j = json::parse(r);
    macqc_free(r);
    return j;
}

static json run_logged(const std::string& label, const json& input) {
    const auto t0 = Clock::now();
    const json j = run(input);
    std::printf("     %-60s %6.1f s%s\n", label.c_str(), since(t0),
                j["error"].is_null() ? "" : (" error: " + j["error"].get<std::string>()).c_str());
    return j;
}

static const json* find_case(const json& cases, const std::string& mol, const std::string& basis, const std::string& method) {
    for (const auto& c : cases)
        if (c["molecule"] == mol && c["basis"] == basis && c["method"] == method) return &c;
    return nullptr;
}

constexpr double kBohr = macqc::kBohrInAngstrom;  // as in the reference files (CODATA 2022)

static double angstrom(const json& g, const char* key) { return g[key].get<double>() * kBohr; }

int main() {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);  // progress visible in CI logs
    // full mode (manual "full" CI job): PSILAB_TEST_MAX_NBF >= 100, as for the other test programs
    const char* max_nbf = std::getenv("PSILAB_TEST_MAX_NBF");
    const bool full = max_nbf && std::atoi(max_nbf) >= 100;
    const std::string data = MACQC_TEST_DATA_DIR;
    const json opt_ref = json::parse(read_file(data + "/optimization_reference.json"))["cases"];
    const json freq_ref = json::parse(read_file(data + "/frequency_reference.json"))["cases"];

    // ---- job types and validation
    CHECK(macqc::canonical_job_type("OPT FREQ") == "opt+freq");
    CHECK(macqc::canonical_job_type("optfreq") == "opt+freq");
    CHECK(macqc::canonical_job_type("Opt+Freq+SP") == "opt+freq+sp");
    CHECK(macqc::canonical_job_type("5") == "opt+freq+sp");
    CHECK(macqc::canonical_job_type("frequencies") == "freq");
    CHECK(macqc::canonical_job_type("sp") == "energy");
    CHECK(macqc::canonical_job_type("scan").empty());
    {
        json bad = run({{"molecule", "H2"}, {"job", "scan"}});
        CHECK(contains(bad["error"].get<std::string>(), "job must be energy, opt, freq, opt+freq or opt+freq+sp"));
        bad = run({{"molecule", "HeH+"}, {"job", "opt"}});
        CHECK(contains(bad["error"].get<std::string>(), "HeH+ (STO-nG) supports only the energy job"));
        bad = run({{"molecule", "H2"}, {"job", "freq"}, {"temperature", -5.0}});
        CHECK(contains(bad["error"].get<std::string>(), "temperature"));
        bad = run({{"molecule", "H2"}, {"job", "opt"}, {"optConvergence", "loose"}});
        CHECK(contains(bad["error"].get<std::string>(), "optConvergence must be normal or tight"));
        bad = run({{"molecule", "H2"}, {"job", "opt+freq+sp"}, {"spMethod", "MP2"}});
        CHECK(contains(bad["error"].get<std::string>(), "spMethod"));
        // HeH+ energy: the 4.0 path, no job output
        const json heh = run({{"molecule", "HeH+"}});
        CHECK(heh["error"].is_null() && !heh["hf"].is_null());
        CHECK(contains(heh["report"].get<std::string>(), "STO-3G FOR ATOMIC NUMBERS"));
    }

    // ---- every job type, RHF/aug-cc-pVDZ, against PySCF
    std::printf("Job types with RHF/aug-cc-pVDZ (PySCF + geomeTRIC, PySCF Hessians and thermochemistry)\n");
    const std::vector<std::string> jobs = {"energy", "opt", "freq", "opt+freq", "opt+freq+sp"};
    struct Mol {
        std::string name, basis;
    };
    for (const Mol& mol : {Mol{"H2", "aug-cc-pVTZ"}, Mol{"LiH", "aug-cc-pVDZ"}, Mol{"H2O", "aug-cc-pVDZ"}}) {
        const json* oref = find_case(opt_ref, mol.name, mol.basis, "RHF");
        const json* fref = find_case(freq_ref, mol.name, mol.basis, "RHF");
        // geometry of the frequency reference (bohr, degrees)
        double fR = 0.0, fangle = 104.52;
        if (fref) {
            const auto& a = (*fref)["atoms"];
            auto pos = [&](int i, int k) { return a[i][1][k].get<double>(); };
            auto dist = [&](int i, int j) {
                double d = 0.0;
                for (int k = 0; k < 3; ++k) d += (pos(i, k) - pos(j, k)) * (pos(i, k) - pos(j, k));
                return std::sqrt(d);
            };
            fR = dist(0, 1);
            if (mol.name == "H2O") {
                double dot = 0.0;
                for (int k = 0; k < 3; ++k) dot += (pos(1, k) - pos(0, k)) * (pos(2, k) - pos(0, k));
                fangle = std::acos(dot / (dist(0, 1) * dist(0, 2))) * 180.0 / 3.14159265358979323846;
            }
        }
        auto compare_frequencies = [&](const std::string& label, const json& r, double tol_cm) {
            const json& f = r["frequencies"];
            const auto& fr = (*fref)["frequencies"];
            CHECK(f["modes"].size() == fr.size());
            double worst = 0.0;
            for (std::size_t k = 0; k < fr.size() && k < f["modes"].size(); ++k)
                worst = std::max(worst, std::fabs(f["modes"][k]["wavenumber"].get<double>() - fr[k].get<double>()));
            check_max("  " + label + ": wavenumbers vs PySCF analytic Hessian (cm-1)", worst, tol_cm);
            const json& t = r["thermochemistry"];
            const json& pt = (*fref)["pyscf_thermo"];
            constexpr double kJmolPerEh = 2625499.639479;  // J/mol per E_h (CODATA 2018)
            check_max("  " + label + ": ZPE vs PySCF (E_h)",
                      std::fabs(t["zeroPointEnergy"].get<double>() / kJmolPerEh - pt["ZPE"].get<double>()), 3e-7);
            check_max("  " + label + ": S vs PySCF (J/(mol K))",
                      std::fabs(t["terms"]["total"]["entropy"].get<double>() - pt["S_tot"].get<double>() * kJmolPerEh), 0.01);
            check_max("  " + label + ": G vs PySCF (E_h)", std::fabs(t["totals"]["gibbs"].get<double>() - pt["G_tot"].get<double>()),
                      3e-6);
        };
        for (const auto& job : jobs) {
            json in = {{"molecule", mol.name}, {"basis", mol.basis}, {"method", "HF"}, {"job", job}, {"optConvergence", "tight"}};
            if (job == "freq" && fref) {
                in["R"] = fR;
                in["angle"] = fangle;
            }
            const json r = run_logged(mol.name + " RHF/" + mol.basis + " " + job, in);
            CHECK(r["error"].is_null());
            CHECK(r["job"]["type"] == job);
            CHECK(r["job"]["termination"] == "normal");
            const std::string out = r["output"].is_string() ? r["output"].get<std::string>() : "";
            CHECK(contains(out, "Normal termination of"));
            CHECK(contains(out, "---- Input"));
            CHECK(contains(out, "Atomic charges (SCF density)"));
            CHECK(contains(out, "Orbital energies"));
            CHECK(r["optimization"].is_null() == !macqc::job_has_optimization(job));
            CHECK(r["frequencies"].is_null() == !macqc::job_has_frequencies(job));
            CHECK(r["thermochemistry"].is_null() == !macqc::job_has_frequencies(job));
            CHECK(r["singlePoint"].is_null() == !macqc::job_has_single_point(job));
            if (macqc::job_has_optimization(job)) {
                const json& o = r["optimization"];
                CHECK(o["converged"].get<bool>());
                CHECK(contains(out, "Stationary point found"));
                CHECK(o["history"].size() == static_cast<std::size_t>(o["steps"].get<int>()));
                CHECK(o["history"][0]["deltaE"].is_null());
                if (oref) {
                    const json& g = o["finalGeometry"];
                    double dr = 0.0;
                    if (mol.name == "H2O") {
                        dr = std::max(std::fabs(angstrom(g, "R1") - (*oref)["R1"].get<double>()),
                                      std::fabs(angstrom(g, "R2") - (*oref)["R2"].get<double>()));
                        check_max("  " + mol.name + " " + job + ": angle vs PySCF + geomeTRIC (deg)",
                                  std::fabs(g["angle"].get<double>() - (*oref)["angle"].get<double>()), 2e-3);
                        CHECK(o["symmetricEnd"].get<bool>());
                    } else {
                        dr = std::fabs(angstrom(g, "R") - (*oref)["R"].get<double>());
                    }
                    check_max("  " + mol.name + " " + job + ": bond lengths vs PySCF + geomeTRIC (angstrom)", dr, 2e-5);
                    check_max("  " + mol.name + " " + job + ": energy vs PySCF + geomeTRIC (E_h)",
                              std::fabs(o["finalEnergy"].get<double>() - (*oref)["energy"].get<double>()), 1e-9);
                }
            }
            if (macqc::job_has_frequencies(job) && fref) {
                CHECK(r["frequencies"]["analytic"].get<bool>());
                CHECK(contains(out, "Thermochemistry (ideal gas, rigid rotor, harmonic oscillators)"));
                CHECK(contains(out, "Zero-point energy"));
                if (job == "freq") {
                    compare_frequencies(mol.name + " freq", r, 0.01);
                } else {
                    // at ΨLab's optimised geometry, if PySCF's reference is at the same minimum
                    const json& g = r["optimization"]["finalGeometry"];
                    if (std::fabs(g["R"].get<double>() - fR) * kBohr < 1e-4) compare_frequencies(mol.name + " " + job, r, 0.05);
                }
            }
            if (job == "opt+freq+sp") {
                CHECK(r["singlePoint"]["method"] == "RHF");
                const json& t = r["thermochemistry"];
                // same method and basis: the single point reproduces the optimised energy
                check_max("  " + mol.name + " opt+freq+sp: single point = optimised energy (E_h)",
                          std::fabs(r["singlePoint"]["energy"].get<double>() - r["optimization"]["finalEnergy"].get<double>()),
                          1e-9);
                check_max("  " + mol.name + " opt+freq+sp: G = E(SP) + thermal correction (E_h)",
                          std::fabs(t["totals"]["gibbs"].get<double>() -
                                    (r["job"]["finalEnergy"].get<double>() + t["corrections"]["gibbs"].get<double>() / 2625499.639479)),
                          1e-10);
            }
        }
    }
    {
        // frequencies at a non-stationary geometry: a warning in the output and the result
        const json r = run({{"molecule", "H2O"}, {"method", "HF"}, {"job", "freq"}});
        CHECK(r["error"].is_null());
        CHECK(contains(r["output"].get<std::string>(), "not a stationary point"));
        bool warned = false;
        for (const auto& w : r["job"]["warnings"])
            if (contains(w.get<std::string>(), "not a stationary point")) warned = true;
        CHECK(warned);
    }

    // ---- opt+freq+sp with a higher-level single point
    {
        std::printf("opt+freq+sp: H2O RHF/aug-cc-pVDZ geometry and frequencies, CCSD(T)/aug-cc-pVTZ single point\n");
        const json r = run_logged("H2O opt+freq+sp HF -> CCSD(T)/aug-cc-pVTZ",
                                  {{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "HF"}, {"job", "opt+freq+sp"},
                                   {"spMethod", "CCSD(T)"}, {"spBasis", "aug-cc-pVTZ"}, {"optConvergence", "tight"}});
        CHECK(r["error"].is_null());
        CHECK(r["singlePoint"]["method"] == "CCSD(T)");
        CHECK(r["job"]["finalEnergyLabel"] == "CCSD(T)/aug-cc-pVTZ");
        const json& g = r["optimization"]["finalGeometry"];
        const json e = run({{"molecule", "H2O"}, {"basis", "aug-cc-pVTZ"}, {"method", "CCSD(T)"}, {"R", g["R"]}, {"angle", g["angle"]}});
        // the energy job converges the amplitudes to 1e-8, the job's single point to 1e-9
        check_max("  single point = CCSD(T)/aug-cc-pVTZ energy job at the optimised geometry (E_h)",
                  std::fabs(r["singlePoint"]["energy"].get<double>() - e["coupledCluster"]["totalEnergy"].get<double>()), 5e-9);
        const json& t = r["thermochemistry"];
        const double esp = r["singlePoint"]["energy"].get<double>();
        double worst = 0.0;
        for (const char* k : {"zeroPoint", "thermalEnergy", "enthalpy", "gibbs", "helmholtz"})
            worst = std::max(worst, std::fabs(t["totals"][k].get<double>() - (esp + t["corrections"][k].get<double>() / 2625499.639479)));
        check_max("  E+ZPE, U, H, G, A = E(CCSD(T)/aug-cc-pVTZ) + corrections of RHF/aug-cc-pVDZ (E_h)", worst, 1e-10);
        // the 4.0 fields describe the final (single-point) calculation
        CHECK(r["coupledCluster"]["method"] == "CCSD(T)" && r["basis"]["name"] == "aug-cc-pVTZ");
        const std::string out = r["output"].get<std::string>();
        CHECK(contains(out, "Single-point energy at the optimised geometry: CCSD(T)/aug-cc-pVTZ"));
    }

    // ---- correlated methods for H2 and LiH: exact identities
    // H2 (2 electrons): opt+freq with every method; CCSD = CCSDT = CCSDTQ = Full CI, (T) = 0.
    // LiH (4 electrons): opt+freq with CCSD, CCSD(T) and Full CI; in the full mode also CCSDT and
    // CCSDTQ (numerical gradients and Hessians from energies) with CCSDTQ = Full CI.
    std::printf("Correlated methods: H2 (2 electrons) and LiH (4 electrons), aug-cc-pVDZ\n");
    {
        const std::vector<std::string> methods = {"CCSD", "CCSD(T)", "CCSDT", "CCSDTQ", "FCI"};
        for (const std::string mol : {"H2", "LiH"}) {
            std::vector<json> res(methods.size());
            std::vector<bool> ran(methods.size(), false);
            // Full CI first: its geometry is the reference of the frequency-only runs
            for (int k = static_cast<int>(methods.size()) - 1; k >= 0; --k) {
                const std::string& m = methods[static_cast<std::size_t>(k)];
                const bool numerical = m == "CCSDT" || m == "CCSDTQ";
                // LiH CCSDT and CCSDTQ (about 90 s per energy in aug-cc-pVDZ): full mode only
                if (mol == "LiH" && numerical && !full) continue;
                json in = {{"molecule", mol}, {"basis", "aug-cc-pVDZ"}, {"method", m}, {"job", "opt+freq"},
                           {"optConvergence", "tight"}};
                const json r = run_logged(mol + " " + m + "/aug-cc-pVDZ " + in["job"].get<std::string>(), in);
                CHECK(r["error"].is_null());
                CHECK(r["frequencies"]["modes"].size() == 1);
                CHECK(r["frequencies"]["analytic"].get<bool>() == false);
                // 2 displaced energies or gradients; LiH from energies also 4 finite-field energies (dipole)
                CHECK(r["frequencies"][numerical ? "energyEvaluations" : "gradientEvaluations"].get<int>() ==
                      (numerical && mol == "LiH" ? 6 : 2));
                if (in["job"] == "opt+freq") {
                    CHECK(r["optimization"]["converged"].get<bool>());
                    CHECK(contains(r["optimization"]["gradient"].get<std::string>(), numerical ? "numerical" : "analytic"));
                }
                res[static_cast<std::size_t>(k)] = r;
                ran[static_cast<std::size_t>(k)] = true;
            }
            const json& fci = res.back();
            const double r_fci = fci["optimization"]["finalGeometry"]["R"].get<double>();
            const double w_fci = fci["frequencies"]["modes"][0]["wavenumber"].get<double>();
            const double e_fci = fci["optimization"]["finalEnergy"].get<double>();
            for (std::size_t k = 0; k + 1 < methods.size(); ++k) {
                if (!ran[k]) continue;
                const bool exact = mol == "H2" || methods[k] == "CCSDTQ";
                if (!exact) continue;
                const bool numerical = methods[k] == "CCSDT" || methods[k] == "CCSDTQ";
                const json& r = res[k];
                const double dw = std::fabs(r["frequencies"]["modes"][0]["wavenumber"].get<double>() - w_fci);
                if (!r["optimization"].is_null()) {
                    const double dr = std::fabs(r["optimization"]["finalGeometry"]["R"].get<double>() - r_fci) * kBohr;
                    const double de = std::fabs(r["optimization"]["finalEnergy"].get<double>() - e_fci);
                    check_max("  " + mol + " " + methods[k] + " = Full CI: R_e (angstrom)", dr, numerical ? 5e-5 : 1e-6);
                    check_max("  " + mol + " " + methods[k] + " = Full CI: energy at R_e (E_h)", de, numerical ? 1e-8 : 1e-9);
                } else {
                    check_max("  " + mol + " " + methods[k] + " = Full CI: energy at the Full CI R_e (E_h)",
                              std::fabs(r["frequencies"]["energy"].get<double>() - e_fci), 1e-9);
                }
                check_max("  " + mol + " " + methods[k] + " = Full CI: wavenumber (cm-1)", dw, numerical ? 0.5 : 0.05);
            }
            if (mol == "LiH") {
                // LiH CCSD(T) and Full CI against PySCF + geomeTRIC
                for (const std::string m : {"CCSD(T)", "FCI"}) {
                    const json* ref = find_case(opt_ref, "LiH", "aug-cc-pVDZ", m);
                    const json& r = res[m == "FCI" ? 4 : 1];
                    if (!ref) continue;
                    check_max("  LiH " + m + ": R_e vs PySCF + geomeTRIC (angstrom)",
                              std::fabs(angstrom(r["optimization"]["finalGeometry"], "R") - (*ref)["R"].get<double>()), 2e-5);
                }
            }
        }
    }

    // ---- H2O with CCSD and CCSD(T)
    std::printf("H2O opt+freq with CCSD and CCSD(T)/aug-cc-pVDZ (analytic gradients, Hessian from 6 gradients)\n");
    for (const std::string m : {"CCSD", "CCSD(T)"}) {
        const json r = run_logged("H2O " + m + "/aug-cc-pVDZ opt+freq",
                                  {{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", m}, {"job", "opt+freq"},
                                   {"optConvergence", "tight"}});
        CHECK(r["error"].is_null());
        const json* ref = find_case(opt_ref, "H2O", "aug-cc-pVDZ", m);
        if (ref) {
            const json& g = r["optimization"]["finalGeometry"];
            check_max("  H2O " + m + ": R(O-H) vs PySCF + geomeTRIC (angstrom)",
                      std::fabs(angstrom(g, "R1") - (*ref)["R1"].get<double>()), 2e-5);
            check_max("  H2O " + m + ": angle vs PySCF + geomeTRIC (deg)", std::fabs(g["angle"].get<double>() - (*ref)["angle"].get<double>()),
                      2e-3);
        }
        CHECK(r["frequencies"]["gradientEvaluations"].get<int>() == 6);
        CHECK(r["frequencies"]["intensities"].get<bool>());
        CHECK(r["frequencies"]["modes"].size() == 3);
        const auto& est = r["job"]["estimate"]["steps"];
        CHECK(est.size() == 2 && est[1]["evaluations"].get<int>() == 6 && est[1]["evaluationKind"] == "gradients");
    }

    // ---- refusals and estimates
    std::printf("Refusals and estimates\n");
    {
        // with 8 GB the CCSDTQ energy of H2O/aug-cc-pVDZ is possible (5.5 GB), its frequencies are not
        json r = run({{"molecule", "H2O"}, {"method", "CCSDTQ"}, {"job", "freq"}, {"fciMemoryMB", 8192}});
        CHECK(!r["error"].is_null() && contains(r["error"].get<std::string>(), "CCSDTQ frequencies for H2O is not possible"));
        CHECK(r["job"]["termination"] == "refused");
        CHECK(contains(r["output"].get<std::string>(), "NOT POSSIBLE"));
        CHECK(contains(r["output"].get<std::string>(), "Error termination"));
        r = run({{"molecule", "H2O"}, {"method", "CCSDTQ"}, {"job", "opt"}, {"estimateOnly", true}, {"fciMemoryMB", 8192}});
        CHECK(contains(r["job"]["estimate"]["refusal"].get<std::string>(), "CCSDTQ optimisation for H2O is not possible"));
        r = run({{"molecule", "H2O"}, {"method", "CCSDTQ"}, {"job", "energy"}, {"estimateOnly", true}, {"fciMemoryMB", 8192}});
        CHECK(r["job"]["estimate"]["refusal"].is_null());
        r = run({{"molecule", "H2O"}, {"method", "CCSDTQ"}, {"job", "energy"}, {"estimateOnly", true}});
        CHECK(contains(r["job"]["estimate"]["refusal"].get<std::string>(), "CCSDTQ is not possible here"));
        r = run({{"molecule", "H2O"}, {"method", "FCI"}, {"job", "opt+freq"}});
        CHECK(!r["error"].is_null() && contains(r["error"].get<std::string>(), "impossible"));
        // job 5 with an impossible single point is refused before anything runs
        r = run({{"molecule", "H2O"}, {"method", "HF"}, {"job", "opt+freq+sp"}, {"spMethod", "FCI"}});
        CHECK(!r["error"].is_null() && contains(r["error"].get<std::string>(), "single-point step is not possible"));
        CHECK(r["optimization"].is_null());
        // numerical Hessians: counts of energies and gradients
        r = run({{"molecule", "H2O"}, {"method", "CCSDT"}, {"job", "freq"}, {"estimateOnly", true}});
        auto st = r["job"]["estimate"]["steps"][0];
        CHECK(st["evaluationKind"] == "energies" && st["evaluations"].get<int>() == 19 && !st["analytic"].get<bool>());
        CHECK(contains(st["how"].get<std::string>(), "9 CCSDT energies") && contains(st["how"].get<std::string>(), "10 finite-field"));
        r = run({{"molecule", "LiH"}, {"method", "CCSD(T)"}, {"job", "freq"}, {"estimateOnly", true}});
        st = r["job"]["estimate"]["steps"][0];
        CHECK(st["evaluationKind"] == "gradients" && st["evaluations"].get<int>() == 2);
        r = run({{"molecule", "H2O"}, {"method", "HF"}, {"job", "opt+freq"}, {"estimateOnly", true}});
        st = r["job"]["estimate"]["steps"][1];
        CHECK(st["analytic"].get<bool>() && st["evaluationKind"] == "Hessians");
        CHECK(r["job"]["estimate"]["totalSeconds"].get<double>() > 0.0);
    }

    // ---- streaming, progress text, cancellation
    std::printf("Streamed output, progress text, cancellation\n");
    {
        macqc_job* job = macqc_job_create();
        std::string result;
        std::atomic<bool> done{false};
        std::thread worker([&] {
            char* r = macqc_run_job(json({{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "HF"}, {"job", "opt+freq"},
                                          {"R", 1.9}, {"angle", 112.0}})
                                        .dump()
                                        .c_str(),
                                    job);
            result = r;
            macqc_free(r);
            done = true;
        });
        std::string streamed;
        long long offset = 0;
        int pieces = 0;
        std::set<std::string> progress;
        while (true) {
            const bool finished = done.load();
            long long total = 0;
            char* t = macqc_job_output(job, offset, &total);
            const std::string piece = t ? t : "";
            macqc_free(t);
            if (!piece.empty()) ++pieces;
            streamed += piece;
            offset += static_cast<long long>(piece.size());
            CHECK(offset == total);
            char* p = macqc_job_progress_text(job);
            if (p && *p) progress.insert(std::string(p).substr(0, 18));
            macqc_free(p);
            if (finished) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
        }
        worker.join();
        const json r = json::parse(result);
        CHECK(r["error"].is_null());
        CHECK(streamed == r["output"].get<std::string>());
        std::printf("     output read in %d pieces while the job ran (%zu bytes)\n", pieces, streamed.size());
        CHECK(pieces >= 5);
        bool step_seen = false;
        for (const auto& p : progress)
            if (p.rfind("Optimisation step", 0) == 0) step_seen = true;
        CHECK(step_seen);
        macqc_job_free(job);
    }
    {
        // cancel during the optimisation (from another thread)
        macqc_job* job = macqc_job_create();
        std::string result;
        std::thread worker([&] {
            char* r = macqc_run_job(json({{"molecule", "H2O"}, {"basis", "aug-cc-pVTZ"}, {"method", "HF"}, {"job", "opt+freq"},
                                          {"R", 1.9}, {"angle", 112.0}})
                                        .dump()
                                        .c_str(),
                                    job);
            result = r;
            macqc_free(r);
        });
        const auto t0 = Clock::now();
        while (since(t0) < 120.0) {
            char* t = macqc_job_output(job, 0, nullptr);
            const bool step2 = t && contains(t, "Optimisation step 2");
            macqc_free(t);
            if (step2) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        const auto tc = Clock::now();
        macqc_job_cancel(job);
        worker.join();
        const double latency = since(tc);
        const json r = json::parse(result);
        CHECK(r["error"] == "calculation cancelled");
        char* t = macqc_job_output(job, 0, nullptr);
        CHECK(t && contains(t, "Calculation cancelled"));
        CHECK(t && !contains(t, "Harmonic frequencies"));
        macqc_free(t);
        check_max("  cancellation latency during an optimisation step (s)", latency, 5.0);
        macqc_job_free(job);
    }

    if (full) {
        std::printf("Full tests: H2O CCSDT/aug-cc-pVDZ frequencies (numerical Hessian from 19 energies)\n");
        const json r = run_logged("H2O CCSDT/aug-cc-pVDZ freq", {{"molecule", "H2O"}, {"method", "CCSDT"}, {"job", "freq"}});
        CHECK(r["error"].is_null());
        CHECK(r["frequencies"]["energyEvaluations"].get<int>() == 19);
    }

    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
