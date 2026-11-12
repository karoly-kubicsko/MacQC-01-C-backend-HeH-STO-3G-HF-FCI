#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>

#include "macqc/calculation.h"
#include "macqc/macqc_api.h"
#include "nlohmann/json.hpp"

namespace macqc {

using nlohmann::json;

namespace {

json matrix_json(const Matrix& m) {
    json rows = json::array();
    for (std::size_t i = 0; i < m.rows(); ++i) {
        json row = json::array();
        for (std::size_t j = 0; j < m.cols(); ++j) row.push_back(m(i, j));
        rows.push_back(row);
    }
    return rows;
}

json input_json(const Input& in) {
    if (!is_heh_input(in) && !canonical_molecule_name(in.molecule).empty())
        return {{"molecule", in.molecule},
                {"basis", in.basis},
                {"R", in.R},
                {"angle", in.angle},
                {"method", in.method},
                {"basisDirectory", in.basisDirectory},
                {"threads", in.threads},
                {"memoryMB", in.memoryMB},
                {"fciMemoryMB", in.fciMemoryMB},
                {"scfMaxIter", in.scfMaxIter},
                {"energyConvergence", in.energyConvergence},
                {"gradientConvergence", in.gradientConvergence},
                {"estimateOnly", in.estimateOnly},
                {"printLevel", in.printLevel}};
    return {{"R", in.R},           {"zeta1", in.zeta1},         {"zeta2", in.zeta2},
            {"ZA", in.ZA},         {"ZB", in.ZB},               {"nGaussians", in.nGaussians},
            {"exactErf", in.exactErf}, {"maxIter", in.maxIter}, {"convergence", in.convergence},
            {"printLevel", in.printLevel}};
}

double get_number(const json& j, const char* key, double fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_number()) throw std::runtime_error(std::string("\"") + key + "\" must be a number");
    return it->get<double>();
}

int get_int(const json& j, const char* key, int fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_number()) throw std::runtime_error(std::string("\"") + key + "\" must be an integer");
    const double v = it->get<double>();
    if (v != std::floor(v) || std::fabs(v) > 1.0e9)
        throw std::runtime_error(std::string("\"") + key + "\" must be an integer");
    return static_cast<int>(v);
}

std::string get_string(const json& j, const char* key, const std::string& fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_string()) throw std::runtime_error(std::string("\"") + key + "\" must be a string");
    return it->get<std::string>();
}

json molecular_json(const MolecularResult& m) {
    json mol;
    mol["name"] = m.molecule.name;
    mol["electrons"] = m.molecule.electrons();
    mol["charge"] = m.molecule.charge;
    mol["R"] = m.input.R;
    mol["angle"] = m.molecule.name == "H2O" ? json(m.input.angle) : json(nullptr);
    mol["nuclearRepulsion"] = m.molecule.nuclear_repulsion();
    json atoms = json::array();
    for (const Atom& a : m.molecule.atoms)
        atoms.push_back({{"symbol", a.symbol}, {"Z", a.Z}, {"position", {a.position[0], a.position[1], a.position[2]}}});
    mol["atoms"] = atoms;

    json basis;
    basis["name"] = m.basis.name;
    basis["functions"] = m.basis.n_functions;
    basis["shells"] = m.basis.shells.size();
    basis["maxL"] = m.basis.max_l;
    json ba = json::array();
    for (const auto& a : m.basis.atoms)
        ba.push_back({{"symbol", a.symbol},
                      {"basis", a.basis},
                      {"contraction", a.contraction},
                      {"primitives", a.primitives},
                      {"functions", a.functions}});
    basis["atoms"] = ba;

    const CostEstimate& e = m.estimate;
    json est = {{"basisFunctions", e.basis_functions},
                {"uniqueIntegrals", e.unique_eris},
                {"integralMemoryMB", e.eri_memory_mb},
                {"integralsStored", e.eri_stored},
                {"threads", e.threads},
                {"rhfSeconds", e.rhf_seconds},
                {"totalSeconds", e.total_seconds},
                {"longCalculation", e.total_seconds > kLongCalculationSeconds}};
    if (e.fci) {
        const FCISize& f = *e.fci;
        est["fci"] = {{"orbitals", f.orbitals},
                      {"alphaElectrons", f.alpha},
                      {"betaElectrons", f.beta},
                      {"alphaStrings", f.alpha_strings},
                      {"betaStrings", f.beta_strings},
                      {"determinants", f.determinants},
                      {"memoryMB", f.memory_bytes / (1024.0 * 1024.0)},
                      {"seconds", f.seconds},
                      {"feasible", f.feasible},
                      {"reason", f.reason.empty() ? json(nullptr) : json(f.reason)}};
    } else {
        est["fci"] = nullptr;
    }

    json out;
    out["molecule"] = mol;
    out["basis"] = basis;
    out["estimate"] = est;
    out["notes"] = m.notes;
    out["warnings"] = m.warnings;
    out["refusal"] = m.refusal.empty() ? json(nullptr) : json(m.refusal);
    if (m.fci) {
        const FCIResultGeneral& f = *m.fci;
        json fc;
        fc["converged"] = f.converged;
        fc["iterations"] = f.iterations;
        fc["totalEnergy"] = f.total_energy;
        fc["electronicEnergy"] = f.electronic_energy;
        fc["correlationEnergy"] = f.correlation_energy;
        fc["S2"] = f.s2;
        fc["hfWeight"] = f.hf_weight;
        fc["orbitals"] = f.orbitals;
        fc["alphaElectrons"] = f.alpha;
        fc["betaElectrons"] = f.beta;
        fc["alphaStrings"] = f.alpha_strings;
        fc["betaStrings"] = f.beta_strings;
        fc["determinants"] = f.determinants;
        json lead = json::array();
        for (const auto& d : f.leading)
            lead.push_back({{"coefficient", d.coefficient}, {"alpha", d.alpha}, {"beta", d.beta},
                            {"excitation", d.excitation}});
        fc["leadingDeterminants"] = lead;
        fc["naturalOccupations"] = f.natural_occupations;
        json hist = json::array();
        for (const auto& it : f.history)
            hist.push_back({{"iteration", it.iteration},
                            {"totalEnergy", it.energy},
                            {"residual", it.residual},
                            {"subspace", it.subspace}});
        fc["history"] = hist;
        fc["transformSeconds"] = f.transform_seconds;
        fc["ciSeconds"] = f.ci_seconds;
        out["fullCI"] = fc;
    } else {
        out["fullCI"] = nullptr;
    }
    if (m.rhf) {
        const RHFResult& h = *m.rhf;
        json r;
        r["converged"] = h.converged;
        r["iterations"] = h.iterations;
        r["totalEnergy"] = h.total_energy;
        r["electronicEnergy"] = h.electronic_energy;
        r["nuclearRepulsion"] = h.nuclear_repulsion;
        r["kineticEnergy"] = h.kinetic_energy;
        r["nuclearAttractionEnergy"] = h.nuclear_attraction_energy;
        r["electronRepulsionEnergy"] = h.two_electron_energy;
        r["virialRatio"] = h.virial_ratio;
        r["orbitalEnergies"] = h.orbital_energies;
        r["occupiedOrbitals"] = h.n_occupied;
        const int no = h.n_occupied, nmo = static_cast<int>(h.orbital_energies.size());
        r["homo"] = no >= 1 ? json(h.orbital_energies[no - 1]) : json(nullptr);
        r["lumo"] = no < nmo ? json(h.orbital_energies[no]) : json(nullptr);
        r["smallestOverlapEigenvalue"] = h.smallest_overlap_eigenvalue;
        r["removedFunctions"] = h.removed_functions;
        r["integralsStored"] = h.eri_stored;
        r["uniqueIntegrals"] = h.unique_integrals;
        r["threads"] = h.threads;
        r["integralSeconds"] = h.integral_seconds;
        r["scfSeconds"] = h.scf_seconds;
        json hist = json::array();
        for (const auto& it : h.history)
            hist.push_back({{"iteration", it.iteration},
                            {"totalEnergy", it.energy},
                            {"deltaE", it.delta_e},
                            {"gradient", it.gradient},
                            {"diis", it.diis}});
        r["scfHistory"] = hist;
        out["rhf"] = r;
    } else {
        out["rhf"] = nullptr;
    }
    return out;
}

bool get_bool(const json& j, const char* key, bool fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_boolean()) throw std::runtime_error(std::string("\"") + key + "\" must be a boolean");
    return it->get<bool>();
}

}  // namespace

Input parse_input_json(const std::string& text) {
    json j;
    try {
        j = json::parse(text);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("invalid JSON input: ") + e.what());
    }
    if (!j.is_object()) throw std::runtime_error("input JSON must be an object");
    Input in;
    in.molecule = get_string(j, "molecule", in.molecule);
    const std::string canonical = canonical_molecule_name(in.molecule);
    if (canonical.empty()) throw std::runtime_error("\"molecule\" must be HeH+, H2, LiH or H2O");
    in.molecule = canonical;
    if (canonical != "HeH+") in.R = 0.0;  // default bond length of the molecule unless given
    in.R = get_number(j, "R", in.R);
    in.zeta1 = get_number(j, "zeta1", in.zeta1);
    in.zeta2 = get_number(j, "zeta2", in.zeta2);
    in.ZA = get_number(j, "ZA", in.ZA);
    in.ZB = get_number(j, "ZB", in.ZB);
    in.nGaussians = get_int(j, "nGaussians", in.nGaussians);
    in.exactErf = get_bool(j, "exactErf", in.exactErf);
    in.maxIter = get_int(j, "maxIter", in.maxIter);
    in.convergence = get_number(j, "convergence", in.convergence);
    in.printLevel = get_int(j, "printLevel", in.printLevel);
    in.basis = get_string(j, "basis", in.basis);
    in.angle = get_number(j, "angle", in.angle);
    in.method = get_string(j, "method", in.method);
    in.basisDirectory = get_string(j, "basisDirectory", in.basisDirectory);
    in.threads = get_int(j, "threads", in.threads);
    in.memoryMB = get_number(j, "memoryMB", in.memoryMB);
    in.fciMemoryMB = get_number(j, "fciMemoryMB", in.fciMemoryMB);
    in.scfMaxIter = get_int(j, "scfMaxIter", in.scfMaxIter);
    in.energyConvergence = get_number(j, "energyConvergence", in.energyConvergence);
    in.gradientConvergence = get_number(j, "gradientConvergence", in.gradientConvergence);
    in.estimateOnly = get_bool(j, "estimateOnly", in.estimateOnly);
    return in;
}

std::string output_to_json(const Output& out, int indent) {
    json j;
    j["input"] = input_json(out.input);
    j["error"] = out.error.empty() ? json(nullptr) : json(out.error);
    if (!is_heh_input(out.input) && !canonical_molecule_name(out.input.molecule).empty()) {
        // General molecules: {"input","error","molecule","basis","estimate","notes",
        // "warnings","refusal","rhf","fullCI","hf":null,"fci":null,"report"}.
        if (out.molecular) {
            const json m = molecular_json(*out.molecular);
            for (auto it = m.begin(); it != m.end(); ++it) j[it.key()] = it.value();
        } else {
            j["molecule"] = nullptr;
            j["basis"] = nullptr;
            j["estimate"] = nullptr;
            j["notes"] = json::array();
            j["warnings"] = json::array();
            j["refusal"] = nullptr;
            j["rhf"] = nullptr;
            j["fullCI"] = nullptr;
        }
        j["hf"] = nullptr;
        j["fci"] = nullptr;
        j["report"] = out.report;
        return j.dump(indent, ' ', false, json::error_handler_t::replace);
    }

    if (out.hf) {
        const HFResult& h = *out.hf;
        json hf;
        hf["converged"] = h.converged;
        hf["iterations"] = h.iterations;
        hf["electronicEnergy"] = h.electronic_energy;
        hf["nuclearRepulsion"] = h.nuclear_repulsion;
        hf["totalEnergy"] = h.total_energy;
        hf["orbitalEnergies"] = h.orbital_energies;
        hf["C"] = matrix_json(h.C);
        hf["P"] = matrix_json(h.P);
        hf["S"] = matrix_json(h.S);
        hf["X"] = matrix_json(h.X);
        hf["H"] = matrix_json(h.H);
        hf["F"] = matrix_json(h.F);
        hf["G"] = matrix_json(h.G);
        hf["PS"] = matrix_json(h.PS);
        const Integrals& v = h.integrals;
        hf["integrals"] = {{"S12", v.S12},     {"T11", v.T11},     {"T12", v.T12},   {"T22", v.T22},
                           {"V11A", v.V11A},   {"V12A", v.V12A},   {"V22A", v.V22A}, {"V11B", v.V11B},
                           {"V12B", v.V12B},   {"V22B", v.V22B},   {"V1111", v.V1111}, {"V2111", v.V2111},
                           {"V2121", v.V2121}, {"V2211", v.V2211}, {"V2221", v.V2221}, {"V2222", v.V2222}};
        json hist = json::array();
        for (const auto& s : h.history)
            hist.push_back({{"iteration", s.iteration},
                            {"electronicEnergy", s.electronic_energy},
                            {"deltaP", s.delta_p}});
        hf["scfHistory"] = hist;
        j["hf"] = hf;
    } else {
        j["hf"] = nullptr;
    }

    if (out.fci) {
        const FCIResult& f = *out.fci;
        json fci;
        fci["totalEnergy"] = f.total_energy;
        fci["electronicEnergy"] = f.electronic_energy;
        fci["correlationEnergy"] = f.correlation_energy;
        json states = json::array();
        for (std::size_t k = 0; k < f.states.size(); ++k)
            states.push_back({{"energy", f.states[k].energy},
                              {"S2", f.states[k].S2},
                              {"spin", f.states[k].spin},
                              {"vector", f.vectors[k]}});
        fci["states"] = states;
        fci["groundStateVector"] = f.ground_state_vector;
        fci["determinants"] = f.determinants;
        json w = json::array();
        for (const auto& s : f.singlet_weights)
            w.push_back({{"configuration", s.configuration}, {"weight", s.weight}});
        fci["singletWeights"] = w;
        j["fci"] = fci;
    } else {
        j["fci"] = nullptr;
    }

    j["report"] = out.report;
    return j.dump(indent, ' ', false, json::error_handler_t::replace);
}

}  // namespace macqc

extern "C" {

char* macqc_run(const char* input_json) {
    std::string result;
    try {
        macqc::Output out;
        try {
            const macqc::Input in = macqc::parse_input_json(input_json ? input_json : "");
            out = macqc::run_calculation(in);
        } catch (const std::exception& e) {
            out = macqc::Output();
            out.error = e.what();
        }
        result = macqc::output_to_json(out);
    } catch (...) {
        result = "{\"error\":\"internal error\",\"hf\":null,\"fci\":null,\"report\":\"\"}";
    }
    char* buf = static_cast<char*>(std::malloc(result.size() + 1));
    if (!buf) return nullptr;
    std::memcpy(buf, result.c_str(), result.size() + 1);
    return buf;
}

void macqc_free(char* p) { std::free(p); }

}  // extern "C"
