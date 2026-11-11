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
    return in;
}

std::string output_to_json(const Output& out, int indent) {
    json j;
    j["input"] = input_json(out.input);
    j["error"] = out.error.empty() ? json(nullptr) : json(out.error);

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
