#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>

#include "macqc/calculation.h"
#include "macqc/job.h"
#include "macqc/jobs.h"
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
                {"printLevel", in.printLevel},
                {"job", in.job},
                {"spMethod", in.spMethod},
                {"spBasis", in.spBasis},
                {"temperature", in.temperature},
                {"pressure", in.pressure},
                {"optConvergence", in.optConvergence},
                {"maxSteps", in.maxSteps}};
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
    if (e.cc) {
        const CCSize& c = *e.cc;
        est["coupledCluster"] = {{"method", c.method},
                                 {"engine", c.engine},
                                 {"occupied", c.occupied},
                                 {"virtuals", c.virtuals},
                                 {"amplitudes", c.amplitudes},
                                 {"memoryMB", c.memory_bytes / (1024.0 * 1024.0)},
                                 {"ladderIntegralsStored", c.ladder_stored},
                                 {"iterations", c.iterations},
                                 {"iterationSeconds", c.iteration_seconds},
                                 {"triplesSeconds", c.triples_seconds},
                                 {"seconds", c.seconds},
                                 {"feasible", c.feasible},
                                 {"reason", c.reason.empty() ? json(nullptr) : json(c.reason)}};
    } else {
        est["coupledCluster"] = nullptr;
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
    if (m.cc) {
        const CCResult& c = *m.cc;
        json cj;
        cj["method"] = c.method;
        cj["engine"] = c.engine;
        cj["converged"] = c.converged;
        cj["iterations"] = c.iterations;
        cj["occupied"] = c.occupied;
        cj["virtuals"] = c.virtuals;
        cj["referenceEnergy"] = c.reference_energy;
        cj["correlationEnergy"] = c.correlation_energy;
        cj["ccEnergy"] = c.cc_energy;
        cj["triplesCorrection"] = c.has_triples ? json(c.triples_correction) : json(nullptr);
        cj["totalEnergy"] = c.total_energy;
        cj["t1Diagnostic"] = c.t1_diagnostic;
        cj["amplitudeCounts"] = c.amplitude_counts;
        cj["amplitudeNorms"] = c.rank_norms;
        json lead = json::array();
        for (const auto& a : c.leading)
            lead.push_back({{"value", a.value},
                            {"rank", a.rank},
                            {"occupied", a.occupied},
                            {"virtuals", a.virtuals},
                            {"spins", a.spins}});
        cj["leadingAmplitudes"] = lead;
        json hist = json::array();
        for (const auto& it : c.history)
            hist.push_back({{"iteration", it.iteration},
                            {"totalEnergy", it.energy},
                            {"correlationEnergy", it.correlation},
                            {"deltaE", it.delta_e},
                            {"residual", it.residual}});
        cj["history"] = hist;
        cj["ladder"] = c.ladder.empty() ? json(nullptr) : json(c.ladder);
        cj["transformSeconds"] = c.transform_seconds;
        cj["ccSeconds"] = c.cc_seconds;
        cj["triplesSeconds"] = c.has_triples ? json(c.triples_seconds) : json(nullptr);
        out["coupledCluster"] = cj;
    } else {
        out["coupledCluster"] = nullptr;
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

json atoms_json(const Molecule& m) {
    json atoms = json::array();
    for (const Atom& a : m.atoms)
        atoms.push_back({{"symbol", a.symbol}, {"Z", a.Z}, {"position", {a.position[0], a.position[1], a.position[2]}}});
    return atoms;
}

double distance(const Atom& a, const Atom& b) {
    double s = 0.0;
    for (int k = 0; k < 3; ++k) s += (a.position[k] - b.position[k]) * (a.position[k] - b.position[k]);
    return std::sqrt(s);
}

// The parameters of the inputs: R (bohr; H2O: both O-H distances and their mean) and the
// H-O-H angle (degrees).
json geometry_parameters_json(const Molecule& m) {
    json g;
    g["atoms"] = atoms_json(m);
    if (m.atoms.size() == 2) {
        g["R"] = distance(m.atoms[0], m.atoms[1]);
        g["angle"] = nullptr;
        return g;
    }
    int o = 0;
    for (std::size_t i = 0; i < m.atoms.size(); ++i)
        if (m.atoms[i].Z == 8) o = static_cast<int>(i);
    std::vector<int> h;
    for (std::size_t i = 0; i < m.atoms.size(); ++i)
        if (static_cast<int>(i) != o) h.push_back(static_cast<int>(i));
    if (h.size() != 2) {
        g["R"] = nullptr;
        g["angle"] = nullptr;
        return g;
    }
    const Atom &O = m.atoms[o], &H1 = m.atoms[h[0]], &H2 = m.atoms[h[1]];
    const double r1 = distance(O, H1), r2 = distance(O, H2);
    double dot = 0.0;
    for (int k = 0; k < 3; ++k) dot += (H1.position[k] - O.position[k]) * (H2.position[k] - O.position[k]);
    const double c = std::max(-1.0, std::min(1.0, dot / (r1 * r2)));
    g["R"] = 0.5 * (r1 + r2);
    g["R1"] = r1;
    g["R2"] = r2;
    g["angle"] = std::acos(c) * 180.0 / 3.14159265358979323846;
    return g;
}

json estimate_json(const JobEstimate& e) {
    json steps = json::array();
    for (const auto& s : e.steps)
        steps.push_back({{"name", s.name},
                         {"method", s.method},
                         {"basis", s.basis},
                         {"how", s.how},
                         {"analytic", s.analytic},
                         {"evaluations", s.evaluations},
                         {"evaluationKind", s.evaluation_kind},
                         {"seconds", s.seconds},
                         {"memoryMB", s.memory_mb}});
    return {{"steps", steps},
            {"totalSeconds", e.total_seconds},
            {"memoryMB", e.memory_mb},
            {"longCalculation", e.total_seconds > kLongCalculationSeconds},
            {"refusal", e.refusal.empty() ? json(nullptr) : json(e.refusal)},
            {"warnings", e.warnings},
            {"notes", e.notes}};
}

json point_json(const SinglePoint& p) {
    json j;
    j["method"] = p.method;
    j["energy"] = p.energy;
    j["converged"] = p.converged;
    j["rhfEnergy"] = p.rhf.total_energy;
    j["correlationEnergy"] = p.energy - p.rhf.total_energy;
    j["geometry"] = geometry_parameters_json(p.molecule);
    j["warnings"] = p.warnings;
    j["seconds"] = p.seconds;
    return j;
}

json optimization_json(const OptimizationResult& o) {
    json j;
    j["converged"] = o.converged;
    j["steps"] = o.steps;
    j["message"] = o.message;
    j["initialHessian"] = o.initial_hessian;
    j["criteria"] = {{"maxForce", o.criteria.max_force},
                     {"rmsForce", o.criteria.rms_force},
                     {"maxDisplacement", o.criteria.max_displacement},
                     {"rmsDisplacement", o.criteria.rms_displacement},
                     {"energy", o.criteria.energy}};
    j["initialEnergy"] = o.initial_energy;
    j["finalEnergy"] = o.final_energy;
    j["initialGeometry"] = geometry_parameters_json(o.initial);
    j["finalGeometry"] = geometry_parameters_json(o.final_geometry);
    auto params = [](const std::vector<GeometricParameter>& v) {
        json a = json::array();
        for (const auto& p : v) a.push_back({{"name", p.name}, {"value", p.value}, {"angle", p.angle}});
        return a;
    };
    j["initialParameters"] = params(o.initial_parameters);
    j["finalParameters"] = params(o.final_parameters);
    j["symmetricStart"] = o.symmetric_start;
    j["symmetricEnd"] = o.symmetric_end;
    j["symmetryDeviation"] = o.symmetry_deviation;
    j["gradient"] = o.history.empty() ? json(nullptr) : json(o.history.front().gradient_description);
    json hist = json::array();
    for (const auto& s : o.history) {
        json h;
        h["step"] = s.step;
        h["energy"] = s.energy;
        h["deltaE"] = s.step == 1 ? json(nullptr) : json(s.delta_e);
        h["geometry"] = geometry_parameters_json(s.geometry);
        json coords = json::array();
        for (std::size_t k = 0; k < s.names.size(); ++k)
            coords.push_back({{"name", s.names[k]},
                              {"definition", k < s.definitions.size() ? s.definitions[k] : ""},
                              {"value", k < s.q.size() ? s.q[k] : 0.0},
                              {"gradient", k < s.gq.size() ? s.gq[k] : 0.0},
                              {"step", k < s.dq.size() ? json(s.dq[k]) : json(nullptr)}});
        h["coordinates"] = coords;
        json conv = json::array();
        for (const auto& c : s.convergence)
            conv.push_back({{"name", c.name},
                            {"value", c.applicable ? json(c.value) : json(nullptr)},
                            {"threshold", c.threshold},
                            {"met", c.met},
                            {"applicable", c.applicable}});
        h["convergence"] = conv;
        h["converged"] = s.converged;
        h["hessianUpdate"] = s.hessian_update;
        h["trustRadius"] = s.trust_radius;
        h["stepNorm"] = s.step_norm;
        h["stepScaled"] = s.step_scaled;
        h["predictedChange"] = s.predicted_change;
        h["ratio"] = s.has_ratio ? json(s.ratio) : json(nullptr);
        h["warnings"] = s.warnings;
        h["seconds"] = s.seconds;
        hist.push_back(h);
    }
    j["history"] = hist;
    j["seconds"] = o.seconds;
    return j;
}

json frequencies_json(const HessianResult& hr, const HarmonicAnalysis& h) {
    json j;
    j["analytic"] = hr.analytic;
    j["description"] = hr.description;
    j["numericalReason"] = hr.numerical_reason.empty() ? json(nullptr) : json(hr.numerical_reason);
    j["gradientEvaluations"] = hr.gradient_evaluations;
    j["energyEvaluations"] = hr.energy_evaluations;
    j["coordinates"] = hr.coordinates;
    j["intensities"] = h.has_intensities;
    j["intensityDescription"] = hr.dipole_description;
    j["dipole"] = hr.has_dipole ? json({hr.dipole[0], hr.dipole[1], hr.dipole[2]}) : json(nullptr);
    j["linear"] = h.linear;
    j["vibrations"] = h.n_vib;
    j["masses"] = h.masses;
    j["isotopes"] = h.isotopes;
    j["geometry"] = geometry_parameters_json(hr.reference.point.molecule);
    j["energy"] = hr.reference.point.energy;
    j["lowFrequencies"] = h.low_frequencies;
    json modes = json::array();
    for (const auto& m : h.modes) {
        json d = json::array();
        for (const auto& v : m.displacement) d.push_back({v[0], v[1], v[2]});
        modes.push_back({{"wavenumber", m.wavenumber},
                         {"imaginary", m.imaginary},
                         {"symmetry", m.symmetry},
                         {"reducedMass", m.reduced_mass},
                         {"forceConstant", m.force_constant},
                         {"irIntensity", m.ir_intensity >= 0.0 ? json(m.ir_intensity) : json(nullptr)},
                         {"displacement", d}});
    }
    j["modes"] = modes;
    j["seconds"] = hr.seconds;
    return j;
}

json thermo_json(const ThermoResult& t) {
    auto term = [](const ThermoTerm& x) {
        return json{{"energy", x.energy}, {"entropy", x.entropy}, {"cv", x.cv}, {"q", x.q}};
    };
    auto arr3 = [](const std::array<double, 3>& a) { return json{a[0], a[1], a[2]}; };
    json j;
    j["temperature"] = t.temperature;
    j["pressure"] = t.pressure;
    j["symmetryNumber"] = t.symmetry_number;
    j["linear"] = t.linear;
    j["mass"] = t.mass;
    j["moments"] = arr3(t.moments);
    j["rotationalConstantsGHz"] = arr3(t.rotational_constants_ghz);
    j["rotationalConstantsCm"] = arr3(t.rotational_constants_cm);
    j["rotationalTemperatures"] = arr3(t.rotational_temperatures);
    j["vibrationalTemperatures"] = t.vibrational_temperatures;
    j["excludedModes"] = t.excluded_modes;
    j["warnings"] = t.warnings;
    // J/mol, J/(mol K): energies of the terms include the zero-point energy (vibration).
    j["terms"] = {{"translation", term(t.translation)},
                  {"rotation", term(t.rotation)},
                  {"vibration", term(t.vibration)},
                  {"electronic", term(t.electronic)},
                  {"total", term(t.total)}};
    j["zeroPointEnergy"] = t.zero_point_energy;  // J/mol
    j["corrections"] = {{"zeroPoint", t.zero_point_energy},
                        {"thermalEnergy", t.thermal_energy},
                        {"enthalpy", t.enthalpy},
                        {"gibbs", t.gibbs},
                        {"helmholtz", t.helmholtz}};  // J/mol
    j["electronicEnergy"] = t.electronic_energy;      // E_h
    j["totals"] = {{"zeroPoint", t.total_zpe},
                   {"thermalEnergy", t.total_u},
                   {"enthalpy", t.total_h},
                   {"gibbs", t.total_g},
                   {"helmholtz", t.total_a}};  // E_h
    return j;
}

json job_json(const JobResult& r, const std::string& type) {
    json j;
    j["type"] = type;
    j["title"] = job_title(type);
    j["termination"] = r.termination;
    j["completed"] = r.completed;
    j["message"] = r.message.empty() ? json(nullptr) : json(r.message);
    j["estimate"] = estimate_json(r.estimate);
    j["warnings"] = r.warnings;
    j["seconds"] = r.seconds;
    if (r.termination == "estimate") {
        j["initialGeometry"] = nullptr;
        j["finalGeometry"] = nullptr;
        j["finalEnergy"] = nullptr;
        j["finalEnergyLabel"] = nullptr;
    } else {
        j["initialGeometry"] = geometry_parameters_json(r.initial_geometry);
        j["finalGeometry"] = geometry_parameters_json(r.final_geometry);
        j["finalEnergy"] = r.final_energy_label.empty() ? json(nullptr) : json(r.final_energy);
        j["finalEnergyLabel"] = r.final_energy_label.empty() ? json(nullptr) : json(r.final_energy_label);
    }
    j["singlePointMethod"] = r.single_point_method.empty() ? json(nullptr) : json(r.single_point_method);
    j["singlePointBasis"] = r.single_point_basis.empty() ? json(nullptr) : json(r.single_point_basis);
    return j;
}

// "job", "optimization", "frequencies", "thermochemistry", "singlePoint", "output".
void add_job_json(json& j, const Output& out) {
    if (!out.job) {
        j["job"] = nullptr;
        j["optimization"] = nullptr;
        j["frequencies"] = nullptr;
        j["thermochemistry"] = nullptr;
        j["singlePoint"] = nullptr;
        j["output"] = nullptr;
        return;
    }
    const JobResult& r = *out.job;
    const std::string type = canonical_job_type(out.input.job).empty() ? "energy" : canonical_job_type(out.input.job);
    j["job"] = job_json(r, type);
    j["optimization"] = r.optimization ? optimization_json(*r.optimization) : json(nullptr);
    j["frequencies"] = (r.hessian && r.harmonic) ? frequencies_json(*r.hessian, *r.harmonic) : json(nullptr);
    j["thermochemistry"] = r.thermo ? thermo_json(*r.thermo) : json(nullptr);
    j["singlePoint"] = r.single_point ? point_json(*r.single_point) : json(nullptr);
    j["output"] = r.termination == "estimate" ? json(nullptr) : json(r.output);
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
    in.job = get_string(j, "job", in.job);
    in.spMethod = get_string(j, "spMethod", in.spMethod);
    in.spBasis = get_string(j, "spBasis", in.spBasis);
    in.temperature = get_number(j, "temperature", in.temperature);
    in.pressure = get_number(j, "pressure", in.pressure);
    in.optConvergence = get_string(j, "optConvergence", in.optConvergence);
    in.maxSteps = get_int(j, "maxSteps", in.maxSteps);
    return in;
}

std::string output_to_json(const Output& out, int indent) {
    json j;
    j["input"] = input_json(out.input);
    j["error"] = out.error.empty() ? json(nullptr) : json(out.error);
    if (!is_heh_input(out.input) && !canonical_molecule_name(out.input.molecule).empty()) {
        // General molecules: {"input","error","molecule","basis","estimate","notes",
        // "warnings","refusal","rhf","fullCI","coupledCluster","job","optimization",
        // "frequencies","thermochemistry","singlePoint","output","hf":null,"fci":null,"report"}.
        // rhf/fullCI/coupledCluster/report describe the final energy calculation of the job.
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
            j["coupledCluster"] = nullptr;
        }
        add_job_json(j, out);
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

struct macqc_job {
    macqc::Job job;
};

extern "C" {

char* macqc_run(const char* input_json) { return macqc_run_job(input_json, nullptr); }

macqc_job* macqc_job_create(void) {
    try {
        return new macqc_job();
    } catch (...) {
        return nullptr;
    }
}

void macqc_job_cancel(macqc_job* job) {
    if (job) job->job.cancelled.store(true);
}

void macqc_job_progress(const macqc_job* job, int* stage, int* iteration, double* energy, double* residual) {
    if (!job) return;
    if (stage) *stage = job->job.stage.load(std::memory_order_relaxed);
    if (iteration) *iteration = job->job.iteration.load(std::memory_order_relaxed);
    if (energy) *energy = job->job.energy.load(std::memory_order_relaxed);
    if (residual) *residual = job->job.residual.load(std::memory_order_relaxed);
}

void macqc_job_free(macqc_job* job) { delete job; }

namespace {
char* copy_string(const std::string& s) {
    char* buf = static_cast<char*>(std::malloc(s.size() + 1));
    if (!buf) return nullptr;
    std::memcpy(buf, s.c_str(), s.size() + 1);
    return buf;
}
}  // namespace

char* macqc_job_output(const macqc_job* job, long long offset, long long* total) {
    if (total) *total = 0;
    if (!job) return copy_string("");
    try {
        std::size_t t = 0;
        const std::string s = job->job.output_since(offset < 0 ? 0 : static_cast<std::size_t>(offset), &t);
        if (total) *total = static_cast<long long>(t);
        return copy_string(s);
    } catch (...) {
        return nullptr;
    }
}

char* macqc_job_progress_text(const macqc_job* job) {
    if (!job) return copy_string("");
    try {
        return copy_string(job->job.progress_text());
    } catch (...) {
        return nullptr;
    }
}

char* macqc_run_job(const char* input_json, macqc_job* job) {
    macqc::JobScope scope(job ? &job->job : nullptr);
    std::string result;
    try {
        macqc::Output out;
        try {
            const macqc::Input in = macqc::parse_input_json(input_json ? input_json : "");
            macqc::OutputSink sink;
            if (job) sink = [job](const std::string& text) { job->job.append_output(text); };
            out = macqc::run_calculation(in, sink);
        } catch (const std::exception& e) {
            out = macqc::Output();
            out.error = e.what();
        }
        result = macqc::output_to_json(out);
        macqc::report_progress(macqc::JobStage::Done);
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
