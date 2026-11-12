#include "psilab_text.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace psilab {

namespace {

std::string format(const char* fmt, ...) {
    char buf[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    return buf;
}

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Keys compare case-insensitively, ignoring '_', '-' and '(' ')' spacing variants.
std::string normalize_key(const std::string& key) {
    std::string out;
    for (char c : lower(key))
        if (c != '_' && c != '-' && c != ' ' && c != '\t') out += c;
    return out;
}

[[noreturn]] void fail(int line, const std::string& message) {
    throw std::runtime_error("line " + std::to_string(line) + ": " + message);
}

double parse_number(const std::string& text, int line, const std::string& key) {
    const std::string t = trim(text);
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (t.empty() || end == t.c_str() || *end != '\0' || !std::isfinite(v))
        fail(line, "'" + key + "' needs a number, got '" + t + "'");
    return v;
}

bool parse_bool(const std::string& text, int line, const std::string& key) {
    const std::string t = lower(trim(text));
    if (t == "yes" || t == "true" || t == "on" || t == "1") return true;
    if (t == "no" || t == "false" || t == "off" || t == "0") return false;
    fail(line, "'" + key + "' must be yes or no, got '" + trim(text) + "'");
}

int parse_basis(const std::string& text, int line) {
    std::string t;
    for (char c : lower(trim(text)))
        if (c != '-' && c != ' ') t += c;
    if (t == "sto1g" || t == "1") return 1;
    if (t == "sto2g" || t == "2") return 2;
    if (t == "sto3g" || t == "3") return 3;
    fail(line, "basis must be STO-1G, STO-2G or STO-3G, got '" + trim(text) + "'");
}

// "1.4632", "1.4632 bohr", "0.7743 angstrom", "0.7743 Å" -> bohr
double parse_length(const std::string& text, int line) {
    std::istringstream in(trim(text));
    std::string number, unit, extra;
    in >> number >> unit >> extra;
    if (!extra.empty()) fail(line, "R has unexpected text '" + extra + "'");
    const double v = parse_number(number, line, "R");
    const std::string u = lower(unit);
    if (u.empty() || u == "bohr" || u == "a0" || u == "au") return v;
    if (u == "angstrom" || u == "ang" || u == "a" || unit == "\xC3\x85" /* Å */ ||
        unit == "\xE2\x84\xAB" /* Å (angstrom sign) */)
        return v / kBohrToAngstrom;
    fail(line, "unknown length unit '" + unit + "' (use bohr or angstrom)");
}

std::string basis_name(int n) { return "STO-" + std::to_string(n) + "G"; }

std::string length_both(double bohr, int decimals) {
    return format("%.*f \xC3\x85 (%.*f bohr)", decimals, bohr * kBohrToAngstrom, decimals, bohr);
}

std::string summary_text(const macqc::Output& out) {
    const macqc::Input& in = out.input;
    std::vector<std::string> lines;
    lines.push_back("\xCE\xA8Lab \xE2\x80\x93 HeH+ (" + basis_name(in.nGaussians) + ")");
    lines.push_back("R = " + length_both(in.R, 4) +
                    format(", zeta(He) = %.4f, zeta(H) = %.4f, %s erf", in.zeta1, in.zeta2,
                           in.exactErf ? "exact" : "Abramowitz-Stegun"));
    if (out.hf) {
        const auto& hf = *out.hf;
        lines.push_back("");
        lines.push_back("Hartree-Fock");
        lines.push_back(format("  converged:            %s (%d iterations)", hf.converged ? "yes" : "no",
                               hf.iterations));
        lines.push_back(format("  electronic energy:    %.8f E_h", hf.electronic_energy));
        lines.push_back(format("  nuclear repulsion:    %.8f E_h", hf.nuclear_repulsion));
        lines.push_back(format("  total energy:         %.8f E_h", hf.total_energy));
        if (hf.orbital_energies.size() >= 2) {
            std::string e;
            for (std::size_t i = 0; i < hf.orbital_energies.size(); ++i)
                e += (i ? ", " : "") + format("%.8f", hf.orbital_energies[i]);
            lines.push_back("  orbital energies:     " + e + " E_h");
        }
    }
    if (out.fci) {
        const auto& fci = *out.fci;
        lines.push_back("");
        lines.push_back("Full CI");
        lines.push_back(format("  total energy:         %.8f E_h", fci.total_energy));
        lines.push_back(format("  correlation energy:   %.8f E_h", fci.correlation_energy));
        for (std::size_t i = 0; i < fci.states.size(); ++i)
            lines.push_back(format("  state %zu: %.8f E_h  <S^2> = %.4f  %s", i + 1, fci.states[i].energy,
                                   fci.states[i].S2, fci.states[i].spin.c_str()));
        for (const auto& w : fci.singlet_weights)
            lines.push_back(format("  weight %s: %.6f", w.configuration.c_str(), w.weight));
    }
    if (!out.error.empty()) {
        lines.push_back("");
        lines.push_back("Error: " + out.error);
    }
    std::string s;
    for (std::size_t i = 0; i < lines.size(); ++i) s += (i ? "\n" : "") + lines[i];
    return s;
}

std::string general_summary(const macqc::Output& out) {
    const macqc::MolecularResult& m = *out.molecular;
    std::vector<std::string> lines;
    lines.push_back("\xCE\xA8Lab \xE2\x80\x93 " + m.molecule.name + ", RHF/" + m.basis.name);
    if (m.molecule.name == "H2O")
        lines.push_back("R(O-H) = " + length_both(m.input.R, 4) + format(", angle H-O-H = %.2f deg", m.input.angle));
    else
        lines.push_back(std::string(m.molecule.name == "LiH" ? "R(Li-H) = " : "R(H-H) = ") + length_both(m.input.R, 4));
    std::string per_atom;
    std::vector<std::string> seen;
    for (const auto& a : m.basis.atoms) {
        if (std::find(seen.begin(), seen.end(), a.symbol) != seen.end()) continue;
        seen.push_back(a.symbol);
        per_atom += (per_atom.empty() ? "" : ", ") + a.symbol + " " + a.contraction;
    }
    lines.push_back(format("%d basis functions (%s)", m.basis.n_functions, per_atom.c_str()));
    for (const auto& n : m.notes) lines.push_back("Note: " + n);
    if (m.rhf) {
        const auto& h = *m.rhf;
        lines.push_back("");
        lines.push_back("Hartree-Fock (RHF)");
        lines.push_back(format("  converged:            %s (%d iterations)", h.converged ? "yes" : "no", h.iterations));
        lines.push_back(format("  total energy:         %.10f E_h", h.total_energy));
        lines.push_back(format("  electronic energy:    %.10f E_h", h.electronic_energy));
        lines.push_back(format("  nuclear repulsion:    %.10f E_h", h.nuclear_repulsion));
        const int no = h.n_occupied;
        if (no >= 1 && no < static_cast<int>(h.orbital_energies.size()))
            lines.push_back(format("  HOMO, LUMO:           %.8f, %.8f E_h", h.orbital_energies[no - 1],
                                   h.orbital_energies[no]));
    }
    for (const auto& w : m.warnings) lines.push_back("Warning: " + w);
    if (!out.error.empty()) {
        lines.push_back("");
        lines.push_back("Error: " + out.error);
    }
    std::string s;
    for (std::size_t i = 0; i < lines.size(); ++i) s += (i ? "\n" : "") + lines[i];
    return s;
}

}  // namespace

Settings parse_input_text(const std::string& text) {
    struct Entry {
        std::string key, value;
        int line = 0;
    };
    // Canonical key for every accepted spelling.
    auto canonical_key = [](const std::string& k) -> std::string {
        if (k == "molecule") return "molecule";
        if (k == "basis" || k == "basisset") return "basis";
        if (k == "r" || k == "bondlength") return "R";
        if (k == "angle" || k == "bondangle" || k == "hoh") return "angle";
        if (k == "method") return "method";
        if (k == "zetahe" || k == "zeta1") return "zeta_He";
        if (k == "zetah" || k == "zeta2") return "zeta_H";
        if (k == "exacterf") return "exact_erf";
        if (k == "scan") return "scan";
        if (k == "scanstep") return "scan_step";
        if (k == "threads") return "threads";
        if (k == "memory" || k == "memorymb") return "memory";
        if (k == "basisdir" || k == "basisdirectory") return "basis_dir";
        return "";
    };
    std::map<std::string, Entry> entries;
    std::istringstream in(text);
    std::string raw;
    int line = 0;
    while (std::getline(in, raw)) {
        ++line;
        if (line == 1 && raw.size() >= 3 && raw.compare(0, 3, "\xEF\xBB\xBF") == 0) raw.erase(0, 3);  // BOM
        const auto hash = raw.find('#');
        const std::string content = trim(hash == std::string::npos ? raw : raw.substr(0, hash));
        if (content.empty()) continue;
        const auto eq = content.find_first_of("=:");
        if (eq == std::string::npos) fail(line, "expected 'key = value', got '" + content + "'");
        const std::string key = trim(content.substr(0, eq));
        const std::string value = trim(content.substr(eq + 1));
        if (value.empty()) fail(line, "'" + key + "' has no value");
        const std::string k = canonical_key(normalize_key(key));
        if (k.empty())
            fail(line, "unknown key '" + key +
                           "' (known: molecule, basis, R, angle, method, zeta_He, zeta_H, exact_erf, scan, "
                           "scan_step, threads, memory, basis_dir)");
        if (entries.count(k)) fail(line, "'" + key + "' is given twice");
        entries[k] = Entry{key, value, line};
    }

    Settings s;
    std::string molecule = "HeH+";
    if (entries.count("molecule")) {
        const Entry& e = entries["molecule"];
        molecule = macqc::canonical_molecule_name(e.value);
        if (molecule.empty()) fail(e.line, "molecule must be HeH+, H2, LiH or H2O, got '" + e.value + "'");
    }
    s.input.molecule = molecule;
    const bool heh = molecule == "HeH+";
    auto only = [&](const char* key, bool allowed, const std::string& who) {
        auto it = entries.find(key);
        if (it != entries.end() && !allowed)
            fail(it->second.line, "'" + it->second.key + "' applies to " + who + " only (molecule = " + molecule + ")");
    };
    only("zeta_He", heh, "HeH+");
    only("zeta_H", heh, "HeH+");
    only("exact_erf", heh, "HeH+");
    only("scan", heh, "HeH+");
    only("scan_step", heh, "HeH+");
    only("angle", molecule == "H2O", "H2O");
    only("method", !heh, "H2, LiH and H2O (HeH+ always runs Hartree-Fock and Full CI)");
    only("threads", !heh, "H2, LiH and H2O");
    only("memory", !heh, "H2, LiH and H2O");
    only("basis_dir", !heh, "H2, LiH and H2O");

    std::vector<std::pair<std::string, Entry>> ordered(entries.begin(), entries.end());
    std::sort(ordered.begin(), ordered.end(),
              [](const auto& a, const auto& b) { return a.second.line < b.second.line; });
    for (const auto& [k, e] : ordered) {
        const int ln = e.line;
        const std::string& value = e.value;
        if (k == "basis") {
            if (heh) {
                s.input.nGaussians = parse_basis(value, ln);
            } else {
                s.input.basis = macqc::canonical_basis_name(value);
                if (s.input.basis.empty()) {
                    std::string list;
                    for (const auto& b : macqc::supported_basis_sets()) list += (list.empty() ? "" : ", ") + b;
                    fail(ln, "basis must be one of " + list + ", got '" + value + "'");
                }
            }
        } else if (k == "R") {
            s.input.R = parse_length(value, ln);
            s.rGiven = true;
        } else if (k == "angle") {
            s.input.angle = parse_number(value, ln, e.key);
        } else if (k == "method") {
            std::string m;
            for (char c : lower(value))
                if (c != '-' && c != ' ') m += c;
            if (m == "rhf" || m == "hf" || m == "hartreefock") {
                s.input.method = "RHF";
            } else if (m == "fci" || m == "fullci") {
                fail(ln, "Full CI for H2, LiH and H2O is not available yet; use method = RHF");
            } else {
                fail(ln, "method must be RHF, got '" + value + "'");
            }
        } else if (k == "zeta_He") {
            s.input.zeta1 = parse_number(value, ln, e.key);
        } else if (k == "zeta_H") {
            s.input.zeta2 = parse_number(value, ln, e.key);
        } else if (k == "exact_erf") {
            s.input.exactErf = parse_bool(value, ln, e.key);
        } else if (k == "scan") {
            s.scan = parse_bool(value, ln, e.key);
        } else if (k == "scan_step") {
            s.scanStep = parse_number(value, ln, e.key);
            if (!(s.scanStep > 0.0) || s.scanStep > kScanEnd - kScanStart)
                fail(ln, format("scan_step must be > 0 and <= %.1f bohr", kScanEnd - kScanStart));
        } else if (k == "threads") {
            const double t = parse_number(value, ln, e.key);
            if (t != std::floor(t) || t < 0 || t > 64) fail(ln, "threads must be an integer from 0 (all) to 64");
            s.input.threads = static_cast<int>(t);
        } else if (k == "memory") {
            std::istringstream ms(value);
            std::string number, unit, extra;
            ms >> number >> unit >> extra;
            if (!extra.empty()) fail(ln, "memory has unexpected text '" + extra + "'");
            double mb = parse_number(number, ln, e.key);
            const std::string u = lower(unit);
            if (u == "gb" || u == "g") mb *= 1024.0;
            else if (!u.empty() && u != "mb" && u != "m") fail(ln, "memory unit must be MB or GB");
            if (mb < 0.0) fail(ln, "memory must not be negative");
            s.input.memoryMB = mb;
        } else if (k == "basis_dir") {
            s.input.basisDirectory = value;
        }
    }
    if (!heh && !s.rGiven) s.input.R = macqc::default_bond_length(molecule);
    const std::string problem = macqc::validate(s.input);
    if (!problem.empty()) throw std::runtime_error(problem);
    return s;
}

std::string example_input() {
    return
        "# \xCE\xA8Lab input file: HeH+ with Hartree-Fock and Full CI\n"
        "# One 'key = value' per line; '#' starts a comment.\n"
        "\n"
        "basis     = STO-3G           # STO-1G, STO-2G or STO-3G\n"
        "R         = 1.4632 bohr      # bond length; unit bohr or angstrom (e.g. 0.7743 angstrom)\n"
        "zeta_He   = 2.0925           # Slater exponent of He\n"
        "zeta_H    = 1.24             # Slater exponent of H\n"
        "exact_erf = no               # no: the original program's erf; yes: the exact erf\n"
        "\n"
        "# Potential curve, R = 0.8 ... 5.0 bohr (uses basis, zetas and exact_erf from above)\n"
        "scan      = yes\n"
        "scan_step = 0.1              # bohr\n";
}

std::vector<double> scan_grid(double step) {
    std::vector<double> values;
    if (!(step > 0.0)) return values;
    const int count = static_cast<int>(std::floor((kScanEnd - kScanStart) / step + 1.0e-9));
    for (int i = 0; i <= count; ++i) values.push_back(kScanStart + i * step);
    if (!values.empty() && values.back() < kScanEnd - 1.0e-9) values.push_back(kScanEnd);
    return values;
}

std::optional<CurveMinimum> curve_minimum(const std::vector<std::pair<double, double>>& samples,
                                          const std::string& method) {
    if (samples.empty()) return std::nullopt;
    std::size_t best = 0;
    for (std::size_t i = 1; i < samples.size(); ++i)
        if (samples[i].second < samples[best].second) best = i;
    CurveMinimum m{method, samples[best].first, samples[best].second, false};
    if (best > 0 && best + 1 < samples.size()) {
        const auto [r0, e0] = samples[best - 1];
        const auto [r1, e1] = samples[best];
        const auto [r2, e2] = samples[best + 1];
        const double h1 = r1 - r0, h2 = r2 - r1;
        const double curvature = e0 - 2.0 * e1 + e2;
        if (std::fabs(h1 - h2) < 1.0e-9 && curvature > 0.0) {
            const double shift = -0.5 * h1 * (e2 - e0) / curvature;
            const double energy = e1 - (e2 - e0) * (e2 - e0) / (8.0 * curvature);
            if (std::fabs(shift) <= h1 && std::isfinite(energy)) {
                m.R = r1 + shift;
                m.energy = energy;
                m.refined = true;
            }
        }
    }
    return m;
}

std::string estimate_text(const Settings& settings) {
    if (macqc::is_heh_input(settings.input)) return "";
    macqc::Input in = settings.input;
    in.estimateOnly = true;
    in.printLevel = 1;
    const macqc::Output out = macqc::run_calculation(in);
    if (!out.molecular) return "Error: " + out.error + "\n";
    const auto& e = out.molecular->estimate;
    std::string t = out.report;
    t += format("Estimated wall time: %.0f s on %d threads (rough estimate)\n", e.rhf_seconds, e.threads);
    return t;
}

std::vector<std::string> estimate_warnings(const Settings& settings) {
    if (macqc::is_heh_input(settings.input)) return {};
    macqc::Input in = settings.input;
    in.estimateOnly = true;
    in.printLevel = 0;
    const macqc::Output out = macqc::run_calculation(in);
    if (!out.molecular) return {};
    return out.molecular->warnings;
}

RunResult run_to_text(const Settings& settings) {
    RunResult result;
    macqc::Input in = settings.input;
    in.printLevel = 2;
    if (!macqc::is_heh_input(in)) {
        macqc::Output out = macqc::run_calculation(in);
        result.output = out;
        if (!out.molecular || !out.molecular->rhf) {
            result.error = out.error.empty() ? "no result" : out.error;
            return result;
        }
        result.error = out.error;
        const std::string rule(78, '=');
        result.text = general_summary(out) + "\n\n" + rule + "\nDETAILED REPORT\n" + rule + "\n" + out.report;
        if (result.text.back() != '\n') result.text += "\n";
        return result;
    }
    macqc::Output out = macqc::run_calculation(in);
    result.output = out;
    if (!out.hf) {
        result.error = out.error.empty() ? "no result" : out.error;
        return result;
    }
    result.error = out.error;  // e.g. SCF not converged: still write what exists

    const std::string rule(78, '=');
    std::string text = summary_text(out);
    if (!out.report.empty())
        text += "\n\n" + rule +
                "\nDETAILED REPORT (Hartree-Fock as in Szabo & Ostlund, Appendix B; then Full CI)\n" +
                rule + "\n" + out.report;

    if (settings.scan) {
        std::vector<CurvePoint> points;
        std::vector<std::string> failures;
        macqc::Input point = settings.input;
        point.printLevel = 0;
        for (double r : scan_grid(settings.scanStep)) {
            point.R = r;
            const macqc::Output o = macqc::run_calculation(point);
            if (o.hf && o.fci && o.error.empty())
                points.push_back({r, o.hf->total_energy, o.fci->total_energy});
            else
                failures.push_back(format("R = %.3f bohr: %s", r, o.error.empty() ? "no result" : o.error.c_str()));
        }
        text += "\n\n" + rule + "\nPOTENTIAL CURVE\n";
        text += basis_name(settings.input.nGaussians) +
                format(", zeta(He) = %.4f, zeta(H) = %.4f, %s erf, step %.4f bohr\n\n", settings.input.zeta1,
                       settings.input.zeta2, settings.input.exactErf ? "exact" : "A&S", settings.scanStep);
        text += "    R (bohr)       R (\xC3\x85)      E_HF (E_h)     E_FCI (E_h)\n";
        std::vector<std::pair<double, double>> hf, fci;
        for (const auto& p : points) {
            text += format("%12.4f  %12.4f  %14.8f  %14.8f\n", p.R, p.R * kBohrToAngstrom, p.hf, p.fci);
            hf.emplace_back(p.R, p.hf);
            fci.emplace_back(p.R, p.fci);
        }
        text += "\n";
        for (const auto& m : {curve_minimum(hf, "HF"), curve_minimum(fci, "FCI")})
            if (m)
                text += m->method + " minimum (" + (m->refined ? "parabolic fit" : "grid point") +
                        "): R = " + length_both(m->R, 4) + format(", E = %.8f E_h\n", m->energy);
        for (const auto& f : failures) text += "failed point: " + f + "\n";
    }
    if (text.empty() || text.back() != '\n') text += "\n";
    result.text = text;
    return result;
}

}  // namespace psilab
