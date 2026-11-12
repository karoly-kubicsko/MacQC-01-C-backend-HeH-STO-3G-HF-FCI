#include "psilab_text.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
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

}  // namespace

Settings parse_input_text(const std::string& text) {
    Settings s;
    std::set<std::string> seen;
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
        const std::string k = normalize_key(key);
        if (!seen.insert(k).second) fail(line, "'" + key + "' is given twice");

        if (k == "basis" || k == "basisset") {
            s.input.nGaussians = parse_basis(value, line);
        } else if (k == "r" || k == "bondlength") {
            s.input.R = parse_length(value, line);
        } else if (k == "zetahe" || k == "zeta1") {
            s.input.zeta1 = parse_number(value, line, key);
        } else if (k == "zetah" || k == "zeta2") {
            s.input.zeta2 = parse_number(value, line, key);
        } else if (k == "exacterf") {
            s.input.exactErf = parse_bool(value, line, key);
        } else if (k == "scan") {
            s.scan = parse_bool(value, line, key);
        } else if (k == "scanstep") {
            s.scanStep = parse_number(value, line, key);
            if (!(s.scanStep > 0.0) || s.scanStep > kScanEnd - kScanStart)
                fail(line, format("scan_step must be > 0 and <= %.1f bohr", kScanEnd - kScanStart));
        } else {
            fail(line, "unknown key '" + key +
                           "' (known: basis, R, zeta_He, zeta_H, exact_erf, scan, scan_step)");
        }
    }
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

RunResult run_to_text(const Settings& settings) {
    RunResult result;
    macqc::Input in = settings.input;
    in.printLevel = 2;
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
