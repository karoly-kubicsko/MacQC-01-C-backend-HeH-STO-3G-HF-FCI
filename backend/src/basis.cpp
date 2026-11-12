#include "macqc/basis.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>

namespace macqc {

namespace {

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> split(const std::string& s) {
    std::istringstream in(s);
    std::vector<std::string> t;
    std::string w;
    while (in >> w) t.push_back(w);
    return t;
}

[[noreturn]] void fail(int line, const std::string& msg) {
    throw std::runtime_error("line " + std::to_string(line) + ": " + msg);
}

// Fortran-style real ("1.5D+01", "0.0297400", "2.0E-3").
bool parse_real(std::string t, double& value) {
    for (char& c : t)
        if (c == 'D' || c == 'd') c = 'E';
    if (t.empty()) return false;
    char* end = nullptr;
    value = std::strtod(t.c_str(), &end);
    return end && *end == '\0' && std::isfinite(value);
}

bool parse_int(const std::string& t, int& value) {
    if (t.empty()) return false;
    char* end = nullptr;
    const long v = std::strtol(t.c_str(), &end, 10);
    if (!end || *end != '\0' || v < -1000000 || v > 1000000) return false;
    value = static_cast<int>(v);
    return true;
}

int shell_l(const std::string& type) {
    static const std::string letters = "SPDFGHIK";
    if (type.size() == 1) {
        const auto k = letters.find(type[0]);
        if (k != std::string::npos) return static_cast<int>(k);
    }
    return -1;
}

double double_factorial(int n) {  // n!! for n >= -1
    double r = 1.0;
    for (int k = n; k > 1; k -= 2) r *= k;
    return r;
}

struct Spec {
    const char* name;
    const char* file;
};

const Spec kBases[] = {
    {"aug-cc-pVDZ", "aug-cc-pvdz.gbs"},   {"aug-cc-pVTZ", "aug-cc-pvtz.gbs"},
    {"aug-cc-pVQZ", "aug-cc-pvqz.gbs"},   {"aug-cc-pV5Z", "aug-cc-pv5z.gbs"},
    {"aug-cc-pCVDZ", "aug-cc-pcvdz.gbs"}, {"aug-cc-pCVTZ", "aug-cc-pcvtz.gbs"},
    {"aug-cc-pCVQZ", "aug-cc-pcvqz.gbs"}, {"aug-cc-pCV5Z", "aug-cc-pcv5z.gbs"},
};

// aug-cc-pCVXZ -> aug-cc-pVXZ (the valence set used for hydrogen), else "".
std::string valence_partner(const std::string& canonical) {
    if (canonical.rfind("aug-cc-pCV", 0) == 0) return "aug-cc-pV" + canonical.substr(10);
    return "";
}

std::string read_text_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read the basis-set file " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

}  // namespace

char angular_letter(int l) {
    static const char letters[] = "spdfghik";
    return (l >= 0 && l < 8) ? letters[l] : '?';
}

BasisFile parse_gaussian94(const std::string& text) {
    BasisFile file;
    std::istringstream in(text);
    std::string raw;
    int line = 0;
    ElementBasis* current = nullptr;
    bool in_header = true;
    while (std::getline(in, raw)) {
        ++line;
        if (line == 1 && raw.size() >= 3 && raw.compare(0, 3, "\xEF\xBB\xBF") == 0) raw.erase(0, 3);
        const std::string s = trim(raw);
        if (s.empty()) continue;
        if (s[0] == '!') {
            if (in_header) file.header += s + "\n";
            continue;
        }
        in_header = false;
        if (s.rfind("****", 0) == 0) {
            current = nullptr;
            continue;
        }
        const auto t = split(s);
        if (!current) {
            // Element line: "<symbol> 0" (a leading '-' is allowed).
            int zero = -1;
            std::string sym = t[0];
            if (!sym.empty() && sym[0] == '-') sym.erase(0, 1);
            if (t.size() != 2 || !parse_int(t[1], zero) || zero != 0)
                fail(line, "expected an element line such as 'O 0', got '" + s + "'");
            const int Z = atomic_number(sym);
            if (Z == 0) fail(line, "unknown or unsupported element '" + t[0] + "'");
            if (file.elements.count(Z)) fail(line, "element " + sym + " is defined twice");
            current = &file.elements[Z];
            current->Z = Z;
            continue;
        }
        // Shell line: "<type> <nprim> <scale>".
        int nprim = 0;
        double scale = 1.0;
        if (t.size() != 3 || !parse_int(t[1], nprim) || !parse_real(t[2], scale))
            fail(line, "expected a shell line such as 'S 3 1.00', got '" + s + "'");
        const std::string type = upper(t[0]);
        const bool sp = (type == "SP" || type == "L");
        const int l = sp ? 0 : shell_l(type);
        if (l < 0) fail(line, "unknown shell type '" + t[0] + "'");
        if (l > kMaxAngularMomentum - 1)
            fail(line, std::string("angular momentum ") + angular_letter(l) +
                           " is beyond the supported maximum (h, l = 5)");
        if (nprim < 1 || nprim > 100) fail(line, "invalid number of primitives " + t[1]);
        if (!(scale > 0.0)) fail(line, "invalid scale factor " + t[2]);
        ShellData a, b;
        a.l = l;
        b.l = 1;
        for (int k = 0; k < nprim; ++k) {
            if (!std::getline(in, raw)) fail(line, "file ends inside a shell");
            ++line;
            const auto p = split(trim(raw));
            double e = 0.0, c1 = 0.0, c2 = 0.0;
            if (p.size() != (sp ? 3u : 2u) || !parse_real(p[0], e) || !parse_real(p[1], c1) ||
                (sp && !parse_real(p[2], c2)))
                fail(line, std::string("expected ") + (sp ? "exponent and two coefficients" : "exponent and coefficient") +
                               ", got '" + trim(raw) + "'");
            if (!(e > 0.0)) fail(line, "exponents must be positive");
            e *= scale * scale;
            a.exponents.push_back(e);
            a.coefficients.push_back(c1);
            if (sp) {
                b.exponents.push_back(e);
                b.coefficients.push_back(c2);
            }
        }
        current->shells.push_back(a);
        if (sp) current->shells.push_back(b);
    }
    if (current) throw std::runtime_error("basis file ends without '****' after the last element");
    if (file.elements.empty()) throw std::runtime_error("no element found in the basis-set file");
    return file;
}

const std::vector<std::string>& supported_basis_sets() {
    static const std::vector<std::string> names = [] {
        std::vector<std::string> v;
        for (const auto& b : kBases) v.push_back(b.name);
        return v;
    }();
    return names;
}

std::string canonical_basis_name(const std::string& name) {
    const std::string n = lower(trim(name));
    for (const auto& b : kBases)
        if (n == lower(b.name)) return b.name;
    return "";
}

std::string basis_file_name(const std::string& canonical_name) {
    for (const auto& b : kBases)
        if (canonical_name == b.name) return b.file;
    return lower(canonical_name) + ".gbs";
}

std::string default_basis_directory() {
    if (const char* env = std::getenv("PSILAB_BASIS_DIR"))
        if (*env) return env;
#ifdef MACQC_BASIS_DIR
    return MACQC_BASIS_DIR;
#else
    return "basis";
#endif
}

std::vector<double> normalized_coefficients(int l, const std::vector<double>& exps,
                                            const std::vector<double>& coefs) {
    const std::size_t n = exps.size();
    std::vector<double> c(n);
    const double pi = 3.14159265358979323846;
    for (std::size_t i = 0; i < n; ++i)
        c[i] = coefs[i] * std::pow(2.0 * exps[i] / pi, 0.75) * std::pow(4.0 * exps[i], 0.5 * l) /
               std::sqrt(double_factorial(2 * l - 1));
    // <phi|phi> with x^l-normalised primitives: sum d_i d_j (2 sqrt(a_i a_j)/(a_i + a_j))^(l + 3/2)
    double norm = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            norm += coefs[i] * coefs[j] *
                    std::pow(2.0 * std::sqrt(exps[i] * exps[j]) / (exps[i] + exps[j]), l + 1.5);
    if (!(norm > 0.0)) throw std::runtime_error("a contracted shell has zero norm");
    const double scale = 1.0 / std::sqrt(norm);
    for (double& x : c) x *= scale;
    return c;
}

BasisSet build_basis(const Molecule& mol, const std::string& basis_name,
                     const std::map<std::string, BasisFile>& files) {
    const std::string name = canonical_basis_name(basis_name);
    if (name.empty()) throw std::runtime_error("unknown basis set '" + basis_name + "'");
    auto main_file = files.find(name);
    if (main_file == files.end()) throw std::runtime_error("basis set " + name + " was not loaded");

    BasisSet basis;
    basis.name = name;
    basis.source = main_file->second.header;
    std::set<std::string> extra_sources;
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) {
        const Atom& atom = mol.atoms[a];
        AtomBasisInfo info;
        info.symbol = atom.symbol;
        info.basis = name;
        const ElementBasis* eb = nullptr;
        auto it = main_file->second.elements.find(atom.Z);
        if (it != main_file->second.elements.end()) {
            eb = &it->second;
        } else {
            const std::string partner = valence_partner(name);
            if (atom.Z == 1 && !partner.empty()) {
                auto pf = files.find(partner);
                if (pf == files.end()) throw std::runtime_error("basis set " + partner + " was not loaded");
                auto pit = pf->second.elements.find(1);
                if (pit == pf->second.elements.end())
                    throw std::runtime_error(partner + " has no functions for H");
                eb = &pit->second;
                info.basis = partner;
                info.note = name + " is not defined for H (hydrogen has no core electrons); the " + partner +
                            " functions are used for H, as is customary.";
                extra_sources.insert(partner);
            } else {
                std::string msg = name + " is not available for " + atom.symbol +
                                  ": the Basis Set Exchange does not define it for this element.";
                if (name == "aug-cc-pCV5Z" && atom.Z == 3)
                    msg += " (The BSE has aug-cc-pCV5Z only for B-Ar.) For " + mol.name +
                           " choose aug-cc-pCVQZ or aug-cc-pV5Z.";
                throw std::runtime_error(msg);
            }
        }
        int count[kMaxAngularMomentum] = {0};
        std::set<double> prims[kMaxAngularMomentum];
        for (const ShellData& sd : eb->shells) {
            Shell sh;
            sh.l = sd.l;
            sh.atom = static_cast<int>(a);
            sh.center = atom.position;
            sh.exponents = sd.exponents;
            sh.coefficients = normalized_coefficients(sd.l, sd.exponents, sd.coefficients);
            sh.first_function = basis.n_functions;
            basis.n_functions += sh.size();
            basis.n_primitives += static_cast<int>(sd.exponents.size()) * sh.size();
            basis.max_l = std::max(basis.max_l, sd.l);
            info.functions += sh.size();
            ++count[sd.l];
            for (double e : sd.exponents) prims[sd.l].insert(e);
            basis.shells.push_back(std::move(sh));
        }
        std::string con = "[", pri = "(";
        for (int l = 0; l < kMaxAngularMomentum; ++l) {
            if (count[l]) con += std::to_string(count[l]) + angular_letter(l);
            if (!prims[l].empty()) pri += std::to_string(prims[l].size()) + angular_letter(l);
        }
        info.contraction = con + "]";
        info.primitives = pri + ")";
        basis.atoms.push_back(info);
    }
    for (const auto& s : extra_sources) basis.source += files.at(s).header;
    return basis;
}

BasisSet build_basis(const Molecule& mol, const std::string& basis_name, const std::string& directory) {
    const std::string name = canonical_basis_name(basis_name);
    if (name.empty()) {
        std::string list;
        for (const auto& b : supported_basis_sets()) list += (list.empty() ? "" : ", ") + b;
        throw std::runtime_error("unknown basis set '" + basis_name + "' (available: " + list + ")");
    }
    std::string dir = directory.empty() ? default_basis_directory() : directory;
    if (!dir.empty() && dir.back() != '/') dir += '/';
    std::map<std::string, BasisFile> files;
    auto load = [&](const std::string& n) {
        const std::string path = dir + basis_file_name(n);
        std::string text;
        try {
            text = read_text_file(path);
        } catch (const std::exception&) {
            throw std::runtime_error("cannot read the basis-set file " + path +
                                     " (set the basis directory, e.g. PSILAB_BASIS_DIR)");
        }
        try {
            files[n] = parse_gaussian94(text);
        } catch (const std::exception& e) {
            throw std::runtime_error(path + ", " + e.what());
        }
    };
    load(name);
    const std::string partner = valence_partner(name);
    bool need_partner = false;
    for (const Atom& at : mol.atoms)
        if (at.Z == 1 && !files[name].elements.count(1)) need_partner = true;
    if (need_partner && !partner.empty()) load(partner);
    return build_basis(mol, name, files);
}

}  // namespace macqc
