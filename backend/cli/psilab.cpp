// ΨLab for the terminal (Linux, macOS): reads a plain-text input file, writes a
// plain-text output file.  Usage: psilab INPUT [OUTPUT]   (see psilab --help)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "psilab_text.h"

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif
#include <sys/stat.h>
#include <unistd.h>

namespace {

void usage(const char* prog) {
    std::fprintf(stdout,
                 "\xCE\xA8Lab - Hartree-Fock and Full CI in the terminal\n\n"
                 "Usage:\n"
                 "  %s INPUT [OUTPUT]   run INPUT, write the results to OUTPUT\n"
                 "                      (default: INPUT with the extension replaced by .out)\n"
                 "  %s --estimate INPUT size of the calculation and time estimate, no calculation\n"
                 "  %s --example        print a commented example input file (HeH+)\n"
                 "  %s --version        print the version\n"
                 "  %s --help           this help\n\n"
                 "Input: one 'key = value' per line, '#' starts a comment. Keys:\n"
                 "  molecule = HeH+ | H2 | LiH | H2O      (default HeH+)\n"
                 "  R = <number> [bohr | angstrom]       bond length (H2O: O-H)\n"
                 "HeH+ (Hartree-Fock and Full CI, Szabo-Ostlund):\n"
                 "  basis = STO-1G | STO-2G | STO-3G      zeta_He = <number>   zeta_H = <number>\n"
                 "  exact_erf = yes | no   scan = yes | no   scan_step = <bohr>  (R = 0.8 ... 5.0 bohr)\n"
                 "H2, LiH, H2O (restricted Hartree-Fock, coupled cluster, Full CI for H2 and LiH):\n"
                 "  basis = aug-cc-pVDZ | aug-cc-pVTZ | aug-cc-pVQZ | aug-cc-pV5Z |\n"
                 "          aug-cc-pCVDZ | aug-cc-pCVTZ | aug-cc-pCVQZ | aug-cc-pCV5Z\n"
                 "  angle = <degrees> (H2O)   threads = <n>\n"
                 "  method = HF | CCSD | CCSD(T) | CCSDT | CCSDTQ | FCI   (all electrons correlated)\n"
                 "  memory = <MB> (two-electron integrals kept in memory up to this size)\n"
                 "  fci_memory = <MB> (memory limit of Full CI and coupled cluster, default 4096)\n"
                 "  basis_dir = <folder with the .gbs files> (or environment variable PSILAB_BASIS_DIR)\n"
                 "  job = energy | opt | freq | opt+freq | opt+freq+sp         (default energy)\n"
                 "        opt: geometry optimisation (redundant internal coordinates, RFO, BFGS)\n"
                 "        freq: harmonic frequencies, IR intensities and thermochemistry\n"
                 "        opt+freq+sp: then a single-point energy with sp_method / sp_basis\n"
                 "  sp_method = <method>   sp_basis = <basis>    (opt+freq+sp; default: method, basis)\n"
                 "  temperature = <K> (298.15)   pressure = <number> [Pa | kPa | bar | atm] (101325 Pa)\n"
                 "  opt_convergence = normal | tight   max_steps = <n> (50)\n\n"
                 "The output file is written while the job runs (tail -f OUTPUT follows it).\n",
                 prog, prog, prog, prog, prog);
}

bool is_directory(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

// <folder of the executable>/../share/psilab/basis, if it exists (installed psilab).
std::string installed_basis_directory() {
    std::string exe;
#if defined(__APPLE__)
    char buf[4096];
    uint32_t size = sizeof buf;
    if (_NSGetExecutablePath(buf, &size) == 0) exe = buf;
#elif defined(__linux__)
    char buf[4096];
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof buf - 1);
    if (n > 0) exe.assign(buf, static_cast<std::size_t>(n));
#endif
    const auto slash = exe.find_last_of('/');
    if (slash == std::string::npos) return "";
    const std::string dir = exe.substr(0, slash) + "/../share/psilab/basis";
    return is_directory(dir) ? dir : "";
}

std::string read_input(const std::string& path, psilab::Settings& settings, int& status) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        std::fprintf(stderr, "error: cannot open %s\n", path.c_str());
        status = 2;
        return "";
    }
    const std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    try {
        settings = psilab::parse_input_text(text);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error in %s: %s\n", path.c_str(), e.what());
        status = 2;
        return "";
    }
    // Basis-set folder: input file > PSILAB_BASIS_DIR > installed copy > source tree.
    const char* env = std::getenv("PSILAB_BASIS_DIR");
    if (settings.input.basisDirectory.empty() && !(env && *env))
        settings.input.basisDirectory = installed_basis_directory();
    status = 0;
    return text;
}

std::string default_output_name(const std::string& input) {
    const auto slash = input.find_last_of('/');
    const auto dot = input.find_last_of('.');
    if (dot != std::string::npos && (slash == std::string::npos || dot > slash) && dot != 0 &&
        input.substr(dot) != ".out")
        return input.substr(0, dot) + ".out";
    return input + ".out";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
        usage(argv[0]);
        return argc < 2 ? 2 : 0;
    }
    if (!std::strcmp(argv[1], "--version")) {
        std::printf("\xCE\xA8Lab %s\n", psilab::kVersion);
        return 0;
    }
    if (!std::strcmp(argv[1], "--example")) {
        std::cout << psilab::example_input();
        return 0;
    }
    if (!std::strcmp(argv[1], "--estimate")) {
        if (argc != 3) {
            std::fprintf(stderr, "usage: %s --estimate INPUT\n", argv[0]);
            return 2;
        }
        psilab::Settings settings;
        int status = 0;
        read_input(argv[2], settings, status);
        if (status) return status;
        const std::string t = psilab::estimate_text(settings);
        if (t.empty()) {
            std::printf("HeH+ (STO-nG): a two-function calculation, it takes milliseconds.\n");
            return 0;
        }
        std::fputs(t.c_str(), stdout);
        return t.rfind("Error: ", 0) == 0 ? 1 : 0;
    }
    if (argc > 3) {
        std::fprintf(stderr, "error: too many arguments (see %s --help)\n", argv[0]);
        return 2;
    }

    const std::string inputPath = argv[1];
    psilab::Settings settings;
    int status = 0;
    read_input(inputPath, settings, status);
    if (status) return status;
    for (const auto& w : psilab::estimate_warnings(settings)) std::fprintf(stderr, "warning: %s\n", w.c_str());

    // The output file is opened at the first piece of output and flushed after every piece, so
    // `tail -f` follows the job (H2, LiH, H2O); HeH+ writes it at the end.
    const std::string outputPath = argc == 3 ? argv[2] : default_output_name(inputPath);
    std::ofstream out;
    bool write_failed = false;
    auto sink = [&](const std::string& text) {
        if (write_failed) return;
        if (!out.is_open()) {
            out.open(outputPath, std::ios::binary | std::ios::trunc);
            if (!out) {
                write_failed = true;
                std::fprintf(stderr, "error: cannot write %s\n", outputPath.c_str());
                return;
            }
        }
        out << text;
        out.flush();
        if (!out) write_failed = true;
    };
    const bool streamed = !macqc::is_heh_input(settings.input);
    const psilab::RunResult result = psilab::run_to_text(settings, streamed ? macqc::OutputSink(sink) : macqc::OutputSink());
    if (result.text.empty()) {
        std::fprintf(stderr, "error: %s\n", result.error.c_str());
        return 1;
    }
    if (!streamed) sink(result.text);
    if (write_failed) {
        std::fprintf(stderr, "error: cannot write %s\n", outputPath.c_str());
        return 2;
    }
    out.close();

    const auto& o = *result.output;
    std::printf("Wrote %s\n", outputPath.c_str());
    if (o.job && o.input.job != "energy") {
        const auto& j = *o.job;
        if (j.optimization)
            std::printf("  optimisation: %s in %d steps, E = %.10f E_h\n",
                        j.optimization->converged ? "converged" : "NOT converged", j.optimization->steps,
                        j.optimization->final_energy);
        if (j.harmonic) {
            std::printf("  frequencies (cm^-1):");
            for (const auto& m : j.harmonic->modes) std::printf(" %.2f%s", std::fabs(m.wavenumber), m.imaginary ? "i" : "");
            std::printf("\n");
        }
        if (j.thermo)
            std::printf("  G(%.2f K) = %.10f E_h, H = %.10f E_h, ZPE = %.10f E_h\n", j.thermo->temperature,
                        j.thermo->total_g, j.thermo->total_h, j.thermo->total_zpe - j.thermo->electronic_energy);
        if (!j.final_energy_label.empty())
            std::printf("  final energy %s: %.10f E_h\n", j.final_energy_label.c_str(), j.final_energy);
    } else if (o.molecular && o.molecular->rhf) {
        std::printf("  RHF total energy: %.10f E_h\n", o.molecular->rhf->total_energy);
        if (o.molecular->cc) std::printf("  %s total energy: %.10f E_h\n", o.molecular->cc->method.c_str(), o.molecular->cc->total_energy);
        if (o.molecular->fci) std::printf("  FCI total energy: %.10f E_h\n", o.molecular->fci->total_energy);
    } else if (o.hf) {
        std::printf("  HF  total energy: %.8f E_h\n", o.hf->total_energy);
        if (o.fci) std::printf("  FCI total energy: %.8f E_h\n", o.fci->total_energy);
    }
    if (!result.error.empty()) {
        const bool failed = o.job && o.job->termination != "normal";
        std::fprintf(stderr, "%s: %s\n", failed ? "error" : "warning", result.error.c_str());
        return 1;
    }
    return 0;
}
