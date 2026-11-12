// ΨLab for the terminal (Linux, macOS): reads a plain-text input file, writes a
// plain-text output file.  Usage: psilab INPUT [OUTPUT]   (see psilab --help)
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "psilab_text.h"

namespace {

void usage(const char* prog) {
    std::fprintf(stdout,
                 "\xCE\xA8Lab - HeH+ Hartree-Fock and Full CI in an STO-nG basis\n\n"
                 "Usage:\n"
                 "  %s INPUT [OUTPUT]   run INPUT, write the results to OUTPUT\n"
                 "                      (default: INPUT with the extension replaced by .out)\n"
                 "  %s --example        print a commented example input file\n"
                 "  %s --help           this help\n\n"
                 "Input: one 'key = value' per line, '#' starts a comment. Keys:\n"
                 "  basis = STO-1G | STO-2G | STO-3G      R = <number> [bohr | angstrom]\n"
                 "  zeta_He = <number>   zeta_H = <number>   exact_erf = yes | no\n"
                 "  scan = yes | no      scan_step = <bohr>  (potential curve, R = 0.8 ... 5.0 bohr)\n",
                 prog, prog, prog);
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
    if (!std::strcmp(argv[1], "--example")) {
        std::cout << psilab::example_input();
        return 0;
    }
    if (argc > 3) {
        std::fprintf(stderr, "error: too many arguments (see %s --help)\n", argv[0]);
        return 2;
    }

    const std::string inputPath = argv[1];
    std::ifstream file(inputPath, std::ios::binary);
    if (!file) {
        std::fprintf(stderr, "error: cannot open %s\n", inputPath.c_str());
        return 2;
    }
    const std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    psilab::Settings settings;
    try {
        settings = psilab::parse_input_text(text);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error in %s: %s\n", inputPath.c_str(), e.what());
        return 2;
    }

    const psilab::RunResult result = psilab::run_to_text(settings);
    if (result.text.empty()) {
        std::fprintf(stderr, "error: %s\n", result.error.c_str());
        return 1;
    }

    const std::string outputPath = argc == 3 ? argv[2] : default_output_name(inputPath);
    std::ofstream out(outputPath, std::ios::binary);
    if (!out || !(out << result.text)) {
        std::fprintf(stderr, "error: cannot write %s\n", outputPath.c_str());
        return 2;
    }
    out.close();

    const auto& o = *result.output;
    std::printf("Wrote %s\n", outputPath.c_str());
    std::printf("  HF  total energy: %.8f E_h\n", o.hf->total_energy);
    if (o.fci) std::printf("  FCI total energy: %.8f E_h\n", o.fci->total_energy);
    if (!result.error.empty()) {
        std::fprintf(stderr, "warning: %s\n", result.error.c_str());
        return 1;
    }
    return 0;
}
