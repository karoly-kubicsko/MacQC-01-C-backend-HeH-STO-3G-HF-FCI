// MacQC command-line tool.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "macqc/calculation.h"
#include "macqc/macqc_api.h"

namespace {

void usage(const char* prog) {
    std::fprintf(stderr,
                 "Usage: %s [options]\n"
                 "  Restricted HF + Full CI for HeH+ with an STO-nG (n = 1..3) s basis.\n\n"
                 "  --R <bohr>            bond length            (default 1.4632)\n"
                 "  --zeta1 <z>           Slater exponent, atom A (default 2.0925)\n"
                 "  --zeta2 <z>           Slater exponent, atom B (default 1.24)\n"
                 "  --ZA <z> --ZB <z>     nuclear charges         (default 2, 1)\n"
                 "  --n <1|2|3>           number of Gaussians     (default 3)\n"
                 "  --exact-erf           use std::erf instead of the Abramowitz-Stegun erf\n"
                 "  --max-iter <n>        maximum SCF iterations  (default 25)\n"
                 "  --conv <x>            density convergence     (default 1e-4)\n"
                 "  --print <0|1|2>       report detail           (default 2)\n"
                 "  --json                print the JSON result instead of the text report\n"
                 "  --input <file|->      read the input JSON (the API format); flags override it\n"
                 "  -h, --help            this help\n",
                 prog);
}

}  // namespace

int main(int argc, char** argv) {
    macqc::Input in;
    bool json = false;
    std::string input_file;

    // Pass 1: find --input so explicit flags can override it.
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], "--input") && i + 1 < argc) input_file = argv[i + 1];
    if (!input_file.empty()) {
        std::string text;
        if (input_file == "-") {
            text.assign(std::istreambuf_iterator<char>(std::cin), std::istreambuf_iterator<char>());
        } else {
            std::ifstream f(input_file);
            if (!f) {
                std::fprintf(stderr, "error: cannot open %s\n", input_file.c_str());
                return 2;
            }
            text.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
        }
        try {
            in = macqc::parse_input_json(text);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            return 2;
        }
    }

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s needs a value\n", name);
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (a == "--json") json = true;
        else if (a == "--exact-erf") in.exactErf = true;
        else if (a == "--input") next("--input");
        else if (a == "--R") in.R = std::atof(next("--R"));
        else if (a == "--zeta1") in.zeta1 = std::atof(next("--zeta1"));
        else if (a == "--zeta2") in.zeta2 = std::atof(next("--zeta2"));
        else if (a == "--ZA") in.ZA = std::atof(next("--ZA"));
        else if (a == "--ZB") in.ZB = std::atof(next("--ZB"));
        else if (a == "--n") in.nGaussians = std::atoi(next("--n"));
        else if (a == "--max-iter") in.maxIter = std::atoi(next("--max-iter"));
        else if (a == "--conv") in.convergence = std::atof(next("--conv"));
        else if (a == "--print") in.printLevel = std::atoi(next("--print"));
        else {
            std::fprintf(stderr, "error: unknown option %s\n", a.c_str());
            usage(argv[0]);
            return 2;
        }
    }

    const macqc::Output out = macqc::run_calculation(in);
    if (json) {
        std::cout << macqc::output_to_json(out, 2) << "\n";
    } else {
        std::cout << out.report;
        if (!out.error.empty()) std::fprintf(stderr, "error: %s\n", out.error.c_str());
    }
    return out.error.empty() ? 0 : 1;
}
