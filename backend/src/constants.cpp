#include "macqc/constants.h"

#include <stdexcept>

namespace macqc {

namespace {
// AME 2020 atomic masses (u) of the most abundant isotopes.
const double kMass[] = {1.007825031898, 4.002603254130, 7.016003434, 9.012183062, 11.009305167,
                        12.0,           14.003074004251, 15.994914619257, 18.998403162067, 19.992440176};
const char* const kLabel[] = {"1H", "4He", "7Li", "9Be", "11B", "12C", "14N", "16O", "19F", "20Ne"};
}  // namespace

double isotope_mass(int Z) {
    if (Z < 1 || Z > 10) throw std::invalid_argument("isotope_mass: Z out of range");
    return kMass[Z - 1];
}

const char* isotope_label(int Z) {
    if (Z < 1 || Z > 10) throw std::invalid_argument("isotope_label: Z out of range");
    return kLabel[Z - 1];
}

}  // namespace macqc
