// Fortran-style number formatting and text-report helpers.
#pragma once

#include <string>

#include "macqc/linalg.h"

namespace macqc {

// printf-style fixed format, e.g. fixed(x, 11, 6) == Fortran F11.6.
std::string fixed(double value, int width, int decimals);
// Right-aligned integer, Fortran Iw.
std::string integer(long value, int width);
// Fortran D format (Dw.d), e.g. fortran_d(x, 18, 10).
std::string fortran_d(double value, int width, int ndigits);

// Append a matrix in the original MATOUT layout (label is padded/cut to 4 chars).
void matout(std::string& out, const Matrix& a, const std::string& label);

}  // namespace macqc
