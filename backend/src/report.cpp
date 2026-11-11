#include "macqc/report.h"

#include <cmath>
#include <cstdio>

namespace macqc {

namespace {
std::string pad_left(const std::string& s, int width) {
    if (static_cast<int>(s.size()) >= width) return s;
    return std::string(static_cast<std::size_t>(width) - s.size(), ' ') + s;
}
}  // namespace

std::string fixed(double value, int width, int decimals) {
    char buf[512];
    std::snprintf(buf, sizeof buf, "%*.*f", width, decimals, value);
    return buf;
}

std::string integer(long value, int width) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%*ld", width, value);
    return buf;
}

std::string fortran_d(double value, int width, int ndigits) {
    double mantissa;
    int exponent;
    if (value == 0.0) {
        mantissa = 0.0;
        exponent = 0;
    } else {
        exponent = static_cast<int>(std::floor(std::log10(std::fabs(value)))) + 1;
        mantissa = value / std::pow(10.0, exponent);
    }

    // Rounding may push the mantissa to 1.000...; renormalise (as the reference does).
    char probe[512];
    std::snprintf(probe, sizeof probe, "%.*f", ndigits, std::fabs(mantissa));
    if (std::fabs(std::strtod(probe, nullptr)) >= 1.0) {
        mantissa /= 10.0;
        exponent += 1;
    }

    char buf[600];
    std::snprintf(buf, sizeof buf, "%.*fD%+03d", ndigits, mantissa, exponent);
    return pad_left(buf, width);
}

void matout(std::string& out, const Matrix& a, const std::string& label_in) {
    std::string label = label_in.substr(0, 4);
    label.resize(4, ' ');
    out += "\n\n\n    THE " + label + " ARRAY\n";

    std::string heading(15, ' ');
    for (std::size_t j = 0; j < a.cols(); ++j)
        heading += std::string(10, ' ') + integer(static_cast<long>(j) + 1, 3) + std::string(6, ' ');
    heading.erase(heading.find_last_not_of(' ') + 1);
    out += heading + "\n";

    for (std::size_t i = 0; i < a.rows(); ++i) {
        out += integer(static_cast<long>(i) + 1, 10) + std::string(5, ' ');
        for (std::size_t j = 0; j < a.cols(); ++j) out += " " + fortran_d(a(i, j), 18, 10);
        out += "\n";
    }
}

}  // namespace macqc
