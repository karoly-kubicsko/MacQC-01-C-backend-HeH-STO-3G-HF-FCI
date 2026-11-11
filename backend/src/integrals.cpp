#include "macqc/integrals.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace macqc {

const double kPi = 3.1415926535898;

double derfother(double arg) {
    const double p = 0.3275911;
    const double a[5] = {0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429};
    const double t = 1.0 / (1.0 + p * arg);
    double tn = t;
    double poly = a[0] * tn;
    for (int i = 1; i < 5; ++i) {
        tn *= t;
        poly += a[i] * tn;
    }
    return 1.0 - poly * std::exp(-arg * arg);
}

double boys_f0(double arg, bool exact_erf) {
    if (arg < 1.0e-6) return 1.0 - arg / 3.0;
    const double e = exact_erf ? std::erf(std::sqrt(arg)) : derfother(std::sqrt(arg));
    return std::sqrt(kPi / arg) * e / 2.0;
}

double overlap_primitive(double a, double b, double rab2) {
    return std::pow(kPi / (a + b), 1.5) * std::exp(-a * b * rab2 / (a + b));
}

double kinetic_primitive(double a, double b, double rab2) {
    return a * b / (a + b) * (3.0 - 2.0 * a * b * rab2 / (a + b)) * std::pow(kPi / (a + b), 1.5) *
           std::exp(-a * b * rab2 / (a + b));
}

double nuclear_attraction_primitive(double a, double b, double rab2, double rcp2, double zc,
                                    bool exact_erf) {
    double value = 2.0 * kPi / (a + b) * boys_f0((a + b) * rcp2, exact_erf);
    value *= std::exp(-a * b * rab2 / (a + b));
    return -value * zc;
}

double two_electron_primitive(double a, double b, double c, double d, double rab2, double rcd2,
                              double rpq2, bool exact_erf) {
    return 2.0 * std::pow(kPi, 2.5) / ((a + b) * (c + d) * std::sqrt(a + b + c + d)) *
           boys_f0((a + b) * (c + d) * rpq2 / (a + b + c + d), exact_erf) *
           std::exp(-a * b * rab2 / (a + b) - c * d * rcd2 / (c + d));
}

Integrals calculate_integrals(const DiatomicBasis& basis, bool exact_erf) {
    const int n = basis.n_gaussians;
    if (n < 1 || n > 3) throw std::invalid_argument("n_gaussians must be 1, 2 or 3");

    // STO-nG expansion coefficients / exponents for zeta = 1 (columns: n = 1, 2, 3).
    static const double coef[3][3] = {
        {1.000000, 0.678914, 0.444635},
        {0.000000, 0.430129, 0.535328},
        {0.000000, 0.000000, 0.154329}};
    static const double expon[3][3] = {
        {0.270950, 0.151623, 0.109818},
        {0.000000, 0.851819, 0.405771},
        {0.000000, 0.000000, 2.227660}};

    const double r = basis.r;
    const double r2 = r * r;
    std::vector<double> a1(n), d1(n), a2(n), d2(n);
    const int col = n - 1;
    for (int i = 0; i < n; ++i) {
        a1[i] = expon[i][col] * (basis.zeta1 * basis.zeta1);
        d1[i] = coef[i][col] * std::pow(2.0 * a1[i] / kPi, 0.75);
        a2[i] = expon[i][col] * (basis.zeta2 * basis.zeta2);
        d2[i] = coef[i][col] * std::pow(2.0 * a2[i] / kPi, 0.75);
    }

    Integrals v;
    const double za = basis.za, zb = basis.zb;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double rap = a2[j] * r / (a1[i] + a2[j]);
            const double rap2 = rap * rap;
            const double rbp2 = (r - rap) * (r - rap);

            v.S12 += overlap_primitive(a1[i], a2[j], r2) * d1[i] * d2[j];
            v.T11 += kinetic_primitive(a1[i], a1[j], 0.0) * d1[i] * d1[j];
            v.T12 += kinetic_primitive(a1[i], a2[j], r2) * d1[i] * d2[j];
            v.T22 += kinetic_primitive(a2[i], a2[j], 0.0) * d2[i] * d2[j];

            v.V11A += nuclear_attraction_primitive(a1[i], a1[j], 0.0, 0.0, za, exact_erf) * d1[i] * d1[j];
            v.V12A += nuclear_attraction_primitive(a1[i], a2[j], r2, rap2, za, exact_erf) * d1[i] * d2[j];
            v.V22A += nuclear_attraction_primitive(a2[i], a2[j], 0.0, r2, za, exact_erf) * d2[i] * d2[j];
            v.V11B += nuclear_attraction_primitive(a1[i], a1[j], 0.0, r2, zb, exact_erf) * d1[i] * d1[j];
            v.V12B += nuclear_attraction_primitive(a1[i], a2[j], r2, rbp2, zb, exact_erf) * d1[i] * d2[j];
            v.V22B += nuclear_attraction_primitive(a2[i], a2[j], 0.0, 0.0, zb, exact_erf) * d2[i] * d2[j];
        }
    }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l) {
                    const double rap = a2[i] * r / (a2[i] + a1[j]);
                    const double raq = a2[k] * r / (a2[k] + a1[l]);
                    const double rbq = r - raq;
                    const double rpq = rap - raq;
                    const double rap2 = rap * rap;
                    const double rbq2 = rbq * rbq;
                    const double rpq2 = rpq * rpq;

                    v.V1111 += two_electron_primitive(a1[i], a1[j], a1[k], a1[l], 0.0, 0.0, 0.0, exact_erf) *
                               d1[i] * d1[j] * d1[k] * d1[l];
                    v.V2111 += two_electron_primitive(a2[i], a1[j], a1[k], a1[l], r2, 0.0, rap2, exact_erf) *
                               d2[i] * d1[j] * d1[k] * d1[l];
                    v.V2121 += two_electron_primitive(a2[i], a1[j], a2[k], a1[l], r2, r2, rpq2, exact_erf) *
                               d2[i] * d1[j] * d2[k] * d1[l];
                    v.V2211 += two_electron_primitive(a2[i], a2[j], a1[k], a1[l], 0.0, 0.0, r2, exact_erf) *
                               d2[i] * d2[j] * d1[k] * d1[l];
                    v.V2221 += two_electron_primitive(a2[i], a2[j], a2[k], a1[l], 0.0, r2, rbq2, exact_erf) *
                               d2[i] * d2[j] * d2[k] * d1[l];
                    v.V2222 += two_electron_primitive(a2[i], a2[j], a2[k], a2[l], 0.0, 0.0, 0.0, exact_erf) *
                               d2[i] * d2[j] * d2[k] * d2[l];
                }
    return v;
}

}  // namespace macqc
