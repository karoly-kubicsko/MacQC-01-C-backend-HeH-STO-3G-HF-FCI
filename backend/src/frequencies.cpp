#include "macqc/frequencies.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "macqc/constants.h"

namespace macqc {

namespace {

using namespace codata2018;

// sqrt(E_h / (a0^2 u)) / (2 pi c), in cm^-1: wavenumber of an eigenvalue of 1 E_h/(bohr^2 u).
double wavenumber_factor() { return std::sqrt(kHartree / (kBohr * kBohr * kAtomicMassUnit)) / (2.0 * kPi * kSpeedOfLight * 100.0); }

// IR intensity (km/mol) of |dmu/dQ|^2 = 1 e^2/u: N_A e^2 / (12 eps0 c^2 u).
double intensity_factor() {
    return kAvogadro * kElementaryCharge * kElementaryCharge /
           (12.0 * kVacuumPermittivity * kSpeedOfLight * kSpeedOfLight * kAtomicMassUnit) / 1000.0;
}

double wavenumber(double lambda) {
    const double w = wavenumber_factor() * std::sqrt(std::fabs(lambda));
    return lambda < 0.0 ? -w : w;
}

// Gram-Schmidt: appends v (normalised) to Q if it is independent; returns true if appended.
bool orthonormalize_into(std::vector<std::vector<double>>& Q, std::vector<double> v, double tol = 1.0e-6) {
    for (int pass = 0; pass < 2; ++pass)
        for (const auto& q : Q) {
            double d = 0.0;
            for (std::size_t i = 0; i < v.size(); ++i) d += q[i] * v[i];
            for (std::size_t i = 0; i < v.size(); ++i) v[i] -= d * q[i];
        }
    double n = 0.0;
    for (double x : v) n += x * x;
    n = std::sqrt(n);
    if (n < tol) return false;
    for (double& x : v) x /= n;
    Q.push_back(v);
    return true;
}

struct Inertia {
    std::array<double, 3> com{};
    std::array<double, 3> moments{};  // ascending, u bohr^2
    Matrix axes;                      // columns
};

Inertia inertia(const Molecule& mol, const std::vector<double>& m) {
    Inertia I;
    double M = 0.0;
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) {
        M += m[a];
        for (std::size_t k = 0; k < 3; ++k) I.com[k] += m[a] * mol.atoms[a].position[k];
    }
    for (double& c : I.com) c /= M;
    Matrix T(3, 3);
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) {
        double r[3], r2 = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            r[k] = mol.atoms[a].position[k] - I.com[k];
            r2 += r[k] * r[k];
        }
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) T(i, j) += m[a] * ((i == j ? r2 : 0.0) - r[i] * r[j]);
    }
    const EigenResult e = symmetric_eigen(T);
    for (std::size_t k = 0; k < 3; ++k) I.moments[k] = std::max(0.0, e.values[k]);
    I.axes = e.vectors;
    return I;
}

}  // namespace

int symmetry_number(const Molecule& mol) {
    if (mol.atoms.size() == 2) return mol.atoms[0].Z == mol.atoms[1].Z ? 2 : 1;
    if (mol.atoms.size() == 3) {
        auto d = [&](int a, int b) {
            double s = 0.0;
            for (std::size_t k = 0; k < 3; ++k) {
                const double x = mol.atoms[static_cast<std::size_t>(a)].position[k] - mol.atoms[static_cast<std::size_t>(b)].position[k];
                s += x * x;
            }
            return std::sqrt(s);
        };
        if (mol.atoms[1].Z == mol.atoms[2].Z && std::fabs(d(0, 1) - d(0, 2)) < 1.0e-4) return 2;
    }
    return 1;
}

double joule_per_mol_to_hartree(double x) { return x / (kHartree * kAvogadro); }

HarmonicAnalysis harmonic_analysis(const Molecule& mol, const Matrix& H,
                                   const std::vector<std::array<double, 3>>* apt) {
    const std::size_t nat = mol.atoms.size(), n3 = 3 * nat;
    if (H.rows() != n3 || H.cols() != n3) throw std::invalid_argument("harmonic_analysis: Hessian has the wrong size");
    HarmonicAnalysis res;
    for (const Atom& a : mol.atoms) {
        res.masses.push_back(isotope_mass(a.Z));
        res.isotopes.push_back(isotope_label(a.Z));
    }
    const Inertia I = inertia(mol, res.masses);
    res.center_of_mass = I.com;
    res.linear = I.moments[0] < 1.0e-6 * std::max(1.0, I.moments[2]);
    // translations and rotations (mass-weighted)
    std::vector<std::vector<double>> TR;
    for (std::size_t k = 0; k < 3; ++k) {
        std::vector<double> t(n3, 0.0);
        for (std::size_t a = 0; a < nat; ++a) t[3 * a + k] = std::sqrt(res.masses[a]);
        orthonormalize_into(TR, t);
    }
    for (std::size_t ax = 0; ax < 3; ++ax) {
        double e[3] = {I.axes(0, ax), I.axes(1, ax), I.axes(2, ax)};
        std::vector<double> r(n3, 0.0);
        for (std::size_t a = 0; a < nat; ++a) {
            double p[3];
            for (std::size_t k = 0; k < 3; ++k) p[k] = mol.atoms[a].position[k] - I.com[k];
            const double c[3] = {e[1] * p[2] - e[2] * p[1], e[2] * p[0] - e[0] * p[2], e[0] * p[1] - e[1] * p[0]};
            for (std::size_t k = 0; k < 3; ++k) r[3 * a + k] = std::sqrt(res.masses[a]) * c[k];
        }
        orthonormalize_into(TR, r);
    }
    const std::size_t ntr = TR.size();
    std::vector<std::vector<double>> Q = TR;  // TR first, then the vibrational complement
    for (std::size_t x = 0; x < n3 && Q.size() < n3; ++x) {
        std::vector<double> u(n3, 0.0);
        u[x] = 1.0;
        orthonormalize_into(Q, u);
    }
    const std::size_t nvib = n3 - ntr;
    res.n_vib = static_cast<int>(nvib);
    // mass-weighted Hessian
    Matrix Hm(n3, n3);
    for (std::size_t i = 0; i < n3; ++i)
        for (std::size_t j = 0; j < n3; ++j)
            Hm(i, j) = 0.5 * (H(i, j) + H(j, i)) / std::sqrt(res.masses[i / 3] * res.masses[j / 3]);
    auto project = [&](std::size_t first, std::size_t count) {
        Matrix P(count, count);
        for (std::size_t a = 0; a < count; ++a)
            for (std::size_t b = 0; b < count; ++b) {
                double s = 0.0;
                const auto& qa = Q[first + a];
                const auto& qb = Q[first + b];
                for (std::size_t i = 0; i < n3; ++i) {
                    if (qa[i] == 0.0) continue;
                    double t = 0.0;
                    for (std::size_t j = 0; j < n3; ++j) t += Hm(i, j) * qb[j];
                    s += qa[i] * t;
                }
                P(a, b) = s;
            }
        return P;
    };
    // residual curvature of the translations and rotations (a quality indicator)
    {
        const EigenResult et = symmetric_eigen(project(0, ntr));
        for (double l : et.values) res.low_frequencies.push_back(wavenumber(l));
    }
    const EigenResult ev = symmetric_eigen(project(ntr, nvib));
    // C2 operation of a C2v triatomic (H1 <-> H2) for the symmetry species
    const bool c2v = nat == 3 && symmetry_number(mol) == 2;
    std::array<double, 3> c2axis{};
    if (c2v) {
        for (std::size_t k = 0; k < 3; ++k) {
            double u1 = 0.0, u2 = 0.0;
            for (std::size_t j = 0; j < 3; ++j) {
                u1 += std::pow(mol.atoms[1].position[j] - mol.atoms[0].position[j], 2);
                u2 += std::pow(mol.atoms[2].position[j] - mol.atoms[0].position[j], 2);
            }
            c2axis[k] = (mol.atoms[1].position[k] - mol.atoms[0].position[k]) / std::sqrt(u1) +
                        (mol.atoms[2].position[k] - mol.atoms[0].position[k]) / std::sqrt(u2);
        }
        const double n = std::sqrt(c2axis[0] * c2axis[0] + c2axis[1] * c2axis[1] + c2axis[2] * c2axis[2]);
        for (double& c : c2axis) c /= n;
    }
    const double ifac = intensity_factor();
    for (std::size_t k = 0; k < nvib; ++k) {
        NormalMode m;
        std::vector<double> l(n3, 0.0);
        for (std::size_t a = 0; a < nvib; ++a)
            for (std::size_t i = 0; i < n3; ++i) l[i] += Q[ntr + a][i] * ev.vectors(a, k);
        m.wavenumber = wavenumber(ev.values[k]);
        m.imaginary = ev.values[k] < 0.0;
        double norm2 = 0.0;
        std::vector<double> x(n3);
        for (std::size_t i = 0; i < n3; ++i) {
            x[i] = l[i] / std::sqrt(res.masses[i / 3]);
            norm2 += x[i] * x[i];
        }
        m.reduced_mass = 1.0 / norm2;
        const double nu = std::fabs(m.wavenumber) * 100.0;  // m^-1
        m.force_constant = std::pow(2.0 * kPi * kSpeedOfLight * nu, 2) * m.reduced_mass * kAtomicMassUnit / 100.0;
        if (m.imaginary) m.force_constant = -m.force_constant;
        m.displacement.assign(nat, {0.0, 0.0, 0.0});
        const double xn = std::sqrt(norm2);
        for (std::size_t a = 0; a < nat; ++a)
            for (std::size_t c = 0; c < 3; ++c) m.displacement[a][c] = x[3 * a + c] / xn;
        // a fixed sign: the largest component positive
        double big = 0.0;
        for (std::size_t i = 0; i < n3; ++i)
            if (std::fabs(x[i]) > std::fabs(big) + 1e-12) big = x[i];
        if (big < 0.0) {
            for (auto& d : m.displacement)
                for (double& c : d) c = -c;
            for (double& v : l) v = -v;
        }
        if (apt) {
            for (std::size_t c = 0; c < 3; ++c) {
                double s = 0.0;
                for (std::size_t i = 0; i < n3; ++i) s += (*apt)[i][c] * l[i] / std::sqrt(res.masses[i / 3]);
                m.dmu_dq[c] = s;
            }
            m.ir_intensity = ifac * (m.dmu_dq[0] * m.dmu_dq[0] + m.dmu_dq[1] * m.dmu_dq[1] + m.dmu_dq[2] * m.dmu_dq[2]);
        }
        if (nat == 2) {
            m.symmetry = mol.atoms[0].Z == mol.atoms[1].Z ? "Sigma_g+" : "Sigma+";
        } else if (c2v) {
            // character of C2: sum_i (R d_i) . d_{pi(i)}, R v = 2 (a.v) a - v, pi swaps atoms 1 and 2
            const std::size_t perm[3] = {0, 2, 1};
            double ch = 0.0;
            for (std::size_t a = 0; a < 3; ++a) {
                const auto& d = m.displacement[a];
                const double ad = c2axis[0] * d[0] + c2axis[1] * d[1] + c2axis[2] * d[2];
                const auto& t = m.displacement[perm[a]];
                for (std::size_t c = 0; c < 3; ++c) ch += (2.0 * ad * c2axis[c] - d[c]) * t[c];
            }
            m.symmetry = ch > 0.5 ? "A1" : (ch < -0.5 ? "B2" : "");
        }
        res.modes.push_back(m);
    }
    res.has_intensities = apt != nullptr;
    return res;
}

ThermoResult thermochemistry(const Molecule& mol, const HarmonicAnalysis& h, double electronic_energy,
                             double T, double P, int sigma, double cutoff) {
    if (!(T > 0.0) || !(P > 0.0)) throw std::invalid_argument("thermochemistry: temperature and pressure must be positive");
    ThermoResult r;
    r.temperature = T;
    r.pressure = P;
    r.symmetry_number = sigma;
    r.linear = h.linear;
    r.electronic_energy = electronic_energy;
    const double R = kGasConstant, kT = kBoltzmann * T;
    for (double m : h.masses) r.mass += m;
    // translation
    {
        const double M = r.mass * kAtomicMassUnit;
        const double qt = std::pow(2.0 * kPi * M * kT / (kPlanck * kPlanck), 1.5) * kT / P;
        r.translation.q = qt;
        r.translation.entropy = R * (std::log(qt) + 2.5);
        r.translation.energy = 1.5 * R * T;
        r.translation.cv = 1.5 * R;
    }
    // rotation
    {
        const Inertia I = inertia(mol, h.masses);
        r.moments = I.moments;
        for (std::size_t k = 0; k < 3; ++k) {
            const double Ikg = I.moments[k] * kAtomicMassUnit * kBohr * kBohr;
            if (Ikg <= 0.0 || (h.linear && k == 0)) continue;
            const double B = kPlanck / (8.0 * kPi * kPi * Ikg);  // Hz
            r.rotational_constants_ghz[k] = B / 1.0e9;
            r.rotational_constants_cm[k] = B / (kSpeedOfLight * 100.0);
            r.rotational_temperatures[k] = kPlanck * B / kBoltzmann;
        }
        if (h.linear) {
            const double th = r.rotational_temperatures[2];
            const double q = T / (sigma * th);
            r.rotation.q = q;
            r.rotation.entropy = R * (std::log(q) + 1.0);
            r.rotation.energy = R * T;
            r.rotation.cv = R;
        } else {
            const double prod = r.rotational_temperatures[0] * r.rotational_temperatures[1] * r.rotational_temperatures[2];
            const double q = std::sqrt(kPi) / sigma * std::pow(T, 1.5) / std::sqrt(prod);
            r.rotation.q = q;
            r.rotation.entropy = R * (std::log(q) + 1.5);
            r.rotation.energy = 1.5 * R * T;
            r.rotation.cv = 1.5 * R;
        }
    }
    // vibration
    r.vibration.q = 1.0;
    for (const NormalMode& m : h.modes) {
        if (m.imaginary) {
            ++r.excluded_modes;
            char buf[160];
            std::snprintf(buf, sizeof buf, "The imaginary frequency %.2fi cm-1 is left out of the vibrational sums.",
                          std::fabs(m.wavenumber));
            r.warnings.push_back(buf);
            continue;
        }
        if (m.wavenumber < cutoff) {
            ++r.excluded_modes;
            char buf[160];
            std::snprintf(buf, sizeof buf,
                          "The frequency %.2f cm-1 (below %.0f cm-1) is left out of the vibrational sums.",
                          m.wavenumber, cutoff);
            r.warnings.push_back(buf);
            continue;
        }
        const double th = kPlanck * kSpeedOfLight * m.wavenumber * 100.0 / kBoltzmann;
        r.vibrational_temperatures.push_back(th);
        const double x = th / T;
        r.zero_point_energy += 0.5 * R * th;
        r.vibration.energy += R * th * (0.5 + 1.0 / std::expm1(x));
        r.vibration.entropy += R * (x / std::expm1(x) - std::log(-std::expm1(-x)));
        r.vibration.cv += R * x * x * std::exp(x) / (std::expm1(x) * std::expm1(x));
        r.vibration.q *= std::exp(-0.5 * x) / (-std::expm1(-x));
    }
    // electronic: closed-shell singlet
    r.electronic.q = 1.0;
    r.total.energy = r.translation.energy + r.rotation.energy + r.vibration.energy + r.electronic.energy;
    r.total.entropy = r.translation.entropy + r.rotation.entropy + r.vibration.entropy + r.electronic.entropy;
    r.total.cv = r.translation.cv + r.rotation.cv + r.vibration.cv + r.electronic.cv;
    r.total.q = r.translation.q * r.rotation.q * r.vibration.q * r.electronic.q;
    r.thermal_energy = r.total.energy;
    r.enthalpy = r.thermal_energy + R * T;
    r.gibbs = r.enthalpy - T * r.total.entropy;
    r.helmholtz = r.thermal_energy - T * r.total.entropy;
    r.total_zpe = electronic_energy + joule_per_mol_to_hartree(r.zero_point_energy);
    r.total_u = electronic_energy + joule_per_mol_to_hartree(r.thermal_energy);
    r.total_h = electronic_energy + joule_per_mol_to_hartree(r.enthalpy);
    r.total_g = electronic_energy + joule_per_mol_to_hartree(r.gibbs);
    r.total_a = electronic_energy + joule_per_mol_to_hartree(r.helmholtz);
    return r;
}

}  // namespace macqc
