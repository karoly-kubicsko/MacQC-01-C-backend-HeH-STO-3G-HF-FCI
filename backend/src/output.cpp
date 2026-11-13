#include "macqc/output.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <ctime>

#include "macqc/constants.h"
#include "macqc/gto_integrals.h"
#include "macqc/rhf_derivatives.h"

namespace macqc {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kWidth = 88;

std::string atom_label(const Molecule& m, std::size_t a) { return m.atoms[a].symbol + std::to_string(a + 1); }

}  // namespace

std::string strprintf(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    va_list copy;
    va_copy(copy, args);
    const int n = std::vsnprintf(nullptr, 0, fmt, copy);
    va_end(copy);
    std::string out(static_cast<std::size_t>(std::max(0, n)), '\0');
    std::vsnprintf(out.data(), out.size() + 1, fmt, args);
    va_end(args);
    return out;
}

void OutputWriter::write(const std::string& text) {
    text_ += text;
    if (sink_) sink_(text);
}

void OutputWriter::section(const std::string& title) {
    std::string s = "\n---- " + title + " ";
    // count UTF-8 code points for the visual width
    std::size_t width = 0;
    for (unsigned char c : s)
        if ((c & 0xC0) != 0x80) ++width;
    if (width < static_cast<std::size_t>(kWidth)) s += std::string(static_cast<std::size_t>(kWidth) - width, '-');
    line(s);
}

void OutputWriter::subsection(const std::string& title) { line("\n  " + title); }

std::string local_time_string() {
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[64];
    std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S", &tm);
    return buf;
}

std::string duration_string(double s) {
    if (s < 60.0) return strprintf("%.1f s", s);
    if (s < 3600.0) return strprintf("%d min %02d s", static_cast<int>(s / 60.0), static_cast<int>(std::fmod(s, 60.0)));
    return strprintf("%d h %02d min", static_cast<int>(s / 3600.0), static_cast<int>(std::fmod(s, 3600.0) / 60.0));
}

std::string geometry_table(const Molecule& mol) {
    std::string s = "    Atom   Z          x (Å)          y (Å)          z (Å)      x (bohr)      y (bohr)      z (bohr)\n";
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) {
        const auto& p = mol.atoms[a].position;
        s += strprintf("    %-5s %2d  %13.8f  %13.8f  %13.8f  %12.8f  %12.8f  %12.8f\n", atom_label(mol, a).c_str(),
                       mol.atoms[a].Z, p[0] * kBohrInAngstrom, p[1] * kBohrInAngstrom, p[2] * kBohrInAngstrom, p[0],
                       p[1], p[2]);
    }
    return s;
}

std::string parameter_lines(const Molecule& mol, const std::string& indent) {
    std::string s;
    for (const GeometricParameter& p : geometric_parameters(mol)) {
        if (p.angle)
            s += indent + strprintf("%-14s = %12.6f deg\n", p.name.c_str(), p.value);
        else
            s += indent + strprintf("%-14s = %12.8f Å   (%.8f bohr)\n", p.name.c_str(), p.value * kBohrInAngstrom, p.value);
    }
    return s;
}

std::string basis_block(const BasisSet& b, const std::vector<std::string>& notes) {
    std::string s = "    " + b.name + ": spherical harmonics, Gaussian94 files from basissetexchange.org\n";
    std::vector<std::string> seen;
    for (const auto& a : b.atoms) {
        if (std::find(seen.begin(), seen.end(), a.symbol) != seen.end()) continue;
        seen.push_back(a.symbol);
        s += strprintf("    %-2s  %-13s %-16s -> %-16s %4d functions\n", a.symbol.c_str(), a.basis.c_str(),
                       a.primitives.c_str(), a.contraction.c_str(), a.functions);
    }
    for (const auto& n : notes) s += "    Note: " + n + "\n";
    s += strprintf("    %d basis functions, %zu shells, highest angular momentum %c (l = %d)\n", b.n_functions,
                   b.shells.size(), angular_letter(b.max_l), b.max_l);
    return s;
}

std::string scf_header() {
    return "    iter        total energy (E_h)        dE (E_h)   orbital gradient\n";
}

std::string scf_iteration_line(const RHFIteration& it) {
    return strprintf("    %4d  %24.12f  %14.3e  %12.3e%s\n", it.iteration, it.energy, it.delta_e, it.gradient,
                     it.diis ? "  DIIS" : "");
}

std::vector<double> mulliken_charges(const Molecule& mol, const BasisSet& b, const Matrix& P, const Matrix& S) {
    std::vector<double> q(mol.atoms.size(), 0.0);
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) q[a] = mol.atoms[a].Z;
    const std::size_t n = P.rows();
    for (const Shell& sh : b.shells)
        for (int f = 0; f < sh.size(); ++f) {
            const std::size_t mu = static_cast<std::size_t>(sh.first_function + f);
            double ps = 0.0;
            for (std::size_t nu = 0; nu < n; ++nu) ps += P(mu, nu) * S(nu, mu);
            q[static_cast<std::size_t>(sh.atom)] -= ps;
        }
    return q;
}

// Löwdin charges: q_A = Z_A - sum_{mu on A} (S^1/2 P S^1/2)_{mu mu}.
std::vector<double> lowdin_charges(const Molecule& mol, const BasisSet& b, const Matrix& P, const Matrix& S) {
    const std::size_t n = S.rows();
    const EigenResult e = symmetric_eigen(S);
    Matrix half(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            double t = 0.0;
            for (std::size_t k = 0; k < n; ++k)
                t += e.vectors(i, k) * std::sqrt(std::max(e.values[k], 0.0)) * e.vectors(j, k);
            half(i, j) = t;
        }
    const Matrix sps = matmul(matmul(half, P), half);
    std::vector<double> q(mol.atoms.size(), 0.0);
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) q[a] = mol.atoms[a].Z;
    for (const Shell& sh : b.shells)
        for (int f = 0; f < sh.size(); ++f) {
            const std::size_t mu = static_cast<std::size_t>(sh.first_function + f);
            q[static_cast<std::size_t>(sh.atom)] -= sps(mu, mu);
        }
    return q;
}

std::string scf_result_block(const RHFResult& h, const Molecule& mol, const BasisSet& b, bool full_orbitals) {
    std::string s;
    s += h.converged ? strprintf("    SCF converged in %d iterations (integrals %.2f s, SCF %.2f s, %d threads)\n",
                                 h.iterations, h.integral_seconds, h.scf_seconds, h.threads)
                     : strprintf("    SCF NOT converged after %d iterations\n", h.iterations);
    s += strprintf("    overlap matrix: smallest eigenvalue %.3e; %s\n", h.smallest_overlap_eigenvalue,
                   h.removed_functions ? strprintf("%d near-linear dependencies removed", h.removed_functions).c_str()
                                       : "no functions removed");
    s += "\n    Energy components (E_h)\n";
    s += strprintf("      total energy (RHF)        %20.10f\n", h.total_energy);
    s += strprintf("      electronic energy         %20.10f\n", h.electronic_energy);
    s += strprintf("      nuclear repulsion         %20.10f\n", h.nuclear_repulsion);
    s += strprintf("      kinetic energy            %20.10f\n", h.kinetic_energy);
    s += strprintf("      nuclear attraction        %20.10f\n", h.nuclear_attraction_energy);
    s += strprintf("      electron repulsion        %20.10f\n", h.two_electron_energy);
    s += strprintf("      virial ratio -V/T         %20.10f\n", h.virial_ratio);
    const int nocc = h.n_occupied;
    const int nall = static_cast<int>(h.orbital_energies.size());
    const int nshow = full_orbitals ? nall : std::min(nall, nocc + 5);
    s += strprintf("\n    Orbital energies (E_h), %d doubly occupied (*)%s\n", nocc,
                   nshow < nall ? strprintf(", the lowest %d virtual orbitals of %d", nshow - nocc, nall - nocc).c_str() : "");
    for (int k = 0; k < nshow; ++k) {
        s += strprintf("  %5d%s %14.8f", k + 1, k < nocc ? "*" : " ", h.orbital_energies[static_cast<std::size_t>(k)]);
        if (k % 4 == 3 || k + 1 == nshow) s += "\n";
    }
    if (nocc >= 1 && nocc < nall)
        s += strprintf("    HOMO %.8f E_h, LUMO %.8f E_h, gap %.8f E_h (%.4f eV)\n", h.orbital_energies[static_cast<std::size_t>(nocc - 1)],
                       h.orbital_energies[static_cast<std::size_t>(nocc)],
                       h.orbital_energies[static_cast<std::size_t>(nocc)] - h.orbital_energies[static_cast<std::size_t>(nocc - 1)],
                       (h.orbital_energies[static_cast<std::size_t>(nocc)] - h.orbital_energies[static_cast<std::size_t>(nocc - 1)]) *
                           codata2018::kHartree / codata2018::kElementaryCharge);
    // Mulliken charges and the SCF dipole moment from the density of the final orbitals
    const std::size_t N = h.C.rows();
    Matrix P(N, N);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) {
            double t = 0.0;
            for (int k = 0; k < nocc; ++k) t += h.C(i, static_cast<std::size_t>(k)) * h.C(j, static_cast<std::size_t>(k));
            P(i, j) = 2.0 * t;
        }
    const std::vector<double> q = mulliken_charges(mol, b, P, h.S);
    const std::vector<double> ql = lowdin_charges(mol, b, P, h.S);
    s += "\n    Atomic charges (SCF density)      Mulliken        Löwdin\n";
    double sum = 0.0, suml = 0.0;
    for (std::size_t a = 0; a < q.size(); ++a) {
        s += strprintf("      %-28s %12.6f  %12.6f\n", atom_label(mol, a).c_str(), q[a], ql[a]);
        sum += q[a];
        suml += ql[a];
    }
    s += strprintf("      %-28s %12.6f  %12.6f\n", "sum", sum + 0.0, suml + 0.0);
    const std::array<double, 3> mu = dipole_moment(mol, b, P);
    const double nm = std::sqrt(mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]);
    s += "\n    Dipole moment (SCF density, origin at 0)\n";
    s += strprintf("      x %12.8f   y %12.8f   z %12.8f   |mu| %12.8f e bohr   (%.6f D)\n", mu[0], mu[1], mu[2], nm,
                   nm / codata2018::kDebyeInAtomicUnits);
    return s;
}

std::string cc_header(const std::string& method) {
    return "    " + method + " iterations\n    iter        total energy (E_h)     correlation (E_h)        dE (E_h)   max |dt|\n";
}

std::string cc_iteration_line(const CCIteration& it) {
    return strprintf("    %4d  %24.12f  %20.12f  %14.3e  %9.2e\n", it.iteration, it.energy, it.correlation, it.delta_e,
                     it.residual);
}

std::string fci_header() {
    return "    Davidson iterations\n    iter        total energy (E_h)        residual   subspace\n";
}

std::string fci_iteration_line(const FCIIteration& it) {
    return strprintf("    %4d  %24.12f  %12.3e  %6d\n", it.iteration, it.energy, it.residual, it.subspace);
}

std::string gradient_block(const Molecule& mol, const std::vector<double>& g, const std::string& how) {
    std::string s = "    Cartesian gradient (E_h/bohr), " + how + "\n";
    s += "    Atom            dE/dx            dE/dy            dE/dz\n";
    double mx = 0.0, rms = 0.0;
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) {
        s += strprintf("    %-5s %16.10f %16.10f %16.10f\n", atom_label(mol, a).c_str(), g[3 * a], g[3 * a + 1],
                       g[3 * a + 2]);
        for (int k = 0; k < 3; ++k) {
            mx = std::max(mx, std::fabs(g[3 * a + static_cast<std::size_t>(k)]));
            rms += g[3 * a + static_cast<std::size_t>(k)] * g[3 * a + static_cast<std::size_t>(k)];
        }
    }
    s += strprintf("    largest component %.3e, RMS %.3e E_h/bohr\n", mx, std::sqrt(rms / static_cast<double>(g.size())));
    return s;
}

std::string optimization_step_block(const OptimizationStep& st, const ConvergenceCriteria& c) {
    std::string s;
    s += "\n    Internal coordinates (bohr, radians; force = -dE/dq)\n";
    s += "    Name  Definition               Value (bohr/rad)    Value (Å/deg)       Force     Step\n";
    for (std::size_t i = 0; i < st.q.size(); ++i) {
        const bool angle = st.names[i] != "R";
        const double shown = angle ? (st.names[i] == "A" ? st.q[i] * 180.0 / kPi : st.q[i]) : st.q[i] * kBohrInAngstrom;
        s += strprintf("    %-4s  %-22s %16.8f %16.8f %11.6f %8.5f\n", (st.names[i] + std::to_string(i + 1)).c_str(),
                       st.definitions[i].c_str(), st.q[i], shown, -st.gq[i], i < st.dq.size() ? st.dq[i] : 0.0);
    }
    s += "    " + st.hessian_update + "\n";
    if (st.has_ratio)
        s += strprintf("    previous step: energy change / predicted = %.3f\n", st.ratio);
    if (!st.converged)
        s += strprintf("    RFO step: lambda = %.6e, |step| = %.6f%s, trust radius %.4f, predicted energy change %.3e "
                       "E_h\n",
                       st.rfo_lambda, st.step_norm, st.step_scaled ? " (scaled to the trust radius)" : "", st.trust_radius,
                       st.predicted_change);
    s += "\n    Convergence             Value      Threshold   Converged?\n";
    for (const auto& item : st.convergence) {
        if (!item.applicable) {
            s += strprintf("    %-22s %10s %12.2e   --\n", item.name.c_str(), "--", item.threshold);
            continue;
        }
        s += strprintf("    %-22s %10.2e %12.2e   %s\n", item.name.c_str(), item.value, item.threshold, item.met ? "YES" : "NO");
    }
    (void)c;
    for (const auto& w : st.warnings) s += "    Warning: " + w + "\n";
    return s;
}

std::string optimization_result_block(const OptimizationResult& r) {
    std::string s;
    if (r.converged) {
        s += "\n    *** Stationary point found after " + std::to_string(r.steps) + " step" + (r.steps == 1 ? "" : "s") +
             " ***\n";
    } else {
        s += "\n    *** " + r.message + " ***\n";
    }
    s += "\n    Step      Energy (E_h)        dE (E_h)    Max force    RMS force    Max step\n";
    for (const auto& st : r.history)
        s += strprintf("    %4d  %18.10f  %14s  %11.3e  %11.3e  %10.3e\n", st.step, st.energy,
                       st.step == 1 ? "--" : strprintf("%.3e", st.delta_e).c_str(), st.convergence[0].value,
                       st.convergence[1].value, st.convergence[2].value);
    s += strprintf("\n    %s geometry, energy %.10f E_h (%s%.10f E_h relative to the start)\n",
                   r.converged ? "Optimised" : "Last", r.final_energy, r.final_energy - r.initial_energy >= 0 ? "+" : "",
                   r.final_energy - r.initial_energy);
    s += geometry_table(r.final_geometry);
    s += "\n    Optimised parameters (start in parentheses)\n";
    for (std::size_t i = 0; i < r.final_parameters.size(); ++i) {
        const auto& p = r.final_parameters[i];
        const double start = i < r.initial_parameters.size() ? r.initial_parameters[i].value : 0.0;
        if (p.angle)
            s += strprintf("    %-14s = %12.6f deg            (%.6f deg)\n", p.name.c_str(), p.value, start);
        else
            s += strprintf("    %-14s = %12.8f Å   %.8f bohr   (%.6f Å)\n", p.name.c_str(), p.value * kBohrInAngstrom,
                           p.value, start * kBohrInAngstrom);
    }
    if (r.final_geometry.atoms.size() == 3) {
        if (r.symmetric_start)
            s += strprintf("    C2v symmetry %s: |R(O1-H2) - R(O1-H3)| = %.2e bohr\n",
                           r.symmetric_end ? "kept" : "NOT kept", r.symmetry_deviation);
        else
            s += strprintf("    Start without C2v symmetry; at the end |R(O1-H2) - R(O1-H3)| = %.2e bohr\n",
                           r.symmetry_deviation);
    }
    if (r.removed_functions_changed)
        s += "    Warning: the number of removed near-linear dependencies changed during the optimisation.\n";
    s += "    Initial Hessian: " + r.initial_hessian + strprintf("; total time %s\n", duration_string(r.seconds).c_str());
    return s;
}

std::string frequency_block(const Molecule& mol, const HessianResult& h, const HarmonicAnalysis& a) {
    std::string s;
    s += "    Hessian: " + h.description + "\n";
    if (!h.analytic && !h.numerical_reason.empty()) s += "    (numerical because: " + h.numerical_reason + ")\n";
    if (!h.coordinates.empty()) {
        s += "    displaced coordinates:\n";
        for (const auto& c : h.coordinates) s += "      " + c + "\n";
    }
    if (h.analytic) s += strprintf("    CPHF: %d iterations%s\n", h.cphf_iterations, h.cphf_converged ? "" : " (NOT converged)");
    s += "    Dipole derivatives: " + (h.has_dipole_derivatives ? h.dipole_description : std::string("not available")) + "\n";
    s += "\n    Isotopes (AME 2020 masses):";
    for (std::size_t i = 0; i < a.masses.size(); ++i) s += strprintf("  %s %.9f u", a.isotopes[i].c_str(), a.masses[i]);
    s += "\n";
    s += strprintf("    %s molecule: %d vibration%s (3N - %d); translations and rotations projected out (Eckart)\n",
                   a.linear ? "Linear" : "Nonlinear", a.n_vib, a.n_vib == 1 ? "" : "s", a.linear ? 5 : 6);
    s += "    Residual curvature of translations and rotations (cm-1):";
    for (double f : a.low_frequencies) s += strprintf(" %.2f", f);
    s += "\n";
    int imag = 0;
    for (const auto& m : a.modes) imag += m.imaginary ? 1 : 0;
    if (imag) s += strprintf("    Warning: %d imaginary frequenc%s (a saddle point, not a minimum)\n", imag, imag == 1 ? "y" : "ies");
    s += "\n    Mode  Symmetry   Wavenumber (cm-1)   Reduced mass (u)   Force const. (mdyn/Å)   IR intensity (km/mol)\n";
    for (std::size_t k = 0; k < a.modes.size(); ++k) {
        const auto& m = a.modes[k];
        const std::string w = m.imaginary ? strprintf("%11.2fi", std::fabs(m.wavenumber)) : strprintf("%12.2f", m.wavenumber);
        s += strprintf("    %4zu  %-9s %18s %18.6f %22.6f %22s\n", k + 1, m.symmetry.c_str(), w.c_str(), m.reduced_mass,
                       m.force_constant, m.ir_intensity >= 0.0 ? strprintf("%.4f", m.ir_intensity).c_str() : "n/a");
    }
    s += "\n    Normal modes: normalised Cartesian displacements\n";
    for (std::size_t k0 = 0; k0 < a.modes.size(); k0 += 3) {
        const std::size_t k1 = std::min(a.modes.size(), k0 + 3);
        std::string head = "    Atom    ", sym = "            ", freq = "            ";
        for (std::size_t k = k0; k < k1; ++k) {
            head += strprintf("  %-27s", strprintf("mode %zu", k + 1).c_str());
            sym += strprintf("  %-27s", a.modes[k].symmetry.c_str());
            freq += strprintf("  %-27s", (a.modes[k].imaginary ? strprintf("%.2fi cm-1", std::fabs(a.modes[k].wavenumber))
                                                              : strprintf("%.2f cm-1", a.modes[k].wavenumber))
                                             .c_str());
        }
        s += head + "\n" + sym + "\n" + freq + "\n            ";
        for (std::size_t k = k0; k < k1; ++k) s += "        x        y        z  ";
        s += "\n";
        for (std::size_t at = 0; at < mol.atoms.size(); ++at) {
            s += strprintf("    %-7s ", atom_label(mol, at).c_str());
            for (std::size_t k = k0; k < k1; ++k)
                s += strprintf("  %8.4f %8.4f %8.4f ", a.modes[k].displacement[at][0], a.modes[k].displacement[at][1],
                               a.modes[k].displacement[at][2]);
            s += "\n";
        }
    }
    s += strprintf("    Hessian time %s\n", duration_string(h.seconds).c_str());
    return s;
}

std::string thermochemistry_block(const ThermoResult& t) {
    using namespace codata2018;
    std::string s;
    s += strprintf("    Temperature %.2f K, pressure %.2f Pa (%.6f atm, %.5f bar)\n", t.temperature, t.pressure,
                   t.pressure / 101325.0, t.pressure / 1.0e5);
    s += strprintf("    Molecular mass %.6f u; symmetry number %d; %s rigid rotor; closed-shell singlet (q_el = 1)\n",
                   t.mass, t.symmetry_number, t.linear ? "linear" : "nonlinear");
    s += strprintf("    Principal moments of inertia (u bohr^2): %.6f  %.6f  %.6f\n", t.moments[0], t.moments[1], t.moments[2]);
    s += "    Rotational constants (GHz):";
    for (int k = 0; k < 3; ++k)
        if (t.rotational_constants_ghz[static_cast<std::size_t>(k)] > 0) s += strprintf("  %.6f", t.rotational_constants_ghz[static_cast<std::size_t>(k)]);
    s += "\n    Rotational constants (cm-1):";
    for (int k = 0; k < 3; ++k)
        if (t.rotational_constants_cm[static_cast<std::size_t>(k)] > 0) s += strprintf("  %.6f", t.rotational_constants_cm[static_cast<std::size_t>(k)]);
    s += "\n    Rotational temperatures (K):";
    for (int k = 0; k < 3; ++k)
        if (t.rotational_temperatures[static_cast<std::size_t>(k)] > 0) s += strprintf("  %.5f", t.rotational_temperatures[static_cast<std::size_t>(k)]);
    s += "\n    Vibrational temperatures (K):";
    for (double v : t.vibrational_temperatures) s += strprintf("  %.2f", v);
    s += "\n";
    for (const auto& w : t.warnings) s += "    Warning: " + w + "\n";
    const double toH = 1.0 / (kHartree * kAvogadro);
    auto row = [&](const char* name, double jmol) {
        return strprintf("    %-34s %14.8f %12.4f %12.4f\n", name, jmol * toH, jmol / 1000.0 / kCalorie, jmol / 1000.0);
    };
    s += "\n                                         E_h          kcal/mol      kJ/mol\n";
    s += row("Zero-point energy", t.zero_point_energy);
    s += row("Thermal correction to U", t.thermal_energy);
    s += row("Thermal correction to H = U + kT", t.enthalpy);
    s += row("Thermal correction to A = U - TS", t.helmholtz);
    s += row("Thermal correction to G = H - TS", t.gibbs);
    s += strprintf("\n    Electronic energy E_el              %20.10f E_h\n", t.electronic_energy);
    s += strprintf("    E_el + ZPE                          %20.10f E_h\n", t.total_zpe);
    s += strprintf("    U = E_el + thermal energy           %20.10f E_h\n", t.total_u);
    s += strprintf("    H = E_el + thermal enthalpy         %20.10f E_h\n", t.total_h);
    s += strprintf("    A = E_el + thermal Helmholtz energy %20.10f E_h\n", t.total_a);
    s += strprintf("    G = E_el + thermal Gibbs energy     %20.10f E_h\n", t.total_g);
    s += "\n                        E (kcal/mol)  Cv (cal/mol K)   S (cal/mol K)   S (J/mol K)    q\n";
    auto trow = [&](const char* name, const ThermoTerm& x) {
        return strprintf("    %-18s %12.4f %15.4f %15.4f %13.4f   %.6e\n", name, x.energy / 1000.0 / kCalorie,
                         x.cv / kCalorie, x.entropy / kCalorie, x.entropy, x.q);
    };
    s += trow("Electronic", t.electronic);
    s += trow("Translational", t.translation);
    s += trow("Rotational", t.rotation);
    s += trow("Vibrational", t.vibration);
    s += trow("Total", t.total);
    s += "    (vibrational energy including the zero-point energy; vibrational q from the bottom of the well)\n";
    return s;
}

}  // namespace macqc
