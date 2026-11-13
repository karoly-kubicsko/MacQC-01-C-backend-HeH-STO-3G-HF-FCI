// The output file of ΨLab 5.0: a complete, human-readable record of a calculation, written
// progressively (streamed to the app's Output view and to psilab's output file) in the spirit
// of a Gaussian .log or ORCA .out file, with ΨLab's own wording and layout.
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/cc.h"
#include "macqc/ci.h"
#include "macqc/frequencies.h"
#include "macqc/gradient.h"
#include "macqc/hessian.h"
#include "macqc/molecule.h"
#include "macqc/optimizer.h"
#include "macqc/rhf.h"

namespace macqc {

// printf into a std::string.
std::string strprintf(const char* fmt, ...);

// Collects the text and passes every piece to an optional sink (streaming).
class OutputWriter {
public:
    using Sink = std::function<void(const std::string&)>;
    explicit OutputWriter(Sink sink = Sink()) : sink_(std::move(sink)) {}
    void write(const std::string& text);
    void line(const std::string& text = "") { write(text + "\n"); }
    // A ruled heading: blank line, "---- Title ----...".
    void section(const std::string& title);
    void subsection(const std::string& title);
    const std::string& text() const { return text_; }

private:
    Sink sink_;
    std::string text_;
};

std::string local_time_string();  // "2026-11-13 20:41:07"
std::string duration_string(double seconds);  // "2.4 s", "3 min 12 s", "1 h 05 min"

// Cartesian coordinates in angstrom and bohr.
std::string geometry_table(const Molecule& mol);
// Bond lengths (angstrom and bohr) and angles (degrees).
std::string parameter_lines(const Molecule& mol, const std::string& indent = "    ");
// Basis set per atom and size.
std::string basis_block(const BasisSet& b, const std::vector<std::string>& notes);

std::string scf_header();
std::string scf_iteration_line(const RHFIteration& it);
// SCF result: convergence, energy components, orbital energies, Mulliken charges, dipole moment.
std::string scf_result_block(const RHFResult& r, const Molecule& mol, const BasisSet& b, bool full_orbitals);
std::string cc_header(const std::string& method);
std::string cc_iteration_line(const CCIteration& it);
std::string fci_header();
std::string fci_iteration_line(const FCIIteration& it);

// Mulliken charges of the atoms for the total density P.
std::vector<double> mulliken_charges(const Molecule& mol, const BasisSet& b, const Matrix& P, const Matrix& S);

// Cartesian gradient table (E_h/bohr) with its largest and RMS component.
std::string gradient_block(const Molecule& mol, const std::vector<double>& g, const std::string& how);

std::string optimization_step_block(const OptimizationStep& s, const ConvergenceCriteria& c);
std::string optimization_result_block(const OptimizationResult& r);

std::string frequency_block(const Molecule& mol, const HessianResult& h, const HarmonicAnalysis& a);
std::string thermochemistry_block(const ThermoResult& t);

}  // namespace macqc
