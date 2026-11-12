//
//  GeneralViews.swift
//  MacQC
//
//  Detail tabs for H₂, LiH and H₂O: summary, Hartree–Fock, Full CI, and the text report
//  (shared with HeH⁺).
//

import SwiftUI
import Charts
import Foundation

// MARK: - Summary

struct GeneralSummaryView: View {
    var model: AppModel
    var out: MacQCOutput

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 18) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                }
                header

                if let rhf = out.rhf {
                    GlassGroup(spacing: 14) {
                        HStack(spacing: 14) {
                            EnergyCard(title: "RHF total energy", value: rhf.totalEnergy, decimals: 10)
                            if let fci = out.fullCI {
                                EnergyCard(title: "Full CI total energy", value: fci.totalEnergy, decimals: 10)
                                EnergyCard(title: "Correlation energy", value: fci.correlationEnergy, decimals: 10)
                            }
                            EnergyCard(title: "Nuclear repulsion", value: rhf.nuclearRepulsion, decimals: 10)
                        }
                    }
                }

                NotesPanel(notes: out.notes ?? [], warnings: out.warnings ?? [])

                if let mol = out.molecule {
                    SectionTitle("Molecule")
                    moleculePanel(mol)
                }
                if let b = out.basis {
                    SectionTitle("Basis set")
                    basisPanel(b)
                }
                if let rhf = out.rhf {
                    SectionTitle("Hartree–Fock")
                    VStack(spacing: 8) {
                        InfoRow("Converged", rhf.converged ? "yes (\(rhf.iterations) iterations)" : "no")
                        if let homo = rhf.homo {
                            InfoRow("HOMO", Fmt.energyWithUnit(homo))
                        }
                        if let lumo = rhf.lumo {
                            InfoRow("LUMO", Fmt.energyWithUnit(lumo))
                        }
                        if let homo = rhf.homo, let lumo = rhf.lumo {
                            InfoRow("HOMO–LUMO gap", Fmt.energyWithUnit(lumo - homo))
                        }
                        InfoRow("Virial ratio −V/T", String(format: "%.8f", rhf.virialRatio))
                        InfoRow("Two-electron integrals",
                                Fmt.count(rhf.uniqueIntegrals) + (rhf.integralsStored ? " (kept in memory)" : " (integral-direct)"))
                        InfoRow("Time", String(format: "%.2f s (integrals %.2f s, SCF %.2f s, %ld threads)",
                                               rhf.integralSeconds + rhf.scfSeconds, rhf.integralSeconds,
                                               rhf.scfSeconds, rhf.threads))
                    }
                    .dataPanel()
                }
                if let f = out.fullCI {
                    SectionTitle("Full CI")
                    VStack(spacing: 8) {
                        InfoRow("Determinants", Fmt.count(f.determinants) + String(
                            format: " (%@ α × %@ β strings)", Fmt.count(f.alphaStrings), Fmt.count(f.betaStrings)))
                        InfoRow("Electrons in orbitals", "\(f.alphaElectrons + f.betaElectrons) in \(f.orbitals)")
                        InfoRow("Converged", f.converged ? "yes (\(f.iterations) Davidson iterations)" : "no")
                        InfoRow("⟨S²⟩", String(format: "%.8f", f.S2))
                        InfoRow("Weight of the RHF determinant |C₀|²", String(format: "%.8f", f.hfWeight))
                        InfoRow("Time", String(format: "%.2f s (integral transformation %.2f s, CI %.2f s)",
                                               f.transformSeconds + f.ciSeconds, f.transformSeconds, f.ciSeconds))
                    }
                    .dataPanel()
                }
            }
            .padding(20)
        }
    }

    private var header: some View {
        VStack(alignment: .leading, spacing: 4) {
            let name = Fmt.formula(out.molecule?.name ?? out.input.molecule)
            let method = out.fullCI != nil ? "RHF and Full CI" : "RHF"
            Text("\(name) · \(out.basis?.name ?? out.input.basis) · \(method)")
                .font(.title2.weight(.semibold))
            if let mol = out.molecule {
                Text(geometryText(mol))
                    .foregroundStyle(.secondary)
            }
        }
    }

    private func geometryText(_ mol: MacQCMoleculeInfo) -> String {
        let label = MoleculeChoice(rawValue: mol.name)?.bondLabel ?? "R"
        var text = "\(label) = " + Fmt.lengthBoth(mol.R)
        if let angle = mol.angle {
            text += String(format: ", angle H–O–H = %.2f°", angle)
        }
        return text
    }

    private func moleculePanel(_ mol: MacQCMoleculeInfo) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 24, verticalSpacing: 5) {
            GridRow {
                Text("Atom").gridColumnAlignment(.leading)
                Text("x (bohr)")
                Text("y (bohr)")
                Text("z (bohr)")
            }
            .font(.subheadline.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider().gridCellUnsizedAxes(.horizontal)
            ForEach(Array(mol.atoms.enumerated()), id: \.offset) { _, atom in
                GridRow {
                    Text(atom.symbol)
                    ForEach(0..<3, id: \.self) { k in
                        Text(String(format: "%.6f", k < atom.position.count ? atom.position[k] : 0))
                    }
                }
            }
            Divider().gridCellUnsizedAxes(.horizontal)
            GridRow {
                Text("\(mol.electrons) electrons, charge \(mol.charge), singlet")
                    .gridCellColumns(2)
                    .gridColumnAlignment(.leading)
                Text("E_nn = " + Fmt.energy10(mol.nuclearRepulsion) + " E_h")
                    .gridCellColumns(2)
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }

    private func basisPanel(_ b: MacQCBasisInfo) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            Grid(alignment: .leading, horizontalSpacing: 24, verticalSpacing: 5) {
                GridRow {
                    Text("Atom")
                    Text("Basis set")
                    Text("Primitives")
                    Text("Contracted")
                    Text("Functions").gridColumnAlignment(.trailing)
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(Fmt.uniqueAtoms(b.atoms), id: \.symbol) { atom in
                    GridRow {
                        Text(atom.symbol)
                        Text(atom.basis)
                        Text(atom.primitives)
                        Text(atom.contraction)
                        Text("\(atom.functions)")
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text("\(b.functions) spherical basis functions in \(b.shells) shells; highest angular momentum l = \(b.maxL).")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }
}

/// Notes (secondary) and warnings (orange) from the backend.
struct NotesPanel: View {
    var notes: [String]
    var warnings: [String]

    var body: some View {
        if !notes.isEmpty || !warnings.isEmpty {
            VStack(alignment: .leading, spacing: 6) {
                ForEach(warnings, id: \.self) { w in
                    Label(w, systemImage: "exclamationmark.triangle.fill")
                        .foregroundStyle(.orange)
                }
                ForEach(notes, id: \.self) { n in
                    Label(n, systemImage: "info.circle")
                        .foregroundStyle(.secondary)
                }
            }
            .font(.callout)
            .frame(maxWidth: .infinity, alignment: .leading)
            .dataPanel()
        }
    }
}

// MARK: - Hartree–Fock

struct GeneralHartreeFockView: View {
    var model: AppModel
    var out: MacQCOutput

    var body: some View {
        if let rhf = out.rhf {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }
                    SectionTitle("Energy components (E_h)")
                    VStack(spacing: 8) {
                        InfoRow("Total energy", Fmt.energy10(rhf.totalEnergy))
                        InfoRow("Electronic energy", Fmt.energy10(rhf.electronicEnergy))
                        InfoRow("Nuclear repulsion", Fmt.energy10(rhf.nuclearRepulsion))
                        InfoRow("Kinetic energy T", Fmt.energy10(rhf.kineticEnergy))
                        InfoRow("Nuclear attraction", Fmt.energy10(rhf.nuclearAttractionEnergy))
                        InfoRow("Electron repulsion", Fmt.energy10(rhf.electronRepulsionEnergy))
                        InfoRow("Virial ratio −V/T", String(format: "%.10f", rhf.virialRatio))
                        InfoRow("Smallest overlap eigenvalue", Fmt.scientific(rhf.smallestOverlapEigenvalue)
                                + (rhf.removedFunctions > 0 ? " (\(rhf.removedFunctions) removed)" : ""))
                    }
                    .dataPanel()

                    SectionTitle("Orbital energies (E_h)")
                    orbitalGrid(rhf)

                    SectionTitle("SCF iterations")
                    scfTable(rhf)
                }
                .padding(20)
            }
        } else {
            NoResultView(model: model)
        }
    }

    private func orbitalGrid(_ rhf: MacQCRHF) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            LazyVGrid(columns: [GridItem(.adaptive(minimum: 190), spacing: 18, alignment: .leading)],
                      alignment: .leading, spacing: 4) {
                ForEach(Array(rhf.orbitalEnergies.enumerated()), id: \.offset) { index, e in
                    HStack {
                        Text("\(index + 1)\(index < rhf.occupiedOrbitals ? " (occ)" : "")")
                            .foregroundStyle(index < rhf.occupiedOrbitals ? Color.primary : Color.secondary)
                        Spacer(minLength: 8)
                        Text(Fmt.energy(e))
                            .monospacedDigit()
                            .textSelection(.enabled)
                    }
                }
            }
            Text("\(rhf.occupiedOrbitals) doubly occupied orbitals, \(rhf.orbitalEnergies.count - rhf.occupiedOrbitals) virtual orbitals.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func scfTable(_ rhf: MacQCRHF) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 28, verticalSpacing: 5) {
            GridRow {
                Text("Iteration").gridColumnAlignment(.leading)
                Text("Total energy (E_h)")
                Text("ΔE (E_h)")
                Text("Orbital gradient")
                Text("DIIS").gridColumnAlignment(.leading)
            }
            .font(.subheadline.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider().gridCellUnsizedAxes(.horizontal)
            ForEach(rhf.scfHistory) { step in
                GridRow {
                    Text("\(step.iteration)")
                    Text(Fmt.energy10(step.totalEnergy))
                    Text(Fmt.scientific(step.deltaE))
                    Text(Fmt.scientific(step.gradient))
                    Text(step.diis ? "yes" : "")
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }
}

// MARK: - Full CI

struct GeneralFullCIView: View {
    var model: AppModel
    var out: MacQCOutput

    var body: some View {
        if let f = out.fullCI {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }
                    SectionTitle("Energies (E_h)")
                    VStack(spacing: 8) {
                        if let rhf = out.rhf {
                            InfoRow("RHF energy", Fmt.energy10(rhf.totalEnergy))
                        }
                        InfoRow("Full CI energy", Fmt.energy10(f.totalEnergy))
                        InfoRow("Correlation energy", Fmt.energy10(f.correlationEnergy))
                        InfoRow("⟨S²⟩", String(format: "%.8f", f.S2))
                        InfoRow("Determinants (M_S = 0)", Fmt.count(f.determinants))
                    }
                    .dataPanel()

                    SectionTitle("Leading determinants")
                    determinantTable(f)

                    SectionTitle("Natural-orbital occupation numbers")
                    occupationChart(f)

                    SectionTitle("Davidson iterations")
                    davidsonTable(f)
                }
                .padding(20)
            }
        } else {
            VStack(spacing: 14) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                        .frame(maxWidth: 560)
                }
                ContentUnavailableView(
                    "No Full CI result",
                    systemImage: "function",
                    description: Text(placeholderText)
                )
            }
            .padding()
        }
    }

    private var placeholderText: String {
        if let refusal = model.refusal, model.method == .fci {
            return refusal
        }
        if out.input.method != "FCI" {
            return "This result is RHF only. Choose Full CI in the sidebar and press Run (H₂ and LiH up to aug-cc-pCVTZ; for H₂O Full CI is impossible)."
        }
        return "Full CI needs a converged Hartree–Fock calculation."
    }

    private func orbitals(_ list: [Int]) -> String {
        list.map(String.init).joined(separator: " ")
    }

    private func determinantTable(_ f: MacQCFullCI) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            Grid(alignment: .trailing, horizontalSpacing: 26, verticalSpacing: 5) {
                GridRow {
                    Text("Coefficient")
                    Text("Weight c²")
                    Text("α orbitals").gridColumnAlignment(.leading)
                    Text("β orbitals").gridColumnAlignment(.leading)
                    Text("Excitation")
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(Array(f.leadingDeterminants.enumerated()), id: \.offset) { _, d in
                    GridRow {
                        Text(String(format: "%+.8f", d.coefficient))
                        Text(Fmt.value8(d.coefficient * d.coefficient))
                        Text(orbitals(d.alpha))
                        Text(orbitals(d.beta))
                        Text(d.excitation == 0 ? "RHF" : "\(d.excitation)")
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text("Occupied molecular orbitals (numbered by RHF orbital energy) of the α and β spin strings.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    /// Points and line of log10(occupation) on a linear axis (a logarithmic axis cannot hold the
    /// zero baseline of bars), strongly occupied orbitals in blue, weakly occupied in orange.
    private func occupationChart(_ f: MacQCFullCI) -> some View {
        let values = Array(f.naturalOccupations.prefix(20))
        let points: [(index: Int, log: Double, strong: Bool)] = values.enumerated().map { item in
            (index: item.offset + 1, log: log10(max(item.element, 1.0e-12)), strong: item.element > 0.5)
        }
        let lowest = (points.map { $0.log }.min() ?? -6).rounded(.down)
        return VStack(alignment: .leading, spacing: 8) {
            Chart {
                ForEach(points, id: \.index) { p in
                    LineMark(x: .value("Natural orbital", p.index), y: .value("log10 occupation", p.log))
                        .foregroundStyle(Color.secondary.opacity(0.5))
                    PointMark(x: .value("Natural orbital", p.index), y: .value("log10 occupation", p.log))
                        .foregroundStyle(p.strong ? Color.blue : Color.orange)
                        .symbolSize(60)
                }
            }
            .chartYScale(domain: min(lowest, -1.0)...0.5)
            .chartXScale(domain: 0...(points.count + 1))
            .chartXAxisLabel("Natural orbital")
            .chartYAxisLabel("log₁₀ occupation")
            .frame(height: 220)
            Text(values.map { String(format: "%.6f", $0) }.joined(separator: "  "))
                .font(.caption.monospacedDigit())
                .textSelection(.enabled)
            Text("Largest \(values.count) of \(f.naturalOccupations.count); they sum to \(f.alphaElectrons + f.betaElectrons) electrons. Blue: occupation above 0.5, orange: weakly occupied (correlation).")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func davidsonTable(_ f: MacQCFullCI) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 28, verticalSpacing: 5) {
            GridRow {
                Text("Iteration").gridColumnAlignment(.leading)
                Text("Total energy (E_h)")
                Text("Residual")
                Text("Subspace")
            }
            .font(.subheadline.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider().gridCellUnsizedAxes(.horizontal)
            ForEach(f.history) { step in
                GridRow {
                    Text("\(step.iteration)")
                    Text(Fmt.energy10(step.totalEnergy))
                    Text(Fmt.scientific(step.residual))
                    Text("\(step.subspace)")
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }
}

// MARK: - Report

/// The backend's text report (FORTRAN-style for HeH⁺, the general report otherwise).
struct ReportView: View {
    var model: AppModel

    var body: some View {
        if let out = model.result, !out.report.isEmpty {
            ScrollView([.vertical, .horizontal]) {
                Text(out.report)
                    .font(.system(.callout, design: .monospaced))
                    .textSelection(.enabled)
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(20)
            }
        } else {
            NoResultView(model: model)
        }
    }
}
