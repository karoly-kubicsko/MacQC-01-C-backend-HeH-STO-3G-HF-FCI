//
//  HartreeFockView.swift
//  MacQC
//
//  Hartree–Fock tab: integrals, S/H/F/C/P matrices, orbital energies, SCF iterations.
//  Tables of numbers use plain (non-glass) panels.
//

import SwiftUI

struct HartreeFockView: View {
    var model: AppModel

    var body: some View {
        if let out = model.result, let hf = out.hf {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }

                    SectionTitle("Orbital energies")
                    orbitalEnergies(hf)

                    SectionTitle("Integrals (atomic units)")
                    integralsTable(hf.integrals)

                    SectionTitle("Matrices")
                    LazyVGrid(columns: [GridItem(.adaptive(minimum: 230), spacing: 16, alignment: .top)],
                              alignment: .leading, spacing: 16) {
                        MatrixPanel(title: "S – overlap", rows: hf.S)
                        MatrixPanel(title: "H – core Hamiltonian", rows: hf.H)
                        MatrixPanel(title: "F – Fock", rows: hf.F)
                        MatrixPanel(title: "C – MO coefficients", rows: hf.C)
                        MatrixPanel(title: "P – density", rows: hf.P)
                    }

                    SectionTitle("SCF iterations")
                    scfTable(hf)
                }
                .padding(20)
            }
        } else {
            NoResultView(model: model)
        }
    }

    private func orbitalEnergies(_ hf: MacQCHF) -> some View {
        VStack(spacing: 8) {
            ForEach(Array(hf.orbitalEnergies.enumerated()), id: \.offset) { index, e in
                InfoRow("ε\(index + 1)  (occupation \(index == 0 ? "2" : "0"))", Fmt.energyWithUnit(e))
            }
        }
        .dataPanel()
    }

    private func integralsTable(_ integrals: MacQCIntegrals) -> some View {
        let items = integrals.namedValues
        return LazyVGrid(
            columns: [GridItem(.adaptive(minimum: 170), spacing: 18, alignment: .leading)],
            alignment: .leading, spacing: 6
        ) {
            ForEach(Array(items.enumerated()), id: \.offset) { _, item in
                HStack {
                    Text(item.name)
                        .foregroundStyle(.secondary)
                    Spacer(minLength: 8)
                    Text(Fmt.value6(item.value))
                        .monospacedDigit()
                        .textSelection(.enabled)
                }
            }
        }
        .dataPanel()
    }

    private func scfTable(_ hf: MacQCHF) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 28, verticalSpacing: 5) {
            GridRow {
                Text("Iteration")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                    .gridColumnAlignment(.leading)
                Text("Electronic energy (E_h)")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                Text("ΔP (RMS)")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
            }

            Divider()
                .gridCellUnsizedAxes(.horizontal)

            ForEach(hf.scfHistory) { step in
                GridRow {
                    Text("\(step.iteration)")
                    Text(Fmt.energy(step.electronicEnergy))
                    Text(Fmt.scientific(step.deltaP))
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }
}

/// A titled 2×2 (or n×m) matrix with monospaced digits.
struct MatrixPanel: View {
    var title: String
    var rows: [[Double]]

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            Text(title)
                .font(.subheadline.weight(.semibold))
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                ForEach(Array(rows.enumerated()), id: \.offset) { _, row in
                    GridRow {
                        ForEach(Array(row.enumerated()), id: \.offset) { _, value in
                            Text(Fmt.value6(value))
                        }
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .dataPanel()
    }
}
