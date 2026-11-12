//
//  FullCIView.swift
//  MacQC
//
//  Full CI tab: determinants and CI vector, singlet configuration weights, all states.
//

import SwiftUI

struct FullCIView: View {
    var model: AppModel

    var body: some View {
        if let out = model.result, out.isGeneral {
            GeneralFullCIView(model: model, out: out)
        } else if let out = model.result, let fci = out.fci {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }

                    SectionTitle("Ground-state CI vector")
                    determinantTable(fci)

                    SectionTitle("Singlet configuration weights (ground state)")
                    weightsPanel(fci)

                    SectionTitle("All M_s = 0 states")
                    statesTable(fci)
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
                    description: Text("Full CI needs a converged Hartree–Fock calculation. Press Run.")
                )
            }
            .padding()
        }
    }

    private func determinantTable(_ fci: MacQCFCI) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 32, verticalSpacing: 5) {
            GridRow {
                Text("Determinant")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                    .gridColumnAlignment(.leading)
                Text("Coefficient c")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                Text("Weight c²")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
            }
            Divider()
                .gridCellUnsizedAxes(.horizontal)
            ForEach(Array(fci.determinants.enumerated()), id: \.offset) { index, label in
                let c = index < fci.groundStateVector.count ? fci.groundStateVector[index] : 0.0
                GridRow {
                    Text(label)
                    Text(Fmt.value8(c))
                    Text(Fmt.value8(c * c))
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }

    private func weightsPanel(_ fci: MacQCFCI) -> some View {
        VStack(spacing: 12) {
            ForEach(Array(fci.singletWeights.enumerated()), id: \.offset) { _, item in
                VStack(alignment: .leading, spacing: 4) {
                    HStack {
                        Text(displayName(item.configuration))
                        Spacer()
                        Text(Fmt.value8(item.weight))
                            .monospacedDigit()
                            .textSelection(.enabled)
                    }
                    ProgressView(value: min(max(item.weight, 0.0), 1.0))
                        .progressViewStyle(.linear)
                }
            }
        }
        .dataPanel()
    }

    private func statesTable(_ fci: MacQCFCI) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 28, verticalSpacing: 5) {
            GridRow {
                Text("State")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                    .gridColumnAlignment(.leading)
                Text("Total energy (E_h)")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                Text("⟨S²⟩")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                Text("Spin")
                    .font(.subheadline.weight(.semibold))
                    .foregroundStyle(.secondary)
                    .gridColumnAlignment(.leading)
            }
            Divider()
                .gridCellUnsizedAxes(.horizontal)
            ForEach(Array(fci.states.enumerated()), id: \.offset) { index, state in
                GridRow {
                    Text("\(index + 1)")
                    Text(Fmt.energy(state.energy))
                    Text(String(format: "%.4f", state.S2))
                    Text(state.spin.capitalized)
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }

    /// "(1sigma)^2" -> "(1σ)²", "(1sigma)^1(2sigma)^1" -> "(1σ)¹(2σ)¹"
    private func displayName(_ configuration: String) -> String {
        configuration
            .replacingOccurrences(of: "sigma", with: "σ")
            .replacingOccurrences(of: ")^1", with: ")¹")
            .replacingOccurrences(of: ")^2", with: ")²")
    }
}
