//
//  SummaryView.swift
//  MacQC
//
//  Summary tab: HF and FCI energies, correlation energy, nuclear repulsion, convergence.
//

import SwiftUI

struct SummaryView: View {
    var model: AppModel

    var body: some View {
        if let out = model.result, let hf = out.hf {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }

                    Text(String(format: "HeH⁺ · STO-%ldG · R = %.4f bohr", out.input.nGaussians, out.input.R))
                        .font(.title2.weight(.semibold))

                    GlassGroup(spacing: 14) {
                        HStack(spacing: 14) {
                            EnergyCard(title: "HF total energy", value: hf.totalEnergy)
                            if let fci = out.fci {
                                EnergyCard(title: "FCI total energy", value: fci.totalEnergy)
                                EnergyCard(title: "Correlation energy", value: fci.correlationEnergy)
                            }
                            EnergyCard(title: "Nuclear repulsion", value: hf.nuclearRepulsion)
                        }
                    }

                    SectionTitle("Hartree–Fock convergence")
                    VStack(spacing: 8) {
                        InfoRow("Converged", hf.converged ? "yes" : "no")
                        InfoRow("SCF iterations", "\(hf.iterations)")
                        if let last = hf.scfHistory.last {
                            InfoRow("Final ΔP (RMS)", Fmt.scientific(last.deltaP))
                        }
                        InfoRow("Convergence threshold", Fmt.scientific(out.input.convergence))
                        InfoRow("HF electronic energy", Fmt.energyWithUnit(hf.electronicEnergy))
                        InfoRow("Error function", out.input.exactErf ? "std::erf (exact)" : "Abramowitz–Stegun (original)")
                        ForEach(Array(hf.orbitalEnergies.enumerated()), id: \.offset) { index, e in
                            InfoRow("Orbital energy ε\(index + 1)", Fmt.energyWithUnit(e))
                        }
                    }
                    .dataPanel()

                    if let fci = out.fci {
                        SectionTitle("Full CI")
                        VStack(spacing: 8) {
                            InfoRow("Ground state", fci.states.first.map { $0.spin.capitalized } ?? "–")
                            InfoRow("FCI electronic energy", Fmt.energyWithUnit(fci.electronicEnergy))
                            InfoRow("Correlation energy", Fmt.energyWithUnit(fci.correlationEnergy))
                            InfoRow("Determinants (M_s = 0)", "\(fci.determinants.count)")
                        }
                        .dataPanel()
                    }
                }
                .padding(20)
            }
        } else {
            NoResultView(model: model)
        }
    }
}

struct EnergyCard: View {
    var title: String
    var value: Double

    var body: some View {
        VStack(alignment: .leading, spacing: 6) {
            Text(title)
                .font(.subheadline)
                .foregroundStyle(.secondary)
            Text(Fmt.energy(value))
                .font(.title3.monospacedDigit().weight(.semibold))
                .textSelection(.enabled)
                .lineLimit(1)
                .minimumScaleFactor(0.7)
            Text("E_h")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .padding(14)
        .glassPanel(cornerRadius: 16)
    }
}

/// A label/value row: label on the left, monospaced-digit value on the right.
struct InfoRow: View {
    var label: String
    var value: String

    init(_ label: String, _ value: String) {
        self.label = label
        self.value = value
    }

    var body: some View {
        HStack(alignment: .firstTextBaseline) {
            Text(label)
                .foregroundStyle(.secondary)
            Spacer(minLength: 20)
            Text(value)
                .monospacedDigit()
                .textSelection(.enabled)
        }
    }
}
