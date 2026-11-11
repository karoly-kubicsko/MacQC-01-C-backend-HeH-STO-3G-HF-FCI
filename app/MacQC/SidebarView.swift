//
//  SidebarView.swift
//  MacQC
//
//  Input sidebar: molecule (fixed HeH+), geometry, exponents, basis, options, Run.
//

import SwiftUI

struct SidebarView: View {
    @Bindable var model: AppModel

    private static let numberFormat = FloatingPointFormatStyle<Double>()
        .precision(.fractionLength(0...6))

    var body: some View {
        Form {
            Section("Molecule") {
                LabeledContent("System", value: "HeH⁺")
                LabeledContent("Electrons", value: "2")
                LabeledContent("Charges (Z_A, Z_B)", value: "2, 1")
            }

            Section("Geometry") {
                LabeledContent("R (bohr)") {
                    TextField("R", value: $model.R, format: Self.numberFormat)
                        .labelsHidden()
                        .multilineTextAlignment(.trailing)
                        .monospacedDigit()
                        .frame(width: 96)
                }
            }

            Section("Basis") {
                Picker("Basis set", selection: $model.nGaussians) {
                    Text("STO-1G").tag(1)
                    Text("STO-2G").tag(2)
                    Text("STO-3G").tag(3)
                }
                LabeledContent("ζ (He)") {
                    TextField("zeta He", value: $model.zeta1, format: Self.numberFormat)
                        .labelsHidden()
                        .multilineTextAlignment(.trailing)
                        .monospacedDigit()
                        .frame(width: 96)
                }
                LabeledContent("ζ (H)") {
                    TextField("zeta H", value: $model.zeta2, format: Self.numberFormat)
                        .labelsHidden()
                        .multilineTextAlignment(.trailing)
                        .monospacedDigit()
                        .frame(width: 96)
                }
            }

            Section("Options") {
                Toggle("Exact erf (std::erf)", isOn: $model.exactErf)
                    .help("Use the C++ standard-library erf instead of the Abramowitz–Stegun approximation of the original FORTRAN program.")
            }

            Section {
                Button("Reset to Szabo–Ostlund defaults") {
                    model.resetToDefaults()
                }
                .secondaryActionStyle()
                .accessibilityLabel("Reset to Szabo–Ostlund defaults")
                .frame(maxWidth: .infinity)
                .disabled(model.isRunning)
            }
        }
        .formStyle(.grouped)
        .safeAreaInset(edge: .bottom) {
            runBar
        }
    }

    private var runBar: some View {
        GlassGroup(spacing: 10) {
            VStack(spacing: 8) {
                Button {
                    Task { await model.run() }
                } label: {
                    HStack(spacing: 8) {
                        if model.isRunning {
                            ProgressView()
                                .controlSize(.small)
                        } else {
                            Image(systemName: "play.fill")
                        }
                        Text(model.isRunning ? "Running…" : "Run")
                            .fontWeight(.semibold)
                    }
                    .frame(maxWidth: .infinity)
                }
                .prominentActionStyle()
                .accessibilityLabel("Run")
                .controlSize(.large)
                .keyboardShortcut("r", modifiers: [.command])
                .disabled(model.isRunning)

                if let message = model.errorMessage {
                    Text(message)
                        .font(.caption)
                        .foregroundStyle(.red)
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .lineLimit(3)
                }
            }
            .padding(12)
            .glassPanel(cornerRadius: 18)
        }
        .padding(.horizontal, 12)
        .padding(.bottom, 10)
    }
}
