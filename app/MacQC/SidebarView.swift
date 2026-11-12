//
//  SidebarView.swift
//  MacQC
//
//  Input sidebar: molecule, geometry, basis set and method, the size of the calculation
//  (with warnings and refusals), appearance, and the Run bar.
//

import SwiftUI

struct SidebarView: View {
    @Bindable var model: AppModel
    @AppStorage(AppearanceChoice.storageKey) private var appearanceRaw: String = AppearanceChoice.system.rawValue

    private static let numberFormat = FloatingPointFormatStyle<Double>()
        .precision(.fractionLength(0...6))

    var body: some View {
        Form {
            Section("Molecule") {
                Picker("Molecule", selection: $model.molecule) {
                    ForEach(MoleculeChoice.allCases) { choice in
                        Text(choice.title).tag(choice)
                    }
                }
                .disabled(model.isRunning || model.isScanning)
                Text(model.molecule.subtitle)
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }

            if model.isGeneral {
                generalSections
            } else {
                hehSections
            }

            Section("Appearance") {
                Picker("Appearance", selection: $appearanceRaw) {
                    ForEach(AppearanceChoice.allCases) { choice in
                        Text(choice.title).tag(choice.rawValue)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
            }
        }
        .formStyle(.grouped)
        .safeAreaInset(edge: .bottom) {
            runBar
        }
        .task(id: model.estimateKey) {
            // Size and time of the calculation, refreshed shortly after an input changes.
            try? await Task.sleep(for: .milliseconds(200))
            guard !Task.isCancelled else { return }
            await model.refreshEstimate()
        }
        .confirmationDialog("Long calculation", isPresented: $model.showLongRunConfirmation) {
            Button("Run") {
                Task { await model.run() }
            }
            Button("Cancel", role: .cancel) {}
        } message: {
            Text(model.longRunMessage + "\n\nYou can cancel it while it runs.")
        }
    }

    // MARK: HeH+ (Szabo–Ostlund)

    @ViewBuilder
    private var hehSections: some View {
        Section("Geometry") {
            LabeledContent("R (Å)") {
                numberField("R angstrom", value: $model.rAngstrom)
            }
            LabeledContent("R (bohr)") {
                numberField("R bohr", value: $model.R)
            }
        }

        Section("Basis") {
            Picker("Basis set", selection: $model.nGaussians) {
                Text("STO-1G").tag(1)
                Text("STO-2G").tag(2)
                Text("STO-3G").tag(3)
            }
            LabeledContent("ζ (He)") {
                numberField("zeta He", value: $model.zeta1)
            }
            LabeledContent("ζ (H)") {
                numberField("zeta H", value: $model.zeta2)
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

    // MARK: H2, LiH, H2O

    @ViewBuilder
    private var generalSections: some View {
        Section("Geometry") {
            LabeledContent("\(model.molecule.bondLabel) (Å)") {
                numberField("bond length angstrom", value: $model.bondLengthAngstrom)
            }
            LabeledContent("\(model.molecule.bondLabel) (bohr)") {
                numberField("bond length bohr", value: $model.bondLength)
            }
            if model.molecule == .h2o {
                LabeledContent("Angle H–O–H (°)") {
                    numberField("angle", value: $model.waterAngle)
                }
            }
            Button("Experimental geometry") {
                model.resetGeometry()
            }
            .secondaryActionStyle()
            .accessibilityLabel("Experimental geometry")
            .frame(maxWidth: .infinity)
            .disabled(model.isRunning)
        }

        Section("Method and basis set") {
            Picker("Method", selection: $model.method) {
                ForEach(MethodChoice.allCases) { choice in
                    Text(choice.title).tag(choice)
                }
            }
            .pickerStyle(.segmented)
            Picker("Basis set", selection: $model.basis) {
                ForEach(BasisSets.all, id: \.self) { name in
                    Text(name).tag(name)
                }
            }
            Text("Spherical Gaussian functions from basissetexchange.org; the core-valence sets (CV) add tight functions for the 1s core.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .disabled(model.isRunning || model.isScanning)

        Section("Size of the calculation") {
            EstimatePanel(model: model)
        }
    }

    private func numberField(_ title: String, value: Binding<Double>) -> some View {
        TextField(title, value: value, format: Self.numberFormat)
            .labelsHidden()
            .multilineTextAlignment(.trailing)
            .monospacedDigit()
            .frame(width: 96)
    }

    // MARK: Run bar

    private var runBar: some View {
        GlassGroup(spacing: 10) {
            VStack(spacing: 8) {
                HStack(spacing: 8) {
                    Button {
                        model.requestRun()
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
                    .disabled(model.isRunning || model.isScanning || model.refusal != nil)

                    if model.isRunning {
                        Button("Cancel") {
                            model.cancelRun()
                        }
                        .secondaryActionStyle()
                        .controlSize(.large)
                        .accessibilityLabel("Cancel calculation")
                        .keyboardShortcut(".", modifiers: [.command])
                    }
                }

                if model.isRunning && model.isGeneral {
                    TimelineView(.periodic(from: .now, by: 1)) { context in
                        Text(model.progress.text + "  ·  " + model.elapsedText(at: context.date))
                            .font(.caption.monospacedDigit())
                            .foregroundStyle(.secondary)
                            .frame(maxWidth: .infinity, alignment: .leading)
                            .lineLimit(2)
                    }
                }

                if let message = model.errorMessage {
                    Text(message)
                        .font(.caption)
                        .foregroundStyle(.red)
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .lineLimit(4)
                } else if model.isGeneral && !model.isRunning, model.refusal != nil {
                    Label("Not possible with these settings (see Size of the calculation).", systemImage: "xmark.octagon.fill")
                        .font(.caption)
                        .foregroundStyle(.red)
                        .frame(maxWidth: .infinity, alignment: .leading)
                } else if model.isGeneral && !model.isRunning && model.isLongCalculation,
                          let e = model.estimate?.estimate {
                    Label("Long calculation: \(Fmt.duration(e.totalSeconds)).", systemImage: "clock.badge.exclamationmark")
                        .font(.caption)
                        .foregroundStyle(.orange)
                        .frame(maxWidth: .infinity, alignment: .leading)
                }
            }
            .padding(12)
            .glassPanel(cornerRadius: 18)
        }
        .padding(.horizontal, 12)
        .padding(.bottom, 10)
    }
}

/// Basis functions, integrals, determinants, memory and time of the selected calculation,
/// with the backend's notes, warnings and refusal.
struct EstimatePanel: View {
    var model: AppModel

    var body: some View {
        if let out = model.estimate {
            VStack(alignment: .leading, spacing: 6) {
                if let b = out.basis {
                    InfoRow("Basis functions", "\(b.functions)")
                    ForEach(Fmt.uniqueAtoms(b.atoms), id: \.symbol) { atom in
                        InfoRow("  \(atom.symbol)", (atom.basis == model.basis ? "" : atom.basis + " ") + atom.contraction)
                    }
                }
                if let e = out.estimate {
                    InfoRow("Two-electron integrals",
                            Fmt.count(e.uniqueIntegrals) + (e.integralsStored ? " (in memory)" : " (direct)"))
                    if let f = e.fci {
                        InfoRow("Determinants", Fmt.count(f.determinants))
                        if f.feasible {
                            InfoRow("Full CI memory", Fmt.memory(f.memoryMB))
                        }
                    }
                    if model.refusal == nil {
                        InfoRow("Estimated time", Fmt.duration(e.totalSeconds) + " (\(e.threads) threads)")
                    }
                }
                ForEach(out.notes ?? [], id: \.self) { note in
                    Text(note)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                ForEach(out.warnings ?? [], id: \.self) { warning in
                    Label(warning, systemImage: "clock.badge.exclamationmark")
                        .font(.caption)
                        .foregroundStyle(.orange)
                }
                if let refusal = model.refusal {
                    Label(refusal, systemImage: "xmark.octagon.fill")
                        .font(.caption)
                        .foregroundStyle(.red)
                        .textSelection(.enabled)
                }
            }
        } else {
            HStack(spacing: 8) {
                ProgressView()
                    .controlSize(.small)
                Text("Estimating…")
                    .foregroundStyle(.secondary)
            }
        }
    }
}
