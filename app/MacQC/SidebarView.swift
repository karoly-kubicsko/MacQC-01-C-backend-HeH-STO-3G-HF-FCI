//
//  SidebarView.swift
//  MacQC
//
//  Input sidebar: molecule, job (ΨLab 5.0), geometry (the start of an optimisation, with the
//  optimised values and a button that copies them), method and basis set, the single point of
//  job 5, thermochemistry and optimisation settings, the size of the calculation (each step of
//  the job, analytic or numerical, with warnings and refusals), appearance, and the Run bar.
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
                jobSection
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

    private var jobSection: some View {
        Section("Job") {
            LabeledContent("Job") {
                Menu {
                    Picker("Job", selection: $model.job) {
                        ForEach(JobChoice.allCases) { choice in
                            Text("\(choice.number)   \(choice.title)").tag(choice)
                        }
                    }
                    .pickerStyle(.inline)
                    .labelsHidden()
                } label: {
                    Text(model.job.shortTitle)
                }
                .fixedSize()
                .accessibilityLabel("Job")
            }
            Text(model.job.title + ": " + model.job.explanation)
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .disabled(model.isRunning || model.isScanning)
    }

    @ViewBuilder
    private var generalSections: some View {
        Section(model.job.hasOptimization ? "Starting geometry" : "Geometry") {
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
            if let g = model.optimizedGeometry {
                VStack(alignment: .leading, spacing: 6) {
                    Text("Optimised: " + optimizedText(g))
                        .font(.caption.monospacedDigit())
                        .textSelection(.enabled)
                    Button("Use optimised geometry") {
                        model.useOptimizedGeometry()
                    }
                    .secondaryActionStyle()
                    .accessibilityLabel("Use optimised geometry")
                    .frame(maxWidth: .infinity)
                    .disabled(!model.canUseOptimizedGeometry || model.isRunning)
                }
            }
        }

        Section("Method and basis set") {
            Picker("Method", selection: $model.method) {
                ForEach(MethodChoice.allCases) { choice in
                    Text(choice.title).tag(choice)
                }
            }
            .pickerStyle(.menu)
            Text(model.method.explanation + (model.method == .hf ? "" : " All electrons are correlated."))
                .font(.caption)
                .foregroundStyle(.secondary)
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

        if model.job.hasSinglePoint {
            Section("Single-point energy") {
                Picker("SP method", selection: $model.spMethod) {
                    Text("Same (\(model.method.title))").tag(MethodChoice?.none)
                    ForEach(MethodChoice.allCases) { choice in
                        Text(choice.title).tag(MethodChoice?.some(choice))
                    }
                }
                .pickerStyle(.menu)
                Picker("SP basis set", selection: $model.spBasis) {
                    Text("Same (\(model.basis))").tag(String?.none)
                    ForEach(BasisSets.all, id: \.self) { name in
                        Text(name).tag(String?.some(name))
                    }
                }
                .pickerStyle(.menu)
                Text("Energy at the optimised geometry; E + ZPE, U, H and G add the thermal corrections of the frequency step.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .disabled(model.isRunning || model.isScanning)
        }

        if model.job.hasFrequencies {
            Section("Thermochemistry") {
                LabeledContent("Temperature (K)") {
                    numberField("temperature", value: $model.temperature)
                }
                LabeledContent("Pressure (Pa)") {
                    numberField("pressure", value: $model.pressure)
                }
                HStack {
                    Button("1 atm") { model.pressure = 101325 }
                    Button("1 bar") { model.pressure = 100000 }
                    Spacer()
                }
                .controlSize(.small)
                Text("Ideal gas, rigid rotor, harmonic oscillators; CODATA 2018 constants, AME 2020 isotope masses.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .disabled(model.isRunning || model.isScanning)
        }

        if model.job.hasOptimization {
            Section("Optimisation") {
                Picker("Convergence", selection: $model.tightConvergence) {
                    Text("Normal").tag(false)
                    Text("Tight").tag(true)
                }
                .pickerStyle(.segmented)
                Stepper(value: $model.maxSteps, in: 5...500, step: 5) {
                    LabeledContent("Maximum steps", value: "\(model.maxSteps)")
                }
                Text(model.tightConvergence
                     ? "Max force 1.5e-5, RMS force 1.0e-5 E_h/bohr; max step 6.0e-5, RMS step 4.0e-5 bohr; energy change 1e-8 E_h."
                     : "Max force 4.5e-4, RMS force 3.0e-4 E_h/bohr; max step 1.8e-3, RMS step 1.2e-3 bohr; energy change 1e-6 E_h.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .disabled(model.isRunning || model.isScanning)
        }

        Section("Size of the calculation") {
            EstimatePanel(model: model)
        }
    }

    private func optimizedText(_ g: (R: Double, angle: Double?, symmetric: Bool)) -> String {
        var text = String(format: "%@ = %.6f Å", model.molecule.bondLabel, g.R * Units.bohrToAngstrom)
        if let a = g.angle {
            text += String(format: ", angle %.4f°", a)
        }
        if !g.symmetric {
            text += " (not C₂ᵥ: mean of both O–H)"
        }
        return text
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
                        Text(model.progressDescription + "  ·  " + model.elapsedText(at: context.date))
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
                          let seconds = model.estimatedSeconds {
                    Label("Long calculation: \(Fmt.duration(seconds)).", systemImage: "clock.badge.exclamationmark")
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
                    if let c = e.coupledCluster {
                        InfoRow("\(MethodChoice.title(of: c.method)) amplitudes", Fmt.count(c.amplitudes))
                        InfoRow("Orbitals (occ./virt.)", "\(c.occupied) / \(c.virtuals)")
                        InfoRow("\(MethodChoice.title(of: c.method)) memory", Fmt.memory(c.memoryMB))
                    }
                    if let j = out.job, j.type != "energy" {
                        ForEach(j.estimate.steps) { step in
                            VStack(alignment: .leading, spacing: 2) {
                                InfoRow(step.name, model.refusal == nil ? Fmt.duration(step.seconds) : "")
                                Text("\(MethodChoice.title(of: step.method))/\(step.basis): \(step.how)")
                                    .font(.caption)
                                    .foregroundStyle(step.analytic ? Color.secondary : Color.orange)
                            }
                        }
                        if model.refusal == nil {
                            InfoRow("Estimated time (job)", Fmt.duration(j.estimate.totalSeconds) + " (\(e.threads) threads)")
                        }
                    } else if model.refusal == nil {
                        InfoRow("Estimated time", Fmt.duration(e.totalSeconds) + " (\(e.threads) threads)")
                    }
                }
                ForEach(out.notes ?? [], id: \.self) { note in
                    Text(note)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                ForEach(warnings(out), id: \.self) { warning in
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

    /// Warnings of the energy and, for jobs, of the whole job (its "long calculation" warning
    /// replaces the one of a single energy).
    private func warnings(_ out: MacQCOutput) -> [String] {
        guard let j = out.job, j.type != "energy" else { return out.warnings ?? [] }
        let own = (out.warnings ?? []).filter { !$0.hasPrefix("Long calculation") }
        return own + j.estimate.warnings
    }
}
