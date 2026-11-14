//
//  JobViews.swift
//  MacQC
//
//  ΨLab 5.0 job results: the Optimisation tab (steps table, energy and force charts, the
//  geometry and internal coordinates of each step, the optimised geometry with a button that
//  copies it into the inputs), the Frequencies tab with its Frequencies and Thermochemistry
//  pages (normal modes with an animated drawing, IR spectrum), and the job panel of the
//  Summary tab.
//

import SwiftUI
import Charts

// MARK: - Optimisation

struct OptimizationView: View {
    @Bindable var model: AppModel
    @State private var selectedStep: Int?

    var body: some View {
        if let out = model.result, let opt = out.optimization {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }
                    header(out, opt)
                    resultPanel(out, opt)

                    SectionTitle("Steps")
                    stepsTable(opt)
                    charts(opt)

                    if let step = shownStep(opt) {
                        SectionTitle("Step \(step.step)")
                        stepDetail(step, opt: opt)
                    }
                }
                .padding(20)
            }
        } else {
            JobPlaceholder(model: model, title: "No optimisation", systemImage: "arrow.down.right.and.arrow.up.left",
                           text: "Choose the job Geometry optimisation (or Optimisation + frequencies) in the sidebar and press Run. The geometry inputs are the starting geometry.")
        }
    }

    private func shownStep(_ opt: MacQCOptimization) -> MacQCOptimizationStep? {
        if let s = selectedStep, let step = opt.history.first(where: { $0.step == s }) { return step }
        return opt.history.last
    }

    private func header(_ out: MacQCOutput, _ opt: MacQCOptimization) -> some View {
        VStack(alignment: .leading, spacing: 4) {
            Text("\(Fmt.formula(out.input.molecule)) · \(MethodChoice.title(of: out.input.method))/\(out.input.basis) · geometry optimisation")
                .font(.title2.weight(.semibold))
            Label(opt.converged ? "Stationary point found after \(opt.steps) steps" : opt.message,
                  systemImage: opt.converged ? "checkmark.seal.fill" : "exclamationmark.triangle.fill")
                .foregroundStyle(opt.converged ? Color.green : Color.orange)
        }
    }

    private func resultPanel(_ out: MacQCOutput, _ opt: MacQCOptimization) -> some View {
        VStack(alignment: .leading, spacing: 10) {
            Grid(alignment: .trailing, horizontalSpacing: 24, verticalSpacing: 5) {
                GridRow {
                    Text("Parameter").gridColumnAlignment(.leading)
                    Text("Start")
                    Text(opt.converged ? "Optimised" : "Last step")
                    Text("Change")
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(Array(opt.finalParameters.enumerated()), id: \.offset) { index, p in
                    let start = index < opt.initialParameters.count ? opt.initialParameters[index].value : p.value
                    GridRow {
                        Text(p.name)
                        Text(Fmt.parameterValue(start, angle: p.angle))
                        Text(Fmt.parameterValue(p.value, angle: p.angle))
                        Text(Fmt.parameterChange(p.value - start, angle: p.angle))
                    }
                }
                Divider().gridCellUnsizedAxes(.horizontal)
                GridRow {
                    Text("Energy (E_h)")
                    Text(Fmt.energy10(opt.initialEnergy))
                    Text(Fmt.energy10(opt.finalEnergy))
                    Text(String(format: "%+.10f", opt.finalEnergy - opt.initialEnergy))
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)

            VStack(alignment: .leading, spacing: 4) {
                if out.input.molecule == "H2O" {
                    Text(opt.symmetricEnd
                         ? String(format: "C₂ᵥ symmetry kept (|R₁ − R₂| = %.1e bohr); all three internal coordinates were optimised without constraints.", opt.symmetryDeviation)
                         : String(format: "The optimised geometry is not C₂ᵥ (|R₁ − R₂| = %.1e bohr).", opt.symmetryDeviation))
                }
                Text("Gradient: " + (opt.gradient ?? "–"))
                Text("Initial Hessian: \(opt.initialHessian); " + String(format: "%.1f s", opt.seconds))
                Text(String(format: "Criteria (%@): max force %.1e, RMS force %.1e E_h/bohr; max step %.1e, RMS step %.1e bohr; energy change %.0e E_h",
                            out.input.optConvergence, opt.criteria.maxForce, opt.criteria.rmsForce,
                            opt.criteria.maxDisplacement, opt.criteria.rmsDisplacement, opt.criteria.energy))
            }
            .font(.caption)
            .foregroundStyle(.secondary)

            if model.optimizedGeometry != nil {
                HStack {
                    Button {
                        model.useOptimizedGeometry()
                    } label: {
                        Label("Use optimised geometry as input", systemImage: "arrow.uturn.backward")
                    }
                    .secondaryActionStyle()
                    .disabled(!model.canUseOptimizedGeometry || model.isRunning)
                    if !model.canUseOptimizedGeometry {
                        Text("The geometry inputs already hold the optimised geometry.")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                    }
                }
            }
        }
        .dataPanel()
    }

    private func stepsTable(_ opt: MacQCOptimization) -> some View {
        VStack(alignment: .leading, spacing: 6) {
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("Step").gridColumnAlignment(.leading)
                    Text("Energy (E_h)")
                    Text("ΔE (E_h)")
                    Text("Max force")
                    Text("RMS force")
                    Text("Max step")
                    Text("RMS step")
                    Text("")
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(opt.history) { s in
                    GridRow {
                        Button {
                            selectedStep = s.step
                        } label: {
                            Text("\(s.step)")
                                .frame(minWidth: 22, alignment: .leading)
                        }
                        .buttonStyle(.link)
                        .accessibilityLabel("Optimisation step \(s.step)")
                        Text(Fmt.energy10(s.energy))
                        Text(s.deltaE.map { String(format: "%+.2e", $0) } ?? "–")
                        criterionText(s, 0)
                        criterionText(s, 1)
                        criterionText(s, 2)
                        criterionText(s, 3)
                        Image(systemName: s.converged ? "checkmark.circle.fill" : "circle")
                            .foregroundStyle(s.converged ? Color.green : Color.secondary)
                    }
                    .fontWeight(shownStep(opt)?.step == s.step ? .semibold : .regular)
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text("Green values meet their threshold. Click a step number to show its geometry, internal coordinates and convergence table below.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func criterionText(_ s: MacQCOptimizationStep, _ index: Int) -> some View {
        let item = index < s.convergence.count ? s.convergence[index] : nil
        return Text(item?.value.map { String(format: "%.2e", $0) } ?? "–")
            .foregroundStyle(item?.met == true ? Color.green : Color.primary)
    }

    private func charts(_ opt: MacQCOptimization) -> some View {
        let final = opt.finalEnergy
        let energies = opt.history.map { (step: $0.step, value: $0.energy - final) }
        let forces = opt.history.compactMap { s in s.criterion(0).map { (step: s.step, value: max($0, 1.0e-12)) } }
        return HStack(alignment: .top, spacing: 16) {
            VStack(alignment: .leading, spacing: 6) {
                Text("Energy relative to the last step (mE_h)")
                    .font(.subheadline.weight(.semibold))
                Chart {
                    ForEach(energies, id: \.step) { p in
                        LineMark(x: .value("Step", p.step), y: .value("E - E_last (mE_h)", p.value * 1000))
                        PointMark(x: .value("Step", p.step), y: .value("E - E_last (mE_h)", p.value * 1000))
                    }
                }
                .chartXAxisLabel("Step")
                .frame(height: 180)
            }
            .dataPanel()
            VStack(alignment: .leading, spacing: 6) {
                Text("Maximum force (E_h/bohr, logarithmic)")
                    .font(.subheadline.weight(.semibold))
                Chart {
                    ForEach(forces, id: \.step) { p in
                        LineMark(x: .value("Step", p.step), y: .value("Max force", p.value))
                            .foregroundStyle(Color.orange)
                        PointMark(x: .value("Step", p.step), y: .value("Max force", p.value))
                            .foregroundStyle(Color.orange)
                    }
                    RuleMark(y: .value("Threshold", opt.criteria.maxForce))
                        .foregroundStyle(Color.green)
                        .lineStyle(StrokeStyle(lineWidth: 1, dash: [4, 3]))
                        .annotation(position: .top, alignment: .trailing) {
                            Text("threshold").font(.caption2).foregroundStyle(.green)
                        }
                }
                .chartYScale(domain: .automatic(includesZero: false), type: .log)
                .chartXAxisLabel("Step")
                .frame(height: 180)
            }
            .dataPanel()
        }
    }

    private func stepDetail(_ s: MacQCOptimizationStep, opt: MacQCOptimization) -> some View {
        VStack(alignment: .leading, spacing: 14) {
            HStack(alignment: .top, spacing: 16) {
                VStack(alignment: .leading, spacing: 6) {
                    Text("Geometry of this step (Å)")
                        .font(.subheadline.weight(.semibold))
                    AtomTable(atoms: s.geometry.atoms)
                    Text(geometryLine(s.geometry))
                        .font(.caption.monospacedDigit())
                        .foregroundStyle(.secondary)
                }
                .dataPanel()
                VStack(alignment: .leading, spacing: 6) {
                    Text("Convergence")
                        .font(.subheadline.weight(.semibold))
                    Grid(alignment: .trailing, horizontalSpacing: 16, verticalSpacing: 4) {
                        GridRow {
                            Text("Item").gridColumnAlignment(.leading)
                            Text("Value")
                            Text("Threshold")
                            Text("Met")
                        }
                        .foregroundStyle(.secondary)
                        ForEach(s.convergence) { c in
                            GridRow {
                                Text(c.name)
                                Text(c.value.map { String(format: "%.2e", $0) } ?? "–")
                                Text(String(format: "%.2e", c.threshold))
                                Text(c.applicable ? (c.met ? "yes" : "no") : "–")
                                    .foregroundStyle(c.met ? Color.green : Color.primary)
                            }
                        }
                    }
                    .monospacedDigit()
                }
                .dataPanel()
            }
            VStack(alignment: .leading, spacing: 6) {
                Text("Internal coordinates (force = −dE/dq; step = change taken from this geometry)")
                    .font(.subheadline.weight(.semibold))
                Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                    GridRow {
                        Text("Name").gridColumnAlignment(.leading)
                        Text("Definition").gridColumnAlignment(.leading)
                        Text("Value")
                        Text("Force")
                        Text("Step")
                    }
                    .foregroundStyle(.secondary)
                    ForEach(s.coordinates) { q in
                        let angle = !q.definition.hasPrefix("R(")
                        GridRow {
                            Text(q.name)
                            Text(q.definition)
                            Text(angle ? String(format: "%.4f°", q.value * 180 / Double.pi)
                                       : String(format: "%.6f Å", q.value * Units.bohrToAngstrom))
                            Text(String(format: "%+.6f", -q.gradient))
                            Text(q.step.map { angle ? String(format: "%+.4f°", $0 * 180 / Double.pi)
                                                    : String(format: "%+.6f Å", $0 * Units.bohrToAngstrom) } ?? "–")
                        }
                    }
                }
                .monospacedDigit()
                .textSelection(.enabled)
                Text(stepCaption(s))
                    .font(.caption)
                    .foregroundStyle(.secondary)
                ForEach(s.warnings, id: \.self) { w in
                    Label(w, systemImage: "exclamationmark.triangle.fill")
                        .font(.caption)
                        .foregroundStyle(.orange)
                }
            }
            .dataPanel()
        }
    }

    private func geometryLine(_ g: MacQCGeometry) -> String {
        if let r1 = g.R1, let r2 = g.R2, let a = g.angle {
            return String(format: "R₁ = %.6f Å, R₂ = %.6f Å, angle = %.4f°", r1 * Units.bohrToAngstrom,
                          r2 * Units.bohrToAngstrom, a)
        }
        if let r = g.R { return String(format: "R = %.6f Å (%.6f bohr)", r * Units.bohrToAngstrom, r) }
        return ""
    }

    private func stepCaption(_ s: MacQCOptimizationStep) -> String {
        var parts = ["Hessian: \(s.hessianUpdate)"]
        if s.stepNorm > 0 {
            parts.append(String(format: "RFO step length %.4f (trust radius %.3f%@), predicted change %.2e E_h",
                                s.stepNorm, s.trustRadius, s.stepScaled ? ", scaled to it" : "", s.predictedChange))
        }
        if let r = s.ratio {
            parts.append(String(format: "previous step: actual/predicted energy change %.3f", r))
        }
        parts.append(String(format: "%.1f s", s.seconds))
        return parts.joined(separator: "; ")
    }
}

/// Cartesian coordinates of a geometry in Å.
struct AtomTable: View {
    var atoms: [MacQCAtom]

    var body: some View {
        Grid(alignment: .trailing, horizontalSpacing: 16, verticalSpacing: 3) {
            GridRow {
                Text("Atom").gridColumnAlignment(.leading)
                Text("x")
                Text("y")
                Text("z")
            }
            .foregroundStyle(.secondary)
            ForEach(Array(atoms.enumerated()), id: \.offset) { index, a in
                GridRow {
                    Text("\(a.symbol)\(index + 1)")
                    ForEach(0..<3, id: \.self) { k in
                        Text(String(format: "%.6f", (k < a.position.count ? a.position[k] : 0) * Units.bohrToAngstrom))
                    }
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
    }
}

/// Placeholder of a job tab without a result: refusal, warning or how to get one.
struct JobPlaceholder: View {
    var model: AppModel
    var title: String
    var systemImage: String
    var text: String

    var body: some View {
        VStack(spacing: 14) {
            if let message = model.errorMessage {
                ErrorBanner(message: message)
                    .frame(maxWidth: 560)
            }
            if model.isGeneral, let refusal = model.refusal {
                StatusBanner(kind: .refusal, text: refusal)
                    .frame(maxWidth: 560)
            }
            ContentUnavailableView(title, systemImage: systemImage,
                                   description: Text(model.isGeneral ? text : "Jobs other than the energy are available for H₂, LiH and H₂O."))
        }
        .padding()
    }
}

// MARK: - Frequencies and thermochemistry

struct FrequenciesView: View {
    @Bindable var model: AppModel

    var body: some View {
        VStack(spacing: 0) {
            Picker("Page", selection: $model.vibrationPage) {
                ForEach(VibrationPage.allCases) { page in
                    Text(page.rawValue).tag(page)
                }
            }
            .pickerStyle(.segmented)
            .labelsHidden()
            .frame(maxWidth: 360)
            .padding(.top, 10)
            .padding(.bottom, 4)

            if let out = model.result, let f = out.frequencies {
                switch model.vibrationPage {
                case .frequencies:
                    FrequencyPage(model: model, out: out, f: f)
                case .thermochemistry:
                    if let t = out.thermochemistry {
                        ThermochemistryPage(model: model, out: out, t: t)
                    } else {
                        Spacer()
                    }
                }
            } else {
                JobPlaceholder(model: model, title: model.vibrationPage == .frequencies ? "No frequencies" : "No thermochemistry",
                               systemImage: "waveform.path",
                               text: "Choose the job Frequencies, Optimisation + frequencies or Optimisation + frequencies + single-point energy in the sidebar and press Run.")
                    .frame(maxHeight: .infinity)
            }
        }
    }
}

struct FrequencyPage: View {
    var model: AppModel
    var out: MacQCOutput
    var f: MacQCFrequencies
    @State private var selectedMode = 0
    @State private var animate = true

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 18) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                }
                VStack(alignment: .leading, spacing: 4) {
                    Text("\(Fmt.formula(out.input.molecule)) · harmonic frequencies · \(MethodChoice.title(of: out.input.method))/\(out.input.basis)")
                        .font(.title2.weight(.semibold))
                    Text(out.input.job == "freq" ? "At the input geometry." : "At the optimised geometry.")
                        .foregroundStyle(.secondary)
                }
                NotesPanel(notes: [], warnings: (out.job?.warnings ?? []) + (out.thermochemistry?.warnings ?? []))

                SectionTitle("How the Hessian was computed")
                VStack(alignment: .leading, spacing: 6) {
                    InfoRow("Hessian", f.analytic ? "analytic" : "numerical")
                    Text(f.description)
                        .font(.callout)
                        .textSelection(.enabled)
                    if let reason = f.numericalReason {
                        Text("Numerical because: " + reason)
                            .font(.caption)
                            .foregroundStyle(.secondary)
                    }
                    if f.gradientEvaluations > 0 {
                        InfoRow("Gradient evaluations", "\(f.gradientEvaluations)")
                    }
                    if f.energyEvaluations > 0 {
                        InfoRow("Energy evaluations", "\(f.energyEvaluations)")
                    }
                    if !f.coordinates.isEmpty {
                        InfoRow("Displaced coordinates", f.coordinates.joined(separator: ", "))
                    }
                    InfoRow("IR intensities", f.intensityDescription.isEmpty ? (f.intensities ? "yes" : "none") : f.intensityDescription)
                    InfoRow("Isotopes", zip(f.isotopes, f.masses).map { String(format: "%@ %.9f u", $0.0, $0.1) }.joined(separator: ", "))
                    InfoRow("Residual curvature (cm⁻¹)", f.lowFrequencies.map { String(format: "%.2f", $0) }.joined(separator: " "))
                }
                .dataPanel()

                SectionTitle("Normal modes")
                modesTable
                if !f.modes.isEmpty {
                    HStack(alignment: .top, spacing: 16) {
                        modeDrawing
                        irSpectrum
                    }
                }
            }
            .padding(20)
        }
    }

    private var modesTable: some View {
        Grid(alignment: .trailing, horizontalSpacing: 22, verticalSpacing: 4) {
            GridRow {
                Text("Mode").gridColumnAlignment(.leading)
                Text("Symmetry").gridColumnAlignment(.leading)
                Text("Wavenumber (cm⁻¹)")
                Text("Reduced mass (u)")
                Text("Force constant (mdyn/Å)")
                Text("IR intensity (km/mol)")
            }
            .font(.subheadline.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider().gridCellUnsizedAxes(.horizontal)
            ForEach(Array(f.modes.enumerated()), id: \.offset) { index, m in
                GridRow {
                    Button {
                        selectedMode = index
                    } label: {
                        Text("\(index + 1)").frame(minWidth: 22, alignment: .leading)
                    }
                    .buttonStyle(.link)
                    .accessibilityLabel("Normal mode \(index + 1)")
                    Text(m.symmetry.isEmpty ? "–" : m.symmetry)
                    Text(Fmt.wavenumber(m))
                        .foregroundStyle(m.imaginary ? Color.red : Color.primary)
                    Text(String(format: "%.6f", m.reducedMass))
                    Text(String(format: "%.6f", m.forceConstant))
                    Text(m.irIntensity.map { String(format: "%.4f", $0) } ?? "–")
                }
                .fontWeight(index == selectedMode ? .semibold : .regular)
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }

    private var modeDrawing: some View {
        let index = min(selectedMode, f.modes.count - 1)
        let mode = f.modes[index]
        return VStack(alignment: .leading, spacing: 8) {
            HStack {
                Text("Mode \(index + 1): \(Fmt.wavenumber(mode)) cm⁻¹ \(mode.symmetry)")
                    .font(.subheadline.weight(.semibold))
                Spacer()
                Toggle("Animate", isOn: $animate)
                    .toggleStyle(.switch)
                    .controlSize(.small)
            }
            NormalModeCanvas(atoms: f.geometry.atoms, displacement: mode.displacement, animate: animate)
                .frame(height: 220)
            Text("Arrows: normalised Cartesian displacements, drawn in the plane of the molecule. Click a mode number in the table to show it.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private var irSpectrum: some View {
        let sticks = f.modes.enumerated().compactMap { item -> (index: Int, w: Double, i: Double)? in
            guard let i = item.element.irIntensity, !item.element.imaginary else { return nil }
            return (item.offset, item.element.wavenumber, i)
        }
        return VStack(alignment: .leading, spacing: 8) {
            Text("IR stick spectrum")
                .font(.subheadline.weight(.semibold))
            if sticks.isEmpty {
                Text(f.intensities ? "No IR-active real modes." : "No IR intensities for this calculation.")
                    .foregroundStyle(.secondary)
                    .frame(maxWidth: .infinity, minHeight: 220)
            } else {
                Chart {
                    ForEach(sticks, id: \.index) { s in
                        RuleMark(x: .value("Wavenumber (cm⁻¹)", s.w), yStart: .value("Intensity", 0), yEnd: .value("Intensity", s.i))
                            .lineStyle(StrokeStyle(lineWidth: 3))
                            .foregroundStyle(s.index == selectedMode ? Color.orange : Color.blue)
                    }
                }
                .chartXScale(domain: 0...(max(4500, (sticks.map { $0.w }.max() ?? 4000) + 300)))
                .chartXAxisLabel("Wavenumber (cm⁻¹)")
                .chartYAxisLabel("km/mol")
                .frame(height: 220)
            }
        }
        .dataPanel()
    }
}

/// Atoms (circles) and the displacement of a normal mode (arrows), projected on the two
/// Cartesian axes along which the molecule extends most; optionally animated.
struct NormalModeCanvas: View {
    var atoms: [MacQCAtom]
    var displacement: [[Double]]
    var animate: Bool

    var body: some View {
        TimelineView(.animation(minimumInterval: 1.0 / 30.0, paused: !animate)) { context in
            let phase = animate ? sin(context.date.timeIntervalSinceReferenceDate * 2 * Double.pi / 1.6) : 0
            Canvas { gc, size in
                draw(gc, size: size, phase: phase)
            }
        }
        .accessibilityLabel("Normal mode drawing")
    }

    private func axes() -> (Int, Int) {
        var spread = [0.0, 0.0, 0.0]
        for k in 0..<3 {
            let values = atoms.map { k < $0.position.count ? $0.position[k] : 0 }
            spread[k] = (values.max() ?? 0) - (values.min() ?? 0)
        }
        let order = [0, 1, 2].sorted { spread[$0] > spread[$1] }
        // horizontal: largest extent; vertical: second (for a linear molecule a zero range)
        return (order[0], order[1])
    }

    private func draw(_ gc: GraphicsContext, size: CGSize, phase: Double) {
        guard !atoms.isEmpty else { return }
        let (h, v) = axes()
        let width = Double(size.width), height = Double(size.height)
        let xs: [Double] = atoms.map { h < $0.position.count ? $0.position[h] : 0 }
        let ys: [Double] = atoms.map { v < $0.position.count ? $0.position[v] : 0 }
        let cx = ((xs.max() ?? 0) + (xs.min() ?? 0)) / 2
        let cy = ((ys.max() ?? 0) + (ys.min() ?? 0)) / 2
        let extent = max((xs.max() ?? 0) - (xs.min() ?? 0), (ys.max() ?? 0) - (ys.min() ?? 0), 1.0)
        let scale: Double = 0.55 * min(width, height * 1.6) / extent
        let arrowScale = 1.6   // bohr per unit of the normalised displacement
        func point(_ x: Double, _ y: Double) -> CGPoint {
            CGPoint(x: width / 2 + (x - cx) * scale, y: height / 2 - (y - cy) * scale)
        }
        func shift(_ i: Int) -> [Double] {
            displacement.count > i && displacement[i].count == 3 ? displacement[i] : [0, 0, 0]
        }
        // bonds between neighbouring atoms (distance below 3.6 bohr)
        for i in 0..<atoms.count {
            for j in (i + 1)..<atoms.count {
                var d2 = 0.0
                for k in 0..<3 {
                    let a = k < atoms[i].position.count ? atoms[i].position[k] : 0
                    let b = k < atoms[j].position.count ? atoms[j].position[k] : 0
                    d2 += (a - b) * (a - b)
                }
                guard d2 < 3.6 * 3.6 else { continue }
                let di = shift(i), dj = shift(j)
                var path = Path()
                path.move(to: point(xs[i] + 0.25 * phase * di[h], ys[i] + 0.25 * phase * di[v]))
                path.addLine(to: point(xs[j] + 0.25 * phase * dj[h], ys[j] + 0.25 * phase * dj[v]))
                gc.stroke(path, with: .color(.secondary), lineWidth: 4)
            }
        }
        for (i, atom) in atoms.enumerated() {
            let d = shift(i)
            let p = point(xs[i] + 0.25 * phase * d[h], ys[i] + 0.25 * phase * d[v])
            let radius: Double = atom.Z == 1 ? 11 : (atom.Z == 3 ? 17 : 15)
            let color: Color = atom.Z == 1 ? Color(white: 0.85) : (atom.Z == 8 ? .red : .purple)
            let rect = CGRect(x: Double(p.x) - radius, y: Double(p.y) - radius, width: 2 * radius, height: 2 * radius)
            gc.fill(Path(ellipseIn: rect), with: .color(color))
            gc.stroke(Path(ellipseIn: rect), with: .color(Color.primary.opacity(0.6)), lineWidth: 1)
            gc.draw(Text(atom.symbol).font(.caption.bold()).foregroundStyle(atom.Z == 1 ? Color.black : Color.white), at: p)
            // arrow of the displacement
            let tip = point(xs[i] + arrowScale * d[h], ys[i] + arrowScale * d[v])
            let base = point(xs[i], ys[i])
            let dx = Double(tip.x - base.x), dy = Double(tip.y - base.y)
            let length = (dx * dx + dy * dy).squareRoot()
            guard length > 4 else { continue }
            let ux = dx / length, uy = dy / length
            let tx = Double(tip.x), ty = Double(tip.y)
            var arrow = Path()
            arrow.move(to: base)
            arrow.addLine(to: tip)
            arrow.move(to: tip)
            arrow.addLine(to: CGPoint(x: tx - 9 * ux + 5 * uy, y: ty - 9 * uy - 5 * ux))
            arrow.move(to: tip)
            arrow.addLine(to: CGPoint(x: tx - 9 * ux - 5 * uy, y: ty - 9 * uy + 5 * ux))
            gc.stroke(arrow, with: .color(.orange), lineWidth: 2)
        }
    }
}

struct ThermochemistryPage: View {
    var model: AppModel
    var out: MacQCOutput
    var t: MacQCThermochemistry

    private func hartree(_ jPerMol: Double) -> Double { jPerMol / Units.joulePerMolPerHartree }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 18) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                }
                VStack(alignment: .leading, spacing: 4) {
                    Text("\(Fmt.formula(out.input.molecule)) · thermochemistry")
                        .font(.title2.weight(.semibold))
                    Text(String(format: "Ideal gas, rigid rotor, harmonic oscillators · T = %.2f K · p = %.2f Pa (%.6f atm, %.5f bar) · symmetry number %ld",
                                t.temperature, t.pressure, t.pressure / 101325, t.pressure / 1.0e5, t.symmetryNumber))
                        .foregroundStyle(.secondary)
                }
                NotesPanel(notes: [], warnings: t.warnings)

                SectionTitle("Energies")
                energyTable
                if let sp = out.singlePoint, let label = out.job?.finalEnergyLabel {
                    Text("The electronic energy is the single-point energy \(label) (\(Fmt.energy10(sp.energy)) E_h); the thermal corrections are from the frequency step (\(MethodChoice.title(of: out.input.method))/\(out.input.basis)).")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }

                SectionTitle("Entropy, heat capacity and thermal energy by contribution")
                contributionTable

                SectionTitle("Rotation and vibration")
                VStack(spacing: 8) {
                    InfoRow("Molecular mass", String(format: "%.6f u", t.mass))
                    InfoRow("Rotor", t.linear ? "linear" : "nonlinear (asymmetric top)")
                    InfoRow("Moments of inertia (u bohr²)", t.moments.map { String(format: "%.6f", $0) }.joined(separator: "  "))
                    InfoRow("Rotational constants (GHz)", t.rotationalConstantsGHz.map { String(format: "%.6f", $0) }.joined(separator: "  "))
                    InfoRow("Rotational constants (cm⁻¹)", t.rotationalConstantsCm.map { String(format: "%.6f", $0) }.joined(separator: "  "))
                    InfoRow("Rotational temperatures (K)", t.rotationalTemperatures.map { String(format: "%.5f", $0) }.joined(separator: "  "))
                    InfoRow("Vibrational temperatures (K)", t.vibrationalTemperatures.map { String(format: "%.2f", $0) }.joined(separator: "  "))
                    if t.excludedModes > 0 {
                        InfoRow("Modes left out", "\(t.excludedModes) (imaginary or below 10 cm⁻¹)")
                    }
                }
                .dataPanel()
                Text("Constants: CODATA 2018; masses of the most abundant isotopes (AME 2020); electronic ground state a closed-shell singlet (q_el = 1).")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .padding(20)
        }
    }

    private var energyTable: some View {
        let rows: [(label: String, value: Double)] = [
            ("Zero-point energy", t.corrections.zeroPoint),
            ("Thermal correction to U", t.corrections.thermalEnergy),
            ("Thermal correction to H = U + kT", t.corrections.enthalpy),
            ("Thermal correction to A = U − TS", t.corrections.helmholtz),
            ("Thermal correction to G = H − TS", t.corrections.gibbs),
        ]
        let totals: [(label: String, value: Double)] = [
            ("E (electronic)", t.electronicEnergy),
            ("E + ZPE", t.totals.zeroPoint),
            ("U", t.totals.thermalEnergy),
            ("H", t.totals.enthalpy),
            ("A", t.totals.helmholtz),
            ("G", t.totals.gibbs),
        ]
        return HStack(alignment: .top, spacing: 16) {
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("Correction").gridColumnAlignment(.leading)
                    Text("E_h")
                    Text("kcal/mol")
                    Text("kJ/mol")
                }
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(rows, id: \.label) { row in
                    GridRow {
                        Text(row.label)
                        Text(String(format: "%.8f", hartree(row.value)))
                        Text(String(format: "%.4f", row.value / 1000 / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.value / 1000))
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            .dataPanel()
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("Total").gridColumnAlignment(.leading)
                    Text("E_h")
                }
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(totals, id: \.label) { row in
                    GridRow {
                        Text(row.label)
                        Text(Fmt.energy10(row.value))
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            .dataPanel()
        }
    }

    private var contributionTable: some View {
        let terms: [(label: String, term: MacQCThermoTerm)] = [
            ("Electronic", t.terms.electronic),
            ("Translational", t.terms.translation),
            ("Rotational", t.terms.rotation),
            ("Vibrational", t.terms.vibration),
            ("Total", t.terms.total),
        ]
        return VStack(alignment: .leading, spacing: 6) {
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("").gridColumnAlignment(.leading)
                    Text("E (kcal/mol)")
                    Text("C_V (cal/(mol K))")
                    Text("S (cal/(mol K))")
                    Text("S (J/(mol K))")
                    Text("q")
                }
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(terms, id: \.label) { row in
                    GridRow {
                        Text(row.label)
                        Text(String(format: "%.4f", row.term.energy / 1000 / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.term.cv / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.term.entropy / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.term.entropy))
                        Text(String(format: "%.6e", row.term.q))
                    }
                    .fontWeight(row.label == "Total" ? .semibold : .regular)
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text("The vibrational energy includes the zero-point energy; the vibrational partition function is counted from the bottom of the well.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }
}

// MARK: - Summary panel of a job

/// Job, final energy, optimised geometry, wavenumbers and thermochemistry (Summary tab).
struct JobSummaryPanel: View {
    var model: AppModel
    var out: MacQCOutput

    var body: some View {
        if let j = out.job, j.type != "energy" {
            VStack(alignment: .leading, spacing: 8) {
                InfoRow("Job", j.title)
                InfoRow("Termination", j.termination == "normal" ? "normal" : (j.message ?? j.termination))
                if let label = j.finalEnergyLabel, let e = j.finalEnergy {
                    InfoRow("Final energy (\(label))", Fmt.energy10(e) + " E_h")
                }
                if let opt = out.optimization {
                    InfoRow("Optimisation", opt.converged ? "converged in \(opt.steps) steps" : "not converged (\(opt.steps) steps)")
                    ForEach(opt.finalParameters) { p in
                        InfoRow("  " + p.name, Fmt.parameterValue(p.value, angle: p.angle))
                    }
                }
                if let f = out.frequencies {
                    InfoRow("Harmonic wavenumbers (cm⁻¹)", f.modes.map { Fmt.wavenumber($0) }.joined(separator: "  "))
                    InfoRow("Hessian", f.analytic ? "analytic" : "numerical (\(f.gradientEvaluations > 0 ? "\(f.gradientEvaluations) gradients" : "\(f.energyEvaluations) energies"))")
                }
                if let t = out.thermochemistry {
                    InfoRow(String(format: "ZPE, H, G at %.2f K", t.temperature),
                            String(format: "%.8f, %@, %@ E_h", t.corrections.zeroPoint / Units.joulePerMolPerHartree,
                                   Fmt.energy10(t.totals.enthalpy), Fmt.energy10(t.totals.gibbs)))
                    InfoRow("Entropy", String(format: "%.4f J/(mol K)", t.terms.total.entropy))
                }
                if model.optimizedGeometry != nil {
                    Button {
                        model.useOptimizedGeometry()
                    } label: {
                        Label("Use optimised geometry as input", systemImage: "arrow.uturn.backward")
                    }
                    .secondaryActionStyle()
                    .disabled(!model.canUseOptimizedGeometry || model.isRunning)
                }
            }
            .dataPanel()
        }
    }
}
