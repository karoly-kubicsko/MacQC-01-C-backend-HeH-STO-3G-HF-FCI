//
//  AppModel.swift
//  MacQC
//
//  Observable application state: sidebar inputs (molecule, job type, method, basis set,
//  thermochemistry and optimisation settings), size estimate, the running job with its streamed
//  output file, last result, potential-curve scan, and the menu actions (export text / JSON,
//  copy summary).
//

import SwiftUI
import AppKit
import UniformTypeIdentifiers

/// The six tabs of the detail area (ΨLab 5.0). Methods, Frequencies and Output have pages.
enum DetailTab: Hashable {
    case summary, methods, optimization, frequencies, potentialCurve, output
}

/// Pages of the Methods tab.
enum MethodsPage: String, CaseIterable, Identifiable {
    case hartreeFock = "Hartree–Fock"
    case coupledCluster = "Coupled cluster"
    case fullCI = "Full CI"
    var id: String { rawValue }
}

/// Pages of the Frequencies tab.
enum VibrationPage: String, CaseIterable, Identifiable {
    case frequencies = "Frequencies"
    case thermochemistry = "Thermochemistry"
    var id: String { rawValue }
}

/// One point of the potential-curve scan: the (R)HF energy and that of the correlated method
/// (Full CI or coupled cluster; missing in an HF-only scan).
struct CurvePoint: Identifiable, Sendable {
    var R: Double
    var hf: Double
    var correlated: Double?
    var id: Double { R }
}

/// Flattened (method, R, energy) sample for Swift Charts.
struct CurveSample: Identifiable {
    let id = UUID()
    var method: String
    var R: Double
    var energy: Double
}

/// Minimum of one curve (refined by a parabola through the three lowest grid points).
struct CurveMinimum: Identifiable {
    var method: String
    var R: Double
    var energy: Double
    var refined: Bool
    var id: String { method }
}

@MainActor
@Observable
final class AppModel {

    // MARK: Sidebar inputs

    static let defaultR = 1.4632
    static let defaultZeta1 = 2.0925
    static let defaultZeta2 = 1.24
    static let defaultGaussians = 3

    var molecule: MoleculeChoice = .heh

    // HeH+ (Szabo–Ostlund)
    var R: Double = AppModel.defaultR             // bohr (the backend's unit)
    var zeta1: Double = AppModel.defaultZeta1   // He
    var zeta2: Double = AppModel.defaultZeta2   // H
    var nGaussians: Int = AppModel.defaultGaussians
    var exactErf: Bool = false

    // H2, LiH, H2O
    var basis: String = "aug-cc-pVDZ"
    var method: MethodChoice = .hf
    var bondLengths: [MoleculeChoice: Double] = [:]   // bohr; missing = experimental value
    var waterAngle: Double = MoleculeChoice.defaultWaterAngle

    // Job (ΨLab 5.0)
    var job: JobChoice = .energy
    var spMethod: MethodChoice? = nil      // single point of job 5; nil = the job's method
    var spBasis: String? = nil             // nil = the job's basis set
    var temperature: Double = 298.15       // K
    var pressure: Double = 101325          // Pa
    var tightConvergence = false
    var maxSteps = 50

    /// The bond length in angstrom; reads and writes `R` (bohr).
    var rAngstrom: Double {
        get { R * Units.bohrToAngstrom }
        set { R = newValue / Units.bohrToAngstrom }
    }

    /// Bond length (bohr) of the selected general molecule.
    var bondLength: Double {
        get { bondLengths[molecule] ?? molecule.defaultBondLength }
        set { bondLengths[molecule] = newValue }
    }

    var bondLengthAngstrom: Double {
        get { bondLength * Units.bohrToAngstrom }
        set { bondLength = newValue / Units.bohrToAngstrom }
    }

    var isGeneral: Bool { molecule != .heh }

    // MARK: Size estimate (general molecules)

    /// Output of an `estimateOnly` call for the current inputs (nil for HeH+).
    var estimate: MacQCOutput?

    /// Changes whenever an input that affects the estimate changes.
    var estimateKey: String {
        guard isGeneral else { return "HeH+" }
        return "\(molecule.rawValue)|\(basis)|\(method.rawValue)|\(bondLength)|\(waterAngle)|\(job.rawValue)|"
            + "\(spMethod?.rawValue ?? "")|\(spBasis ?? "")|\(temperature)|\(pressure)|\(tightConvergence)|\(maxSteps)"
    }

    /// The backend refuses the calculation (reason), or the inputs are invalid.
    var refusal: String? {
        guard isGeneral, let e = estimate else { return nil }
        if let r = e.job?.estimate.refusal { return r }
        if let r = e.refusal { return r }
        if e.estimate == nil, let message = e.error { return message }
        return nil
    }

    /// Refusal of a single energy with these settings (potential curves run energies only).
    var energyRefusal: String? {
        guard isGeneral, let e = estimate else { return nil }
        if let r = e.refusal { return r }
        if e.estimate == nil, let message = e.error { return message }
        return nil
    }

    /// The whole job is estimated to take longer than 30 s.
    var isLongCalculation: Bool {
        estimate?.job?.estimate.longCalculation ?? estimate?.estimate?.longCalculation ?? false
    }

    /// Estimated time of the whole job (seconds).
    var estimatedSeconds: Double? {
        estimate?.job?.estimate.totalSeconds ?? estimate?.estimate?.totalSeconds
    }

    /// The displayed result was computed with the current molecule, job, basis set, method and
    /// starting geometry.
    var resultMatchesInputs: Bool {
        guard let input = result?.input, isGeneral, input.isGeneral else { return false }
        return input.molecule == molecule.rawValue && input.basis == basis && input.method == method.rawValue
            && input.job == job.rawValue
            && abs(input.R - bondLength) < 1e-9 && (molecule != .h2o || abs(input.angle - waterAngle) < 1e-9)
    }

    /// The single-point method and basis set of job 5 (the job's own when not chosen).
    var effectiveSPMethod: MethodChoice { spMethod ?? method }
    var effectiveSPBasis: String { spBasis ?? basis }

    // MARK: Results and the running job

    var selectedTab: DetailTab = .summary
    var methodsPage: MethodsPage = .hartreeFock
    var vibrationPage: VibrationPage = .frequencies
    var isRunning = false
    var result: MacQCOutput?
    var errorMessage: String?
    var progress = MacQCProgress()
    /// The backend's description of the current step ("Optimisation step 4 · E = … · max force …").
    var progressLine = ""
    var runStarted: Date?
    var showLongRunConfirmation = false

    /// The output file of the latest run (streamed while it runs; also after a cancelled or
    /// failed run). `outputGeneration` changes when the text is replaced rather than extended.
    var outputText = ""
    var outputGeneration = 0

    /// Starts the job at launch even when it is estimated to be long (launch argument
    /// -psilabAutoRun YES, used by the UI smoke test).
    var autoRunLongJobs = false

    @ObservationIgnored private var runningJob: MacQCJob?

    // MARK: Potential curve

    /// Step choices in angstrom (converted to bohr for the scan grid).
    static let scanStepChoices: [Double] = [0.025, 0.05, 0.075, 0.1, 0.15]

    var scanStep: Double = 0.05   // angstrom
    var isScanning = false
    var scanProgress: Double = 0
    var scanMessage: String?
    var curvePoints: [CurvePoint] = []
    var curveInputDescription: String?   // parameters the scan was computed with
    var curveHFLabel = "HF"               // "HF" for HeH+, "RHF" for the other molecules
    var curveCorrelatedLabel = "FCI"      // "FCI", "CCSD", "CCSD(T)", "CCSDT" or "CCSDTQ"
    var showLongScanConfirmation = false

    @ObservationIgnored private var scanTask: Task<Void, Never>?
    @ObservationIgnored private var scanJob: MacQCJob?

    // MARK: Start-up

    init() {
        // Launch arguments (also used by the UI test): -psilabMolecule H2O -psilabBasis aug-cc-pVTZ
        // -psilabMethod FCI (or HF, CCSD, CCSD(T), CCSDT, CCSDTQ)
        let defaults = UserDefaults.standard
        if let m = defaults.string(forKey: "psilabMolecule"),
           let choice = MoleculeChoice.allCases.first(where: { $0.rawValue.lowercased() == m.lowercased() }) {
            molecule = choice
        }
        if let b = defaults.string(forKey: "psilabBasis"),
           let name = BasisSets.all.first(where: { $0.lowercased() == b.lowercased() }) {
            basis = name
        }
        if let m = defaults.string(forKey: "psilabMethod"),
           let choice = MethodChoice(name: m) {
            method = choice
        }
        // ΨLab 5.0: -psilabJob opt+freq -psilabSPMethod "CCSD(T)" -psilabSPBasis aug-cc-pVTZ
        // -psilabR 1.0 (angstrom) -psilabAngle 110 -psilabTab output -psilabAutoRun YES
        if let j = defaults.string(forKey: "psilabJob"), let choice = JobChoice(name: j) {
            job = choice
        }
        if let m = defaults.string(forKey: "psilabSPMethod"), let choice = MethodChoice(name: m) {
            spMethod = choice
        }
        if let b = defaults.string(forKey: "psilabSPBasis"),
           let name = BasisSets.all.first(where: { $0.lowercased() == b.lowercased() }) {
            spBasis = name
        }
        if defaults.object(forKey: "psilabR") != nil, molecule != .heh {
            let r = defaults.double(forKey: "psilabR")
            if r > 0 { bondLengths[molecule] = r / Units.bohrToAngstrom }
        }
        if defaults.object(forKey: "psilabAngle") != nil {
            let a = defaults.double(forKey: "psilabAngle")
            if a > 0 { waterAngle = a }
        }
        if let t = defaults.string(forKey: "psilabTab") {
            switch t.lowercased() {
            case "methods": selectedTab = .methods
            case "optimisation", "optimization": selectedTab = .optimization
            case "frequencies": selectedTab = .frequencies
            case "thermochemistry":
                selectedTab = .frequencies
                vibrationPage = .thermochemistry
            case "curve", "potential curve": selectedTab = .potentialCurve
            case "output": selectedTab = .output
            default: selectedTab = .summary
            }
        }
        autoRunLongJobs = defaults.bool(forKey: "psilabAutoRun")
    }

    // MARK: Inputs

    func currentInput(printLevel: Int = 2) -> MacQCInput {
        var input = MacQCInput()
        input.printLevel = printLevel
        if isGeneral {
            input.molecule = molecule.rawValue
            input.R = bondLength
            input.angle = waterAngle
            input.basis = basis
            input.method = method.rawValue
            input.basisDirectory = BasisSets.directory
            input.threads = 0
            input.memoryMB = MachineLimits.integralMemoryMB
            input.fciMemoryMB = MachineLimits.fciMemoryMB
            input.job = job.rawValue
            if job.hasSinglePoint {
                input.spMethod = spMethod?.rawValue ?? ""
                input.spBasis = spBasis ?? ""
            }
            input.temperature = temperature
            input.pressure = pressure
            input.optConvergence = tightConvergence ? "tight" : "normal"
            input.maxSteps = maxSteps
        } else {
            input.molecule = MoleculeChoice.heh.rawValue
            input.R = R
            input.zeta1 = zeta1
            input.zeta2 = zeta2
            input.ZA = 2
            input.ZB = 1
            input.nGaussians = nGaussians
            input.exactErf = exactErf
        }
        return input
    }

    func resetToDefaults() {
        R = AppModel.defaultR
        zeta1 = AppModel.defaultZeta1
        zeta2 = AppModel.defaultZeta2
        nGaussians = AppModel.defaultGaussians
        exactErf = false
    }

    /// Experimental equilibrium geometry of the selected molecule.
    func resetGeometry() {
        bondLengths[molecule] = nil
        if molecule == .h2o { waterAngle = MoleculeChoice.defaultWaterAngle }
    }

    /// Bond length (bohr; H₂O the mean of both O–H distances) and angle of the optimised
    /// geometry of the displayed result, if it is for the selected molecule and converged.
    var optimizedGeometry: (R: Double, angle: Double?, symmetric: Bool)? {
        guard isGeneral, let out = result, out.input.molecule == molecule.rawValue,
              let opt = out.optimization, opt.converged, let r = opt.finalGeometry.R else { return nil }
        return (r, opt.finalGeometry.angle, opt.symmetricEnd || out.input.molecule != "H2O")
    }

    /// The optimised geometry differs from the inputs (so copying it changes something).
    var canUseOptimizedGeometry: Bool {
        guard let g = optimizedGeometry else { return false }
        if abs(g.R - bondLength) > 1e-9 { return true }
        if let a = g.angle, abs(a - waterAngle) > 1e-9 { return true }
        return false
    }

    /// Copies the optimised geometry into the geometry inputs (the start of the next job).
    func useOptimizedGeometry() {
        guard let g = optimizedGeometry else { return }
        bondLengths[molecule] = g.R
        if molecule == .h2o, let a = g.angle { waterAngle = a }
    }

    // MARK: Estimate

    /// Asks the backend for the size, time and memory of the current input (no calculation).
    func refreshEstimate() async {
        guard isGeneral else {
            estimate = nil
            return
        }
        let key = estimateKey
        var input = currentInput(printLevel: 0)
        input.estimateOnly = true
        do {
            let output = try await MacQCEngine.run(input)
            if key == estimateKey {   // inputs unchanged meanwhile
                estimate = output
            }
        } catch {
            if key == estimateKey {
                estimate = nil
            }
        }
    }

    /// Text of the estimate for a confirmation dialog.
    var longRunMessage: String {
        if let j = estimate?.job, !j.estimate.warnings.isEmpty {
            return j.estimate.warnings.joined(separator: "\n")
        }
        return estimate?.warnings?.joined(separator: "\n") ?? "This calculation may take a long time."
    }

    // MARK: Single calculation

    /// Run button: asks for confirmation first when the estimate warns of a long calculation.
    func requestRun() {
        guard !isRunning else { return }
        if isGeneral && refusal != nil { return }
        if isGeneral && isLongCalculation {
            showLongRunConfirmation = true
        } else {
            Task { await run() }
        }
    }

    func run() async {
        guard !isRunning else { return }
        isRunning = true
        errorMessage = nil
        progress = MacQCProgress()
        progressLine = ""
        runStarted = Date()
        let input = currentInput(printLevel: 2)
        if input.isGeneral {
            outputText = ""
            outputGeneration += 1
        }
        let handle = MacQCJob()
        runningJob = handle
        // Poll the backend's progress and the new part of the output file while the job runs.
        let poller = Task { [weak self] in
            var offset = 0
            while !Task.isCancelled {
                let p = handle.progress()
                let line = handle.progressText()
                let piece = handle.output(from: offset)
                offset = max(offset, piece.total)
                guard let self else { return }
                self.progress = p
                self.progressLine = line
                if !piece.text.isEmpty { self.outputText += piece.text }
                try? await Task.sleep(for: .milliseconds(250))
            }
        }
        defer {
            poller.cancel()
            isRunning = false
            runningJob = nil
            runStarted = nil
            progressLine = ""
        }

        do {
            let output = try await MacQCEngine.run(input, job: handle)
            // The complete output file (also of a cancelled or failed job).
            if input.isGeneral {
                let full = output.output ?? handle.output(from: 0).text
                if !full.isEmpty && full != outputText {
                    outputText = full
                    outputGeneration += 1
                }
            }
            if output.error == "calculation cancelled" {
                errorMessage = "Calculation cancelled."   // keep the previous result visible
                return
            }
            result = output
            errorMessage = output.error
        } catch {
            // Keep the previous result visible; show why the new run failed.
            errorMessage = error.localizedDescription
        }
    }

    /// The progress line under the Run button.
    var progressDescription: String {
        progressLine.isEmpty ? progress.text : progressLine
    }

    func cancelRun() {
        runningJob?.cancel()
    }

    /// Elapsed time of the running calculation, e.g. "1:05".
    func elapsedText(at now: Date) -> String {
        guard let start = runStarted else { return "" }
        let s = Int(now.timeIntervalSince(start))
        return String(format: "%d:%02d", s / 60, s % 60)
    }

    // MARK: Potential-curve scan

    var scanRange: ClosedRange<Double> { molecule.scanRange }

    /// Grid of R values over the molecule's scan range (inclusive) with the selected step.
    func scanGrid() -> [Double] {
        let step = scanStep / Units.bohrToAngstrom   // grid is in bohr
        let start = scanRange.lowerBound, end = scanRange.upperBound
        guard step > 0 else { return [] }
        var values: [Double] = []
        let count = Int(((end - start) / step + 1.0e-9).rounded(.down))
        for i in 0...count {
            values.append(start + Double(i) * step)
        }
        if let last = values.last, last < end - 1.0e-9 {
            values.append(end)
        }
        return values
    }

    /// Estimated time of a scan (general molecules), from the single-point estimate plus
    /// about 0.2 s per point for setting up each calculation (basis files, start-up).
    var scanSeconds: Double {
        guard isGeneral, let e = estimate?.estimate else { return 0 }
        return Double(scanGrid().count) * (e.totalSeconds + 0.2)
    }

    /// Run scan button: confirmation first when the scan would take more than 30 s.
    func requestScan() {
        guard !isScanning else { return }
        if isGeneral && energyRefusal != nil {
            scanMessage = energyRefusal
            return
        }
        if isGeneral && scanSeconds > 30 {
            showLongScanConfirmation = true
        } else {
            startScan()
        }
    }

    func startScan() {
        guard !isScanning else { return }
        scanTask = Task { await self.runScan() }
    }

    func cancelScan() {
        scanTask?.cancel()
        scanJob?.cancel()
    }

    private func runScan() async {
        isScanning = true
        scanProgress = 0
        scanMessage = nil
        defer {
            isScanning = false
            scanTask = nil
            scanJob = nil
        }

        let grid = scanGrid()
        var points: [CurvePoint] = []
        var failures: [String] = []
        var template = currentInput(printLevel: 0)
        template.job = "energy"   // a curve is a series of single-point energies, whatever the job
        let general = template.isGeneral

        for (index, r) in grid.enumerated() {
            if Task.isCancelled {
                scanMessage = "Scan cancelled."
                break
            }
            var input = template
            input.R = r
            let job = MacQCJob()
            scanJob = job
            do {
                let output = try await MacQCEngine.run(input, job: job)
                if output.error == "calculation cancelled" {
                    scanMessage = "Scan cancelled."
                    break
                }
                if general, let rhf = output.rhf, output.error == nil {
                    points.append(CurvePoint(R: r, hf: rhf.totalEnergy, correlated: output.correlatedEnergy))
                } else if !general, let hf = output.hf, let fci = output.fci, output.error == nil {
                    points.append(CurvePoint(R: r, hf: hf.totalEnergy, correlated: fci.totalEnergy))
                } else {
                    failures.append(String(format: "R = %.3f: %@", r, output.error ?? "no result"))
                }
            } catch {
                failures.append(String(format: "R = %.3f: %@", r, error.localizedDescription))
            }
            scanProgress = Double(index + 1) / Double(grid.count)
        }

        curvePoints = points
        curveHFLabel = general ? "RHF" : "HF"
        curveCorrelatedLabel = general ? (template.method == "FCI" ? "FCI" : template.method) : "FCI"
        if general {
            let m = MoleculeChoice(rawValue: template.molecule)?.title ?? template.molecule
            curveInputDescription = String(
                format: "%@ · %@ · %@%@", m, template.basis,
                template.method == "RHF" ? "RHF" : "RHF and " + MethodChoice.title(of: template.method),
                template.molecule == "H2O" ? String(format: " · angle H–O–H %.2f°", template.angle) : "")
        } else {
            curveInputDescription = String(
                format: "STO-%ldG, ζ(He) = %.4f, ζ(H) = %.4f, %@ erf",
                template.nGaussians, template.zeta1, template.zeta2,
                template.exactErf ? "exact" : "A&S")
        }
        if !failures.isEmpty {
            scanMessage = "\(failures.count) point(s) failed: " + failures.prefix(3).joined(separator: "; ")
        }
    }

    var curveSamples: [CurveSample] {
        var samples: [CurveSample] = []
        for p in curvePoints {
            samples.append(CurveSample(method: curveHFLabel, R: p.R, energy: p.hf))
        }
        for p in curvePoints {
            if let e = p.correlated {
                samples.append(CurveSample(method: curveCorrelatedLabel, R: p.R, energy: e))
            }
        }
        return samples
    }

    var curveMinima: [CurveMinimum] {
        var minima: [CurveMinimum] = []
        if let m = AppModel.minimum(of: curvePoints.map { ($0.R, $0.hf) }, method: curveHFLabel) {
            minima.append(m)
        }
        let correlatedSamples = curvePoints.compactMap { p in p.correlated.map { (p.R, $0) } }
        if let m = AppModel.minimum(of: correlatedSamples, method: curveCorrelatedLabel) {
            minima.append(m)
        }
        return minima
    }

    /// Lowest grid point; if it is interior and the grid is equally spaced, refined by the
    /// vertex of the parabola through it and its two neighbours.
    static func minimum(of samples: [(Double, Double)], method: String) -> CurveMinimum? {
        guard !samples.isEmpty else { return nil }
        var best = 0
        for i in samples.indices where samples[i].1 < samples[best].1 {
            best = i
        }
        var result = CurveMinimum(method: method, R: samples[best].0, energy: samples[best].1,
                                  refined: false)
        if best > 0 && best < samples.count - 1 {
            let (r0, e0) = samples[best - 1]
            let (r1, e1) = samples[best]
            let (r2, e2) = samples[best + 1]
            let h1 = r1 - r0
            let h2 = r2 - r1
            let curvature = e0 - 2.0 * e1 + e2
            if abs(h1 - h2) < 1.0e-9 && curvature > 0 {
                let shift = -0.5 * h1 * (e2 - e0) / curvature
                let energy = e1 - (e2 - e0) * (e2 - e0) / (8.0 * curvature)
                if abs(shift) <= h1 && energy.isFinite {
                    result.R = r1 + shift
                    result.energy = energy
                    result.refined = true
                }
            }
        }
        return result
    }

    // MARK: Summary text / export

    func summaryText() -> String? {
        guard let out = result else { return nil }
        if out.isGeneral { return generalSummaryText(out) }
        guard let hf = out.hf else { return nil }
        var lines: [String] = []
        lines.append("ΨLab – HeH+ (STO-\(out.input.nGaussians)G)")
        lines.append(String(format: "R = %@, zeta(He) = %.4f, zeta(H) = %.4f, %@ erf",
                            Fmt.lengthBoth(out.input.R), out.input.zeta1, out.input.zeta2,
                            out.input.exactErf ? "exact" : "Abramowitz-Stegun"))
        lines.append("")
        lines.append("Hartree-Fock")
        lines.append("  converged:            \(hf.converged ? "yes" : "no") (\(hf.iterations) iterations)")
        lines.append("  electronic energy:    " + Fmt.energyWithUnit(hf.electronicEnergy))
        lines.append("  nuclear repulsion:    " + Fmt.energyWithUnit(hf.nuclearRepulsion))
        lines.append("  total energy:         " + Fmt.energyWithUnit(hf.totalEnergy))
        if hf.orbitalEnergies.count >= 2 {
            lines.append("  orbital energies:     " + hf.orbitalEnergies.map { Fmt.energy($0) }.joined(separator: ", ") + " E_h")
        }
        if let fci = out.fci {
            lines.append("")
            lines.append("Full CI")
            lines.append("  total energy:         " + Fmt.energyWithUnit(fci.totalEnergy))
            lines.append("  correlation energy:   " + Fmt.energyWithUnit(fci.correlationEnergy))
            for (i, s) in fci.states.enumerated() {
                lines.append(String(format: "  state %ld: %@ E_h  <S^2> = %.4f  %@",
                                    i + 1, Fmt.energy(s.energy), s.S2, s.spin))
            }
            for w in fci.singletWeights {
                lines.append(String(format: "  weight %@: %.6f", w.configuration, w.weight))
            }
        }
        if let message = out.error {
            lines.append("")
            lines.append("Error: \(message)")
        }
        return lines.joined(separator: "\n")
    }

    private func generalSummaryText(_ out: MacQCOutput) -> String? {
        guard let mol = out.molecule, let rhf = out.rhf else { return nil }
        var lines: [String] = []
        lines.append("ΨLab – \(mol.name), \(out.input.method)/\(out.basis?.name ?? out.input.basis)")
        var geometry = String(format: "%@ = %@", MoleculeChoice(rawValue: mol.name)?.bondLabel ?? "R",
                              Fmt.lengthBoth(mol.R))
        if let angle = mol.angle {
            geometry += String(format: ", angle H-O-H = %.2f deg", angle)
        }
        lines.append(geometry)
        if let b = out.basis {
            let atoms = Fmt.uniqueAtoms(b.atoms).map { "\($0.symbol) \($0.contraction)" }.joined(separator: ", ")
            lines.append("\(b.functions) basis functions (\(atoms))")
        }
        for note in out.notes ?? [] {
            lines.append("Note: \(note)")
        }
        lines.append("")
        lines.append("Hartree-Fock (RHF)")
        lines.append("  converged:            \(rhf.converged ? "yes" : "no") (\(rhf.iterations) iterations)")
        lines.append("  total energy:         " + Fmt.energy10(rhf.totalEnergy) + " E_h")
        lines.append("  electronic energy:    " + Fmt.energy10(rhf.electronicEnergy) + " E_h")
        lines.append("  nuclear repulsion:    " + Fmt.energy10(rhf.nuclearRepulsion) + " E_h")
        if let homo = rhf.homo, let lumo = rhf.lumo {
            lines.append(String(format: "  HOMO, LUMO:           %.8f, %.8f E_h", homo, lumo))
        }
        if let f = out.fullCI {
            lines.append("")
            lines.append("Full CI")
            lines.append(String(format: "  determinants:         %@ (%ld orbitals, %ld electrons)",
                                Fmt.count(f.determinants), f.orbitals, f.alphaElectrons + f.betaElectrons))
            lines.append("  converged:            \(f.converged ? "yes" : "no") (\(f.iterations) Davidson iterations)")
            lines.append("  total energy:         " + Fmt.energy10(f.totalEnergy) + " E_h")
            lines.append("  correlation energy:   " + Fmt.energy10(f.correlationEnergy) + " E_h")
            lines.append(String(format: "  <S^2>:                %.6f", f.S2))
            lines.append(String(format: "  RHF determinant:      |C_0|^2 = %.6f", f.hfWeight))
        }
        if let c = out.coupledCluster {
            lines.append("")
            lines.append("Coupled cluster \(c.method)")
            lines.append("  converged:            \(c.converged ? "yes" : "no") (\(c.iterations) iterations)")
            if let t = c.triplesCorrection {
                lines.append("  CCSD energy:          " + Fmt.energy10(c.ccEnergy) + " E_h")
                lines.append("  (T) correction:       " + Fmt.energy10(t) + " E_h")
            }
            lines.append("  total energy:         " + Fmt.energy10(c.totalEnergy) + " E_h")
            lines.append("  correlation energy:   " + Fmt.energy10(c.totalCorrelation) + " E_h")
            lines.append(String(format: "  T1 diagnostic:        %.6f", c.t1Diagnostic))
        }
        if let j = out.job, j.type != "energy" {
            lines.append("")
            lines.append("Job: \(j.title)")
            if let label = j.finalEnergyLabel, let e = j.finalEnergy {
                lines.append("  final energy (\(label)): " + Fmt.energy10(e) + " E_h")
            }
            if let opt = out.optimization {
                lines.append("  optimisation:         \(opt.converged ? "converged" : "NOT converged") in \(opt.steps) steps")
                for p in opt.finalParameters {
                    lines.append("  " + Fmt.parameter(p))
                }
            }
            if let f = out.frequencies {
                let list = f.modes.map { Fmt.wavenumber($0) }.joined(separator: ", ")
                lines.append("  harmonic wavenumbers: " + list + " cm-1 (" + (f.analytic ? "analytic Hessian" : "numerical Hessian") + ")")
            }
            if let t = out.thermochemistry {
                lines.append(String(format: "  thermochemistry at %.2f K, %.0f Pa:", t.temperature, t.pressure))
                lines.append("    ZPE " + String(format: "%.10f E_h", t.corrections.zeroPoint / Units.joulePerMolPerHartree))
                lines.append("    E + ZPE " + Fmt.energy10(t.totals.zeroPoint) + " E_h, H " + Fmt.energy10(t.totals.enthalpy)
                             + " E_h, G " + Fmt.energy10(t.totals.gibbs) + " E_h")
                lines.append(String(format: "    S %.4f J/(mol K)", t.terms.total.entropy))
            }
        }
        for warning in out.warnings ?? [] {
            lines.append("Warning: \(warning)")
        }
        if let message = out.error {
            lines.append("")
            lines.append("Error: \(message)")
        }
        return lines.joined(separator: "\n")
    }

    func copySummary() {
        guard let text = summaryText() else { return }
        let pasteboard = NSPasteboard.general
        pasteboard.clearContents()
        pasteboard.setString(text, forType: .string)
    }

    /// Plain-text export. H₂, LiH, H₂O: the output file of the job, exactly as in the Output
    /// tab and in psilab's output file. HeH⁺: summary, the FORTRAN-style report and, if a scan
    /// was run, the potential-curve table with its minima.
    func fullTextReport() -> String? {
        if let out = result, out.isGeneral, let file = out.output, !file.isEmpty {
            return file
        }
        guard let out = result, let summary = summaryText() else { return nil }
        let rule = String(repeating: "=", count: 78)
        var parts: [String] = [summary]
        if !out.report.isEmpty {
            let title = out.isGeneral
                ? "DETAILED REPORT"
                : "DETAILED REPORT (Hartree-Fock as in Szabo & Ostlund, Appendix B; then Full CI)"
            parts.append(rule + "\n" + title + "\n" + rule + "\n" + out.report)
        }
        if let curve = curveTableText() {
            parts.append(rule + "\n" + curve)
        }
        return parts.joined(separator: "\n\n") + "\n"
    }

    /// The potential-curve scan as a fixed-width table, or nil if no scan has been run.
    func curveTableText() -> String? {
        guard !curvePoints.isEmpty else { return nil }
        var lines: [String] = []
        lines.append("POTENTIAL CURVE")
        if let description = curveInputDescription {
            lines.append(description)
        }
        lines.append("")
        let hfHead = curveHFLabel == "HF" ? " E_HF (E_h)" : "E_RHF (E_h)"
        let hasCorrelated = curvePoints.contains { $0.correlated != nil }
        let corrHead = String(repeating: " ", count: max(0, 11 - curveCorrelatedLabel.count)) + "E_\(curveCorrelatedLabel) (E_h)"
        lines.append("    R (bohr)       R (Å)     \(hfHead)" + (hasCorrelated ? "  " + corrHead : ""))
        for p in curvePoints {
            var line = String(format: "%12.4f  %12.4f  %14.8f", p.R, p.R * Units.bohrToAngstrom, p.hf)
            if let e = p.correlated {
                line += String(format: "  %14.8f", e)
            }
            lines.append(line)
        }
        lines.append("")
        for m in curveMinima {
            lines.append(String(format: "%@ minimum (%@): R = %@, E = %.8f E_h",
                                m.method, m.refined ? "parabolic fit" : "grid point",
                                Fmt.lengthBoth(m.R), m.energy))
        }
        return lines.joined(separator: "\n")
    }

    /// File name stem for the exports, e.g. "PsiLab-HeH+-STO3G-R1.4632" or "PsiLab-H2O-CCSD(T)-aug-cc-pVDZ".
    private func exportName(_ out: MacQCOutput) -> String {
        if out.isGeneral {
            let method = out.input.method
            let job = out.input.job == "energy" ? "" : "-" + out.input.job
            return String(format: "PsiLab-%@-%@-%@%@-R%.4f", out.input.molecule, method, out.input.basis, job, out.input.R)
        }
        return String(format: "PsiLab-HeH+-STO%ldG-R%.4f", out.input.nGaussians, out.input.R)
    }

    func exportText() {
        guard let out = result, let text = fullTextReport() else { return }

        let panel = NSSavePanel()
        panel.title = "Export Results as Text"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = exportName(out) + (out.isGeneral ? ".out" : ".txt")
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try text.write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    func exportJSON() {
        guard let out = result else { return }

        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        let data: Data
        do {
            data = try encoder.encode(out)
        } catch {
            errorMessage = "Could not encode the results: \(error.localizedDescription)"
            return
        }

        let panel = NSSavePanel()
        panel.title = "Export Results as JSON"
        panel.allowedContentTypes = [UTType.json]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = exportName(out) + ".json"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try data.write(to: url, options: .atomic)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    var canExportOrCopy: Bool { result?.hf != nil || result?.rhf != nil || (result?.output?.isEmpty == false) }

    /// File > Export Potential Curve as Text… (the scan table with its minima).
    func exportCurve() {
        guard let text = curveTableText() else { return }
        let panel = NSSavePanel()
        panel.title = "Export Potential Curve as Text"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = "PsiLab-curve-\(molecule.rawValue).txt"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try (text + "\n").write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    /// Selects a tab (and its page) from the View menu.
    func show(_ tab: DetailTab, methods: MethodsPage? = nil, vibration: VibrationPage? = nil) {
        if let m = methods { methodsPage = m }
        if let v = vibration { vibrationPage = v }
        selectedTab = tab
    }
}
