//
//  AppModel.swift
//  MacQC
//
//  Observable application state: sidebar inputs, size estimate, the running job, last result,
//  potential-curve scan, and the menu actions (export text / JSON, copy summary).
//

import SwiftUI
import AppKit
import UniformTypeIdentifiers

enum DetailTab: Hashable {
    case summary, hartreeFock, fullCI, potentialCurve, report
}

/// One point of the potential-curve scan (FCI is missing in an RHF-only scan).
struct CurvePoint: Identifiable, Sendable {
    var R: Double
    var hf: Double
    var fci: Double?
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
    var method: MethodChoice = .rhf
    var bondLengths: [MoleculeChoice: Double] = [:]   // bohr; missing = experimental value
    var waterAngle: Double = MoleculeChoice.defaultWaterAngle

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
        return "\(molecule.rawValue)|\(basis)|\(method.rawValue)|\(bondLength)|\(waterAngle)"
    }

    /// The backend refuses the calculation (reason), or the inputs are invalid.
    var refusal: String? {
        guard isGeneral, let e = estimate else { return nil }
        if let r = e.refusal { return r }
        if e.estimate == nil, let message = e.error { return message }
        return nil
    }

    var isLongCalculation: Bool { estimate?.estimate?.longCalculation ?? false }

    // MARK: Results and the running job

    var selectedTab: DetailTab = .summary
    var isRunning = false
    var result: MacQCOutput?
    var errorMessage: String?
    var progress = MacQCProgress()
    var runStarted: Date?
    var showLongRunConfirmation = false

    @ObservationIgnored private var job: MacQCJob?

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
    var showLongScanConfirmation = false

    @ObservationIgnored private var scanTask: Task<Void, Never>?
    @ObservationIgnored private var scanJob: MacQCJob?

    // MARK: Start-up

    init() {
        // Launch arguments (also used by the UI test): -psilabMolecule H2O -psilabBasis aug-cc-pVTZ
        // -psilabMethod FCI
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
           let choice = MethodChoice.allCases.first(where: { $0.rawValue.lowercased() == m.lowercased() }) {
            method = choice
        }
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
        estimate?.warnings?.joined(separator: "\n") ?? "This calculation may take a long time."
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
        runStarted = Date()
        let job = MacQCJob()
        self.job = job
        // Poll the backend's progress while the calculation runs.
        let poller = Task { [weak self] in
            while !Task.isCancelled {
                let p = job.progress()
                self?.progress = p
                try? await Task.sleep(for: .milliseconds(250))
            }
        }
        defer {
            poller.cancel()
            isRunning = false
            self.job = nil
            runStarted = nil
        }

        do {
            let output = try await MacQCEngine.run(currentInput(printLevel: 2), job: job)
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

    func cancelRun() {
        job?.cancel()
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
        if isGeneral && refusal != nil {
            scanMessage = refusal
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
        let template = currentInput(printLevel: 0)
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
                    points.append(CurvePoint(R: r, hf: rhf.totalEnergy, fci: output.fullCI?.totalEnergy))
                } else if !general, let hf = output.hf, let fci = output.fci, output.error == nil {
                    points.append(CurvePoint(R: r, hf: hf.totalEnergy, fci: fci.totalEnergy))
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
        if general {
            let m = MoleculeChoice(rawValue: template.molecule)?.title ?? template.molecule
            curveInputDescription = String(
                format: "%@ · %@ · %@%@", m, template.basis,
                template.method == "FCI" ? "RHF and Full CI" : "RHF",
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
            if let fci = p.fci {
                samples.append(CurveSample(method: "FCI", R: p.R, energy: fci))
            }
        }
        return samples
    }

    var curveMinima: [CurveMinimum] {
        var minima: [CurveMinimum] = []
        if let m = AppModel.minimum(of: curvePoints.map { ($0.R, $0.hf) }, method: curveHFLabel) {
            minima.append(m)
        }
        let fciSamples = curvePoints.compactMap { p in p.fci.map { (p.R, $0) } }
        if let m = AppModel.minimum(of: fciSamples, method: "FCI") {
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
        let method = out.fullCI != nil ? "FCI" : "RHF"
        lines.append("ΨLab – \(mol.name), \(method)/\(out.basis?.name ?? out.input.basis)")
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

    /// Plain-text export: summary, the backend's full report and, if a scan was run, the
    /// potential-curve table with its minima.
    func fullTextReport() -> String? {
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
        let hasFCI = curvePoints.contains { $0.fci != nil }
        lines.append("    R (bohr)       R (Å)     \(hfHead)" + (hasFCI ? "     E_FCI (E_h)" : ""))
        for p in curvePoints {
            var line = String(format: "%12.4f  %12.4f  %14.8f", p.R, p.R * Units.bohrToAngstrom, p.hf)
            if let fci = p.fci {
                line += String(format: "  %14.8f", fci)
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

    /// File name stem for the exports, e.g. "PsiLab-HeH+-STO3G-R1.4632" or "PsiLab-H2O-FCI-aug-cc-pVDZ".
    private func exportName(_ out: MacQCOutput) -> String {
        if out.isGeneral {
            let method = out.fullCI != nil ? "FCI" : "RHF"
            return String(format: "PsiLab-%@-%@-%@-R%.4f", out.input.molecule, method, out.input.basis, out.input.R)
        }
        return String(format: "PsiLab-HeH+-STO%ldG-R%.4f", out.input.nGaussians, out.input.R)
    }

    func exportText() {
        guard let out = result, let text = fullTextReport() else { return }

        let panel = NSSavePanel()
        panel.title = "Export Results as Text"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = exportName(out) + ".txt"
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

    var canExportOrCopy: Bool { result?.hf != nil || result?.rhf != nil }
}
