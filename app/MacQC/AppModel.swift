//
//  AppModel.swift
//  MacQC
//
//  Observable application state: sidebar inputs, last result, potential-curve scan,
//  and the menu actions (export JSON, copy summary).
//

import SwiftUI
import AppKit
import UniformTypeIdentifiers

enum DetailTab: Hashable {
    case summary, hartreeFock, fullCI, potentialCurve
}

/// One point of the potential-curve scan.
struct CurvePoint: Identifiable, Sendable {
    var R: Double
    var hf: Double
    var fci: Double
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

    var R: Double = AppModel.defaultR
    var zeta1: Double = AppModel.defaultZeta1   // He
    var zeta2: Double = AppModel.defaultZeta2   // H
    var nGaussians: Int = AppModel.defaultGaussians
    var exactErf: Bool = false

    // MARK: Results

    var selectedTab: DetailTab = .summary
    var isRunning = false
    var result: MacQCOutput?
    var errorMessage: String?

    // MARK: Potential curve

    static let scanStart = 0.8
    static let scanEnd = 5.0
    static let scanStepChoices: [Double] = [0.05, 0.1, 0.15, 0.2, 0.3]

    var scanStep: Double = 0.1
    var isScanning = false
    var scanProgress: Double = 0
    var scanMessage: String?
    var curvePoints: [CurvePoint] = []
    var curveInputDescription: String?   // parameters the scan was computed with

    @ObservationIgnored private var scanTask: Task<Void, Never>?

    // MARK: Inputs

    func currentInput(printLevel: Int = 2) -> MacQCInput {
        var input = MacQCInput()
        input.R = R
        input.zeta1 = zeta1
        input.zeta2 = zeta2
        input.ZA = 2
        input.ZB = 1
        input.nGaussians = nGaussians
        input.exactErf = exactErf
        input.printLevel = printLevel
        return input
    }

    func resetToDefaults() {
        R = AppModel.defaultR
        zeta1 = AppModel.defaultZeta1
        zeta2 = AppModel.defaultZeta2
        nGaussians = AppModel.defaultGaussians
        exactErf = false
    }

    // MARK: Single calculation

    func run() async {
        guard !isRunning else { return }
        isRunning = true
        errorMessage = nil
        defer { isRunning = false }

        do {
            let output = try await MacQCEngine.run(currentInput(printLevel: 2))
            result = output
            errorMessage = output.error
        } catch {
            // Keep the previous result visible; show why the new run failed.
            errorMessage = error.localizedDescription
        }
    }

    // MARK: Potential-curve scan

    /// Grid of R values from 0.8 to 5.0 bohr (inclusive) with the selected step.
    func scanGrid() -> [Double] {
        let step = scanStep
        guard step > 0 else { return [] }
        var values: [Double] = []
        let count = Int(((AppModel.scanEnd - AppModel.scanStart) / step + 1.0e-9).rounded(.down))
        for i in 0...count {
            values.append(AppModel.scanStart + Double(i) * step)
        }
        if let last = values.last, last < AppModel.scanEnd - 1.0e-9 {
            values.append(AppModel.scanEnd)
        }
        return values
    }

    func startScan() {
        guard !isScanning else { return }
        scanTask = Task { await self.runScan() }
    }

    func cancelScan() {
        scanTask?.cancel()
    }

    private func runScan() async {
        isScanning = true
        scanProgress = 0
        scanMessage = nil
        defer {
            isScanning = false
            scanTask = nil
        }

        let grid = scanGrid()
        var points: [CurvePoint] = []
        var failures: [String] = []
        let template = currentInput(printLevel: 0)

        for (index, r) in grid.enumerated() {
            if Task.isCancelled {
                scanMessage = "Scan cancelled."
                break
            }
            var input = template
            input.R = r
            do {
                let output = try await MacQCEngine.run(input)
                if let hf = output.hf, let fci = output.fci, output.error == nil {
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
        curveInputDescription = String(
            format: "STO-%ldG, ζ(He) = %.4f, ζ(H) = %.4f, %@ erf",
            template.nGaussians, template.zeta1, template.zeta2,
            template.exactErf ? "exact" : "A&S")
        if !failures.isEmpty {
            scanMessage = "\(failures.count) point(s) failed: " + failures.prefix(3).joined(separator: "; ")
        }
    }

    var curveSamples: [CurveSample] {
        var samples: [CurveSample] = []
        for p in curvePoints {
            samples.append(CurveSample(method: "HF", R: p.R, energy: p.hf))
        }
        for p in curvePoints {
            samples.append(CurveSample(method: "FCI", R: p.R, energy: p.fci))
        }
        return samples
    }

    var curveMinima: [CurveMinimum] {
        var minima: [CurveMinimum] = []
        if let m = AppModel.minimum(of: curvePoints.map { ($0.R, $0.hf) }, method: "HF") {
            minima.append(m)
        }
        if let m = AppModel.minimum(of: curvePoints.map { ($0.R, $0.fci) }, method: "FCI") {
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
        guard let out = result, let hf = out.hf else { return nil }
        var lines: [String] = []
        lines.append("ΨLab – HeH+ (STO-\(out.input.nGaussians)G)")
        lines.append(String(format: "R = %.4f bohr, zeta(He) = %.4f, zeta(H) = %.4f, %@ erf",
                            out.input.R, out.input.zeta1, out.input.zeta2,
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

    func copySummary() {
        guard let text = summaryText() else { return }
        let pasteboard = NSPasteboard.general
        pasteboard.clearContents()
        pasteboard.setString(text, forType: .string)
    }

    /// Plain-text export: summary, the backend's full FORTRAN-style report (HF + Full CI) and,
    /// if a scan was run, the potential-curve table with its minima.
    func fullTextReport() -> String? {
        guard let out = result, let summary = summaryText() else { return nil }
        let rule = String(repeating: "=", count: 78)
        var parts: [String] = [summary]
        if !out.report.isEmpty {
            parts.append(rule + "\nDETAILED REPORT (Hartree-Fock as in Szabo & Ostlund, Appendix B; then Full CI)\n"
                         + rule + "\n" + out.report)
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
        lines.append("    R (bohr)      E_HF (E_h)     E_FCI (E_h)")
        for p in curvePoints {
            lines.append(String(format: "%12.4f  %14.8f  %14.8f", p.R, p.hf, p.fci))
        }
        lines.append("")
        for m in curveMinima {
            lines.append(String(format: "%@ minimum (%@): R = %.4f bohr, E = %.8f E_h",
                                m.method, m.refined ? "parabolic fit" : "grid point", m.R, m.energy))
        }
        return lines.joined(separator: "\n")
    }

    func exportText() {
        guard let out = result, let text = fullTextReport() else { return }

        let panel = NSSavePanel()
        panel.title = "Export Results as Text"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = String(format: "PsiLab-HeH+-STO%ldG-R%.4f.txt",
                                            out.input.nGaussians, out.input.R)
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
        panel.nameFieldStringValue = String(format: "PsiLab-HeH+-STO%ldG-R%.4f.json",
                                            out.input.nGaussians, out.input.R)
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try data.write(to: url, options: .atomic)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    var canExportOrCopy: Bool { result?.hf != nil }
}
