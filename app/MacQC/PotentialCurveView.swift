//
//  PotentialCurveView.swift
//  MacQC
//
//  Potential-curve tab: scan the bond length (HeH⁺: 0.8 … 5.0 bohr; H₂, LiH, H₂O: their own
//  ranges, H₂O as a symmetric stretch), plot the HF/RHF and Full CI total energies with
//  Swift Charts and mark the minimum of each curve.
//

import SwiftUI
import Charts

struct PotentialCurveView: View {
    @Bindable var model: AppModel

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 16) {
                controls

                if let message = model.scanMessage {
                    ErrorBanner(message: message)
                }

                if model.curvePoints.isEmpty {
                    if !model.isScanning {
                        ContentUnavailableView(
                            "No potential curve yet",
                            systemImage: "chart.xyaxis.line",
                            description: Text(emptyText)
                        )
                        .frame(maxWidth: .infinity, minHeight: 260)
                    }
                } else {
                    chart
                    minimaPanel
                }
            }
            .padding(20)
        }
    }

    // MARK: Controls

    private var controls: some View {
        GlassGroup(spacing: 12) {
            HStack(spacing: 14) {
                Picker("Step (Å)", selection: $model.scanStep) {
                    ForEach(AppModel.scanStepChoices, id: \.self) { step in
                        Text(String(format: "%.3f", step)).tag(step)
                    }
                }
                .pickerStyle(.menu)
                .frame(width: 190)
                .disabled(model.isScanning)

                Button {
                    model.requestScan()
                } label: {
                    Label("Run scan", systemImage: "play.fill")
                }
                .prominentActionStyle()
                .accessibilityLabel("Run scan")
                .disabled(model.isScanning)

                if model.isScanning {
                    Button("Stop") {
                        model.cancelScan()
                    }
                    .secondaryActionStyle()
                    .accessibilityLabel("Stop scan")

                    ProgressView(value: model.scanProgress)
                        .progressViewStyle(.linear)
                        .frame(width: 160)
                }

                Spacer()

                Text(pointsText)
                    .foregroundStyle(.secondary)
                    .monospacedDigit()
            }
            .padding(12)
            .glassPanel(cornerRadius: 18)
        }
        .confirmationDialog("Long scan", isPresented: $model.showLongScanConfirmation) {
            Button("Run scan") {
                model.startScan()
            }
            Button("Cancel", role: .cancel) {}
        } message: {
            Text("The scan has \(model.scanGrid().count) points and is estimated to take \(Fmt.duration(model.scanSeconds)). You can stop it while it runs.")
        }
    }

    private var pointsText: String {
        let n = model.scanGrid().count
        if model.isGeneral && model.scanSeconds > 0 {
            return "\(n) points · \(Fmt.duration(model.scanSeconds))"
        }
        return "\(n) points"
    }

    private var emptyText: String {
        let range = model.scanRange
        let rangeText = String(format: "R = %.2f … %.2f Å (%.2f … %.2f bohr)",
                               range.lowerBound * Units.bohrToAngstrom, range.upperBound * Units.bohrToAngstrom,
                               range.lowerBound, range.upperBound)
        if model.isGeneral {
            let method = model.method == .fci ? "RHF and Full CI" : "RHF"
            let what = model.molecule == .h2o ? "symmetric stretch, \(rangeText) at the chosen angle" : rangeText
            return "Press Run scan to compute \(method) total energies of \(model.molecule.title) (\(what)) in \(model.basis)."
        }
        return "Press Run scan to compute HF and FCI total energies for \(rangeText) with the basis and exponents from the sidebar."
    }

    // MARK: Chart

    private var chart: some View {
        VStack(alignment: .leading, spacing: 8) {
            if let description = model.curveInputDescription {
                Text(description)
                    .font(.subheadline)
                    .foregroundStyle(.secondary)
            }

            Chart {
                ForEach(model.curveSamples) { sample in
                    LineMark(
                        x: .value("R (Å)", sample.R * Units.bohrToAngstrom),
                        y: .value("Total energy (E_h)", sample.energy),
                        series: .value("Method", sample.method)
                    )
                    .foregroundStyle(by: .value("Method", sample.method))
                    .interpolationMethod(.catmullRom)
                }

                ForEach(model.curveMinima) { minimum in
                    PointMark(
                        x: .value("R (Å)", minimum.R * Units.bohrToAngstrom),
                        y: .value("Total energy (E_h)", minimum.energy)
                    )
                    .foregroundStyle(by: .value("Method", minimum.method))
                    .symbol(by: .value("Method", minimum.method))
                    .symbolSize(110)
                }
            }
            .chartForegroundStyleScale([model.curveHFLabel: Color.orange, "FCI": Color.blue])
            // Same shapes as the minimum markers, so the legend matches the plot.
            .chartSymbolScale([model.curveHFLabel: BasicChartSymbolShape.diamond, "FCI": BasicChartSymbolShape.circle])
            .chartXScale(domain: xDomain)
            .chartYScale(domain: yDomain)
            .chartXAxisLabel("R (Å)")
            .chartYAxisLabel("Total energy (E_h)")
            .chartLegend(position: .top, alignment: .trailing)
            .chartOverlay { proxy in
                // One readout box in an empty corner of the plot area, so the labels can
                // neither overlap each other nor the curves.
                GeometryReader { geometry in
                    if let plotFrame = proxy.plotFrame, !model.curveMinima.isEmpty {
                        let frame = geometry[plotFrame]
                        minimaReadout
                            .padding(10)
                            .frame(width: frame.width, height: frame.height, alignment: readoutAlignment)
                            .offset(x: frame.minX, y: frame.minY)
                            .allowsHitTesting(false)
                    }
                }
            }
            .frame(minHeight: 380)
        }
        .dataPanel()
    }

    // MARK: Minima markers and readout

    private func color(for method: String) -> Color {
        method == "FCI" ? .blue : .orange
    }

    /// Bottom-right normally (the curves level off high up at large R); bottom-left if a
    /// minimum lies in the right half of the scan.
    private var readoutAlignment: Alignment {
        let middle = (xDomain.lowerBound + xDomain.upperBound) / 2
        let anyOnRight = model.curveMinima.contains { $0.R * Units.bohrToAngstrom > middle }
        return anyOnRight ? .bottomLeading : .bottomTrailing
    }

    private var minimaReadout: some View {
        VStack(alignment: .leading, spacing: 4) {
            ForEach(model.curveMinima) { minimum in
                HStack(spacing: 6) {
                    Image(systemName: minimum.method == "FCI" ? "circle.fill" : "diamond.fill")
                        .font(.caption2)
                        .foregroundStyle(color(for: minimum.method))
                    Text(String(format: "%@ min: R = %@, E = %.6f E_h",
                                minimum.method, Fmt.lengthBoth(minimum.R, decimals: 3), minimum.energy))
                        .font(.caption.monospacedDigit())
                }
            }
        }
        .padding(.horizontal, 8)
        .padding(.vertical, 6)
        .background(.regularMaterial, in: RoundedRectangle(cornerRadius: 8))
        .fixedSize()
    }

    // MARK: Axis ranges

    /// Swift Charts would otherwise include 0 on the energy axis and squash the curves.
    private var xDomain: ClosedRange<Double> {
        let values = model.curvePoints.map { $0.R * Units.bohrToAngstrom }   // axis in Å
        guard let lo = values.min(), let hi = values.max(), lo < hi else {
            return (model.scanRange.lowerBound * Units.bohrToAngstrom)...(model.scanRange.upperBound * Units.bohrToAngstrom)
        }
        return lo...hi
    }

    private var yDomain: ClosedRange<Double> {
        let values = model.curveSamples.map { $0.energy }
        guard let lo = values.min(), let hi = values.max(), lo < hi else { return -3.0 ... -2.0 }
        let pad = max((hi - lo) * 0.12, 0.005)
        return (lo - pad)...(hi + pad)
    }

    // MARK: Minima

    private var minimaPanel: some View {
        VStack(alignment: .leading, spacing: 8) {
            SectionTitle("Minima")
            ForEach(model.curveMinima) { minimum in
                InfoRow(
                    "\(minimum.method) (\(minimum.refined ? "parabolic fit" : "grid point"))",
                    String(format: "R = %@,  E = %@ E_h", Fmt.lengthBoth(minimum.R), Fmt.energy(minimum.energy))
                )
            }
            Text("The minimum is the lowest scanned point, refined by a parabola through it and its two neighbours when it is not at the edge of the scan range.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }
}
