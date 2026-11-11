//
//  PotentialCurveView.swift
//  MacQC
//
//  Potential-curve tab: scan R from 0.8 to 5.0 bohr, plot HF and FCI total energies with
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
                            description: Text("Press Run scan to compute HF and FCI total energies for R = 0.8 … 5.0 bohr with the basis and exponents from the sidebar.")
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
                Picker("Step (bohr)", selection: $model.scanStep) {
                    ForEach(AppModel.scanStepChoices, id: \.self) { step in
                        Text(String(format: "%.2f", step)).tag(step)
                    }
                }
                .pickerStyle(.menu)
                .frame(width: 190)
                .disabled(model.isScanning)

                Button {
                    model.startScan()
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

                Text("\(model.scanGrid().count) points")
                    .foregroundStyle(.secondary)
                    .monospacedDigit()
            }
            .padding(12)
            .glassPanel(cornerRadius: 18)
        }
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
                        x: .value("R (bohr)", sample.R),
                        y: .value("Total energy (E_h)", sample.energy),
                        series: .value("Method", sample.method)
                    )
                    .foregroundStyle(by: .value("Method", sample.method))
                    .interpolationMethod(.catmullRom)
                }

                ForEach(model.curveMinima) { minimum in
                    PointMark(
                        x: .value("R (bohr)", minimum.R),
                        y: .value("Total energy (E_h)", minimum.energy)
                    )
                    .foregroundStyle(by: .value("Method", minimum.method))
                    .symbol(.diamond)
                    .symbolSize(140)
                    .annotation(position: minimum.method == "HF" ? .bottom : .top, spacing: 6) {
                        Text(String(format: "%@ min: R = %.3f, E = %.6f", minimum.method, minimum.R, minimum.energy))
                            .font(.caption2.monospacedDigit())
                            .padding(.horizontal, 6)
                            .padding(.vertical, 2)
                            .background(.regularMaterial, in: RoundedRectangle(cornerRadius: 6))
                    }
                }
            }
            .chartForegroundStyleScale(["HF": Color.blue, "FCI": Color.orange])
            .chartXScale(domain: xDomain)
            .chartYScale(domain: yDomain)
            .chartXAxisLabel("R (bohr)")
            .chartYAxisLabel("Total energy (E_h)")
            .chartLegend(position: .top, alignment: .trailing)
            .frame(minHeight: 380)
        }
        .dataPanel()
    }

    // MARK: Axis ranges

    /// Swift Charts would otherwise include 0 on the energy axis and squash the curves.
    private var xDomain: ClosedRange<Double> {
        let values = model.curvePoints.map { $0.R }
        guard let lo = values.min(), let hi = values.max(), lo < hi else { return 0.8...5.0 }
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
                    String(format: "R = %.4f bohr,  E = %@ E_h", minimum.R, Fmt.energy(minimum.energy))
                )
            }
            Text("The minimum is the lowest scanned point, refined by a parabola through it and its two neighbours when it is not at the edge of the scan range.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }
}
