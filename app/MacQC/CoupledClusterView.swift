//
//  CoupledClusterView.swift
//  MacQC
//
//  Coupled-cluster tab (H₂, LiH, H₂O): energies, (T) correction, amplitude statistics, the
//  largest amplitudes and the convergence of the iterations.
//

import SwiftUI
import Charts

struct CoupledClusterView: View {
    var model: AppModel

    var body: some View {
        if let out = model.result, out.isGeneral, let cc = out.coupledCluster {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }
                    SectionTitle("Energies (E_h)")
                    energies(cc)

                    SectionTitle("Amplitudes")
                    amplitudes(cc)

                    SectionTitle("Largest amplitudes")
                    amplitudeTable(cc)

                    SectionTitle("Iterations")
                    convergenceChart(cc)
                    iterationTable(cc)
                }
                .padding(20)
            }
        } else {
            VStack(spacing: 14) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                        .frame(maxWidth: 560)
                }
                ContentUnavailableView(
                    "No coupled-cluster result",
                    systemImage: "point.3.connected.trianglepath.dotted",
                    description: Text(placeholderText)
                )
            }
            .padding()
        }
    }

    private var placeholderText: String {
        if !model.isGeneral {
            return "Coupled cluster is available for H₂, LiH and H₂O. HeH⁺ runs Hartree–Fock and Full CI."
        }
        if model.method.isCoupledCluster, let refusal = model.refusal {
            return refusal
        }
        if let out = model.result, out.isGeneral, !(MethodChoice(rawValue: out.input.method)?.isCoupledCluster ?? false) {
            return "This result is \(out.methodTitle). Choose CCSD, CCSD(T), CCSDT or CCSDTQ in the sidebar and press Run."
        }
        return "Choose CCSD, CCSD(T), CCSDT or CCSDTQ in the sidebar and press Run."
    }

    private func energies(_ c: MacQCCoupledCluster) -> some View {
        let base = c.triplesCorrection != nil ? "CCSD" : c.method
        return VStack(spacing: 8) {
            InfoRow("RHF energy", Fmt.energy10(c.referenceEnergy))
            InfoRow("\(base) correlation energy", Fmt.energy10(c.correlationEnergy))
            InfoRow("\(base) energy", Fmt.energy10(c.ccEnergy))
            if let t = c.triplesCorrection {
                InfoRow("(T) correction", Fmt.energy10(t))
                InfoRow("CCSD(T) correlation energy", Fmt.energy10(c.totalCorrelation))
                InfoRow("CCSD(T) energy", Fmt.energy10(c.totalEnergy))
            }
            InfoRow("T1 diagnostic", String(format: "%.6f", c.t1Diagnostic))
        }
        .dataPanel()
    }

    private func amplitudes(_ c: MacQCCoupledCluster) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            ForEach(1..<max(c.amplitudeCounts.count, 1), id: \.self) { rank in
                InfoRow("T\(rank)", Fmt.count(c.amplitudeCounts[rank]) + String(
                    format: " amplitudes, |T%ld| = %.6f", rank, rank < c.amplitudeNorms.count ? c.amplitudeNorms[rank] : 0))
            }
            InfoRow("Orbitals", "\(c.occupied) occupied, \(c.virtuals) virtual (all electrons correlated)")
            InfoRow("Implementation", c.engine)
            if let ladder = c.ladder {
                InfoRow("(vv|vv) ladder", ladder)
            }
            Text(c.engine.contains("spin-orbital")
                 ? "Amplitudes are stored once per set of distinct spin orbitals (antisymmetric, blocked by spin)."
                 : "Closed-shell amplitudes over spatial orbitals; T2 is the αβ amplitude t_ij^ab.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func orbitals(_ list: [Int]) -> String {
        list.map(String.init).joined(separator: " ")
    }

    private func amplitudeTable(_ c: MacQCCoupledCluster) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            Grid(alignment: .trailing, horizontalSpacing: 26, verticalSpacing: 5) {
                GridRow {
                    Text("Rank").gridColumnAlignment(.leading)
                    Text("Amplitude")
                    Text("Occupied").gridColumnAlignment(.leading)
                    Text("Virtual").gridColumnAlignment(.leading)
                    Text("Spins").gridColumnAlignment(.leading)
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(Array(c.leadingAmplitudes.enumerated()), id: \.offset) { _, a in
                    GridRow {
                        Text("T\(a.rank)")
                        Text(String(format: "%+.8f", a.value))
                        Text(orbitals(a.occupied))
                        Text(orbitals(a.virtuals))
                        Text(a.spins.map { $0 == "a" ? "α" : "β" }.joined(separator: " "))
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text("Molecular orbitals numbered by RHF orbital energy; the spins are those of the occupied (and the matching virtual) orbitals. Spin-flipped copies of an amplitude are listed once.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    /// log10 of the largest amplitude change per iteration on a linear axis.
    private func convergenceChart(_ c: MacQCCoupledCluster) -> some View {
        let points: [(iteration: Int, log: Double)] = c.history.map {
            (iteration: $0.iteration, log: log10(max($0.residual, 1.0e-16)))
        }
        let lowest = (points.map { $0.log }.min() ?? -8).rounded(.down)
        let highest = (points.map { $0.log }.max() ?? 0).rounded(.up)
        return VStack(alignment: .leading, spacing: 8) {
            Chart {
                ForEach(points, id: \.iteration) { p in
                    LineMark(x: .value("Iteration", p.iteration), y: .value("log10 max |Δt|", p.log))
                        .foregroundStyle(Color.blue.opacity(0.6))
                    PointMark(x: .value("Iteration", p.iteration), y: .value("log10 max |Δt|", p.log))
                        .foregroundStyle(Color.blue)
                        .symbolSize(40)
                }
            }
            .chartYScale(domain: min(lowest, -8.0)...max(highest, 0.0))
            .chartXAxisLabel("Iteration")
            .chartYAxisLabel("log₁₀ max |Δt|")
            .frame(height: 200)
            Text(String(format: "%@ in %ld iterations (%.2f s), largest amplitude change below 1e-8 and |ΔE| below 1e-10 E_h; DIIS extrapolation.",
                        c.converged ? "Converged" : "Not converged", c.iterations, c.ccSeconds))
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func iterationTable(_ c: MacQCCoupledCluster) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 28, verticalSpacing: 5) {
            GridRow {
                Text("Iteration").gridColumnAlignment(.leading)
                Text("Total energy (E_h)")
                Text("ΔE (E_h)")
                Text("max |Δt|")
            }
            .font(.subheadline.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider().gridCellUnsizedAxes(.horizontal)
            ForEach(c.history) { step in
                GridRow {
                    Text("\(step.iteration)")
                    Text(Fmt.energy10(step.totalEnergy))
                    Text(Fmt.scientific(step.deltaE))
                    Text(Fmt.scientific(step.residual))
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }
}

/// Coupled-cluster rows of the summary tab.
struct CoupledClusterSummaryRows: View {
    var cc: MacQCCoupledCluster

    var body: some View {
        VStack(spacing: 8) {
            InfoRow("Converged", cc.converged ? "yes (\(cc.iterations) iterations)" : "no")
            if let t = cc.triplesCorrection {
                InfoRow("CCSD energy", Fmt.energyWithUnit(cc.ccEnergy))
                InfoRow("(T) correction", Fmt.energyWithUnit(t))
            }
            InfoRow("\(cc.method) energy", Fmt.energyWithUnit(cc.totalEnergy))
            InfoRow("Correlation energy", Fmt.energyWithUnit(cc.totalCorrelation))
            InfoRow("T1 diagnostic", String(format: "%.6f", cc.t1Diagnostic))
            InfoRow("Amplitudes", Fmt.count(cc.amplitudeCounts.reduce(0, +)) + " (\(cc.occupied) occupied, \(cc.virtuals) virtual orbitals)")
            InfoRow("Time", String(format: "%.2f s (integral transformation %.2f s, iterations %.2f s%@)",
                                   cc.transformSeconds + cc.ccSeconds + (cc.triplesSeconds ?? 0), cc.transformSeconds,
                                   cc.ccSeconds, cc.triplesSeconds.map { String(format: ", (T) %.2f s", $0) } ?? ""))
        }
    }
}
