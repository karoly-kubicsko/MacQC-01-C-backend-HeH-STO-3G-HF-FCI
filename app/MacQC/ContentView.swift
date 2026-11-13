//
//  ContentView.swift
//  MacQC
//
//  NavigationSplitView: input sidebar + tabbed detail area.
//

import SwiftUI

struct ContentView: View {
    @Bindable var model: AppModel
    @AppStorage(AppearanceChoice.storageKey) private var appearanceRaw: String = AppearanceChoice.system.rawValue

    var body: some View {
        NavigationSplitView {
            SidebarView(model: model)
                .navigationSplitViewColumnWidth(min: 260, ideal: 290, max: 360)
        } detail: {
            DetailView(model: model)
        }
        .navigationTitle("ΨLab")
        .preferredColorScheme((AppearanceChoice(rawValue: appearanceRaw) ?? .system).colorScheme)
        .task {
            // Show a result right away (the Szabo–Ostlund defaults, or the molecule given as a
            // launch argument). Large general calculations wait for Run.
            if model.result == nil {
                if model.isGeneral {
                    await model.refreshEstimate()
                    if model.refusal == nil && (!model.isLongCalculation || model.autoRunLongJobs) {
                        await model.run()
                    }
                } else {
                    await model.run()
                }
            }
        }
    }
}

/// Six tabs (they fit next to the sidebar at the minimum window width of 1000 pt): Summary,
/// Methods (Hartree–Fock, coupled cluster, Full CI), Optimisation, Frequencies (with the
/// thermochemistry page), Potential curve, Output.
struct DetailView: View {
    @Bindable var model: AppModel

    var body: some View {
        TabView(selection: $model.selectedTab) {
            SummaryView(model: model)
                .tabItem { Label("Summary", systemImage: "list.bullet.rectangle") }
                .tag(DetailTab.summary)

            MethodsView(model: model)
                .tabItem { Label("Methods", systemImage: "square.grid.2x2") }
                .tag(DetailTab.methods)

            OptimizationView(model: model)
                .tabItem { Label("Optimisation", systemImage: "arrow.down.right.and.arrow.up.left") }
                .tag(DetailTab.optimization)

            FrequenciesView(model: model)
                .tabItem { Label("Frequencies", systemImage: "waveform.path") }
                .tag(DetailTab.frequencies)

            PotentialCurveView(model: model)
                .tabItem { Label("Curve", systemImage: "chart.xyaxis.line") }
                .tag(DetailTab.potentialCurve)

            OutputView(model: model)
                .tabItem { Label("Output", systemImage: "doc.plaintext") }
                .tag(DetailTab.output)
        }
        .padding(.top, 4)
    }
}

// MARK: - Shared small views

/// Banner in the detail area before a result exists: why the current settings cannot be calculated (red), or that the
/// calculation will take long (orange).
struct StatusBanner: View {
    enum Kind {
        case refusal, warning
    }

    var kind: Kind
    var text: String

    var body: some View {
        let color: Color = kind == .refusal ? .red : .orange
        Label {
            VStack(alignment: .leading, spacing: 2) {
                Text(kind == .refusal ? "Not possible with these settings" : "Long calculation")
                    .font(.headline)
                Text(text)
                    .textSelection(.enabled)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
        } icon: {
            Image(systemName: kind == .refusal ? "xmark.octagon.fill" : "clock.badge.exclamationmark")
                .foregroundStyle(color)
        }
        .padding(12)
        .background(color.opacity(0.12), in: RoundedRectangle(cornerRadius: 10))
        .overlay(RoundedRectangle(cornerRadius: 10).strokeBorder(color.opacity(0.35), lineWidth: 1))
    }
}

/// Red banner for errors returned by the backend (or by the wrapper).
struct ErrorBanner: View {
    var message: String

    var body: some View {
        Label {
            Text(message)
                .textSelection(.enabled)
                .frame(maxWidth: .infinity, alignment: .leading)
        } icon: {
            Image(systemName: "exclamationmark.triangle.fill")
                .foregroundStyle(.red)
        }
        .padding(12)
        .background(Color.red.opacity(0.12), in: RoundedRectangle(cornerRadius: 10))
        .overlay(RoundedRectangle(cornerRadius: 10).strokeBorder(Color.red.opacity(0.35), lineWidth: 1))
    }
}

/// Section heading used inside the detail tabs.
struct SectionTitle: View {
    var text: String

    init(_ text: String) {
        self.text = text
    }

    var body: some View {
        Text(text)
            .font(.title3.weight(.semibold))
    }
}

/// Placeholder when no result is available yet.
struct NoResultView: View {
    var model: AppModel

    var body: some View {
        VStack(spacing: 14) {
            if let message = model.errorMessage {
                ErrorBanner(message: message)
                    .frame(maxWidth: 560)
            }
            if model.isGeneral, let refusal = model.refusal {
                StatusBanner(kind: .refusal, text: refusal)
                    .frame(maxWidth: 560)
            } else if model.isGeneral && model.isLongCalculation {
                StatusBanner(kind: .warning, text: model.longRunMessage + " Run asks for confirmation; a running calculation can be cancelled.")
                    .frame(maxWidth: 560)
            }
            ContentUnavailableView(
                "No results yet",
                systemImage: "atom",
                description: Text(model.isGeneral
                                  ? "Choose the job, method and basis set in the sidebar and press Run."
                                  : "Press Run to start a Hartree–Fock and Full CI calculation.")
            )
        }
        .padding()
    }
}
