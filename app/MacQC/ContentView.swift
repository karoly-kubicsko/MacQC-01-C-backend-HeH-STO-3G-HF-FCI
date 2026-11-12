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
                    if model.refusal == nil && !model.isLongCalculation {
                        await model.run()
                    }
                } else {
                    await model.run()
                }
            }
        }
    }
}

struct DetailView: View {
    @Bindable var model: AppModel

    var body: some View {
        VStack(spacing: 0) {
            if model.isGeneral && !model.isRunning, let refusal = model.refusal {
                StatusBanner(kind: .refusal, text: refusal)
                    .padding([.horizontal, .top], 14)
            } else if model.isGeneral && !model.isRunning && model.isLongCalculation && !model.resultMatchesInputs,
                      let warning = model.estimate?.warnings?.first {
                StatusBanner(kind: .warning, text: warning + " Run asks for confirmation; a running calculation can be cancelled.")
                    .padding([.horizontal, .top], 14)
            }
            tabs
        }
    }

    private var tabs: some View {
        TabView(selection: $model.selectedTab) {
            SummaryView(model: model)
                .tabItem { Label("Summary", systemImage: "list.bullet.rectangle") }
                .tag(DetailTab.summary)

            HartreeFockView(model: model)
                .tabItem { Label("Hartree–Fock", systemImage: "square.grid.2x2") }
                .tag(DetailTab.hartreeFock)

            FullCIView(model: model)
                .tabItem { Label("Full CI", systemImage: "function") }
                .tag(DetailTab.fullCI)

            PotentialCurveView(model: model)
                .tabItem { Label("Potential curve", systemImage: "chart.xyaxis.line") }
                .tag(DetailTab.potentialCurve)

            ReportView(model: model)
                .tabItem { Label("Report", systemImage: "doc.plaintext") }
                .tag(DetailTab.report)
        }
        .padding(.top, 4)
    }
}

// MARK: - Shared small views

/// Banner above the tabs: why the current settings cannot be calculated (red), or that the
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
                    .fixedSize(horizontal: false, vertical: true)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
        } icon: {
            Image(systemName: kind == .refusal ? "xmark.octagon.fill" : "clock.badge.exclamationmark")
                .foregroundStyle(color)
        }
        .padding(12)
        .background(color.opacity(0.12), in: RoundedRectangle(cornerRadius: 10))
        .overlay(RoundedRectangle(cornerRadius: 10).strokeBorder(color.opacity(0.35), lineWidth: 1))
        .accessibilityElement(children: .combine)
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
            ContentUnavailableView(
                "No results yet",
                systemImage: "atom",
                description: Text(model.isGeneral
                                  ? "Choose the basis set and method in the sidebar and press Run."
                                  : "Press Run to start a Hartree–Fock and Full CI calculation.")
            )
        }
        .padding()
    }
}
