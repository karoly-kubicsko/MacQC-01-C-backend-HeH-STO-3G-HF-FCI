//
//  ContentView.swift
//  MacQC
//
//  NavigationSplitView: input sidebar + tabbed detail area.
//

import SwiftUI

struct ContentView: View {
    @Bindable var model: AppModel

    var body: some View {
        NavigationSplitView {
            SidebarView(model: model)
                .navigationSplitViewColumnWidth(min: 260, ideal: 290, max: 360)
        } detail: {
            DetailView(model: model)
        }
        .navigationTitle("MacQC")
        .task {
            // Show a result right away with the Szabo–Ostlund defaults.
            if model.result == nil {
                await model.run()
            }
        }
    }
}

struct DetailView: View {
    @Bindable var model: AppModel

    var body: some View {
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
        }
        .padding(.top, 4)
    }
}

// MARK: - Shared small views

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
                description: Text("Press Run to start a Hartree–Fock and Full CI calculation.")
            )
        }
        .padding()
    }
}
