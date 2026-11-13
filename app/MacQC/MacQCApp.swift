//
//  MacQCApp.swift
//  MacQC
//
//  App entry point and menu commands.
//

import SwiftUI

@main
struct MacQCApp: App {
    @State private var model = AppModel()

    var body: some Scene {
        WindowGroup {
            ContentView(model: model)
                .frame(minWidth: 1000, minHeight: 640)
        }
        .defaultSize(width: 1180, height: 760)
        .commands {
            // File > Export Results as Text… / as JSON… (replace the unused "New Window" item)
            CommandGroup(replacing: .newItem) {
                Button("Export Results as Text…") {
                    model.exportText()
                }
                .keyboardShortcut("e", modifiers: [.command, .shift])
                .disabled(!model.canExportOrCopy)

                Button("Export Results as JSON…") {
                    model.exportJSON()
                }
                .keyboardShortcut("e", modifiers: [.command, .shift, .option])
                .disabled(!model.canExportOrCopy)
            }

            // View > Summary … Report (⌘1 … ⌘5)
            CommandGroup(before: .toolbar) {
                Button("Summary") { model.selectedTab = .summary }
                    .keyboardShortcut("1", modifiers: [.command])
                Button("Hartree–Fock") { model.selectedTab = .hartreeFock }
                    .keyboardShortcut("2", modifiers: [.command])
                Button("Coupled Cluster") { model.selectedTab = .coupledCluster }
                    .keyboardShortcut("3", modifiers: [.command])
                Button("Full CI") { model.selectedTab = .fullCI }
                    .keyboardShortcut("4", modifiers: [.command])
                Button("Potential Curve") { model.selectedTab = .potentialCurve }
                    .keyboardShortcut("5", modifiers: [.command])
                Button("Report") { model.selectedTab = .report }
                    .keyboardShortcut("6", modifiers: [.command])
                Divider()
            }

            // Edit > Copy Summary as Text
            CommandGroup(after: .pasteboard) {
                Button("Copy Summary as Text") {
                    model.copySummary()
                }
                .keyboardShortcut("c", modifiers: [.command, .shift])
                .disabled(!model.canExportOrCopy)
            }
        }
    }
}
