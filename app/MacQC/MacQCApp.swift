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
            // File > Export Results as Text… (H2, LiH, H2O: the output file) / as JSON… /
            // Potential Curve as Text… (replace the unused "New Window" item)
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

                Button("Export Potential Curve as Text…") {
                    model.exportCurve()
                }
                .disabled(model.curvePoints.isEmpty)
            }

            // View > the tabs and their pages (⌘1 … ⌘9)
            CommandGroup(before: .toolbar) {
                Button("Summary") { model.show(.summary) }
                    .keyboardShortcut("1", modifiers: [.command])
                Button("Hartree–Fock") { model.show(.methods, methods: .hartreeFock) }
                    .keyboardShortcut("2", modifiers: [.command])
                Button("Coupled Cluster") { model.show(.methods, methods: .coupledCluster) }
                    .keyboardShortcut("3", modifiers: [.command])
                Button("Full CI") { model.show(.methods, methods: .fullCI) }
                    .keyboardShortcut("4", modifiers: [.command])
                Button("Optimisation") { model.show(.optimization) }
                    .keyboardShortcut("5", modifiers: [.command])
                Button("Frequencies") { model.show(.frequencies, vibration: .frequencies) }
                    .keyboardShortcut("6", modifiers: [.command])
                Button("Thermochemistry") { model.show(.frequencies, vibration: .thermochemistry) }
                    .keyboardShortcut("7", modifiers: [.command])
                Button("Potential Curve") { model.show(.potentialCurve) }
                    .keyboardShortcut("8", modifiers: [.command])
                Button("Output") { model.show(.output) }
                    .keyboardShortcut("9", modifiers: [.command])
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
