//
//  MethodsView.swift
//  MacQC
//
//  Methods tab: the detailed results of Hartree–Fock, coupled cluster and Full CI, one page each
//  (the tabs of versions 3.0 and 4.0, combined so that the window keeps six tabs).
//

import SwiftUI

struct MethodsView: View {
    @Bindable var model: AppModel

    var body: some View {
        VStack(spacing: 0) {
            Picker("Method", selection: $model.methodsPage) {
                ForEach(MethodsPage.allCases) { page in
                    Text(page.rawValue).tag(page)
                }
            }
            .pickerStyle(.segmented)
            .labelsHidden()
            .frame(maxWidth: 420)
            .padding(.top, 10)
            .padding(.bottom, 4)
            .accessibilityLabel("Method page")

            switch model.methodsPage {
            case .hartreeFock:
                HartreeFockView(model: model)
            case .coupledCluster:
                CoupledClusterView(model: model)
            case .fullCI:
                FullCIView(model: model)
            }
        }
    }
}
