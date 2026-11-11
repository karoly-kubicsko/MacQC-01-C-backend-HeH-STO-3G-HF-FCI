//
//  AppearanceChoice.swift
//  MacQC
//
//  In-app appearance override: follow the system (default), or force light / dark.
//  Stored in UserDefaults under the key `appearanceChoice` (so it can also be set on the
//  command line: open ΨLab.app --args -appearanceChoice dark).
//

import SwiftUI

enum AppearanceChoice: String, CaseIterable, Identifiable {
    case system
    case light
    case dark

    static let storageKey = "appearanceChoice"

    var id: String { rawValue }

    var title: String {
        switch self {
        case .system: return "System"
        case .light: return "Light"
        case .dark: return "Dark"
        }
    }

    /// `nil` means "follow the system appearance".
    var colorScheme: ColorScheme? {
        switch self {
        case .system: return nil
        case .light: return .light
        case .dark: return .dark
        }
    }
}
