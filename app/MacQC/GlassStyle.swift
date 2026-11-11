//
//  GlassStyle.swift
//  MacQC
//
//  Liquid Glass on macOS 26 and later, with material / standard-style fallbacks on
//  macOS 14-15. Glass is used for controls and panels, never for tables of numbers.
//
//  NOTE: this file needs the macOS 26 SDK (Xcode 26) to compile, because it references
//  glassEffect, GlassEffectContainer, .glass and .glassProminent. The `#available`
//  checks only protect the *runtime*.
//

import SwiftUI
import AppKit

extension View {
    /// A glass panel on macOS 26+, a regular-material rounded rectangle before that.
    @ViewBuilder
    func glassPanel(cornerRadius: CGFloat = 16) -> some View {
        if #available(macOS 26, *) {
            self.glassEffect(.regular, in: .rect(cornerRadius: cornerRadius))
        } else {
            self.background(.regularMaterial, in: RoundedRectangle(cornerRadius: cornerRadius))
        }
    }

    /// Prominent action button (the "Run" button).
    @ViewBuilder
    func prominentActionStyle() -> some View {
        if #available(macOS 26, *) {
            self.buttonStyle(.glassProminent)
        } else {
            self.buttonStyle(.borderedProminent)
        }
    }

    /// Secondary button.
    @ViewBuilder
    func secondaryActionStyle() -> some View {
        if #available(macOS 26, *) {
            self.buttonStyle(.glass)
        } else {
            self.buttonStyle(.bordered)
        }
    }

    /// Plain, non-glass container for tables of numbers (readable in light and dark mode).
    func dataPanel(cornerRadius: CGFloat = 12) -> some View {
        self
            .padding(14)
            .background(Color(nsColor: .controlBackgroundColor),
                        in: RoundedRectangle(cornerRadius: cornerRadius))
            .overlay(
                RoundedRectangle(cornerRadius: cornerRadius)
                    .strokeBorder(Color.secondary.opacity(0.18), lineWidth: 1)
            )
    }
}

/// Groups glass views so they blend / morph on macOS 26; a plain pass-through earlier.
struct GlassGroup<Content: View>: View {
    var spacing: CGFloat = 12
    @ViewBuilder var content: () -> Content

    var body: some View {
        if #available(macOS 26, *) {
            GlassEffectContainer(spacing: spacing) {
                content()
            }
        } else {
            content()
        }
    }
}
