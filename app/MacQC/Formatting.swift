//
//  Formatting.swift
//  MacQC
//
//  Number formatting helpers (fixed decimals, "E_h" units).
//

import Foundation

nonisolated enum Fmt {
    /// Energy with 8 decimals, e.g. "-2.86066199". Unit not included.
    static func energy(_ x: Double) -> String {
        String(format: "%.8f", x)
    }

    /// Energy with 8 decimals and unit, e.g. "-2.86066199 E_h".
    static func energyWithUnit(_ x: Double) -> String {
        String(format: "%.8f E_h", x)
    }

    /// Matrix / integral entries with 6 decimals (as in the text report).
    static func value6(_ x: Double) -> String {
        String(format: "%.6f", x)
    }

    static func value8(_ x: Double) -> String {
        String(format: "%.8f", x)
    }

    static func scientific(_ x: Double) -> String {
        String(format: "%.3e", x)
    }

    static func bohr(_ x: Double) -> String {
        String(format: "%.4f bohr", x)
    }

    /// Length given in bohr, shown as "0.7743 Å (1.4632 bohr)".
    static func lengthBoth(_ bohr: Double, decimals: Int = 4) -> String {
        let digits = Int32(decimals)  // "%.*f" takes a C int
        return String(format: "%.*f Å (%.*f bohr)", digits, bohr * Units.bohrToAngstrom, digits, bohr)
    }
}

nonisolated enum Units {
    /// 1 bohr in angstrom (CODATA 2022).
    static let bohrToAngstrom = 0.529177210544
}
