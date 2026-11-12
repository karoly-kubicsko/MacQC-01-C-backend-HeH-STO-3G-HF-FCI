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

    /// Energy with 10 decimals (general molecules: agrees with the backend's report).
    static func energy10(_ x: Double) -> String {
        String(format: "%.10f", x)
    }

    /// Integer count with thousands separators ("11,029,041"), or "5.62e+11" when very large.
    static func count(_ x: Double) -> String {
        if x >= 1.0e10 { return String(format: "%.2e", x) }
        let formatter = NumberFormatter()
        formatter.numberStyle = .decimal
        formatter.groupingSeparator = ","
        formatter.usesGroupingSeparator = true
        formatter.maximumFractionDigits = 0
        return formatter.string(from: NSNumber(value: x.rounded())) ?? String(format: "%.0f", x)
    }

    /// Seconds as "less than a second", "about 25 s", "about 4 min", "about 1.5 h".
    static func duration(_ seconds: Double) -> String {
        if seconds < 1 { return "less than a second" }
        if seconds < 90 { return String(format: "about %.0f s", seconds) }
        if seconds < 5400 { return String(format: "about %.0f min", seconds / 60) }
        return String(format: "about %.1f h", seconds / 3600)
    }

    /// Memory in MB as "43 MB" or "2.8 GB".
    static func memory(_ megabytes: Double) -> String {
        if megabytes < 1024 { return String(format: "%.0f MB", max(megabytes, 1)) }
        return String(format: "%.1f GB", megabytes / 1024)
    }

    /// One entry per element (atoms of the same element share the basis set).
    static func uniqueAtoms(_ atoms: [MacQCAtomBasis]) -> [MacQCAtomBasis] {
        var seen = Set<String>()
        return atoms.filter { seen.insert($0.symbol).inserted }
    }

    /// Chemical formula with subscripts: "H2O" -> "H₂O".
    static func formula(_ name: String) -> String {
        let subscripts: [Character: Character] = ["0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄",
                                                  "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉"]
        return String(name.map { subscripts[$0] ?? $0 })
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
