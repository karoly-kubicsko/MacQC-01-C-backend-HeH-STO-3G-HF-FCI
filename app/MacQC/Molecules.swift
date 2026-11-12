//
//  Molecules.swift
//  MacQC
//
//  The molecules, methods and basis sets offered in the sidebar (ΨLab 3.0).
//

import Foundation

/// HeH⁺ is the Szabo–Ostlund STO-nG calculation of versions 1 and 2; H₂, LiH and H₂O use the
/// correlation-consistent Gaussian basis sets.
nonisolated enum MoleculeChoice: String, CaseIterable, Identifiable, Sendable {
    case heh = "HeH+"
    case h2 = "H2"
    case lih = "LiH"
    case h2o = "H2O"

    var id: String { rawValue }

    /// Name in the user interface.
    var title: String {
        switch self {
        case .heh: return "HeH⁺"
        case .h2: return "H₂"
        case .lih: return "LiH"
        case .h2o: return "H₂O"
        }
    }

    var electrons: Int {
        switch self {
        case .heh, .h2: return 2
        case .lih: return 4
        case .h2o: return 10
        }
    }

    /// Label of the bond length, e.g. "R(O–H)".
    var bondLabel: String {
        switch self {
        case .heh: return "R"
        case .h2: return "R(H–H)"
        case .lih: return "R(Li–H)"
        case .h2o: return "R(O–H)"
        }
    }

    /// Experimental equilibrium bond length in bohr (H₂ 0.74144 Å, LiH 1.5949 Å, H₂O 0.9572 Å);
    /// HeH⁺: the Szabo–Ostlund value 1.4632 bohr.
    var defaultBondLength: Double {
        switch self {
        case .heh: return 1.4632
        case .h2: return 0.74144 / Units.bohrToAngstrom
        case .lih: return 1.5949 / Units.bohrToAngstrom
        case .h2o: return 0.9572 / Units.bohrToAngstrom
        }
    }

    static let defaultWaterAngle = 104.52

    /// Potential-curve range in bohr (H₂O: symmetric stretch at the chosen angle).
    var scanRange: ClosedRange<Double> {
        switch self {
        case .heh: return 0.8...5.0
        case .h2: return (0.40 / Units.bohrToAngstrom)...(3.00 / Units.bohrToAngstrom)
        case .lih: return (1.00 / Units.bohrToAngstrom)...(4.00 / Units.bohrToAngstrom)
        case .h2o: return (0.70 / Units.bohrToAngstrom)...(2.00 / Units.bohrToAngstrom)
        }
    }

    /// Short note shown under the molecule picker.
    var subtitle: String {
        switch self {
        case .heh: return "Szabo–Ostlund, STO-nG, HF and Full CI"
        case .h2: return "2 electrons, closed shell"
        case .lih: return "4 electrons, closed shell"
        case .h2o: return "10 electrons, closed shell (Full CI impossible)"
        }
    }
}

nonisolated enum MethodChoice: String, CaseIterable, Identifiable, Sendable {
    case rhf = "RHF"
    case fci = "FCI"

    var id: String { rawValue }

    var title: String {
        switch self {
        case .rhf: return "RHF"
        case .fci: return "Full CI"
        }
    }
}

nonisolated enum BasisSets {
    /// The Gaussian94 files in the app bundle's basis folder (from basissetexchange.org).
    static let all = [
        "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z",
        "aug-cc-pCVDZ", "aug-cc-pCVTZ", "aug-cc-pCVQZ", "aug-cc-pCV5Z",
    ]

    /// Folder with the .gbs files inside the app bundle ("" if missing; the backend then reports
    /// that it cannot read the file).
    static var directory: String {
        Bundle.main.resourceURL?.appendingPathComponent("basis", isDirectory: true).path ?? ""
    }
}

nonisolated enum MachineLimits {
    static var physicalMemoryMB: Double {
        Double(ProcessInfo.processInfo.physicalMemory) / 1_048_576
    }

    /// Two-electron integrals are kept in memory up to a quarter of the RAM (at most 8 GB).
    static var integralMemoryMB: Double {
        min(max(physicalMemoryMB / 4, 512), 8192)
    }

    /// Memory limit of a Full CI: half of the RAM (at least 2 GB, at most 16 GB).
    static var fciMemoryMB: Double {
        min(max(physicalMemoryMB / 2, 2048), 16384)
    }
}
