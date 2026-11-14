//
//  Molecules.swift
//  MacQC
//
//  The molecules, job types, methods and basis sets offered in the sidebar (ΨLab 5.0).
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

/// Methods for H₂, LiH and H₂O, in the order of the method menu. The raw value is the
/// backend's method name; all correlated methods start from RHF and correlate all electrons.
nonisolated enum MethodChoice: String, CaseIterable, Identifiable, Sendable {
    case hf = "RHF"
    case ccsd = "CCSD"
    case ccsdT = "CCSD(T)"
    case ccsdt = "CCSDT"
    case ccsdtq = "CCSDTQ"
    case fci = "FCI"

    var id: String { rawValue }

    /// Name in the user interface.
    var title: String {
        switch self {
        case .hf: return "HF"
        case .ccsd: return "CCSD"
        case .ccsdT: return "CCSD(T)"
        case .ccsdt: return "CCSDT"
        case .ccsdtq: return "CCSDTQ"
        case .fci: return "Full CI"
        }
    }

    var isCoupledCluster: Bool {
        switch self {
        case .ccsd, .ccsdT, .ccsdt, .ccsdtq: return true
        case .hf, .fci: return false
        }
    }

    /// One-line description under the method menu.
    var explanation: String {
        switch self {
        case .hf: return "Restricted Hartree–Fock (no electron correlation)."
        case .ccsd: return "Coupled cluster with single and double excitations."
        case .ccsdT: return "CCSD plus the perturbative triples correction (T)."
        case .ccsdt: return "Coupled cluster with single, double and triple excitations."
        case .ccsdtq: return "Coupled cluster up to quadruple excitations."
        case .fci: return "Full configuration interaction: exact in the basis set."
        }
    }

    /// Accepts the raw value or the title (launch arguments), in any case.
    init?(name: String) {
        let n = name.lowercased()
        guard let m = MethodChoice.allCases.first(where: { $0.rawValue.lowercased() == n || $0.title.lowercased() == n })
        else { return nil }
        self = m
    }

    /// Title of a backend method name ("RHF" -> "HF").
    static func title(of backendName: String) -> String {
        MethodChoice(rawValue: backendName)?.title ?? backendName
    }
}

/// Job types of ΨLab 5.0 (H₂, LiH, H₂O), in the order of the job menu. The raw value is the
/// backend's job name.
nonisolated enum JobChoice: String, CaseIterable, Identifiable, Sendable {
    case energy = "energy"
    case opt = "opt"
    case freq = "freq"
    case optFreq = "opt+freq"
    case optFreqSP = "opt+freq+sp"

    var id: String { rawValue }

    /// Full name (job menu, output file).
    var title: String {
        switch self {
        case .energy: return "Energy"
        case .opt: return "Geometry optimisation"
        case .freq: return "Frequencies"
        case .optFreq: return "Optimisation + frequencies"
        case .optFreqSP: return "Optimisation + frequencies + single-point energy"
        }
    }

    /// Short name for the closed menu in the sidebar and for the result headers.
    var shortTitle: String {
        switch self {
        case .energy: return "Energy"
        case .opt: return "Optimisation"
        case .freq: return "Frequencies"
        case .optFreq: return "Opt + freq"
        case .optFreqSP: return "Opt + freq + SP"
        }
    }

    var number: Int { (JobChoice.allCases.firstIndex(of: self) ?? 0) + 1 }

    var hasOptimization: Bool { self == .opt || self == .optFreq || self == .optFreqSP }
    var hasFrequencies: Bool { self == .freq || self == .optFreq || self == .optFreqSP }
    var hasSinglePoint: Bool { self == .optFreqSP }

    /// One-line description under the job menu.
    var explanation: String {
        switch self {
        case .energy:
            return "Single-point energy at the input geometry."
        case .opt:
            return "Minimum-energy geometry from the input geometry (redundant internal coordinates, RFO steps, BFGS)."
        case .freq:
            return "Harmonic frequencies, IR intensities and thermochemistry at the input geometry (a warning if it is not a stationary point)."
        case .optFreq:
            return "Optimisation, then frequencies and thermochemistry at the optimised geometry."
        case .optFreqSP:
            return "Optimisation and frequencies, then a single-point energy with its own method and basis set; the final energies add the thermal corrections."
        }
    }

    /// Accepts the raw value, the title or the number (launch arguments), in any case.
    init?(name: String) {
        let n = name.lowercased()
        guard let j = JobChoice.allCases.first(where: {
            $0.rawValue == n || $0.title.lowercased() == n || "\($0.number)" == n
        }) else { return nil }
        self = j
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

    /// Memory limit of Full CI and coupled cluster: 60 % of the RAM (at least 2 GB, at most
    /// 48 GB).
    static var fciMemoryMB: Double {
        min(max(physicalMemoryMB * 0.6, 2048), 49152)
    }
}
