//
//  MacQCEngine.swift
//  MacQC
//
//  Swift wrapper around the C API of the C++ backend (macqc_run / macqc_free)
//  and Codable models that mirror the JSON format documented in README.md.
//
//  All model types are marked `nonisolated` so they work unchanged whether or not
//  the project uses "Default Actor Isolation = MainActor" (the Xcode 26 template
//  default) and can be decoded on a background thread.
//

import Foundation

// MARK: - Input

/// Input JSON of `macqc_run`.
nonisolated struct MacQCInput: Codable, Sendable, Equatable {
    var R: Double = 1.4632
    var zeta1: Double = 2.0925
    var zeta2: Double = 1.24
    var ZA: Double = 2
    var ZB: Double = 1
    var nGaussians: Int = 3
    var exactErf: Bool = false
    var maxIter: Int = 25
    var convergence: Double = 1.0e-4
    var printLevel: Int = 2
}

// MARK: - Output

/// The 16 unique integrals, named exactly as in the FORTRAN program.
nonisolated struct MacQCIntegrals: Codable, Sendable, Equatable {
    var S12: Double
    var T11: Double
    var T12: Double
    var T22: Double
    var V11A: Double
    var V12A: Double
    var V22A: Double
    var V11B: Double
    var V12B: Double
    var V22B: Double
    var V1111: Double
    var V2111: Double
    var V2121: Double
    var V2211: Double
    var V2221: Double
    var V2222: Double

    /// Name/value pairs in the order used by the printed report.
    var namedValues: [(name: String, value: Double)] {
        [
            ("S12", S12), ("T11", T11), ("T12", T12), ("T22", T22),
            ("V11A", V11A), ("V12A", V12A), ("V22A", V22A),
            ("V11B", V11B), ("V12B", V12B), ("V22B", V22B),
            ("V1111", V1111), ("V2111", V2111), ("V2121", V2121),
            ("V2211", V2211), ("V2221", V2221), ("V2222", V2222),
        ]
    }
}

nonisolated struct MacQCSCFIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var electronicEnergy: Double
    var deltaP: Double

    var id: Int { iteration }
}

nonisolated struct MacQCHF: Codable, Sendable, Equatable {
    var converged: Bool
    var iterations: Int
    var electronicEnergy: Double
    var nuclearRepulsion: Double
    var totalEnergy: Double
    var orbitalEnergies: [Double]
    var C: [[Double]]
    var P: [[Double]]
    var S: [[Double]]
    var X: [[Double]]
    var H: [[Double]]
    var F: [[Double]]
    var G: [[Double]]
    var PS: [[Double]]
    var integrals: MacQCIntegrals
    var scfHistory: [MacQCSCFIteration]
}

nonisolated struct MacQCFCIState: Codable, Sendable, Equatable {
    var energy: Double
    var S2: Double
    var spin: String
    var vector: [Double]
}

nonisolated struct MacQCSingletWeight: Codable, Sendable, Equatable {
    var configuration: String
    var weight: Double
}

nonisolated struct MacQCFCI: Codable, Sendable, Equatable {
    var totalEnergy: Double
    var electronicEnergy: Double
    var correlationEnergy: Double
    var states: [MacQCFCIState]
    var groundStateVector: [Double]
    var determinants: [String]
    var singletWeights: [MacQCSingletWeight]
}

/// Output JSON of `macqc_run`. `hf` / `fci` are null when the backend reports an error
/// before (or instead of) producing them.
nonisolated struct MacQCOutput: Codable, Sendable, Equatable {
    var input: MacQCInput
    var error: String?
    var hf: MacQCHF?
    var fci: MacQCFCI?
    var report: String
}

// MARK: - Engine

nonisolated enum MacQCEngineError: LocalizedError {
    case encodingFailed
    case backendReturnedNull
    case decodingFailed(String)

    var errorDescription: String? {
        switch self {
        case .encodingFailed:
            return "The input could not be encoded as JSON."
        case .backendReturnedNull:
            return "The backend returned no result (out of memory?)."
        case .decodingFailed(let detail):
            return "The backend output could not be decoded: \(detail)"
        }
    }
}

nonisolated enum MacQCEngine {

    /// Synchronous call into the C++ backend. Blocks the calling thread, so do not call it
    /// from the main thread for long scans; use `run(_:)` instead.
    static func runSync(_ input: MacQCInput) throws -> MacQCOutput {
        let encoder = JSONEncoder()
        let inputData = try encoder.encode(input)
        guard let inputJSON = String(data: inputData, encoding: .utf8) else {
            throw MacQCEngineError.encodingFailed
        }

        // macqc_run never throws and returns a malloc'ed, NUL-terminated UTF-8 string.
        guard let raw = macqc_run(inputJSON) else {
            throw MacQCEngineError.backendReturnedNull
        }
        defer { macqc_free(raw) }
        let outputText = String(cString: raw)

        do {
            return try JSONDecoder().decode(MacQCOutput.self, from: Data(outputText.utf8))
        } catch {
            throw MacQCEngineError.decodingFailed(String(describing: error))
        }
    }

    /// Runs the calculation on a background thread.
    static func run(_ input: MacQCInput) async throws -> MacQCOutput {
        let task = Task.detached(priority: .userInitiated) {
            try MacQCEngine.runSync(input)
        }
        return try await task.value
    }
}
