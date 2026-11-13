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

/// Input JSON of `macqc_run`. `molecule` "HeH+" selects the Szabo–Ostlund STO-nG calculation
/// (keys R … printLevel); "H2", "LiH", "H2O" the general molecules (R, basis, angle, method, …).
/// The backend echoes only the keys of the selected path, so every key is optional when decoding.
nonisolated struct MacQCInput: Codable, Sendable, Equatable {
    var molecule: String = "HeH+"
    var R: Double = 1.4632          // bohr (H2O: O–H)
    var zeta1: Double = 2.0925
    var zeta2: Double = 1.24
    var ZA: Double = 2
    var ZB: Double = 1
    var nGaussians: Int = 3
    var exactErf: Bool = false
    var maxIter: Int = 25
    var convergence: Double = 1.0e-4
    var printLevel: Int = 2
    // General molecules
    var basis: String = "aug-cc-pVDZ"
    var angle: Double = 104.52      // degrees, H2O
    var method: String = "RHF"      // "RHF", "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ" or "FCI"
    var basisDirectory: String = ""
    var threads: Int = 0
    var memoryMB: Double = 1024
    var fciMemoryMB: Double = 4096
    var estimateOnly: Bool = false
    // ΨLab 5.0 job types (general molecules)
    var job: String = "energy"      // "energy", "opt", "freq", "opt+freq", "opt+freq+sp"
    var spMethod: String = ""       // opt+freq+sp ("" = method)
    var spBasis: String = ""        // opt+freq+sp ("" = basis)
    var temperature: Double = 298.15
    var pressure: Double = 101325
    var optConvergence: String = "normal"
    var maxSteps: Int = 50

    init() {}

    private enum CodingKeys: String, CodingKey {
        case molecule, R, zeta1, zeta2, ZA, ZB, nGaussians, exactErf, maxIter, convergence, printLevel
        case basis, angle, method, basisDirectory, threads, memoryMB, fciMemoryMB, estimateOnly
        case job, spMethod, spBasis, temperature, pressure, optConvergence, maxSteps
    }

    init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        let d = MacQCInput()
        molecule = try c.decodeIfPresent(String.self, forKey: .molecule) ?? d.molecule
        R = try c.decodeIfPresent(Double.self, forKey: .R) ?? d.R
        zeta1 = try c.decodeIfPresent(Double.self, forKey: .zeta1) ?? d.zeta1
        zeta2 = try c.decodeIfPresent(Double.self, forKey: .zeta2) ?? d.zeta2
        ZA = try c.decodeIfPresent(Double.self, forKey: .ZA) ?? d.ZA
        ZB = try c.decodeIfPresent(Double.self, forKey: .ZB) ?? d.ZB
        nGaussians = try c.decodeIfPresent(Int.self, forKey: .nGaussians) ?? d.nGaussians
        exactErf = try c.decodeIfPresent(Bool.self, forKey: .exactErf) ?? d.exactErf
        maxIter = try c.decodeIfPresent(Int.self, forKey: .maxIter) ?? d.maxIter
        convergence = try c.decodeIfPresent(Double.self, forKey: .convergence) ?? d.convergence
        printLevel = try c.decodeIfPresent(Int.self, forKey: .printLevel) ?? d.printLevel
        basis = try c.decodeIfPresent(String.self, forKey: .basis) ?? d.basis
        angle = try c.decodeIfPresent(Double.self, forKey: .angle) ?? d.angle
        method = try c.decodeIfPresent(String.self, forKey: .method) ?? d.method
        basisDirectory = try c.decodeIfPresent(String.self, forKey: .basisDirectory) ?? d.basisDirectory
        threads = try c.decodeIfPresent(Int.self, forKey: .threads) ?? d.threads
        memoryMB = try c.decodeIfPresent(Double.self, forKey: .memoryMB) ?? d.memoryMB
        fciMemoryMB = try c.decodeIfPresent(Double.self, forKey: .fciMemoryMB) ?? d.fciMemoryMB
        estimateOnly = try c.decodeIfPresent(Bool.self, forKey: .estimateOnly) ?? d.estimateOnly
        job = try c.decodeIfPresent(String.self, forKey: .job) ?? d.job
        spMethod = try c.decodeIfPresent(String.self, forKey: .spMethod) ?? d.spMethod
        spBasis = try c.decodeIfPresent(String.self, forKey: .spBasis) ?? d.spBasis
        temperature = try c.decodeIfPresent(Double.self, forKey: .temperature) ?? d.temperature
        pressure = try c.decodeIfPresent(Double.self, forKey: .pressure) ?? d.pressure
        optConvergence = try c.decodeIfPresent(String.self, forKey: .optConvergence) ?? d.optConvergence
        maxSteps = try c.decodeIfPresent(Int.self, forKey: .maxSteps) ?? d.maxSteps
    }

    /// True for H2, LiH and H2O (the general Gaussian-basis path).
    var isGeneral: Bool { molecule != "HeH+" }
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

// MARK: - Output of the general molecules (H2, LiH, H2O)

nonisolated struct MacQCAtom: Codable, Sendable, Equatable {
    var symbol: String
    var Z: Int
    var position: [Double]   // bohr
}

nonisolated struct MacQCMoleculeInfo: Codable, Sendable, Equatable {
    var name: String
    var electrons: Int
    var charge: Int
    var R: Double            // bohr
    var angle: Double?       // degrees (H2O)
    var nuclearRepulsion: Double
    var atoms: [MacQCAtom]
}

nonisolated struct MacQCAtomBasis: Codable, Sendable, Equatable {
    var symbol: String
    var basis: String
    var contraction: String
    var primitives: String
    var functions: Int
}

nonisolated struct MacQCBasisInfo: Codable, Sendable, Equatable {
    var name: String
    var functions: Int
    var shells: Int
    var maxL: Int
    var atoms: [MacQCAtomBasis]
}

nonisolated struct MacQCFCIEstimate: Codable, Sendable, Equatable {
    var orbitals: Int
    var alphaElectrons: Int
    var betaElectrons: Int
    var alphaStrings: Double
    var betaStrings: Double
    var determinants: Double
    var memoryMB: Double
    var seconds: Double
    var feasible: Bool
    var reason: String?
}

nonisolated struct MacQCEstimate: Codable, Sendable, Equatable {
    var basisFunctions: Int
    var uniqueIntegrals: Double
    var integralMemoryMB: Double
    var integralsStored: Bool
    var threads: Int
    var rhfSeconds: Double
    var totalSeconds: Double
    var longCalculation: Bool
    var fci: MacQCFCIEstimate?
    var coupledCluster: MacQCCCEstimate?
}

nonisolated struct MacQCCCEstimate: Codable, Sendable, Equatable {
    var method: String
    var engine: String
    var occupied: Int
    var virtuals: Int
    var amplitudes: Double
    var memoryMB: Double
    var ladderIntegralsStored: Bool
    var iterations: Int
    var iterationSeconds: Double
    var triplesSeconds: Double
    var seconds: Double
    var feasible: Bool
    var reason: String?
}

nonisolated struct MacQCRHFIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var totalEnergy: Double
    var deltaE: Double
    var gradient: Double
    var diis: Bool

    var id: Int { iteration }
}

nonisolated struct MacQCRHF: Codable, Sendable, Equatable {
    var converged: Bool
    var iterations: Int
    var totalEnergy: Double
    var electronicEnergy: Double
    var nuclearRepulsion: Double
    var kineticEnergy: Double
    var nuclearAttractionEnergy: Double
    var electronRepulsionEnergy: Double
    var virialRatio: Double
    var orbitalEnergies: [Double]
    var occupiedOrbitals: Int
    var homo: Double?
    var lumo: Double?
    var smallestOverlapEigenvalue: Double
    var removedFunctions: Int
    var integralsStored: Bool
    var uniqueIntegrals: Double
    var threads: Int
    var integralSeconds: Double
    var scfSeconds: Double
    var scfHistory: [MacQCRHFIteration]
}

nonisolated struct MacQCDeterminant: Codable, Sendable, Equatable {
    var coefficient: Double
    var alpha: [Int]
    var beta: [Int]
    var excitation: Int
}

nonisolated struct MacQCDavidsonIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var totalEnergy: Double
    var residual: Double
    var subspace: Int

    var id: Int { iteration }
}

nonisolated struct MacQCFullCI: Codable, Sendable, Equatable {
    var converged: Bool
    var iterations: Int
    var totalEnergy: Double
    var electronicEnergy: Double
    var correlationEnergy: Double
    var S2: Double
    var hfWeight: Double
    var orbitals: Int
    var alphaElectrons: Int
    var betaElectrons: Int
    var alphaStrings: Double
    var betaStrings: Double
    var determinants: Double
    var leadingDeterminants: [MacQCDeterminant]
    var naturalOccupations: [Double]
    var history: [MacQCDavidsonIteration]
    var transformSeconds: Double
    var ciSeconds: Double
}

nonisolated struct MacQCCCAmplitude: Codable, Sendable, Equatable {
    var value: Double
    var rank: Int
    var occupied: [Int]
    var virtuals: [Int]
    var spins: String
}

nonisolated struct MacQCCCIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var totalEnergy: Double
    var correlationEnergy: Double
    var deltaE: Double
    var residual: Double

    var id: Int { iteration }
}

/// Coupled-cluster result (CCSD, CCSD(T), CCSDT, CCSDTQ).
nonisolated struct MacQCCoupledCluster: Codable, Sendable, Equatable {
    var method: String
    var engine: String
    var converged: Bool
    var iterations: Int
    var occupied: Int
    var virtuals: Int
    var referenceEnergy: Double
    var correlationEnergy: Double      // without (T)
    var ccEnergy: Double               // reference + correlation
    var triplesCorrection: Double?     // CCSD(T) only
    var totalEnergy: Double
    var t1Diagnostic: Double
    var amplitudeCounts: [Double]      // index = rank
    var amplitudeNorms: [Double]
    var leadingAmplitudes: [MacQCCCAmplitude]
    var history: [MacQCCCIteration]
    var ladder: String?
    var transformSeconds: Double
    var ccSeconds: Double
    var triplesSeconds: Double?

    /// Correlation energy including (T).
    var totalCorrelation: Double { totalEnergy - referenceEnergy }
}

// MARK: - Jobs of ΨLab 5.0 (optimisation, frequencies, thermochemistry)

/// A geometry: atoms (bohr) and the parameters of the inputs (R in bohr; H₂O: both O–H
/// distances, their mean and the H–O–H angle in degrees).
nonisolated struct MacQCGeometry: Codable, Sendable, Equatable {
    var atoms: [MacQCAtom]
    var R: Double?
    var R1: Double?
    var R2: Double?
    var angle: Double?
}

nonisolated struct MacQCJobStep: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var method: String
    var basis: String
    var how: String
    var analytic: Bool
    var evaluations: Int
    var evaluationKind: String      // "steps", "gradients", "energies", "Hessians"
    var seconds: Double
    var memoryMB: Double

    var id: String { name }
}

nonisolated struct MacQCJobEstimate: Codable, Sendable, Equatable {
    var steps: [MacQCJobStep]
    var totalSeconds: Double
    var memoryMB: Double
    var longCalculation: Bool
    var refusal: String?
    var warnings: [String]
    var notes: [String]
}

nonisolated struct MacQCJobInfo: Codable, Sendable, Equatable {
    var type: String
    var title: String
    var termination: String          // "normal", "error", "refused", "cancelled", "estimate"
    var completed: Bool
    var message: String?
    var estimate: MacQCJobEstimate
    var warnings: [String]
    var seconds: Double
    var initialGeometry: MacQCGeometry?
    var finalGeometry: MacQCGeometry?
    var finalEnergy: Double?
    var finalEnergyLabel: String?
    var singlePointMethod: String?
    var singlePointBasis: String?
}

nonisolated struct MacQCParameter: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var value: Double                // bohr or degrees
    var angle: Bool

    var id: String { name }
}

nonisolated struct MacQCInternalCoordinate: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var definition: String
    var value: Double                // bohr or radians
    var gradient: Double             // E_h/bohr or E_h/rad
    var step: Double?

    var id: String { name }
}

nonisolated struct MacQCConvergenceItem: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var value: Double?
    var threshold: Double
    var met: Bool
    var applicable: Bool

    var id: String { name }
}

nonisolated struct MacQCOptimizationStep: Codable, Sendable, Equatable, Identifiable {
    var step: Int
    var energy: Double
    var deltaE: Double?
    var geometry: MacQCGeometry
    var coordinates: [MacQCInternalCoordinate]
    var convergence: [MacQCConvergenceItem]
    var converged: Bool
    var hessianUpdate: String
    var trustRadius: Double
    var stepNorm: Double
    var stepScaled: Bool
    var predictedChange: Double
    var ratio: Double?
    var warnings: [String]
    var seconds: Double

    var id: Int { step }

    /// Value of a convergence item by position: 0 max force, 1 RMS force, 2 max step, 3 RMS step.
    func criterion(_ index: Int) -> Double? {
        index < convergence.count ? convergence[index].value : nil
    }
}

nonisolated struct MacQCCriteria: Codable, Sendable, Equatable {
    var maxForce: Double
    var rmsForce: Double
    var maxDisplacement: Double
    var rmsDisplacement: Double
    var energy: Double
}

nonisolated struct MacQCOptimization: Codable, Sendable, Equatable {
    var converged: Bool
    var steps: Int
    var message: String
    var initialHessian: String
    var gradient: String?
    var criteria: MacQCCriteria
    var initialEnergy: Double
    var finalEnergy: Double
    var initialGeometry: MacQCGeometry
    var finalGeometry: MacQCGeometry
    var initialParameters: [MacQCParameter]
    var finalParameters: [MacQCParameter]
    var symmetricStart: Bool
    var symmetricEnd: Bool
    var symmetryDeviation: Double
    var history: [MacQCOptimizationStep]
    var seconds: Double
}

nonisolated struct MacQCNormalMode: Codable, Sendable, Equatable {
    var wavenumber: Double           // cm⁻¹, negative = imaginary
    var imaginary: Bool
    var symmetry: String
    var reducedMass: Double          // u
    var forceConstant: Double        // mdyn/Å
    var irIntensity: Double?         // km/mol
    var displacement: [[Double]]     // normalised Cartesian displacement per atom
}

nonisolated struct MacQCFrequencies: Codable, Sendable, Equatable {
    var analytic: Bool
    var description: String
    var numericalReason: String?
    var gradientEvaluations: Int
    var energyEvaluations: Int
    var coordinates: [String]
    var intensities: Bool
    var intensityDescription: String
    var dipole: [Double]?
    var linear: Bool
    var vibrations: Int
    var masses: [Double]
    var isotopes: [String]
    var geometry: MacQCGeometry
    var energy: Double
    var lowFrequencies: [Double]
    var modes: [MacQCNormalMode]
    var seconds: Double
}

nonisolated struct MacQCThermoTerm: Codable, Sendable, Equatable {
    var energy: Double               // J/mol
    var entropy: Double              // J/(mol K)
    var cv: Double                   // J/(mol K)
    var q: Double
}

nonisolated struct MacQCThermoTerms: Codable, Sendable, Equatable {
    var translation: MacQCThermoTerm
    var rotation: MacQCThermoTerm
    var vibration: MacQCThermoTerm
    var electronic: MacQCThermoTerm
    var total: MacQCThermoTerm
}

/// Thermal corrections (J/mol) or totals (E_h).
nonisolated struct MacQCThermoValues: Codable, Sendable, Equatable {
    var zeroPoint: Double
    var thermalEnergy: Double
    var enthalpy: Double
    var gibbs: Double
    var helmholtz: Double
}

nonisolated struct MacQCThermochemistry: Codable, Sendable, Equatable {
    var temperature: Double
    var pressure: Double
    var symmetryNumber: Int
    var linear: Bool
    var mass: Double
    var moments: [Double]                    // u bohr²
    var rotationalConstantsGHz: [Double]
    var rotationalConstantsCm: [Double]
    var rotationalTemperatures: [Double]
    var vibrationalTemperatures: [Double]
    var excludedModes: Int
    var warnings: [String]
    var terms: MacQCThermoTerms
    var zeroPointEnergy: Double              // J/mol
    var corrections: MacQCThermoValues       // J/mol
    var electronicEnergy: Double             // E_h
    var totals: MacQCThermoValues            // E_h
}

nonisolated struct MacQCSinglePoint: Codable, Sendable, Equatable {
    var method: String
    var energy: Double
    var converged: Bool
    var rhfEnergy: Double
    var correlationEnergy: Double
    var geometry: MacQCGeometry
    var warnings: [String]
    var seconds: Double
}

/// Output JSON of `macqc_run`. HeH+: `hf` / `fci` (null when the backend reports an error
/// before producing them). H2, LiH, H2O: `molecule`, `basis`, `estimate`, `notes`,
/// `warnings`, `refusal`, `rhf`, `fullCI`, `coupledCluster`.
nonisolated struct MacQCOutput: Codable, Sendable, Equatable {
    var input: MacQCInput
    var error: String?
    var hf: MacQCHF?
    var fci: MacQCFCI?
    var report: String
    var molecule: MacQCMoleculeInfo?
    var basis: MacQCBasisInfo?
    var estimate: MacQCEstimate?
    var notes: [String]?
    var warnings: [String]?
    var refusal: String?
    var rhf: MacQCRHF?
    var fullCI: MacQCFullCI?
    var coupledCluster: MacQCCoupledCluster?
    // ΨLab 5.0 (H2, LiH, H2O): the job and its output file; rhf / fullCI / coupledCluster /
    // molecule / basis describe the final energy calculation of the job.
    var job: MacQCJobInfo?
    var optimization: MacQCOptimization?
    var frequencies: MacQCFrequencies?
    var thermochemistry: MacQCThermochemistry?
    var singlePoint: MacQCSinglePoint?
    var output: String?

    /// Energy of the correlated method of this result (Full CI or coupled cluster).
    var correlatedEnergy: Double? { coupledCluster?.totalEnergy ?? fullCI?.totalEnergy }

    /// Title of the method of this result: "HF", "CCSD(T)", "Full CI", …
    var methodTitle: String { MethodChoice.title(of: input.method) }

    /// True for H2, LiH and H2O.
    var isGeneral: Bool { input.isGeneral }
}

// MARK: - Jobs (cancellation and progress)

/// Progress of a running calculation, as reported by the backend.
nonisolated struct MacQCProgress: Sendable, Equatable {
    var stage: Int = 0        // 0 idle, 1 integrals, 2 SCF, 3 transformation, 4 Full CI, 5 done,
                              // 6 coupled cluster, 7 (T) correction (iteration = percent);
                              // 5.0: 8 gradient, 9 Hessian, 10 CPHF, 11 Lambda, 12 optimisation,
                              // 13 numerical Hessian
    var iteration: Int = 0
    var energy: Double = 0
    var residual: Double = 0

    var text: String {
        switch stage {
        case 1:
            return "Integrals…"
        case 2:
            return String(format: "SCF iteration %ld · E = %.8f E_h · gradient %.1e", iteration, energy, residual)
        case 3:
            return "Transforming the integrals to the RHF orbitals…"
        case 4:
            if iteration == 0 { return "Full CI: setting up the string Hamiltonian…" }
            return String(format: "Full CI · Davidson iteration %ld · E = %.8f E_h · residual %.1e",
                          iteration, energy, residual)
        case 5:
            return "Done"
        case 6:
            if iteration == 0 { return "Coupled cluster: setting up…" }
            return String(format: "Coupled cluster iteration %ld · E = %.8f E_h · max |Δt| %.1e",
                          iteration, energy, residual)
        case 7:
            return String(format: "(T) correction · %ld %% done", iteration)
        case 8:
            return "Gradient…"
        case 9:
            return "Hessian…"
        case 10:
            return String(format: "Orbital response · iteration %ld · residual %.1e", iteration, residual)
        case 11:
            return String(format: "Lambda equations · iteration %ld · max change %.1e", iteration, residual)
        case 12:
            return String(format: "Optimisation step %ld", iteration)
        case 13:
            return String(format: "Numerical Hessian · displacement %ld", iteration)
        default:
            return "Starting…"
        }
    }
}

/// A backend job handle: lets another thread cancel a running calculation and read its progress.
nonisolated final class MacQCJob: @unchecked Sendable {
    private let handle: OpaquePointer?

    init() {
        handle = macqc_job_create()
    }

    deinit {
        macqc_job_free(handle)
    }

    func cancel() {
        macqc_job_cancel(handle)
    }

    func progress() -> MacQCProgress {
        var stage: Int32 = 0
        var iteration: Int32 = 0
        var energy = 0.0
        var residual = 0.0
        macqc_job_progress(handle, &stage, &iteration, &energy, &residual)
        return MacQCProgress(stage: Int(stage), iteration: Int(iteration), energy: energy, residual: residual)
    }

    /// The output file written so far from byte `offset` on, and the total length in bytes
    /// (the next call passes that total as its offset).
    func output(from offset: Int) -> (text: String, total: Int) {
        var total: Int64 = 0
        guard let raw = macqc_job_output(handle, Int64(offset), &total) else { return ("", offset) }
        defer { macqc_free(raw) }
        return (String(cString: raw), Int(total))
    }

    /// One line about the current step of the job ("Optimisation step 4 · E = … · max force …").
    func progressText() -> String {
        guard let raw = macqc_job_progress_text(handle) else { return "" }
        defer { macqc_free(raw) }
        return String(cString: raw)
    }

    fileprivate var pointer: OpaquePointer? { handle }
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
    static func runSync(_ input: MacQCInput, job: MacQCJob? = nil) throws -> MacQCOutput {
        let encoder = JSONEncoder()
        let inputData = try encoder.encode(input)
        guard let inputJSON = String(data: inputData, encoding: .utf8) else {
            throw MacQCEngineError.encodingFailed
        }

        // macqc_run_job never throws and returns a malloc'ed, NUL-terminated UTF-8 string.
        guard let raw = macqc_run_job(inputJSON, job?.pointer) else {
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
    static func run(_ input: MacQCInput, job: MacQCJob? = nil) async throws -> MacQCOutput {
        let task = Task.detached(priority: .userInitiated) {
            try MacQCEngine.runSync(input, job: job)
        }
        return try await task.value
    }
}
