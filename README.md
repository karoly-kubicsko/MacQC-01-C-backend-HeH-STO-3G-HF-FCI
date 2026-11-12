# ΨLab — a native macOS quantum-chemistry app

**Current version: 2.0** — see [Versions](#versions). Version 3.0 (H₂, LiH and H₂O in
correlation-consistent basis sets) is being built in stages; its backend and terminal version
are described in [General molecules](#general-molecules-h2-lih-h2o-version-30-in-progress).

ΨLab computes the HeH⁺ molecular ion in a minimal STO-nG basis (n = 1, 2, 3): a restricted
Hartree–Fock calculation and then Full CI on top of it. It shows the energies, integrals and
matrices, the CI vector and the spin of every state, and scans the potential-energy curve from
0.8 to 5.0 bohr with the HF and FCI minima marked. SwiftUI with Liquid Glass on macOS 26; runs
on macOS 14 and later, on Apple Silicon and Intel Macs.

![App icon](app/icon/AppIcon-1024.png)

## Using the app

- **Download**: *Actions* ▸ *app* ▸ the newest run on `main` ▸ *Artifacts* ▸ **ΨLab-app**. Unzip
  `ΨLab-app.zip`, then the `ΨLab.app.zip` inside it. The app is not notarized: right-click
  `ΨLab.app` ▸ *Open* once, or run `xattr -dr com.apple.quarantine` on it.
- **Sidebar**: bond length R in Å and, below it, in bohr (both editable), Slater exponents ζ(He) and ζ(H), basis STO-1G/2G/3G, the
  exact-erf option and the appearance (System / Light / Dark). **Run** (⌘R) calculates;
  *Reset to Szabo–Ostlund defaults* restores R = 1.4632 bohr, ζ = 2.0925 / 1.24, STO-3G.
- **Tabs**: *Summary* (HF and FCI energies, correlation energy, convergence), *Hartree–Fock*
  (integrals, S/H/F/C/P matrices, orbital energies, SCF iterations), *Full CI* (determinants, CI
  vector, singlet weights, all states with ⟨S²⟩) and *Potential curve* (scan with a selectable step in Å).
- **Menus**: *File ▸ Export Results as Text…* (⇧⌘E) saves a `.txt` file with the summary, the full
  calculation report (every SCF iteration, the Full CI section) and, after a scan, the potential-curve
  table; *File ▸ Export Results as JSON…* (⌥⇧⌘E) saves the same results as JSON;
  *Edit ▸ Copy Summary as Text* (⇧⌘C) copies the short summary.
- **Build from source**: see [`app/SETUP.md`](app/SETUP.md) (Xcode 26, or XcodeGen with
  `xcode/project.yml`).

## Versions

### 2.0
- **Terminal version for Linux** (`psilab`, also runs on macOS): reads a plain-text input file
  (`basis`, `R` in bohr or Å, `zeta_He`, `zeta_H`, `exact_erf`, `scan`, `scan_step`) and writes a
  plain-text output file; installs as `psilab input.inp output.txt`.
- **Bond length in Å and bohr**: the sidebar shows R (Å) and below it R (bohr), both editable;
  the summary, the minima and the exported text give both units.
- **Export Results as Text** (⇧⌘E): summary, full Hartree–Fock and Full CI report and the
  potential-curve table; JSON export moved to ⌥⇧⌘E.
- **App icon** and the new name **ΨLab** (formerly MacQC).
- Potential curve: FCI blue circles, HF orange diamonds; the minima are listed in a readout box
  that no longer overlaps the curves.

### 1.0
- Native macOS app (SwiftUI, Liquid Glass on macOS 26, runs on macOS 14+ on Apple Silicon and Intel).
- HeH⁺ in an STO-1G/2G/3G basis: restricted Hartree–Fock (a C++ port of the Szabo & Ostlund
  program, output identical to it) and Full CI with ⟨S²⟩ spin labels and singlet weights.
- Tabs Summary, Hartree–Fock, Full CI and Potential curve (scan R = 0.8 … 5.0 bohr with minima).
- Light, Dark and System appearance; Export Results as JSON; Copy Summary as Text.

## Terminal version (Linux)

The same calculation runs without the graphical interface as the terminal program `psilab`, on
Linux (and macOS). It reads a plain-text input file and writes a plain-text output file.

Build and install once (needs a C++17 compiler and CMake):

```bash
cmake -S . -B build && cmake --build build
sudo cmake --install build                       # installs /usr/local/bin/psilab
```

Without administrator rights, install into your home folder instead with
`cmake --install build --prefix ~/.local` (then `~/.local/bin` must be on your `PATH`). After that
`psilab` runs from any folder:

```bash
psilab input.inp output.txt       # read input.inp, write the results to output.txt
psilab input.inp                  # output name defaults to input.out
psilab --example > input.inp      # start from a commented example input
```

The input file has one `key = value` per line; `#` starts a comment; keys are case-insensitive:

```text
basis     = STO-3G           # STO-1G, STO-2G or STO-3G
R         = 1.4632 bohr      # bond length; unit bohr (default) or angstrom, e.g. 0.7743 angstrom
zeta_He   = 2.0925           # Slater exponent of He
zeta_H    = 1.24             # Slater exponent of H
exact_erf = no               # no: the original program's erf; yes: the exact erf
scan      = yes              # potential curve, R = 0.8 ... 5.0 bohr
scan_step = 0.1              # bohr
```

Missing keys take the values shown (the Szabo–Ostlund case, without scan unless `scan = yes`).

For H₂, LiH and H₂O (restricted Hartree–Fock, see
[General molecules](#general-molecules-h2-lih-h2o-version-30-in-progress)) the input names the
molecule and a basis set; `R` and `angle` default to the experimental equilibrium geometry:

```text
molecule = H2O               # H2, LiH or H2O (default HeH+)
basis    = aug-cc-pVTZ       # aug-cc-pVXZ or aug-cc-pCVXZ, X = D, T, Q, 5
R        = 0.9572 angstrom   # O-H (H2: H-H, LiH: Li-H); default 0.74144 / 1.5949 / 0.9572 Å
angle    = 104.52            # H-O-H in degrees (H2O only)
method   = RHF
threads  = 0                 # 0: all processor cores
memory   = 1024 MB           # two-electron integrals are kept in memory up to this size
```

`psilab --estimate input.inp` prints the size of the calculation (basis functions, two-electron
integrals, memory) and a rough time estimate without computing; a calculation estimated to take
longer than 30 s prints a warning first. Calculations that cannot be done are refused with an
explanation (for example LiH in aug-cc-pCV5Z, which is not defined for Li). The basis-set files are
found in `basis_dir` (input key), `$PSILAB_BASIS_DIR`, the installed copy
(`<prefix>/share/psilab/basis`, put there by `cmake --install`) or the source tree.
The output file contains the same text as the app's *File ▸ Export Results as Text*: the summary,
the full Hartree–Fock and Full CI report and, with `scan = yes`, the potential-curve table (R in
bohr and Å) with the HF and FCI minima. Mistakes in the input are reported with their line number;
the exit status is 0 on success, 1 if the calculation failed (e.g. no SCF convergence) and 2 for
input or file errors.

## The C++ backend

The calculations are done by a C++17 backend (standard library plus the vendored
`nlohmann/json` single header, MIT) that the app compiles directly and calls through a small C
API. The rest of this README documents the backend; its code and API keep the original working
name `macqc`.

### What it computes

For HeH⁺ (two electrons, one s-type STO-nG contracted function per atom, n = 1, 2, 3):

1. **Restricted Hartree–Fock** – a port of the Szabo & Ostlund (*Modern Quantum Chemistry*,
   Appendix B) FORTRAN program, via `reference/hartree_fock_sto3g.py`. The routine structure is
   kept (`calculate_integrals`, `collect_matrices`, `form_g`, `diag_2x2`, `scf`; primitive
   overlap/kinetic/nuclear-attraction/two-electron integrals; Boys `F0`; `derfother`), as are the
   details that make the printed output identical to the FORTRAN program: `PI = 3.1415926535898`,
   the Abramowitz–Stegun `erf` approximation (optionally `std::erf`), the analytic 2×2
   diagonalisation with its `E2 <= E1` swap rule, and the convergence test
   `sqrt(Σ(P − Pold)² / 4) < 1e-4` (max. 25 iterations, `P = 0` start).
2. **Full CI** in the converged HF MO basis: `h` and `(pq|rs)` are transformed to the MO basis, the
   4 determinants with M_s = 0 in 4 spin orbitals are coupled with Slater–Condon rules, the
   Hamiltonian is diagonalised with a built-in Jacobi solver, ⟨S²⟩ is evaluated for every
   eigenvector (singlet/triplet labels), and the ground state's singlet configuration weights
   `(1σ)²`, `(1σ)¹(2σ)¹`, `(2σ)²` are reported. The correlation energy is
   `E_FCI − E_HF(total)`, with `E_HF` the converged SCF energy as printed.

Reference values (R = 1.4632 bohr, STO-3G, ζ = 2.0925 / 1.24): E_HF,total = −2.860661992 Eh
(6 iterations), E_FCI = −2.880711848 Eh, E_corr = −0.020049857 Eh.

### General molecules: H2, LiH, H2O (version 3.0, in progress)

ΨLab 3.0 adds H₂, LiH and H₂O (closed shell, restricted Hartree–Fock) in the
correlation-consistent basis sets aug-cc-pVXZ and aug-cc-pCVXZ (X = D, T, Q, 5), up to h functions
(l = 5) and up to 341 basis functions (H₂O in aug-cc-pCV5Z). Stage 1 (this part) is the integral
engine and RHF; Full CI and the app's user interface follow.

- **Basis sets**: the Gaussian94 files in [`basis/`](basis/README.md), downloaded from the Basis
  Set Exchange and read at run time (`basis.cpp`, G94 parser with `D` exponents, SP shells, scale
  factors). Spherical harmonics; primitive and contraction normalisation. aug-cc-pCVXZ is not
  defined for H, so aug-cc-pVXZ is used for H (and reported); aug-cc-pCV5Z does not exist for Li,
  so LiH/aug-cc-pCV5Z is refused.
- **Integrals** (`gto_integrals.cpp`): overlap and kinetic energy by Obara–Saika recurrences;
  nuclear attraction and electron repulsion by the Obara–Saika vertical recurrence per primitive
  quartet, contraction, then the Head-Gordon–Pople horizontal recurrence and the transformation to
  real solid harmonics (Helgaker, Jørgensen, Olsen, Eqs. 6.4.47–50). Boys function from a
  tabulated 7-term Taylor expansion with downward recursion (asymptotic form above T = 45),
  relative error ≈ 1e-15. Schwarz screening; 8-fold permutational symmetry.
- **RHF** (`rhf.cpp`): canonical orthogonalisation (overlap eigenvalues below 1e-7 are removed and
  reported), core-Hamiltonian guess, DIIS (8 vectors), convergence |ΔE| < 1e-10 Eₕ and
  max |X†(FPS − SPF)X| < 1e-7. The unique two-electron integrals are kept in memory when they fit
  into 1 GB (default), otherwise every Fock matrix is built integral-direct. Integrals and Fock
  builds run on all processor cores (`std::thread`). Householder/QL eigensolver.
- **Driver** (`molecular.cpp`): input validation, basis assignment, a size and time estimate
  before the calculation (warning above 30 s), the text report (geometry, basis per atom, SCF
  iterations, energy components, virial ratio, orbital energies).

**Verification** (`backend/tests/test_rhf.cpp`, run by CTest):

- RHF total energies agree with [PySCF](https://pyscf.org) 2.15 — an independent program reading
  the same `.gbs` files with its own parser, at identical geometries — to better than 1e-9 Eₕ in
  all 23 cases (`reference/pyscf_rhf_reference.py` → `backend/tests/data/rhf_reference.json`;
  CI runs the cases up to 130 functions, the workflow input *full* runs all). Timings with 4
  threads on a 2.1 GHz x86-64 machine: H₂O aug-cc-pVQZ 2.9 s (integrals in memory), aug-cc-pCVQZ
  29 s, aug-cc-pV5Z 2.4 min, aug-cc-pCV5Z 3.9 min (integral-direct); LiH aug-cc-pV5Z 34 s; all
  other cases below 2 s.
- Integral invariants (traces and norms of S, T, V; norm of the ERI tensor; Σ(pp|qq), Σ(pq|pq))
  agree with PySCF to 1e-10, for a test basis with every l up to 5 and for H₂O/aug-cc-pV(D,T)Z;
  individual shell quartets up to (hh|hh) agree to 1e-13.
- Published Hartree–Fock limits from numerical (basis-set-free) HF calculations: H₂ at
  R = 1.4 bohr, −1.133629571 Eₕ (Kołos & Roothaan, Rev. Mod. Phys. 32, 219 (1960) give −1.13363);
  LiH at R = 3.015 bohr, −7.987352 Eₕ (Laaksonen, Pyykkö & Sundholm, Comput. Phys. Rep. 4, 313
  (1986)). The aug-cc-pVXZ energies decrease monotonically towards the limit from above:
  H₂ 4.8e-3, 6.0e-4, 1.6e-4, 1.9e-5 Eₕ and LiH 3.2e-3, 5.6e-4, 1.4e-4, 2.0e-5 Eₕ above it for
  X = D, T, Q, 5.

RHF energies at the default geometries (E_h; identical to PySCF in every printed digit, with the same
threshold 1e-7 for removing near-linear dependencies — PySCF's own default, 1e-6, would remove one
function combination in LiH/aug-cc-pCVQZ, whose smallest overlap eigenvalue is 8.5e-7):

| Basis set | H₂ (R = 0.74144 Å) | LiH (R = 1.5949 Å) | H₂O (0.9572 Å, 104.52°) |
|---|---|---|---|
| aug-cc-pVDZ | -1.1287937367 (18) | -7.9841394069 (32) | -76.0414279605 (41) |
| aug-cc-pVTZ | -1.1330212852 (46) | -7.9867905250 (69) | -76.0606132999 (92) |
| aug-cc-pVQZ | -1.1334669066 (92) | -7.9872114378 (126) | -76.0660008181 (172) |
| aug-cc-pV5Z | -1.1336044053 (160) | -7.9873310001 (207) | -76.0673205023 (287) |
| aug-cc-pCVDZ | -1.1287937367 (18) | -7.9844061601 (36) | -76.0417253915 (45) |
| aug-cc-pCVTZ | -1.1330212852 (46) | -7.9869238626 (82) | -76.0608386487 (105) |
| aug-cc-pCVQZ | -1.1334669066 (92) | -7.9872421887 (155) | -76.0660945771 (201) |
| aug-cc-pCV5Z | -1.1336044053 (160) | not defined for Li | -76.0673328556 (341) |

(number of basis functions in parentheses; for H₂ the aug-cc-pCVXZ rows equal aug-cc-pVXZ, because
H takes the aug-cc-pVXZ functions)

### Layout

```
backend/include/macqc/   public headers (macqc_api.h is the C API)
backend/src/             linalg, integrals, scf, fci, report, calculation, api (HeH+, STO-nG);
                         molecule, basis, gto_integrals, rhf, molecular (H2, LiH, H2O)
basis/                   Gaussian94 basis-set files from basissetexchange.org (+ README)
backend/cli/main.cpp     developer command-line tool `macqc` (flags / JSON)
backend/cli/psilab*      terminal version `psilab` (text input file -> text output file)
examples/                example input files for `psilab` (heh.inp, h2o.inp)
backend/tests/           test program (CTest), expected reports in tests/data/
reference/               the verified Python reference (HeH+); PySCF reference scripts (H2, LiH, H2O)
third_party/nlohmann/    json.hpp (MIT)
CMakeLists.txt           library + CLI + tests (for testing; Xcode compiles the sources directly)
```

### Build, test and run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

```bash
./build/macqc                     # text report (default: HeH+, STO-3G, print level 2)
./build/macqc --n 2 --R 1.6 --print 1
./build/macqc --json              # JSON instead of the text report
./build/macqc --input in.json     # read the API input JSON ("-" = stdin); flags override it
./build/macqc --help
```

Options: `--R --zeta1 --zeta2 --ZA --ZB --n --exact-erf --max-iter --conv --print --json --input`.
General molecules go through `--input`, e.g.
`echo '{"molecule":"H2O","basis":"aug-cc-pVTZ"}' | ./build/macqc --input -`.
The exit status is 0 on success and 1 (message on stderr) if the input is invalid or the SCF did
not converge.

The expected HF reports in `backend/tests/data/` are generated by the reference script
(`python3 backend/tests/gen_expected.py`, needs NumPy) and must never be edited by hand; the C++
text report is compared with them byte for byte (print levels 0, 1, 2; STO-1G/2G/3G; R = 2.5).

### Text report

The HF part reproduces the output of `reference/hartree_fock_sto3g.py` exactly (headings, spacing,
MATOUT layout, Fortran `D` numbers). The Full CI section is appended after it:
determinants, (level 2: MO `h`, CI Hamiltonian and MO integrals), all eigenvalues with ⟨S²⟩ and
spin label, the ground-state CI vector, singlet weights, FCI energies and the correlation energy.
Print level 0 is silent (empty report, use the JSON fields), level 1 is the summary, level 2 prints
every SCF iteration.

### C API

```c
char* macqc_run(const char* input_json);   // never throws; returns malloc'ed JSON
void  macqc_free(char* p);                 // release the string
```

#### Input JSON

```json
{"R":1.4632,"zeta1":2.0925,"zeta2":1.24,"ZA":2,"ZB":1,"nGaussians":3,
 "exactErf":false,"maxIter":25,"convergence":1e-4,"printLevel":2}
```

All keys are optional (the values above are the defaults). Validation: `R`, `zeta1`, `zeta2`,
`ZA`, `ZB`, `convergence` > 0; `nGaussians` ∈ {1,2,3}; `maxIter` ≥ 1; `printLevel` ∈ {0,1,2}.
Invalid or malformed input yields an error JSON; the call never crashes.

General molecules are selected with `"molecule"` (`"HeH+"` is the default and keeps everything
above unchanged):

```json
{"molecule":"H2O","basis":"aug-cc-pVTZ","R":1.808846,"angle":104.52,"method":"RHF",
 "basisDirectory":"/path/to/basis","threads":0,"memoryMB":1024,"scfMaxIter":100,
 "energyConvergence":1e-10,"gradientConvergence":1e-7,"estimateOnly":false,"printLevel":2}
```

`R` is in bohr (H₂O: O–H) and defaults to the experimental equilibrium bond length; `angle`
(degrees, H₂O) defaults to 104.52; `basisDirectory` defaults to `$PSILAB_BASIS_DIR` or the
compiled-in source folder. With `"estimateOnly":true` only the molecule, the basis set, the size
and time estimate and the warnings are returned.

#### Output JSON

```text
{
  "input":  { ...the (defaulted) input echoed... },
  "error":  null | "message",
  "hf": {
    "converged", "iterations", "electronicEnergy", "nuclearRepulsion", "totalEnergy",
    "orbitalEnergies": [e1, e2],
    "C", "P", "S", "X", "H", "F", "G", "PS":   2x2 matrices as [[row0],[row1]],
    "integrals": {"S12","T11","T12","T22","V11A","V12A","V22A","V11B","V12B","V22B",
                  "V1111","V2111","V2121","V2211","V2221","V2222"},
    "scfHistory": [{"iteration","electronicEnergy","deltaP"}, ...]
  },
  "fci": {
    "totalEnergy", "electronicEnergy", "correlationEnergy",
    "states": [{"energy","S2","spin","vector"}, ...],      // ascending, total energies
    "groundStateVector": [c1..c4],
    "determinants": ["|1α 1β>", "|1α 2β>", "|2α 1β>", "|2α 2β>"],
    "singletWeights": [{"configuration":"(1sigma)^2","weight":…}, …]   // ground state
  },
  "report": "FORTRAN-style text report"
}
```

On invalid input `hf` and `fci` are `null`. If the SCF does not converge, `hf` is present with
`converged:false`, `fci` is `null` and `error` is set.

For general molecules the output is

```text
{
  "input": {...}, "error": null | "message",
  "molecule": {"name","electrons","charge","R","angle","nuclearRepulsion",
               "atoms":[{"symbol","Z","position":[x,y,z]}]},                    // bohr
  "basis":    {"name","functions","shells","maxL",
               "atoms":[{"symbol","basis","contraction","primitives","functions"}]},
  "estimate": {"basisFunctions","uniqueIntegrals","integralMemoryMB","integralsStored",
               "threads","rhfSeconds","longCalculation"},
  "notes":    ["aug-cc-pCVDZ is not defined for H ..."],
  "warnings": ["Long calculation: ... Estimated time about 6 min on 8 threads."],
  "rhf": {"converged","iterations","totalEnergy","electronicEnergy","nuclearRepulsion",
          "kineticEnergy","nuclearAttractionEnergy","electronRepulsionEnergy","virialRatio",
          "orbitalEnergies","occupiedOrbitals","homo","lumo","smallestOverlapEigenvalue",
          "removedFunctions","integralsStored","uniqueIntegrals","threads",
          "integralSeconds","scfSeconds",
          "scfHistory":[{"iteration","totalEnergy","deltaE","gradient","diis"}]} | null,
  "hf": null, "fci": null,
  "report": "text report"
}
```

A refused calculation (e.g. LiH/aug-cc-pCV5Z) has `error` set and `rhf` `null`.

### Adding to Xcode

How the backend is wired into the ΨLab Xcode target (the full step-by-step guide for the app is
[`app/SETUP.md`](app/SETUP.md); `xcode/project.yml` encodes the same settings for XcodeGen):

1. Copy or reference this repository's folders in the Xcode project (*File ▸ Add Files to …*,
   "Create folder references" or groups; do **not** tick "Copy items" if you want to keep a single
   source of truth).
2. In the app target's **Build Phases ▸ Compile Sources** add exactly these files:
   `backend/src/linalg.cpp`, `integrals.cpp`, `report.cpp`, `scf.cpp`, `fci.cpp`,
   `molecule.cpp`, `basis.cpp`, `gto_integrals.cpp`, `rhf.cpp`, `molecular.cpp`,
   `calculation.cpp`, `api.cpp`. Do **not** add `backend/cli/main.cpp` or anything in
   `backend/tests/` (they have their own `main`). Headers (`backend/include/macqc/*.h` and
   `third_party/nlohmann/json.hpp`) need no compile phase but should be in the project.
3. **Build Settings ▸ Header Search Paths** (`HEADER_SEARCH_PATHS`): add
   `$(SRCROOT)/backend/include` and `$(SRCROOT)/third_party` (adjust the prefix to where the folders
   live relative to the `.xcodeproj`).
4. **Build Settings ▸ C++ Language Dialect** (`CLANG_CXX_LANGUAGE_STANDARD`): `C++17`
   (or `GNU++17`). C++ standard library: `libc++` (default).
5. Recommended: **Other C++ Flags** `-ffp-contract=off`. Apple clang contracts `a*b+c` into fused
   multiply-adds on Apple Silicon by default, which can change the last bits of the results
   relative to the reference output (the CMake build and CI use this flag).
6. **Build Settings ▸ Objective-C Bridging Header** (`SWIFT_OBJC_BRIDGING_HEADER`): point to a
   header such as `MacQC-Bridging-Header.h` containing

   ```c
   #include "macqc/macqc_api.h"
   ```

   The header is plain C (`extern "C"` guarded), so Swift imports `macqc_run`/`macqc_free`
   directly; no Objective-C++ wrapper is needed.
7. Call it from Swift:

   ```swift
   import Foundation

   func runMacQC(inputJSON: String) -> Data? {
       guard let raw = macqc_run(inputJSON) else { return nil }   // nil only if out of memory
       defer { macqc_free(raw) }
       return Data(bytes: raw, count: strlen(raw))                // UTF-8 JSON
   }

   // let data = runMacQC(inputJSON: #"{"R":1.4632,"nGaussians":3}"#)
   // let result = try JSONSerialization.jsonObject(with: data!)  // or Codable
   ```

   Run it off the main thread for scans; the library keeps no global state, so concurrent calls
   are safe.

## App sources

The SwiftUI sources of ΨLab are in `app/MacQC/`, the app icon in `app/MacQC/Assets.xcassets`
(regenerated by `app/icon/make_icon.py`, master image `app/icon/AppIcon-1024.png`), and the
XcodeGen project spec in `xcode/project.yml`. See [`app/SETUP.md`](app/SETUP.md) for building.

## Continuous integration

- `.github/workflows/backend.yml` builds the C++ backend and runs its tests on `ubuntu-latest`
  (Python 3 + NumPy regenerate the expected reports and verify they are unchanged). A second job
  installs PySCF and recomputes the stored integral and RHF reference values; run manually with
  the input *full*, a third job checks every RHF case up to 341 basis functions.
- `.github/workflows/app.yml` builds the universal ΨLab app on a `macos-26` runner (XcodeGen +
  Xcode), runs a UI smoke test (launch, all tabs, potential-curve scan, menus, light and dark
  appearance, screenshots) and publishes the **ΨLab-app** download.
