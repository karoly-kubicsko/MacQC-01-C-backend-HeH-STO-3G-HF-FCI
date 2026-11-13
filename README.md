# ΨLab — a native macOS quantum-chemistry app

**Current version: 3.1** — see [Versions](#versions).

ΨLab computes small molecules with restricted Hartree–Fock, coupled cluster and Full CI:

- **H₂, LiH and H₂O** in the correlation-consistent basis sets aug-cc-pVXZ and aug-cc-pCVXZ
  (X = D, T, Q, 5; spherical Gaussians up to h functions, from basissetexchange.org) with the
  methods **HF, CCSD, CCSD(T), CCSDT, CCSDTQ and Full CI** (all electrons correlated): RHF and
  CCSD(T) for all of them, CCSDT and CCSDTQ where they fit into memory, Full CI for H₂ and LiH
  where it is feasible. Before a calculation ΨLab shows its size (basis functions, two-electron
  integrals, amplitudes, determinants), memory and estimated time; it warns about long
  calculations and refuses impossible ones with an explanation (e.g. Full CI for H₂O:
  5.6·10¹¹ determinants already in aug-cc-pVDZ).
- **HeH⁺** in a minimal STO-nG basis (n = 1, 2, 3), the Szabo–Ostlund calculation: Hartree–Fock
  identical to the book's program, Full CI with all states and their spin.

It shows energies, orbital energies, SCF, coupled-cluster and Davidson iterations, the largest
cluster amplitudes, the leading CI determinants and natural-orbital occupations, the full text report, and scans potential-energy curves with the
minima marked. SwiftUI with Liquid Glass on macOS 26; runs on macOS 14 and later, on Apple
Silicon and Intel Macs. The same calculations run in the terminal on Linux and macOS
([`psilab`](#terminal-version-linux)).

![App icon](app/icon/AppIcon-1024.png)

## Using the app

- **Download**: *Actions* ▸ *app* ▸ the newest run on `main` ▸ *Artifacts* ▸ **ΨLab-app**. Unzip
  `ΨLab-app.zip`, then the `ΨLab.app.zip` inside it. The app is not notarized: right-click
  `ΨLab.app` ▸ *Open* once, or run `xattr -dr com.apple.quarantine` on it.
- **Sidebar**: the **molecule** (HeH⁺, H₂, LiH, H₂O) and, for H₂, LiH and H₂O:
  - the geometry: bond length in Å and in bohr (both editable) and, for H₂O, the H–O–H angle;
    *Experimental geometry* restores 0.74144 Å (H₂), 1.5949 Å (LiH), 0.9572 Å and 104.52° (H₂O);
  - the **method** — HF, CCSD, CCSD(T), CCSDT, CCSDTQ or Full CI (a menu, in this order) — and
    the **basis set** (aug-cc-pVDZ … aug-cc-pV5Z, aug-cc-pCVDZ … aug-cc-pCV5Z);
  - *Size of the calculation*, updated as you type: basis functions per atom, two-electron
    integrals (kept in memory or integral-direct), coupled-cluster amplitudes and memory, Full CI
    determinants and memory, estimated time.
    Notes (e.g. H takes aug-cc-pVXZ in the core-valence sets) in grey, warnings in orange, and a
    red explanation when the calculation is not possible (then **Run** is disabled); the run bar and,
    before the first result, the detail area repeat the refusal or the long-calculation warning.

  For HeH⁺: R in Å and bohr, Slater exponents ζ(He) and ζ(H), STO-1G/2G/3G, the exact-erf option;
  *Reset to Szabo–Ostlund defaults* restores R = 1.4632 bohr, ζ = 2.0925 / 1.24, STO-3G. Below:
  the appearance (System / Light / Dark).
- **Run** (⌘R) calculates on all processor cores. A calculation estimated to take longer than
  30 s asks for confirmation first; while it runs, the sidebar shows the progress (SCF iteration,
  integral transformation, coupled-cluster iteration with energy and largest amplitude change,
  progress of the (T) correction, Davidson iteration with energy and residual) and the elapsed time, and
  **Cancel** (⌘.) stops it.
- **Tabs**: *Summary* (energies, geometry, basis set per atom, convergence, timings),
  *Hartree–Fock* (energy components and virial ratio, orbital energies, SCF iterations; for HeH⁺
  also the integrals and the S/H/F/C/P matrices), *Coupled cluster* (energies, (T) correction,
  T1 diagnostic, amplitude counts and norms, largest amplitudes, convergence of the iterations),
  *Full CI* (energies, ⟨S²⟩, leading determinants,
  natural-orbital occupations, Davidson iterations; for HeH⁺ the CI vector, singlet weights and all
  states), *Potential curve* (bond-length scan with a selectable step in Å — H₂O as a symmetric
  stretch — with the RHF minimum and that of the selected correlated method; long scans ask for
  confirmation and can be stopped) and *Report* (the complete text report); *View ▸ Summary …
  Report* (⌘1 … ⌘6) switches tabs.
- **Menus**: *File ▸ Export Results as Text…* (⇧⌘E) saves a `.txt` file with the summary, the full
  calculation report (every SCF iteration, the coupled-cluster or Full CI section) and, after a
  scan, the potential-curve
  table; *File ▸ Export Results as JSON…* (⌥⇧⌘E) saves the same results as JSON;
  *Edit ▸ Copy Summary as Text* (⇧⌘C) copies the short summary.
- **Build from source**: see [`app/SETUP.md`](app/SETUP.md) (Xcode 26, or XcodeGen with
  `xcode/project.yml`).

## Versions

### 3.1
- **Coupled cluster**: CCSD, CCSD(T), CCSDT and CCSDTQ for H₂, LiH and H₂O, all electrons
  correlated, on the RHF orbitals. The method menu offers HF, CCSD, CCSD(T), CCSDT, CCSDTQ and
  Full CI.
  - CCSD and CCSD(T): closed-shell, spin-adapted equations; the (vv|vv) ladder runs over AO
    integrals (in memory or recomputed), so H₂O CCSD(T) works up to aug-cc-pCV5Z.
  - CCSDT and CCSDTQ: a general-order spin-orbital engine whose 183 diagrams are generated and
    checked by `backend/tools/cc_diagrams.py`.
- Verified: CCSD and CCSD(T) agree with PySCF to 1e-10 Eₕ (H₂, LiH and H₂O, seven basis sets);
  CCSD, CCSDT and CCSDTQ agree with brute-force coupled cluster in the determinant space to 1e-12;
  CCSDTQ reproduces Full CI for LiH (four electrons) and every method reproduces Full CI for H₂.
- Estimates of amplitudes, memory and time; refusals above the memory limit (60 % of the RAM in
  the app) or three days, warnings for long calculations, notes where a method is exact.
- New *Coupled cluster* tab; potential curves with any method; matrix products through Apple's
  Accelerate in the app.
- Terminal version: `method = HF | CCSD | CCSD(T) | CCSDT | CCSDTQ | FCI`.

### 3.0
- **H₂, LiH and H₂O** in eight correlation-consistent basis sets (aug-cc-pVXZ, aug-cc-pCVXZ,
  X = D … 5), read from Gaussian94 files from basissetexchange.org (`basis/`): a general Gaussian
  integral engine up to h functions, restricted Hartree–Fock (DIIS, integrals in memory or
  integral-direct, multithreaded) and a determinant-based Full CI (H₂ in all basis sets, LiH up to
  aug-cc-pCVTZ with 11 million determinants).
- Verified against PySCF (all 23 RHF energies and every Full CI energy identical to 10 decimals),
  against dense-matrix Full CI, and against published Hartree–Fock limits and exact energies.
- Size, memory and time estimates before running; warnings for long calculations; impossible
  ones are refused with an explanation (Full CI for H₂O, LiH from aug-cc-pVQZ on, LiH in
  aug-cc-pCV5Z which is not defined for Li).
- App: molecule, method and basis-set pickers, geometry inputs (bond length; H₂O bond length and
  angle), live size estimate, progress and cancellation of running calculations, new Summary,
  Hartree–Fock and Full CI views for the new molecules, a Report tab, potential curves for all
  molecules.
- Terminal version: `molecule`, `basis`, `angle`, `method = RHF | FCI`, `threads`, `memory`,
  `fci_memory`, `basis_dir` and `psilab --estimate`.

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
[General molecules](#general-molecules-h2-lih-h2o-version-30)) the input names the
molecule and a basis set; `R` and `angle` default to the experimental equilibrium geometry:

```text
molecule = H2O               # H2, LiH or H2O (default HeH+)
basis    = aug-cc-pVTZ       # aug-cc-pVXZ or aug-cc-pCVXZ, X = D, T, Q, 5
R        = 0.9572 angstrom   # O-H (H2: H-H, LiH: Li-H); default 0.74144 / 1.5949 / 0.9572 Å
angle    = 104.52            # H-O-H in degrees (H2O only)
method   = RHF               # HF (RHF), CCSD, CCSD(T), CCSDT, CCSDTQ, or FCI (H2 and LiH up to aug-cc-pCVTZ)
threads  = 0                 # 0: all processor cores
memory   = 1024 MB           # two-electron integrals are kept in memory up to this size
fci_memory = 4 GB            # memory limit of Full CI and coupled cluster
```

`psilab --estimate input.inp` prints the size of the calculation (basis functions, two-electron
integrals, memory) and a rough time estimate without computing; a calculation estimated to take
longer than 30 s prints a warning first. Calculations that cannot be done are refused with an
explanation (for example LiH in aug-cc-pCV5Z, which is not defined for Li, or Full CI for H₂O,
which would need 5.6·10¹¹ determinants already in aug-cc-pVDZ). The basis-set files are
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

### General molecules: H2, LiH, H2O (version 3.0)

ΨLab 3.0 adds H₂, LiH and H₂O (closed shell) in the correlation-consistent basis sets aug-cc-pVXZ
and aug-cc-pCVXZ (X = D, T, Q, 5), up to h functions (l = 5) and up to 341 basis functions (H₂O in
aug-cc-pCV5Z), with restricted Hartree–Fock and Full CI, in the app and in `psilab`.

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

- **Full CI** (`mo_integrals.cpp`, `ci.cpp`): all electrons correlated, determinants with
  M_S = 0, lowest singlet state. The two-electron integrals are transformed to the RHF orbitals
  integral-direct (shell quartets → (μν|rs) → (pq|rs), matrix products). The CI vector is a matrix
  C(I_α, I_β) over occupation strings; σ = HC is built from the one-spin string Hamiltonian F
  (σ = CF + FC) and, for the opposite-spin part, the (N−1)-electron resolution of the identity:
  D_jl(K_α,K_β) = ⟨K_α|a_j|J_α⟩⟨K_β|a_l|J_β⟩C(J_α,J_β), E = (ij|kl)·D as one matrix product
  (cost n⁶ for two electrons per spin, instead of n⁸), then the scatter back with a⁺_i a⁺_k.
  Davidson method with diagonal preconditioner; the trial vectors are kept symmetric (C = Cᵀ), which
  removes triplet components. Output: energy, correlation energy, ⟨S²⟩ (= |S₊Ψ|²), weight of the
  RHF determinant, leading determinants, natural-orbital occupation numbers.
- **Limits**: Full CI is refused, with an explanation, above 2·10⁷ determinants (160 MB per CI
  vector) or above the memory limit (`fciMemoryMB`, default 4 GB). Determinants per case:

  | Basis set | H₂ | LiH | H₂O |
  |---|---|---|---|
  | aug-cc-pVDZ | 324 | 246,016 | 5.6e11 — refused |
  | aug-cc-pVTZ | 2,116 | 5,503,716 | 2.4e15 — refused |
  | aug-cc-pVQZ | 8,464 | 62,015,625 — refused | 1.4e18 — refused |
  | aug-cc-pV5Z | 25,600 | 454,585,041 — refused | 2.5e20 — refused |
  | aug-cc-pCVDZ | 324 | 396,900 | 1.5e12 — refused |
  | aug-cc-pCVTZ | 2,116 | 11,029,041 | 9.3e15 — refused |
  | aug-cc-pCVQZ | 8,464 | 142,444,225 — refused | 6.8e18 — refused |
  | aug-cc-pCV5Z | 25,600 | – (basis not defined for Li) | 1.4e21 — refused |

  Before a calculation ΨLab estimates its memory and time and warns above 30 s.

**Verification** (`backend/tests/test_rhf.cpp`, `test_fci.cpp`, run by CTest):

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

- Full CI energies agree with PySCF to every printed digit (10 decimals; `reference/pyscf_fci_reference.py`
  → `backend/tests/data/fci_reference.json`): H₂ in aug-cc-pV{D,T,Q,5}Z (reference: PySCF CCSD, which
  is exact for two electrons) and LiH in aug-cc-pVDZ, aug-cc-pCVDZ, aug-cc-pVTZ (PySCF's determinant
  Full CI, `fci.direct_spin0`). LiH/aug-cc-pCVTZ (11 million determinants) runs through the same
  code; it was not recomputed with PySCF, whose solver needed 2.2 h already for aug-cc-pVTZ (ΨLab:
  57 s). The string-driven σ algorithm is also compared with
  a dense Hamiltonian matrix built from spin-orbital Slater–Condon rules for 1, 2 and 5 electrons
  per spin (H₂, LiH and H₂O in truncated orbital spaces): agreement to 1e-12 Eₕ; ⟨S²⟩ = 0 to 1e-8.
- Exact non-relativistic energies: Full CI lies above them and approaches them with the basis set,
  e.g. H₂ at R = 1.4 bohr (exact −1.1744757142 Eₕ, Pachucki, Phys. Rev. A 82, 032509 (2010)):
  9.9e-3, 1.8e-3, 6.1e-4 Eₕ above it in aug-cc-pV{D,T,Q}Z. For LiH (exact ≈ −8.07055 Eₕ at
  3.015 bohr, explicitly correlated Gaussians) the core-valence set aug-cc-pCVDZ recovers most of the
  Li 1s correlation that aug-cc-pVDZ misses (1.9e-2 instead of 4.9e-2 Eₕ above the exact energy).

Full CI energies at the default geometries (E_h):

| Basis set | H₂ | LiH |
|---|---|---|
| aug-cc-pVDZ | −1.1646244795 | −8.0212084789 |
| aug-cc-pVTZ | −1.1726340136 | −8.0386964311 |
| aug-cc-pVQZ | −1.1738672462 | refused |
| aug-cc-pV5Z | −1.1742521870 | refused |
| aug-cc-pCVDZ | (= aug-cc-pVDZ) | −8.0518951298 |
| aug-cc-pCVTZ | (= aug-cc-pVTZ) | −8.0652388024 |

Timings with 4 threads (2.1 GHz x86-64; RHF, integral transformation and Davidson): H₂/aug-cc-pV5Z
21 s (1.6 GB), LiH/aug-cc-pVTZ 57 s (1.5 GB), LiH/aug-cc-pCVTZ 2.7 min (2.8 GB); all other Full CI
cases take a few seconds. The built-in estimates are 20–40 % above these times and within 5 % of
the memory.

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

### Coupled cluster (version 3.1)

ΨLab 3.1 adds CCSD, CCSD(T), CCSDT and CCSDTQ for H₂, LiH and H₂O on the canonical RHF orbitals,
all electrons correlated (no frozen core); the Fock matrix is used in full.

- **CCSD and CCSD(T)** (`cc_closed_shell.cpp`): closed-shell, spin-adapted CCSD equations
  (Stanton–Gauss intermediates in the form of PySCF's `rccsd`) written as matrix products. The
  integrals with at least one occupied index come from a transformation that keeps one occupied
  index in the ket (memory o·n·N(N+1)/2). The particle–particle ladder Σ_cd (ac|bd) τ_ij^cd runs in
  the AO basis: τ_ij is back-transformed with the virtual orbitals, split into its symmetric and
  antisymmetric parts, contracted with the AO integrals like an exchange matrix (o² "densities"
  at once) and transformed back. The AO integrals are kept in memory when they fit, otherwise
  recomputed with Schwarz screening, so no (vv|vv) integrals are stored. (T):
  E = Σ_ijk Σ_abc (4W_abc + W_bca + W_cab)(V_abc − V_cba)/(3D_ijk^abc), looping over i ≥ j ≥ k
  with matrix products (cost o³v⁴).
- **CCSDT and CCSDTQ** (`cc_general.cpp`): a general-order spin-orbital engine. Amplitudes t(A;I)
  are stored once per pair of sorted sets of spin orbitals, blocked by the number of α orbitals
  and addressed with the combinatorial number system (6.8·10⁷ for CCSDTQ of H₂O/aug-cc-pVDZ
  instead of 2.7·10¹¹ in a full spin-orbital array). The residuals ⟨μ|(H e^T)_C|0⟩ are the 183
  connected diagrams up to quadruples from `backend/tools/cc_diagrams.py`, which enumerates them,
  orders the cluster vertices of each for the lowest cost and determines every sign by fitting to
  a brute-force evaluation of ⟨μ|e^−T H e^T|0⟩ in a determinant space (exact to 1e-13). Each
  diagram is evaluated from the Hamiltonian vertex outwards, one cluster vertex at a time, as one
  matrix product per spin case (gather, product, signed scatter); diagrams that begin alike share
  their intermediates. The engine also runs CCSD (used for cross-checks).
- Both engines: Jacobi updates with DIIS, converged when the largest amplitude change is below
  1e-8 and |ΔE| below 1e-10 Eₕ. The report lists the energies, the (T) correction, the T1
  diagnostic, the amplitude norms per rank and the largest amplitudes.
- **Limits** (`cc_estimate.cpp`): memory and time are estimated before running; a calculation is
  refused above `fciMemoryMB` (60 % of the RAM in the app, at most 48 GB) or above three days.
  CCSDTQ for H₂O needs about 5.5 GB in aug-cc-pVDZ and is refused in larger basis sets; CCSDT
  for H₂O runs up to aug-cc-pVTZ (about 4 GB, estimated 77 min on 4 threads). For two electrons
  every method equals Full CI; for four electrons (LiH) CCSDTQ equals Full CI. ΨLab adds a note
  in both cases.

Verification:

- CCSD and CCSD(T) agree with PySCF (`reference/pyscf_cc_reference.py`: same basis files,
  geometries and linear-dependency threshold) within 3e-10 Eₕ for H₂ (aug-cc-pVDZ, aug-cc-pVTZ),
  LiH (aug-cc-pVDZ, aug-cc-pCVDZ) and H₂O (aug-cc-pVDZ, aug-cc-pCVDZ, aug-cc-pVTZ), with the AO
  integrals of the ladder in memory and recomputed.
- CCSD, CCSDT and CCSDTQ from both engines agree within 1e-12 Eₕ with brute-force coupled cluster
  in the determinant space (`backend/tests/gen_cc_reference.py`, no coupled-cluster equations
  coded): random integrals with a non-diagonal Fock matrix (6 electrons in 7 orbitals), H₂O/STO-3G
  and LiH/6-31G.
- Exactness: CCSD, CCSDT and CCSDTQ equal Full CI for H₂, and CCSDTQ equals Full CI for
  LiH/aug-cc-pVDZ (−8.0212084789 Eₕ, 246,016 determinants).

Energies at the default geometries (Eₕ, aug-cc-pVDZ):

| Method | H₂ | LiH | H₂O |
| --- | --- | --- | --- |
| RHF | −1.1287937367 | −7.9841394069 | −76.0414279605 |
| CCSD | −1.1646244795 | −8.0211854592 | −76.2707727838 |
| CCSD(T) | −1.1646244795 | −8.0212069372 | −76.2760053269 |
| CCSDT | −1.1646244795 | −8.0212097170 | −76.2762366785 |
| CCSDTQ | −1.1646244795 | −8.0212084789 | (about 5.5 GB) |
| Full CI | −1.1646244795 | −8.0212084789 | not possible |

Timings with 4 threads (2.1 GHz x86-64, portable matrix products; the app uses Accelerate):
H₂O CCSD(T) in aug-cc-pVDZ 0.5 s, aug-cc-pVTZ 6 s, aug-cc-pVQZ 51 s; H₂O CCSDT/aug-cc-pVDZ
2.6 min (20 iterations); LiH CCSDTQ/aug-cc-pVDZ 1.6 min.

### Layout

```
backend/include/macqc/   public headers (macqc_api.h is the C API)
backend/src/             linalg, integrals, scf, fci, report, calculation, api (HeH+, STO-nG);
                         molecule, basis, gto_integrals, rhf, molecular (H2, LiH, H2O);
                         gemm, mo_integrals, ci (Full CI for H2, LiH); job (cancel/progress);
                         cc_closed_shell (CCSD, (T)), cc_general (CCSDT, CCSDTQ), cc_estimate
backend/tools/           cc_diagrams.py: generates and verifies the coupled-cluster diagrams
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

// Jobs: cancellation and progress of long calculations (thread-safe).
macqc_job* macqc_job_create(void);
char* macqc_run_job(const char* input_json, macqc_job* job);   // like macqc_run
void  macqc_job_cancel(macqc_job* job);    // from another thread; the run then returns
                                           // "error": "calculation cancelled"
void  macqc_job_progress(const macqc_job* job, int* stage, int* iteration,
                         double* energy, double* residual);
      // stage 0 idle, 1 integrals, 2 SCF, 3 integral transformation, 4 Full CI, 5 done,
      // 6 coupled-cluster iterations (residual = largest amplitude change), 7 (T) correction
      // (iteration = percent done)
void  macqc_job_free(macqc_job* job);      // after macqc_run_job has returned
```

The app polls `macqc_job_progress` four times a second for the progress line of the sidebar.

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
 "basisDirectory":"/path/to/basis","threads":0,"memoryMB":1024,"fciMemoryMB":4096,
 "scfMaxIter":100,"energyConvergence":1e-10,"gradientConvergence":1e-7,"estimateOnly":false,
 "printLevel":2}
```

`R` is in bohr (H₂O: O–H) and defaults to the experimental equilibrium bond length; `angle`
(degrees, H₂O) defaults to 104.52; `method` is `"RHF"` (also `"HF"`), `"CCSD"`, `"CCSD(T)"`,
`"CCSDT"`, `"CCSDTQ"` or `"FCI"` (the correlated methods follow RHF); `fciMemoryMB` is the
memory limit of Full CI and coupled cluster;
`basisDirectory` defaults to `$PSILAB_BASIS_DIR` or the compiled-in source folder. With
`"estimateOnly":true` only the molecule, the basis set, the size and time estimate, the warnings
and (for a Full CI or coupled-cluster calculation that is not possible) the refusal are
returned, without an error.

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
               "threads","rhfSeconds","totalSeconds","longCalculation",
               "fci": {"orbitals","alphaElectrons","betaElectrons","alphaStrings","betaStrings",
                       "determinants","memoryMB","seconds","feasible","reason"} | null,
               "coupledCluster": {"method","engine","occupied","virtuals","amplitudes","memoryMB",
                       "ladderIntegralsStored","iterations","iterationSeconds","triplesSeconds",
                       "seconds","feasible","reason"} | null},
  "notes":    ["aug-cc-pCVDZ is not defined for H ..."],
  "warnings": ["Long calculation: ... Estimated time about 6 min on 8 threads."],
  "refusal":  null | "Full CI is not possible here: ...",
  "rhf": {"converged","iterations","totalEnergy","electronicEnergy","nuclearRepulsion",
          "kineticEnergy","nuclearAttractionEnergy","electronRepulsionEnergy","virialRatio",
          "orbitalEnergies","occupiedOrbitals","homo","lumo","smallestOverlapEigenvalue",
          "removedFunctions","integralsStored","uniqueIntegrals","threads",
          "integralSeconds","scfSeconds",
          "scfHistory":[{"iteration","totalEnergy","deltaE","gradient","diis"}]} | null,
  "fullCI": {"converged","iterations","totalEnergy","electronicEnergy","correlationEnergy","S2",
             "hfWeight","orbitals","alphaElectrons","betaElectrons","alphaStrings","betaStrings",
             "determinants","naturalOccupations","transformSeconds","ciSeconds",
             "leadingDeterminants":[{"coefficient","alpha":[MOs],"beta":[MOs],"excitation"}],
             "history":[{"iteration","totalEnergy","residual","subspace"}]} | null,
  "coupledCluster": {"method","engine","converged","iterations","occupied","virtuals",
             "referenceEnergy","correlationEnergy","ccEnergy","triplesCorrection" | null,
             "totalEnergy","t1Diagnostic","amplitudeCounts":[by rank],"amplitudeNorms":[by rank],
             "leadingAmplitudes":[{"value","rank","occupied":[MOs],"virtuals":[MOs],"spins"}],
             "history":[{"iteration","totalEnergy","correlationEnergy","deltaE","residual"}],
             "ladder","transformSeconds","ccSeconds","triplesSeconds"} | null,
  "hf": null, "fci": null,
  "report": "text report"
}
```

A refused calculation (e.g. LiH/aug-cc-pCV5Z, Full CI for H₂O, CCSDTQ for H₂O/aug-cc-pVTZ) has
`error` set and `rhf`, `fullCI` and `coupledCluster` `null`; for Full CI and coupled cluster the
reason is also in `refusal`. `correlationEnergy` and `ccEnergy` exclude (T); `totalEnergy`
includes it.

### Adding to Xcode

How the backend is wired into the ΨLab Xcode target (the full step-by-step guide for the app is
[`app/SETUP.md`](app/SETUP.md); `xcode/project.yml` encodes the same settings for XcodeGen):

1. Copy or reference this repository's folders in the Xcode project (*File ▸ Add Files to …*,
   "Create folder references" or groups; do **not** tick "Copy items" if you want to keep a single
   source of truth).
2. In the app target's **Build Phases ▸ Compile Sources** add exactly these files:
   `backend/src/linalg.cpp`, `integrals.cpp`, `report.cpp`, `scf.cpp`, `fci.cpp`,
   `molecule.cpp`, `basis.cpp`, `gto_integrals.cpp`, `rhf.cpp`, `molecular.cpp`, `gemm.cpp`,
   `mo_integrals.cpp`, `ci.cpp`, `job.cpp`, `calculation.cpp`, `api.cpp`; add the folder `basis/`
   as a folder reference to *Copy Bundle Resources* (the app passes its path as `basisDirectory`). Do **not** add `backend/cli/main.cpp` or anything in
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
  Xcode), checks that the basis-set files are inside the bundle, runs a UI smoke test (launch, all
  tabs, potential-curve scans, menus, light and dark appearance; H₂O RHF, LiH Full CI, an H₂ curve,
  the refused Full CI of H₂O, a long-calculation warning, H₂O CCSD(T), the open method menu, LiH
  CCSDT and the refused CCSDTQ of H₂O/aug-cc-pVTZ, selected with launch arguments
  `-psilabMolecule`, `-psilabBasis`, `-psilabMethod`; screenshots) and publishes the **ΨLab-app**
  download.
