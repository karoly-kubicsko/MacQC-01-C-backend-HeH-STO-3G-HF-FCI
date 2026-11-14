# ΨLab — a native macOS quantum-chemistry app

**Current version: 5.1** — see [Versions](#versions).

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

For H₂, LiH and H₂O it also **optimises geometries** (redundant internal coordinates, analytic
gradients for HF, CCSD, CCSD(T) and Full CI), computes **harmonic frequencies, IR intensities and
thermochemistry** (analytic RHF Hessian; numerical Hessians from analytic gradients or energies
for the correlated methods) and runs the job types *energy*, *optimisation*, *frequencies*,
*optimisation + frequencies* and *optimisation + frequencies + single-point energy*, each with a
complete output file that is written while the job runs.

It shows energies, orbital energies, SCF, coupled-cluster and Davidson iterations, the largest
cluster amplitudes, the leading CI determinants and natural-orbital occupations, every
optimisation step with a drawing of the molecule (played step by step), every normal mode as an
animation, thermochemical tables and the output file, and scans
potential-energy curves with the minima marked. SwiftUI with Liquid Glass on macOS 26; runs on macOS 14 and later, on Apple
Silicon and Intel Macs. The same calculations run in the terminal on Linux and macOS
([`psilab`](#terminal-version-linux-and-macos)).

![App icon](app/icon/AppIcon-1024.png)

## Using the app

- **Download**: the repository's **Releases** page has every version, 1.0 to the current one
  (`PsiLab-<version>.app.zip`, built by CI from `main`); the newest build of a branch or pull
  request is under *Actions* ▸ *app* ▸ the run ▸ *Artifacts* ▸ **ΨLab-app** (kept 90 days; it
  contains `ΨLab.app.zip`). Unzip it; the app is not notarized: right-click `ΨLab.app` ▸ *Open*
  once, or run `xattr -dr com.apple.quarantine` on it.
- **Sidebar**: the **molecule** (HeH⁺, H₂, LiH, H₂O) and, for H₂, LiH and H₂O:
  - the **job** (a menu, in this order): 1 *Energy*, 2 *Geometry optimisation*, 3 *Frequencies*
    (at the input geometry; a warning if it is not a stationary point), 4 *Optimisation +
    frequencies* (frequencies at the optimised geometry), 5 *Optimisation + frequencies +
    single-point energy* (the single point with its own method and basis set, by default the
    same; the final energies add the thermal corrections). Jobs 4 and 5 stop with a message if the
    optimisation does not converge;
  - the geometry: bond length in Å and in bohr (both editable) and, for H₂O, the H–O–H angle —
    the *starting* geometry of an optimisation; *Experimental geometry* restores 0.74144 Å (H₂),
    1.5949 Å (LiH), 0.9572 Å and 104.52° (H₂O); after an optimisation the optimised values are
    shown with **Use optimised geometry**, which copies them into the inputs;
  - the **method** — HF, CCSD, CCSD(T), CCSDT, CCSDTQ or Full CI (a menu, in this order) — and
    the **basis set** (aug-cc-pVDZ … aug-cc-pV5Z, aug-cc-pCVDZ … aug-cc-pCV5Z);
  - job 5: the method and basis set of the single point; frequency jobs: temperature (K) and
    pressure (Pa; buttons for 1 atm and 1 bar); optimisations: normal or tight convergence and
    the maximum number of steps;
  - *Size of the calculation*, updated as you type: basis functions per atom, two-electron
    integrals (kept in memory or integral-direct), coupled-cluster amplitudes and memory, Full CI
    determinants and memory, and for jobs each step — how it is computed (analytic, or numerical
    with the number of energies or gradients, in orange) and its time — and the estimated time of
    the whole job.
    Notes (e.g. H takes aug-cc-pVXZ in the core-valence sets) in grey, warnings in orange, and a
    red explanation when the calculation is not possible (then **Run** is disabled); the run bar and,
    before the first result, the detail area repeat the refusal or the long-calculation warning.

  For HeH⁺: R in Å and bohr, Slater exponents ζ(He) and ζ(H), STO-1G/2G/3G, the exact-erf option;
  *Reset to Szabo–Ostlund defaults* restores R = 1.4632 bohr, ζ = 2.0925 / 1.24, STO-3G. Below:
  the appearance (System / Light / Dark).
- **Run** (⌘R) calculates on all processor cores. A job estimated to take longer than 30 s asks
  for confirmation first; while it runs, the sidebar shows one progress line (for example
  "Optimisation step 4 · E = -76.0418435426 E_h · max force 1.10e-04", "Hessian: displacement 3
  of 6 (S2 +)", "Lambda iteration 7 · max change 3.1e-06", SCF, coupled-cluster and Davidson
  iterations) and the elapsed time, and **Cancel** (⌘.) stops it (the output file says so).
- **Tabs** (six, so that they fit next to the sidebar at the minimum window width of 1000 pt):
  - *Summary*: energies, the job (final energy, optimised parameters, wavenumbers, ZPE, H, G, S),
    geometry, basis set per atom, convergence, timings;
  - *Methods*, with the pages *Hartree–Fock* (energy components and virial ratio, orbital
    energies, SCF iterations; for HeH⁺ also the integrals and the S/H/F/C/P matrices), *Coupled
    cluster* (energies, (T) correction, T1 diagnostic, amplitude counts and norms, largest
    amplitudes, convergence) and *Full CI* (energies, ⟨S²⟩, leading determinants, natural-orbital
    occupations, Davidson iterations; for HeH⁺ the CI vector, singlet weights and all states);
  - *Optimisation*: the **molecule of each step**, drawn in its plane to the same scale for all
    steps, with the bond lengths in Å and, for H₂O, the H–O–H angle (the starting geometry dashed),
    next to the step's energy, ΔE, maximum force and geometry; the buttons and the slider go to the
    first, previous, next or last step, **Play steps** shows the steps one after the other (0.7 s
    each, the last held 1.6 s, then from the start); then the start and optimised parameters with
    their change, whether H₂O kept C₂ᵥ, the gradient method, the steps table (energy, ΔE, maximum
    and RMS force and step, converged; a click on a step number shows that step), the energy and
    maximum-force charts (logarithmic, with the threshold; the shown step is marked) and, for the
    shown step, its Cartesian geometry, internal coordinates (value, force, step) and convergence
    table; **Use optimised geometry as input**;
  - *Frequencies*, with the pages *Frequencies* (**every normal mode as an animation**: the atoms
    vibrate about the equilibrium geometry (dashed) along the mode, at a speed proportional to its
    wavenumber, with *Play*/*Pause*, the amplitude (at 1 the largest atomic displacement is a
    quarter of the shortest bond) and arrows to the turning points; the table of modes with
    symmetry, wavenumber (imaginary in red), reduced mass, force constant and IR intensity; the IR
    stick spectrum, the selected mode in orange; how the Hessian was computed — analytic, or
    numerical with the displaced coordinates and the number of gradients or energies —, isotopes)
    and
    *Thermochemistry* (ZPE and the thermal corrections to U, H, A, G in Eₕ, kcal/mol and kJ/mol;
    E + ZPE, U, H, A, G; E, C_V and S by contribution in kcal/mol, cal/(mol K) and J/(mol K);
    moments of inertia, rotational constants and temperatures, vibrational temperatures);
  - *Curve*: the potential-curve scan (single-point energies whatever the job; a selectable step
    in Å — H₂O as a symmetric stretch — with the RHF minimum and that of the selected correlated
    method; long scans ask for confirmation and can be stopped);
  - *Output*: the output file, streamed while the job runs (it follows the end while *Follow
    output* is on), selectable, with a search field (Return: next match; arrows: next and
    previous); for HeH⁺ the FORTRAN-style report.

  *View* switches tabs and pages: Summary ⌘1, Hartree–Fock ⌘2, Coupled Cluster ⌘3, Full CI ⌘4,
  Optimisation ⌘5, Frequencies ⌘6, Thermochemistry ⌘7, Potential Curve ⌘8, Output ⌘9.
- **Menus**: *File ▸ Export Results as Text…* (⇧⌘E) saves, for H₂, LiH and H₂O, the output file
  of the job — exactly the text of the *Output* tab and of `psilab` — and for HeH⁺ the summary, the
  FORTRAN-style report and, after a scan, the potential-curve table; *File ▸ Export Results as
  JSON…* (⌥⇧⌘E) saves the results as JSON; *File ▸ Export Potential Curve as Text…* saves the scan
  table with its minima; *Edit ▸ Copy Summary as Text* (⇧⌘C) copies the short summary.
- **Build from source**: see [`app/SETUP.md`](app/SETUP.md) (Xcode 26, or XcodeGen with
  `xcode/project.yml`).

## Versions

### 5.1
The molecule of each geometry-optimisation step is drawn in the *Optimisation* tab (bond lengths,
the H–O–H angle, the starting geometry dashed), with buttons, a slider and **Play steps**; every
normal mode is animated in the *Frequencies* tab (speed proportional to the wavenumber, adjustable
amplitude, arrows). The calculations and their results are those of 5.0; the output file and
`psilab --version` say 5.1.

### 5.0
Geometry optimisation, harmonic frequencies, IR intensities and thermochemistry, five job types,
analytic gradients (RHF, CCSD, CCSD(T), Full CI) and the analytic RHF Hessian, and a complete
output file written while the job runs — see
[Version 5.0 (next version)](#version-50-next-version) at the end of this file.

### 4.0
- **Coupled cluster**: CCSD, CCSD(T), CCSDT and CCSDTQ for H₂, LiH and H₂O, all electrons
  correlated, on the RHF orbitals; the method menu offers HF, CCSD, CCSD(T), CCSDT, CCSDTQ and
  Full CI. CCSD and CCSD(T): closed-shell, spin-adapted equations with the (vv|vv) ladder over AO
  integrals (H₂O CCSD(T) up to aug-cc-pCV5Z); CCSDT and CCSDTQ: a general-order spin-orbital
  engine whose 183 diagrams are generated and checked by `backend/tools/cc_diagrams.py`.
- Verified: CCSD and CCSD(T) agree with PySCF within 3e-10 Eₕ (H₂, LiH and H₂O, seven cases);
  CCSD, CCSDT and CCSDTQ agree with brute-force coupled cluster in the determinant space to 1e-12;
  CCSDTQ reproduces Full CI for LiH (four electrons) and every method reproduces Full CI for H₂.
- Estimates of amplitudes, memory and time; refusals above the memory limit (60 % of the RAM in
  the app) or three days; a *CC* tab; potential curves with any method; Accelerate in the app.
- Terminal version: `method = HF | CCSD | CCSD(T) | CCSDT | CCSDTQ | FCI`. The energies of all
  methods for H₂, LiH and H₂O/aug-cc-pVDZ are listed under
  [Coupled cluster (version 4.0)](#coupled-cluster-version-40).

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

## Terminal version (Linux and macOS)

The same calculation runs without the graphical interface as the terminal program `psilab`, on
Linux and macOS (Apple Silicon and Intel). It reads a plain-text input file and writes a
plain-text output file.

Build and install once (needs a C++17 compiler and CMake; on a Mac: `xcode-select --install` for
Apple's compiler and CMake, e.g. `brew install cmake`):

```bash
cmake -S . -B build && cmake --build build
sudo cmake --install build                       # installs /usr/local/bin/psilab
```

On macOS the matrix products go through Apple's Accelerate, as in the app, which is much faster
for coupled cluster, Full CI and the integral transformations than the portable blocked GEMM
(about 5 GFlop/s); `-DMACQC_USE_ACCELERATE=OFF` keeps the portable code. The output file names the
code path in its header (*Matrices*). Both builds are tested in CI (Linux and macOS).

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
job      = energy            # energy | opt | freq | opt+freq | opt+freq+sp (version 5.0)
sp_method = CCSD(T)          # opt+freq+sp: method of the final single point (default: method)
sp_basis = aug-cc-pVTZ       # opt+freq+sp: its basis set (default: basis)
temperature = 298.15 K       # frequency jobs: thermochemistry temperature
pressure = 101325 Pa         # and pressure (Pa, kPa, MPa, bar, atm or Torr)
opt_convergence = normal     # normal | tight (optimisation criteria, see below)
max_steps = 50               # optimisation steps
```

The job keys and the output file are described under
[Job types and the output file](#job-types-and-the-output-file-version-50). For H₂, LiH and H₂O
`psilab` writes the output file while the job runs and flushes it after every line block, so
`tail -f water.out` follows SCF, coupled-cluster, Lambda, CPHF and optimisation iterations as they
happen; the terminal shows the final numbers (optimised energy and number of steps, wavenumbers,
G, H and ZPE, the final energy) when the job ends.

`psilab --estimate input.inp` prints the size of the calculation (basis functions, two-electron
integrals, memory) and a rough time estimate without computing; a calculation estimated to take
longer than 30 s prints a warning first. Calculations that cannot be done are refused with an
explanation (for example LiH in aug-cc-pCV5Z, which is not defined for Li, or Full CI for H₂O,
which would need 5.6·10¹¹ determinants already in aug-cc-pVDZ). The basis-set files are
found in `basis_dir` (input key), `$PSILAB_BASIS_DIR`, the installed copy
(`<prefix>/share/psilab/basis`, put there by `cmake --install`) or the source tree.
The output file contains the same text as the app's *File ▸ Export Results as Text*: for HeH⁺ the
summary, the full Hartree–Fock and Full CI report and, with `scan = yes`, the potential-curve table
(R in bohr and Å) with the HF and FCI minima; for H₂, LiH and H₂O the output file of the job
(version 5.0, see below). Mistakes in the input are reported with their line number;
the exit status is 0 on success, 1 if the calculation failed (e.g. no SCF convergence, an
optimisation that did not converge or a refused job) and 2 for input or file errors.

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

### Coupled cluster (version 4.0)

ΨLab 4.0 adds CCSD, CCSD(T), CCSDT and CCSDTQ for H₂, LiH and H₂O on the canonical RHF orbitals,
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

### Analytic derivatives (version 5.0)

ΨLab 5.0 differentiates the energy with respect to the nuclear coordinates. The first stage is
the integral machinery and the analytic RHF gradient and Hessian (`gto_derivatives.cpp`,
`rhf_derivatives.cpp`).

- **Derivative integrals**: the derivative of a contracted spherical Gaussian with respect to its
  centre, d/dA_k [x^c exp(−a r²)] = 2a G(c + 1_k) − c_k G(c − 1_k), is written as "component
  shells" (angular momentum l + 1 with the contraction coefficients times 2a, and l − 1), each with
  an output transform whose rows are the 3(2l+1) derivatives of the spherical functions; second
  derivatives use l + 2, l and l − 2 components with 6(2l+1) rows. The unchanged Obara–Saika /
  Head-Gordon–Pople engines evaluate them (shells up to l = 7 for second derivatives of h
  functions). Translational invariance gives the derivative with respect to the last centre:
  B for overlap and kinetic energy, the nucleus for nuclear attraction, D for electron repulsion
  (second derivatives from the AA, BB, CC, AB, AC, BC blocks).
- **RHF gradient**: dE/dx = Σ P h^x + Σ Γ (μν|λσ)^x − Σ W S^x + dV_nn/dx with the
  energy-weighted density W = 2 Σ_i ε_i C_i C_iᵀ, over the unique shell quartets (8-fold symmetry),
  multithreaded. When near-linearly dependent combinations are removed (canonical
  orthogonalisation), the kept orbital space rotates with the geometry; ΨLab adds the exact term
  2 Σ_{l removed, k kept} (e_lᵀ S^x e_k)(e_lᵀ F P S e_k) / (s_k (s_k − s_l)) (e, s: eigenvectors and
  eigenvalues of S), so the gradient is the derivative of the energy that is actually computed.
  PySCF omits this term.
- **RHF Hessian**: second-derivative integrals contracted with P, Γ and W, the explicit Fock
  derivatives F^(x) = h^x + G^x(P), and the coupled-perturbed HF equations for all 3N coordinates
  at once (Krylov subspace solver, orbital-energy preconditioner, one Fock build per coordinate and
  iteration). The same response gives the **dipole derivatives** (atomic polar tensors) for IR
  intensities. Without the second-order term of the orbital-space rotation, the analytic Hessian
  is exact only when no functions are removed; with removed functions, ΨLab differentiates the
  analytic gradients numerically instead (see the frequency section).

Verification (`backend/tests/test_gradients.cpp`, run by CTest; PySCF values from
`reference/pyscf_gradient_reference.py` → `backend/tests/data/gradient_reference.json`):

- Integral derivatives against 4-point central finite differences of the integrals, for a test
  basis with s … h shells on O, H and Li (62 functions): first derivatives of S, h = T + V
  (including the operator derivative), the dipole integrals and Σ Γ (μν|λσ) to 1e-9 relative
  (matrices of dS/dx, dh/dx and G^x(D) element by element); second derivatives to 1e-9 (S),
  1e-8 (h, rounding floor of the finite differences) and 3e-9 (ERI) relative.
- RHF gradient against finite differences of the energy, at distorted geometries without
  symmetry: H₂O (aug-cc-pVDZ, aug-cc-pCVDZ), LiH (aug-cc-pCVDZ), H₂ (aug-cc-pVTZ), all within
  2e-10 Eₕ/bohr. With two removed combinations (H₂O/aug-cc-pVDZ, threshold raised to 1e-2) the
  gradient including the rotation term agrees within 1.4e-10, while the term itself is
  3.2e-4 Eₕ/bohr.
- Against PySCF 2.15 (analytic gradient, analytic CPHF Hessian, dipole moment): energies within
  1e-9 Eₕ, gradients within 1.2e-10 Eₕ/bohr (H₂O aug-cc-pVDZ … aug-cc-pV5Z, LiH aug-cc-pVDZ and
  aug-cc-pCVTZ, H₂ aug-cc-pVTZ and aug-cc-pV5Z; with removed functions, ΨLab's gradient minus the
  rotation term), Hessians within 6e-10 Eₕ/bohr² (H₂O aug-cc-pVDZ and aug-cc-pVTZ, LiH
  aug-cc-pVDZ, H₂ aug-cc-pVTZ), dipole moments within 2.5e-9 e·bohr.
- RHF Hessian against finite differences of the analytic gradient (H₂O/aug-cc-pVDZ, distorted):
  1e-10 Eₕ/bohr²; dipole derivatives against finite differences of the dipole moment: 5e-10;
  translational sum rules of the Hessian and of the dipole derivatives: 1e-11.

Timings with 4 threads (2.1 GHz x86-64): H₂O/aug-cc-pVTZ gradient 1.2 s, Hessian 7 s;
H₂O/aug-cc-pVQZ gradient 10 s, Hessian 73 s (61 s derivative integrals, 12 s CPHF).

#### Correlated gradients (stage 2)

The gradients of CCSD, CCSD(T) and Full CI are analytic; CCSDT and CCSDTQ are differentiated
numerically (`gradient.cpp`, `correlated_gradient.cpp`, `cc_gradient.cpp`, `ci.cpp`).

- **Relaxed densities** (`correlated_gradient.cpp`): a method supplies its one- and
  two-particle response densities D and G in the RHF orbitals (E = Σ D h + Σ G (pq|rs) + V_nn,
  every Fock-matrix dependence written through h and (pq|rs)). ΨLab forms the orbital Lagrangian
  X_tp = 2 Σ h_tq D_qp + 4 Σ (tq|rs) G_pqrs, solves one Z-vector (CPHF) equation for the
  occupied–virtual orbital response, back-transforms G to the AO basis (stored once per
  8-fold-symmetric index class) and contracts everything with the derivative integrals of stage 1.
  The same relaxed density gives the relaxed dipole moment. The methods are invariant to
  occupied–occupied and virtual–virtual rotations, which the program checks (X_pq − X_qp below
  1e-11 in those blocks).
- **CCSD and CCSD(T)** (`cc_gradient.cpp`): closed-shell Λ equations (Stanton–Gauss
  intermediates in PySCF's formulation, DIIS), with the (vv|vv) ladder through the same
  AO-driven operator as the amplitude equations, so no (vv|vv) integrals are stored. The
  spin-traced densities are accumulated block by block into the packed two-particle density
  (the vvvv block slice by slice). CCSD(T): the Λ equations of the CCSD(T) Lagrangian get the (T)
  source terms, and the densities the (T) contributions, with triples denominators from the full
  Fock matrix (so the Lagrangian stays rotation-invariant); these terms are evaluated per ordered
  triple of occupied orbitals with v³ blocks (no o³v³ arrays).
- **Full CI** (`ci.cpp`): the two-particle density from the CI vector through (N−1)-electron
  intermediates A_xy(K) = ⟨K| a_y,β a_x,α |Ψ⟩ (opposite spins) and (N−2)-electron ones (same spin),
  each one matrix product. Full CI is invariant to all orbital rotations, so the Z-vector
  vanishes (checked).
- **Numerical gradients**: CCSDT and CCSDTQ, and correlated methods when near-linear
  dependencies were removed or the analytic gradient would not fit into the memory limit, use
  central differences of energies (step 0.005 bohr / rad) along internal coordinates: R for H₂
  and LiH (2 energies); for H₂O the symmetric stretch and the bend when the molecule is C₂ᵥ
  (4 energies, the antisymmetric derivative vanishes by symmetry), otherwise both O–H bonds and
  the angle (6 energies). The result says that it is numerical and why. The SCF of every displaced
  geometry starts from the orbitals of the reference geometry.

Verification (`backend/tests/test_correlated_gradients.cpp`; PySCF values from
`reference/pyscf_correlated_gradient_reference.py` →
`backend/tests/data/correlated_gradient_reference.json`):

- Against 4-point finite differences of the energy: CCSD and CCSD(T) for H₂O/aug-cc-pVDZ without
  symmetry within 1.7e-10 and 8.6e-11 Eₕ/bohr; LiH/aug-cc-pVDZ CCSD, CCSD(T) and Full CI within
  6e-11; H₂/aug-cc-pVTZ Full CI within 4e-10; in the *full* CI job also H₂O/aug-cc-pVTZ CCSD(T)
  within 2.7e-10 and LiH/aug-cc-pCVDZ Full CI within 1e-10.
- Against PySCF 2.15 (`pyscf.grad.ccsd`, `pyscf.grad.ccsd_t`, CASCI with all orbitals active):
  energies within 5e-12 Eₕ, gradients within 7.5e-10 Eₕ/bohr (CCSD, CCSD(T) of H₂O and LiH; CCSD of
  H₂) and 4.9e-9 (Full CI of LiH). For H₂/aug-cc-pVTZ, PySCF's CASCI gradient differs from PySCF's
  own CCSD gradient (exact for two electrons, energies equal to 1e-12) by 2.5e-7; ΨLab's Full CI
  gradient agrees with the CCSD one within 3e-10.
- The densities reproduce the energies within 3e-12 Eₕ; the blocked (T) terms agree with the
  full-array formulation within 2e-18; numerical gradients agree with the analytic ones within
  1e-5 Eₕ/bohr; the numerical CCSDT gradient of H₂ equals the analytic Full CI one within 8e-6.

Timings with 4 threads: H₂O CCSD(T) gradient in aug-cc-pVDZ 1.5 s, aug-cc-pVTZ 24 s; LiH Full CI
gradient in aug-cc-pVTZ 89 s.

### Geometry optimisation (version 5.0)

`internal_coordinates.cpp`, `optimizer.cpp`: minimisation in redundant internal coordinates
(Pulay & Fogarasi, J. Chem. Phys. 96, 2856 (1992); Peng, Ayala, Schlegel & Frisch, J. Comput.
Chem. 17, 49 (1996)).

- **Coordinates**: bonds between atoms closer than 1.3 × the sum of their covalent radii (separate
  fragments are joined through their closest atoms), bends at every atom with two bonds, and two
  linear-bend components in perpendicular planes when an angle exceeds 175° (back to a bend below
  170°). Wilson B matrix (analytic for stretches and bends), G = BBᵀ with its generalised inverse,
  internal gradient g_q = G⁻ B g_x, projector P = G G⁻.
- **Step**: rational-function optimisation (lowest eigenvector of the augmented Hessian of the
  projected gradient and Hessian P H P + 1000 (1 − P)), scaled to the trust radius (0.3 at the
  start, adapted between 0.001 and 0.5 from the ratio of the actual to the predicted energy
  change); iterative back-transformation x ← x + Bᵀ G⁻ (q_target − q(x)) to 1e-10 bohr.
- **Hessian**: initial Hessian analytic (the RHF Hessian of stage 1 in internal coordinates,
  including the gradient term Σ g_q ∂²q/∂x², when the method is RHF and the basis has at most 100
  functions), otherwise Lindh's model Hessian (Chem. Phys. Lett. 241, 423 (1995)); Schlegel's is
  available. Negative eigenvalues of an analytic initial Hessian are replaced by their magnitude;
  if the start is a saddle point with no gradient along the negative-curvature direction (e.g. H₂O
  at exactly 180°), the first step follows that direction. BFGS updates (skipped, and reported,
  when the curvature condition fails).
- **Convergence** (as in Gaussian): maximum and RMS internal force (normal 4.5e-4 and 3.0e-4
  Eₕ/bohr, tight 1.5e-5 and 1.0e-5) and maximum and RMS displacement of the next step (normal
  1.8e-3 and 1.2e-3 bohr, tight 6.0e-5 and 4.0e-5), plus the energy change since the previous
  step (normal 1e-6, tight 1e-8 Eₕ); converged when all are met, or when both forces are 100 times
  below their thresholds. Each step reports which criteria are met. Default at most 50 steps; an
  optimisation that does not converge stops with a message, and no frequency or single-point step
  follows (jobs 4 and 5).
- **Symmetry**: H₂O is optimised in all three internal coordinates without constraints; the report
  says whether the C₂ᵥ start stays C₂ᵥ (|R(O–H1) − R(O–H2)|, 1e-14 bohr in practice). A change of the
  number of removed near-linear dependencies along the path is reported (the energy surface is
  then not continuous).

Verification (`backend/tests/test_optimizer.cpp`; references from
`reference/pyscf_optimization_reference.py` → `backend/tests/data/optimization_reference.json`):

- Tight optimisations from the default geometries against PySCF analytic gradients minimised by
  geomeTRIC (an independent optimiser): H₂O RHF/aug-cc-pVDZ and aug-cc-pVTZ, CCSD/aug-cc-pVDZ,
  CCSD(T)/aug-cc-pVDZ; LiH RHF, CCSD(T) and Full CI/aug-cc-pVDZ; H₂ Full CI/aug-cc-pVTZ. Bond
  lengths agree within 8e-6 Å, angles within 6e-4°, energies within 1e-10 Eₕ. (PySCF's own
  gradient scanner for CCSD(T) passes the CCSD Λ multipliers to the (T) gradient and returns the
  CCSD energy; the reference script uses a scanner with the CCSD(T) Λ equations instead — with
  PySCF's scanner the H₂O CCSD(T) optimum is 1.9e-3 Å and 0.12° off and 6.7e-6 Eₕ higher.)
- Difficult starts reach the same RHF minimum: an asymmetric start (1.05 Å, 0.90 Å, 120°) ends
  C₂ᵥ; 179°; exactly 180° (a saddle point, left along the bend, C₂ᵥ kept); 0.8 Å and 70°.
- Wilson B matrices (also linear bends) against finite differences to 1e-10; back-transformation
  residual below 1e-15; a CCSDT optimisation of H₂ with numerical gradients reaches the analytic
  Full CI optimum within 7e-6 Å; the step limit stops with a message.

### Harmonic frequencies and thermochemistry (version 5.0)

`hessian.cpp`, `frequencies.cpp`, `constants.cpp`.

- **Hessian**: RHF analytic (stage 1: second-derivative integrals and CPHF, with analytic dipole
  derivatives). CCSD, CCSD(T) and Full CI: central differences of analytic gradients along internal
  coordinates (step 0.005 bohr / rad); for C₂ᵥ H₂O the symmetry-adapted S1 (A1, symmetric stretch),
  S2 (A1, bend) and S3 (B2, antisymmetric stretch), 6 gradients; for H₂ and LiH the bond, 2
  gradients. The Cartesian Hessian is Bᵀ H_q B + Σ g_q ∂²q/∂x², exact also away from a stationary
  point. The displaced geometries are superimposed on the reference (mass-weighted, Horn's
  quaternion method, i.e. the Eckart frame), so the relaxed dipole moments of the gradients give
  the dipole derivatives. CCSDT and CCSDTQ: second differences of energies along the same
  coordinates (mixed A1 × B2 terms vanish by symmetry; C₂ᵥ H₂O 9 energies, diatomics 2), dipole
  moments by finite fields μ = −dE/dF at the displaced geometries (10 and 4 more energies). With
  removed near-linear dependencies, RHF also differentiates its analytic gradients numerically.
  The output and the app say which was used. The SCF accepts a uniform electric field for this
  (and for the tests of the relaxed dipoles).
- **Normal modes**: masses of the most abundant isotopes (AME 2020: ¹H 1.007825031898 u,
  ⁷Li 7.016003434 u, ¹⁶O 15.994914619257 u); mass-weighted Hessian in the orthonormal complement of
  the translations and of the rotations about the principal axes through the centre of mass
  (Eckart conditions; 5 for linear molecules); harmonic wavenumbers (negative = imaginary),
  normalised Cartesian displacements, reduced masses, force constants (mdyn/Å), IR intensities
  (km/mol) from dμ/dQ, symmetry species for C₂ᵥ H₂O (A1, B2) and diatomics; the residual
  curvature of the translations and rotations is printed as a check.
- **Thermochemistry** (CODATA 2018): ideal-gas translation at the chosen pressure (Sackur–Tetrode),
  rigid rotor (linear for H₂ and LiH, nonlinear for H₂O; σ(H₂) = 2, σ(LiH) = 1, σ(H₂O) = 2), harmonic
  oscillators with the zero-point energy, closed-shell singlet (q_el = 1). Reported: zero-point
  energy; thermal corrections and totals of U, H = U + k_BT, S (total and its electronic,
  translational, rotational and vibrational parts), A = U − TS, G = H − TS; C_V; moments of
  inertia, rotational constants (GHz, cm⁻¹), rotational and vibrational temperatures; in hartree,
  kcal/mol, kJ/mol and cal/(mol K) (also J/(mol K)). Imaginary frequencies and those below 10 cm⁻¹
  are left out of the vibrational sums, with a warning.

Verification (`backend/tests/test_frequencies.cpp`; `reference/pyscf_frequency_reference.py` →
`backend/tests/data/frequency_reference.json`):

- RHF wavenumbers, reduced masses and force constants against PySCF's analytic Hessian and
  harmonic analysis (given the same AME 2020 masses; by default PySCF uses standard atomic
  weights such as 1.008 for H, which shifts the wavenumbers of H₂O by 0.4 cm⁻¹) at the
  PySCF/geomeTRIC geometries of H₂O (aug-cc-pVDZ, aug-cc-pVTZ), LiH and H₂: within 2e-3 cm⁻¹,
  1e-6 u and 1e-5 mdyn/Å.
- Thermochemistry against a hand evaluation of the formulas with CODATA 2018 (in the reference
  script, independent of both programs): ZPE, U, H, G, A within 1e-9 Eₕ, entropies and C_V within
  5e-5 J/(mol K); against PySCF's `thermo` (older CODATA constants) within 3e-6 Eₕ and
  5e-3 J/(mol K).
- Numerical Hessians: finite differences of analytic RHF gradients reproduce the analytic
  wavenumbers within 0.05 cm⁻¹ and the IR intensities within 6e-4 km/mol; for CCSD(T) (H₂O) and
  Full CI (LiH) the gradient-based and the energy-based Hessians agree within 0.03 cm⁻¹, and the IR
  intensities from relaxed dipoles and from finite-field dipoles within 0.007 km/mol; relaxed
  dipole moments equal −dE/dF within 1e-9 e·bohr.
- At linear H₂O the degenerate bend is imaginary (twice) and left out of the thermochemistry with
  two warnings.

### Job types and the output file (version 5.0)

`jobs.cpp`, `output.cpp`. Five job types for H₂, LiH and H₂O (HeH⁺ keeps its 4.0 energy path):

| Job | What runs | Geometry of each step |
|---|---|---|
| 1 `energy` | single-point energy, as in 4.0 | input |
| 2 `opt` | geometry optimisation | input → optimised |
| 3 `freq` | Hessian, normal modes, IR intensities, thermochemistry; a warning if the input geometry is not stationary (largest gradient component above 4.5e-4 Eₕ/bohr) | input |
| 4 `opt+freq` | optimisation, then frequencies | frequencies at the optimised geometry |
| 5 `opt+freq+sp` | optimisation and frequencies with `method`/`basis`, then a single-point energy with `sp_method`/`sp_basis`; final U, H, G, A = single-point energy + thermal corrections of step 2 | optimised |

If the optimisation of job 4 or 5 does not converge within `max_steps`, the job stops after it
with an error termination that says why ("the frequency and single-point steps were not run:
they need a stationary point"). Every method can be used where its energy is possible: analytic
gradients for RHF, CCSD, CCSD(T) and Full CI; numerical gradients (central differences of
energies along internal coordinates) for CCSDT and CCSDTQ; analytic RHF Hessians; numerical
Hessians from analytic gradients (CCSD, CCSD(T), Full CI) or from energies (CCSDT, CCSDTQ). The
estimate lists each step with how it is computed, the number of energies, gradients or steps,
the time and the memory; jobs that would take more than three days are refused, and so are
CCSDTQ optimisations and frequencies of H₂O (no analytic CCSDTQ gradient: each step would need
CCSDTQ energies of hours each and more than 5 GB; the CCSDTQ energy itself remains possible with
enough memory, also as the single point of job 5). A job whose single point is impossible is
refused before anything runs.

The **output file** is written while the job runs and is identical in the app's *Output* view, in
*File ▸ Export Results as Text* and in `psilab`'s output file. In order:

1. header: program and version, start time, job type, threads; the input (keys as understood, for
   `psilab` also the input file as given);
2. molecule: electrons, charge, Cartesian coordinates in Å and bohr, bond lengths and angle,
   nuclear repulsion; basis set per element (primitives → contraction, functions), number of
   functions and shells, storage of the two-electron integrals;
3. estimates: size of Full CI or coupled cluster, each step of the job, total time, warnings or the
   refusal;
4. every energy calculation as it happens: SCF iterations (energy, change, orbital gradient,
   DIIS); CC iterations (energy, correlation energy, change, largest amplitude change) or Davidson
   iterations; Lambda iterations; CPHF / Z-vector iterations;
5. at the energy job's geometry or the optimised geometry: energy components (kinetic, nuclear
   attraction, electron repulsion, virial ratio), all orbital energies with HOMO, LUMO and gap,
   Mulliken and Löwdin charges, dipole moment (SCF density; relaxed density for correlated
   gradients), the coupled-cluster or Full CI report (T1 diagnostic, leading amplitudes or
   determinants, natural occupations);
6. optimisation, every step: the Cartesian geometry it uses, the energy and its change, the
   Cartesian gradient and how it was computed, the table of internal coordinates (name,
   definition, value in bohr/rad and Å/deg, force, step), the Hessian (initial or update), the
   RFO step (λ, length, trust radius, predicted change, ratio of actual to predicted change), and
   the convergence table (item, value, threshold, converged?); then "Stationary point found", a
   step summary, the optimised geometry and parameters (start in parentheses) and whether H₂O kept
   C₂ᵥ;
7. frequencies: the geometry used, how the Hessian and the dipole derivatives were obtained (and
   each displaced calculation of a numerical Hessian), isotopes and masses, the residual curvature
   of translations and rotations, the table of modes (symmetry, wavenumber, reduced mass, force
   constant, IR intensity) and the normalised Cartesian displacements;
8. thermochemistry: T, p, mass, symmetry number, moments of inertia, rotational constants and
   temperatures, vibrational temperatures; ZPE and the thermal corrections to U, H, A, G in Eₕ,
   kcal/mol and kJ/mol; E + ZPE, U, H, A, G; E, C_V and S of the electronic, translational,
   rotational and vibrational parts in kcal/mol, cal/(mol K) and J/(mol K), with partition
   functions;
9. single point of job 5 with its own SCF and correlated iterations;
10. summary (final geometry, final energy and its label, convergence, wavenumbers, thermochemistry)
    and "Normal termination" with the total time, or "Error termination" with the reason
    (refusal, non-converged optimisation, failed SCF) or "Calculation cancelled".

Warnings (non-stationary geometry, removed near-linear dependencies, imaginary frequencies left out
of the thermochemistry, rises of the energy during the optimisation) appear where they occur and
again in the summary.

Verification (`backend/tests/test_jobs.cpp`, `backend/tests/test_psilab.cpp`):

- every job type with RHF for H₂ (aug-cc-pVTZ), LiH and H₂O (aug-cc-pVDZ), tight convergence,
  through the C API: optimised bond lengths and angle against PySCF + geomeTRIC within
  2e-5 Å and 2e-3°, energies within 1e-9 Eₕ; frequencies against PySCF's analytic Hessian at the
  same geometry within 0.01 cm⁻¹, at ΨLab's own optimised geometry within 0.05 cm⁻¹; ZPE, S and G
  against PySCF's `thermo`;
- job 5 (H₂O, RHF/aug-cc-pVDZ geometry and frequencies, CCSD(T)/aug-cc-pVTZ single point): the
  single point equals an energy job at the optimised geometry (8e-10 Eₕ; the energy job converges
  the amplitudes to 1e-8, the job to 1e-9), and E + ZPE,
  U, H, G, A equal the single-point energy plus the thermal corrections (1e-10 Eₕ);
- correlated methods (aug-cc-pVDZ): `opt+freq` with CCSD, CCSD(T), CCSDT, CCSDTQ and Full CI for
  H₂, where CCSD = CCSDT = CCSDTQ = Full CI exactly (two electrons): R_e, the wavenumber and the
  energy agree (analytic CCSD against analytic Full CI within 1e-6 Å and 0.05 cm⁻¹; numerical
  CCSDT and CCSDTQ within 5e-5 Å and 0.5 cm⁻¹); LiH `opt+freq` with CCSD, CCSD(T) and Full CI (in the
  *full* CI job also CCSDT and CCSDTQ, about 90 s per energy, with CCSDTQ = Full CI for four
  electrons); LiH CCSD(T) and Full CI and H₂O
  CCSD and CCSD(T) geometries against PySCF + geomeTRIC;
- refusals (CCSDTQ frequencies and optimisation of H₂O, Full CI of H₂O, an impossible single
  point in job 5), the counts of the numerical Hessians in the estimates (H₂O CCSDT: 9 energies +
  10 finite-field energies; CCSD(T): 6 gradients; diatomics: 2);
- the streamed output: read with `macqc_job_output` while the job runs, it is identical to the
  final `output`; progress text; cancellation during an optimisation step (the run returns
  "calculation cancelled" and the output says so);
- `psilab`: all new keys and their errors, the streamed text equals the output file.

### Layout

```
backend/include/macqc/   public headers (macqc_api.h is the C API)
backend/src/             linalg, integrals, scf, fci, report, calculation, api (HeH+, STO-nG);
                         molecule, basis, gto_integrals, rhf, molecular (H2, LiH, H2O);
                         gemm, mo_integrals, ci (Full CI for H2, LiH); job (cancel/progress);
                         cc_closed_shell (CCSD, (T)), cc_general (CCSDT, CCSDTQ), cc_estimate;
                         gto_derivatives (derivative integrals), rhf_derivatives (RHF gradient,
                         Hessian, CPHF, dipole derivatives); tensor (einsum), cc_gradient (CCSD
                         and CCSD(T) Lambda equations and densities), correlated_gradient
                         (relaxed densities, Z-vector), gradient (energies and gradients of
                         every method at any geometry); internal_coordinates, optimizer
                         (redundant internal coordinates, RFO, BFGS); hessian (analytic and
                         numerical Hessians, dipole derivatives), frequencies (normal modes,
                         thermochemistry), constants (CODATA 2018, AME 2020)
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
      // (iteration = percent done); 5.0: 8 gradient, 9 Hessian, 10 CPHF, 11 Lambda equations,
      // 12 geometry optimisation, 13 numerical Hessian
char* macqc_job_output(const macqc_job* job, long long offset, long long* total);
      // 5.0: the output file written so far, from byte `offset` on (free with macqc_free);
      // *total = its full length, so a reader polls with offset = the previous total
char* macqc_job_progress_text(const macqc_job* job);
      // 5.0: one line about the current step (free with macqc_free), e.g.
      // "Optimisation step 4 · E = -76.0418435426 E_h · max force 1.10e-04",
      // "Hessian: displacement 3 of 6 (S2 -)", "Orbital response iteration 7 · residual 3.1e-06"
void  macqc_job_free(macqc_job* job);      // after macqc_run_job has returned
```

The app polls `macqc_job_progress` four times a second for the progress line of the sidebar.
The output buffer of a job is a string guarded by a mutex inside the job handle: the calculation
appends each piece of the output file as it is produced and `macqc_job_output` copies the part a
reader has not seen yet, so the app's *Output* view follows the job without blocking it.

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

Version 5.0 adds the job keys (H₂, LiH, H₂O; HeH⁺ only runs the energy job):

```json
{"job":"opt+freq+sp","spMethod":"CCSD(T)","spBasis":"aug-cc-pVTZ",
 "temperature":298.15,"pressure":101325,"optConvergence":"normal","maxSteps":50}
```

`job` is `"energy"` (default), `"opt"`, `"freq"`, `"opt+freq"` or `"opt+freq+sp"` (other
spellings such as `"OPT FREQ"` or `"5"` are accepted); `spMethod`/`spBasis` default to `method`
and `basis`; `temperature` in K (0 < T ≤ 10⁵), `pressure` in Pa; `optConvergence` `"normal"` or
`"tight"`; `maxSteps` 1–1000. With `"estimateOnly":true` the result's `job.estimate` lists the
steps of the job (how each is computed, how many energies, gradients or steps, time, memory) and
the refusal, if any.

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

Version 5.0 adds, for H₂, LiH and H₂O (`null` where the job has no such step):

```text
  "job": {"type","title","termination":"normal"|"error"|"refused"|"estimate",
          "completed","message","seconds","warnings",
          "estimate": {"steps":[{"name","method","basis","how","analytic","evaluations",
                                 "evaluationKind":"steps"|"gradients"|"energies"|"Hessians",
                                 "seconds","memoryMB"}],
                       "totalSeconds","memoryMB","longCalculation","refusal","warnings","notes"},
          "initialGeometry","finalGeometry",           // {"atoms","R","R1","R2","angle"}, bohr, deg
          "finalEnergy","finalEnergyLabel",            // e.g. -76.33, "CCSD(T)/aug-cc-pVTZ"
          "singlePointMethod","singlePointBasis"},
  "optimization": {"converged","steps","message","initialHessian","gradient",
          "criteria":{"maxForce","rmsForce","maxDisplacement","rmsDisplacement","energy"},
          "initialEnergy","finalEnergy","initialGeometry","finalGeometry",
          "initialParameters","finalParameters":[{"name","value","angle"}],   // bohr or degrees
          "symmetricStart","symmetricEnd","symmetryDeviation","seconds",
          "history":[{"step","energy","deltaE","geometry",
                      "coordinates":[{"name","definition","value","gradient","step"}],
                      "convergence":[{"name","value","threshold","met","applicable"}],
                      "converged","hessianUpdate","trustRadius","stepNorm","stepScaled",
                      "predictedChange","ratio","warnings","seconds"}]},
  "frequencies": {"analytic","description","numericalReason","gradientEvaluations",
          "energyEvaluations","coordinates","intensities","intensityDescription","dipole",
          "linear","vibrations","masses","isotopes","geometry","energy","lowFrequencies",
          "modes":[{"wavenumber","imaginary","symmetry","reducedMass","forceConstant",
                    "irIntensity","displacement":[[x,y,z] per atom]}],"seconds"},
  "thermochemistry": {"temperature","pressure","symmetryNumber","linear","mass","moments",
          "rotationalConstantsGHz","rotationalConstantsCm","rotationalTemperatures",
          "vibrationalTemperatures","excludedModes","warnings",
          "terms":{"translation"|"rotation"|"vibration"|"electronic"|"total":
                   {"energy","entropy","cv","q"}},                 // J/mol, J/(mol K)
          "zeroPointEnergy","corrections":{"zeroPoint","thermalEnergy","enthalpy","gibbs",
          "helmholtz"},                                             // J/mol
          "electronicEnergy","totals":{...same keys...}},           // E_h
  "singlePoint": {"method","energy","converged","rhfEnergy","correlationEnergy","geometry",
          "warnings","seconds"},
  "output": "the complete output file"
```

`molecule`, `basis`, `rhf`, `fullCI`, `coupledCluster` and `report` then describe the final
energy calculation of the job (the optimised geometry; for `opt+freq+sp` the single point). In
`thermochemistry.totals` the energy is the job's final energy (for `opt+freq+sp` the single point)
plus the thermal corrections of the frequency step.

### Adding to Xcode

How the backend is wired into the ΨLab Xcode target (the full step-by-step guide for the app is
[`app/SETUP.md`](app/SETUP.md); `xcode/project.yml` encodes the same settings for XcodeGen):

1. Copy or reference this repository's folders in the Xcode project (*File ▸ Add Files to …*,
   "Create folder references" or groups; do **not** tick "Copy items" if you want to keep a single
   source of truth).
2. In the app target's **Build Phases ▸ Compile Sources** add exactly these files:
   `backend/src/linalg.cpp`, `integrals.cpp`, `report.cpp`, `scf.cpp`, `fci.cpp`,
   `molecule.cpp`, `basis.cpp`, `gto_integrals.cpp`, `rhf.cpp`, `molecular.cpp`, `gemm.cpp`,
   `mo_integrals.cpp`, `ci.cpp`, `cc_general.cpp`, `cc_closed_shell.cpp`, `cc_estimate.cpp`,
   `gto_derivatives.cpp`, `rhf_derivatives.cpp`, `tensor.cpp`, `correlated_gradient.cpp`,
   `cc_gradient.cpp`, `gradient.cpp`, `internal_coordinates.cpp`, `optimizer.cpp`, `constants.cpp`,
   `frequencies.cpp`, `hessian.cpp`, `output.cpp`, `jobs.cpp`, `job.cpp`,
   `calculation.cpp`, `api.cpp`; add the folder `basis/`
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
  the input *full*, a third job checks every RHF case up to 341 basis functions. A macOS job
  (`macos-26`, Apple clang, Accelerate) builds the backend and `psilab`, runs all tests, installs
  `psilab` and runs an optimisation + frequencies + single-point job with it, with the basis-set
  files taken from the installed copy.
- `.github/workflows/app.yml` builds the universal ΨLab app on a `macos-26` runner (XcodeGen +
  Xcode), checks that the basis-set files are inside the bundle, runs a UI smoke test (launch, all
  tabs, potential-curve scans, menus, light and dark appearance; H₂O RHF, LiH Full CI, an H₂ curve,
  the refused Full CI of H₂O, a long-calculation warning, H₂O CCSD(T), the open method menu, LiH
  CCSDT and the refused CCSDTQ of H₂O/aug-cc-pVTZ; for 5.0 the open job menu, an H₂O RHF
  optimisation, LiH CCSD(T) frequencies and thermochemistry, H₂O optimisation + frequencies with
  a CCSD(T)/aug-cc-pVTZ single point, the Output tab while a job runs, Cancel and the refused
  CCSDTQ frequencies of H₂O; for 5.1 the molecule of the optimisation steps, still and while
  *Play steps* runs, and the animated normal modes of LiH and H₂O (two captures 0.3 s apart must
  differ; a few screenshots are printed into the log as small JPEGs), selected with launch arguments `-psilabMolecule`, `-psilabBasis`,
  `-psilabMethod`, `-psilabJob`, `-psilabSPMethod`, `-psilabSPBasis`, `-psilabR` (Å),
  `-psilabAngle`, `-psilabTab`, `-psilabAutoRun`; screenshots on the runner's 1024-pt-wide screen,
  where the window has its minimum width of 1000 pt) and publishes the **ΨLab-app** download.
  Every build of `main` is also published as a GitHub Release, one per app version (tag
  `v<version>`, asset `PsiLab-<version>.app.zip`; a later build of the same version replaces it).
- `.github/workflows/releases.yml` turned the last builds of versions 1.0 to 4.0 (CI artifacts)
  into the releases `v1.0` … `v4.0`; it runs when the file changes on `main` or by hand and leaves
  existing releases alone.

## Version 5.0 (next version)

ΨLab 5.0 is the next version of the application, after 4.0. It adds geometry optimisation,
harmonic frequencies, IR intensities and thermochemistry for H₂, LiH and H₂O, five job types and a
complete output file; HeH⁺ (Szabo–Ostlund) is unchanged (its Hartree–Fock text output is still
identical to `reference/hartree_fock_sto3g.py`).

**Job types** (app: the sidebar's *Job* menu; `psilab`: `job = …`; C API: `"job"`):
1 *energy*; 2 *opt* (geometry optimisation); 3 *freq* (frequencies and thermochemistry at the
input geometry, with a warning if it is not stationary); 4 *opt+freq* (frequencies at the
optimised geometry); 5 *opt+freq+sp* (then a single point with `sp_method` / `sp_basis`, final
U, H, G, A = single-point energy + thermal corrections). Jobs 4 and 5 stop with a message if the
optimisation does not converge. Every step prints the geometry it uses.

**What is analytic and what is numerical**

| Method | Gradient | Hessian | Dipole derivatives (IR) |
|---|---|---|---|
| HF (RHF) | analytic (derivative integrals, energy-weighted density) | analytic (second-derivative integrals, CPHF) | analytic (CPHF) |
| CCSD, CCSD(T) | analytic (Λ equations, relaxed densities, Z-vector) | central differences of analytic gradients along symmetry-adapted internal coordinates (H₂O 6, diatomics 2 gradients) | relaxed dipoles of those gradients |
| Full CI | analytic (relaxed densities, Z-vector) | as CCSD(T) | as CCSD(T) |
| CCSDT, CCSDTQ | numerical (central differences of energies along internal coordinates: H₂O 4, diatomics 2 energies) | second differences of energies (H₂O 9, diatomics 2) | finite fields μ = −dE/dF (H₂O 10, LiH 4 more energies) |

Correlated methods also differentiate numerically when near-linear dependencies were removed from
the basis or the analytic gradient would not fit into the memory limit, and RHF uses numerical
second derivatives of its analytic gradients then; the output file, the app and the JSON say
which was used and why. CCSDTQ optimisations and frequencies of H₂O are refused (hours per energy
and more than 5 GB); a single point with CCSDTQ in job 5 remains possible where the energy is.

**Verification** (each stage has its own test program, run by `ctest` and in CI; references from
PySCF 2.15 and geomeTRIC by the scripts in `reference/`):

| Stage | What | Result |
|---|---|---|
| 1 | first and second derivative integrals vs finite differences | 7e-10 relative (first); S 1e-11, h 2e-9, ERI 5e-9 (second) |
| 1 | RHF gradient vs finite differences / PySCF (aug-cc-pVDZ … aug-cc-pV5Z) | 2e-10 / 1.2e-10 Eₕ/bohr |
| 1 | RHF Hessian vs PySCF / finite differences; dipole derivatives | 6e-10 / 1e-10 Eₕ/bohr²; 5e-10 e |
| 2 | CCSD, CCSD(T), Full CI gradients vs finite differences | ≤ 1.7e-10 Eₕ/bohr |
| 2 | vs PySCF (`grad.ccsd`, `grad.ccsd_t`, CASCI) | ≤ 7.5e-10 (CC), 4.9e-9 (Full CI) |
| 3 | optimised geometries vs PySCF + geomeTRIC (8 cases) | ≤ 8e-6 Å, 6e-4°, 1e-10 Eₕ |
| 4 | RHF wavenumbers vs PySCF analytic Hessian | ≤ 2e-3 cm⁻¹ |
| 4 | thermochemistry vs the formulas evaluated by hand (CODATA 2018) / PySCF `thermo` | 1e-9 Eₕ, 5e-5 J/(mol K) / 3e-6 Eₕ |
| 4 | numerical vs analytic Hessians; relaxed dipoles vs −dE/dF | ≤ 0.05 cm⁻¹; 1e-9 e bohr |
| 5 | all job types through the C API vs PySCF; H₂: CCSD = CCSDT = CCSDTQ = Full CI; LiH: CCSDTQ = Full CI (*full* job) | 2e-5 Å, 0.01–0.05 cm⁻¹; 1.5e-9 Å (analytic), 7.4e-6 Å and 0.15 cm⁻¹ (numerical); 6e-12 Eₕ and 2.4e-3 cm⁻¹ |
| 5 | streamed output = final output file; cancellation; refusals; estimates | identical; < 1 ms during an optimisation step |

**Comparison with experiment.** Harmonic constants from optimisations and frequencies with the
largest feasible calculations (all electrons correlated):

| Quantity | ΨLab | Experiment |
|---|---|---|
| H₂ rₑ, Full CI/aug-cc-pV5Z | 0.741602 Å | 0.74144 Å |
| H₂ ωₑ, Full CI/aug-cc-pV5Z | 4403.42 cm⁻¹ | 4401.21 cm⁻¹ |
| LiH rₑ, Full CI/aug-cc-pCVTZ | 1.601196 Å | 1.5949 Å |
| LiH ωₑ, Full CI/aug-cc-pCVTZ | 1393.36 cm⁻¹ | 1405.65 cm⁻¹ |
| H₂O rₑ, θₑ, CCSD(T)/aug-cc-pVTZ | 0.958852 Å, 104.342° | 0.95782 Å, 104.48° |
| H₂O ω₁ (A₁), ω₂ (A₁), ω₃ (B₂), CCSD(T)/aug-cc-pVTZ | 3838.53, 1641.38, 3942.48 cm⁻¹ | 3832.2, 1648.5, 3942.5 cm⁻¹ |
| same, CCSD(T)/aug-cc-pCVTZ | 3813.59, 1645.52, 3920.97 cm⁻¹ (rₑ 0.960573 Å, θₑ 104.289°) | |
| S°(H₂O, 298.15 K, 1 bar), CCSD(T)/aug-cc-pVTZ | 188.69 J/(mol K) | 188.83 J/(mol K) |

Sources: H₂O harmonic frequencies from the NIST Computational Chemistry Comparison and Benchmark
Database (experimental data); H₂O equilibrium structure from Császár et al., J. Chem. Phys. 122,
214305 (2005) (rₑ = 0.95782 Å, θₑ = 104.48(5)°; the default geometry of ΨLab, 0.9572 Å and
104.52°, is the older, widely quoted value); S°(H₂O) from the NIST-JANAF tables (188.834) and the
CODATA key values for thermodynamics (188.835 ± 0.010); LiH and H₂ spectroscopic constants from
Huber & Herzberg, *Constants of Diatomic Molecules* (1979). The H₂ values (Full CI in a nearly
complete basis) differ by 2 cm⁻¹ and 2e-4 Å — what the Born–Oppenheimer approximation leaves
out. For LiH and H₂O the differences are those of the basis sets (aug-cc-pCVTZ and aug-cc-pVTZ are
not complete; for H₂O also the missing core-correlating functions in aug-cc-pVTZ) and, for
CCSD(T), of the method. The rigid-rotor, harmonic-oscillator entropy of H₂O is 0.14 J/(mol K)
below the tabulated value, the size of the anharmonic and rotation–vibration contributions it
leaves out.

**Terminal version on macOS**: `psilab` builds with CMake on macOS too, with Apple's Accelerate for
the matrix products; CI builds and tests it on Linux and on macOS (all tests, plus an installed
`psilab` running job 5).

**Not in this version**: analytic gradients of CCSDT and CCSDTQ (they need general-order Λ
equations; the gradients are numerical), analytic Hessians of the correlated methods (numerical
from analytic gradients), analytic correlated gradients with removed near-linear dependencies
(numerical then), transition-state searches, anharmonic frequencies, open-shell molecules.

**Usage**: choose the job in the sidebar's **Job** menu, the method and basis set (and for job 5
the single-point method and basis set), set the starting geometry, temperature, pressure and
convergence, and press **Run**; the *Optimisation*, *Frequencies* and *Output* tabs fill as the
job runs. In the terminal: `job = opt+freq` (see [Terminal version](#terminal-version-linux-and-macos));
`tail -f` follows the output file.
