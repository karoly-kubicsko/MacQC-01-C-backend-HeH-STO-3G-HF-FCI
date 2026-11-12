#!/usr/bin/env python3
"""Independent RHF reference energies for ΨLab 3.0, computed with PySCF.

PySCF (https://pyscf.org) reads the very same Gaussian94 basis files in basis/ with its
own parser (pyscf.gto.basis.parse_gaussian); the geometries are built exactly as in
backend/src/molecule.cpp (coordinates in bohr, 1 bohr = 0.529177210544 Å). For H in the
core-valence sets aug-cc-pCVXZ (not defined for H) the aug-cc-pVXZ functions are used,
as in ΨLab. Spherical harmonics, tight convergence (energy 1e-12), no screening loss, and
overlap eigenvalues below 1e-7 removed as in ΨLab.

Usage:
  python3 reference/pyscf_rhf_reference.py                 # all cases -> backend/tests/data/rhf_reference.json
  python3 reference/pyscf_rhf_reference.py --max-nbf 120   # only cases with <= 120 functions
  python3 reference/pyscf_rhf_reference.py --check         # recompute and compare with the stored file

Needs: pip install pyscf
"""
import argparse
import json
import math
import pathlib
import sys

import pyscf
from pyscf import gto, scf
from pyscf.gto.basis import parse_gaussian

# Same treatment of near-linear dependencies as ΨLab: overlap eigenvalues below 1e-7 are
# removed (PySCF's own default threshold is 1e-6; LiH/aug-cc-pCVQZ has 8.5e-7, so with the
# default PySCF would drop one function and give an energy 3e-6 E_h higher).
LINDEP_THRESHOLD = 1e-7
scf.hf.remove_overlap_zero_eigenvalue = True
scf.hf.overlap_zero_eigenvalue_threshold = LINDEP_THRESHOLD

ROOT = pathlib.Path(__file__).resolve().parents[1]
BASIS_DIR = ROOT / "basis"
OUT = ROOT / "backend" / "tests" / "data" / "rhf_reference.json"

BOHR_IN_ANGSTROM = 0.529177210544  # CODATA 2022, as in macqc/molecule.h

BASES = ["aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z",
         "aug-cc-pCVDZ", "aug-cc-pCVTZ", "aug-cc-pCVQZ", "aug-cc-pCV5Z"]

# Default geometries of ΨLab (experimental equilibrium structures).
GEOMETRIES = {
    "H2": dict(R=0.74144 / BOHR_IN_ANGSTROM),
    "LiH": dict(R=1.5949 / BOHR_IN_ANGSTROM),
    "H2O": dict(R=0.9572 / BOHR_IN_ANGSTROM, angle=104.52),
}


def atoms(name, R, angle=None):
    if name == "H2":
        return [("H", (0.0, 0.0, -0.5 * R)), ("H", (0.0, 0.0, 0.5 * R))]
    if name == "LiH":
        return [("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, R))]
    if name == "H2O":
        half = 0.5 * angle * 3.14159265358979323846 / 180.0
        x, z = R * math.sin(half), R * math.cos(half)
        return [("O", (0.0, 0.0, 0.0)), ("H", (x, 0.0, z)), ("H", (-x, 0.0, z))]
    raise ValueError(name)


def defines(path, symbol):
    """True if the Gaussian94 file has an element line '<symbol> 0'."""
    for line in path.read_text().splitlines():
        t = line.split()
        if len(t) == 2 and t[1] == "0" and t[0].lstrip("-").lower() == symbol.lower():
            return True
    return False


def element_basis(basis, symbol):
    path = BASIS_DIR / (basis.lower() + ".gbs")
    if defines(path, symbol):
        return parse_gaussian.load(str(path), symbol), basis
    if symbol == "H" and basis.startswith("aug-cc-pCV"):
        partner = "aug-cc-pV" + basis[len("aug-cc-pCV"):]
        path = BASIS_DIR / (partner.lower() + ".gbs")
        if defines(path, "H"):
            return parse_gaussian.load(str(path), "H"), partner
    return None, None


def run(name, basis, threads=None):
    geo = GEOMETRIES[name]
    at = atoms(name, geo["R"], geo.get("angle"))
    bas = {}
    for sym, _ in at:
        b, used = element_basis(basis, sym)
        if b is None:
            return None  # basis set not defined for this element
        bas[sym] = b
    mol = gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-8
    mf.direct_scf_tol = 1e-14
    mf.max_cycle = 200
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"PySCF did not converge for {name}/{basis}")
    s_min = float(min(__import__("numpy").linalg.eigvalsh(mol.intor("int1e_ovlp"))))
    return dict(molecule=name, basis=basis, R_bohr=geo["R"], angle=geo.get("angle"),
                nbf=mol.nao, energy=float(e), nuclear_repulsion=float(mol.energy_nuc()),
                homo=float(mf.mo_energy[mol.nelectron // 2 - 1]),
                lumo=float(mf.mo_energy[mol.nelectron // 2]),
                smallest_overlap_eigenvalue=s_min)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-nbf", type=int, default=10**9)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--only", nargs="*", help="molecule/basis pairs, e.g. H2O/aug-cc-pVDZ")
    ap.add_argument("--tolerance", type=float, default=2e-9)
    args = ap.parse_args()

    stored = json.loads(OUT.read_text()) if OUT.exists() else None
    results = []
    worst = 0.0
    for name in GEOMETRIES:
        for basis in BASES:
            if args.only and f"{name}/{basis}" not in args.only:
                continue
            # rough size check before building
            r = None
            try:
                r = run_if_small(name, basis, args.max_nbf)
            except RuntimeError as exc:
                print(exc, file=sys.stderr)
                return 1
            if r is None:
                continue
            results.append(r)
            line = f"{name:4s} {basis:13s} nbf={r['nbf']:4d}  E_RHF = {r['energy']:.10f}"
            if args.check and stored:
                ref = next((c for c in stored["cases"] if c["molecule"] == name and c["basis"] == basis), None)
                if ref:
                    d = abs(ref["energy"] - r["energy"])
                    worst = max(worst, d)
                    line += f"   stored {ref['energy']:.10f}  |diff| = {d:.1e}"
            print(line, flush=True)
    if args.check:
        print(f"largest difference from the stored values: {worst:.1e} E_h")
        return 0 if worst <= args.tolerance else 1
    doc = {
        "description": "RHF total energies (E_h) from PySCF with the basis files in basis/; "
                       "generated by reference/pyscf_rhf_reference.py",
        "program": f"PySCF {pyscf.__version__}",
        "settings": "spherical harmonics, conv_tol 1e-12, conv_tol_grad 1e-8, direct_scf_tol 1e-14, "
                    f"overlap eigenvalues below {LINDEP_THRESHOLD:g} removed",
        "bohr_in_angstrom": BOHR_IN_ANGSTROM,
        "cases": results,
    }
    if stored and (args.max_nbf < 10**9 or args.only):  # keep the other cases computed earlier
        have = {(c["molecule"], c["basis"]) for c in results}
        doc["cases"] += [c for c in stored["cases"] if (c["molecule"], c["basis"]) not in have]
        order = {(m, b): i for i, (m, b) in enumerate((m, b) for m in GEOMETRIES for b in BASES)}
        doc["cases"].sort(key=lambda c: order[(c["molecule"], c["basis"])])
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


def run_if_small(name, basis, max_nbf):
    geo = GEOMETRIES[name]
    at = atoms(name, geo["R"], geo.get("angle"))
    bas = {}
    for sym, _ in at:
        b, _ = element_basis(basis, sym)
        if b is None:
            return None
        bas[sym] = b
    probe = gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, verbose=0)
    if probe.nao > max_nbf:
        return None
    return run(name, basis)


if __name__ == "__main__":
    sys.exit(main())
