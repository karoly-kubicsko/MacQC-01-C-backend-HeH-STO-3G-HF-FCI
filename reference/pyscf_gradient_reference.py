#!/usr/bin/env python3
"""Independent RHF gradients, Hessians and dipole moments for ΨLab 5.0, computed with PySCF.

Same basis files and linear-dependency threshold (1e-7) as reference/pyscf_rhf_reference.py.
The geometries are deliberately away from equilibrium and, for H2O, without symmetry, so
that every Cartesian component of the gradient is non-zero; they are given in bohr and stored
in the output file, which the C++ tests read (backend/tests/test_gradients.cpp).

PySCF: scf.RHF, its analytic gradient (mf.Gradients()) and analytic Hessian (mf.Hessian(),
coupled-perturbed HF), and mf.dip_moment() about the origin.

One case removes near-linear dependencies on purpose (threshold 1e-2 instead of 1e-7, two
combinations removed in H2O/aug-cc-pVDZ). PySCF's gradient formula does not contain the term
from the rotation of the kept orbital space with the geometry; ΨLab reports that term
separately ("lindep"), so the test compares PySCF with ΨLab's gradient minus that term and,
separately, ΨLab's full gradient with finite differences of the energy.

Usage:
  python3 reference/pyscf_gradient_reference.py              # -> backend/tests/data/gradient_reference.json
  python3 reference/pyscf_gradient_reference.py --check      # recompute and compare with the stored values
  python3 reference/pyscf_gradient_reference.py --max-nbf 120 --check
"""
import argparse
import json
import math
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402  (sets the 1e-7 linear-dependency threshold)
import numpy  # noqa: E402
import pyscf  # noqa: E402
from pyscf import gto, scf  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "gradient_reference.json"


def water_distorted():
    """H2O near 1.85 bohr / 101 degrees, displaced so that no symmetry is left."""
    R, half = 1.85, 0.5 * 101.0 * math.pi / 180.0
    x, z = R * math.sin(half), R * math.cos(half)
    return [("O", (0.0, 0.03, 0.0)), ("H", (x + 0.05, 0.0, z)), ("H", (-x, 0.0, z - 0.07))]


GEOMETRIES = {
    "H2O-distorted": water_distorted(),
    "LiH-3.10": [("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 3.10))],
    "H2-1.30": [("H", (0.0, 0.0, -0.65)), ("H", (0.0, 0.0, 0.65))],
}

# (geometry, basis, hessian?, linear-dependency threshold)
CASES = [
    ("H2O-distorted", "aug-cc-pVDZ", True, 1e-7),
    ("H2O-distorted", "aug-cc-pCVDZ", False, 1e-7),
    ("H2O-distorted", "aug-cc-pVTZ", True, 1e-7),
    ("H2O-distorted", "aug-cc-pVQZ", False, 1e-7),
    ("H2O-distorted", "aug-cc-pV5Z", False, 1e-7),
    ("H2O-distorted", "aug-cc-pVDZ", False, 1e-2),
    ("LiH-3.10", "aug-cc-pVDZ", True, 1e-7),
    ("LiH-3.10", "aug-cc-pCVTZ", False, 1e-7),
    ("H2-1.30", "aug-cc-pVTZ", True, 1e-7),
    ("H2-1.30", "aug-cc-pV5Z", False, 1e-7),
]


def build(geometry, basis):
    at = GEOMETRIES[geometry]
    bas = {s: rr.element_basis(basis, s)[0] for s, _ in at}
    return gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, verbose=0)


def run(geometry, basis, hessian, lindep):
    scf.hf.overlap_zero_eigenvalue_threshold = lindep
    try:
        mol = build(geometry, basis)
        t0 = time.time()
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.conv_tol_grad = 1e-9
        mf.direct_scf_tol = 1e-14
        mf.max_cycle = 200
        mf.kernel()
        if not mf.converged:
            raise RuntimeError(f"PySCF RHF did not converge for {geometry}/{basis}")
        g = mf.Gradients().kernel()
        out = dict(geometry=geometry, basis=basis, lindep_threshold=lindep,
                   atoms=[[s, list(xyz)] for s, xyz in GEOMETRIES[geometry]],
                   nbf=mol.nao, orbitals=int(mf.mo_coeff.shape[1]), energy=float(mf.e_tot),
                   gradient=[float(v) for v in g.ravel()],
                   dipole=[float(v) for v in mf.dip_moment(unit="AU", verbose=0)])
        if hessian:
            hs = mf.Hessian()
            hs.conv_tol = 1e-10
            h = hs.kernel()  # (natm, natm, 3, 3)
            n = mol.natm
            out["hessian"] = [[float(h[i // 3, j // 3, i % 3, j % 3]) for j in range(3 * n)] for i in range(3 * n)]
        out["seconds"] = round(time.time() - t0, 2)
        return out
    finally:
        scf.hf.overlap_zero_eigenvalue_threshold = rr.LINDEP_THRESHOLD


def nbf_of(geometry, basis):
    return build(geometry, basis).nao


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--max-nbf", type=int, default=10**9)
    ap.add_argument("--tolerance", type=float, default=1e-8)
    args = ap.parse_args()
    stored = json.loads(OUT.read_text()) if OUT.exists() else None
    results, worst = [], 0.0
    for geometry, basis, hessian, lindep in CASES:
        if nbf_of(geometry, basis) > args.max_nbf:
            continue
        r = run(geometry, basis, hessian, lindep)
        results.append(r)
        line = f"{geometry:14s} {basis:13s} lindep {lindep:g}  nbf {r['nbf']:4d}  E = {r['energy']:.10f}  ({r['seconds']} s)"
        if args.check and stored:
            ref = next((c for c in stored["cases"] if c["geometry"] == geometry and c["basis"] == basis
                        and c["lindep_threshold"] == lindep), None)
            if ref:
                d = max(abs(a - b) for a, b in zip(ref["gradient"], r["gradient"]))
                d = max(d, abs(ref["energy"] - r["energy"]))
                if hessian:
                    d = max(d, float(numpy.max(numpy.abs(numpy.array(ref["hessian"]) - numpy.array(r["hessian"])))))
                worst = max(worst, d)
                line += f"  |diff| {d:.1e}"
        print(line, flush=True)
    if args.check:
        print(f"largest difference from the stored values: {worst:.1e}")
        return 0 if worst <= args.tolerance else 1
    doc = {
        "description": "RHF energies (E_h), gradients (E_h/bohr), Hessians (E_h/bohr^2) and dipole moments "
                       "(e bohr, origin 0) from PySCF; generated by reference/pyscf_gradient_reference.py",
        "program": f"PySCF {pyscf.__version__}",
        "settings": "spherical harmonics, conv_tol 1e-12, conv_tol_grad 1e-9, direct_scf_tol 1e-14, CPHF conv_tol "
                    "1e-10; overlap eigenvalues below lindep_threshold removed (canonical orthogonalisation)",
        "cases": results,
    }
    if stored and args.max_nbf < 10**9:
        have = {(c["geometry"], c["basis"], c["lindep_threshold"]) for c in results}
        doc["cases"] += [c for c in stored["cases"]
                         if (c["geometry"], c["basis"], c["lindep_threshold"]) not in have]
        order = {(g, b, t): i for i, (g, b, _, t) in enumerate(CASES)}
        doc["cases"].sort(key=lambda c: order[(c["geometry"], c["basis"], c["lindep_threshold"])])
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
