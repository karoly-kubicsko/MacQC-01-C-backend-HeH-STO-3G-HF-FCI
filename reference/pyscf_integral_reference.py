#!/usr/bin/env python3
"""Orientation-independent integral checks for ΨLab's Gaussian integral engine, from PySCF.

For a few molecule/basis combinations, and for a small test basis that contains every
angular momentum up to h (l = 5) with contracted and uncontracted shells, this script
writes quantities that do not depend on the order or sign convention of the spherical
functions within a shell:
  trace and Frobenius norm of S, T, V;  Frobenius norm of the ERI tensor (pq|rs);
  sum_pq (pp|qq)  and  sum_pq (pq|pq).
Output: backend/tests/data/integral_reference.json      Needs: pip install pyscf
With --check the values are recomputed and compared with the stored file instead.
"""
import argparse
import json
import math
import pathlib

import numpy as np
import pyscf
from pyscf import gto
from pyscf.gto.basis import parse_gaussian

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "backend" / "tests" / "data" / "integral_reference.json"
BOHR_IN_ANGSTROM = 0.529177210544

# Test basis with l = 0..5 (exponent, coefficient pairs); atoms in bohr.
TEST_ATOMS = [("O", (0.0, 0.0, 0.0)), ("H", (1.43, 0.0, 1.11)), ("Li", (-1.43, 0.3, 1.11))]
TEST_SHELLS = [
    (0, 0, [(130.7, 0.154), (23.8, 0.535), (6.44, 0.444)]),
    (0, 1, [(5.03, 0.156), (1.17, 0.607)]),
    (0, 2, [(1.2, 1.0)]),
    (0, 3, [(0.9, 1.0)]),
    (0, 4, [(0.8, 1.0)]),
    (0, 5, [(1.5, 0.6), (0.5, 0.5)]),
    (1, 0, [(3.4, 0.15), (0.62, 0.9)]),
    (1, 1, [(0.7, 1.0)]),
    (1, 2, [(0.9, 1.0)]),
    (2, 0, [(0.3, 1.0)]),
    (2, 3, [(1.1, 1.0)]),
    (2, 4, [(0.6, 1.0)]),
]


def invariants(mol):
    S, T, V = mol.intor("int1e_ovlp"), mol.intor("int1e_kin"), mol.intor("int1e_nuc")
    eri = mol.intor("int2e", aosym="s1")
    n = mol.nao
    return dict(
        nbf=n,
        trace_S=float(np.trace(S)), norm_S=float(np.linalg.norm(S)),
        trace_T=float(np.trace(T)), norm_T=float(np.linalg.norm(T)),
        trace_V=float(np.trace(V)), norm_V=float(np.linalg.norm(V)),
        norm_ERI=float(np.linalg.norm(eri)),
        coulomb_sum=float(np.einsum("ppqq->", eri)),
        exchange_sum=float(np.einsum("pqpq->", eri)),
    )


def water(basis):
    R = 0.9572 / BOHR_IN_ANGSTROM
    half = 0.5 * 104.52 * 3.14159265358979323846 / 180.0
    x, z = R * math.sin(half), R * math.cos(half)
    path = ROOT / "basis" / (basis.lower() + ".gbs")
    bas = {s: parse_gaussian.load(str(path), s) for s in ("O", "H")}
    return gto.M(atom=[["O", (0, 0, 0)], ["H", (x, 0, z)], ["H", (-x, 0, z)]], basis=bas,
                 unit="Bohr", cart=False, verbose=0)


def test_molecule():
    atom = [[f"{s}{i}", xyz] for i, (s, xyz) in enumerate(TEST_ATOMS)]
    basis = {f"{s}{i}": [[l] + [[e, c] for e, c in prims] for a, l, prims in TEST_SHELLS if a == i]
             for i, (s, _) in enumerate(TEST_ATOMS)}
    return gto.M(atom=atom, basis=basis, unit="Bohr", cart=False, verbose=0, spin=None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    cases = []
    r = invariants(test_molecule())
    r.update(name="test-basis-l5", atoms=[[s, list(x)] for s, x in TEST_ATOMS],
             shells=[[a, l, [list(p) for p in prims]] for a, l, prims in TEST_SHELLS])
    cases.append(r)
    for basis in ("aug-cc-pVDZ", "aug-cc-pVTZ"):
        r = invariants(water(basis))
        r.update(name=f"H2O/{basis}")
        cases.append(r)
    if args.check:
        stored = {c["name"]: c for c in json.loads(OUT.read_text())["cases"]}
        worst = 0.0
        for c in cases:
            for k, v in c.items():
                if isinstance(v, float):
                    worst = max(worst, abs(v - stored[c["name"]][k]) / max(1.0, abs(v)))
        print(f"largest relative difference from the stored values: {worst:.1e}")
        return 0 if worst < 1e-11 else 1
    OUT.write_text(json.dumps({"program": f"PySCF {pyscf.__version__}",
                               "generated_by": "reference/pyscf_integral_reference.py",
                               "cases": cases}, indent=1) + "\n")
    print(f"wrote {OUT}")
    for c in cases:
        print(c["name"], c["nbf"], c["norm_ERI"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
