#!/usr/bin/env python3
"""Independent Full CI reference energies for ΨLab 3.0, computed with PySCF.

Same molecules, geometries, basis files and linear-dependency threshold as
reference/pyscf_rhf_reference.py; all electrons are correlated (no frozen core), and the
CI uses the canonical RHF orbitals (the FCI energy does not depend on that choice).
LiH: PySCF's determinant Full CI solver fci.direct_spin0 (singlet ground state).
H2: coupled cluster with singles and doubles (CCSD), which is exact - identical to Full CI -
for two electrons; PySCF's determinant solver is impractically slow with 160 orbitals.

Usage:
  python3 reference/pyscf_fci_reference.py                   # -> backend/tests/data/fci_reference.json
  python3 reference/pyscf_fci_reference.py --check --max-det 300000 --max-orbitals 100
"""
import argparse
import json
import time
import pathlib
import sys
from math import comb

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402  (sets the 1e-7 linear-dependency threshold)
from pyscf import cc, fci, gto, scf  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "fci_reference.json"

# Cases ΨLab allows (H2O and the largest LiH bases are refused).
CASES = [("H2", b) for b in ("aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z")] + \
        [("LiH", b) for b in ("aug-cc-pVDZ", "aug-cc-pCVDZ", "aug-cc-pVTZ", "aug-cc-pCVTZ")]


def build(name, basis, R=None):
    geo = dict(rr.GEOMETRIES[name])
    if R is not None:
        geo["R"] = R
    at = rr.atoms(name, geo["R"], geo.get("angle"))
    bas = {s: rr.element_basis(basis, s)[0] for s, _ in at}
    return geo, gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, verbose=0)


def run(name, basis, R=None):
    geo, mol = build(name, basis, R)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-8
    mf.direct_scf_tol = 1e-14
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"RHF did not converge for {name}/{basis}")
    norb = mf.mo_coeff.shape[1]
    na = mol.nelectron // 2
    if mol.nelectron == 2:
        ccsd = cc.RCCSD(mf)
        ccsd.conv_tol = 1e-13
        ccsd.conv_tol_normt = 1e-10
        ccsd.max_cycle = 200
        ccsd.kernel()
        if not ccsd.converged:
            raise RuntimeError(f"CCSD did not converge for {name}/{basis}")
        e = mf.e_tot + ccsd.e_corr
        t1, t2 = ccsd.t1, ccsd.t2
        # |C_0|^2 of the normalised wave function (1 + T1 + T2 + T1^2/2)|0> for two electrons:
        # C_ia^ab-type amplitudes: c_ij^ab = t2 + t1 t1, singles c_i^a = t1 (each spin).
        import numpy as np
        c2 = t2[0, 0] + np.einsum("a,b->ab", t1[0], t1[0])
        norm2 = 1.0 + 2.0 * np.dot(t1[0], t1[0]) + np.sum(c2 * c2)
        return dict(molecule=name, basis=basis, R_bohr=geo["R"], orbitals=norb, determinants=comb(norb, na) ** 2,
                    rhf_energy=float(mf.e_tot), fci_energy=float(e), correlation_energy=float(e - mf.e_tot),
                    s2=0.0, hf_weight=float(1.0 / norm2), method="CCSD (exact for two electrons)")
    solver = fci.direct_spin0.FCI(mol)
    solver.conv_tol = 1e-12
    solver.max_cycle = 200
    solver.max_space = 12
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = pyscf.ao2mo.kernel(mol, mf.mo_coeff)
    e, civec = solver.kernel(h1, eri, norb, (na, na), ecore=mol.energy_nuc())
    # PySCF's spin_square supports fewer than 64 orbitals only.
    ss = solver.spin_square(civec, norb, (na, na))[0] if norb < 64 else None
    return dict(molecule=name, basis=basis, R_bohr=geo["R"], orbitals=norb,
                determinants=comb(norb, na) ** 2, rhf_energy=float(mf.e_tot),
                fci_energy=float(e), correlation_energy=float(e - mf.e_tot), s2=None if ss is None else float(ss),
                hf_weight=float(civec[0, 0] ** 2), method="FCI (fci.direct_spin0)")


import pyscf  # noqa: E402
import pyscf.ao2mo  # noqa: E402,F401


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-det", type=float, default=1e12)
    ap.add_argument("--max-orbitals", type=int, default=10**6)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--tolerance", type=float, default=2e-8)
    args = ap.parse_args()
    stored = json.loads(OUT.read_text()) if OUT.exists() else None
    results, worst = [], 0.0
    for name, basis in CASES:
        _, mol = build(name, basis)
        if comb(mol.nao, mol.nelectron // 2) ** 2 > args.max_det or mol.nao > args.max_orbitals:
            continue
        t0 = time.time()
        r = run(name, basis)
        r["pyscf_seconds"] = round(time.time() - t0, 1)
        results.append(r)
        line = f"{name:4s} {basis:13s} n={r['orbitals']:4d} det={r['determinants']:>10d}  E_FCI = {r['fci_energy']:.10f}"
        if args.check and stored:
            ref = next((c for c in stored["cases"] if c["molecule"] == name and c["basis"] == basis), None)
            if ref:
                d = abs(ref["fci_energy"] - r["fci_energy"])
                worst = max(worst, d)
                line += f"  stored {ref['fci_energy']:.10f}  |diff| = {d:.1e}"
        print(line, flush=True)
        if not args.check:
            save(results, stored)
    if args.check:
        print(f"largest difference from the stored values: {worst:.1e} E_h")
        return 0 if worst <= args.tolerance else 1
    print(f"wrote {OUT}")
    return 0


def save(results, stored):
    """Writes the computed cases, keeping stored ones that were not recomputed."""
    cases = list(results)
    if stored:
        have = {(c["molecule"], c["basis"]) for c in cases}
        cases += [c for c in stored["cases"] if (c["molecule"], c["basis"]) not in have]
    order = {k: i for i, k in enumerate(CASES)}
    cases.sort(key=lambda c: order[(c["molecule"], c["basis"])])
    OUT.write_text(json.dumps({
        "description": "Full CI total energies (E_h), all electrons, from PySCF with the basis files in basis/; "
                       "generated by reference/pyscf_fci_reference.py",
        "program": f"PySCF {pyscf.__version__}",
        "settings": "RHF orbitals (conv_tol 1e-12); LiH: fci.direct_spin0, conv_tol 1e-12; H2: RCCSD "
                    "(exact for two electrons), conv_tol 1e-13, conv_tol_normt 1e-10; overlap eigenvalues "
                    "below 1e-7 removed",
        "cases": cases}, indent=1) + "\n")


if __name__ == "__main__":
    sys.exit(main())
