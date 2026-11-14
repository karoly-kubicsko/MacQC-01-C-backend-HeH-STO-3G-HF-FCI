#!/usr/bin/env python3
"""Optimised geometries for ΨLab 5.0 with PySCF analytic gradients and geomeTRIC.

PySCF computes the energies and analytic gradients (RHF, CCSD, CCSD(T) through
pyscf.grad.ccsd_t, Full CI as CASCI with every orbital active), geomeTRIC (an independent
optimiser in translation-rotation-internal coordinates) minimises with tight thresholds
(gradient max 1.5e-6, RMS 1e-6 E_h/bohr; displacement 1.5e-5 / 1e-5 bohr; energy 1e-10 E_h).
Same basis files, all electrons correlated, starting from ΨLab's default (experimental)
geometries. The results (bond lengths in angstrom, H-O-H angle in degrees, energies) are
stored for backend/tests/test_optimizer.cpp.

Usage (needs: pip install pyscf geometric):
  python3 reference/pyscf_optimization_reference.py           # -> backend/tests/data/optimization_reference.json
  python3 reference/pyscf_optimization_reference.py --check   # recompute and compare
"""
import argparse
import json
import math
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402
import numpy  # noqa: E402
import pyscf  # noqa: E402
from pyscf import cc, gto, lib, mcscf, scf  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402
from pyscf.grad import ccsd_t as ccsd_t_grad  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "optimization_reference.json"

CASES = [
    ("H2O", "aug-cc-pVDZ", "RHF"),
    ("H2O", "aug-cc-pVTZ", "RHF"),
    ("H2O", "aug-cc-pVDZ", "CCSD"),
    ("H2O", "aug-cc-pVDZ", "CCSD(T)"),
    ("LiH", "aug-cc-pVDZ", "RHF"),
    ("LiH", "aug-cc-pVDZ", "CCSD(T)"),
    ("LiH", "aug-cc-pVDZ", "FCI"),
    ("H2", "aug-cc-pVTZ", "FCI"),
]

CONV = dict(convergence_energy=1e-10, convergence_grms=1e-6, convergence_gmax=1.5e-6,
            convergence_drms=1e-5, convergence_dmax=1.5e-5)


def build(name, basis):
    geo = rr.GEOMETRIES[name]
    at = rr.atoms(name, geo["R"], geo.get("angle"))
    bas = {s: rr.element_basis(basis, s)[0] for s, _ in at}
    return gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, verbose=0)


def parameters(mol):
    x = mol.atom_coords(unit="Bohr")
    bohr = rr.BOHR_IN_ANGSTROM
    if mol.natm == 2:
        return dict(R=float(numpy.linalg.norm(x[1] - x[0]) * bohr))
    u, v = x[1] - x[0], x[2] - x[0]
    ang = math.degrees(math.acos(numpy.dot(u, v) / (numpy.linalg.norm(u) * numpy.linalg.norm(v))))
    return dict(R1=float(numpy.linalg.norm(u) * bohr), R2=float(numpy.linalg.norm(v) * bohr), angle=float(ang))


class CCSDTScanner(lib.GradScanner):
    """CCSD(T) energy and gradient for geomeTRIC.

    PySCF's own scanner of pyscf.grad.ccsd_t.Gradients (inherited from the CCSD gradient
    scanner) solves the CCSD Lambda equations, passes those multipliers to the CCSD(T)
    gradient and returns the CCSD energy, so geomeTRIC would minimise the CCSD energy with an
    inconsistent gradient. This scanner returns the CCSD(T) energy and the gradient with the
    CCSD(T) Lambda equations (ccsd_t.Gradients(...).kernel() without multipliers)."""

    def __init__(self, g):
        lib.GradScanner.__init__(self, g)
        self._ok = False

    def __call__(self, mol_or_geom, **kwargs):
        mol = mol_or_geom if isinstance(mol_or_geom, gto.MoleBase) else self.mol.set_geom_(mol_or_geom, inplace=False)
        self.mol = mol
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.conv_tol_grad = 1e-9
        mf.direct_scf_tol = 1e-14
        mf.kernel()
        mc = cc.CCSD(mf)
        mc.conv_tol = 1e-12
        mc.conv_tol_normt = 1e-9
        mc.max_cycle = 300
        mc.kernel()
        et = mc.ccsd_t()
        g = ccsd_t_grad.Gradients(mc).kernel()
        self._ok = bool(mf.converged and mc.converged)
        return mc.e_tot + et, g

    @property
    def converged(self):
        return self._ok


def solver(mol, method):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.direct_scf_tol = 1e-14
    if method == "RHF":
        return mf
    mf.kernel()
    if method in ("CCSD", "CCSD(T)"):
        mc = cc.CCSD(mf)
        mc.conv_tol = 1e-12
        mc.conv_tol_normt = 1e-9
        mc.max_cycle = 300
        if method == "CCSD":
            return mc
        mc.kernel()
        return CCSDTScanner(ccsd_t_grad.Gradients(mc))
    mcas = mcscf.CASCI(mf, mol.nao, mol.nelectron)
    mcas.fcisolver.conv_tol = 1e-12
    return mcas


def energy(mol, method):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    if method == "RHF":
        return mf.e_tot
    if method in ("CCSD", "CCSD(T)"):
        mc = cc.CCSD(mf)
        mc.conv_tol = 1e-12
        mc.conv_tol_normt = 1e-9
        mc.kernel()
        return mc.e_tot + (mc.ccsd_t() if method == "CCSD(T)" else 0.0)
    mcas = mcscf.CASCI(mf, mol.nao, mol.nelectron)
    mcas.fcisolver.conv_tol = 1e-12
    return mcas.kernel()[0]


def run(name, basis, method):
    t0 = time.time()
    mol = build(name, basis)
    mol_eq = optimize(solver(mol, method), maxsteps=100, **CONV)
    out = dict(molecule=name, basis=basis, method=method, energy=float(energy(mol_eq, method)))
    out.update(parameters(mol_eq))
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--tolerance", type=float, default=2e-6)
    ap.add_argument("--only", nargs="*", help="cases as molecule/basis/method, e.g. H2O/aug-cc-pVDZ/CCSD(T)")
    args = ap.parse_args()
    stored = json.loads(OUT.read_text()) if OUT.exists() else None
    results, worst = [], 0.0
    for name, basis, method in CASES:
        if args.only and f"{name}/{basis}/{method}" not in args.only:
            if stored:
                old = next((c for c in stored["cases"]
                            if (c["molecule"], c["basis"], c["method"]) == (name, basis, method)), None)
                if old and not args.check:
                    results.append(old)
            continue
        r = run(name, basis, method)
        results.append(r)
        geo = ", ".join(f"{k} = {r[k]:.6f}" for k in ("R", "R1", "R2", "angle") if k in r)
        line = f"{name:4s} {basis:12s} {method:8s} E = {r['energy']:.10f}  {geo}  ({r['seconds']} s)"
        if args.check and stored:
            ref = next((c for c in stored["cases"] if (c["molecule"], c["basis"], c["method"]) == (name, basis, method)), None)
            if ref:
                d = max(abs(ref[k] - r[k]) for k in ("R", "R1", "R2", "angle") if k in r)
                worst = max(worst, d)
                line += f"  |diff| {d:.1e}"
        print(line, flush=True)
    if args.check:
        print(f"largest difference from the stored values: {worst:.1e}")
        return 0 if worst <= args.tolerance else 1
    doc = {
        "description": "Optimised geometries (angstrom, degrees) and energies (E_h) from PySCF analytic gradients and "
                       "the geomeTRIC optimiser; generated by reference/pyscf_optimization_reference.py",
        "program": f"PySCF {pyscf.__version__}, geomeTRIC",
        "settings": "all electrons correlated; geomeTRIC thresholds " + json.dumps(CONV),
        "bohr_in_angstrom": rr.BOHR_IN_ANGSTROM,
        "cases": results,
    }
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
