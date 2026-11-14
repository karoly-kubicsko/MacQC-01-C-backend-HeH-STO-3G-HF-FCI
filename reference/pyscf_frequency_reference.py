#!/usr/bin/env python3
"""Harmonic frequencies and thermochemistry references for ΨLab 5.0.

At the RHF geometries optimised with PySCF + geomeTRIC (backend/tests/data/optimization_reference.json;
H2 at the Full CI/aug-cc-pVTZ geometry), PySCF computes the analytic RHF Hessian, its harmonic
analysis (pyscf.hessian.thermo.harmonic_analysis: mass-weighted Hessian, translations and
rotations projected out) and its thermochemistry (pyscf.hessian.thermo.thermo) at 298.15 K and
101325 Pa. PySCF is given the same AME 2020 isotope masses through mol.nucprop (by default its harmonic
analysis and thermochemistry use standard atomic weights, e.g. 1.008 for H); its physical
constants are older CODATA values, which differ from CODATA 2018 in the 7th-9th digit.

Independently, the thermochemistry is evaluated here "by hand" from the textbook formulas with
CODATA 2018 constants and AME 2020 masses (ideal-gas translation, rigid rotor, harmonic
oscillators, singlet ground state), from PySCF's frequencies and the geometry.

Usage:
  python3 reference/pyscf_frequency_reference.py           # -> backend/tests/data/frequency_reference.json
  python3 reference/pyscf_frequency_reference.py --check   # recompute and compare
"""
import argparse
import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402
import numpy  # noqa: E402
import pyscf  # noqa: E402
from pyscf import gto, scf  # noqa: E402
from pyscf.hessian import thermo  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "frequency_reference.json"
OPT = rr.ROOT / "backend" / "tests" / "data" / "optimization_reference.json"

# CODATA 2018
C = 299792458.0
H = 6.62607015e-34
KB = 1.380649e-23
NA = 6.02214076e23
EH = 4.3597447222071e-18
A0 = 5.29177210903e-11
U = 1.66053906660e-27
R = NA * KB
# AME 2020
MASS = {"H": 1.007825031898, "Li": 7.016003434, "O": 15.994914619257}

CASES = [("H2O", "aug-cc-pVDZ", "RHF"), ("H2O", "aug-cc-pVTZ", "RHF"), ("LiH", "aug-cc-pVDZ", "RHF"),
         ("H2", "aug-cc-pVTZ", "FCI")]


def geometry(case):
    """Atoms (bohr) from the optimised parameters; H2O in the xz plane as in ΨLab."""
    b = rr.BOHR_IN_ANGSTROM
    if case["molecule"] == "H2O":
        r, ang = case["R1"] / b, math.radians(case["angle"])
        x, z = r * math.sin(ang / 2), r * math.cos(ang / 2)
        return [("O", (0.0, 0.0, 0.0)), ("H", (x, 0.0, z)), ("H", (-x, 0.0, z))]
    r = case["R"] / b
    if case["molecule"] == "LiH":
        return [("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, r))]
    return [("H", (0.0, 0.0, -0.5 * r)), ("H", (0.0, 0.0, 0.5 * r))]


def hand_thermo(atoms, freqs_cm, sigma, T=298.15, P=101325.0, e_elec=0.0):
    m = numpy.array([MASS[s] for s, _ in atoms])
    x = numpy.array([p for _, p in atoms])
    com = (m[:, None] * x).sum(0) / m.sum()
    y = x - com
    I = sum(mi * (numpy.dot(r, r) * numpy.eye(3) - numpy.outer(r, r)) for mi, r in zip(m, y))
    mom = numpy.linalg.eigvalsh(I)  # u bohr^2
    linear = mom[0] < 1e-6 * mom[2]
    M = m.sum() * U
    kT = KB * T
    qt = (2 * math.pi * M * kT / H**2) ** 1.5 * kT / P
    St, Et, Cvt = R * (math.log(qt) + 2.5), 1.5 * R * T, 1.5 * R
    theta = [H * H / (8 * math.pi**2 * mi * U * A0 * A0 * KB) for mi in mom if mi > 1e-6 * mom[2]]
    if linear:
        qr = T / (sigma * theta[-1])
        Sr, Er, Cvr = R * (math.log(qr) + 1), R * T, R
    else:
        qr = math.sqrt(math.pi) / sigma * T**1.5 / math.sqrt(theta[0] * theta[1] * theta[2])
        Sr, Er, Cvr = R * (math.log(qr) + 1.5), 1.5 * R * T, 1.5 * R
    zpe = Ev = Sv = Cvv = 0.0
    for f in freqs_cm:
        th = H * C * f * 100 / KB
        u = th / T
        zpe += 0.5 * R * th
        Ev += R * th * (0.5 + 1 / math.expm1(u))
        Sv += R * (u / math.expm1(u) - math.log(-math.expm1(-u)))
        Cvv += R * u * u * math.exp(u) / math.expm1(u) ** 2
    to_h = 1.0 / (EH * NA)
    U_ = Et + Er + Ev
    S = St + Sr + Sv
    return dict(ZPE=zpe * to_h, U=U_ * to_h, H=(U_ + R * T) * to_h, G=(U_ + R * T - T * S) * to_h,
                A=(U_ - T * S) * to_h, S=S, S_trans=St, S_rot=Sr, S_vib=Sv, Cv=Cvt + Cvr + Cvv,
                moments=[float(v) for v in mom])


def run(case):
    atoms = geometry(case)
    bas = {s: rr.element_basis(case["basis"], s)[0] for s, _ in atoms}
    nucprop = {s: {"mass": MASS[s]} for s, _ in atoms}
    mol = gto.M(atom=[[s, p] for s, p in atoms], basis=bas, unit="Bohr", verbose=0)
    mol.nucprop = nucprop
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.kernel()
    hs = mf.Hessian()
    hs.conv_tol = 1e-10
    h = hs.kernel()
    ha = thermo.harmonic_analysis(mol, h)
    freqs = [float(f) for f in numpy.real(ha["freq_wavenumber"])]
    sigma = {"H2O": 2, "LiH": 1, "H2": 2}[case["molecule"]]
    molsym = gto.M(atom=[[s, p] for s, p in atoms], basis=bas, unit="Bohr", verbose=0, symmetry=True)
    molsym.nucprop = nucprop
    mf.mol = molsym
    t = thermo.thermo(mf, ha["freq_au"], 298.15, 101325)
    pt = {k: float(numpy.real(v[0])) for k, v in t.items()
          if isinstance(v, tuple) and k in ("ZPE", "E_tot", "H_tot", "G_tot", "S_tot", "Cv_tot",
                                             "S_trans", "S_rot", "S_vib", "E0")}
    return dict(molecule=case["molecule"], basis=case["basis"], method="RHF",
                atoms=[[s, list(p)] for s, p in atoms], energy=float(mf.e_tot),
                frequencies=freqs, reduced_masses=[float(v) for v in ha["reduced_mass"]],
                force_constants=[float(v) for v in ha["force_const_dyne"]],
                pyscf_thermo=pt, symmetry_number=sigma,
                hand_thermo=hand_thermo(atoms, freqs, sigma, e_elec=float(mf.e_tot)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--tolerance", type=float, default=1e-6)
    args = ap.parse_args()
    opt = json.loads(OPT.read_text())
    stored = json.loads(OUT.read_text()) if OUT.exists() else None
    results, worst = [], 0.0
    for name, basis, method in CASES:
        case = next(c for c in opt["cases"] if (c["molecule"], c["basis"], c["method"]) == (name, basis, method))
        case = dict(case, basis=basis if method == "RHF" else "aug-cc-pVTZ")
        r = run(case)
        results.append(r)
        line = f"{name:4s} {r['basis']:12s} RHF  frequencies " + " ".join(f"{f:.3f}" for f in r["frequencies"])
        if args.check and stored:
            ref = next((c for c in stored["cases"] if (c["molecule"], c["basis"]) == (name, r["basis"])), None)
            if ref:
                d = max(abs(a - b) for a, b in zip(ref["frequencies"], r["frequencies"]))
                worst = max(worst, d)
                line += f"  |diff| {d:.1e}"
        print(line, flush=True)
    if args.check:
        print(f"largest difference from the stored values: {worst:.1e} cm-1")
        return 0 if worst <= 1e-3 else 1
    doc = {
        "description": "RHF harmonic frequencies (cm-1), reduced masses (u), force constants (mdyn/A) and "
                       "thermochemistry at 298.15 K, 101325 Pa from PySCF (analytic Hessian, "
                       "pyscf.hessian.thermo), and the thermochemistry evaluated by hand with CODATA 2018 and "
                       "AME 2020; generated by reference/pyscf_frequency_reference.py",
        "program": f"PySCF {pyscf.__version__}",
        "units": "pyscf_thermo and hand_thermo: energies in E_h (E_tot, H_tot, G_tot include the electronic "
                 "energy in pyscf_thermo; U, H, G, A, ZPE in hand_thermo are corrections), entropies and Cv "
                 "in E_h/K (PySCF) or J/(mol K) (hand)",
        "cases": results,
    }
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
