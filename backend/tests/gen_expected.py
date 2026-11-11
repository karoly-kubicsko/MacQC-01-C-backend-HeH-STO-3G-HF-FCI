#!/usr/bin/env python3
"""Regenerate the expected HF reports in backend/tests/data/ from the reference script.

Usage: python3 backend/tests/gen_expected.py
"""
import io
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "reference"))
from hartree_fock_sto3g import hfcalc  # noqa: E402

DATA = pathlib.Path(__file__).resolve().parent / "data"

# name -> hfcalc keyword arguments
CASES = {
    "hf_sto3g_level0": dict(iop=0, n=3),
    "hf_sto3g_level1": dict(iop=1, n=3),
    "hf_sto3g_level2": dict(iop=2, n=3),
    "hf_sto2g_level2": dict(iop=2, n=2),
    "hf_sto1g_level2": dict(iop=2, n=1),
    "hf_sto3g_R2.5_level2": dict(iop=2, n=3, r=2.5),
}

DATA.mkdir(parents=True, exist_ok=True)
for name, kwargs in CASES.items():
    buf = io.StringIO()
    hfcalc(out=buf, **kwargs)
    (DATA / f"{name}.txt").write_text(buf.getvalue(), newline="\n")
    print(f"wrote {name}.txt ({len(buf.getvalue())} bytes)")
