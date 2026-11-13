#!/usr/bin/env python3
"""Reference energies for the coupled-cluster tests (backend/tests/data/cc_reference.json).

Usage: python3 backend/tests/gen_cc_reference.py

Each case is a set of spatial MO integrals (h, (pq|rs)) with a closed-shell reference.
CCSD, CCSDT and CCSDTQ are solved by brute force: the cluster operator is a matrix in
the M_S = 0 determinant space and the residuals <mu| exp(-T) H exp(T) |0> are formed
by matrix-vector products, so no coupled-cluster equations are coded here. Full CI is
the lowest eigenvalue of H in the same space. Cases:

  random-3o4v  random integrals with the permutational symmetry of real orbitals
               (3 doubly occupied, 4 virtual orbitals; non-diagonal Fock matrix);
  h2o-sto3g    H2O, STO-3G (PySCF integrals), 5 occupied, 2 virtual orbitals;
  lih-631g     LiH, 6-31G (PySCF integrals), 2 occupied, 9 virtual orbitals
               (CCSDTQ = Full CI for four electrons).

The PySCF CCSD and CCSD(T) energies of the molecular cases are stored as a check of the
brute-force solver and of the closed-shell implementation.
"""
import itertools
import json
import pathlib

import numpy as np

DATA = pathlib.Path(__file__).resolve().parent / "data"


def spin_orbital_hamiltonian(h, eri, nocc):
    """Spin-orbital h and antisymmetrized <pq||rs>; order: alpha orbitals then beta."""
    n = h.shape[0]
    N = 2 * n
    spatial = np.array(list(range(n)) * 2)
    spin = np.array([0] * n + [1] * n)
    hs = np.zeros((N, N))
    for p in range(N):
        for q in range(N):
            if spin[p] == spin[q]:
                hs[p, q] = h[spatial[p], spatial[q]]
    # <pq|rs> = (pr|qs) delta(sp,sr) delta(sq,ss)
    g = eri[np.ix_(spatial, spatial, spatial, spatial)]  # (pr|qs) indexed [p,r,q,s]
    same = (spin[:, None] == spin[None, :]).astype(float)
    phys = g.transpose(0, 2, 1, 3) * same[:, None, :, None] * same[None, :, None, :]
    return hs, phys - phys.transpose(0, 1, 3, 2), spin


class Space:
    def __init__(self, n, nocc):
        self.n, self.nocc = n, nocc
        alpha = [sum(1 << i for i in c) for c in itertools.combinations(range(n), nocc)]
        self.dets = [(a, b) for a in alpha for b in alpha]
        self.index = {d: k for k, d in enumerate(self.dets)}
        ref = (1 << nocc) - 1
        self.ref = self.index[(ref, ref)]

    def to_bits(self, det):
        a, b = det
        return a | (b << self.n)

    def from_bits(self, bits):
        mask = (1 << self.n) - 1
        return (bits & mask, bits >> self.n)


def apply_ops(ops, bits):
    sign = 1
    for orb, creation in reversed(ops):
        bit = 1 << orb
        if creation == bool(bits & bit):
            return None
        if bin(bits & (bit - 1)).count("1") % 2:
            sign = -sign
        bits ^= bit
    return bits, sign


def build_hamiltonian(space, hs, V):
    N = 2 * space.n
    dim = len(space.dets)
    H = np.zeros((dim, dim))
    for col, det in enumerate(space.dets):
        bits = space.to_bits(det)
        occ = [i for i in range(N) if bits >> i & 1]
        for q in occ:
            for p in range(N):
                if hs[p, q] != 0.0:
                    r = apply_ops([(p, True), (q, False)], bits)
                    if r:
                        H[space.index[space.from_bits(r[0])], col] += r[1] * hs[p, q]
        for r_, s in itertools.combinations(occ, 2):
            for p, q in itertools.combinations(range(N), 2):
                if V[p, q, r_, s] != 0.0:
                    res = apply_ops([(p, True), (q, True), (s, False), (r_, False)], bits)
                    if res:
                        H[space.index[space.from_bits(res[0])], col] += res[1] * V[p, q, r_, s]
    return H


def excitations(n, nocc, max_rank):
    """M_S-conserving excitations: list of (A, I) spin-orbital tuples (sorted)."""
    occ_a = list(range(nocc))
    vir_a = list(range(nocc, n))
    occ_b = [n + i for i in occ_a]
    vir_b = [n + a for a in vir_a]
    out = []
    for rank in range(1, max_rank + 1):
        for ka in range(rank + 1):
            kb = rank - ka
            for Ia in itertools.combinations(occ_a, ka):
                for Ib in itertools.combinations(occ_b, kb):
                    for Aa in itertools.combinations(vir_a, ka):
                        for Ab in itertools.combinations(vir_b, kb):
                            out.append((Aa + Ab, Ia + Ib))
    return out


def solve_cc(space, H, hs, V, max_rank, tol=1e-11, max_iter=200):
    n, nocc = space.n, space.nocc
    dim = len(space.dets)
    exc = excitations(n, nocc, max_rank)
    ref_bits = space.to_bits(space.dets[space.ref])
    # Operator structure and target determinants
    struct = []
    target = []
    for A, I in exc:
        ops = [(a, True) for a in A] + [(i, False) for i in reversed(I)]
        rows, cols, signs = [], [], []
        for col, det in enumerate(space.dets):
            r = apply_ops(ops, space.to_bits(det))
            if r:
                rows.append(space.index[space.from_bits(r[0])])
                cols.append(col)
                signs.append(r[1])
        struct.append((np.array(rows), np.array(cols), np.array(signs, dtype=float)))
        bits, sign = apply_ops(ops, ref_bits)
        target.append((space.index[space.from_bits(bits)], sign))
    # Fock diagonal for the denominators
    occ_so = [i for i in range(2 * n) if ref_bits >> i & 1]
    f = hs + np.einsum("piqi->pq", V[:, occ_so][:, :, :, occ_so])
    eps = np.diag(f)
    denom = np.array([eps[list(A)].sum() - eps[list(I)].sum() for A, I in exc])
    e_ref = H[space.ref, space.ref]
    HN = H - e_ref * np.eye(dim)
    t = np.zeros(len(exc))
    diis_t, diis_e = [], []

    def residual(t):
        T = np.zeros((dim, dim))
        for (rows, cols, signs), amp in zip(struct, t):
            if amp != 0.0:
                np.add.at(T, (rows, cols), signs * amp)
        psi = np.zeros(dim)
        psi[space.ref] = 1.0

        def expm(vec, s):
            out, term = vec.copy(), vec.copy()
            for k in range(1, 2 * nocc + 1):
                term = s * (T @ term) / k
                out += term
            return out

        phi = expm(HN @ expm(psi, 1.0), -1.0)
        return phi[space.ref], np.array([sign * phi[idx] for idx, sign in target])

    energy = 0.0
    for it in range(max_iter):
        e_corr, R = residual(t)
        t_new = t - R / denom
        diis_t.append(t_new)
        diis_e.append(t_new - t)
        if len(diis_t) > 8:
            diis_t.pop(0)
            diis_e.pop(0)
        m = len(diis_t)
        B = -np.ones((m + 1, m + 1))
        B[m, m] = 0.0
        for i in range(m):
            for j in range(m):
                B[i, j] = diis_e[i] @ diis_e[j]
        rhs = np.zeros(m + 1)
        rhs[m] = -1.0
        c = np.linalg.solve(B, rhs)[:m]
        t = sum(ci * ti for ci, ti in zip(c, diis_t))
        if np.max(np.abs(R)) < tol and abs(e_corr - energy) < 1e-13:
            return e_ref, e_corr, it + 1
        energy = e_corr
    raise RuntimeError("brute-force CC did not converge")


def case_from_integrals(name, h, eri, nocc, enuc, extra=None):
    n = h.shape[0]
    hs, V, _ = spin_orbital_hamiltonian(h, eri, nocc)
    space = Space(n, nocc)
    H = build_hamiltonian(space, hs, V)
    e_fci = np.linalg.eigvalsh(H)[0]
    # (pq|rs) packed: p >= q, r >= s, pq >= rs (row-major over pq, then rs)
    packed = [float(eri[p, q, r, s]) for p in range(n) for q in range(p + 1)
              for r in range(p + 1) for s in range(r + 1) if p * (p + 1) // 2 + q >= r * (r + 1) // 2 + s]
    out = {"name": name, "n": n, "nocc": nocc, "nuclear_repulsion": enuc,
           "h": h.ravel().tolist(), "eri_packed": packed, "energies": {}}
    for label, rank in (("CCSD", 2), ("CCSDT", 3), ("CCSDTQ", 4)):
        if rank > 2 * nocc:
            continue
        e_ref, e_corr, iters = solve_cc(space, H, hs, V, rank)
        out["energies"][label] = e_corr
        out["reference_energy"] = e_ref + enuc
        print(f"{name:12s} {label:7s} E_corr = {e_corr:.12f}  ({iters} iterations)")
    out["energies"]["FCI"] = e_fci - out["reference_energy"] + enuc
    print(f"{name:12s} FCI     E_corr = {out['energies']['FCI']:.12f}")
    if extra:
        out.update(extra)
    return out


def random_case(seed=11, nocc=3, nvir=4):
    rng = np.random.default_rng(seed)
    n = nocc + nvir
    h = rng.normal(size=(n, n)) * 0.3
    h = h + h.T + np.diag([-3.0 + 0.5 * p if p < nocc else 1.0 + 0.4 * p for p in range(n)])
    eri = rng.normal(size=(n, n, n, n)) * 0.03
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)
    for p in range(n):
        for q in range(n):
            eri[p, p, q, q] += 0.4  # Coulomb-like diagonal, keeps the reference stable
    return case_from_integrals("random-3o4v", h, eri, nocc, 0.0)


def molecular_case(name, atom, basis):
    from pyscf import gto, scf, cc, ao2mo
    from pyscf.cc import ccsd_t
    mol = gto.M(atom=atom, basis=basis, unit="Bohr", verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    C = mf.mo_coeff
    h = C.T @ mf.get_hcore() @ C
    n = C.shape[1]
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), n)
    mc = cc.CCSD(mf).run(conv_tol=1e-12, conv_tol_normt=1e-9)
    et = ccsd_t.kernel(mc, mc.ao2mo(), verbose=0)
    return case_from_integrals(name, h, eri, mol.nelectron // 2, mol.energy_nuc(),
                               {"pyscf": {"CCSD": mc.e_corr, "CCSD(T)": mc.e_corr + et}})


def main():
    cases = [
        random_case(),
        molecular_case("h2o-sto3g", "O 0 0 0; H 0 1.430429 1.107157; H 0 -1.430429 1.107157", "sto-3g"),
        molecular_case("lih-631g", "Li 0 0 0; H 0 0 3.015", "6-31g"),
    ]
    for c in cases:
        if "pyscf" in c:
            d = c["energies"]["CCSD"] - c["pyscf"]["CCSD"]
            print(f"{c['name']}: brute-force CCSD - PySCF CCSD = {d:.2e}")
            assert abs(d) < 1e-9
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "cc_reference.json").write_text(json.dumps({"cases": cases}, indent=1) + "\n")
    print("wrote cc_reference.json")


if __name__ == "__main__":
    main()
