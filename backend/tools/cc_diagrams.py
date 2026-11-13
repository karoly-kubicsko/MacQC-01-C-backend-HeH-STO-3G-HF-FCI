#!/usr/bin/env python3
"""Spin-orbital coupled-cluster equations up to CCSDTQ for the ΨLab C++ engine.

Usage:  python3 backend/tools/cc_diagrams.py            (verify and write backend/include/macqc/cc_diagrams.inc)
        python3 backend/tools/cc_diagrams.py --check    (verify only)

Theory
------
Residuals  R_m(A;I) = <Phi_I^A| (H_N e^T)_C |0>,  m = 0 (energy) .. 4, with the
normal-ordered Hamiltonian H_N = F_N + V_N,

    F_N = sum_pq f_pq {p+ q},            V_N = 1/4 sum_pqrs <pq||rs> {p+ q+ s r},
    T_n = sum_{a1<..<an, i1<..<in} t(a1..an; i1..in) a1+ .. an+ in .. i1.

Every connected term is one Hamiltonian vertex (F or V) and k <= 4 cluster vertices
T_n, each joined to the Hamiltonian by cp particle and ch hole lines. A diagram is
fixed by the vertex type, x = number of external particle lines leaving the vertex
(virtual creation slots), y = number of particle lines entering it (virtual
annihilation slots) and the list of (n, cp, ch). The engine evaluates a diagram
from the Hamiltonian outwards, one cluster vertex at a time, with packed
antisymmetric tensors:

  * index layout: vertex creation slots [x external particles, contracted holes in
    cluster order], annihilation slots [contracted particles in cluster order,
    external holes]; cluster vertex particles [contracted, external], holes
    [contracted, external]; target particles [vertex, T_1, T_2, ...], holes alike;
  * contraction = sum over sorted index sets; external groups are merged with
    the sign of the shuffle (sorting the concatenation).

In this convention the coefficient of a diagram is s / n_v, where n_v is the number
of permutations of identical cluster vertices and s = +-1. The signs are not
derived by hand: for every diagram the unrestricted tensor expression is evaluated
with numpy and the signs are fitted to a brute-force evaluation of
<mu| e^{-T} H e^{T} |0> in a determinant space (random Fock matrix, integrals and
amplitudes); the fit must be exact to rounding, otherwise the script stops.

The cluster vertices of each diagram are ordered to minimise the cost (operation
count and intermediate size for 10 occupied and 80 virtual spin orbitals).
Diagonal Fock elements f_ii and f_aa are not in the list: the engine adds
(sum f_aa - sum f_ii) t(A;I) itself.
"""
import itertools
import math
import pathlib
import sys

import numpy as np

MAX_RANK = 4
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "backend" / "include" / "macqc" / "cc_diagrams.inc"


# ----------------------------------------------------------------------------------
# Diagram enumeration
# ----------------------------------------------------------------------------------

class Diagram:
    def __init__(self, vertex, x, y, clusters):
        self.vertex = vertex          # 'F' or 'V'
        self.x = x                    # external particle lines of the vertex
        self.y = y                    # particle lines into the vertex (contracted)
        self.clusters = list(clusters)  # ordered list of (n, cp, ch)
        c = 1 if vertex == 'F' else 2
        self.c = c
        self.cp = y                   # contracted particle lines of the vertex
        self.ch = c - x               # contracted hole lines of the vertex
        self.ext_p = x + sum(n - cp for n, cp, ch in self.clusters)
        self.ext_h = (c - y) + sum(n - ch for n, cp, ch in self.clusters)
        self.sign = 0

    @property
    def target(self):
        return self.ext_p

    @property
    def symmetry(self):
        counts = {}
        for t in self.clusters:
            counts[t] = counts.get(t, 0) + 1
        s = 1
        for k in counts.values():
            s *= math.factorial(k)
        return s

    def key(self):
        return (self.target, self.vertex, self.x, self.y, tuple(sorted(self.clusters)))

    def __repr__(self):
        ts = " ".join(f"T{n}({cp},{ch})" for n, cp, ch in self.clusters)
        return f"R{self.target} <- {self.vertex}[x={self.x},y={self.y}] {ts}"


def enumerate_diagrams(max_rank=MAX_RANK):
    out = []
    for vertex, c in (('F', 1), ('V', 2)):
        for x in range(c + 1):
            for y in range(c + 1):
                cp_h, ch_h = y, c - x
                types = [(n, cp, ch) for n in range(1, max_rank + 1)
                         for cp in range(min(n, cp_h) + 1) for ch in range(min(n, ch_h) + 1) if cp + ch >= 1]
                for k in range(0, cp_h + ch_h + 1):
                    for combo in itertools.combinations_with_replacement(types, k):
                        if sum(t[1] for t in combo) != cp_h or sum(t[2] for t in combo) != ch_h:
                            continue
                        d = Diagram(vertex, x, y, combo)
                        if d.ext_p != d.ext_h or d.ext_p > max_rank:
                            continue
                        out.append(d)
    return out


# ----------------------------------------------------------------------------------
# Cost model and ordering of the cluster vertices
# ----------------------------------------------------------------------------------

def comb(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def order_cost(d, order, no=10, nv=80):
    """Operation count and largest intermediate for evaluating the clusters in `order`."""
    op, oh = d.cp, d.ch
    ep, eh = d.x, d.c - d.y
    flops = 0.0
    peak = 0.0
    for (n, cp, ch) in order:
        op_rest, oh_rest = op - cp, oh - ch
        rows = comb(nv, op_rest) * comb(no, oh_rest) * comb(nv, ep) * comb(no, eh)
        k = comb(nv, cp) * comb(no, ch)
        cols = comb(nv, n - cp) * comb(no, n - ch)
        flops += rows * k * cols
        op, oh = op_rest, oh_rest
        ep, eh = ep + n - cp, eh + n - ch
        size = comb(nv, op) * comb(no, oh) * comb(nv, ep) * comb(no, eh)
        peak = max(peak, size)
    return flops, peak


def choose_order(d):
    best = None
    for perm in set(itertools.permutations(d.clusters)):
        flops, peak = order_cost(d, perm)
        key = (flops + 4.0 * peak, perm)
        if best is None or key < best:
            best = key
    d.clusters = list(best[1])


# ----------------------------------------------------------------------------------
# Unrestricted tensor expression of a diagram (numpy, dense, no spin)
# ----------------------------------------------------------------------------------

def perm_sign(p):
    p = list(p)
    s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def sorted_tuples(m, no, nv):
    return [(A, I) for A in itertools.combinations(range(nv), m) for I in itertools.combinations(range(no), m)]


_PERMS = {}


def antisymmetrized_values(z, m, no, nv):
    """Full antisymmetriser over the m particle and m hole axes of z (no 1/m! factors),
    evaluated at all sorted index tuples (A, I); returns a vector."""
    if m == 0:
        return np.array([float(z)])
    if m not in _PERMS:
        perms = [(p, perm_sign(p)) for p in itertools.permutations(range(m))]
        _PERMS[m] = perms
    perms = _PERMS[m]
    tuples = sorted_tuples(m, no, nv)
    out = np.zeros(len(tuples))
    for pa, sa in perms:
        for pi, si in perms:
            idx = tuple(np.array([[A[k] for k in pa] for A, I in tuples]).T) + \
                  tuple(np.array([[I[k] for k in pi] for A, I in tuples]).T)
            out += sa * si * z[idx]
    return out


def evaluate_unrestricted(d, f, V, t, no, nv):
    """A[ sum_internal H T_1 .. T_k ] with the layout of the module docstring, at the
    sorted index tuples of the target (vector in the order of sorted_tuples)."""
    letters = iter("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
    o, v = slice(0, no), slice(no, no + nv)
    ext_p_vertex = [next(letters) for _ in range(d.x)]
    ext_h_vertex = [next(letters) for _ in range(d.c - d.y)]
    cont_p = [[next(letters) for _ in range(cp)] for (n, cp, ch) in d.clusters]
    cont_h = [[next(letters) for _ in range(ch)] for (n, cp, ch) in d.clusters]
    ext_p = [[next(letters) for _ in range(n - cp)] for (n, cp, ch) in d.clusters]
    ext_h = [[next(letters) for _ in range(n - ch)] for (n, cp, ch) in d.clusters]
    cre = ext_p_vertex + sum(cont_h, [])
    ann = sum(cont_p, []) + ext_h_vertex
    space = lambda idx: v if idx in ext_p_vertex or any(idx in g for g in cont_p) else o  # noqa: E731
    sl = tuple(space(i) for i in cre + ann)
    H = (f if d.vertex == 'F' else V)[sl]
    operands = [H]
    subs = ["".join(cre + ann)]
    for j, (n, cp, ch) in enumerate(d.clusters):
        operands.append(t[n])
        subs.append("".join(cont_p[j] + ext_p[j] + cont_h[j] + ext_h[j]))
    out_p = ext_p_vertex + sum(ext_p, [])
    out_h = ext_h_vertex + sum(ext_h, [])
    z = np.einsum(",".join(subs) + "->" + "".join(out_p + out_h), *operands, optimize=True)
    return antisymmetrized_values(z, len(out_p), no, nv)


def packing_factor(d):
    """K_D = product of factorials of all line groups (unrestricted / packed)."""
    k = math.factorial(d.x) * math.factorial(d.c - d.y)
    for n, cp, ch in d.clusters:
        k *= math.factorial(cp) * math.factorial(ch) * math.factorial(n - cp) * math.factorial(n - ch)
    return k


# ----------------------------------------------------------------------------------
# Brute force: <mu| e^{-T} H_N e^{T} |0> in the determinant space
# ----------------------------------------------------------------------------------

class Determinants:
    def __init__(self, no, nv):
        self.no, self.nv = no, nv
        self.N = no + nv
        self.dets = [sum(1 << i for i in c) for c in itertools.combinations(range(self.N), no)]
        self.index = {d: i for i, d in enumerate(self.dets)}
        self.ref = self.index[(1 << no) - 1]

    @staticmethod
    def apply(ops, det):
        """Apply operators [(orbital, is_creation), ...] right to left; returns (det, sign) or None."""
        sign = 1
        for orb, creation in reversed(ops):
            bit = 1 << orb
            occupied = bool(det & bit)
            if creation == occupied:
                return None
            if bin(det & (bit - 1)).count("1") % 2:
                sign = -sign
            det ^= bit
        return det, sign

    def excitations(self, max_rank):
        """List of (rank, A, I) with sorted A (virtual) and I (occupied), and their phases."""
        out = []
        for n in range(1, max_rank + 1):
            for A in itertools.combinations(range(self.no, self.N), n):
                for I in itertools.combinations(range(self.no), n):
                    out.append((n, A, I))
        return out

    def excitation_ops(self, A, I):
        return [(a, True) for a in A] + [(i, False) for i in reversed(I)]


def brute_force_residuals(dets, f, V, amplitudes, max_rank):
    """R_m(A;I) for all excitations up to max_rank (dict), H_N = F_N + V_N with Fock f."""
    no, N = dets.no, dets.N
    occ = range(no)
    # Plain-operator Hamiltonian with h chosen so that its Fock matrix is f.
    h = f - np.einsum("piqi->pq", V[:, :no, :, :no])
    ndet = len(dets.dets)
    Hm = np.zeros((ndet, ndet))
    for col, det in enumerate(dets.dets):
        for p in range(N):
            for q in range(N):
                if h[p, q] != 0.0:
                    r = dets.apply([(p, True), (q, False)], det)
                    if r:
                        Hm[dets.index[r[0]], col] += r[1] * h[p, q]
        occs = [i for i in range(N) if det >> i & 1]
        for r_, s in itertools.combinations(occs, 2):
            for p, q in itertools.combinations(range(N), 2):
                res = dets.apply([(p, True), (q, True), (s, False), (r_, False)], det)
                if res:
                    Hm[dets.index[res[0]], col] += res[1] * V[p, q, r_, s]
    e_ref = Hm[dets.ref, dets.ref]
    Hm -= e_ref * np.eye(ndet)
    Tm = np.zeros((ndet, ndet))
    for (n, A, I), amp in amplitudes.items():
        ops = dets.excitation_ops(A, I)
        for col, det in enumerate(dets.dets):
            r = dets.apply(ops, det)
            if r:
                Tm[dets.index[r[0]], col] += r[1] * amp
    psi = np.zeros(ndet)
    psi[dets.ref] = 1.0
    def expm_apply(M, vec, s):
        out, term = vec.copy(), vec.copy()
        for k in range(1, no + 2):
            term = s * (M @ term) / k
            out += term
        return out
    phi = expm_apply(Tm, Hm @ expm_apply(Tm, psi, 1.0), -1.0)
    res = {0: phi[dets.ref]}
    for n in range(1, max_rank + 1):
        for A in itertools.combinations(range(no, N), n):
            for I in itertools.combinations(occ, n):
                det, sign = dets.apply(dets.excitation_ops(A, I), (1 << no) - 1)
                res[(n, A, I)] = sign * phi[dets.index[det]]
    return res


def random_problem(rng, no, nv, max_rank):
    N = no + nv
    f = rng.normal(size=(N, N)) * 0.3
    f = f + f.T
    V = rng.normal(size=(N, N, N, N)) * 0.2
    V = V - V.transpose(1, 0, 2, 3)
    V = V - V.transpose(0, 1, 3, 2)
    V = V + V.transpose(2, 3, 0, 1)   # real Hermitian symmetry (not needed, but physical)
    amplitudes = {}
    t = {}
    for n in range(1, max_rank + 1):
        dense = np.zeros((nv,) * n + (no,) * n)
        for A in itertools.combinations(range(nv), n):
            for I in itertools.combinations(range(no), n):
                val = rng.normal() * (0.5 / n)
                amplitudes[(n, tuple(a + no for a in A), I)] = val
                for pa in itertools.permutations(range(n)):
                    for pi in itertools.permutations(range(n)):
                        dense[tuple(A[k] for k in pa) + tuple(I[k] for k in pi)] = perm_sign(pa) * perm_sign(pi) * val
        t[n] = dense
    return f, V, t, amplitudes


def fit_signs(diagrams, no=6, nv=7, samples=3, seed=7):
    rng = np.random.default_rng(seed)
    dets = Determinants(no, nv)
    rows = {m: [] for m in range(MAX_RANK + 1)}
    rhs = {m: [] for m in range(MAX_RANK + 1)}
    for _ in range(samples):
        f, V, t, amplitudes = random_problem(rng, no, nv, MAX_RANK)
        brute = brute_force_residuals(dets, f, V, amplitudes, MAX_RANK)
        # Generic F terms use f without its occupied-occupied and virtual-virtual diagonal.
        f_off = f.copy()
        for p in range(no + nv):
            f_off[p, p] = 0.0
        values = {}
        for d in diagrams:
            values[id(d)] = evaluate_unrestricted(d, f_off, V, t, no, nv) / (d.symmetry * packing_factor(d))
        diag_f = np.diag(f)
        for m in range(MAX_RANK + 1):
            ds = [d for d in diagrams if d.target == m]
            cols = np.array([values[id(d)] for d in ds]).T
            if m == 0:
                rows[0].extend(cols.tolist())
                rhs[0].append(brute[0])
                continue
            rows[m].extend(cols.tolist())
            for A, I in sorted_tuples(m, no, nv):
                amp = amplitudes[(m, tuple(a + no for a in A), I)]
                denom = sum(diag_f[a + no] for a in A) - sum(diag_f[i] for i in I)
                rhs[m].append(brute[(m, tuple(a + no for a in A), I)] - denom * amp)
    worst = 0.0
    for m in range(MAX_RANK + 1):
        ds = [d for d in diagrams if d.target == m]
        M = np.array(rows[m])
        b = np.array(rhs[m])
        coef, *_ = np.linalg.lstsq(M, b, rcond=None)
        for d, c in zip(ds, coef):
            if abs(abs(c) - 1.0) > 1e-8:
                raise SystemExit(f"coefficient {c:.6f} is not +-1 for {d}")
            d.sign = 1 if c > 0 else -1
        err = np.max(np.abs(M @ np.array([d.sign for d in ds], dtype=float) - b))
        worst = max(worst, err)
        print(f"R{m}: {len(ds):3d} diagrams, {M.shape[0]:5d} equations, max error {err:.2e}")
    if worst > 1e-10:
        raise SystemExit("the diagram expansion does not reproduce the brute-force residuals")


# ----------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------

def write_include(diagrams):
    lines = [
        "// Generated by backend/tools/cc_diagrams.py -- do not edit.",
        "// Connected coupled-cluster diagrams up to CCSDTQ (spin orbitals), verified against",
        "// a brute-force evaluation of <mu|exp(-T) H exp(T)|0>. Fields: target excitation level,",
        "// vertex ('F' Fock, 'V' two-electron), x, y, coefficient (+-1/n_v), number of cluster",
        "// vertices, then (rank, contracted particle lines, contracted hole lines) per vertex in",
        "// evaluation order.",
    ]
    for d in sorted(diagrams, key=lambda d: (d.target, d.vertex, d.x, d.y, len(d.clusters), d.clusters)):
        cl = ", ".join(f"{{{n}, {cp}, {ch}}}" for n, cp, ch in d.clusters)
        cl = cl if cl else ""
        coef = f"{d.sign}.0 / {d.symmetry}.0" if d.symmetry != 1 else f"{d.sign}.0"
        lines.append(f"{{{d.target}, '{d.vertex}', {d.x}, {d.y}, {coef}, {len(d.clusters)}, {{{cl}}}}},")
    OUTPUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)} ({len(diagrams)} diagrams)")


def main():
    diagrams = enumerate_diagrams()
    for d in diagrams:
        choose_order(d)
    print(f"{len(diagrams)} diagrams up to CCSDTQ")
    fit_signs(diagrams)
    if "--check" not in sys.argv:
        write_include(diagrams)


if __name__ == "__main__":
    main()
