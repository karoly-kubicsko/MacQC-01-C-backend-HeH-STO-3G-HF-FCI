"""
Minimal-basis STO-nG restricted Hartree--Fock calculation for a
closed-shell, two-electron diatomic molecule.

Python translation of the Fortran program szabo.f, based on Appendix B of
Szabo and Ostlund, Modern Quantum Chemistry.

The default run reproduces the printed report of the original Fortran program
for HeH+ in an STO-3G minimal basis, including Fortran-style D exponents.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan, cos, exp, floor, log10, sin, sqrt
from typing import TextIO
import sys

import numpy as np


# The original Fortran program uses this DATA constant, not a compiler/library pi.
PI = 3.1415926535898


@dataclass
class SCFResult:
    """Container for the final Hartree--Fock result."""

    electronic_energy: float
    total_energy: float
    iterations: int
    density: np.ndarray
    coefficients: np.ndarray
    orbital_energies: np.ndarray
    overlap: np.ndarray
    core_hamiltonian: np.ndarray
    fock: np.ndarray
    g_matrix: np.ndarray
    mulliken_ps: np.ndarray


def fortran_d(value: float, width: int, ndigits: int) -> str:
    """Format a float like a Fortran D descriptor, e.g. D18.10 or D20.12."""
    if value == 0.0:
        mantissa = 0.0
        exponent = 0
    else:
        exponent = floor(log10(abs(value))) + 1
        mantissa = value / (10.0**exponent)

    # Handle rare cases where rounding pushes the mantissa to 1.000...
    rounded = round(abs(mantissa), ndigits)
    if rounded >= 1.0:
        mantissa /= 10.0
        exponent += 1

    number = f"{mantissa:.{ndigits}f}D{exponent:+03d}"
    return f"{number:>{width}}"


def derfother(arg: float) -> float:
    """
    Error-function approximation used by the original Fortran program.

    It is the Abramowitz--Stegun rational approximation, retained here so that
    the printed numerical values match the Fortran output exactly.
    """
    p = 0.3275911
    a = [0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429]
    t = 1.0 / (1.0 + p * arg)
    tn = t
    poly = a[0] * tn
    for coeff in a[1:]:
        tn *= t
        poly += coeff * tn
    return 1.0 - poly * exp(-arg * arg)


class TwoElectronDiatomicHF:
    """
    STO-nG RHF code for a two-electron diatomic molecule.

    This is an algorithm-preserving Python translation of the original Fortran
    code. Matrix indexing is zero-based internally, but the printed report uses
    the original one-based labels.
    """

    def __init__(
        self,
        n: int = 3,
        r: float = 1.4632,
        zeta1: float = 2.0925,
        zeta2: float = 1.24,
        za: float = 2.0,
        zb: float = 1.0,
        iop: int = 2,
        out: TextIO | None = None,
    ) -> None:
        if n not in (1, 2, 3):
            raise ValueError("n must be 1, 2, or 3 for STO-nG.")

        self.n = n
        self.r = float(r)
        self.zeta1 = float(zeta1)
        self.zeta2 = float(zeta2)
        self.za = float(za)
        self.zb = float(zb)
        self.iop = int(iop)
        self.out = sys.stdout if out is None else out

        self.Smat = np.zeros((2, 2), dtype=float)
        self.X = np.zeros((2, 2), dtype=float)
        self.XT = np.zeros((2, 2), dtype=float)
        self.H = np.zeros((2, 2), dtype=float)
        self.F = np.zeros((2, 2), dtype=float)
        self.G = np.zeros((2, 2), dtype=float)
        self.C = np.zeros((2, 2), dtype=float)
        self.Fprime = np.zeros((2, 2), dtype=float)
        self.Cprime = np.zeros((2, 2), dtype=float)
        self.P = np.zeros((2, 2), dtype=float)
        self.oldP = np.zeros((2, 2), dtype=float)
        self.TT = np.zeros((2, 2, 2, 2), dtype=float)
        self.E = np.zeros((2, 2), dtype=float)

        self.S12 = 0.0
        self.T11 = 0.0
        self.T12 = 0.0
        self.T22 = 0.0
        self.V11A = 0.0
        self.V12A = 0.0
        self.V22A = 0.0
        self.V11B = 0.0
        self.V12B = 0.0
        self.V22B = 0.0
        self.V1111 = 0.0
        self.V2111 = 0.0
        self.V2121 = 0.0
        self.V2211 = 0.0
        self.V2221 = 0.0
        self.V2222 = 0.0

    def write(self, text: str) -> None:
        self.out.write(text)

    @staticmethod
    def f0(arg: float) -> float:
        """Boys function F_0(arg), matching the original Fortran implementation."""
        if arg < 1.0e-6:
            return 1.0 - arg / 3.0
        return sqrt(PI / arg) * derfother(sqrt(arg)) / 2.0

    @staticmethod
    def overlap_primitive(a: float, b: float, rab2: float) -> float:
        return (PI / (a + b)) ** 1.5 * exp(-a * b * rab2 / (a + b))

    @staticmethod
    def kinetic_primitive(a: float, b: float, rab2: float) -> float:
        return (
            a
            * b
            / (a + b)
            * (3.0 - 2.0 * a * b * rab2 / (a + b))
            * (PI / (a + b)) ** 1.5
            * exp(-a * b * rab2 / (a + b))
        )

    @classmethod
    def nuclear_attraction_primitive(
        cls, a: float, b: float, rab2: float, rcp2: float, zc: float
    ) -> float:
        value = 2.0 * PI / (a + b) * cls.f0((a + b) * rcp2)
        value *= exp(-a * b * rab2 / (a + b))
        return -value * zc

    @classmethod
    def two_electron_primitive(
        cls,
        a: float,
        b: float,
        c: float,
        d: float,
        rab2: float,
        rcd2: float,
        rpq2: float,
    ) -> float:
        return (
            2.0
            * (PI**2.5)
            / ((a + b) * (c + d) * sqrt(a + b + c + d))
            * cls.f0((a + b) * (c + d) * rpq2 / (a + b + c + d))
            * exp(-a * b * rab2 / (a + b) - c * d * rcd2 / (c + d))
        )

    def calculate_integrals(self) -> None:
        """Calculate all one- and two-electron integrals."""
        coef = np.array(
            [
                [1.000000, 0.678914, 0.444635],
                [0.000000, 0.430129, 0.535328],
                [0.000000, 0.000000, 0.154329],
            ],
            dtype=float,
        )
        expon = np.array(
            [
                [0.270950, 0.151623, 0.109818],
                [0.000000, 0.851819, 0.405771],
                [0.000000, 0.000000, 2.227660],
            ],
            dtype=float,
        )

        r2 = self.r * self.r
        a1 = np.zeros(self.n, dtype=float)
        d1 = np.zeros(self.n, dtype=float)
        a2 = np.zeros(self.n, dtype=float)
        d2 = np.zeros(self.n, dtype=float)

        col = self.n - 1
        for i in range(self.n):
            a1[i] = expon[i, col] * (self.zeta1**2)
            d1[i] = coef[i, col] * ((2.0 * a1[i] / PI) ** 0.75)
            a2[i] = expon[i, col] * (self.zeta2**2)
            d2[i] = coef[i, col] * ((2.0 * a2[i] / PI) ** 0.75)

        self.S12 = self.T11 = self.T12 = self.T22 = 0.0
        self.V11A = self.V12A = self.V22A = 0.0
        self.V11B = self.V12B = self.V22B = 0.0
        self.V1111 = self.V2111 = self.V2121 = 0.0
        self.V2211 = self.V2221 = self.V2222 = 0.0

        for i in range(self.n):
            for j in range(self.n):
                rap = a2[j] * self.r / (a1[i] + a2[j])
                rap2 = rap**2
                rbp2 = (self.r - rap) ** 2

                self.S12 += self.overlap_primitive(a1[i], a2[j], r2) * d1[i] * d2[j]
                self.T11 += self.kinetic_primitive(a1[i], a1[j], 0.0) * d1[i] * d1[j]
                self.T12 += self.kinetic_primitive(a1[i], a2[j], r2) * d1[i] * d2[j]
                self.T22 += self.kinetic_primitive(a2[i], a2[j], 0.0) * d2[i] * d2[j]

                self.V11A += self.nuclear_attraction_primitive(
                    a1[i], a1[j], 0.0, 0.0, self.za
                ) * d1[i] * d1[j]
                self.V12A += self.nuclear_attraction_primitive(
                    a1[i], a2[j], r2, rap2, self.za
                ) * d1[i] * d2[j]
                self.V22A += self.nuclear_attraction_primitive(
                    a2[i], a2[j], 0.0, r2, self.za
                ) * d2[i] * d2[j]
                self.V11B += self.nuclear_attraction_primitive(
                    a1[i], a1[j], 0.0, r2, self.zb
                ) * d1[i] * d1[j]
                self.V12B += self.nuclear_attraction_primitive(
                    a1[i], a2[j], r2, rbp2, self.zb
                ) * d1[i] * d2[j]
                self.V22B += self.nuclear_attraction_primitive(
                    a2[i], a2[j], 0.0, 0.0, self.zb
                ) * d2[i] * d2[j]

        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    for l in range(self.n):
                        rap = a2[i] * self.r / (a2[i] + a1[j])
                        rbp = self.r - rap
                        raq = a2[k] * self.r / (a2[k] + a1[l])
                        rbq = self.r - raq
                        rpq = rap - raq

                        rap2 = rap * rap
                        rbq2 = rbq * rbq
                        rpq2 = rpq * rpq

                        self.V1111 += self.two_electron_primitive(
                            a1[i], a1[j], a1[k], a1[l], 0.0, 0.0, 0.0
                        ) * d1[i] * d1[j] * d1[k] * d1[l]
                        self.V2111 += self.two_electron_primitive(
                            a2[i], a1[j], a1[k], a1[l], r2, 0.0, rap2
                        ) * d2[i] * d1[j] * d1[k] * d1[l]
                        self.V2121 += self.two_electron_primitive(
                            a2[i], a1[j], a2[k], a1[l], r2, r2, rpq2
                        ) * d2[i] * d1[j] * d2[k] * d1[l]
                        self.V2211 += self.two_electron_primitive(
                            a2[i], a2[j], a1[k], a1[l], 0.0, 0.0, r2
                        ) * d2[i] * d2[j] * d1[k] * d1[l]
                        self.V2221 += self.two_electron_primitive(
                            a2[i], a2[j], a2[k], a1[l], 0.0, r2, rbq2
                        ) * d2[i] * d2[j] * d2[k] * d1[l]
                        self.V2222 += self.two_electron_primitive(
                            a2[i], a2[j], a2[k], a2[l], 0.0, 0.0, 0.0
                        ) * d2[i] * d2[j] * d2[k] * d2[l]

        if self.iop != 0:
            self.write(
                f"   STO-{self.n}G FOR ATOMIC NUMBERS {self.za:5.2f} AND {self.zb:5.2f}\n\n\n"
            )
            self.write("   R          ZETA1      ZETA2      S12        T11\n\n")
            self.write(
                f"{self.r:11.6f}{self.zeta1:11.6f}{self.zeta2:11.6f}"
                f"{self.S12:11.6f}{self.T11:11.6f}\n\n\n"
            )
            self.write("   T12        T22        V11A       V12A       V22A\n\n")
            self.write(
                f"{self.T12:11.6f}{self.T22:11.6f}{self.V11A:11.6f}"
                f"{self.V12A:11.6f}{self.V22A:11.6f}\n\n\n"
            )
            self.write("   V11B       V12B       V22B       V1111      V2111\n\n")
            self.write(
                f"{self.V11B:11.6f}{self.V12B:11.6f}{self.V22B:11.6f}"
                f"{self.V1111:11.6f}{self.V2111:11.6f}\n\n\n"
            )
            self.write("   V2121      V2211      V2221      V2222\n\n")
            self.write(
                f"{self.V2121:11.6f}{self.V2211:11.6f}"
                f"{self.V2221:11.6f}{self.V2222:11.6f}\n"
            )

    def collect_matrices(self) -> None:
        """Assemble S, H, X, X^T, and the two-electron integral tensor."""
        self.H[0, 0] = self.T11 + self.V11A + self.V11B
        self.H[0, 1] = self.T12 + self.V12A + self.V12B
        self.H[1, 0] = self.H[0, 1]
        self.H[1, 1] = self.T22 + self.V22A + self.V22B

        self.Smat[0, 0] = 1.0
        self.Smat[0, 1] = self.S12
        self.Smat[1, 0] = self.S12
        self.Smat[1, 1] = 1.0

        self.X[0, 0] = 1.0 / sqrt(2.0 * (1.0 + self.S12))
        self.X[1, 0] = self.X[0, 0]
        self.X[0, 1] = 1.0 / sqrt(2.0 * (1.0 - self.S12))
        self.X[1, 1] = -self.X[0, 1]
        self.XT[:, :] = self.X.T

        self.TT[:, :, :, :] = 0.0
        self.TT[0, 0, 0, 0] = self.V1111
        self.TT[1, 0, 0, 0] = self.V2111
        self.TT[0, 1, 0, 0] = self.V2111
        self.TT[0, 0, 1, 0] = self.V2111
        self.TT[0, 0, 0, 1] = self.V2111
        self.TT[1, 0, 1, 0] = self.V2121
        self.TT[0, 1, 1, 0] = self.V2121
        self.TT[1, 0, 0, 1] = self.V2121
        self.TT[0, 1, 0, 1] = self.V2121
        self.TT[1, 1, 0, 0] = self.V2211
        self.TT[0, 0, 1, 1] = self.V2211
        self.TT[1, 1, 1, 0] = self.V2221
        self.TT[1, 1, 0, 1] = self.V2221
        self.TT[1, 0, 1, 1] = self.V2221
        self.TT[0, 1, 1, 1] = self.V2221
        self.TT[1, 1, 1, 1] = self.V2222

        if self.iop != 0:
            self.matout(self.Smat, "S   ")
            self.matout(self.X, "X   ")
            self.matout(self.H, "H   ")
            self.write("\n\n\n")
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for l in range(2):
                            self.write(
                                f"   ({i + 1:2d}{j + 1:2d}{k + 1:2d}{l + 1:2d} )"
                                f"{self.TT[i, j, k, l]:10.6f}\n"
                            )

    def form_g(self) -> None:
        """Build the two-electron part of the Fock matrix from the density matrix."""
        self.G[:, :] = 0.0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        self.G[i, j] += self.P[k, l] * (
                            self.TT[i, j, k, l] - 0.5 * self.TT[i, l, k, j]
                        )

    @staticmethod
    def diag_2x2(f: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Diagonalize a real symmetric 2x2 matrix using the original formula."""
        c = np.zeros((2, 2), dtype=float)
        e = np.zeros((2, 2), dtype=float)

        if abs(f[0, 0] - f[1, 1]) <= 1.0e-20:
            theta = PI / 4.0
        else:
            theta = 0.5 * atan(2.0 * f[0, 1] / (f[0, 0] - f[1, 1]))

        c[0, 0] = cos(theta)
        c[1, 0] = sin(theta)
        c[0, 1] = sin(theta)
        c[1, 1] = -cos(theta)

        e[0, 0] = (
            f[0, 0] * cos(theta) ** 2
            + f[1, 1] * sin(theta) ** 2
            + f[0, 1] * sin(2.0 * theta)
        )
        e[1, 1] = (
            f[1, 1] * cos(theta) ** 2
            + f[0, 0] * sin(theta) ** 2
            - f[0, 1] * sin(2.0 * theta)
        )
        e[1, 0] = 0.0
        e[0, 1] = 0.0

        if e[1, 1] <= e[0, 0]:
            e[0, 0], e[1, 1] = e[1, 1], e[0, 0]
            c[:, [0, 1]] = c[:, [1, 0]]

        return c, e

    def matout(self, a: np.ndarray, label: str) -> None:
        """Print a matrix using the original Fortran MATOUT layout."""
        label = label[:4].ljust(4)
        self.write(f"\n\n\n    THE {label} ARRAY\n")
        heading = " " * 15 + "".join(
            " " * 10 + f"{j + 1:3d}" + " " * 6 for j in range(a.shape[1])
        )
        self.write(heading.rstrip() + "\n")
        for i in range(a.shape[0]):
            values = "".join(" " + fortran_d(float(a[i, j]), 18, 10) for j in range(a.shape[1]))
            self.write(f"{i + 1:10d}" + " " * 5 + values + "\n")

    def scf(self, crit: float = 1.0e-4, maxit: int = 25) -> SCFResult:
        """Perform the SCF iterations and return the final result."""
        self.P[:, :] = 0.0
        if self.iop >= 2:
            self.matout(self.P, "P   ")

        electronic_energy = 0.0
        iterations = 0

        for iterations in range(1, maxit + 1):
            if self.iop >= 2:
                self.write(f"\n    START OF ITERATION NUMBER = {iterations:2d}\n")

            self.form_g()
            if self.iop >= 2:
                self.matout(self.G, "G   ")

            self.F[:, :] = self.H + self.G
            electronic_energy = 0.5 * float(np.sum(self.P * (self.H + self.F)))

            if self.iop >= 2:
                self.matout(self.F, "F   ")
                self.write(f"\n\n\n    ELECTRONIC ENERGY = {fortran_d(electronic_energy, 20, 12)}\n")

            temp_g = self.F @ self.X
            self.Fprime[:, :] = self.XT @ temp_g
            self.Cprime[:, :], self.E[:, :] = self.diag_2x2(self.Fprime)
            self.C[:, :] = self.X @ self.Cprime

            self.oldP[:, :] = self.P
            self.P[:, :] = 0.0
            for i in range(2):
                for j in range(2):
                    self.P[i, j] = 2.0 * self.C[i, 0] * self.C[j, 0]

            if self.iop >= 2:
                self.matout(self.Fprime, "F'  ")
                self.matout(self.Cprime, "C'  ")
                self.matout(self.E, "E   ")
                self.matout(self.C, "C   ")
                self.matout(self.P, "P   ")

            delta = sqrt(float(np.sum((self.P - self.oldP) ** 2)) / 4.0)
            if self.iop != 0:
                self.write(f"\n    DELTA(CONVERGENCE OF DENSITY MATRIX) = {delta:10.6f}\n\n")

            if delta < crit:
                break
        else:
            self.write("    NO CONVERGENCE IN SCF\n")
            raise RuntimeError("No convergence in SCF")

        total_energy = electronic_energy + self.za * self.zb / self.r

        if self.iop != 0:
            self.write("\n\n    CALCULATION CONVERGED\n\n")
            self.write(f"    ELECTRONIC ENERGY = {fortran_d(electronic_energy, 20, 12)}\n\n")
            self.write(f"    TOTAL ENERGY =      {fortran_d(total_energy, 20, 12)}\n")

        if self.iop == 1:
            self.matout(self.G, "G   ")
            self.matout(self.F, "F   ")
            self.matout(self.E, "E   ")
            self.matout(self.C, "C   ")
            self.matout(self.P, "P   ")

        mulliken_ps = self.P @ self.Smat
        if self.iop != 0:
            self.matout(mulliken_ps, "PS  ")

        return SCFResult(
            electronic_energy=electronic_energy,
            total_energy=total_energy,
            iterations=iterations,
            density=self.P.copy(),
            coefficients=self.C.copy(),
            orbital_energies=np.diag(self.E).copy(),
            overlap=self.Smat.copy(),
            core_hamiltonian=self.H.copy(),
            fock=self.F.copy(),
            g_matrix=self.G.copy(),
            mulliken_ps=mulliken_ps.copy(),
        )

    def run(self) -> SCFResult:
        """Run integral evaluation, matrix assembly, and SCF."""
        self.calculate_integrals()
        self.collect_matrices()
        return self.scf()


def hfcalc(
    iop: int = 2,
    n: int = 3,
    r: float = 1.4632,
    zeta1: float = 2.0925,
    zeta2: float = 1.24,
    za: float = 2.0,
    zb: float = 1.0,
    out: TextIO | None = None,
) -> SCFResult:
    """Convenience function matching the Fortran HFCALC argument list."""
    calculation = TwoElectronDiatomicHF(
        iop=iop, n=n, r=r, zeta1=zeta1, zeta2=zeta2, za=za, zb=zb, out=out
    )
    return calculation.run()


if __name__ == "__main__":
    hfcalc(iop=2, n=3, r=1.4632, zeta1=2.0925, zeta2=1.24, za=2.0, zb=1.0)
