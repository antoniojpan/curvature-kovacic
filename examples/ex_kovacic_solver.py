#!/usr/bin/env python3
"""Diagnostics for the self-contained Kovacic solver.

This script runs Cases 1--3 and reports the Case 4 obstruction on a few
classical equations.  It is a quick way to inspect the solver independently
of the nonlinear algorithms.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from curvature_kovacic import KovacicSolver, is_zero, x


def main() -> None:
    solver = KovacicSolver(x)
    examples = [
        ("y'' + y = 0", sp.Integer(1)),
        ("y'' - y = 0", sp.Integer(-1)),
        ("Euler-Cauchy, kappa = -2/x^2", -2 / x**2),
        ("Kovacic Example 1, kappa = 3/(16x^2) - 1/x",
         sp.Rational(3, 16) / x**2 - 1 / x),
        ("Case 2, kappa = (3 + 8x)/(16x^2)",
         (3 + 8 * x) / (16 * x**2)),
        ("Case 3, kappa = (x - 4)/(4x(x - 1)^2)",
         (x - 4) / (4 * x * (x - 1) ** 2)),
        ("Airy, kappa = x", x),
        ("Bessel nu=1, kappa = 1 - 3/(4x^2)",
         1 - sp.Rational(3, 4) / x**2),
    ]

    for label, kappa in examples:
        result = solver.solve(kappa)
        print("=" * 78)
        print(f"{label}")
        print(f"  status : {result.status.value}")
        print(f"  case   : {result.case}")
        if result.delta is not None:
            print(f"  delta  : {result.delta}")
            residual = sp.diff(result.delta, x, 2) + kappa * result.delta
            print(f"  L(delta) = 0 ? {is_zero(residual, [{x: 1}])}")
        if result.message:
            print(f"  message: {result.message}")


if __name__ == "__main__":
    main()
