#!/usr/bin/env python3
"""Euler--Cauchy example (``ex:euler_cauchy``).

Equation:

    u' = u/(2x) + (3/(2x)) sqrt(u^2 - 4x),   x > 0.

The curvature is ``kappa = -2/x^2`` and ``delta = x^2`` solves
``L(delta) = 0``.  The general solution is

    u(x;s) = s x^2 + 1/(s x).

This script runs ``alg:delta`` and ``alg:kovacic`` on the equation.  See
``ex_particular_solution.py`` for ``alg:particular`` on a simpler affine
Euler--Cauchy equation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from curvature_kovacic import (
    is_zero,
    kovacic_integration,
    points_from_family,
    solve_from_linear_solution,
    u,
    v,
    x,
)


def main() -> None:
    phi = (
        u / (2 * x)
        + sp.Rational(3, 2) / x * sp.sqrt(u**2 - 4 * x)
    )
    kappa = -2 / x**2
    x0, v0 = 2, 5
    x_values = [sp.Rational(5, 2), 3]
    points = [
        {x: 2, u: 5},
        {x: sp.Rational(5, 2), u: 8},
        {x: 3, u: 11},
    ]

    print("=" * 78)
    print("A1: alg:delta with delta = x^2")
    print("=" * 78)
    result_a1 = solve_from_linear_solution(
        phi, x**2, x0=x0, v_value=v0, points=points
    )
    _report(result_a1, phi, kappa, x0, v0, x_values)

    print("\n" + "=" * 78)
    print("A3: alg:kovacic on the nonlinear phi")
    print("=" * 78)
    result_a3 = kovacic_integration(
        phi, x0=x0, v_value=v0, points=points
    )
    _report(result_a3, phi, kappa, x0, v0, x_values)

    print("\nAll checks passed.")


def _report(result, phi, kappa, x0, v0, x_values) -> None:
    print(f"status : {result.status}")
    print(f"kappa  : {sp.simplify(result.kappa)}")
    print(f"c      : {sp.simplify(result.c)}")
    print(f"delta  : {result.delta}")
    print(f"I      : {sp.simplify(result.I)}")
    print(f"family : {sp.simplify(result.family)}")

    assert result.status == "ok", result.message
    assert is_zero(result.kappa - kappa, points_from_family(result.family, x_values, v0))
    assert is_zero(result.c, points_from_family(result.family, x_values, v0))

    pts = points_from_family(result.family, x_values, v0)
    assert is_zero(sp.diff(result.I, u) - result.S, pts)
    assert is_zero(sp.diff(result.I, x) + result.S * phi, pts)
    assert is_zero(result.family.subs({x: x0, v: v0}) - v0, [{x: x0, v: v0}])
    assert is_zero(
        sp.diff(result.family, x).subs({x: x0, v: v0})
        - phi.subs({x: x0, u: v0}),
        [{x: x0, v: v0}],
    )
    for xv in x_values:
        uv = result.family.subs({x: xv, v: v0})
        residual = sp.diff(result.family, x).subs({x: xv, v: v0}) - phi.subs(
            {x: xv, u: uv}
        )
        assert abs(complex(sp.N(residual, 30))) < 1e-6
    I_on = result.I.subs(u, result.family.subs({v: v0}))
    for xv in x_values:
        assert abs(complex(sp.N(sp.diff(I_on, x).subs({x: xv}), 30))) < 1e-6


if __name__ == "__main__":
    main()
