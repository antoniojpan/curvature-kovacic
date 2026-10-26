#!/usr/bin/env python3
"""Circle equation (``ex:nonlinear_proj``, equation ``eq:example_sqrt``).

Equation:

    u' = (1/2) sqrt(1 - u^2).

The curvature is ``kappa = 1/4``.  The paper uses
``delta(x) = sin(x/2)`` and the family
``u(x;s) = sin((x-s)/2)``; here we use the independent solution
``delta(x) = cos(x/2)`` based at ``x0 = 0``.
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
    phi = sp.sqrt(1 - u**2) / 2
    kappa = sp.Rational(1, 4)
    x0, v0 = 0, 0
    x_values = [sp.Rational(1, 2), 1]
    # Generic points used for the identities; they do not lie on the
    # solution family based at (0, 0).
    points = [{x: 1, u: 0}, {x: 1, u: sp.Rational(1, 2)}]

    print("=" * 78)
    print("A1: alg:delta with delta = cos(x/2)")
    print("=" * 78)
    result_a1 = solve_from_linear_solution(
        phi, sp.cos(x / 2), x0=x0, v_value=v0, points=points
    )
    _report(result_a1, phi, kappa, x0, v0, x_values)

    print("\n" + "=" * 78)
    print("A3: alg:kovacic on the nonlinear phi")
    print("=" * 78)
    result_a3 = kovacic_integration(phi, x0=x0, v_value=v0, points=points)
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
    pts = points_from_family(result.family, x_values, v0)
    assert is_zero(result.kappa - kappa, pts)
    assert is_zero(result.c, pts)
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


if __name__ == "__main__":
    main()
