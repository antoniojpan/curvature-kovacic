#!/usr/bin/env python3
"""``alg:particular`` on an affine Euler--Cauchy equation.

Equation:

    u' = 2u/x + x^3.

The curvature is ``kappa = -2/x^2``, the coefficient is ``c = 5x^2`` and
``f(x) = x^4/2`` is a particular solution.  The linear solution built by the
algorithm is ``delta_f(x) = exp(int phi_u(t, f(t)) dt)``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from curvature_kovacic import (
    is_zero,
    points_from_family,
    solve_from_particular,
    u,
    v,
    x,
)


def main() -> None:
    phi = 2 * u / x + x**3
    f = x**4 / 2
    kappa = -2 / x**2
    c = 5 * x**2
    x0, v0 = 1, 2
    x_values = [sp.Rational(3, 2), 2]

    result = solve_from_particular(
        phi, f, x0=x0, v_value=v0, points=[{x: 1, u: 3}, {x: 2, u: 4}]
    )
    print(f"status : {result.status.value}")
    print(f"kappa  : {sp.simplify(result.kappa)}")
    print(f"c      : {sp.simplify(result.c)}")
    print(f"delta_f: {result.delta}")
    print(f"I      : {sp.simplify(result.I)}")
    print(f"family : {sp.simplify(result.family)}")

    assert result.status == "ok", result.message
    pts = points_from_family(result.family, x_values, v0)
    assert is_zero(result.kappa - kappa, pts)
    assert is_zero(result.c - c, pts)
    assert is_zero(result.delta.subs(x, x0) - 1, [{x: x0}])
    assert is_zero(
        sp.diff(result.delta, x, 2) + kappa * result.delta, pts
    )
    assert is_zero(sp.diff(result.I, u) - result.S, pts)
    assert is_zero(sp.diff(result.I, x) + result.S * phi, pts)
    assert is_zero(result.family.subs({x: x0, v: v0}) - v0, [{x: x0, v: v0}])
    for xv in x_values:
        uv = result.family.subs({x: xv, v: v0})
        residual = sp.diff(result.family, x).subs({x: xv, v: v0}) - phi.subs(
            {x: xv, u: uv}
        )
        assert abs(complex(sp.N(residual, 30))) < 1e-6

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
