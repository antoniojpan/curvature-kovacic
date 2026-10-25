#!/usr/bin/env python3
"""Kovacic's Example 1, Section 4.2 (``ex:galois_positive``).

Equation:

    u' = u/(4x) + sqrt(u^2 - 4 sqrt(x))/sqrt(x),   x > 0,
    u > 2 x^{1/4}.

The paper displays the Liouvillian solution

    delta(x) = x^{1/4} exp(2 sqrt(x))

of ``L(y) = y'' + (3/(16x^2) - 1/x) y = 0``, the first integral

    I(x,u) = (1/2) x^{-1/4} exp(-2 sqrt(x)) (u + sqrt(u^2 - 4 sqrt(x)))

and the family

    f(x;s) = s x^{1/4} exp(2 sqrt(x))
             + s^{-1} x^{1/4} exp(-2 sqrt(x)).

This script runs the paper's pipeline (``alg:kovacic`` on the genuine
nonlinear ``phi``) and checks the displayed formulas up to the standard
freedom in the choice of first integral.
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
    rational_curvature,
    solve_from_linear_solution,
    u,
    v,
    x,
)


def main() -> None:
    phi = u / (4 * x) + sp.sqrt(u**2 - 4 * sp.sqrt(x)) / sp.sqrt(x)
    kappa = sp.Rational(3, 16) / x**2 - 1 / x

    delta_paper = x ** sp.Rational(1, 4) * sp.exp(2 * sp.sqrt(x))
    I_paper = (
        sp.Rational(1, 2)
        * x ** sp.Rational(-1, 4)
        * sp.exp(-2 * sp.sqrt(x))
        * (u + sp.sqrt(u**2 - 4 * sp.sqrt(x)))
    )
    s_paper = sp.Rational(1, 2) * sp.exp(-2) * (v + sp.sqrt(v**2 - 4))
    f_paper = (
        s_paper * x ** sp.Rational(1, 4) * sp.exp(2 * sp.sqrt(x))
        + (1 / s_paper) * x ** sp.Rational(1, 4) * sp.exp(-2 * sp.sqrt(x))
    )

    points = [
        {x: 1, u: 3},
        {x: sp.Rational(3, 2), u: 4},
        {x: 2, u: 5},
    ]

    print("=" * 78)
    print("A1: alg:delta with the delta displayed in the paper")
    print("=" * 78)
    result_a1 = solve_from_linear_solution(
        phi, delta_paper, x0=1, v_value=3, points=points
    )
    print(f"status : {result_a1.status.value}")
    print(f"kappa  : {rational_curvature(result_a1.kappa, points)}")
    print(f"c      : {sp.simplify(result_a1.c)}")
    print(f"delta  : {result_a1.delta}")
    print(f"I      : {sp.simplify(result_a1.I)}")
    print(f"family : {sp.simplify(result_a1.family)}")

    assert result_a1.status == "ok", result_a1.message
    assert is_zero(result_a1.kappa - kappa, points)
    assert is_zero(result_a1.c, points)

    pts_a1 = points_from_family(result_a1.family, [sp.Rational(3, 2), 2], 3)
    # With delta_+ = x^{1/4} e^{2 sqrt(x)}, alg:delta returns 2 / I_paper.
    assert is_zero(result_a1.I * I_paper - 2, pts_a1)
    assert is_zero(
        sp.simplify(result_a1.family - f_paper),
        [{x: sp.Rational(3, 2), v: 3}, {x: 2, v: 3}],
    )
    assert is_zero(
        sp.diff(I_paper, x) + sp.diff(I_paper, u) * phi, pts_a1
    )

    print("\n" + "=" * 78)
    print("A3: alg:kovacic on the nonlinear phi of the paper")
    print("=" * 78)
    result_a3 = kovacic_integration(
        phi, x0=1, v_value=3, points=points
    )
    print(f"status : {result_a3.status.value}")
    print(f"kappa  : {rational_curvature(result_a3.kappa, points)}")
    print(f"c      : {sp.simplify(result_a3.c)}")
    print(f"delta  : {result_a3.delta}")
    print(f"I      : {sp.simplify(result_a3.I)}")
    print(f"family : {sp.simplify(result_a3.family)}")

    assert result_a3.status == "ok", result_a3.message
    assert is_zero(result_a3.kappa - kappa, points)
    assert is_zero(result_a3.c, points)

    pts_a3 = points_from_family(result_a3.family, [sp.Rational(3, 2), 2], 3)
    # The solver's branch is delta_- = x^{1/4} e^{-2 sqrt(x)}, so the
    # algorithm returns -2 I_paper.
    assert is_zero(result_a3.I + 2 * I_paper, pts_a3)
    assert is_zero(
        sp.simplify(result_a3.family - f_paper),
        [{x: sp.Rational(3, 2), v: 3}, {x: 2, v: 3}],
    )

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
