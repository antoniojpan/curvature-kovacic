"""Checks tied to the formulas displayed in the paper."""

import runpy
from pathlib import Path

import sympy as sp

from curvature_kovacic import (
    KovacicSolver,
    kovacic_integration,
    points_from_family,
    solve_from_linear_solution,
    u,
    v,
    x,
)

from conftest import assert_zero

PHI_GALOIS = u / (4 * x) + sp.sqrt(u**2 - 4 * sp.sqrt(x)) / sp.sqrt(x)
KAPPA_GALOIS = sp.Rational(3, 16) / x**2 - 1 / x
DELTA_GALOIS = x ** sp.Rational(1, 4) * sp.exp(2 * sp.sqrt(x))
# First integral printed in the paper, attached to delta_+.
I_GALOIS_PLUS = (
    x ** sp.Rational(-1, 4)
    * sp.exp(2 * sp.sqrt(x))
    * (u - sp.sqrt(u**2 - 4 * sp.sqrt(x)))
)
# Independent-solution counterpart; I_+ * I_- = 4.
I_GALOIS_MINUS = (
    x ** sp.Rational(-1, 4)
    * sp.exp(-2 * sp.sqrt(x))
    * (u + sp.sqrt(u**2 - 4 * sp.sqrt(x)))
)
S_GALOIS = sp.Rational(1, 2) * sp.exp(-2) * (v + sp.sqrt(v**2 - 4))
F_GALOIS = (
    S_GALOIS * x ** sp.Rational(1, 4) * sp.exp(2 * sp.sqrt(x))
    + (1 / S_GALOIS) * x ** sp.Rational(1, 4) * sp.exp(-2 * sp.sqrt(x))
)
GALOIS_POINTS = [
    {x: 1, u: 3},
    {x: sp.Rational(3, 2), u: 4},
    {x: 2, u: 5},
]


def test_galois_positive_printed_I_and_family():
    result = solve_from_linear_solution(
        PHI_GALOIS, DELTA_GALOIS, x0=1, v_value=3, points=GALOIS_POINTS
    )
    assert result.status == "ok", result.message
    pts = points_from_family(result.family, [sp.Rational(3, 2), 2], 3)

    # alg:delta with delta_+ returns exactly the printed I_+.
    assert_zero(result.I - I_GALOIS_PLUS, pts)
    # The independent solution gives I_-, with I_+ I_- = 4.
    assert_zero(result.I * I_GALOIS_MINUS - 4, pts)
    assert_zero(
        sp.simplify(result.family - F_GALOIS),
        [{x: sp.Rational(3, 2), v: 3}, {x: 2, v: 3}],
    )
    # The printed I_+ is itself a first integral.
    assert_zero(
        sp.diff(I_GALOIS_PLUS, x) + sp.diff(I_GALOIS_PLUS, u) * PHI_GALOIS,
        pts,
    )


def test_galois_positive_kovacic_is_case_two():
    result = KovacicSolver(x).solve(KAPPA_GALOIS)
    assert result.status == "solved"
    assert result.case == 2


def test_galois_positive_algorithm_3():
    result = kovacic_integration(
        PHI_GALOIS, x0=1, v_value=3, points=GALOIS_POINTS
    )
    assert result.status == "ok", result.message
    assert_zero(result.kappa - KAPPA_GALOIS, GALOIS_POINTS)
    assert_zero(result.c, GALOIS_POINTS)

    pts = points_from_family(result.family, [sp.Rational(3, 2), 2], 3)
    # The solver's branch is delta_- = x^{1/4} e^{-2 sqrt(x)},
    # hence alg:delta returns -I_-.
    assert_zero(result.I + I_GALOIS_MINUS, pts)
    assert_zero(
        sp.simplify(result.family - F_GALOIS),
        [{x: sp.Rational(3, 2), v: 3}, {x: 2, v: 3}],
    )


def test_circle_family_is_the_printed_solution():
    phi = sp.sqrt(1 - u**2) / 2
    # Generic points, not lying on the solution family: the zero set of
    # S must be detected from the identity, not from the base point.
    off_family = [{x: 1, u: 0}, {x: 1, u: sp.Rational(1, 2)}]
    result = solve_from_linear_solution(
        phi, sp.cos(x / 2), x0=0, v_value=0, points=off_family
    )
    assert result.status == "ok", result.message
    assert result.branch == "general"
    assert_zero(
        sp.simplify(result.family.subs(v, 0) - sp.sin(x / 2)),
        [{x: sp.Rational(1, 2)}, {x: 1}],
    )


def test_nomaple_script_reproduces_the_printed_example():
    """Run all exact checks in the companion script for ``ex:nomaple``."""
    script = Path(__file__).resolve().parents[1] / "examples" / "ex_nomaple.py"
    namespace = runpy.run_path(str(script))
    namespace["main"]()
