"""Shared test helpers."""

from __future__ import annotations

import sympy as sp

from curvature_kovacic import is_zero, points_from_family, u, v, x


def assert_zero(expr: sp.Expr, points=None, tol: float = 1e-7) -> None:
    simplified = sp.simplify(expr)
    assert is_zero(simplified, points, tol), (
        f"expected zero, got {sp.sstr(simplified)}"
    )


def check_standard_result(result, phi, kappa, c, x0, v0, x_values,
                          tol: float = 1e-7) -> None:
    """Check the identities common to A1, A2 and A3."""
    assert result.status == "ok", result.message
    pts = points_from_family(result.family, x_values, v0)

    assert_zero(result.kappa - kappa, pts, tol)
    assert_zero(result.c - c, pts, tol)

    # dI = S (du - phi dx)
    assert_zero(sp.diff(result.I, u) - result.S, pts, tol)
    assert_zero(sp.diff(result.I, x) + result.S * phi, pts, tol)

    # Initial conditions of the family.
    assert_zero(result.family.subs({x: x0, v: v0}) - v0,
                [{x: x0, v: v0}], tol)
    assert_zero(
        sp.diff(result.family, x).subs({x: x0, v: v0})
        - phi.subs({x: x0, u: v0}),
        [{x: x0, v: v0}],
        tol,
    )

    # The family solves the ODE at the sample points.
    for xv in x_values:
        uv = result.family.subs({x: xv, v: v0})
        residual = sp.diff(result.family, x).subs({x: xv, v: v0}) - phi.subs(
            {x: xv, u: uv}
        )
        value = complex(sp.N(residual, 30))
        assert abs(value) < 1e-6, f"family does not solve the ODE at x={xv}: {value}"

    # The first integral is constant along the family.
    I_on = result.I.subs(u, result.family.subs(v, v0))
    for xv in x_values:
        value = complex(sp.N(sp.diff(I_on, x).subs(x, xv), 30))
        assert abs(value) < 1e-6, f"I is not constant at x={xv}: {value}"


def check_particular_result(result, phi, f, kappa, c, x0, v0, x_values,
                            tol: float = 1e-7) -> None:
    check_standard_result(result, phi, kappa, c, x0, v0, x_values, tol)
    pts = points_from_family(result.family, x_values, v0)
    assert_zero(sp.diff(result.delta, x, 2) + result.kappa * result.delta,
                pts, tol)
    assert_zero(result.delta.subs(x, x0) - 1, [{x: x0}], tol)
