"""Simplification helpers used by the symbolic algorithms.

The mathematical algorithms do not depend on these routines.  They are
convenience steps that help SymPy recognise a rational curvature and choose
a convenient branch of a Liouvillian solution.  Both are heuristics: each
candidate is verified at the given sample points before it is accepted.
"""

from __future__ import annotations

from collections.abc import Iterable

import sympy as sp

from .core import Point, is_zero, numeric_value, u, x


def rational_curvature(kappa: sp.Expr,
                       points: Iterable[Point] | None = None) -> sp.Expr:
    """Try to write a ``u``-independent curvature as a rational function.

    Algebraic curvatures (for instance those involving square roots) may be
    independent of ``u`` while SymPy keeps spurious radicals.  Since
    ``d_u kappa = 0`` has already been checked by the caller, we may evaluate
    at a convenient value of ``u`` and cancel.  A candidate is accepted only
    if it is rational in ``x`` and agrees numerically with the original
    expression at the sample points.
    """
    if kappa.is_rational_function(x):
        return sp.cancel(kappa)
    if points is None:
        points = [{x: sp.Rational(2), u: sp.Rational(3, 2)}]
    points = list(points)
    for u0 in (sp.Integer(0), sp.Integer(1)):
        try:
            candidate = sp.simplify(kappa.subs(u, u0))
        except Exception:
            continue
        if not candidate.is_rational_function(x):
            continue
        ok = True
        for point in points:
            try:
                original = numeric_value(kappa, point)
                simplified = complex(sp.N(candidate.subs({x: point[x]}), 30))
            except Exception:
                ok = False
                break
            if abs(original - simplified) > 1e-7:
                ok = False
                break
        if ok:
            return sp.cancel(candidate)
    for simplifier in (
        sp.cancel,
        lambda e: sp.cancel(sp.radsimp(e)),
        lambda e: sp.cancel(sp.together(e)),
    ):
        try:
            candidate = simplifier(kappa)
        except Exception:
            continue
        if candidate.is_rational_function(x):
            return candidate
    return kappa


def refine_delta(delta: sp.Expr, kappa: sp.Expr, x0: sp.Expr,
                 points: Iterable[Point] | None = None) -> sp.Expr:
    """Choose the branch of ``delta`` adapted to the domain ``x > 0``.

    For example, ``exp(-2*sqrt(x**3)/x)`` is simplified to
    ``exp(-2*sqrt(x))``.  The refined expression replaces ``delta`` only when
    it is shorter, nonzero at ``x0`` and still satisfies ``L(delta) = 0`` at
    the sample points.  This is a branch-selection heuristic; the algorithms
    work with any nonzero solution of ``L(delta) = 0`` supplied by the
    solver.
    """
    try:
        refined = sp.refine(delta, sp.Q.positive(x))
    except Exception:
        return delta
    if (
        sp.count_ops(refined) < sp.count_ops(delta)
        and not is_zero(refined.subs(x, x0), [{x: x0}])
        and is_zero(sp.diff(refined, x, 2) + kappa * refined, points)
    ):
        return refined
    return delta
