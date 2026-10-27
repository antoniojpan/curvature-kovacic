"""Core geometric objects for the curvature--Kovacic algorithms.

The conventions match the paper:

    A = d/dx + phi(x,u) d/du
    K(x,u) = -d/du(A(phi))
    L = d^2/dx^2 + kappa(x)

When ``K(x,u) = kappa(x)`` is independent of ``u``, one has

    A(phi) = -kappa(x) u + c(x).

The functions in this module use the real symbols :data:`x`, :data:`u` and
:data:`v`.  The first two are the independent and dependent variables of the
scalar ODE; ``v`` is the parameter of the solution family.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from functools import lru_cache

import sympy as sp

#: Independent variable.
x = sp.Symbol("x", real=True)
#: Dependent variable.
u = sp.Symbol("u", real=True)
#: Parameter of the one-parameter family of solutions.
v = sp.Symbol("v", real=True)

Point = Mapping[sp.Symbol, sp.Expr]


def safe_simplify(expr: sp.Expr) -> sp.Expr:
    """Call ``sp.simplify`` but return the input if SymPy raises.

    Large radical expressions occasionally make SymPy's polynomial gcd
    fail.  In that case the caller can still fall back to numerical checks.
    """
    try:
        return sp.simplify(expr)
    except Exception:
        return expr


def A(f: sp.Expr, phi: sp.Expr) -> sp.Expr:
    """Apply the vector field ``A = d/dx + phi d/du`` to ``f``."""
    return sp.diff(f, x) + phi * sp.diff(f, u)


@lru_cache(maxsize=16)
def divergence(phi: sp.Expr) -> sp.Expr:
    """Divergence along solutions, ``A(phi) = phi_x + phi phi_u``, simplified.

    It is simplified once, before any differentiation in ``u``.  For an
    equation written with nested radicals this is what keeps the curvature
    computation tractable: the identity ``A(phi) = -kappa(x) u + c(x)`` is a
    simplification of ``A(phi)`` itself, not of its ``u``-derivative.  The
    result is cached because ``curvature`` and ``c_coefficient`` both need
    it, and the simplification can be expensive.
    """
    return safe_simplify(A(phi, phi))


def curvature(phi: sp.Expr) -> sp.Expr:
    """Gauss curvature ``K(x,u) = -d/du(A(phi))`` of the associated surface."""
    return safe_simplify(-sp.diff(divergence(phi), u))


def c_coefficient(phi: sp.Expr, kappa: sp.Expr) -> sp.Expr:
    """The ``u``-independent coefficient ``c(x) = A(phi) + kappa u``."""
    return safe_simplify(divergence(phi) + kappa * u)


def integral_from(expr: sp.Expr, x0: sp.Expr) -> sp.Expr:
    """Return the definite integral ``int_{x0}^{x} expr(t) dt``."""
    t = sp.Dummy("t")
    integrand = safe_simplify(expr.subs(x, t))
    return sp.integrate(integrand, (t, x0, x))


def numeric_value(expr: sp.Expr, subs: Point) -> complex:
    """Evaluate a scalar expression numerically as a complex double."""
    return complex(sp.N(expr.subs(subs), 30))


def is_zero(expr: sp.Expr, points: Iterable[Point] | None = None,
            tol: float = 1e-7) -> bool:
    """Test whether ``expr`` is zero.

    Symbolic simplification is attempted first.  If that is inconclusive,
    the expression is evaluated numerically at ``points`` (or at one default
    point).  The sample points must assign values to every free symbol that
    occurs in ``expr``.  This is a pragmatic verification helper, not a
    proof: see the README for the limitations of the numerical fallback.
    """
    try:
        expr = sp.simplify(expr)
    except Exception:
        # Some large radical expressions make SymPy's polynomial gcd fail;
        # the numerical fallback below is still available.
        pass
    if expr == 0:
        return True
    for simplifier in (sp.cancel, sp.trigsimp, sp.radsimp, sp.expand, sp.factor):
        try:
            if simplifier(expr) == 0:
                return True
        except Exception:
            pass
    if points is None:
        points = [{x: sp.Rational(2), u: sp.Rational(3, 2)}]
    for point in points:
        try:
            value = numeric_value(expr, point)
        except Exception:
            return False
        if abs(value) > tol:
            return False
    return True


def points_from_family(family: sp.Expr, x_values: Iterable[sp.Expr],
                       v_value: sp.Expr) -> list[Point]:
    """Sample points on a solution family ``u(x;v)``.

    The returned dictionaries are used to evaluate identities in ``(x, u)``
    at points that lie on the returned solution.
    """
    pts: list[Point] = []
    for xv in x_values:
        uv = family.subs({x: xv, v: v_value})
        pts.append({x: xv, u: uv})
    return pts
