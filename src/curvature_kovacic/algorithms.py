"""The three constructive integration algorithms of the paper.

``solve_from_linear_solution`` implements ``alg:delta``,
``solve_from_particular`` implements ``alg:particular`` and
``kovacic_integration`` implements ``alg:kovacic``.  The latter combines
Kovacic's algorithm with ``alg:delta``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from enum import Enum

import sympy as sp

from .core import (
    Point,
    c_coefficient,
    curvature,
    integral_from,
    is_zero,
    safe_simplify,
    u,
    v,
    x,
)
from .kovacic import KovacicResult, KovacicSolver
from .results import AlgorithmResult, AlgorithmStatus
from .simplification import rational_curvature, refine_delta

_KOVACIC_SOLVER = KovacicSolver(x)
KovacicCallable = Callable[[sp.Expr], KovacicResult]


def _status_value(status: object) -> str:
    if isinstance(status, Enum):
        return str(status.value)
    return str(status)


def full_kovacic(kappa: sp.Expr) -> KovacicResult:
    """Run the package's default Kovacic solver on ``kappa``."""
    return _KOVACIC_SOLVER.solve(kappa)


def solve_from_linear_solution(phi: sp.Expr, delta: sp.Expr, x0: sp.Expr,
                               v_value: sp.Expr = 0,
                               points: Iterable[Point] | None = None
                               ) -> AlgorithmResult:
    """Algorithm ``alg:delta``: integrate from a solution of ``L(delta)=0``.

    Parameters
    ----------
    phi:
        Right-hand side of ``u' = phi(x,u)``.
    delta:
        Proposed nonzero solution of ``L(delta) = delta'' + kappa delta = 0``.
    x0:
        Base point, chosen with ``delta(x0) != 0``.
    v_value:
        Initial value ``u(x0) = v_value`` used to select a member of the
        returned family.
    points:
        Sample points for the symbolic/numerical zero checks.
    """
    K = curvature(phi)
    if not is_zero(sp.diff(K, u), points):
        return AlgorithmResult(
            AlgorithmStatus.OUTSIDE_CLASS,
            kappa=K,
            message="d_u K != 0",
        )

    kappa = rational_curvature(K, points)
    c = c_coefficient(phi, kappa)

    if safe_simplify(delta) == 0:
        return AlgorithmResult(
            AlgorithmStatus.INVALID_AUXILIARY_SOLUTION,
            kappa=kappa,
            c=c,
            message="delta is identically zero",
        )

    L_delta = safe_simplify(sp.diff(delta, x, 2) + kappa * delta)
    if not is_zero(L_delta, points):
        return AlgorithmResult(
            AlgorithmStatus.INVALID_AUXILIARY_SOLUTION,
            kappa=kappa,
            c=c,
            delta=delta,
            message="L(delta) != 0",
        )

    if is_zero(delta.subs(x, x0), [{x: x0}]):
        return AlgorithmResult(
            AlgorithmStatus.INVALID_AUXILIARY_SOLUTION,
            kappa=kappa,
            c=c,
            delta=delta,
            message="delta(x0) == 0: choose another base point",
        )

    S = safe_simplify(sp.diff(delta, x) - sp.diff(phi, u) * delta)
    Q = integral_from(delta * c, x0)

    # Degenerate case S == 0: the ODE is affine in u.
    if is_zero(S, points):
        beta = safe_simplify(phi - (sp.diff(delta, x) / delta) * u)
        if not is_zero(sp.diff(beta, u), points):
            return AlgorithmResult(
                AlgorithmStatus.INVALID_AUXILIARY_SOLUTION,
                kappa=kappa,
                c=c,
                delta=delta,
                S=S,
                message="S == 0 but phi is not affine in u",
            )
        I = safe_simplify(u / delta - integral_from(beta / delta, x0))
        family = safe_simplify(
            delta
            * (
                v / delta.subs(x, x0)
                + integral_from(beta / delta, x0)
            )
        )
        return AlgorithmResult(
            AlgorithmStatus.OK,
            kappa=kappa,
            c=c,
            delta=delta,
            S=safe_simplify(1 / delta),
            S_raw=S,
            Q=Q,
            I=I,
            family=family,
            branch="degenerate",
        )

    # General case S != 0.
    I = safe_simplify(sp.diff(delta, x) * u - delta * phi + Q)

    delta2 = safe_simplify(delta * integral_from(1 / delta**2, x0))
    S2 = safe_simplify(sp.diff(delta2, x) - sp.diff(phi, u) * delta2)
    Q2 = integral_from(delta2 * c, x0)
    I2 = safe_simplify(sp.diff(delta2, x) * u - delta2 * phi + Q2)

    # If S vanishes at the selected base point, use the independent
    # solution for the local implicit solution (its multiplier is nonzero).
    S_at_base = S.subs({x: x0, u: v_value})
    use_delta2 = is_zero(S_at_base)
    I_regular = I2 if use_delta2 else I
    S_regular = S2 if use_delta2 else S

    u_p = safe_simplify(-delta * integral_from(delta2 * c, x0) + delta2 * Q)
    family = safe_simplify(
        u_p
        + (v / delta.subs(x, x0)) * delta
        + (
            delta.subs(x, x0) * phi.subs({x: x0, u: v})
            - sp.diff(delta, x).subs(x, x0) * v
        )
        * delta2
    )

    return AlgorithmResult(
        AlgorithmStatus.OK,
        kappa=kappa,
        c=c,
        delta=delta,
        S=S_regular,
        S_raw=S,
        Q=Q,
        I=I_regular,
        I_alt=I,
        family=family,
        delta2=delta2,
        branch="general",
    )


def solve_from_particular(phi: sp.Expr, f: sp.Expr, x0: sp.Expr,
                          v_value: sp.Expr = 0,
                          points: Iterable[Point] | None = None
                          ) -> AlgorithmResult:
    """Algorithm ``alg:particular``: integrate from a particular solution.

    The linear solution required by ``alg:delta`` is obtained as
    ``delta_f = exp(int phi_u(t, f(t)) dt)``.
    """
    K = curvature(phi)
    if not is_zero(sp.diff(K, u), points):
        return AlgorithmResult(
            AlgorithmStatus.OUTSIDE_CLASS,
            kappa=K,
            message="d_u K != 0",
        )

    kappa = rational_curvature(K, points)
    c = c_coefficient(phi, kappa)

    if not is_zero(sp.diff(f, x) - phi.subs(u, f), points):
        return AlgorithmResult(
            AlgorithmStatus.INVALID_PARTICULAR_SOLUTION,
            kappa=kappa,
            c=c,
            message="f' != phi(x,f)",
        )

    p = safe_simplify(sp.diff(phi, u).subs(u, f))
    delta_f = safe_simplify(sp.exp(integral_from(p, x0)))

    if not is_zero(sp.diff(delta_f, x, 2) + kappa * delta_f, points):
        return AlgorithmResult(
            AlgorithmStatus.INVALID_AUXILIARY_SOLUTION,
            kappa=kappa,
            c=c,
            delta=delta_f,
            message="delta_f does not solve L(delta_f)=0",
        )

    base = solve_from_linear_solution(phi, delta_f, x0, v_value, points)
    if base.status != AlgorithmStatus.OK:
        base.delta = delta_f
        return base

    delta_f2 = safe_simplify(delta_f * integral_from(1 / delta_f**2, x0))
    c_at = sp.diff(phi, u).subs({x: x0, u: f.subs(x, x0)})
    family = safe_simplify(
        f
        + (v - f.subs(x, x0)) * delta_f
        + (
            phi.subs({x: x0, u: v})
            - sp.diff(f, x).subs(x, x0)
            - c_at * (v - f.subs(x, x0))
        )
        * delta_f2
    )

    base.delta = delta_f
    base.delta2 = delta_f2
    base.family = family
    base.branch = "particular"
    return base


def kovacic_integration(phi: sp.Expr, x0: sp.Expr, v_value: sp.Expr = 0,
                        points: Iterable[Point] | None = None,
                        kovacic: KovacicSolver | KovacicCallable | None = None
                        ) -> AlgorithmResult:
    """Algorithm ``alg:kovacic``: Kovacic's algorithm followed by A1.

    ``kovacic`` may be a :class:`KovacicSolver`, a callable returning a
    :class:`KovacicResult`, or ``None`` to use the package solver.  A solver
    limitation is reported as ``SOLVER_LIMITED``; only a genuine Case 4
    outcome is reported as ``OBSTRUCTION``.
    """
    K = curvature(phi)
    if not is_zero(sp.diff(K, u), points):
        return AlgorithmResult(
            AlgorithmStatus.OUTSIDE_CLASS,
            kappa=K,
            message="d_u K != 0",
        )

    kappa = rational_curvature(K, points)
    c = c_coefficient(phi, kappa)

    if not kappa.is_rational_function(x):
        return AlgorithmResult(
            AlgorithmStatus.NOT_RATIONAL,
            kappa=kappa,
            c=c,
            message="Kovacic's algorithm requires kappa in C(x)",
        )

    solver = _KOVACIC_SOLVER if kovacic is None else kovacic
    if isinstance(solver, KovacicSolver):
        kovacic_result = solver.solve(kappa)
    else:
        kovacic_result = solver(kappa)

    status = _status_value(kovacic_result.status)
    if status in ("solved", "success"):
        delta = kovacic_result.delta
        if delta is None:
            return AlgorithmResult(
                AlgorithmStatus.SOLVER_ERROR,
                kappa=kappa,
                c=c,
                message="solver returned SOLVED without a delta",
            )
        delta = refine_delta(delta, kappa, x0, points)
        result = solve_from_linear_solution(phi, delta, x0, v_value, points)
        result.message = "kovacic_success"
        return result
    if status in ("obstruction",):
        return AlgorithmResult(
            AlgorithmStatus.OBSTRUCTION,
            kappa=kappa,
            c=c,
            message="at most one Liouvillian solution over C(x)",
        )
    if status in ("solver_limited", "not_found"):
        return AlgorithmResult(
            AlgorithmStatus.SOLVER_LIMITED,
            kappa=kappa,
            c=c,
            message=kovacic_result.message
            or "provided Kovacic solver returned no solution",
        )
    return AlgorithmResult(
        AlgorithmStatus.SOLVER_ERROR,
        kappa=kappa,
        c=c,
        message=kovacic_result.message,
    )
