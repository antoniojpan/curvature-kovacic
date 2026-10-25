"""A self-contained implementation of Kovacic's algorithm.

The solver treats ``y'' + kappa(x) y = 0`` with ``kappa`` rational and
implements Cases 1, 2 and 3 of Kovacic (1986).  When none of the three cases
applies, or when all applicable cases yield no solution, the equation has no
nonzero Liouvillian solution over ``C(x)``: this is the classical Case 4
obstruction.

The implementation follows the classical algorithm and was written with
reference to the exposition and public Maple package of N. M. Abbasi
(arXiv:2211.00804).  It is an independent SymPy implementation, not a
line-by-line translation.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from enum import Enum

import sympy as sp

from .core import x as default_x


class KovacicStatus(str, Enum):
    """Outcome of Kovacic's algorithm."""

    SOLVED = "solved"
    OBSTRUCTION = "obstruction"
    SOLVER_LIMITED = "solver_limited"


@dataclass
class KovacicResult:
    """Result of :meth:`KovacicSolver.solve`.

    ``status == SOLVER_LIMITED`` means that a case was applicable but an
    internal symbolic subproblem (polynomial system, root extraction or
    quadrature) could not be completed.  It is deliberately kept distinct
    from ``OBSTRUCTION``, which is a mathematical conclusion.
    """

    status: KovacicStatus
    delta: sp.Expr | None = None
    case: int | None = None
    message: str = ""


class _NotApplicable(Exception):
    """The necessary conditions of a Kovacic case fail."""


class _SolverLimited(Exception):
    """The case applies, but an internal symbolic subproblem failed."""


def _solve_identity(expr: sp.Expr, unknowns: tuple[sp.Symbol, ...],
                    x: sp.Symbol) -> dict[sp.Symbol, sp.Expr] | None:
    """Solve a polynomial identity in ``x`` for ``unknowns``.

    Returns ``None`` when the system has no solution, and raises
    :class:`_SolverLimited` when the system cannot be formed or solved.
    """
    num, _ = sp.fraction(sp.cancel(expr))
    try:
        poly = sp.Poly(sp.expand(num), x)
    except sp.PolynomialError as exc:
        raise _SolverLimited("could not form a polynomial identity") from exc
    equations = poly.coeffs()
    try:
        solutions = sp.solve(equations, unknowns, dict=True)
    except Exception as exc:
        raise _SolverLimited("could not solve the polynomial system") from exc
    if not solutions:
        return None
    return solutions[0]


def _exp_of_integral(root: sp.Expr, x: sp.Symbol) -> sp.Expr:
    """Return ``exp(int root dx)`` or report that the quadrature failed."""
    antiderivative = sp.integrate(root, x)
    if antiderivative.has(sp.Integral):
        raise _SolverLimited("could not integrate the Riccati solution")
    return sp.simplify(sp.exp(antiderivative))


class KovacicSolver:
    """Kovacic's algorithm for ``y'' + kappa(x) y = 0``.

    Parameters
    ----------
    x:
        Independent variable.  Defaults to the real symbol ``x`` used by the
        rest of the package.
    """

    def __init__(self, x: sp.Symbol | None = None):
        self.x = default_x if x is None else x

    # ------------------------------------------------------------------
    # Rational-function utilities
    # ------------------------------------------------------------------
    def cancel(self, expr: sp.Expr) -> sp.Expr:
        return sp.cancel(sp.together(expr))

    def num_den(self, r: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
        return sp.fraction(self.cancel(r))

    def poles_and_orders(self, r: sp.Expr) -> list[tuple[sp.Expr, int]]:
        x = self.x
        _, den = self.num_den(r)
        if den == 1:
            return []
        factors = sp.factor_list(den, x)[1]
        poles: list[tuple[sp.Expr, int]] = []
        for factor, multiplicity in factors:
            roots = sp.solve(factor, x)
            for root in roots:
                poles.append((sp.simplify(root), multiplicity))
        return poles

    def order_at_infinity(self, r: sp.Expr) -> sp.Expr:
        x = self.x
        num, den = self.num_den(r)
        if num == 0:
            return -sp.oo
        return sp.degree(den, x) - sp.degree(num, x)

    def laurent_coeff(self, f: sp.Expr, c: sp.Expr, m: int, n: int) -> sp.Expr:
        """Coefficient of ``(x-c)^(-n)`` in the Laurent expansion at ``c``.

        ``m`` is the pole order of ``f`` at ``c`` (or its order at infinity).
        """
        x = self.x
        if n > m:
            return sp.Integer(0)
        if c is sp.oo or c == sp.oo:
            y = sp.Dummy("y", positive=True)
            fy = f.subs(x, 1 / y)
            expr = y**m * fy
            value = sp.limit(sp.diff(expr, y, m - n), y, 0) / sp.factorial(m - n)
        else:
            expr = (x - c) ** m * f
            value = sp.limit(sp.diff(expr, x, m - n), x, c) / sp.factorial(m - n)
        return sp.simplify(value)

    def coeff_at_infinity(self, f: sp.Expr, n: int) -> sp.Expr:
        """Coefficient of ``x^n`` in the polynomial part of ``f`` at infinity."""
        x = self.x
        f = self.cancel(f)
        num, den = sp.fraction(f)
        quotient, _ = sp.div(sp.Poly(num, x), sp.Poly(den, x))
        return sp.expand(quotient.as_expr()).coeff(x, n)

    def is_zero(self, expr: sp.Expr) -> bool:
        expr = sp.simplify(expr)
        if expr == 0:
            return True
        try:
            return sp.cancel(sp.together(expr)) == 0
        except Exception:
            return False

    # ------------------------------------------------------------------
    # Case 1
    # ------------------------------------------------------------------
    def case1(self, r: sp.Expr) -> sp.Expr | None:
        x = self.x
        poles = self.poles_and_orders(r)
        order_inf = self.order_at_infinity(r)
        orders = [multiplicity for _, multiplicity in poles]
        high_orders = [multiplicity for multiplicity in orders if multiplicity > 2]
        if any(multiplicity % 2 == 1 for multiplicity in high_orders):
            raise _NotApplicable("odd pole of order > 2")
        if not (
            (order_inf < 0 and order_inf % 2 == 0)
            or order_inf == 0
            or order_inf > 1
        ):
            raise _NotApplicable("order at infinity not allowed in Case 1")

        gamma: list[tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]] = []
        for c, multiplicity in poles:
            if multiplicity == 1:
                sqrt_r, alpha_plus, alpha_minus = (
                    sp.Integer(0),
                    sp.Integer(1),
                    sp.Integer(1),
                )
            elif multiplicity == 2:
                b = self.laurent_coeff(r, c, 2, 2)
                disc = sp.sqrt(1 + 4 * b)
                alpha_plus = sp.Rational(1, 2) + sp.Rational(1, 2) * disc
                alpha_minus = sp.Rational(1, 2) - sp.Rational(1, 2) * disc
                sqrt_r = sp.Integer(0)
            else:
                half = multiplicity // 2
                sqrt_r = sp.Integer(0)
                a = None
                for i in range(2, half + 1):
                    ai = self.laurent_coeff(sp.sqrt(r), c, half, i)
                    sqrt_r += ai / (x - c) ** i
                    if i == half:
                        a = ai
                b = self.laurent_coeff(r, c, multiplicity, half + 1) - (
                    self.laurent_coeff(sp.sqrt(r), c, half, half + 1)
                )
                alpha_plus = sp.Rational(1, 2) * (b / a + half)
                alpha_minus = sp.Rational(1, 2) * (-b / a + half)
            gamma.append(
                (
                    c,
                    sp.simplify(sqrt_r),
                    sp.simplify(alpha_plus),
                    sp.simplify(alpha_minus),
                )
            )

        if order_inf > 2:
            sqrt_r_inf, alpha_plus_inf, alpha_minus_inf = (
                sp.Integer(0),
                sp.Integer(0),
                sp.Integer(1),
            )
        elif order_inf == 2:
            num, den = self.num_den(r)
            b = sp.LC(num, x) / sp.LC(den, x) if num != 0 else 0
            disc = sp.sqrt(1 + 4 * b)
            alpha_plus_inf = sp.Rational(1, 2) + sp.Rational(1, 2) * disc
            alpha_minus_inf = sp.Rational(1, 2) - sp.Rational(1, 2) * disc
            sqrt_r_inf = sp.Integer(0)
        else:
            half = (-order_inf) // 2
            sqrt_r_inf = sp.Integer(0)
            a = None
            for i in range(0, half + 1):
                ai = self.laurent_coeff(sp.sqrt(r), sp.oo, half, i)
                sqrt_r_inf += ai * x**i
                if i == half:
                    a = ai
            b_r = self.coeff_at_infinity(r, half - 1)
            b_sq = self.coeff_at_infinity(sqrt_r_inf**2, half - 1)
            b = b_r - b_sq
            alpha_plus_inf = sp.Rational(1, 2) * (b / a - half)
            alpha_minus_inf = sp.Rational(1, 2) * (-b / a - half)

        sqrt_r_inf = sp.simplify(sqrt_r_inf)
        alpha_plus_inf = sp.simplify(alpha_plus_inf)
        alpha_minus_inf = sp.simplify(alpha_minus_inf)

        n = len(gamma)
        for signs in itertools.product([1, -1], repeat=n):
            alphas = [
                g[2] if sign == 1 else g[3]
                for g, sign in zip(gamma, signs)
            ]
            for inf_sign, alpha_inf in [
                (1, alpha_plus_inf),
                (-1, alpha_minus_inf),
            ]:
                d = sp.simplify(
                    sum(
                        (alpha_plus if sign == 1 else -alpha_minus)
                        for (_, _, alpha_plus, alpha_minus), sign in zip(
                            gamma, signs
                        )
                    )
                    + alpha_inf
                )
                if not (d.is_integer and d >= 0):
                    continue
                w = sp.Integer(0)
                for (c, sqrt_r, alpha_plus, alpha_minus), sign, alpha in zip(
                    gamma, signs, alphas
                ):
                    w += sign * sqrt_r + alpha / (x - c)
                w += sqrt_r_inf if inf_sign == 1 else -sqrt_r_inf
                w = sp.simplify(w)
                solution = self.case1_step3(r, int(d), w)
                if solution is not None:
                    return solution
        return None

    def case1_step3(self, r: sp.Expr, d: int, w: sp.Expr) -> sp.Expr | None:
        x = self.x
        if d == 0:
            p = sp.Integer(1)
            expr = sp.simplify(
                sp.diff(p, x, 2)
                + 2 * w * sp.diff(p, x)
                + (sp.diff(w, x) + w**2 - r) * p
            )
            if not self.is_zero(expr):
                return None
        else:
            a = sp.symbols(f"a0:{d}")
            p = x**d + sum(a[i] * x**i for i in range(d))
            expr = sp.simplify(
                sp.diff(p, x, 2)
                + 2 * w * sp.diff(p, x)
                + (sp.diff(w, x) + w**2 - r) * p
            )
            num, _ = sp.fraction(self.cancel(expr))
            try:
                poly = sp.Poly(sp.expand(num), x)
            except sp.PolynomialError as exc:
                raise _SolverLimited("could not form a polynomial identity") from exc
            equations = poly.coeffs()
            try:
                solutions = sp.solve(equations, a, dict=True)
            except Exception as exc:
                raise _SolverLimited("could not solve the polynomial system") from exc
            if not solutions:
                return None
            solution = solutions[0]
            p = p.subs(solution)
            for ai in a:
                if ai not in solution:
                    p = p.subs(ai, 1)
            p = sp.expand(p)

        W = sp.simplify(sp.diff(p, x) / p + w)
        if not self.is_zero(sp.diff(W, x) + W**2 - r):
            return None
        return _exp_of_integral(w, x)

    # ------------------------------------------------------------------
    # Case 2
    # ------------------------------------------------------------------
    def case2(self, r: sp.Expr) -> sp.Expr | None:
        x = self.x
        poles = self.poles_and_orders(r)
        order_inf = self.order_at_infinity(r)
        orders = [multiplicity for _, multiplicity in poles]
        high_orders = [multiplicity for multiplicity in orders if multiplicity > 2]
        if not poles or not (
            any(multiplicity % 2 == 1 for multiplicity in high_orders)
            or 2 in orders
        ):
            raise _NotApplicable("Case 2 conditions fail")

        gamma: list[tuple[sp.Expr, set[sp.Expr]]] = []
        for c, multiplicity in poles:
            if multiplicity == 1:
                exponents = {sp.Integer(4)}
            elif multiplicity == 2:
                b = self.laurent_coeff(r, c, 2, 2)
                disc = sp.sqrt(1 + 4 * b)
                values = [2, 2 + 2 * disc, 2 - 2 * disc]
                exponents = {
                    sp.simplify(value)
                    for value in values
                    if sp.simplify(value).is_integer
                }
            else:
                exponents = {sp.Integer(multiplicity)}
            gamma.append((c, exponents))

        if order_inf > 2:
            infinity_exponents = {sp.Integer(0), sp.Integer(2), sp.Integer(4)}
        elif order_inf == 2:
            num, den = self.num_den(r)
            b = sp.LC(num, x) / sp.LC(den, x) if num != 0 else 0
            disc = sp.sqrt(1 + 4 * b)
            values = [2, 2 + 2 * disc, 2 - 2 * disc]
            infinity_exponents = {
                sp.simplify(value)
                for value in values
                if sp.simplify(value).is_integer
            }
        else:
            infinity_exponents = {sp.Integer(order_inf)}

        lists = [sorted(exponents, key=str) for _, exponents in gamma]
        combinations = itertools.product(*lists) if lists else [()]
        for combination in combinations:
            for infinity_exponent in infinity_exponents:
                d = sp.Rational(1, 2) * (
                    infinity_exponent - sum(combination)
                )
                if not (d.is_integer and d >= 0):
                    continue
                theta = sp.Rational(1, 2) * sum(
                    exponent / (x - c)
                    for (c, _), exponent in zip(gamma, combination)
                )
                theta = sp.simplify(theta)
                solution = self.case2_step3(r, int(d), theta)
                if solution is not None:
                    return solution
        return None

    def case2_step3(self, r: sp.Expr, d: int, theta: sp.Expr) -> sp.Expr | None:
        x = self.x
        if d == 0:
            p = sp.Integer(1)
            expr = sp.simplify(
                sp.diff(theta, x, 2)
                + 3 * theta * sp.diff(theta, x)
                + theta**3
                - 4 * r * theta
                - 2 * sp.diff(r, x)
            )
            if not self.is_zero(expr):
                return None
        else:
            a = sp.symbols(f"b0:{d}")
            p = x**d + sum(a[i] * x**i for i in range(d))
            expr = sp.simplify(
                sp.diff(p, x, 3)
                + 3 * theta * sp.diff(p, x, 2)
                + (3 * theta**2 + 3 * sp.diff(theta, x) - 4 * r)
                * sp.diff(p, x)
                + (
                    sp.diff(theta, x, 2)
                    + 3 * theta * sp.diff(theta, x)
                    + theta**3
                    - 4 * r * theta
                    - 2 * sp.diff(r, x)
                )
                * p
            )
            solution = _solve_identity(expr, a, x)
            if solution is None:
                return None
            p = p.subs(solution)
            for ai in a:
                if ai not in solution:
                    p = p.subs(ai, 1)
            p = sp.expand(p)

        phi = sp.simplify(theta + sp.diff(p, x) / p)
        w = sp.Symbol("w")
        equation = (
            w**2
            - phi * w
            + sp.Rational(1, 2) * sp.diff(phi, x)
            + sp.Rational(1, 2) * phi**2
            - r
        )
        try:
            roots = sp.solve(equation, w)
        except Exception as exc:
            raise _SolverLimited("could not solve the Riccati equation") from exc
        if not roots:
            return None
        for root in roots:
            root = sp.simplify(root)
            if self.is_zero(sp.diff(root, x) + root**2 - r):
                return _exp_of_integral(root, x)
        return None

    # ------------------------------------------------------------------
    # Case 3
    # ------------------------------------------------------------------
    def case3(self, r: sp.Expr) -> sp.Expr | None:
        x = self.x
        poles = self.poles_and_orders(r)
        order_inf = self.order_at_infinity(r)
        if any(multiplicity > 2 for _, multiplicity in poles) or order_inf < 2:
            raise _NotApplicable("Case 3 conditions fail")

        for n in (4, 6, 12):
            gamma: list[tuple[sp.Expr, set[sp.Expr]]] = []
            for c, multiplicity in poles:
                if multiplicity == 1:
                    exponents = {sp.Integer(12)}
                elif multiplicity == 2:
                    b = self.laurent_coeff(r, c, 2, 2)
                    disc = sp.sqrt(1 + 4 * b)
                    values = [
                        6 + sp.Rational(12 * k, n) * disc
                        for k in range(-n // 2, n // 2 + 1)
                    ]
                    exponents = {
                        sp.simplify(value)
                        for value in values
                        if sp.simplify(value).is_integer
                    }
                else:
                    exponents = set()
                gamma.append((c, exponents))

            num, den = self.num_den(r)
            b = sp.LC(num, x) / sp.LC(den, x) if num != 0 else 0
            disc = sp.sqrt(1 + 4 * b)
            values = [
                6 + sp.Rational(12 * k, n) * disc
                for k in range(-n // 2, n // 2 + 1)
            ]
            infinity_exponents = {
                sp.simplify(value)
                for value in values
                if sp.simplify(value).is_integer
            }

            lists = [sorted(exponents, key=str) for _, exponents in gamma]
            combinations = itertools.product(*lists) if lists else [()]
            for combination in combinations:
                for infinity_exponent in infinity_exponents:
                    d = sp.Rational(n, 12) * (
                        infinity_exponent - sum(combination)
                    )
                    if not (d.is_integer and d >= 0):
                        continue
                    theta = sp.Rational(n, 12) * sum(
                        exponent / (x - c)
                        for (c, _), exponent in zip(gamma, combination)
                    )
                    theta = sp.simplify(theta)
                    S = sp.prod([x - c for c, _ in poles])
                    solution = self.case3_step3(
                        r, n, int(d), theta, sp.expand(S)
                    )
                    if solution is not None:
                        return solution
        return None

    def case3_step3(self, r: sp.Expr, n: int, d: int, theta: sp.Expr,
                    S: sp.Expr) -> sp.Expr | None:
        x = self.x
        if d == 0:
            p = sp.Integer(1)
        else:
            a = sp.symbols(f"c0:{d}")
            p = x**d + sum(a[i] * x**i for i in range(d))

        P: dict[int, sp.Expr] = {n: -p}
        for i in range(n, -1, -1):
            P_next = P.get(i + 1, sp.Integer(0))
            if i == n:
                P[i - 1] = -S * sp.diff(P[i], x) - S * theta * P[i]
            else:
                P[i - 1] = (
                    -S * sp.diff(P[i], x)
                    + ((n - i) * sp.diff(S, x) - S * theta) * P[i]
                    - (n - i) * (i + 1) * S**2 * r * P_next
                )
            P[i - 1] = sp.simplify(P[i - 1])

        if d > 0:
            solution = _solve_identity(P[-1], a, x)
            if solution is None:
                return None
            p = p.subs(solution)
            for ai in a:
                if ai not in solution:
                    p = p.subs(ai, 1)
            p = sp.expand(p)
            P = {n: -p}
            for i in range(n, -1, -1):
                P_next = P.get(i + 1, sp.Integer(0))
                if i == n:
                    P[i - 1] = -S * sp.diff(P[i], x) - S * theta * P[i]
                else:
                    P[i - 1] = (
                        -S * sp.diff(P[i], x)
                        + ((n - i) * sp.diff(S, x) - S * theta) * P[i]
                        - (n - i) * (i + 1) * S**2 * r * P_next
                    )
                P[i - 1] = sp.simplify(P[i - 1])

        w = sp.Symbol("w")
        omega_equation = sum(
            S**i * P[i] / sp.factorial(n - i) * w**i
            for i in range(0, n + 1)
        )
        omega_equation = sp.simplify(omega_equation)
        try:
            roots = sp.solve(omega_equation, w)
        except Exception as exc:
            raise _SolverLimited("could not solve the algebraic equation") from exc
        if not roots:
            return None
        for root in roots:
            if root.has(sp.RootOf):
                try:
                    root = sp.allvalues(root)[0]
                except Exception as exc:
                    raise _SolverLimited("could not extract a radical root") from exc
            root = sp.simplify(root)
            if self.is_zero(sp.diff(root, x) + root**2 - r):
                return _exp_of_integral(root, x)
        return None

    # ------------------------------------------------------------------
    # Driver
    # ------------------------------------------------------------------
    def solve(self, kappa: sp.Expr) -> KovacicResult:
        """Run Cases 1, 2 and 3 on ``y'' + kappa y = 0``.

        The returned status is:

        * ``SOLVED`` with a nonzero Liouvillian solution ``delta``;
        * ``OBSTRUCTION`` when no case applies or no applicable case yields a
          solution (Kovacic's Case 4);
        * ``SOLVER_LIMITED`` when a case applies but an internal symbolic
          computation could not be completed.
        """
        kappa = sp.simplify(kappa)
        if not kappa.is_rational_function(self.x):
            raise ValueError("Kovacic's algorithm requires kappa in C(x)")
        if kappa == 0:
            return KovacicResult(
                KovacicStatus.SOLVED,
                sp.Integer(1),
                case=1,
                message="kappa = 0",
            )
        r = sp.simplify(-kappa)
        for case in (1, 2, 3):
            try:
                solution = getattr(self, f"case{case}")(r)
            except _NotApplicable:
                continue
            except _SolverLimited as exc:
                return KovacicResult(
                    KovacicStatus.SOLVER_LIMITED,
                    case=case,
                    message=str(exc),
                )
            except Exception as exc:  # unexpected internal failure
                return KovacicResult(
                    KovacicStatus.SOLVER_LIMITED,
                    case=case,
                    message=f"{type(exc).__name__}: {exc}",
                )
            if solution is not None:
                return KovacicResult(
                    KovacicStatus.SOLVED,
                    sp.simplify(solution),
                    case=case,
                )
        return KovacicResult(
            KovacicStatus.OBSTRUCTION,
            message="no solution in Cases 1-3",
        )
