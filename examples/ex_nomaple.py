#!/usr/bin/env python3
"""Paper example ``ex:nomaple``: the ODE on which Maple's dsolve fails.

The right-hand side is used exactly as printed, with its nested radicals:

    R = ((x u + 3x^3 - 2x^9
          + sqrt(x^2 u^2 + 2x^4 (3 - 2x^6) u + 4 - 3x^6))/2)^(1/3)
    S = R + (x^6 - 1)/R - x^3
    u' = 9 x S^2 - u/x,         D = {0 < x < 1, u > 0}.

The algebraic element is deliberately *not* replaced by an abstract
symbol: that rewrites the example into a different and much easier problem
(on the RootOf presentation Maple's dsolve returns the general solution in
0.1 s, whereas on this input it does not return at all), which would make
the comparison recorded in the paper meaningless.

The curvature is cheap anyway, because ``A(phi) = phi_x + phi phi_u``
simplifies to ``2u/x^2`` *before* the derivative in ``u`` is taken:
``kappa = -2/x^2`` and ``c = 0``.  With ``delta = 1/x`` the construction of
``alg:delta`` returns expressions in the nested radicals (too long to be
informative), equal to the closed forms printed in the paper:

    I = -9 S^2,     f(x;s) = 3 s^2 x^2 + (s^3 + 3 s)/x,   s > 0.

Those two identities are what this script checks.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from curvature_kovacic import kovacic_integration, u, v, x


def main() -> None:
    R = (
        (
            x * u + 3 * x**3 - 2 * x**9
            + sp.sqrt(
                x**2 * u**2
                + 2 * x**4 * (3 - 2 * x**6) * u
                + 4 - 3 * x**6
            )
        ) / 2
    ) ** sp.Rational(1, 3)
    S = R + (x**6 - 1) / R - x**3
    phi = 9 * x * S**2 - u / x
    x0 = sp.Rational(1, 2)

    result = kovacic_integration(phi, x0=x0, v_value=1)

    print(f"status : {result.status.value} ({result.message})")
    print(f"kappa  : {result.kappa}")
    print(f"c      : {result.c}")
    print(f"delta  : {result.delta}")
    print(
        "raw output: "
        f"S = {sp.count_ops(result.S)} ops, "
        f"I = {sp.count_ops(result.I)} ops, "
        f"family = {sp.count_ops(result.family)} ops (nested radicals)"
    )

    assert result.status == "ok", result.message
    assert sp.simplify(result.kappa + 2 / x**2) == 0
    assert sp.simplify(result.c) == 0
    assert sp.simplify(result.delta - 1 / x) == 0

    def close(expr, subs, tol=1e-15) -> bool:
        return abs(complex(sp.N(expr.subs(subs), 30))) < tol

    print("\nclosed forms printed in the paper:")

    # (a) the first integral is -9 S^2
    ok_I = all(
        close(result.I + 9 * S**2, {x: xv, u: uv})
        for xv, uv in [(sp.Rational(1, 4), sp.Rational(3, 2)),
                       (sp.Rational(3, 4), sp.Rational(5, 4))]
    )
    print(f"  I = -9 S^2                              : {ok_I}")

    # (b) the family 3 s^2 x^2 + (s^3+3s)/x solves the ODE
    s = sp.Symbol("s", positive=True)
    f_paper = 3 * s**2 * x**2 + (s**3 + 3 * s) / x
    pairs = [(sp.Rational(1, 4), sp.Integer(1)),
             (sp.Rational(3, 4), sp.Integer(2))]
    ok_ode = all(
        close(
            sp.diff(f_paper, x) - phi,
            {x: xv, u: f_paper.subs({x: xv, s: sv}), s: sv},
        )
        for xv, sv in pairs
    )
    print(f"  3 s^2 x^2 + (s^3+3s)/x solves the ODE    : {ok_ode}")

    # (c) it is the family returned by the algorithm, with v = u(x0)
    fam_s = result.family.subs(v, f_paper.subs(x, x0))
    ok_fam = all(
        close(fam_s - f_paper, {x: xv, s: sv}) for xv, sv in pairs
    )
    print(f"  Algorithm 1 family = paper family        : {ok_fam}")

    assert ok_I and ok_ode and ok_fam
    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
