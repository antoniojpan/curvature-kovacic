#!/usr/bin/env python3
"""Apply Algorithm 3 to the ODE in paper example ``ex:nomaple``."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from curvature_kovacic import kovacic_integration, u, v, x


class CardanoS(sp.Function):
    """Symbolic representation of the algebraic branch in the ODE."""

    nargs = 2

    def fdiff(self, argindex: int = 1) -> sp.Expr:
        xx, uu = self.args
        denominator = 3 * (self**2 + 2 * xx**3 * self + 1)
        if argindex == 1:
            return (uu - 9 * xx**2 * self**2) / denominator
        if argindex == 2:
            return xx / denominator
        raise sp.ArgumentIndexError(self, argindex)


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
    cardano_expression = R + (x**6 - 1) / R - x**3
    phi = 9 * x * cardano_expression**2 - u / x

    S = CardanoS(x, u)
    phi_for_algorithm = phi.xreplace({cardano_expression: S})
    assert phi_for_algorithm == 9 * x * S**2 - u / x
    x0 = sp.Rational(1, 2)

    result = kovacic_integration(phi_for_algorithm, x0=x0, v_value=1)
    assert result.status == "ok", result.message
    assert sp.simplify(result.kappa + 2 / x**2) == 0
    assert sp.simplify(result.c) == 0
    assert sp.simplify(result.delta - 1 / x) == 0
    assert sp.simplify(result.I + 9 * S**2) == 0
    assert result.I.xreplace({S: cardano_expression}) == -9 * cardano_expression**2

    s = sp.Symbol("s", positive=True)
    v_from_s = (s**3 + 3 * x0**3 * s**2 + 3 * s) / x0
    family = result.family.subs(CardanoS(x0, v), s).subs(v, v_from_s)
    f_paper = 3 * s**2 * x**2 + (s**3 + 3 * s) / x
    assert sp.simplify(family - f_paper) == 0

    print(f"status : {result.status.value}")
    print(f"kappa  : {result.kappa}")
    print(f"c      : {result.c}")
    print(f"delta  : {result.delta}")
    R_symbol = sp.Symbol("R")
    print(f"I      : {result.I.xreplace({S: R_symbol + (x**6 - 1)/R_symbol - x**3})}")
    print(f"family : {f_paper},  s > 0")
    print("Maple 2015.0 dsolve: NULL after 4.290 CPU seconds (paper's benchmark).")


if __name__ == "__main__":
    main()
