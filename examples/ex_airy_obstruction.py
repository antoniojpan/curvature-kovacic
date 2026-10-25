#!/usr/bin/env python3
"""Airy obstruction (``ex:airy_obstruction``).

The linear equation is ``u' = (y'/y) u`` where ``y`` solves the Airy
equation ``y'' + x y = 0``.  Its curvature is ``kappa = x``.  Over
``C(x)`` the Airy equation has no nonzero Liouvillian solution, so
``alg:kovacic`` must report an obstruction.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from curvature_kovacic import (
    is_zero,
    kovacic_integration,
    solve_from_linear_solution,
    u,
    x,
)


def main() -> None:
    delta = sp.airyai(-x)
    phi = sp.diff(delta, x) / delta * u
    kappa = x

    print("=" * 78)
    print("A1: alg:delta with delta = Ai(-x)")
    print("=" * 78)
    result_a1 = solve_from_linear_solution(
        phi, delta, x0=0, v_value=1, points=[{x: 1, u: 1}]
    )
    print(f"status : {result_a1.status}")
    print(f"kappa  : {sp.simplify(result_a1.kappa)}")
    print(f"c      : {sp.simplify(result_a1.c)}")
    print(f"delta  : {result_a1.delta}")
    print(f"I      : {sp.simplify(result_a1.I)}")
    print(f"family : {sp.simplify(result_a1.family)}")

    assert result_a1.status == "ok", result_a1.message
    assert is_zero(result_a1.kappa - kappa, [{x: 1}])
    assert is_zero(result_a1.c, [{x: 1}])

    print("\n" + "=" * 78)
    print("A3: alg:kovacic on the Airy equation")
    print("=" * 78)
    result_a3 = kovacic_integration(
        phi, x0=0, v_value=1, points=[{x: 1, u: 1}]
    )
    print(f"status  : {result_a3.status}")
    print(f"kappa   : {sp.simplify(result_a3.kappa)}")
    print(f"c       : {sp.simplify(result_a3.c)}")
    print(f"message : {result_a3.message}")

    assert result_a3.status == "obstruction", result_a3.message
    assert is_zero(result_a3.kappa - kappa, [{x: 1}])
    assert is_zero(result_a3.c, [{x: 1}])

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
