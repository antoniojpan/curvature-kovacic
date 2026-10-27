"""Curvature--Kovacic algorithms for scalar first-order ODEs.

Companion code for

    A. J. Pan-Collantes, J. A. Alvarez-Garcia,
    "From curvature to Kovacic: a geometric approach to integrability of
    scalar ODEs".

The public entry points are the three algorithms

* :func:`solve_from_linear_solution` (``alg:delta``),
* :func:`solve_from_particular` (``alg:particular``),
* :func:`kovacic_integration` (``alg:kovacic``),

and the self-contained :class:`KovacicSolver`.
"""

from .algorithms import (
    full_kovacic,
    kovacic_integration,
    solve_from_linear_solution,
    solve_from_particular,
)
from .core import (
    A,
    c_coefficient,
    curvature,
    divergence,
    is_zero,
    points_from_family,
    safe_simplify,
    u,
    v,
    x,
)
from .kovacic import KovacicResult, KovacicSolver, KovacicStatus
from .results import AlgorithmResult, AlgorithmStatus
from .simplification import rational_curvature, refine_delta

__version__ = "1.0.0"

#: Paper-oriented aliases.
integrate_from_linear_solution = solve_from_linear_solution
integrate_from_particular = solve_from_particular
integrate_with_kovacic = kovacic_integration

__all__ = [
    "A",
    "AlgorithmResult",
    "AlgorithmStatus",
    "KovacicResult",
    "KovacicSolver",
    "KovacicStatus",
    "c_coefficient",
    "curvature",
    "divergence",
    "full_kovacic",
    "integrate_from_linear_solution",
    "integrate_from_particular",
    "integrate_with_kovacic",
    "is_zero",
    "kovacic_integration",
    "points_from_family",
    "rational_curvature",
    "refine_delta",
    "safe_simplify",
    "solve_from_linear_solution",
    "solve_from_particular",
    "u",
    "v",
    "x",
]
