"""Result containers and status codes for the integration algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import sympy as sp


class AlgorithmStatus(str, Enum):
    """Status of an integration algorithm.

    The values are strings so that results remain easy to print and to
    compare with plain strings.  ``OK`` means that a first integral and a
    local family were produced.
    """

    OK = "ok"
    OUTSIDE_CLASS = "outside_class"
    INVALID_AUXILIARY_SOLUTION = "invalid_auxiliary_solution"
    INVALID_PARTICULAR_SOLUTION = "invalid_particular_solution"
    INVALID_DELTA_F = "invalid_delta_f"
    NOT_RATIONAL = "not_rational"
    OBSTRUCTION = "obstruction"
    SOLVER_LIMITED = "solver_limited"
    SOLVER_ERROR = "solver_error"


@dataclass
class AlgorithmResult:
    """Output of one of the three integration algorithms.

    Attributes
    ----------
    status:
        Outcome, see :class:`AlgorithmStatus`.
    kappa, c:
        Curvature and the coefficient in ``L(u) = c`` when they were computed.
    delta:
        Auxiliary solution of ``L(delta) = 0`` used by the construction.
    S, S_raw:
        Integrating factor.  ``S`` is the multiplier used for the implicit
        first integral; ``S_raw`` is the multiplier obtained directly from
        ``delta`` before the possible switch to the independent solution.
    Q, I, I_alt:
        Quadrature and first integrals; ``I_alt`` is the first integral built
        from ``delta`` when ``I`` uses the independent solution instead.
    family:
        Local solution family ``u(x;v)``.
    delta2:
        The independent solution ``delta * int(dx/delta**2)``, when computed.
    branch:
        ``"general"``, ``"degenerate"`` or ``"particular"``.
    message:
        Human-readable diagnostic.
    """

    status: AlgorithmStatus
    kappa: sp.Expr | None = None
    c: sp.Expr | None = None
    delta: sp.Expr | None = None
    S: sp.Expr | None = None
    S_raw: sp.Expr | None = None
    Q: sp.Expr | None = None
    I: sp.Expr | None = None
    I_alt: sp.Expr | None = None
    family: sp.Expr | None = None
    delta2: sp.Expr | None = None
    branch: str = ""
    message: str = ""
