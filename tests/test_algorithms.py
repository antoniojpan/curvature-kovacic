import pytest
import sympy as sp

from curvature_kovacic import (
    KovacicSolver,
    kovacic_integration,
    solve_from_linear_solution,
    solve_from_particular,
    u,
    x,
)

from conftest import (
    assert_zero,
    check_particular_result,
    check_standard_result,
)

X_VALUES = [sp.Rational(3, 2), 2]

PHI_EULER = u / (2 * x) + sp.Rational(3, 2) / x * sp.sqrt(u**2 - 4 * x)
PHI_CIRCLE = sp.sqrt(1 - u**2) / 2
PHI_CIRCLE_FORCING = 1 + sp.sqrt(1 - (u - x) ** 2) / 2
PHI_HYPERBOLA = (
    u / (2 * x)
    + sp.Rational(5, 2) * x**2
    + sp.Rational(3, 2) / x * sp.sqrt((u - x**3) ** 2 - 4 * x)
)
PHI_ZERO = (1 + 2 * u * x - sp.sqrt(1 + 4 * u * x)) / (2 * x**2)
PHI_GALOIS = u / (4 * x) + sp.sqrt(u**2 - 4 * sp.sqrt(x)) / sp.sqrt(x)
KAPPA_GALOIS = sp.Rational(3, 16) / x**2 - 1 / x
DELTA_AIRY = sp.airyai(-x)
PHI_AIRY = sp.diff(DELTA_AIRY, x) / DELTA_AIRY * u
GALOIS_POINTS = [
    {x: 1, u: 3},
    {x: sp.Rational(3, 2), u: 4},
    {x: 2, u: 5},
]

A1_CASES = [
    dict(
        name="affine Euler-Cauchy (c=0)",
        phi=2 * u / x,
        delta=x**2,
        kappa=-2 / x**2,
        c=0,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
    dict(
        name="affine Euler-Cauchy (c=x^3)",
        phi=2 * u / x + x**3,
        delta=x**2,
        kappa=-2 / x**2,
        c=5 * x**2,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
    dict(
        name="affine zero curvature",
        phi=u / x + x,
        delta=x,
        kappa=0,
        c=2,
        x0=1,
        v0=1,
        xs=X_VALUES,
    ),
    dict(
        name="affine kappa=-2/(x^2+1)",
        phi=2 * x * u / (x**2 + 1),
        delta=x**2 + 1,
        kappa=-2 / (x**2 + 1),
        c=0,
        x0=0,
        v0=1,
        xs=[sp.Rational(1, 2), 1],
    ),
    dict(
        name="nonlinear Euler-Cauchy",
        phi=PHI_EULER,
        delta=x**2,
        kappa=-2 / x**2,
        c=0,
        x0=2,
        v0=5,
        xs=[sp.Rational(5, 2), 3],
    ),
    dict(
        name="nonlinear circle",
        phi=PHI_CIRCLE,
        delta=sp.cos(x / 2),
        kappa=sp.Rational(1, 4),
        c=0,
        x0=0,
        v0=0,
        xs=[sp.Rational(1, 2), 1],
    ),
    dict(
        name="nonlinear circle with forcing",
        phi=PHI_CIRCLE_FORCING,
        delta=sp.cos(x / 2),
        kappa=sp.Rational(1, 4),
        c=x / 4,
        x0=0,
        v0=0,
        xs=[sp.Rational(1, 2), 1],
    ),
    dict(
        name="nonlinear hyperbola with forcing",
        phi=PHI_HYPERBOLA,
        delta=x**2,
        kappa=-2 / x**2,
        c=4 * x,
        x0=2,
        v0=12,
        xs=[sp.Rational(5, 2), 3],
    ),
    dict(
        name="nonlinear kappa=0",
        phi=PHI_ZERO,
        delta=sp.Integer(1),
        kappa=0,
        c=0,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
    dict(
        name="Kovacic Example 1 (paper delta)",
        phi=PHI_GALOIS,
        delta=x ** sp.Rational(1, 4) * sp.exp(2 * sp.sqrt(x)),
        kappa=KAPPA_GALOIS,
        c=0,
        x0=1,
        v0=3,
        xs=X_VALUES,
        points=GALOIS_POINTS,
    ),
    dict(
        name="Airy linear",
        phi=PHI_AIRY,
        delta=DELTA_AIRY,
        kappa=x,
        c=0,
        x0=0,
        v0=1,
        xs=[sp.Rational(1, 2), 1],
        points=[{x: 1, u: 1}],
    ),
]


@pytest.mark.parametrize("case", A1_CASES, ids=lambda case: case["name"])
def test_algorithm_1(case):
    result = solve_from_linear_solution(
        case["phi"],
        case["delta"],
        x0=case["x0"],
        v_value=case["v0"],
        points=case.get("points"),
    )
    check_standard_result(
        result,
        case["phi"],
        case["kappa"],
        case["c"],
        case["x0"],
        case["v0"],
        case["xs"],
    )


A2_CASES = [
    dict(
        name="affine Euler-Cauchy (c=0)",
        phi=2 * u / x,
        f=x**2,
        kappa=-2 / x**2,
        c=0,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
    dict(
        name="affine Euler-Cauchy (c=x^3)",
        phi=2 * u / x + x**3,
        f=x**4 / 2,
        kappa=-2 / x**2,
        c=5 * x**2,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
    dict(
        name="affine kappa=-2/(x^2+1)",
        phi=2 * x * u / (x**2 + 1),
        f=x**2 + 1,
        kappa=-2 / (x**2 + 1),
        c=0,
        x0=0,
        v0=1,
        xs=[sp.Rational(1, 2), 1],
    ),
    dict(
        name="affine trigonometric",
        phi=-sp.tan(x) * u,
        f=sp.cos(x),
        kappa=1,
        c=0,
        x0=0,
        v0=1,
        xs=[sp.Rational(1, 2), 1],
    ),
    dict(
        name="nonlinear kappa=0",
        phi=PHI_ZERO,
        f=sp.Integer(0),
        kappa=0,
        c=0,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
]


@pytest.mark.parametrize("case", A2_CASES, ids=lambda case: case["name"])
def test_algorithm_2(case):
    result = solve_from_particular(
        case["phi"],
        case["f"],
        x0=case["x0"],
        v_value=case["v0"],
        points=case.get("points"),
    )
    check_particular_result(
        result,
        case["phi"],
        case["f"],
        case["kappa"],
        case["c"],
        case["x0"],
        case["v0"],
        case["xs"],
    )


def _case2_phi():
    kappa = (3 + 8 * x) / (16 * x**2)
    result = KovacicSolver(x).solve(kappa)
    assert result.delta is not None
    return sp.simplify(sp.diff(result.delta, x) / result.delta * u)


A3_CASES = [
    *[
        dict(
            name=f"linear Euler-Cauchy n={n}",
            phi=n * u / x,
            kappa=-n * (n - 1) / x**2,
            c=0,
            x0=1,
            v0=2,
            xs=X_VALUES,
        )
        for n in (2, 3, 4, 5)
    ],
    dict(
        name="affine Euler-Cauchy with forcing",
        phi=2 * u / x + x**3,
        kappa=-2 / x**2,
        c=5 * x**2,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
    dict(
        name="nonlinear Euler-Cauchy",
        phi=PHI_EULER,
        kappa=-2 / x**2,
        c=0,
        x0=2,
        v0=5,
        xs=[sp.Rational(5, 2), 3],
    ),
    dict(
        name="nonlinear hyperbola with forcing",
        phi=PHI_HYPERBOLA,
        kappa=-2 / x**2,
        c=4 * x,
        x0=2,
        v0=12,
        xs=[sp.Rational(5, 2), 3],
    ),
    dict(
        name="nonlinear circle",
        phi=PHI_CIRCLE,
        kappa=sp.Rational(1, 4),
        c=0,
        x0=0,
        v0=0,
        xs=[sp.Rational(1, 2), 1],
    ),
    dict(
        name="Kovacic Example 1 (nonlinear presentation)",
        phi=PHI_GALOIS,
        kappa=KAPPA_GALOIS,
        c=0,
        x0=1,
        v0=3,
        xs=X_VALUES,
        points=GALOIS_POINTS,
    ),
    dict(
        name="Kovacic Case 2 example 2",
        phi=_case2_phi(),
        kappa=(3 + 8 * x) / (16 * x**2),
        c=0,
        x0=1,
        v0=2,
        xs=X_VALUES,
    ),
]


@pytest.mark.parametrize("case", A3_CASES, ids=lambda case: case["name"])
def test_algorithm_3(case):
    result = kovacic_integration(
        case["phi"],
        x0=case["x0"],
        v_value=case["v0"],
        points=case.get("points"),
    )
    check_standard_result(
        result,
        case["phi"],
        case["kappa"],
        case["c"],
        case["x0"],
        case["v0"],
        case["xs"],
    )


def test_algorithm_3_airy_obstruction():
    result = kovacic_integration(
        PHI_AIRY, x0=0, v_value=1, points=[{x: 1, u: 1}]
    )
    assert result.status == "obstruction", result.message
    assert_zero(result.kappa - x, [{x: 1}])
    assert_zero(result.c, [{x: 1}])


def test_algorithm_3_solver_limited_is_not_obstruction(monkeypatch):
    from curvature_kovacic.kovacic import _SolverLimited

    def boom(_r):
        raise _SolverLimited("simulated failure")

    solver = KovacicSolver(x)
    monkeypatch.setattr(solver, "case1", boom)
    result = kovacic_integration(
        PHI_AIRY, x0=0, v_value=1, kovacic=solver, points=[{x: 1, u: 1}]
    )
    assert result.status == "solver_limited", result.message


def test_algorithm_1_rejects_invalid_delta():
    result = solve_from_linear_solution(
        2 * u / x, x**3, x0=1, v_value=2
    )
    assert result.status == "invalid_auxiliary_solution"


def test_algorithm_2_rejects_invalid_particular_solution():
    result = solve_from_particular(
        2 * u / x, x**3, x0=1, v_value=2
    )
    assert result.status == "invalid_particular_solution"


def test_algorithm_outside_class():
    result = solve_from_linear_solution(
        u**2, x**2, x0=1, v_value=1
    )
    assert result.status == "outside_class"
