import pytest
import sympy as sp

from curvature_kovacic import KovacicSolver, KovacicStatus, is_zero, x
from curvature_kovacic.kovacic import _SolverLimited

PAPER_KAPPA = sp.Rational(3, 16) / x**2 - 1 / x
CASE2_KAPPA = (3 + 8 * x) / (16 * x**2)
CASE3_KAPPA = (x - 4) / (4 * x * (x - 1) ** 2)


@pytest.mark.parametrize(
    "kappa",
    [
        sp.Integer(1),
        sp.Integer(-1),
        -2 / x**2,
        PAPER_KAPPA,
        CASE2_KAPPA,
        CASE3_KAPPA,
    ],
    ids=["y''+y", "y''-y", "Euler-2", "paper-Case2", "Case2-b", "Case3"],
)
def test_solved_cases(kappa):
    result = KovacicSolver(x).solve(kappa)
    assert result.status == KovacicStatus.SOLVED
    assert result.delta is not None
    assert sp.simplify(result.delta) != 0
    assert is_zero(
        sp.diff(result.delta, x, 2) + kappa * result.delta,
        [{x: 2}, {x: 3}],
    )


def test_paper_example_is_case_two():
    result = KovacicSolver(x).solve(PAPER_KAPPA)
    assert result.status == KovacicStatus.SOLVED
    assert result.case == 2


def test_zero_curvature():
    result = KovacicSolver(x).solve(sp.Integer(0))
    assert result.status == KovacicStatus.SOLVED
    assert result.case == 1
    assert result.delta == 1


@pytest.mark.parametrize(
    "kappa",
    [x, 1 - sp.Rational(3, 4) / x**2],
    ids=["Airy", "Bessel-nu1"],
)
def test_obstructions(kappa):
    result = KovacicSolver(x).solve(kappa)
    assert result.status == KovacicStatus.OBSTRUCTION
    assert result.delta is None


def test_case3_direct():
    solver = KovacicSolver(x)
    delta = solver.case3(-CASE3_KAPPA)
    assert delta is not None
    assert is_zero(
        sp.diff(delta, x, 2) + CASE3_KAPPA * delta,
        [{x: 2}, {x: 3}],
    )


def test_non_rational_curvature_raises():
    with pytest.raises(ValueError):
        KovacicSolver(x).solve(sp.sqrt(x))


def test_solver_limited_is_not_an_obstruction(monkeypatch):
    solver = KovacicSolver(x)

    def boom(_r):
        raise _SolverLimited("simulated failure")

    monkeypatch.setattr(solver, "case1", boom)
    result = solver.solve(sp.Integer(1))
    assert result.status == KovacicStatus.SOLVER_LIMITED
    assert result.case == 1
    assert "simulated failure" in result.message
