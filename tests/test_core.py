import sympy as sp

from curvature_kovacic import (
    c_coefficient,
    curvature,
    is_zero,
    points_from_family,
    rational_curvature,
    u,
    v,
    x,
)

from conftest import assert_zero


def galois_phi():
    return u / (4 * x) + sp.sqrt(u**2 - 4 * sp.sqrt(x)) / sp.sqrt(x)


def test_curvature_euler_cauchy():
    phi = u / (2 * x) + sp.Rational(3, 2) / x * sp.sqrt(u**2 - 4 * x)
    assert_zero(curvature(phi) + 2 / x**2,
                [{x: 2, u: 5}, {x: 3, u: 8}])


def test_curvature_galois_example():
    phi = galois_phi()
    kappa = sp.Rational(3, 16) / x**2 - 1 / x
    K = curvature(phi)
    points = [{x: 1, u: 3}, {x: 2, u: 5}]
    assert is_zero(sp.diff(K, u), points)
    assert_zero(rational_curvature(K, points) - kappa, points)


def test_c_coefficient_euler_cauchy():
    phi = 2 * u / x + x**3
    kappa = -2 / x**2
    assert_zero(c_coefficient(phi, kappa) - 5 * x**2,
                [{x: 1, u: 2}, {x: 2, u: 3}])


def test_is_zero_true():
    assert is_zero(sp.sin(x) ** 2 + sp.cos(x) ** 2 - 1)


def test_is_zero_false():
    assert not is_zero(x - 1)


def test_points_from_family():
    family = x**2 + v
    pts = points_from_family(family, [sp.Integer(2)], sp.Integer(3))
    assert pts == [{x: sp.Integer(2), u: sp.Integer(7)}]
