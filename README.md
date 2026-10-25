# curvature-kovacic

Companion code for

> A. J. Pan-Collantes, J. A. Alvarez-Garcia,
> *From curvature to Kovacic: a geometric approach to integrability of
> scalar ODEs*.

The package implements the three constructive algorithms of the paper for
scalar first-order equations

$$
u' = \phi(x,u)
$$

whose Gauss curvature depends only on the independent variable:

$$
\mathcal K(x,u) = \kappa(x).
$$

For such an equation, the flow divergence along solutions satisfies a
Riccati equation which linearizes to

$$
L(y) = y'' + \kappa(x)y = 0,
$$

and every solution of the nonlinear equation satisfies
$L(u) = c(x)$.  A nonzero solution $\delta$ of $L(\delta)=0$ produces an
integrating factor and a first integral by quadratures.  When $\kappa$ is
rational, Kovacic's algorithm supplies such a $\delta$.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

SymPy is pinned to `1.14.0` for reproducibility; see `requirements.txt`.
The package requires Python 3.10 or newer.

## Quick start

```python
import sympy as sp

from curvature_kovacic import solve_from_linear_solution, u, x

phi = u / (4 * x) + sp.sqrt(u**2 - 4 * sp.sqrt(x)) / sp.sqrt(x)
delta = x ** sp.Rational(1, 4) * sp.exp(2 * sp.sqrt(x))

result = solve_from_linear_solution(phi, delta, x0=1, v_value=3)

print(result.status)    # ok
print(result.kappa)     # 3/(16*x**2) - 1/x
print(result.c)         # 0
print(result.I)         # first integral
print(result.family)    # local family u(x;v)
```

The three public entry points are:

| Paper algorithm | Function |
|---|---|
| `alg:delta` | `solve_from_linear_solution(phi, delta, x0, v_value)` |
| `alg:particular` | `solve_from_particular(phi, f, x0, v_value)` |
| `alg:kovacic` | `kovacic_integration(phi, x0, v_value)` |

The paper-oriented aliases `integrate_from_linear_solution`,
`integrate_from_particular` and `integrate_with_kovacic` are also available.

## Examples

| Script | Content |
|---|---|
| `examples/ex_galois_positive.py` | Kovacic's Example 1, Section 4.2 (`ex:galois_positive`) |
| `examples/ex_euler_cauchy.py` | Euler--Cauchy family (`ex:euler_cauchy`) |
| `examples/ex_particular_solution.py` | `alg:particular` on an affine Euler--Cauchy equation |
| `examples/ex_circle.py` | Circle equation `u' = sqrt(1-u^2)/2` (`ex:example_sqrt`) |
| `examples/ex_airy_obstruction.py` | Airy obstruction (`ex:airy_obstruction`) |

Run one with, for example,

```bash
python examples/ex_galois_positive.py
```

Each script prints the curvature, the coefficient $c$, the auxiliary
solution, the first integral and the local family, and verifies the
identities used in the paper.

## Verification

```bash
pytest -q
```

The test suite covers the algorithms, the Kovacic solver, the paper
examples and the regression examples.  The identities are checked by exact
symbolic simplification where possible and by high-precision numerical
evaluation at sample points otherwise.  The numerical fallback is a
verification aid, not a proof; see `is_zero` in
`src/curvature_kovacic/core.py`.

## Kovacic status semantics

`KovacicSolver.solve(kappa)` returns a `KovacicResult` whose status is one of

- `SOLVED`: a nonzero Liouvillian solution $\delta$ was found.
- `OBSTRUCTION`: no Case 1--3 solution exists, so $L(y)=0$ has no nonzero
  Liouvillian solution over $\mathbb C(x)$ (Kovacic's Case 4).
- `SOLVER_LIMITED`: a case was applicable but an internal symbolic
  subproblem could not be completed.  This is **not** a mathematical
  obstruction.

Keeping `SOLVER_LIMITED` separate from `OBSTRUCTION` is deliberate: a
computer-algebra failure must not be reported as a non-integrability
theorem.

## Scope and limitations

- The algorithms are local.  The base point `x0` must satisfy
  $\delta(x_0)\ne0$, and `v_value` must lie in the domain of $\phi$.
- The symbolic simplification in `rational_curvature` and `refine_delta`
  are heuristics that help SymPy recognise a rational curvature and choose
  the branch adapted to $x>0$.  They are not part of the mathematical
  method and every accepted candidate is checked.
- The Kovacic solver implements Cases 1--3.  Its polynomial systems are
  solved with SymPy; when SymPy cannot complete a subproblem the result is
  `SOLVER_LIMITED`, not an obstruction.
- Numerical checks use the real symbols `x`, `u` and `v`; expressions with
  branch cuts are evaluated on the branch selected by SymPy and should be
  read with the paper's domain restrictions in mind.

## Citation

See `CITATION.cff`.  If you use this software, please cite the paper and
this repository.

## License and attribution

The code is released under the MIT License; see `LICENSE`.  The Kovacic
solver follows Kovacic (1986) and was written with reference to the
exposition and public Maple package of N. M. Abbasi (arXiv:2211.00804); see
`NOTICE`.
