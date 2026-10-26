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
| `examples/ex_circle.py` | Circle equation `u' = sqrt(1-u^2)/2` (`ex:nonlinear_proj`, equation `eq:example_sqrt`) |
| `examples/ex_airy_obstruction.py` | Airy obstruction (`ex:airy_obstruction`) |
| `examples/ex_nomaple.py` | Algorithm 3 on the algebraic radical equation (`ex:nomaple`) |

Maple counterparts of the three integration algorithms and all examples are
provided in [`maple/`](maple/). They implement the same diagnostics
(`outside_class`, `invalid_auxiliary_solution`, `invalid_particular_solution`,
`not_rational`, `obstruction`, `solver_limited`), use Maple's
`DEtools[kovacicsols]` for Kovacic's algorithm, and verify themselves with
one `[pass]`/`[FAIL]` line per identity. They were tested with Maple 2015.0;
`maple/transcript.txt` is the recorded output of `sh maple/run_all.sh`
(121 checks, all passing). See `maple/README.md` for details.

Run one with, for example,

```bash
python examples/ex_galois_positive.py
```

Each script prints the curvature, the coefficient $c$, the auxiliary
solution, the first integral and the local family, and verifies the
identities used in the paper.

## Relation to the paper

The three procedures follow §8 literally; three presentation choices are
worth knowing when comparing printed output with the paper.

**Family parameter.** The procedures return `u(x;v)` with `v = u(x0)`,
whereas the examples in the paper write `f(x;s)` with their own parameter.
For the two instances where the paper uses `s`:

- `ex:euler_cauchy`: the paper's `f(x;s) = s x^2 + s^{-1} x^{-1}` and the
  returned family coincide for `s = (v + sqrt(v^2 - 8))/8` at `x0 = 2`
  (checked in `examples/ex_euler_cauchy.py`);
- `ex:galois_positive`: they coincide for
  `s = exp(-2)(v + sqrt(v^2 - 4))/2` at `x0 = 1`.

**First integrals are not unique.** Multiplying `delta` by a constant
rescales `I` and `S`; replacing `delta` by the independent solution changes
`I` nonlinearly; and any nonconstant function of a first integral is again a
first integral. In `ex:galois_positive` the paper prints the first integral
attached to `delta = x^{1/4} exp(-2 sqrt(x))`, while `alg:delta` applied to
the displayed `x^{1/4} exp(2 sqrt(x))` returns its reciprocal up to the
factor 4 (`I_+ * I = 4`). `examples/ex_galois_positive.py` checks both
statements. Likewise `examples/ex_euler_cauchy.py` uses `delta = x^2`,
whereas the continuation of that example in the paper derives
`delta = x^2 - x^{-1}` from the particular solution with `s = 1`; both are
nonzero solutions of the same `L`, as the construction allows.

**Kovacic, step 4.** The implementation always takes the conservative
branch: it returns the implicit solution and the family, but makes no claim
about Liouvillianity of the solutions, which is what `alg:kovacic` asks for
when it is not known that `c` is Liouvillian over `C(x)`.

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
  $\delta(x_0)\ne0$; violating this is reported as
  `invalid_auxiliary_solution`, not as a silent failure.  `v_value` must
  lie in the domain of $\phi$, and membership of $(x,\tilde f(x))$ in the
  domain — which `alg:particular` step 2 asks to check — cannot be decided
  symbolically and remains the caller's responsibility.
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
