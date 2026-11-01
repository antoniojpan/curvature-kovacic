# Maple implementation

Maple counterpart of the constructive integration algorithms of the paper
(§7): `alg:delta`, `alg:particular` and `alg:kovacic`.

| File | Paper | Content |
|---|---|---|
| `curvature_kovacic.mpl` | §7 | the library: the three algorithms, the curvature and coefficient helpers, and the verification procedure |
| `ex_3_4_euler_cauchy.mpl` | Example 3.4 (`ex:euler_cauchy`) | $\kappa=-2/x^2$ |
| `ex_5_6_circle.mpl` | Example 5.6 (`ex:nonlinear_proj`) | equation `eq:example_sqrt`, $u'=\sqrt{1-u^2}/2$, $\kappa=1/4$ |
| `ex_6_3_galois_positive.mpl` | Example 6.3 (`ex:galois_positive`) | Kovacic's Example 1 of §4.2 |
| `ex_6_5_airy_obstruction.mpl` | Example 6.5 (`ex:airy_obstruction`) | Algorithm 1 integrates, Algorithm 3 must report the obstruction |
| `ex_6_7_nomaple.mpl` | Example 6.7 (`ex:nomaple`) | the algebraic equation on which Maple's `dsolve` returns NULL, with its nested radicals exactly as printed |
| `alg_2_particular_solution.mpl` | Algorithm 7.2 | `alg:particular` on an affine Euler–Cauchy equation and on a nonlinear member of the class (general branch, where $\delta_{f2}$ enters) |
| `alg_3_kovacic_solver.mpl` | Algorithm 7.3 | diagnostics for `DEtools[kovacicsols]` itself |
| `alg_diagnostics.mpl` | §7 | the rejection branches: outside the class, invalid auxiliary or particular solution, non-rational curvature, obstruction |
| `run_all.sh` | — | runs every example in turn |
| `transcript.txt` | — | recorded output of `run_all.sh` on Maple 2015.0 |

## Requirements and testing

Only the Maple standard library is used; Kovacic's algorithm is taken from
`DEtools[kovacicsols]`. No package needs to be installed.

The scripts were tested with **Maple 2015.0** (X86 64 LINUX, Build ID
1022128). `transcript.txt` is the recorded output: 141 checks, all passing,
including the identity $dI=S_\delta\,(du-\phi\,dx)$,
$L(\delta)=L(\delta_2)=L(\delta_{f2})=0$, the
initial conditions of the family, the ODE along the family, the
constancy of $I$ along the family and $L(u_p)=c$.

## Usage

From this directory:

```sh
sh run_all.sh                       # every example, with its checks
maple -q < ex_5_6_circle.mpl       # one example
```

or set `MAPLE` to your binary: `MAPLE=/opt/maple/bin/maple sh run_all.sh`.
The example scripts `read "curvature_kovacic.mpl"` from the current
directory, so run them from here (or `currentdir` to this folder first).

Every script ends with one `[pass]`/`[FAIL]` line per identity, and the
library procedure `CK_check` performs the standard verification of a result:

```maple
read "curvature_kovacic.mpl":
phi := sqrt(1-u^2)/2:
points := [{x = 1/2, u = 0}, {x = 1, u = 1/2}]:
out := CK_alg_delta(phi, cos(x/2), 0, 0, points):
CK_show(out);
CK_check(out, phi, 0, 0, [1/2, 1], points);
```

## Interface

| Paper algorithm | Maple procedure | Python function |
|---|---|---|
| `alg:delta` | `CK_alg_delta(phi, delta, x0, v0)` | `solve_from_linear_solution` |
| `alg:particular` | `CK_alg_particular(phi, f, x0, v0)` | `solve_from_particular` |
| `alg:kovacic` | `CK_alg_kovacic(phi, x0, v0)` | `kovacic_integration` |

The independent variable, the dependent variable, the family parameter and
the integration variable are the global names `x`, `u`, `v` and `t`, as in
the paper.

Every procedure takes an optional last argument `points`: the sample set used
by the numerical zero test (`CK_alg_delta(phi, delta, x0, v0, points)`,
`CK_alg_particular(phi, f, x0, v0, points)`,
`CK_alg_kovacic(phi, x0, v0, points)`,
`CK_check(out, phi, x0, v0, xpts, points)`, `CK_zero(e, points)`). It is the
counterpart of `points` in the Python package, and it must lie inside the
domain of the equation (see the note on sample points below). The default is
`CK_sample`, a set meant for equations with no prescribed domain.

Each procedure returns a table with string keys, printed by `CK_show`.
The `status` entry is one of

| status | meaning |
|---|---|
| `ok` | a first integral and a solution family were produced |
| `outside_class` | $\partial_u\mathcal K\ne0$ |
| `invalid_auxiliary_solution` | $\delta\equiv0$, $L(\delta)\ne0$, or $\delta(x_0)=0$ |
| `invalid_particular_solution` | $\tilde f'\ne\phi(x,\tilde f)$ |
| `not_rational` | $\kappa\notin\mathbb C(x)$: Kovacic's algorithm does not apply |
| `obstruction` | no nonzero Liouvillian solution over $\mathbb C(x)$ (Kovacic Case 4) |
| `solver_limited` | `kovacicsols` did not complete or its answer is not Liouvillian — **not** a mathematical obstruction |
| `solver_error` | `kovacicsols` answered but the returned solution failed the check |

Keeping `solver_limited` apart from `obstruction` is deliberate, and matches
the Python package: a computer-algebra failure must not be reported as a
non-integrability theorem.

## Implementation notes

* **Curvature.** `CK_curvature` simplifies $-\partial_u(\phi_x+\phi\phi_u)$;
  when a spurious dependence on $u$ survives, a `u=0` reduction is attempted
  and accepted only after a numerical test of $u$-independence. This is the
  counterpart of `rational_curvature` in the Python package, where SymPy
  needs it even for Example 6.3.
* **Quadratures.** `CK_int` computes the *indefinite* integral and
  substitutes the limits by hand. Maple's definite `int` may return a
  rewrite on a complex branch (for instance
  $\int_0^x\sec^2(t/2)\,dt$ as $-2i/\sqrt{1-\csc^2(x/2)}$, which is
  undefined at $x=0$); the indefinite form avoids this.
* **Degenerate branch.** If $S_\delta=\delta'-\phi_u\delta$ vanishes
  identically, step 4 of `alg:delta` is used
  ($I=u/\delta-\int\beta/\delta$, $\phi=(\delta'/\delta)u+\beta$). This branch
  is what integrates the Airy example, where the nonlinear equation is
  linear.
* **Algorithm 3 and Case 4.** For rational $\kappa$, an empty list from
  `kovacicsols` is reported as Case 4 (`obstruction`). This is verified on the
  examples above with Maple 2015.0, but it rests on the solver being
  exhaustive: if a version of the CAS missed a Liouvillian solution, the empty
  list would be a false obstruction. The status is therefore only as strong as
  the solver. The mathematical statement of `alg:kovacic` does not depend on
  this — it asks for a nonzero Liouvillian solution, and the two obstructions
  quoted in the paper (Airy, Bessel) are verified independently.
  In the versions that have Kovacic's "pullback" presentation, `kovacicsols`
  presents the answers for projective groups $A_4,S_4,A_5$ as `hypergeom`
  expressions, which the library's non-Liouvillian filter would reject. The
  setting `_Env_kovacicsols_pullback := false` — documented in
  `?DEtools[kovacicsols]` — forces the classical presentation; on Maple 2015
  it is inert, because the presentation does not exist there. Independently
  of that, the library rejects any candidate carrying a non-Liouvillian head,
  checks that each remaining candidate solves $L(\delta)=0$, and reports
  `solver_limited` if nothing survives. The solver call is wrapped in
  `timelimit(120, ...)`.
* **`Example 6.7`.** What makes the input tractable here is the order of
  operations in `CK_curvature`: $A(\phi)$ is simplified first (it becomes
  $2u/x^2$), and only then differentiated in $u$; `radnormal` and
  `evala(Normal(...))` handle the residual radicals.
  The closed forms $I=-9S^2$ and $3s^2x^2+(s^3+3s)/x$ are then verified
  against the raw output.
* **Sample points.** `CK_zero` simplifies symbolically first and evaluates
  numerically only as a last resort. Such a numerical decision is meaningful
  only where the equation and its solutions are real: the default `CK_sample`
  has $x>1$, which is outside the domain of Example 6.7 ($0<x<1$, $u>0$);
  there the radicand $x^2u^2+2x^4(3-2x^6)u+4-3x^6$ is negative and Maple
  evaluates a complex branch, so a true identity such as $I=-9S^2$ is worth
  $10^{-7}$, not $0$, at those points. Every example whose domain excludes
  the default therefore passes its own `points`, inside the domain; for
  Example 6.7 these are the points $x=1/2,1/4,3/4$ of the family
  $3s^2x^2+(s^3+3s)/x$ with $s=1$.
* **Normalization.** Different choices of the linear solution $\delta$ give
  different normalizations of $I$ and $S$, and different bases of the
  solution space in Algorithm 3. The implicit solution is therefore defined
  up to a nonzero constant factor; the family is not, since it is pinned by
  $f(x_0)=v$ and $f'(x_0)=\phi(x_0,v)$.
