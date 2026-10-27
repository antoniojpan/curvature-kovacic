# Maple implementation

Maple counterpart of the constructive integration algorithms of the paper
(§8): `alg:delta`, `alg:particular` and `alg:kovacic`.

| File | Content |
|---|---|
| `curvature_kovacic.mpl` | the library: the three algorithms, the curvature and coefficient helpers, and the verification procedure |
| `ex_circle.mpl` | `ex:nonlinear_proj` (equation `eq:example_sqrt`), $u'=\sqrt{1-u^2}/2$, $\kappa=1/4$ |
| `ex_euler_cauchy.mpl` | `ex:euler_cauchy`, $\kappa=-2/x^2$ |
| `ex_galois_positive.mpl` | `ex:galois_positive`, Kovacic's Example 1 of §4.2 |
| `ex_airy_obstruction.mpl` | `ex:airy_obstruction`: Algorithm 1 integrates, Algorithm 3 must report the obstruction |
| `ex_particular_solution.mpl` | `alg:particular` on an affine Euler–Cauchy equation |
| `ex_nomaple.mpl` | `ex:nomaple`, the algebraic equation on which Maple's `dsolve` returns NULL, with its nested radicals exactly as printed |
| `ex_kovacic_solver.mpl` | diagnostics for `DEtools[kovacicsols]` itself |
| `ex_diagnostics.mpl` | the rejection branches: outside the class, invalid auxiliary or particular solution, non-rational curvature, obstruction |
| `run_all.sh` | runs every example in turn |
| `transcript.txt` | recorded output of `run_all.sh` on Maple 2015.0 |

## Requirements and testing

Only the Maple standard library is used; Kovacic's algorithm is taken from
`DEtools[kovacicsols]`. No package needs to be installed.

The scripts were tested with **Maple 2015.0** (X86 64 LINUX, Build ID
1022128). `transcript.txt` is the recorded output: 120 checks, all passing,
including the identity $dI=S_\delta\,(du-\phi\,dx)$, $L(\delta)=0$, the
initial conditions of the family, the ODE along the family, the
constancy of $I$ along the family and $L(u_p)=c$.

## Usage

From this directory:

```sh
sh run_all.sh                       # every example, with its checks
maple -q < ex_circle.mpl            # one example
```

or set `MAPLE` to your binary: `MAPLE=/opt/maple/bin/maple sh run_all.sh`.
The example scripts `read "curvature_kovacic.mpl"` from the current
directory, so run them from here (or `currentdir` to this folder first).

Every script ends with one `[pass]`/`[FAIL]` line per identity, and the
library procedure `CK_check` performs the standard verification of a result:

```maple
read "curvature_kovacic.mpl":
phi := sqrt(1-u^2)/2:
out := CK_alg_delta(phi, cos(x/2), 0, 0):
CK_show(out);
CK_check(out, phi, 0, 0, [1/2, 1]);
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
  needs it even for `ex:galois_positive`.
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
* **Algorithm 3 and Case 4.** For rational $\kappa$, `kovacicsols` is
  complete in Maple 2015, so the empty list means Case 4. In later versions
  the default "pullback" presentation returns `hypergeom` answers for
  projective groups $A_4,S_4,A_5$, which are *not* Liouvillian; the library
  tries to switch that off (`_Env_kovacicsols_pullback := false`), rejects
  any candidate carrying a non-Liouvillian head, checks that each remaining
  candidate solves $L(\delta)=0$, and reports `solver_limited` if nothing
  survives. The solver call is wrapped in `timelimit(120, ...)`.
* **`ex_nomaple`.** The nested radicals are used exactly as printed; the
  algebraic element is *not* replaced by a `RootOf`. That replacement would
  turn the example into a different and much easier problem: on the `RootOf`
  presentation Maple's `dsolve` returns the general solution in about 0.1 s,
  whereas on the printed input it does not return at all, so the benchmark
  recorded in the paper would be lost. What makes the printed input
  tractable here is the order of operations in `CK_curvature`: $A(\phi)$ is
  simplified first (it becomes $2u/x^2$), and only then differentiated in
  $u$; `radnormal` and `evala(Normal(...))` handle the residual radicals.
  The closed forms $I=-9S^2$ and $3s^2x^2+(s^3+3s)/x$ are then verified
  against the raw output.
* **Normalization.** Different choices of the linear solution $\delta$ give
  different normalizations of $I$ and $S$, and different bases of the
  solution space in Algorithm 3. The implicit solution is therefore defined
  up to a nonzero constant factor; the family is not, since it is pinned by
  $f(x_0)=v$ and $f'(x_0)=\phi(x_0,v)$.
