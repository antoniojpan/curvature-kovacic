##############################################################################
# ex_particular_solution.mpl -- Algorithm 2 in two regimes
#
# (a) affine Euler-Cauchy ODE (degenerate branch, S_delta = 0)
#
#   u' = 2u/x + x^3,   kappa = -2/x^2,   c = 5 x^2,
#   f(x) = x^4/2 is a particular solution,
#   delta_f = exp(int phi_u(t,f(t)) dt) = x^2.
#
# (b) nonlinear member of the class (general branch, S_delta <> 0)
#
#   u' = sqrt(1-u^2)/2,   kappa = 1/4,   c = 0,
#   f(x) = sin(x/2) is a particular solution,
#   delta_f = cos(x/2) and delta_f2 = delta_f int dt/delta_f^2 = 2 sin(x/2).
#
#   Here the delta_f2 term of the family does not vanish, so this block
#   would fail if delta_f2 were computed without the factor delta_f.
#
# Run from this directory:  maple -q < ex_particular_solution.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

phi := 2*u/x + x^3:
f := x^4/2:
x0 := 1: v0 := 2:

out := CK_alg_particular(phi, f, x0, v0):
CK_show(out):
CK_check(out, phi, x0, v0, [3/2, 2]):

x := 'x': u := 'u': v := 'v':
assume(x > -Pi, x < Pi):

phi := sqrt(1-u^2)/2:
f := sin(x/2):
x0 := Pi/3: v0 := 3/4:

out := CK_alg_particular(phi, f, x0, v0):
CK_show(out):
CK_check(out, phi, x0, v0, [1, 3/2]):

NULL:
