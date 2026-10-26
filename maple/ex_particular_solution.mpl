##############################################################################
# ex_particular_solution.mpl -- Algorithm 2 on an affine Euler-Cauchy ODE
#
#   u' = 2u/x + x^3,   kappa = -2/x^2,   c = 5 x^2,
#   f(x) = x^4/2 is a particular solution,
#   delta_f = exp(int phi_u(t,f(t)) dt) = x^2.
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

NULL:
