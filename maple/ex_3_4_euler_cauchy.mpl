##############################################################################
# ex_3_4_euler_cauchy.mpl -- Example 3.4 of the paper (ex:euler_cauchy)
#
#   u' = u/(2x) + (3/(2x)) sqrt(u^2 - 4x),   kappa = -2/x^2,
#   delta = x^2,  c = 0.
#
# Run from this directory:  maple -q < ex_3_4_euler_cauchy.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

phi := u/(2*x) + 3/(2*x)*sqrt(u^2-4*x):
x0 := 2: v0 := 5:

# Numerical zero test inside the domain of the equation: x > 0, u^2 > 4x.
points := [{x = 5/2, u = 5}, {x = 3, u = 5}, {x = 2, u = 4}]:

printf("=== A1: alg:delta with delta = x^2 ===\n"):
out1 := CK_alg_delta(phi, x^2, x0, v0, points):
CK_show(out1):
CK_check(out1, phi, x0, v0, [5/2, 3], points):

printf("\n=== A3: alg:kovacic ===\n"):
out3 := CK_alg_kovacic(phi, x0, v0, points):
CK_show(out3):
CK_check(out3, phi, x0, v0, [5/2, 3], points):

NULL:
