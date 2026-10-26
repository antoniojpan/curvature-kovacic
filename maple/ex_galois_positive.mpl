##############################################################################
# ex_galois_positive.mpl -- paper example ex:galois_positive
#
#   u' = u/(4x) + sqrt(u^2 - 4 sqrt(x))/sqrt(x),   x > 0,
#   kappa = 3/(16 x^2) - 1/x,
#   delta = x^(1/4) exp(2 sqrt(x)),   c = 0.
#
# Run from this directory:  maple -q < ex_galois_positive.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

phi := u/(4*x) + sqrt(u^2-4*sqrt(x))/sqrt(x):
x0 := 1: v0 := 3:

printf("=== A1: alg:delta with delta = x^(1/4) exp(2 sqrt(x)) ===\n"):
out1 := CK_alg_delta(phi, x^(1/4)*exp(2*sqrt(x)), x0, v0):
CK_show(out1):
CK_check(out1, phi, x0, v0, [3/2, 2]):

printf("\n=== A3: alg:kovacic (Kovacic's Example 1, Section 4.2) ===\n"):
out3 := CK_alg_kovacic(phi, x0, v0):
CK_show(out3):
CK_check(out3, phi, x0, v0, [3/2, 2]):

NULL:
