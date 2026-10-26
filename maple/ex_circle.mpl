##############################################################################
# ex_circle.mpl -- paper example ex:nonlinear_proj (equation eq:example_sqrt)
#
#   u' = sqrt(1 - u^2)/2,   kappa = 1/4,
#   delta = cos(x/2),  family  u(x) = v cos(x/2) + sqrt(1-v^2) sin(x/2).
#
# Run from this directory:  maple -q < ex_circle.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

phi := sqrt(1-u^2)/2:
x0 := 0: v0 := 0:

printf("=== A1: alg:delta with delta = cos(x/2) ===\n"):
out1 := CK_alg_delta(phi, cos(x/2), x0, v0):
CK_show(out1):
CK_check(out1, phi, x0, v0, [1/2, 1]):

printf("\n=== A3: alg:kovacic ===\n"):
out3 := CK_alg_kovacic(phi, x0, v0):
CK_show(out3):
CK_check(out3, phi, x0, v0, [1/2, 1]):

NULL:
