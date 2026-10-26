##############################################################################
# ex_airy_obstruction.mpl -- paper example ex:airy_obstruction
#
#   The Airy equation  y'' + x y = 0  has no nonzero Liouvillian solution
#   over C(x); the nonlinear equation  u' = (y'/y) u  has curvature
#   kappa = x and is affine, so Algorithm 1 still integrates it.
#
# Run from this directory:  maple -q < ex_airy_obstruction.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

delta := AiryAi(-x):
phi := diff(delta, x)/delta*u:
x0 := 0: v0 := 1:

printf("=== A1: alg:delta with delta = Ai(-x) ===\n"):
out1 := CK_alg_delta(phi, delta, x0, v0):
CK_show(out1):
CK_check(out1, phi, x0, v0, [1, 2]):

printf("\n=== A3: alg:kovacic (must report the obstruction) ===\n"):
out3 := CK_alg_kovacic(phi, x0, v0):
CK_show(out3):
CK_report("status is an obstruction", out3["status"] = "obstruction"):

NULL:
