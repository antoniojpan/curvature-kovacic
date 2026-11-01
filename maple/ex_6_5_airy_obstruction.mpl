##############################################################################
# ex_6_5_airy_obstruction.mpl -- Example 6.5 of the paper (ex:airy_obstruction)
#
#   The Airy equation  y'' + x y = 0  has no nonzero Liouvillian solution
#   over C(x); the nonlinear equation  u' = (y'/y) u  has curvature
#   kappa = x and is affine, so Algorithm 1 still integrates it.
#
# Run from this directory:  maple -q < ex_6_5_airy_obstruction.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

delta := AiryAi(-x):
phi := diff(delta, x)/delta*u:
x0 := 0: v0 := 1:

# Numerical zero test away from the real zeros of Ai(-x) (-2.3381, -4.0879, ...).
points := [{x = 1, u = 1}, {x = 3/2, u = 2}, {x = 2, u = 1}]:

printf("=== A1: alg:delta with delta = Ai(-x) ===\n"):
out1 := CK_alg_delta(phi, delta, x0, v0, points):
CK_show(out1):
CK_check(out1, phi, x0, v0, [1, 2], points):

printf("\n=== A3: alg:kovacic (must report the obstruction) ===\n"):
out3 := CK_alg_kovacic(phi, x0, v0, points):
CK_show(out3):
CK_report("status is an obstruction", out3["status"] = "obstruction"):

NULL:
