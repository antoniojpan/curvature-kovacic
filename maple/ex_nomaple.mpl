##############################################################################
# ex_nomaple.mpl -- paper example ex:nomaple
#
#   The algebraic ODE used in the paper as the benchmark on which Maple's
#   dsolve returns NULL.  The right-hand side is used exactly as printed,
#   with its nested radicals:
#
#       R = ( (x u + 3 x^3 - 2 x^9
#              + sqrt(x^2 u^2 + 2 x^4 (3 - 2 x^6) u + 4 - 3 x^6)) / 2 )^(1/3)
#       S = R + (x^6 - 1)/R - x^3
#       u' = 9 x S^2 - u/x,        D = {0 < x < 1, u > 0}.
#
#   The algebraic element is deliberately NOT replaced by a RootOf: that
#   rewrites the example into a different and much easier problem (on the
#   RootOf presentation Maple's dsolve returns the general solution in
#   0.1 seconds, whereas on this input it returns NULL in about 1.6
#   seconds on Maple 2015.0), which would make the comparison recorded in
#   the paper meaningless.
#
#   The curvature is nevertheless cheap: A(phi) = phi_x + phi phi_u
#   simplifies to 2u/x^2 before the u-derivative is taken, whence
#   kappa = -2/x^2 and c = 0.  With delta = 1/x,
#
#       I = -9 S^2,     f(x;s) = 3 s^2 x^2 + (s^3 + 3 s)/x,   s > 0.
#
#   The raw output of the algorithm is written in the nested radicals and is
#   too long to be informative, so only its size is reported; the closed
#   forms printed in the paper are verified instead.
#
# Run from this directory:  maple -q < ex_nomaple.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

R := ((x*u + 3*x^3 - 2*x^9
        + sqrt(x^2*u^2 + 2*x^4*(3 - 2*x^6)*u + 4 - 3*x^6))/2)^(1/3):
S := R + (x^6 - 1)/R - x^3:
phi := 9*x*S^2 - u/x:
x0 := 1/2: v0 := 1:

printf("=== A(phi), curvature and coefficient, with the explicit radicals ===\n"):
printf("  A(phi) = %a\n", CK_divergence(phi)):
printf("  kappa  = %a\n", CK_curvature(phi)):
printf("  c      = %a\n", CK_c(phi, CK_curvature(phi))):

printf("\n=== A3: alg:kovacic on the printed input ===\n"):
out := CK_alg_kovacic(phi, x0, v0):
printf("  status = %a   (%a)\n", out["status"], out["message"]):
printf("  delta  = %a\n", out["delta"]):
if out["status"] = "ok" then
    printf("  raw output (nested radicals): S = %a chars, I = %a chars, family = %a chars\n",
        length(out["S"]), length(out["I"]), length(out["family"])):
end if:
CK_check(out, phi, x0, v0, [1/4, 3/4]):

printf("\n=== closed forms printed in the paper ===\n"):
printf("  I      = -9 S^2\n"):
CK_report("I = -9 S^2", CK_zero(out["I"] + 9*S^2)):

# The paper's family, on the domain 0 < x < 1.
printf("  family = 3 s^2 x^2 + (s^3 + 3 s)/x\n"):
fpaper := 3*s^2*x^2 + (s^3 + 3*s)/x:
okp := true:
for sv in [1, 2, 3/2] do
    for xv in [1/4, 3/4] do
        uv := subs({s = sv, x = xv}, fpaper):
        okp := okp and abs(evalf(subs({x = xv, u = uv, s = sv}, diff(fpaper, x) - phi), 20)) < 1e-10:
    end do:
end do:
CK_report("3 s^2 x^2 + (s^3 + 3 s)/x solves the ODE", okp):

# v = u(x0) on the paper's family; the two families must then coincide.
v_of_s := subs(x = x0, fpaper):
fam_s := subs(v = v_of_s, out["family"]):
okp := true:
for sv in [1, 2, 3/2] do
    for xv in [1/4, 3/4] do
        okp := okp and abs(evalf(subs({s = sv, x = xv}, fam_s - fpaper), 20)) < 1e-10:
    end do:
end do:
CK_report("Algorithm 1 family = paper family 3 s^2 x^2 + (s^3+3s)/x", okp):

NULL:
