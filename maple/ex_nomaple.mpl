##############################################################################
# ex_nomaple.mpl -- paper example ex:nomaple
#
#   The algebraic ODE used in the paper as the benchmark on which Maple's
#   dsolve returns NULL.  The right-hand side is written as
#
#       phi = 9 x S^2 - u/x,
#
#   where S is the algebraic element defined by
#
#       S^3 + 3 x^3 S^2 + 3 S - x u = 0.
#
#   Using RootOf for S lets Maple differentiate phi through the implicit
#   function theorem (the Python package uses a custom Function for the same
#   purpose).  The curvature is kappa = -2/x^2, c = 0, and the paper's
#   family is  u = 3 s^2 x^2 + (s^3 + 3 s)/x,  which corresponds to the
#   branch S = s.
#
# Run from this directory:  maple -q < ex_nomaple.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

S := RootOf(_Z^3 + 3*x^3*_Z^2 + 3*_Z - x*u, _Z):
phi := 9*x*S^2 - u/x:
x0 := 1/2: v0 := 1:

printf("=== A3: alg:kovacic on the algebraic ODE ===\n"):
out := CK_alg_kovacic(phi, x0, v0):
CK_show(out):
CK_check(out, phi, x0, v0, [1, 3/2]):

printf("\n=== the family printed in the paper ===\n"):
fpaper := 3*s^2*x^2 + (s^3+3*s)/x:
CK_report("3s^2x^2+(s^3+3s)/x solves the ODE on the branch S = s",
    CK_zero(subs({S = s, u = fpaper}, diff(fpaper, x) - phi))):
CK_report("the branch S = s corresponds to that family",
    CK_zero(subs({S = s, u = fpaper}, S^3 + 3*x^3*S^2 + 3*S - x*u))):
CK_report("the first integral equals -9 s^2 on the branch S = s",
    CK_zero(subs(S = s, out["I"]) + 9*s^2)):

# The parameter of Algorithm 1 is v = u(x0); on the branch S = s, u(x0) is
# given by fpaper(x0), so the two families have to coincide at every x.
v_of_s := subs(x = x0, fpaper):
fam_s := subs(v = v_of_s, out["family"]):
okp := true:
for sv in [1, 2] do
    for xv in [1, 3/2] do
        okp := okp and abs(evalf(subs({s = sv, x = xv}, fam_s - fpaper), 20)) < 1e-10:
    end do:
end do:
CK_report("Algorithm 1 family = paper family for s = 1, 2", okp):

NULL:
