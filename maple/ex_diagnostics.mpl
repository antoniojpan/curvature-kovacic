##############################################################################
# ex_diagnostics.mpl -- the rejection branches of the three algorithms
#
#   A companion implementation is only informative if it can also say "no".
#   This script exercises every diagnostic status of curvature_kovacic.mpl.
#
# Run from this directory:  maple -q < ex_diagnostics.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

CK_expect := proc(lab, out, expected)
    local st, okflag;
    if type(out, table) then st := out["status"]; else st := out; end if:
    okflag := evalb(st = expected):
    printf("%-46s -> %-28s ", lab, st):
    CK_report("as expected", okflag):
end proc:

printf("=== outside the class: d_u K <> 0 ===\n"):
CK_expect("alg:delta, phi = u^2", CK_alg_delta(u^2, 1, 1, 1), "outside_class"):
CK_expect("alg:particular, phi = u^2", CK_alg_particular(u^2, x, 1, 1), "outside_class"):
CK_expect("alg:kovacic, phi = u^2", CK_alg_kovacic(u^2, 1, 1), "outside_class"):

printf("\n=== invalid auxiliary solution ===\n"):
# kappa = 1/4, but sin(x) does not solve y'' + y/4 = 0
CK_expect("alg:delta, delta = sin(x) is not a solution",
    CK_alg_delta(sqrt(1-u^2)/2, sin(x), 0, 0), "invalid_auxiliary_solution"):
# cos(x/2) solves the linear equation but vanishes at x0 = Pi
CK_expect("alg:delta, delta(x0) = 0",
    CK_alg_delta(sqrt(1-u^2)/2, cos(x/2), Pi, 0), "invalid_auxiliary_solution"):
CK_expect("alg:delta, delta = 0",
    CK_alg_delta(sqrt(1-u^2)/2, 0, 0, 0), "invalid_auxiliary_solution"):

printf("\n=== invalid particular solution ===\n"):
# f = x^3 does not solve u' = 2u/x + x^3
CK_expect("alg:particular, f = x^3",
    CK_alg_particular(2*u/x + x^3, x^3, 1, 2), "invalid_particular_solution"):

printf("\n=== curvature not rational: Kovacic does not apply ===\n"):
# phi = u sqrt(x) has kappa = -(x + 1/(2 sqrt(x)))
CK_expect("alg:kovacic, kappa not in C(x)",
    CK_alg_kovacic(u*sqrt(x), 1, 1), "not_rational"):

printf("\n=== obstruction (Maple 2015: kovacicsols is complete) ===\n"):
CK_expect("alg:kovacic, Airy obstruction (kappa = x)",
    CK_alg_kovacic(diff(AiryAi(-x), x)/AiryAi(-x)*u, 0, 1), "obstruction"):

NULL:
