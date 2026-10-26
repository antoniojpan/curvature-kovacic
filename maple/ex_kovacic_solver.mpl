##############################################################################
# ex_kovacic_solver.mpl -- diagnostics for DEtools[kovacicsols] itself
#
#   This script does not use the nonlinear algorithms: it only records what
#   Maple's Kovacic implementation returns for the inputs used in the paper
#   and in the Python diagnostics, and checks that every candidate really
#   solves  y'' + kappa y = 0.
#
# Run from this directory:  maple -q < ex_kovacic_solver.mpl
##############################################################################

restart:
read "curvature_kovacic.mpl":

y := 'y':

examples := [
    ["y'' + y = 0",                    1,                        "solved"],
    ["y'' - y = 0",                   -1,                        "solved"],
    ["Euler-Cauchy",                  -2/x^2,                    "solved"],
    ["Kovacic Example 1",              3/(16*x^2) - 1/x,         "solved"],
    ["Case 2",                        (3 + 8*x)/(16*x^2),        "solved"],
    ["Case 3",                        (x - 4)/(4*x*(x - 1)^2),   "solved"],
    ["Airy obstruction",               x,                        "obstruction"],
    ["Bessel nu = 1",                  1 - 3/(4*x^2),            "obstruction"]
]:

for item in examples do
    lab := item[1]:
    kappa := item[2]:
    expected := item[3]:
    sols := CK_kovacic_sols(kappa):
    if sols = FAIL then
        actual := "solver_limited":
    elif nops(sols) = 0 then
        actual := "obstruction":
    else
        actual := "solved":
    end if:
    printf("\n%-22s kappa = %a\n", lab, kappa):
    printf("  expected %-12s actual %-12s\n", expected, actual):
    okflag := evalb(actual = expected):
    CK_report("Kovacic status as expected", okflag):
    for cand in sols do
        CK_report("candidate solves L(y) = 0",
            CK_zero(diff(cand, x, x) + kappa*cand)):
    end do:
end do:

NULL:
