##############################################################################
# curvature_kovacic.mpl
#
# Maple counterpart of the constructive integration algorithms of
#
#   A. J. Pan-Collantes, J. A. Alvarez-Garcia,
#   "From curvature to Kovacic: a geometric approach to integrability of
#    scalar ODEs",
#
# Section 8 (algorithms alg:delta, alg:particular, alg:kovacic).
#
# Load with:  read "curvature_kovacic.mpl":
#
# Conventions (as in the paper):
#
#   u' = phi(x,u)
#   A       = d/dx + phi d/du
#   K(x,u)  = -d/du( A(phi) )                  (Gauss curvature)
#   K = kappa(x)  =>  A(phi) = -kappa(x) u + c(x)
#   L       = d^2/dx^2 + kappa(x)
#
# The independent variable, the dependent variable, the parameter of the
# solution family and the integration variable are the global names x, u, v
# and t.  They are reset below so that a previous assignment in the session
# cannot corrupt a computation.
#
# Every procedure returns a table with string keys; use CK_show to print it
# and CK_check to verify it (see the bottom of this file).  All status
# strings are:
#
#   "ok"                        a first integral and a family were produced
#   "outside_class"             d_u K is not identically zero
#   "invalid_auxiliary_solution" delta is zero, or L(delta) <> 0, or delta(x0)=0
#   "invalid_particular_solution" f' <> phi(x,f)
#   "not_rational"              kappa is not in C(x); Kovacic does not apply
#   "obstruction"               no nonzero Liouvillian solution over C(x)
#   "solver_limited"            kovacicsols did not complete or its answer is
#                               not Liouvillian: NOT a mathematical obstruction
#   "solver_error"              kovacicsols answered but the check failed
#
# Tested with Maple 2015.0 (X86 64 LINUX, Build ID 1022128).
##############################################################################

x := 'x': u := 'u': v := 'v': t := 't':

# Try to disable Kovacic's "pullback" presentation (Maple 2018 and later
# return hypergeometric answers for projective groups A4, S4, A5, which are
# NOT Liouvillian).  In older versions this name is simply unused.
try
    _Env_kovacicsols_pullback := false:
catch:
end try:

# Non-Liouvillian function heads that a solver may return for Kovacic Case 4.
CK_nonLiouvillian := [hypergeom, KummerM, KummerU, WhittakerM, WhittakerW,
    BesselJ, BesselY, BesselI, BesselK, AiryAi, AiryBi, HeunG, Heun, hankelh1,
    hankelh2, LaguerreL, LegendreP, ChebyshevT, ChebyshevU, JacobiP]:

# Sample points used by the numerical zero test.
CK_sample := [{x = 2, u = 3/2}, {x = 3/2, u = 6/5}, {x = 5/2, u = 7/3}]:

##############################################################################
# Helpers
##############################################################################

# CK_zero(e): true when e is identically zero.  Symbolic attempts first;
# as a last resort the expression is evaluated numerically at CK_sample.
CK_zero := proc(e)
    local s, p, val, ok;
    if e = 0 then return true end if:
    s := simplify(e):
    if s = 0 then return true end if:
    try
        s := normal(radnormal(s)):
        if s = 0 then return true end if:
    catch:
    end try:
    ok := true:
    for p in CK_sample do
        try
            val := evalf(subs(p, s), 20):
            if not type(val, numeric) or abs(val) > 1e-10 then
                ok := false: break:
            end if:
        catch:
            ok := false: break:
        end try:
    end do:
    return ok and CK_sample <> [];
end proc:

# CK_report(name, cond): print one verification line and return cond.
CK_report := proc(name, cond)
    if cond then
        printf("  [pass] %s\n", name):
    else
        printf("  [FAIL] %s\n", name):
    end if:
    return cond:
end proc:

# CK_show(T): print the entries of a result table.
CK_show := proc(T)
    local k;
    if not type(T, table) then
        printf("%a\n", T):
        return NULL:
    end if:
    for k in sort(map(convert, [indices(T, nolist)], string)) do
        printf("  %-14s %a\n", k, T[k]):
    end do:
    return NULL:
end proc:

# CK_int(g, x0, x): integral of g(t) from t = x0 to t = x, computed from the
# indefinite integral so that Maple does not rewrite a definite integral on a
# complex branch (this happens, e.g., for csc); the limits are substituted by
# hand and the result is simplified.  Returns an inert Int when the
# antiderivative is not available.
CK_int := proc(g, x0, xto)
    local F;
    F := int(subs(x = t, g), t):
    if has(F, int) or has(F, Int) then
        return Int(subs(x = t, g), t = x0 .. xto):
    end if:
    return simplify(subs(t = xto, F) - subs(t = x0, F)):
end proc:

# CK_is_u_free(K): numeric test that K does not depend on u.
CK_is_u_free := proc(K)
    local ref, q, val;
    try
        ref := evalf(subs({x = 2, u = 0}, K), 25):
    catch:
        return false:
    end try:
    if not type(ref, numeric) then return false end if:
    for q in [1, 3/2, 6/5] do
        try
            val := evalf(subs({x = 2, u = q}, K), 25):
        catch:
            return false:
        end try:
        if not type(val, numeric) or abs(val - ref) > 1e-12 then
            return false:
        end if:
    end do:
    return true:
end proc:

##############################################################################
# Geometric coefficients
##############################################################################

# CK_divergence(phi): A(phi) = phi_x + phi phi_u, simplified once, before
# any differentiation in u.  For equations written with nested radicals this
# is what keeps the curvature computation tractable: the identity
# A(phi) = -kappa u + c is a simplification of A(phi) itself.
CK_divergence := proc(phi)
    local Aphi;
    Aphi := simplify(diff(phi, x) + phi*diff(phi, u)):
    if has(Aphi, u) then
        try
            Aphi := radnormal(Aphi):
        catch:
        end try:
    end if:
    return Aphi:
end proc:

# CK_curvature(phi): K = -d_u( phi_x + phi phi_u ).
# Normalisations tried, in order: simplify, radnormal, evala(Normal), and
# finally the u -> 0 reduction (accepted only after a numerical test of
# u-independence).  A curvature that genuinely depends on u is returned
# unchanged, so that has(K,u) detects it.
CK_curvature := proc(phi)
    local K, cand;
    K := simplify(-diff(CK_divergence(phi), u)):
    if not has(K, u) then return K end if:
    try
        cand := radnormal(K):
        if not has(cand, u) then return cand end if:
    catch:
    end try:
    try
        cand := evala(Normal(K)):
        if not has(cand, u) then return cand end if:
    catch:
    end try:
    if not CK_is_u_free(K) then return K end if:
    cand := simplify(subs(u = 0, K)):
    if has(cand, u) then return K end if:
    return cand:
end proc:

# CK_c(phi, kappa): c(x) = A(phi) + kappa u.
CK_c := proc(phi, kappa)
    return simplify(CK_divergence(phi) + kappa*u):
end proc:

##############################################################################
# Algorithm 1 (alg:delta)
##############################################################################

CK_alg_delta := proc(phi, delta, x0, v0)
    local kappa, cc, S, Q, delta2, S2, Q2, F1, F2, up, fam, Sreg, Freg,
        beta, deg, Fdeg, famdeg;

    # (1) curvature and class check
    kappa := CK_curvature(phi):
    if has(kappa, u) then
        return table(["status" = "outside_class", "kappa" = kappa,
            "message" = "d_u K is not identically zero"]):
    end if:
    cc := CK_c(phi, kappa):

    # (2) validate the auxiliary solution and the base point
    if delta = 0 then
        return table(["status" = "invalid_auxiliary_solution",
            "kappa" = kappa, "c" = cc, "delta" = delta,
            "message" = "delta is identically zero"]):
    end if:
    if not CK_zero(diff(delta, x, x) + kappa*delta) then
        return table(["status" = "invalid_auxiliary_solution",
            "kappa" = kappa, "c" = cc, "delta" = delta,
            "message" = "L(delta) is not zero"]):
    end if:
    if CK_zero(subs(x = x0, delta)) then
        return table(["status" = "invalid_auxiliary_solution",
            "kappa" = kappa, "c" = cc, "delta" = delta,
            "message" = "delta(x0) = 0: choose another base point"]):
    end if:

    # (3) multiplier; the degenerate case S = 0 means phi is affine in u
    S := simplify(diff(delta, x) - diff(phi, u)*delta):
    if CK_zero(S) then
        beta := simplify(phi - (diff(delta, x)/delta)*u):
        if not CK_zero(diff(beta, u)) then
            return table(["status" = "invalid_auxiliary_solution",
                "kappa" = kappa, "c" = cc, "delta" = delta, "S" = S,
                "message" = "S = 0 but phi is not affine in u"]):
        end if:
        deg := CK_int(beta/delta, x0, x):
        Fdeg := simplify(u/delta - deg):
        famdeg := simplify(delta*(v/subs(x = x0, delta) + deg)):
        return table(["status" = "ok", "branch" = "degenerate",
            "kappa" = kappa, "c" = cc, "delta" = delta, "S" = 1/delta,
            "I" = Fdeg, "family" = famdeg, "beta" = beta, "S_raw" = S]):
    end if:

    # (4) quadratures and the two first integrals
    Q := CK_int(delta*cc, x0, x):
    delta2 := simplify(delta*CK_int(1/delta^2, x0, x)):
    S2 := simplify(diff(delta2, x) - diff(phi, u)*delta2):
    Q2 := CK_int(delta2*cc, x0, x):
    F1 := simplify(diff(delta, x)*u - delta*phi + Q):
    F2 := simplify(diff(delta2, x)*u - delta2*phi + Q2):
    up := simplify(-delta*Q2 + delta2*Q):
    fam := simplify(up + (v/subs(x = x0, delta))*delta
        + (subs(x = x0, delta)*subs(u = v, subs(x = x0, phi))
            - subs(x = x0, diff(delta, x))*v)*delta2):

    # (5) if S vanishes at the base point, use the independent solution
    if CK_zero(subs(u = v0, subs(x = x0, S))) then
        Sreg := S2: Freg := F2:
    else
        Sreg := S: Freg := F1:
    end if:

    return table(["status" = "ok", "branch" = "general",
        "kappa" = kappa, "c" = cc, "delta" = delta, "S" = Sreg, "I" = Freg,
        "Q" = Q, "delta2" = delta2, "u_p" = up, "family" = fam,
        "S_raw" = S, "I_raw" = F1, "I_alt" = F2, "S_alt" = S2]):
end proc:

##############################################################################
# Algorithm 2 (alg:particular)
##############################################################################

CK_alg_particular := proc(phi, f, x0, v0)
    local kappa, cc, phi_u, delta_f, delta_f2, out, fam, f0, fu0, pu0;

    # (1) curvature and class check
    kappa := CK_curvature(phi):
    if has(kappa, u) then
        return table(["status" = "outside_class", "kappa" = kappa,
            "message" = "d_u K is not identically zero"]):
    end if:
    cc := CK_c(phi, kappa):

    # (2) validate the particular solution
    if not CK_zero(diff(f, x) - subs(u = f, phi)) then
        return table(["status" = "invalid_particular_solution",
            "kappa" = kappa, "c" = cc,
            "message" = "f' is not phi(x,f)"]):
    end if:

    # (3) linear solution from the divergence along f
    phi_u := diff(phi, u):
    delta_f := simplify(exp(CK_int(subs(u = f, phi_u), x0, x))):
    delta_f2 := simplify(CK_int(1/delta_f^2, x0, x)):
    if not CK_zero(diff(delta_f, x, x) + kappa*delta_f) then
        return table(["status" = "invalid_auxiliary_solution",
            "kappa" = kappa, "c" = cc, "delta" = delta_f,
            "message" = "L(delta_f) is not zero"]):
    end if:

    # (4) implicit solution from Algorithm 1
    out := CK_alg_delta(phi, delta_f, x0, v0):
    if out["status"] <> "ok" then
        return out:
    end if:

    # (5) the family printed in the paper (sequential substitutions!)
    f0 := subs(x = x0, f):
    fu0 := subs(x = x0, diff(f, x)):
    pu0 := subs(u = f0, subs(x = x0, phi_u)):
    fam := simplify(f + (v - f0)*delta_f
        + (subs(u = v, subs(x = x0, phi)) - fu0 - pu0*(v - f0))*delta_f2):

    out["algorithm"] := "particular":
    out["delta_f"] := delta_f:
    out["delta_f2"] := delta_f2:
    out["family_alg1"] := out["family"]:
    out["family"] := fam:
    out["family_matches"] := CK_zero(simplify(fam - out["family_alg1"])):
    return out:
end proc:

##############################################################################
# Algorithm 3 (alg:kovacic)
##############################################################################

# CK_kovacic_sols(kappa): list of candidates returned by DEtools[kovacicsols],
# or FAIL when the call did not complete.
CK_kovacic_sols := proc(kappa)
    local y, sols;
    y := 'y':
    try
        sols := timelimit(120,
            DEtools[kovacicsols](diff(y(x), x, x) + kappa*y(x) = 0, y(x))):
    catch:
        return FAIL:
    end try:
    if not type(sols, list) then
        if sols = NULL then return [] end if:
        return [sols]:
    end if:
    return sols:
end proc:

CK_alg_kovacic := proc(phi, x0, v0)
    local kappa, cc, sols, cand, out, tried, nliv;

    # (1) curvature and class check
    kappa := CK_curvature(phi):
    if has(kappa, u) then
        return table(["status" = "outside_class", "kappa" = kappa,
            "message" = "d_u K is not identically zero"]):
    end if:
    cc := CK_c(phi, kappa):

    # (2) Kovacic's algorithm requires a rational coefficient
    if not type(kappa, ratpoly(anything, x)) then
        return table(["status" = "not_rational", "kappa" = kappa, "c" = cc,
            "message" = "Kovacic's algorithm does not apply over C(x)"]):
    end if:

    # (3) run Kovacic
    sols := CK_kovacic_sols(kappa):
    if sols = FAIL then
        return table(["status" = "solver_limited", "kappa" = kappa, "c" = cc,
            "message" = "kovacicsols did not complete within the time limit"]):
    end if:
    if nops(sols) = 0 then
        # For rational coefficients kovacicsols is complete: Case 4.
        return table(["status" = "obstruction", "kappa" = kappa, "c" = cc,
            "message" = "the nonlinear ODE has at most one Liouvillian solution over C(x)"]):
    end if:

    # (4) keep the Liouvillian candidates that solve L(delta) = 0, and run
    #     Algorithm 1 with the first one that is nonzero at the base point
    tried := []:
    nliv := 0:
    out := FAIL:
    for cand in sols do
        if has(cand, CK_nonLiouvillian) then next end if:
        nliv := nliv + 1:
        if not CK_zero(diff(cand, x, x) + kappa*cand) then next end if:
        tried := [op(tried), cand]:
        out := CK_alg_delta(phi, cand, x0, v0):
        if out["status"] = "ok" then
            out["message"] := "kovacic_success":
            return out:
        end if:
    end do:
    if nliv = 0 then
        # the solver answered with non-Liouvillian expressions only
        return table(["status" = "solver_limited", "kappa" = kappa, "c" = cc,
            "solutions" = sols,
            "message" = "no Liouvillian solution recovered: not an obstruction"]):
    end if:
    if tried = [] then
        # the solver answered, but none of its answers solves L(delta) = 0
        return table(["status" = "solver_error", "kappa" = kappa, "c" = cc,
            "solutions" = sols, "message" = "kovacicsols answers do not solve L(delta) = 0"]):
    end if:
    # every candidate failed Algorithm 1 (typically delta(x0) = 0)
    out["solutions_tried"] := tried:
    out["message"] := "kovacicsols succeeded but Algorithm 1 rejected every candidate":
    return out:
end proc:

##############################################################################
# Verification of a result table
##############################################################################

# CK_check(out, phi, x0, v0, xpts): check the identities used in the paper.
# xpts is a list of x-values at which the family is tested numerically.
CK_check := proc(out, phi, x0, v0, xpts)
    local ok, S, F, Fx, xv, uv, res, Ivals, ref, kk, st;
    ok := true:

    if type(out, table) then
        st := out["status"]:
    else
        st := out:
    end if:
    if not type(out, table) or st <> "ok" then
        printf("  [FAIL] status = %a\n", st):
        return false:
    end if:

    # multiplier and first integral: dI = S (du - phi dx)
    S := out["S"]:
    F := out["I"]:
    ok := CK_report("dI/du = S", CK_zero(diff(F, u) - S)) and ok:
    ok := CK_report("dI/dx + S*phi = 0", CK_zero(diff(F, x) + S*phi)) and ok:
    ok := CK_report("L(delta) = 0",
        CK_zero(diff(out["delta"], x, x) + out["kappa"]*out["delta"])) and ok:

    # explicit family: initial conditions and the ODE itself
    Fx := out["family"]:
    ok := CK_report("family(x0) = v", CK_zero(subs(x = x0, Fx) - v)) and ok:
    ok := CK_report("family'(x0) = phi(x0,v)",
        CK_zero(subs(x = x0, diff(Fx, x)) - subs(u = v, subs(x = x0, phi)))) and ok:
    Ivals := []:
    for xv in xpts do
        uv := evalf(subs({x = xv, v = v0}, Fx), 25):
        res := evalf(subs({x = xv, u = uv, v = v0}, diff(Fx, x) - phi), 25):
        ok := CK_report(sprintf("family solves the ODE at x = %a", xv),
            type(res, numeric) and abs(res) < 1e-8) and ok:
        Ivals := [op(Ivals), evalf(subs({u = uv, x = xv, v = v0}, F), 25)]:
    end do:
    ref := Ivals[1]:
    ok := CK_report("I is constant along the family",
        evalb(max(seq(abs(q - ref), q in Ivals)) < 1e-8)) and ok:

    # particular solution produced by the quadratures
    if has([indices(out, nolist)], "u_p") then
        kk := out["u_p"]:
        ok := CK_report("L(u_p) = c",
            CK_zero(diff(kk, x, x) + out["kappa"]*kk - out["c"])) and ok:
        ok := CK_report("u_p(x0) = u_p'(x0) = 0",
            CK_zero(subs(x = x0, kk)) and CK_zero(subs(x = x0, diff(kk, x)))) and ok:
    end if:

    if has([indices(out, nolist)], "family_matches") then
        ok := CK_report("paper family = Algorithm 1 family",
            out["family_matches"]) and ok:
    end if:

    if ok then
        printf("  ALL CHECKS PASSED\n"):
    else
        printf("  SOME CHECKS FAILED\n"):
    end if:
    return ok:
end proc:

NULL:
