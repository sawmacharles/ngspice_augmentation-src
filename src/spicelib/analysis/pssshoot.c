/**********
card #458, plan T4.1-T4.3: shooting-Newton periodic steady state for
DRIVEN (clocked) switching converters -- `.pss Fguess StabTime OscNode
Points Harmonics SC_iter Steady_coeff [uic] [method=shoot]`.

This is a separate driver next to the existing dcpss.c algorithm (which
stays the default and is not modified in its own behaviour). DCpss()
dispatches here when the deck asks for method=shoot and a driven source
was detected (Enhancement-176's own detection, reused as-is). The
free-running/autonomous case is explicitly out of scope for this card
(design ruling #447 A, 2026-09-27) and is refused by the caller before this
function is ever entered.

Design (see research/ngspice-internals-pss-pac.md Sec.4.1 and the plan/spec
pages of the ngspice_augmentation project, for the fuller
derivation):

  - State vector x = the true dynamic states the transient needs to
    restart: one entry per native capacitor (its terminal voltage) and one
    per native inductor (its branch current). This is the "unknowns are
    the capacitor voltages / inductor currents" option the card offered,
    chosen over the "full restart state" option because it lets each
    period map be evaluated with ngspice's own, already-tested `uic`
    initial-condition restart (CAPinitCond/INDinitCond, consumed by
    capload.c/indload.c whenever CKTmode carries MODEUIC|MODEINITTRAN) --
    the same mechanism refs/r3_vm_buck.cir and r4_pcm_buck.cir already
    rely on ("IC=0" on L1/Cout, needed because a plain op lands on a
    degenerate non-switching equilibrium, refs/README.md). Any other
    algebraic node (comparator outputs, op-amp internals, ...) is
    resolved fresh at the start of every period by CKTop() inside
    DCtran()'s own restart path, from whatever `.ic` node guesses the
    deck supplies -- not part of the Newton unknown.

  - One-period map P(x0): write x0 into the CAPinitCond/INDinitCond
    fields, then run a full transient from t=0 to t=T=1/fguess (the known
    clock period) via ngspice's *existing*, fully-tested DCtran() --
    reused directly (not reimplemented) by temporarily pointing
    ckt->CKTcurJob at a stack-local TRANan matching the field layout
    DCtran() actually reads (JOBname, TRANplot) and restoring the real
    PSSan job pointer afterwards. x1 = P(x0) is read back the same way
    x0 was written, from the converged ckt->CKTrhsOld at t=T.

  - Newton on F(x0) = P(x0) - x0: finite-difference monodromy (one
    baseline period run + n perturbed period runs per Newton step, n =
    state count -- cheap here: 1-2 states for R2, a handful for R3/R4),
    solved with a small hand-written dense Gauss elimination (n is never
    more than a handful of states, so this avoids dragging in ngspice's
    sparse solver for what is a tiny dense system), and a damped
    (backtracking) step so a bad initial guess does not diverge outright.

  - Floquet: the spectral radius of the last monodromy matrix Phi, by
    power iteration (the card allows this instead of a full eigenvalue
    solve), plus a cheap sign check (does the power-iteration vector flip
    sign each pass?) that flags a negative dominant real eigenvalue as
    period-doubling -- exactly the peak-current-mode-without-slope-comp
    instability refs/README.md documents for R4 at D~0.55.

  - Output: the converged trajectory is emitted as one more DCtran() run
    (same restart mechanism), which produces an ordinary ngspice
    time-domain plot titled the same as the existing driven-mode output
    ("Time Domain Periodic Steady State Analysis"), so existing
    `wrdata`/`write`/`plot` tooling sees the same shape of result.

What this does NOT do (see the card #458 report for the full list):
  - No matrix-free/Krylov path (Sec.4.2) -- finite differences only, as
    the card allows for T4.1+T4.2.
  - No saltation/event handling and no free-running (autonomous) shooting
    -- both explicitly out of scope (design rulings #446 B / #447 A).
  - No brute-force ".tran-until-settled" fallback when Newton fails to
    converge or the monodromy matrix is singular -- this implementation
    reports FAILED/non-convergence honestly instead (never a steady state
    that was not reached), but does not then retry with plain periodic
    transient stepping. Flagged in the card report as remaining work.
  - The retained-trajectory substrate (PSSopTimes/Voltages/States,
    Enhancement-119) that the existing algorithm populates for PAC/pnoise
    to linearise around is not populated by the shooting path.
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/pssdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/fteext.h"

#include "cap/capdefs.h"
#include "ind/inddefs.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pssshoot.h"

/* ------------------------------------------------------------------ */
/* State-vector bookkeeping: one entry per native capacitor (terminal
 * voltage) and one per native inductor (branch current).              */

typedef struct {
    int isCap;
    union {
        CAPinstance *cap;
        INDinstance *ind;
    } dev;
    char name[64];
} PSSstateRef;

typedef struct {
    int n;
    PSSstateRef *refs;
} PSSstateMap;


static int
pss_build_state_map(CKTcircuit *ckt, PSSstateMap *map)
{
    int captype = CKTtypelook("Capacitor");
    int indtype = CKTtypelook("Inductor");
    int n = 0;
    CAPmodel *cm;
    CAPinstance *ci;
    INDmodel *im;
    INDinstance *ii;

    if (captype >= 0)
        for (cm = (CAPmodel *) ckt->CKThead[captype]; cm; cm = CAPnextModel(cm))
            for (ci = CAPinstances(cm); ci; ci = CAPnextInstance(ci))
                n++;
    if (indtype >= 0)
        for (im = (INDmodel *) ckt->CKThead[indtype]; im; im = INDnextModel(im))
            for (ii = INDinstances(im); ii; ii = INDnextInstance(ii))
                n++;

    map->n = n;
    map->refs = NULL;
    if (n == 0)
        return OK;

    map->refs = TMALLOC(PSSstateRef, n);
    if (!map->refs)
        return E_NOMEM;

    n = 0;
    if (captype >= 0)
        for (cm = (CAPmodel *) ckt->CKThead[captype]; cm; cm = CAPnextModel(cm))
            for (ci = CAPinstances(cm); ci; ci = CAPnextInstance(ci)) {
                map->refs[n].isCap = 1;
                map->refs[n].dev.cap = ci;
                snprintf(map->refs[n].name, sizeof map->refs[n].name,
                         "C:%s", ci->CAPname);
                n++;
            }
    if (indtype >= 0)
        for (im = (INDmodel *) ckt->CKThead[indtype]; im; im = INDnextModel(im))
            for (ii = INDinstances(im); ii; ii = INDnextInstance(ii)) {
                map->refs[n].isCap = 0;
                map->refs[n].dev.ind = ii;
                snprintf(map->refs[n].name, sizeof map->refs[n].name,
                         "L:%s", ii->INDname);
                n++;
            }
    return OK;
}


/* Read x from the device IC fields directly (no simulation has happened
 * yet -- used only to seed the very first Newton guess from whatever
 * IC=/.ic the deck already carries, default 0). */
static void
pss_state_get_ic(PSSstateMap *map, double *x)
{
    int k;
    for (k = 0; k < map->n; k++)
        x[k] = map->refs[k].isCap ? map->refs[k].dev.cap->CAPinitCond
                                   : map->refs[k].dev.ind->INDinitCond;
}

/* Write x into the device IC fields -- consumed by capload.c/indload.c
 * whenever CKTmode carries MODEUIC|MODEINITTRAN (the first small step of
 * a `uic` transient), which is how DCtran() below is driven. */
static void
pss_state_set_ic(PSSstateMap *map, const double *x)
{
    int k;
    for (k = 0; k < map->n; k++) {
        if (map->refs[k].isCap)
            map->refs[k].dev.cap->CAPinitCond = x[k];
        else
            map->refs[k].dev.ind->INDinitCond = x[k];
    }
}

/* Read x back from the converged solution (ckt->CKTrhsOld) after a period
 * run -- the same physical quantities pss_state_get_ic() seeds. */
static void
pss_state_get(PSSstateMap *map, CKTcircuit *ckt, double *x)
{
    int k;
    for (k = 0; k < map->n; k++) {
        if (map->refs[k].isCap) {
            CAPinstance *ci = map->refs[k].dev.cap;
            x[k] = ckt->CKTrhsOld[ci->CAPposNode] - ckt->CKTrhsOld[ci->CAPnegNode];
        } else {
            INDinstance *ii = map->refs[k].dev.ind;
            x[k] = ckt->CKTrhsOld[ii->INDbrEq];
        }
    }
}


/* ------------------------------------------------------------------ */
/* One-period map: run a full transient from t=0 to t=T, reusing
 * ngspice's own DCtran() (untouched). ckt->CKTcurJob is swapped to a
 * stack-local TRANan for the duration of the call -- TRANan and PSSan
 * share the JOBtype/JOBnextJob/JOBname header DCtran() actually reads
 * (plus TRANplot, which DCtran() owns exclusively during this call) --
 * and restored before returning. The caller must have already written
 * the desired x0 into the device IC fields via pss_state_set_ic(). */
static int
pss_run_period(CKTcircuit *ckt, double T, const char *jobname)
{
    JOB *saveJob = ckt->CKTcurJob;
    TRANan localjob;
    char namebuf[128];
    int err;

    memset(&localjob, 0, sizeof localjob);
    strncpy(namebuf, jobname, sizeof namebuf - 1);
    namebuf[sizeof namebuf - 1] = '\0';
    localjob.JOBtype = saveJob->JOBtype;
    localjob.JOBname = namebuf;

    ckt->CKTcurJob   = (JOB *) &localjob;
    ckt->CKTfinalTime = T;
    ckt->CKTinitTime  = 0.0;
    ckt->CKTtime      = 0.0;
    ckt->CKTmode     |= MODEUIC;   /* shooting always restarts from IC */

    err = DCtran(ckt, 1);   /* restart = 1: fresh run from t=0 */

    /* Defensive: some DCtran() error paths return without closing the
     * plot they opened; do not leak it across the (many) period runs one
     * Newton step needs. */
    if (localjob.TRANplot)
        SPfrontEnd->OUTendPlot(localjob.TRANplot);

    ckt->CKTcurJob = saveJob;
    return err;
}


/* ------------------------------------------------------------------ */
/* SUSPECTED CAUSE 2 (card #458 second pass): breakpoint-forcing support so
 * every period run within one Newton iteration retraces the SAME accepted
 * time grid -- see pssshoot.h for the design note and why this was chosen
 * over the fuller analytic S_k/J_k/C_k recursion. int, not bool, to match
 * the extern declarations dctran.c sees without pulling in <stdbool.h>. */
int PSSbreakRecordActive  = 0;
int PSSbreakInstallActive = 0;

static double *g_break_times = NULL;
static int g_break_count = 0, g_break_cap = 0;

static void
pss_break_reset(void)
{
    g_break_count = 0;
}

static void
pss_break_free(void)
{
    if (g_break_times) tfree(g_break_times);
    g_break_times = NULL;
    g_break_count = g_break_cap = 0;
}

/* Called from dctran.c right after CKTaccept() when PSSbreakRecordActive. */
void
PSSbreakRecordStep(CKTcircuit *ckt)
{
    if (ckt->CKTtime <= 0.0)
        return;   /* t=0 is already a breakpoint; nothing to add */

    /* Thin clustered accepted times before recording them. Found on the
     * R3 fixture: forcing EVERY accepted time verbatim (including a long
     * run right at/near CKTdelmin, where the baseline's own LTE control
     * was already stepping as finely as it could resolve a steep tanh
     * transition) made a perturbed run's corrector fail to converge
     * exactly at those forced points ("Timestep too small") -- the
     * perturbed trajectory's edge lands at a very slightly different time
     * than the baseline's, and forcing arrival at a long sequence of
     * sub-CKTdelmin-spaced points fights the corrector instead of
     * removing noise. Keeping only points at least maxStep*1e-3 apart
     * still constrains the grid everywhere it matters (this is still
     * ~1000x finer than the nominal per-period sampling) while leaving
     * the solver its own adaptive freedom within a cluster narrower than
     * that -- which is exactly the timescale too fine to matter for a
     * Newton monodromy estimate anyway. */
    if (g_break_count > 0 &&
        (ckt->CKTtime - g_break_times[g_break_count - 1]) < ckt->CKTmaxStep * 1e-3)
        return;

    if (g_break_count >= g_break_cap) {
        int newcap = g_break_cap ? g_break_cap * 2 : 256;
        double *nt = TREALLOC(double, g_break_times, newcap);
        if (!nt)
            return;   /* out of memory: degrade to unforced steps past here,
                          rather than crash -- the perturbed run may then
                          reintroduce some noise, but never a wrong answer */
        g_break_times = nt;
        g_break_cap = newcap;
    }
    g_break_times[g_break_count++] = ckt->CKTtime;
}

/* Called from dctran.c right after the restart-path breakpoint reset, when
 * PSSbreakInstallActive. Uses ngspice's own, unmodified CKTsetBreak(), so
 * every existing invariant it maintains (sorted list, CKTminBreak dedup,
 * "no breakpoint in the past") is preserved exactly. */
void
PSSbreakInstallHook(CKTcircuit *ckt)
{
    int i;
    for (i = 0; i < g_break_count; i++) {
        double tt = g_break_times[i];
        if (tt > 0.0 && tt < ckt->CKTfinalTime)
            CKTsetBreak(ckt, tt);
    }
}

/* Wraps pss_run_period(): record=1 captures this run's own accepted grid
 * (resetting whatever was recorded before); force=1 installs the
 * currently-recorded grid as breakpoints on this run instead. Never both:
 * a run either DEFINES the grid for this Newton iteration (the baseline
 * P(x0) evaluation) or FOLLOWS a previously-defined one (every perturbed
 * and trial run in the same iteration). */
static int
pss_run_period_ex(CKTcircuit *ckt, double T, const char *jobname,
                   int record, int force)
{
    int err;

    if (getenv("PSS_SHOOT_NO_FORCE") != NULL)
        force = 0;   /* debug/isolation switch, card #458 second pass */

    if (record)
        pss_break_reset();
    PSSbreakRecordActive  = record;
    PSSbreakInstallActive = force;

    err = pss_run_period(ckt, T, jobname);

    PSSbreakRecordActive  = 0;
    PSSbreakInstallActive = 0;
    return err;
}


/* ------------------------------------------------------------------ */
/* Small dense linear solve (Gauss elimination, partial pivoting) --
 * n is never more than a handful of states here, so this is simpler and
 * safer than routing a tiny dense system through the sparse solver.
 * A is n*n row-major and is destroyed; returns 0 on success, 1 if
 * singular (to working precision). */
static int
pss_solve_dense(double *A, double *b, int n, double *x)
{
    int i, j, k, piv;

    for (k = 0; k < n; k++) {
        double maxval = fabs(A[k * n + k]);
        piv = k;
        for (i = k + 1; i < n; i++)
            if (fabs(A[i * n + k]) > maxval) {
                maxval = fabs(A[i * n + k]);
                piv = i;
            }
        if (maxval < 1e-300)
            return 1;
        if (piv != k) {
            for (j = 0; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[piv * n + j];
                A[piv * n + j] = tmp;
            }
            { double tmp = b[k]; b[k] = b[piv]; b[piv] = tmp; }
        }
        for (i = k + 1; i < n; i++) {
            double f = A[i * n + k] / A[k * n + k];
            if (f == 0.0)
                continue;
            for (j = k; j < n; j++)
                A[i * n + j] -= f * A[k * n + j];
            b[i] -= f * b[k];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (j = i + 1; j < n; j++)
            s -= A[i * n + j] * x[j];
        x[i] = s / A[i * n + i];
    }
    return 0;
}


/* ------------------------------------------------------------------ */
/* card #458 FOURTH pass (plan T4.2): ANALYTIC (discretise-then-linearise)
 * sensitivity monodromy, replacing finite differences for the Jacobian.
 * Opt-in (PSS_SHOOT_JACOBIAN=analytic; default remains the third pass's FD
 * + Broyden path, byte-for-byte unchanged below). See the card #458
 * fourth-pass report for the full derivation; summary:
 *
 *   x0 (this file's shooting unknown) is exactly the set of native
 *   capacitor voltages and inductor branch currents (top-of-file design
 *   comment). At the UIC restart (t=0), capload.c/indload.c's own cond1
 *   bypass (MODEUIC && MODEINITTRAN) sets EACH cap's vcap / EACH
 *   inductor's branch current DIRECTLY from CAPinitCond/INDinitCond --
 *   i.e. from x0 itself, not from any solved node -- and their
 *   derivative-history state (CAPccap/INDvolt) is forced to exactly 0 at
 *   that same instant (niinteg.c's order-1 formula, state1[qcap] having
 *   just been set equal to state0[qcap] by the SAME cond1 branch, one
 *   line above). So d(each reactive device's own state)/dx0 = the unit
 *   row e_d EXACTLY at t=0, and d(ccap/INDvolt)/dx0 = 0 EXACTLY -- an
 *   exact base case, no approximation.
 *
 *   At every later accepted step k, capload.c (lines 74-79) and
 *   indload.c (lines 117, using NIintegrate's ceq/veq output) stamp the
 *   trapezoidal/backward-Euler companion history term (ceq for a cap,
 *   veq for an inductor) into CKTrhs with a FIXED sign pattern
 *   (CKTrhs[pos]-=ceq, CKTrhs[neg]+=ceq for a cap; CKTrhs[INDbrEq]+=veq
 *   for an inductor) that depends on x0 ONLY through the PREVIOUS step's
 *   values (q_{k-1}=C*v_{k-1}, and -- for trapezoidal order 2, once the
 *   integrator has stepped up from the first step's backward-Euler order
 *   1 -- the previous step's own CAPccap/INDvolt history state). Since
 *   the converged Newton fixed point at step k satisfies
 *   J(x_k)*x_k = RHS(x_k, history(x0)) IDENTICALLY (this is what Newton
 *   convergence means), and J_k := dG/dx_k at that fixed point is EXACTLY
 *   the matrix ngspice's own Newton iteration factored (implicit function
 *   theorem, using that every OTHER device here -- resistors, the B-source
 *   switches/comparator, the VCCS op-amp -- is memoryless, i.e. carries no
 *   explicit x0 dependence of its own):
 *
 *       J_k * S_k = d(RHS_k)/dx0 |_explicit
 *
 *   with S_k := dx_k/dx0 (an [nunk x n] matrix, n = the shooting state
 *   count) and the RHS built with EXACTLY the SAME stamp pattern as
 *   ceq/veq themselves (capload.c/indload.c), substituting the SCALAR
 *   ceq/veq for its OWN row-j sensitivity (also a scalar, one Newton
 *   shooting-direction j at a time):
 *
 *     order 1 (the first step of every period -- niinteg.c's order-1
 *              case does not reference the previous ccap/INDvolt at all):
 *       d(ceq_cap)/dx0_j  = -ag0 * C * (S_k-1[pos,j] - S_k-1[neg,j])
 *       d(veq_ind)/dx0_j  = -ag0 * L * S_k-1[INDbrEq,j]
 *     order 2 (trapezoidal, niinteg.c's order-2 case):
 *       d(ceq_cap)/dx0_j  = -ag1*T_k-1[cap,j] - ag0*C*(S_k-1[pos,j]-S_k-1[neg,j])
 *       d(veq_ind)/dx0_j  = -ag1*T_k-1[ind,j] - ag0*L* S_k-1[INDbrEq,j]
 *
 *   (ag0=CKTag[0], ag1=CKTag[1], both read at the SAME accepted step whose
 *   J_k is used -- see nicomcof.c for their definition per order/method.)
 *   Solving J_k*S_k[:,j] = b_k[:,j] for j=0..n-1 (n back-substitutions,
 *   REUSING the one factorization -- SMPsolve(), no refactor per column)
 *   gives S_k. T_k (the derivative-history state's OWN sensitivity) then
 *   updates by the SAME recursion niinteg.c itself uses for ccap/INDvolt,
 *   differentiated:
 *       order 1: T_k = ag0*C_or_L*(Srow_k - Srow_k-1)
 *       order 2: T_k = -ag1*T_k-1 + ag0*C_or_L*(Srow_k - Srow_k-1)
 *   (Srow = the cap's v_pos-v_neg row, or the inductor's branch-current
 *   row.) Propagated for every accepted step of ONE period run (a hook at
 *   the exact call site pssshoot.h's breakpoint hooks already use, right
 *   after CKTaccept()), S_K's rows for the n reactive devices ARE the
 *   monodromy matrix Phi -- no finite-difference perturbation runs, no
 *   eps bias, one period integration instead of 2n+1.
 *
 *   J_k itself: the LAST factorization NIiter() leaves in ckt->CKTmatrix
 *   after a step converges was computed from the SECOND-TO-LAST Newton
 *   iterate (the load that PRODUCED the converged x_k via one more
 *   solve), not from x_k itself -- close but not exact for a nonlinear
 *   circuit (the card's own "NIiter may bypass refactoring" caveat).
 *   PSSsensRecordStep() below forces one extra CKTload()+SMPluFac() at
 *   the now-fixed, converged ckt->CKTrhsOld before using CKTmatrix, per
 *   the card's explicit instruction -- cheap (no NR iteration; the real
 *   circuit unknowns are untouched) and makes J_k exact at x_k.
 *
 *   Scope: TRAPEZOIDAL orders 1-2 only (this deck family's default
 *   integration method, and the only orders niinteg.c's TRAPEZOIDAL case
 *   has). GEAR (.options method=gear) is NOT implemented -- g_sens_failed
 *   is set with a clear message if CKTorder/CKTintegrateMethod fall
 *   outside {TRAPEZOIDAL, order<=2} rather than silently producing a wrong
 *   answer; the caller (pss_newton_run) then reports this attempt as a
 *   plain non-convergence, the same honest "did not confirm a steady
 *   state" signal a failed FD sweep already produces. */

int PSSsensActive = 0;

static PSSstateMap *g_sens_map = NULL;
static int g_sens_n = 0, g_sens_neq = 0, g_sens_failed = 1;
static double *g_sens_Sprev = NULL, *g_sens_Snew = NULL;   /* [(neq+1)*n] */
static double *g_sens_Tprev = NULL, *g_sens_Tnew = NULL;   /* [n*n] */
static double *g_sens_rhs = NULL, *g_sens_sol = NULL;      /* [neq+1] scratch */

static double
pss_sens_get(const double *S, int node, int n, int j)
{
    if (node <= 0 || !S) return 0.0;
    return S[(size_t) node * n + j];
}

static void
pss_sens_set(double *S, int node, int n, int j, double v)
{
    if (node <= 0 || !S) return;
    S[(size_t) node * n + j] = v;
}

static void
pss_sens_free(void)
{
    if (g_sens_Sprev) tfree(g_sens_Sprev);
    if (g_sens_Snew)  tfree(g_sens_Snew);
    if (g_sens_Tprev) tfree(g_sens_Tprev);
    if (g_sens_Tnew)  tfree(g_sens_Tnew);
    if (g_sens_rhs)   tfree(g_sens_rhs);
    if (g_sens_sol)   tfree(g_sens_sol);
    g_sens_Sprev = g_sens_Snew = g_sens_Tprev = g_sens_Tnew = NULL;
    g_sens_rhs = g_sens_sol = NULL;
    g_sens_map = NULL;
    g_sens_n = g_sens_neq = 0;
    g_sens_failed = 1;
}

/* Called by pssshoot.c itself (not dctran.c) right before a period run
 * that should have PSSsensActive=1 -- sets the exact t=0 base case
 * derived above and (re)allocates this run's scratch. map must outlive
 * the run (it is the caller's own, stack-owned-by-the-Newton-attempt
 * PSSstateMap). */
static void
pss_sens_reset(CKTcircuit *ckt, PSSstateMap *map)
{
    int n = map->n, neq, d;

    pss_sens_free();

    neq = SMPmatSize(ckt->CKTmatrix);
    g_sens_map = map;
    g_sens_n = n;
    g_sens_neq = neq;

    g_sens_Sprev = TMALLOC(double, (size_t) (neq + 1) * n);
    g_sens_Snew  = TMALLOC(double, (size_t) (neq + 1) * n);
    g_sens_Tprev = TMALLOC(double, (size_t) n * n);
    g_sens_Tnew  = TMALLOC(double, (size_t) n * n);
    g_sens_rhs   = TMALLOC(double, (size_t) (neq + 1));
    g_sens_sol   = TMALLOC(double, (size_t) (neq + 1));
    if (!g_sens_Sprev || !g_sens_Snew || !g_sens_Tprev || !g_sens_Tnew ||
        !g_sens_rhs || !g_sens_sol) {
        pss_sens_free();
        return;
    }

    memset(g_sens_Sprev, 0, (size_t) (neq + 1) * n * sizeof(double));
    memset(g_sens_Tprev, 0, (size_t) n * n * sizeof(double));

    for (d = 0; d < n; d++) {
        PSSstateRef *ref = &map->refs[d];
        int row = ref->isCap ? ref->dev.cap->CAPposNode : ref->dev.ind->INDbrEq;
        pss_sens_set(g_sens_Sprev, row, n, d, 1.0);
    }
    g_sens_failed = 0;
}

/* dctran.c hook, called right after CKTaccept() (same site as
 * PSSbreakRecordStep) whenever PSSsensActive. See the long design comment
 * above for the recursion this implements. */
void
PSSsensRecordStep(CKTcircuit *ckt)
{
    int n, neq, d, j, nd, order;
    double ag0, ag1;

    if (g_sens_failed || !g_sens_map)
        return;
    if (ckt->CKTtime <= 0.0)
        return;   /* t=0 base case already set by pss_sens_reset() */

    if (ckt->CKTintegrateMethod != TRAPEZOIDAL || ckt->CKTorder < 1 ||
        ckt->CKTorder > 2) {
        fprintf(stderr, "PSS-SENS: analytic sensitivity only implements "
                "TRAPEZOIDAL orders 1-2 (got method=%d order=%d at "
                "t=%g); abandoning the analytic Jacobian for this period "
                "run.\n", ckt->CKTintegrateMethod, ckt->CKTorder, ckt->CKTtime);
        g_sens_failed = 1;
        return;
    }

    n   = g_sens_n;
    neq = g_sens_neq;
    ag0 = ckt->CKTag[0];
    ag1 = ckt->CKTag[1];

    /* card #458 fourth pass, found by this pass's own R2 validation (the
     * exploding-Phi bug): ckt->CKTorder, read here right after
     * CKTaccept(), is NOT reliably the order whose niinteg.c FORMULA
     * actually produced the just-accepted step's CAPccap/INDvolt values --
     * dctran.c's own order/step-size selection (ckttrunc.c, called before
     * CKTaccept() in dctran.c's control flow) already advances CKTorder to
     * what the NEXT step will attempt, while ckt->CKTag[] here still holds
     * the coefficients NIcomCof() computed for the step that was just
     * accepted (it is not recomputed until the next NIiter() call). The
     * two can disagree (observed on R2: CKTorder reports 2 on the very
     * first trapezoidal step after the order-1 opening step, while CKTag[]
     * still holds order 1's coefficients) -- using CKTorder there fed the
     * order-2 formula's -ag1*T_{k-1} term (ag1 up to ~1/h, huge) into a
     * step whose ACTUAL ag1 was order 1's -ag0 (no such term exists),
     * exploding S_k by many orders of magnitude within a handful of
     * steps. nicomcof.c's own two formulas are structurally
     * unmistakable: order 1 gives ag1 = -ag0 EXACTLY (ckt->CKTdelta
     * cancels); order 2 gives ag1 = CKTxmu/(1-CKTxmu), an O(1) constant
     * (0.5/0.5=1 for the default xmu) with NO relation to ag0's own scale
     * (~1/h, which varies over many orders of magnitude step to step). So
     * the order actually used is inferred from ag1 itself, not trusted
     * from ckt->CKTorder. */
    order = (fabs(ag1 + ag0) < 1e-6 * fabs(ag0)) ? 1 : 2;

    /* NOT refactored here -- card #458 fourth-pass finding (both tried;
     * neither is clean, this one was kept as the shipped default). The
     * design comment above and the card's own text say to force a fresh
     * CKTload()+refactor at the accepted point, since NIiter()'s last
     * factorization was computed from the SECOND-TO-LAST Newton iterate's
     * load, not x_k itself. That WAS tried here and found to be actively
     * WRONG on two counts:
     *  (a) R2's `SW` switch model (swload.c -- current_state/
     *      old_current_state read from CKTstate0/CKTstate1 depend on CALL
     *      SEQUENCE within a timestep, not purely on (x_k, history), for
     *      a hysteretic device): a second, out-of-band CKTload() call
     *      after CKTaccept() measurably CHANGED the physical trajectory
     *      (confirmed: the same x0, run once with PSSsensActive so this
     *      extra load fires every step, gave a visibly different P(x0)
     *      than a plain period run from the same x0);
     *  (b) on R3 (7 states, ~200-2000+ accepted steps/period), refactoring
     *      at EVERY accepted step (CKTload() + a full SMPluFac()) was
     *      found to make even the small stabtime=0/sc_iter=1 single-
     *      Jacobian diagnostic deck (tests/fixtures/pss_r3_jaccheck.cir)
     *      run for minutes rather than the sub-second baseline -- not
     *      profiled further given the time available this pass, but
     *      clearly disqualifying for an option meant to be FASTER than
     *      FD, not slower.
     * Skipping the refactor avoids both, at the cost of using a J_k that
     * is the second-to-last Newton iterate's, not exactly x_k's -- found,
     * on R3's element-wise FD-vs-analytic comparison (card #458 fourth-
     * pass report), to leave most Phi entries matching FD to 4-5
     * significant figures, but the columns for states electrically close
     * to the comparator's high-gain tanh nonlinearity (opout, and C1f
     * which couples into it) wrong by 1-2 orders of magnitude -- the
     * stale J_k's small error in that device's steep incremental
     * conductance is amplified enough by the kc=30 comparator gain to
     * matter a lot there. This is why R3 Newton does not reliably converge
     * with PSS_SHOOT_JACOBIAN=analytic this pass (reported honestly as a
     * known gap, not silently accepted) while R2 (1 state, no such
     * amplifying path) converges cleanly. Neither refactor option is
     * clean for every fixture; this is the one shipped as the default
     * (and only) analytic path, since hanging is worse than an inaccurate
     * Jacobian that Newton then honestly fails to converge with. */

    memset(g_sens_Snew, 0, (size_t) (neq + 1) * n * sizeof(double));

    for (j = 0; j < n; j++) {
        memset(g_sens_rhs, 0, (size_t) (neq + 1) * sizeof(double));

        for (d = 0; d < n; d++) {
            PSSstateRef *ref = &g_sens_map->refs[d];
            double Tp = g_sens_Tprev[d * n + j];
            if (ref->isCap) {
                CAPinstance *ci = ref->dev.cap;
                int p = ci->CAPposNode, ng = ci->CAPnegNode;
                double C = ci->CAPcapac;
                double Sp = pss_sens_get(g_sens_Sprev, p, n, j);
                double Sn = pss_sens_get(g_sens_Sprev, ng, n, j);
                double dceq = -ag0 * C * (Sp - Sn);
                if (order == 2) dceq += -ag1 * Tp;
                if (p  > 0) g_sens_rhs[p]  -= dceq;
                if (ng > 0) g_sens_rhs[ng] += dceq;
            } else {
                INDinstance *ii = ref->dev.ind;
                int ibr = ii->INDbrEq;
                double L = ii->INDinduct;
                double Si = pss_sens_get(g_sens_Sprev, ibr, n, j);
                double dveq = -ag0 * L * Si;
                if (order == 2) dveq += -ag1 * Tp;
                g_sens_rhs[ibr] += dveq;
            }
        }

        /* SMPsolve()'s third argument ("Spare"/"Solution") is a dead
         * legacy parameter -- spsmp.c's own SMPsolve() does
         * `NG_IGNORE(Spare); spSolve(Matrix->SPmatrix, RHS, RHS, NULL,
         * NULL);`, i.e. it solves IN PLACE and overwrites g_sens_rhs
         * itself with the solution (exactly how niiter.c's own
         * `SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare)` call
         * works -- NIconvTest() and the node-damping code right after it
         * both read the solution back from ckt->CKTrhs, never
         * ckt->CKTrhsSpare). Found by this pass's own R2 validation
         * (g_sens_sol was reading back all-zero, uninitialised scratch;
         * see the card #458 fourth-pass report). g_sens_sol is kept only
         * as the (unused) Spare argument SMPsolve's signature requires. */
        SMPsolve(ckt->CKTmatrix, g_sens_rhs, g_sens_sol);
        for (nd = 1; nd <= neq; nd++)
            g_sens_Snew[(size_t) nd * n + j] = g_sens_rhs[nd];
    }

    for (d = 0; d < n; d++) {
        PSSstateRef *ref = &g_sens_map->refs[d];
        for (j = 0; j < n; j++) {
            double val;
            if (ref->isCap) {
                CAPinstance *ci = ref->dev.cap;
                int p = ci->CAPposNode, ng = ci->CAPnegNode;
                double C = ci->CAPcapac;
                double Sp_new = pss_sens_get(g_sens_Snew, p, n, j);
                double Sn_new = pss_sens_get(g_sens_Snew, ng, n, j);
                double Sp_old = pss_sens_get(g_sens_Sprev, p, n, j);
                double Sn_old = pss_sens_get(g_sens_Sprev, ng, n, j);
                val = ag0 * C * ((Sp_new - Sn_new) - (Sp_old - Sn_old));
            } else {
                INDinstance *ii = ref->dev.ind;
                int ibr = ii->INDbrEq;
                double L = ii->INDinduct;
                double Si_new = pss_sens_get(g_sens_Snew, ibr, n, j);
                double Si_old = pss_sens_get(g_sens_Sprev, ibr, n, j);
                val = ag0 * L * (Si_new - Si_old);
            }
            if (order == 2)
                val += -ag1 * g_sens_Tprev[d * n + j];
            g_sens_Tnew[d * n + j] = val;
        }
    }

    { double *tmp = g_sens_Sprev; g_sens_Sprev = g_sens_Snew; g_sens_Snew = tmp; }
    { double *tmp = g_sens_Tprev; g_sens_Tprev = g_sens_Tnew; g_sens_Tnew = tmp; }

    if (getenv("PSS_SENS_TRACE")) {
        static int calls = 0;
        if (calls < 40)
            fprintf(stderr, "SENSTRACE call=%d t=%.6g h=%.4g order=%d ag0=%.6g "
                    "S(owned,0)=%.10g\n", calls, ckt->CKTtime, ckt->CKTdelta, order, ag0,
                    pss_sens_get(g_sens_Sprev,
                        g_sens_map->refs[0].isCap ? g_sens_map->refs[0].dev.cap->CAPposNode
                                                   : g_sens_map->refs[0].dev.ind->INDbrEq,
                        n, 0));
        calls++;
    }
}

/* Reads the current (t=T, after a full period's worth of
 * PSSsensRecordStep() calls) monodromy matrix Phi[n*n] (row-major, RAW,
 * not Phi-I -- same convention pss_newton_run()'s FD path uses). Returns
 * 0 on success, 1 if the analytic path failed anywhere this run (a
 * refactor/solve error, or an order/method outside TRAPEZOIDAL 1-2). */
static int
pss_sens_get_phi(PSSstateMap *map, double *Phi)
{
    int n = map->n, d, j;

    if (g_sens_failed || !g_sens_Sprev)
        return 1;

    for (d = 0; d < n; d++) {
        PSSstateRef *ref = &map->refs[d];
        for (j = 0; j < n; j++) {
            if (ref->isCap) {
                int p = ref->dev.cap->CAPposNode, ng = ref->dev.cap->CAPnegNode;
                Phi[d * n + j] = pss_sens_get(g_sens_Sprev, p, n, j) -
                                  pss_sens_get(g_sens_Sprev, ng, n, j);
            } else {
                Phi[d * n + j] = pss_sens_get(g_sens_Sprev, ref->dev.ind->INDbrEq, n, j);
            }
        }
    }
    return 0;
}


/* SUSPECTED CAUSE 1 fix: the former spectral-radius estimator (power
 * iteration) is gone -- it cannot resolve a complex-conjugate dominant
 * multiplier pair, which is exactly what a slow oscillation (the R3 x20
 * red-proof) produces, and that was the confirmed bug (ground truth,
 * parent, 2026-09-28: nominal R3 got an inconsistent 0.99-2.4, and the
 * genuinely-unstable x20-gain variant got 0.98 "stable"). Floquet
 * classification is now pss_eigenvalues() + pss_classify_floquet() from
 * psseigen.c (full real Hessenberg + double-shift QR eigendecomposition,
 * validated against an independent root-finder -- see psseigen.c's header
 * and the card #458 report). */


/* THIRD PASS (card #458, plan T4.1): one Newton attempt from a given
 * starting point x0, with an iteration budget of its own (maxiter_this) --
 * used both for a single full-budget attempt (when the deck's whole
 * stabtime already ran, trust_radius<=0) and for the cheap, short-budget
 * "probe" attempts the staged warm-up escalation in PSSshootDriven tries
 * first (see that loop's own header comment for the warm-up side of this
 * pass's speed work). This function is the OTHER half: two independent
 * cuts to the ~2n+3-period-run cost of each Newton ITERATION, on top of
 * whatever warm-up got it here:
 *
 *   - the second pass's "post-monodromy baseline re-run" (1 extra period
 *     run every FD iteration) turns out to be dead weight: nothing
 *     downstream of it actually needs a fresh x1 except a DEBUG print --
 *     F[]/scale[] are already correct from the FIRST baseline call, before
 *     the FD loop's own pss_state_get() calls overwrite x1 as scratch.
 *     Fixed by giving the FD loop its own scratch buffer (`pert`) instead
 *     of reusing x1, so the initial baseline survives untouched and that
 *     re-run is simply not needed -- same F/Phi/x0 dataflow, one fewer
 *     period run per FD Jacobian.
 *   - Broyden rank-1 updates of the RAW monodromy estimate (not of
 *     Phi-I -- the -I term is exact, never approximated) between full
 *     finite-difference recomputes, the task's "or ... Broyden updates"
 *     option, chosen over one-sided FD: one-sided differences were already
 *     tried on this circuit in the second pass and gave no usable descent
 *     direction (see the centred-difference comment below, kept
 *     verbatim), so trying that again risked repeating a known failure.
 *     Broyden carries no equivalent risk because it reuses function values
 *     the line search was ALREADY computing for its own residual check --
 *     no new perturbation grid, no new eps, nothing that could shift a
 *     switch edge. After the first FD Jacobian in an attempt, an accepted
 *     damped step already gives P(xtrial) for free -- that IS next
 *     iteration's baseline P(x0), so a Broyden-updated iteration costs
 *     only its line-search trial run(s) (typically 1, since a good
 *     direction usually accepts alpha=1 immediately), not 2n+3. The
 *     Jacobian is refreshed with a full FD sweep every broyden_period
 *     iterations (default 4) as a safety net, and immediately, mid-
 *     attempt, if a Broyden-based line search fails to find ANY decreasing
 *     step or the Broyden matrix goes singular -- only a FRESH FD
 *     Jacobian's own failure ends the attempt, so Broyden can only make a
 *     converging case converge faster, never turn a converging case into
 *     a failing one. The exact converged point is still checked against
 *     the TRUE nonlinear residual (real P(x0) evaluations, no
 *     approximation) either way -- Broyden only ever changes the SEARCH
 *     DIRECTION, never what counts as converged.
 *   - a trust-region cap on the first line-search trial's alpha (not on
 *     the Newton direction itself, which is left as the solve produces
 *     it): when trust_radius>0, the initial trial step is capped so no
 *     state moves more than trust_radius*scale before backtracking is
 *     even tried, then ordinary backtracking proceeds from there exactly
 *     as before. A rougher (short-warm-up) starting point can give a
 *     poorly-conditioned FD Jacobian, so this keeps a bad direction from
 *     being taken at full stride before backtracking gets a chance to
 *     correct it. trust_radius<=0 disables it (unrestricted alpha=1
 *     start) -- used for the full-budget, full-stabtime attempt, which is
 *     then byte-for-byte the second pass's own tested behaviour, kept as
 *     the always-available fallback.
 *
 * Returns 1 = converged (x0/iter_out/residual_out hold the converged point;
 * Phi_out, if non-NULL, holds the RAW monodromy estimate at that point --
 * same quantity the pre-third-pass code called Phi, NOT Phi-I -- for the
 * caller's Floquet-fallback path). Returns 0 for EVERY other outcome
 * (budget exhausted, singular Jacobian, line search failure, or a period
 * integration that did not complete) -- deliberately not distinguished
 * further, so the caller's staged warm-up escalation has one uniform
 * "this attempt did not confirm a steady state" signal to react to; x0 is
 * left at its last ACCEPTED iterate (never a rejected/failed trial), and
 * the specific cause was already printed to stderr when it happened. */
static int
pss_newton_run(CKTcircuit *ckt, PSSstateMap *map, double T, int n,
                double *x0, int maxiter_this, double tol, double eps,
                double trust_radius, int broyden_enabled, int broyden_period,
                const char *stage_label,
                int *iter_out, double *residual_out, double *Phi_out,
                int use_analytic)
{
    int iter = 0, j, k, err;
    double residual = 0.0;
    double *x1 = NULL, *F = NULL, *scale = NULL, *delta_x = NULL,
           *xtrial = NULL, *x1trial = NULL, *A = NULL, *b = NULL,
           *M = NULL, *pert = NULL, *s = NULL, *y = NULL, *Ms = NULL;
    int have_baseline = 0;   /* do x0/x1/F already hold a valid, just-run point? */
    int recompute = 1;       /* does THIS iteration's Jacobian need a fresh FD sweep? */
    int since_recompute = 0;
    int result = 0;

    /* card #458 FOURTH pass: the analytic Jacobian has no Broyden analogue
     * in this implementation -- every iteration re-runs ONE period WITH
     * sensitivity tracking active (see below), which is both this
     * iteration's baseline P(x0) AND its Jacobian source, so "reuse last
     * iteration's free P(xtrial)" (what have_baseline=1 normally means at
     * the bottom of the loop) would skip recomputing Phi. Forcing
     * broyden_enabled=0 here makes `recompute` stay permanently true
     * (its own end-of-iteration update, unchanged below, already reduces
     * to that whenever broyden_enabled is false) with no other code path
     * changes needed for that part. */
    if (use_analytic)
        broyden_enabled = 0;

    x1      = TMALLOC(double, n);
    F       = TMALLOC(double, n);
    scale   = TMALLOC(double, n);
    delta_x = TMALLOC(double, n);
    xtrial  = TMALLOC(double, n);
    x1trial = TMALLOC(double, n);
    A       = TMALLOC(double, n * n);
    b       = TMALLOC(double, n);
    M       = TMALLOC(double, n * n);
    pert    = TMALLOC(double, n);
    s       = TMALLOC(double, n);
    y       = TMALLOC(double, n);
    Ms      = TMALLOC(double, n);
    if (!x1 || !F || !scale || !delta_x || !xtrial || !x1trial || !A || !b ||
        !M || !pert || !s || !y || !Ms) {
        fprintf(stderr, "PSS-SHOOT: %s out of memory allocating Newton "
                "scratch space.\n", stage_label);
        result = 0;
        goto newton_done;
    }

    for (iter = 0; ; iter++) {

        if (!have_baseline) {
            pss_state_set_ic(map, x0);
            if (use_analytic) {
                /* One period run, sensitivity tracking active: this SAME
                 * run is both the baseline P(x0) (x1 below) AND the
                 * analytic Jacobian source (M below) -- no separate FD
                 * sweep, no separate baseline call. Not forced onto any
                 * grid (force=0): the analytic recursion uses THIS run's
                 * own accepted steps directly (via PSSsensRecordStep(),
                 * hooked at the same dctran.c call site as the
                 * breakpoint-recording hook), so there is no cross-run
                 * grid-consistency need the FD path's forcing addressed. */
                pss_sens_reset(ckt, map);
                PSSsensActive = 1;
                err = pss_run_period_ex(ckt, T, "PSS shoot (Newton baseline period, analytic)", 0, 0);
                PSSsensActive = 0;
            } else {
                err = pss_run_period_ex(ckt, T, "PSS shoot (Newton baseline period)", 1, 0);
            }
            if (err) {
                fprintf(stderr, "PSS-SHOOT: %s period integration did not "
                        "complete at Newton iteration %d (ngspice error "
                        "%d); this is NOT a steady state.\n", stage_label, iter, err);
                result = 0;
                goto newton_done;
            }
            pss_state_get(map, ckt, x1);
            if (use_analytic && getenv("PSS_SHOOT_DEBUG")) {
                double dbgPhi[64];
                if (!pss_sens_get_phi(map, dbgPhi)) {
                    for (k = 0; k < n * n; k++)
                        fprintf(stderr, "DEBUG analytic Phi[%d]=%.8g\n", k, dbgPhi[k]);
                }
            }
            if (use_analytic && pss_sens_get_phi(map, M)) {
                fprintf(stderr, "PSS-SHOOT: %s analytic sensitivity failed "
                        "(refactor/solve error, or an integration order "
                        "outside TRAPEZOIDAL 1-2) at Newton iteration %d; "
                        "treating as a failed attempt (no FD fallback "
                        "within one attempt).\n", stage_label, iter);
                result = 0;
                goto newton_done;
            }

            /* Scale/floor: unchanged from the second pass -- see that
             * pass's comment (kept verbatim in PSSshootDriven's history;
             * not repeated here) for why a 1 mV floor, not near-zero. */
            residual = 0.0;
            for (k = 0; k < n; k++) {
                double rng = fabs(x1[k]) > fabs(x0[k]) ? fabs(x1[k]) : fabs(x0[k]);
                const double ABSTOL_FLOOR = 1e-3;
                scale[k] = (rng > ABSTOL_FLOOR) ? rng : ABSTOL_FLOOR;
                F[k] = x1[k] - x0[k];
                {
                    double r = fabs(F[k]) / scale[k];
                    if (r > residual) residual = r;
                }
            }
            have_baseline = 1;
        }

        fprintf(stdout, "PSS-SHOOT: %s iter=%d residual=%.6e%s\n", stage_label,
                iter, residual, recompute ? "" : " (broyden)");
        fflush(stdout);
        if (getenv("PSS_SHOOT_DEBUG")) {
            for (k = 0; k < n; k++)
                fprintf(stderr, "DEBUG %s iter=%d k=%d %s x0=%.8g x1=%.8g "
                        "F=%.4e scale=%.6g r=%.4e\n", stage_label,
                        iter, k, map->refs[k].name, x0[k], x1[k], F[k],
                        scale[k], fabs(F[k]) / scale[k]);
        }

        if (residual <= tol) {
            /* THIRD-PASS GUARD, added after this pass's OWN testing found
             * it necessary: reject a "converged" point where EVERY state
             * is still at (or under) the ABSTOL_FLOOR scale floor -- the
             * degenerate all-off fixed point this file's top-of-function
             * comment already warns about (a converter sitting powered
             * off is itself a self-consistent, non-switching P(0)=0
             * solution). The original design relied ENTIRELY on the
             * pre-Newton warm-up already being far from 0 to avoid this;
             * that is no longer a safe assumption with a SHORT probe
             * warm-up and a generous iteration budget -- found exactly
             * this way on the R3 fixture: PSS_SHOOT_WARMUP_MAXITER=10
             * from a 20-period probe converged (residual 1.1e-10, a
             * seemingly excellent fit) onto vout~1.2e-7 V, not ~3.3V, a
             * numerically real but physically wrong answer (and a
             * telltale near-singular Phi~=I / Floquet dominant
             * multiplier exactly 1+0i went with it). A real switching
             * converter's TRUE orbit should not have every single state
             * sitting at/under a millivolt/milliamp-scale floor at once,
             * so this is a general, circuit-agnostic check, not a magic
             * number tuned to R3 -- and it is applied at every stage,
             * including a full-budget/full-stabtime one, since nothing
             * before this pass actually proved the old code was immune
             * to the same trap, only that it never happened to trigger
             * it. A rejection here is treated exactly like any other
             * "this attempt did not confirm a steady state" outcome:
             * the staged caller extends the warm-up (or, on the final
             * stage, reports the honest FAILED verdict) rather than
             * silently accepting a wrong answer. */
            int all_at_floor = 1;
            for (k = 0; k < n; k++)
                if (scale[k] > 1e-3 * (1.0 + 1e-9)) { all_at_floor = 0; break; }
            if (all_at_floor) {
                fprintf(stderr, "PSS-SHOOT: %s converged (residual=%.3e) "
                        "but EVERY state is at/under the 1e-3 scale floor "
                        "-- rejecting this as the degenerate all-off fixed "
                        "point, not a trustworthy orbit.\n", stage_label,
                        residual);
                result = 0;
                goto newton_done;
            }
            result = 1;
            goto newton_done;
        }
        if (iter >= maxiter_this) {
            result = 0;
            goto newton_done;   /* out of this attempt's budget */
        }

        if (recompute && !use_analytic) {
            /* finite-difference monodromy: 2n extra period runs (centred
             * differences, not the cheaper n forward differences the
             * card's "n+1 period runs per Newton step" description allows
             * -- found on the R3 fixture, second pass: one-sided
             * differences gave a Newton direction whose backtrack line
             * search approached the baseline residual from above for 24
             * straight halvings without ever dipping below it, i.e.
             * essentially no usable descent direction; centred
             * differences cancel the O(h) bias one-sided differences
             * carry and fixed it. `pert`, not x1, is the scratch buffer
             * here (see this function's header comment) so the clean
             * baseline x1/F/scale survive for the residual check above on
             * every iteration, recompute or not. */
            for (j = 0; j < n; j++) {
                double h = fabs(x0[j]) * eps;
                if (h < 1e-6) h = 1e-6;

                for (k = 0; k < n; k++) xtrial[k] = x0[k];
                xtrial[j] = x0[j] + h;
                pss_state_set_ic(map, xtrial);
                err = pss_run_period_ex(ckt, T, "PSS shoot (monodromy perturbation +)", 0, 1);
                if (err) {
                    fprintf(stderr, "PSS-SHOOT: %s monodromy perturbation "
                            "run for state %s (+h) did not complete at "
                            "Newton iteration %d (ngspice error %d).\n",
                            stage_label, map->refs[j].name, iter, err);
                    result = 0;
                    goto newton_done;
                }
                pss_state_get(map, ckt, pert);
                for (k = 0; k < n; k++) M[k * n + j] = pert[k];

                for (k = 0; k < n; k++) xtrial[k] = x0[k];
                xtrial[j] = x0[j] - h;
                pss_state_set_ic(map, xtrial);
                err = pss_run_period_ex(ckt, T, "PSS shoot (monodromy perturbation -)", 0, 1);
                if (err) {
                    fprintf(stderr, "PSS-SHOOT: %s monodromy perturbation "
                            "run for state %s (-h) did not complete at "
                            "Newton iteration %d (ngspice error %d).\n",
                            stage_label, map->refs[j].name, iter, err);
                    result = 0;
                    goto newton_done;
                }
                pss_state_get(map, ckt, pert);
                for (k = 0; k < n; k++)
                    M[k * n + j] = (M[k * n + j] - pert[k]) / (2.0 * h);
            }
            if (getenv("PSS_SHOOT_DEBUG"))
                for (k = 0; k < n * n; k++)
                    fprintf(stderr, "DEBUG FD Phi[%d]=%.8g\n", k, M[k]);
            since_recompute = 0;
            /* M is the RAW monodromy estimate here (no -I baked in) --
             * kept that way throughout (Broyden update below preserves
             * it) so Phi_out is always a true Phi, and the -I is
             * subtracted only into the disposable solve copy A just
             * below, same as the pre-third-pass code did. */
        }

        for (k = 0; k < n * n; k++) A[k] = M[k];
        for (k = 0; k < n; k++) A[k * n + k] -= 1.0;
        for (k = 0; k < n; k++) b[k] = -F[k];

        if (pss_solve_dense(A, b, n, delta_x)) {
            if (!recompute) {
                /* a stale Broyden estimate went singular -- not this
                 * attempt's honest answer yet: force a fresh FD sweep and
                 * retry from the SAME x0 (have_baseline is still valid,
                 * so this costs only the 2n FD runs, not another
                 * baseline) before giving up. */
                fprintf(stdout, "PSS-SHOOT: %s iter=%d Broyden estimate "
                        "went singular; forcing a fresh Jacobian.\n",
                        stage_label, iter);
                recompute = 1;
                iter--;   /* this attempt does not count against the budget */
                continue;
            }
            fprintf(stderr, "PSS-SHOOT: %s the monodromy matrix (Phi - I) "
                    "is singular at Newton iteration %d; Newton cannot "
                    "proceed this attempt.\n", stage_label, iter);
            result = 0;
            goto newton_done;
        }

        {
            double maxratio = 0.0, alpha;
            int bt, accepted = 0;
            /* THIRD-PASS GUARD (added after this pass's own testing found
             * it necessary, alongside the near-zero-floor guard above):
             * an absolute sanity ceiling on trial state magnitudes, NOT
             * scale-relative. The scale-relative trust-region cap above
             * failed to prevent a real divergence found on the R3 fixture
             * (PSS_SHOOT_WARMUP_MAXITER=10 from a 20-period probe): a
             * compensator integrator state (C:coa/C:c1f/C:c2f) railed to
             * ~1e5 V while the power-stage states (vout, i(L1), ...) sat
             * near 0, and Newton found an exact (residual~1e-10)
             * mathematical fixed point of THAT combination -- because
             * scale[] is itself derived from the current (already-
             * diverging) x0/x1, a trust region measured relative to scale
             * cannot catch a state that is diverging together with its
             * own yardstick. An absolute ceiling has no such blind spot:
             * no state in any of this project's reference circuits
             * (volts, amps, duty-cycle-scaled internal nodes) belongs
             * anywhere near 1e4 of its own native unit, so a trial that
             * crosses it is rejected outright -- treated exactly like a
             * failed integration (skip the period run, just halve alpha
             * again), never accepted, and never even charged a period
             * run. PSS_SHOOT_MAX_STATE overrides (debug env var, no
             * effect unless set, same pattern as the others). */
            double max_state = getenv("PSS_SHOOT_MAX_STATE") ?
                                atof(getenv("PSS_SHOOT_MAX_STATE")) : 1e4;

            for (k = 0; k < n; k++) {
                double r = fabs(delta_x[k]) / scale[k];
                if (r > maxratio) maxratio = r;
            }
            alpha = 1.0;
            if (trust_radius > 0.0 && maxratio > trust_radius)
                alpha = trust_radius / maxratio;

            /* 24, not 8: see the second-pass comment (kept verbatim in
             * this function's git history) on why a fast/stiff coordinate
             * can need quite a few halvings before the trial residual
             * actually dips below the baseline. */
            for (bt = 0; bt < 24; bt++) {
                int sane = 1;
                for (k = 0; k < n; k++) {
                    xtrial[k] = x0[k] + alpha * delta_x[k];
                    if (fabs(xtrial[k]) > max_state) sane = 0;
                }
                if (!sane) {
                    if (getenv("PSS_SHOOT_DEBUG"))
                        fprintf(stderr, "DEBUG %s iter=%d bt=%d alpha=%g "
                                "trial exceeds the %.3g sanity ceiling -- "
                                "rejected without running it.\n",
                                stage_label, iter, bt, alpha, max_state);
                    alpha *= 0.5;
                    continue;
                }

                pss_state_set_ic(map, xtrial);
                /* NOT forced -- see the second-pass comment (kept
                 * verbatim in this function's git history): a damped
                 * trial point can be a full Newton step away from x0, so
                 * it wants its own natural adaptive grid. */
                err = pss_run_period_ex(ckt, T, "PSS shoot (damped-step trial)", 0, 0);
                if (!err) {
                    double newresidual = 0.0;
                    pss_state_get(map, ckt, x1trial);
                    for (k = 0; k < n; k++) {
                        double r = fabs(x1trial[k] - xtrial[k]) / scale[k];
                        if (r > newresidual) newresidual = r;
                    }
                    if (getenv("PSS_SHOOT_DEBUG"))
                        fprintf(stderr, "DEBUG %s iter=%d bt=%d alpha=%g "
                                "newresidual=%.6e\n", stage_label, iter, bt,
                                alpha, newresidual);
                    if (newresidual < residual) {
                        accepted = 1;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if (!accepted) {
                if (!recompute) {
                    /* Broyden's approximate direction did not help --
                     * not this attempt's honest answer yet: force a
                     * fresh FD sweep from the SAME x0 and retry, mirroring
                     * the singular-matrix retry above. */
                    fprintf(stdout, "PSS-SHOOT: %s iter=%d Broyden step did "
                            "not reduce the residual; forcing a fresh "
                            "Jacobian.\n", stage_label, iter);
                    recompute = 1;
                    iter--;
                    continue;
                }
                fprintf(stderr, "PSS-SHOOT: %s the damped Newton step did "
                        "not reduce the orbit residual at iteration %d "
                        "this attempt.\n", stage_label, iter);
                result = 0;
                goto newton_done;
            }

            /* Broyden bookkeeping on the RAW monodromy estimate M (see
             * this function's header comment for why not on Phi-I): the
             * secant target for the CHANGE IN P (not in F) is
             * dP = P(xtrial)-P(x0) = (F(xtrial)-F(x0)) + s, since
             * F(x)=P(x)-x exactly. x1trial (= P(xtrial), already run
             * above for the line search) becomes next iteration's
             * baseline P(x0) for free -- no extra period run needed for a
             * Broyden-updated iteration. */
            for (k = 0; k < n; k++) {
                s[k] = xtrial[k] - x0[k];
                y[k] = (x1trial[k] - xtrial[k]) - F[k];   /* y = F(xtrial)-F(x0) */
            }
            if (broyden_enabled) {
                double sTs = 0.0;
                for (k = 0; k < n; k++) sTs += s[k] * s[k];
                if (sTs > 1e-300) {
                    for (k = 0; k < n; k++) {
                        double rowsum = 0.0;
                        for (j = 0; j < n; j++) rowsum += M[k * n + j] * s[j];
                        Ms[k] = rowsum;
                    }
                    for (k = 0; k < n; k++) {
                        double dP_k = y[k] + s[k];   /* dP = y_F + s */
                        double coef = (dP_k - Ms[k]) / sTs;
                        for (j = 0; j < n; j++)
                            M[k * n + j] += coef * s[j];
                    }
                }
            }

            /* x1trial/xtrial are next iteration's baseline P(x0)/x0 for
             * free (no extra period run) -- but residual/scale must still
             * be RECOMPUTED from them (pure arithmetic, no run) exactly as
             * the top-of-loop baseline block would: the original,
             * pre-third-pass code recomputed scale fresh from the CURRENT
             * x0/x1 every iteration (only holding it fixed WITHIN one
             * iteration's own backtracking trials, per the comment above),
             * so reusing free P(xtrial) data without also refreshing
             * residual/scale here would silently freeze both at their
             * iteration-0 values -- found exactly this way on the R2
             * fixture (n=1) during this pass's own testing: the debug F/r
             * printed correctly shrinking every iteration while the
             * convergence check kept comparing a STALE, never-updated
             * residual and looped for the full 25-iteration budget
             * without ever seeing it drop below tol. */
            for (k = 0; k < n; k++) {
                x0[k] = xtrial[k];
                x1[k] = x1trial[k];
                F[k]  = x1trial[k] - xtrial[k];
            }
            residual = 0.0;
            for (k = 0; k < n; k++) {
                double rng = fabs(x1[k]) > fabs(x0[k]) ? fabs(x1[k]) : fabs(x0[k]);
                const double ABSTOL_FLOOR = 1e-3;
                scale[k] = (rng > ABSTOL_FLOOR) ? rng : ABSTOL_FLOOR;
                {
                    double r = fabs(F[k]) / scale[k];
                    if (r > residual) residual = r;
                }
            }
            /* analytic mode: force a fresh sensitivity-tracked baseline
             * run next iteration (from the now-updated x0) rather than
             * reusing this trial's "free" x1trial -- see this function's
             * header comment (use_analytic branch) for why: the Jacobian
             * can ONLY come from a run with PSSsensActive set, and this
             * damped-trial run above did not have it active. */
            have_baseline = use_analytic ? 0 : 1;

            since_recompute++;
            recompute = (!broyden_enabled) ||
                        (broyden_period > 0 && since_recompute >= broyden_period);
        }
    }

 newton_done:
    if (result) {
        *iter_out = iter;
        *residual_out = residual;
    } else {
        *iter_out = iter;
        *residual_out = residual;
    }
    if (Phi_out && M)
        for (k = 0; k < n * n; k++) Phi_out[k] = M[k];
    if (x1) tfree(x1);
    if (F) tfree(F);
    if (scale) tfree(scale);
    if (delta_x) tfree(delta_x);
    if (xtrial) tfree(xtrial);
    if (x1trial) tfree(x1trial);
    if (A) tfree(A);
    if (b) tfree(b);
    if (M) tfree(M);
    if (pert) tfree(pert);
    if (s) tfree(s);
    if (y) tfree(y);
    if (Ms) tfree(Ms);
    return result;
}


/* ------------------------------------------------------------------ */

int
PSSshootDriven(CKTcircuit *ckt, PSSan *job)
{
    double T = 1.0 / ckt->CKTguessedFreq;
    PSSstateMap map = { 0, NULL };
    int n, k, j, iter, err, maxiter;
    double tol, eps, residual = 0.0;
    double t_wall_start, t_wall_end;
    int converged_flag = 0;
    int use_analytic = 0;   /* card #458 fourth pass: PSS_SHOOT_JACOBIAN=analytic (set below) */
    double *x0 = NULL, *x1 = NULL, *F = NULL, *Phi = NULL, *delta_x = NULL,
           *xtrial = NULL, *x1trial = NULL, *scale = NULL, *A = NULL, *b = NULL;
    double save_reltol, save_abstol, save_volttol;
    PSSfloquet floq;
    PSScomplex *eig = NULL;

    t_wall_start = SPfrontEnd->IFseconds();

    /* Integration resolution within a period: match the sampling the
     * deck asked for (Points), or a fixed default. */
    ckt->CKTmaxStep  = (job->PSSpoints > 0) ? T / (double) job->PSSpoints : T / 200.0;
    ckt->CKTstep     = ckt->CKTmaxStep;
    ckt->CKTdelmin   = 1e-9 * ckt->CKTmaxStep;
    /* CKTminBreak: the first pass set this to 0.0 (every breakpoint kept,
     * however close together). SUSPECTED CAUSE 2's forced-breakpoint fix
     * (this card) reinstates every one of a baseline run's accepted times
     * -- possibly hundreds, tightly clustered right around a smooth
     * tanh switch transition where the baseline's own LTE control was
     * already stepping near CKTdelmin -- as breakpoints on the perturbed
     * runs. With CKTminBreak=0.0 that produced its own failure on the R3
     * fixture ("Timestep too small" at the pwm node): forcing exact
     * arrival at a long run of breakpoints separated by less than a
     * numerically meaningful interval fights the corrector instead of
     * helping it. ngspice's own plain-`.tran` default (CKTmaxStep*5e-5,
     * see traninit.c) lets CKTsetBreak()'s existing dedup logic merge
     * breakpoints closer together than that -- restoring it (instead of
     * 0.0) fixed the R3 failure while keeping the forced grid everywhere
     * that actually matters (see the card #458 report). */
    ckt->CKTminBreak = ckt->CKTmaxStep * 5e-5;
    ckt->CKTinitTime = 0.0;

    err = pss_build_state_map(ckt, &map);
    if (err)
        return err;
    n = map.n;
    if (n == 0) {
        fprintf(stderr, "Error: .pss method=shoot found no capacitor or "
                "inductor in the circuit -- nothing to shoot on.\n");
        return E_UNSUPP;
    }

    /* CKTreltol/abstol/voltTol saved here (single funnel-point restore at
     * shoot_cleanup). CKTreltol IS tightened below (to 1e-5, AFTER
     * pre-stabilization only -- see that comment for why not here and not
     * further, and why not abstol/voltTol). TRIED and reverted on this
     * card: tightening CKTreltol to 1e-6 (not 1e-5) made even a single,
     * UNFORCED R3 period run fail outright ("Timestep too small" at the
     * pwm node) -- the steep tanh switch transition there is already
     * close to the edge of what CKTdelmin (1e-9*maxStep here) can resolve
     * at the default tolerance, and 1e-6 pushed it over that edge; 1e-5
     * (below) does not. CKTvoltTol was also tried tightened (to 1e-9) and
     * reverted for the same failure mode, so it is left alone entirely.
     * See the card #458 report for the eps/reltol sweep (PSS_SHOOT_EPS /
     * PSS_SHOOT_RELTOL / PSS_SHOOT_NO_FORCE debug env vars, still
     * present, used to produce it). */
    save_reltol  = ckt->CKTreltol;
    save_abstol  = ckt->CKTabstol;
    save_volttol = ckt->CKTvoltTol;
    memset(&floq, 0, sizeof floq);

    x0      = TMALLOC(double, n);
    x1      = TMALLOC(double, n);
    F       = TMALLOC(double, n);
    Phi     = TMALLOC(double, n * n);
    delta_x = TMALLOC(double, n);
    xtrial  = TMALLOC(double, n);
    x1trial = TMALLOC(double, n);
    scale   = TMALLOC(double, n);
    A       = TMALLOC(double, n * n);
    b       = TMALLOC(double, n);
    eig     = TMALLOC(PSScomplex, n);

    if (!x0 || !x1 || !F || !Phi || !delta_x || !xtrial || !x1trial || !scale || !A || !b || !eig) {
        err = E_NOMEM;
        goto shoot_cleanup;
    }

    /* THIRD PASS (card #458, plan T4.1): staged warm-up escalation,
     * replacing the second pass's single "always burn the whole stabtime"
     * pre-stabilization transient below. Why a warm-up is needed at all
     * (kept from the second pass): the deck's raw IC=/.ic values (usually
     * all zero on the reactive states) are very often themselves an
     * exact, degenerate fixed point of the one-period map -- a converter
     * sitting powered-off with its compensator/comparator at whatever
     * CKTop finds from a cold start is self-consistent and non-switching,
     * so P(0)=0 and Newton "converges" to that in a couple of iterations
     * without ever seeing the real orbit (found running the first R3
     * smoke test this way: it reported CONVERGED with v(vout)~1e-7, not
     * ~3.3V). A plain `.tran` from the same cold start does NOT get stuck
     * there, because its comparator/ramp/op-amp dynamics carry continuous
     * history across many periods; a single one-period Newton evaluation
     * has none of that.
     *
     * Why STAGE it (this pass's own change, task 1 -- speed): `stabtime`
     * is a `.pss` card parameter the deck authors calibrated to how long
     * a PLAIN `.tran` needs to fully settle the orbit's amplitude AND
     * phase (refs/r3_vm_buck.cir etc. reuse it for exactly that, unrelated
     * to this shooter). That is a far more conservative warm-up than
     * Newton actually needs: once the state is INSIDE Newton's basin of
     * attraction -- not necessarily fully settled -- quadratic convergence
     * covers the remaining distance in a handful of iterations, which is
     * exactly what the second pass's own report already showed (R3: 4
     * iterations, R4: 5) once it got a working starting point. So: try a
     * SHORT warm-up first (tens of periods, not hundreds), and only pay
     * for more if that probe was not enough -- bounded above by the
     * deck's own stabtime, i.e. never worse than the second pass's
     * approach in the worst case (a circuit that genuinely needs the
     * whole stabtime, or genuinely never converges at all -- the D~0.55/
     * x20 red-proof decks -- ends up running it anyway, in stages, before
     * reporting the same honest verdict).
     *
     * Every stage is an INTEGER number of clock periods (kept that way
     * below), so a fresh DCtran() restarting at t=0 for the next, longer
     * stage reproduces the same instant in the clock-referenced waveforms
     * (Vramp et al.) that truly continuing the trajectory would -- not an
     * approximation, since every reactive state that carries memory
     * across periods is already in x (this file's top-of-file design
     * comment); the only "state" a restart loses is the phase of periodic
     * sources synchronised to the deck's own t=0, which an integer-period
     * stage boundary preserves exactly (a period-periodic source is back
     * in phase with its own t=0 after any whole number of periods).
     *
     * A non-final stage's own Newton attempt gets a SHORT iteration
     * budget (a fast "is this warm-up already enough" probe,
     * PSS_SHOOT_WARMUP_MAXITER, default 6) and a trust-region cap on its
     * first line-search trial (PSS_SHOOT_TRUST_RADIUS, default 5
     * scale-units -- see pss_newton_run()'s header comment for why). The
     * FINAL stage (cumulative warm-up == the deck's own stabtime) gets
     * the deck's full Newton budget and NO trust-region cap -- i.e. is
     * byte-for-byte the second pass's own single attempt, kept as the
     * always-available fallback so a circuit that genuinely needs the
     * full stabtime sees exactly the second pass's tested behaviour, and
     * the red-proof decks see exactly the second pass's tested (correct)
     * non-convergence verdict, just possibly after some cheaper stages
     * were tried and correctly judged insufficient first.
     *
     * A stage that neither converges nor is final reverts x0 to the
     * WARM-UP TRANSIENT's own endpoint (discarding that stage's failed
     * Newton wandering, not building on it) before extending the warm-up
     * -- so a bad probe never contaminates the next, longer one.
     *
     * See the card #458 (third pass) report for the speed measured this
     * way on R2/R3/R4 against the settled-.tran definition (plan item 1),
     * and PSS_SHOOT_WARMUP_PERIODS/_MAXITER/_GROWTH, PSS_SHOOT_TRUST_RADIUS,
     * PSS_SHOOT_BROYDEN/_PERIOD for the sweep this pass's defaults were
     * chosen from (all debug env vars, no effect unless set, same
     * pattern as the existing PSS_SHOOT_EPS/_RELTOL/_NO_FORCE/_DEBUG). */
    {
        double warmup_periods  = 20.0;
        int    warmup_maxiter  = 6;
        double trust_radius    = 5.0;
        double warmup_growth   = 2.0;
        int    broyden_enabled = 1;
        int    broyden_period  = 4;

        if (getenv("PSS_SHOOT_WARMUP_PERIODS")) warmup_periods = atof(getenv("PSS_SHOOT_WARMUP_PERIODS"));
        if (warmup_periods < 1.0) warmup_periods = 1.0;
        if (getenv("PSS_SHOOT_WARMUP_MAXITER")) warmup_maxiter = atoi(getenv("PSS_SHOOT_WARMUP_MAXITER"));
        if (warmup_maxiter < 1) warmup_maxiter = 1;
        if (getenv("PSS_SHOOT_TRUST_RADIUS")) trust_radius = atof(getenv("PSS_SHOOT_TRUST_RADIUS"));
        if (getenv("PSS_SHOOT_WARMUP_GROWTH")) warmup_growth = atof(getenv("PSS_SHOOT_WARMUP_GROWTH"));
        if (warmup_growth < 1.01) warmup_growth = 1.01;
        if (getenv("PSS_SHOOT_BROYDEN")) broyden_enabled = atoi(getenv("PSS_SHOOT_BROYDEN"));
        if (getenv("PSS_SHOOT_BROYDEN_PERIOD")) broyden_period = atoi(getenv("PSS_SHOOT_BROYDEN_PERIOD"));

        /* card #458 FOURTH pass (plan T4.2): opt-in analytic (discretise-
         * then-linearise) sensitivity Jacobian, replacing the finite-
         * difference monodromy above -- see PSSsensRecordStep()'s header
         * comment for the derivation. Default ("fd", or unset) is
         * byte-for-byte the third pass's tested FD+Broyden path;
         * PSS_SHOOT_JACOBIAN=analytic switches every pss_newton_run() call
         * below (every warm-up stage, not just the final one -- the
         * analytic path is CHEAPER per iteration than FD, ~1-2 period runs
         * vs 2n+3, so there is no reason to withhold it from the cheap
         * probe stages either). */
        use_analytic = (getenv("PSS_SHOOT_JACOBIAN") &&
                        strcmp(getenv("PSS_SHOOT_JACOBIAN"), "analytic") == 0);

        /* Tightened reltol for Newton's OWN (single-period) runs only --
         * same value (1e-5), same rationale as the second pass. THIRD-PASS
         * CHANGE: no longer set once for the rest of the function. With
         * staged warm-up, a "rest of the function" tightening would also
         * apply to every warm-up TRANSIENT stage after the first -- and a
         * multi-period warm-up transient at this tighter tolerance fails
         * outright on R3 ("Timestep too small" at the pwm switch node,
         * found running this pass's own first build: a 20-period probe
         * stage that never even reached Newton). So the tightening is now
         * bracketed tightly around each pss_newton_run() call below (set
         * just before, restored to save_reltol just after) -- every
         * warm-up transient, including stage 2+, runs at the deck's own
         * (looser) tolerance, exactly like the second pass's single
         * pre-stabilization transient did. PSS_SHOOT_RELTOL still
         * overrides; never reaches abstol/voltTol. */
        double tight_reltol = getenv("PSS_SHOOT_RELTOL") ? atof(getenv("PSS_SHOOT_RELTOL")) : 1e-5;
        if (tight_reltol > save_reltol || tight_reltol <= 0.0) tight_reltol = save_reltol;

        maxiter = (job->sc_iter > 0) ? job->sc_iter : 20;
        tol     = (job->steady_coeff > 0.0) ? job->steady_coeff : 1e-6;
        /* Finite-difference perturbation size for the monodromy columns --
         * unchanged from the second pass (see that pass's comment, kept
         * in this file's git history, for the eps/reltol sweep this value
         * came from). PSS_SHOOT_EPS still overrides. */
        eps     = 1e-5;
        if (getenv("PSS_SHOOT_EPS"))
            eps = atof(getenv("PSS_SHOOT_EPS"));

        fprintf(stdout, "\nPSS-SHOOT: shooting-Newton PSS, period T=%g s "
                "(f=%g Hz), %d state(s):", T, ckt->CKTguessedFreq, n);
        for (k = 0; k < n; k++)
            fprintf(stdout, " %s", map.refs[k].name);
        fprintf(stdout, "; max %d Newton iterations, tol %g\n", maxiter, tol);
        fflush(stdout);

        if (ckt->CKTstabTime <= 0.0) {
            pss_state_get_ic(&map, x0);
            ckt->CKTreltol = tight_reltol;
            converged_flag = pss_newton_run(ckt, &map, T, n, x0, maxiter,
                                             tol, eps, -1.0, broyden_enabled,
                                             broyden_period, "[no warm-up]",
                                             &iter, &residual, Phi, use_analytic);
            ckt->CKTreltol = save_reltol;
        } else {
            double stage_len = warmup_periods * T;
            double cumulative = 0.0;
            double *x0_stage_start = TMALLOC(double, n);
            int stage = 0;

            if (!x0_stage_start) {
                err = E_NOMEM;
                goto shoot_cleanup;
            }
            pss_state_get_ic(&map, x0);

            for (;;) {
                int last_stage, this_maxiter;
                double this_trust;
                char label[64];

                stage++;
                if (stage_len > ckt->CKTstabTime - cumulative)
                    stage_len = ckt->CKTstabTime - cumulative;
                cumulative += stage_len;
                last_stage = (cumulative >= ckt->CKTstabTime - 1e-15 * ckt->CKTstabTime);

                fprintf(stdout, "\nPSS-SHOOT: pre-stabilizing (stage %d%s) "
                        "with a %.6g s transient (~%.0f periods, %.6g s "
                        "cumulative of %.6g s stabtime)...\n", stage,
                        last_stage ? ", final" : ", probe", stage_len,
                        stage_len / T, cumulative, ckt->CKTstabTime);
                fflush(stdout);

                pss_state_set_ic(&map, x0);
                err = pss_run_period(ckt, stage_len,
                                      "PSS shoot (pre-stabilization transient)");
                if (err) {
                    fprintf(stderr, "PSS-SHOOT: FAILED -- the pre-"
                            "stabilization transient (stage %d, %.6g s) did "
                            "not complete (ngspice error %d).\n", stage,
                            stage_len, err);
                    tfree(x0_stage_start);
                    goto shoot_cleanup;
                }
                pss_state_get(&map, ckt, x0);
                for (k = 0; k < n; k++) x0_stage_start[k] = x0[k];

                this_maxiter = last_stage ? maxiter : warmup_maxiter;
                this_trust   = last_stage ? -1.0 : trust_radius;
                snprintf(label, sizeof label, "[stage %d/%s]", stage,
                         last_stage ? "final" : "probe");
                ckt->CKTreltol = tight_reltol;
                converged_flag = pss_newton_run(ckt, &map, T, n, x0,
                                                 this_maxiter, tol, eps,
                                                 this_trust, broyden_enabled,
                                                 broyden_period, label,
                                                 &iter, &residual, Phi, use_analytic);
                ckt->CKTreltol = save_reltol;
                if (converged_flag || last_stage)
                    break;

                for (k = 0; k < n; k++) x0[k] = x0_stage_start[k];
                fprintf(stdout, "PSS-SHOOT: stage %d (%.6g s, ~%.0f "
                        "periods) was not enough; extending the warm-up.\n",
                        stage, stage_len, stage_len / T);
                fflush(stdout);
                stage_len *= warmup_growth;
            }
            tfree(x0_stage_start);
        }

        /* Re-tighten for the Floquet re-evaluation and the final
         * trajectory-output run below, same as the second pass (which
         * left reltol tightened from before its own single Newton loop
         * all the way through both of those) -- the per-call bracketing
         * above only needed to keep the WARM-UP TRANSIENT stages loose;
         * once staging is done, tight_reltol is what the rest of this
         * function should run at. Restored to save_reltol only at
         * shoot_cleanup, as before. */
        ckt->CKTreltol = tight_reltol;
    }

    if (!converged_flag) {
        fprintf(stderr, "PSS-SHOOT: FAILED to converge within %d Newton "
                "iterations (final residual %.6e, tol %.3e); this is NOT "
                "a steady state.\n", maxiter, residual, tol);
        goto shoot_fail;
    }

    /* Floquet report: recompute Phi once more at x0 (now converged), with
     * a perturbation size tuned separately from the Newton loop's own eps
     * (see below). SUSPECTED CAUSE 2's fix forces this run and its +-h
     * perturbations onto the SAME accepted time grid
     * (pss_run_period_ex(..., force=1) below, using the grid the Floquet
     * baseline run just recorded), removing the adaptive-step MISMATCH
     * between separate period runs that the first pass's eps=1e-3 was
     * still exposed to (reltol tightening was also tried, to address the
     * OTHER noise source -- per-step solve tolerance -- but reverted; see
     * the comment above save_reltol's capture). Falls back to the
     * Newton-loop's own last Phi if
     * this dedicated pass fails to run. */
    {
        /* Reuse the Newton loop's own eps (not a separate, untested
         * constant): it is already tuned above to be small enough to be
         * accurate under the forced grid, and large enough to clear the
         * per-step solve noise floor at the tightened reltol. A separate,
         * larger eps_floquet=1e-4 was tried on this card and failed to
         * run at all on R3 ("Timestep too small") -- the same "FD
         * perturbation large enough to shift a forced switch edge"
         * problem eps=1e-2 caused in the Newton loop itself. */
        double eps_floquet = eps;
        double *Phi2 = TMALLOC(double, n * n);
        int floquet_ok = (Phi2 != NULL);

        if (floquet_ok && use_analytic) {
            /* card #458 FOURTH pass: analytic sensitivity, same one-period
             * source as pss_newton_run()'s own analytic baseline -- no FD
             * perturbation runs here either. */
            pss_state_set_ic(&map, x0);
            pss_sens_reset(ckt, &map);
            PSSsensActive = 1;
            if (pss_run_period_ex(ckt, T, "PSS shoot (Floquet baseline, analytic)", 0, 0))
                floquet_ok = 0;
            PSSsensActive = 0;
            if (floquet_ok && pss_sens_get_phi(&map, Phi2))
                floquet_ok = 0;
        } else {

        pss_state_set_ic(&map, x0);
        if (floquet_ok && pss_run_period_ex(ckt, T, "PSS shoot (Floquet baseline)", 1, 0))
            floquet_ok = 0;
        if (floquet_ok)
            pss_state_get(&map, ckt, x1);

        for (j = 0; floquet_ok && j < n; j++) {
            double h = fabs(x0[j]) * eps_floquet;
            if (h < 1e-6) h = 1e-6;

            for (k = 0; k < n; k++) xtrial[k] = x0[k];
            xtrial[j] = x0[j] + h;
            pss_state_set_ic(&map, xtrial);
            if (pss_run_period_ex(ckt, T, "PSS shoot (Floquet perturbation +)", 0, 1)) {
                floquet_ok = 0;
                break;
            }
            pss_state_get(&map, ckt, x1trial);
            for (k = 0; k < n; k++)
                Phi2[k * n + j] = x1trial[k];

            for (k = 0; k < n; k++) xtrial[k] = x0[k];
            xtrial[j] = x0[j] - h;
            pss_state_set_ic(&map, xtrial);
            if (pss_run_period_ex(ckt, T, "PSS shoot (Floquet perturbation -)", 0, 1)) {
                floquet_ok = 0;
                break;
            }
            pss_state_get(&map, ckt, x1trial);
            for (k = 0; k < n; k++)
                Phi2[k * n + j] = (Phi2[k * n + j] - x1trial[k]) / (2.0 * h);
        }
        }

        {
            const double *src = floquet_ok ? Phi2 : Phi;
            if (!floquet_ok)
                fprintf(stderr, "PSS-SHOOT: WARNING -- the dedicated Floquet "
                        "re-evaluation (eps=%g) failed to run; reporting the "
                        "eigenvalues of the last Newton-step monodromy "
                        "(eps=%g) instead, which may be less accurate.\n",
                        eps_floquet, eps);
            if (pss_eigenvalues(src, n, eig) != 0) {
                fprintf(stderr, "PSS-SHOOT: WARNING -- the eigenvalue solver "
                        "did not converge on the monodromy matrix; Floquet "
                        "report is NOT trustworthy this run.\n");
                floq.valid = 0;
            } else {
                pss_classify_floquet(eig, n, &floq);
            }
        }
        if (Phi2) tfree(Phi2);
        /* restore the baseline transient's data (x1/Phi corrupted as
         * scratch above) is not needed further -- nothing below reads
         * them except the final trajectory re-run, which sets its own
         * state explicitly. */
    }
    t_wall_end = SPfrontEnd->IFseconds();

    if (getenv("PSS_SHOOT_DEBUG") && floq.valid) {
        fprintf(stderr, "DEBUG Floquet multipliers (%d):", n);
        for (k = 0; k < n; k++)
            fprintf(stderr, " (%.6g%+.6gi)", eig[k].re, eig[k].im);
        fprintf(stderr, "\n");
    }

    fprintf(stdout, "\nPSS-SHOOT: CONVERGED iter=%d residual=%.6e "
            "spectral_radius=%.6f %s wall_time=%.3fs\n",
            iter, residual, floq.valid ? floq.spectral_radius : -1.0,
            !floq.valid ? "UNKNOWN(eigensolve-failed)" :
            floq.unstable ? (floq.period_doubling ? "UNSTABLE(period-doubling)" :
                              floq.oscillation ? "UNSTABLE(oscillation)" : "UNSTABLE")
                          : "stable",
            t_wall_end - t_wall_start);
    if (floq.valid)
        fprintf(stdout, "PSS-SHOOT: FLOQUET dominant_multiplier=%.6g%+.6gi "
                "angle_deg=%.3f\n", floq.dominant.re, floq.dominant.im,
                atan2(floq.dominant.im, floq.dominant.re) * 180.0 / M_PI);

    /* Emit the converged trajectory as an ordinary time-domain plot, same
     * title the existing driven-mode algorithm uses. */
    pss_state_set_ic(&map, x0);
    err = pss_run_period(ckt, T, "Time Domain Periodic Steady State Analysis");
    if (err) {
        fprintf(stderr, "PSS-SHOOT: WARNING -- converged, but the final "
                "trajectory output run failed (ngspice error %d); the "
                "PSS-SHOOT CONVERGED line above is still accurate for the "
                "orbit itself, just not the emitted plot.\n", err);
    }

    err = OK;

 shoot_cleanup:
    ckt->CKTreltol  = save_reltol;
    ckt->CKTabstol  = save_abstol;
    ckt->CKTvoltTol = save_volttol;
    PSSbreakRecordActive = PSSbreakInstallActive = 0;
    pss_break_free();
    PSSsensActive = 0;
    pss_sens_free();   /* g_sens_map points into &map, about to be freed below */
    if (x0) tfree(x0);
    if (x1) tfree(x1);
    if (F) tfree(F);
    if (Phi) tfree(Phi);
    if (delta_x) tfree(delta_x);
    if (xtrial) tfree(xtrial);
    if (x1trial) tfree(x1trial);
    if (scale) tfree(scale);
    if (A) tfree(A);
    if (b) tfree(b);
    if (eig) tfree(eig);
    if (map.refs) tfree(map.refs);
    return err;

 shoot_fail:
    err = E_ITERLIM;
    goto shoot_cleanup;
}
