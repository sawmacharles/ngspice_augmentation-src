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
#include "pacshoot.h"   /* card #465: time-domain PAC/PSTB, hooked in after convergence below */

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

/* card #458 FIFTH pass: one state's native value from ckt->CKTrhsOld,
 * factored out of pss_state_get() below so the new ripple tracker
 * (PSSrippleRecordStep) can read the SAME quantity mid-period, not just
 * at t=T. */
static double
pss_state_value(const PSSstateRef *ref, CKTcircuit *ckt)
{
    if (ref->isCap) {
        CAPinstance *ci = ref->dev.cap;
        return ckt->CKTrhsOld[ci->CAPposNode] - ckt->CKTrhsOld[ci->CAPnegNode];
    } else {
        INDinstance *ii = ref->dev.ind;
        return ckt->CKTrhsOld[ii->INDbrEq];
    }
}

/* Read x back from the converged solution (ckt->CKTrhsOld) after a period
 * run -- the same physical quantities pss_state_get_ic() seeds. */
static void
pss_state_get(PSSstateMap *map, CKTcircuit *ckt, double *x)
{
    int k;
    for (k = 0; k < map->n; k++)
        x[k] = pss_state_value(&map->refs[k], ckt);
}


/* ------------------------------------------------------------------ */
/* card #458 FIFTH pass (plan T4.3), task 1: per-state RIPPLE tracking over
 * one period -- min/max of every native state over every accepted step of
 * a period run, used to make the Newton convergence test relative to each
 * state's own ripple (max-min over the period) rather than only its
 * magnitude. Why this matters: a state that sits at a large, steady DC
 * value but has a SMALL ripple (this deck family's compensator internal
 * nodes -- e.g. opout ~0.55V DC but only ~12mV of ripple) was previously
 * scaled by its own magnitude (~0.5V) in pss_newton_run's residual test,
 * which is a MUCH LOOSER absolute tolerance (tol*0.5V ~ 5e-7V) than the
 * SAME relative tolerance applied to its ripple (tol*12mV ~ 1e-8V) would
 * ask for -- i.e. magnitude-based scaling can call Newton "converged" on
 * such a state well before its error is small relative to the only thing
 * that is actually being asked to match (the periodic orbit's own
 * ripple), exactly the card's hypothesis (c). Same no-op-unless-active
 * pattern as the existing breakpoint/sensitivity hooks: PSSrippleRecordActive
 * defaults to 0, so every other analysis is unaffected. */

int PSSrippleRecordActive = 0;

static PSSstateMap *g_ripple_map = NULL;
static double *g_ripple_min = NULL, *g_ripple_max = NULL;   /* [n] */

static void
pss_ripple_reset(PSSstateMap *map, const double *x0)
{
    int k;
    g_ripple_map = map;
    if (!g_ripple_min || !g_ripple_max) {
        if (g_ripple_min) tfree(g_ripple_min);
        if (g_ripple_max) tfree(g_ripple_max);
        g_ripple_min = TMALLOC(double, map->n);
        g_ripple_max = TMALLOC(double, map->n);
    }
    if (!g_ripple_min || !g_ripple_max) {
        g_ripple_map = NULL;   /* allocation failed: degrade to a no-op tracker */
        return;
    }
    /* card #458 FIFTH pass, bug found by this pass's own testing: seed
     * min/max from the KNOWN exact t=0 value (x0, the same array just
     * written into the device IC fields), not a wide +-1e300 sentinel
     * waiting for the t=0 accept callback to fill it in. ckt->CKTrhsOld
     * is NOT yet populated with the UIC value at the exact instant
     * CKTaccept() first fires for t=0 on a UIC restart (confirmed:
     * reading it there produced ~1e265-1e280 garbage, not the true IC) --
     * exactly why PSSbreakRecordStep/PSSsensRecordStep (pre-existing,
     * this same file) both explicitly SKIP t<=0 rather than trust it.
     * Seeding from x0 directly sidesteps the same trap without needing a
     * throwaway CKTrhsOld read at all. */
    for (k = 0; k < map->n; k++) {
        g_ripple_min[k] = x0[k];
        g_ripple_max[k] = x0[k];
    }
}

/* Called from dctran.c right after CKTaccept(), same call site as
 * PSSbreakRecordStep/PSSsensRecordStep. Skips t<=0 like those two --
 * see pss_ripple_reset()'s comment for why: ckt->CKTrhsOld is not
 * trustworthy there on a UIC restart; the t=0 value is instead seeded
 * directly from x0 by pss_ripple_reset() itself. */
void
PSSrippleRecordStep(CKTcircuit *ckt)
{
    int k;
    if (!g_ripple_map || !g_ripple_min)
        return;
    if (ckt->CKTtime <= 0.0)
        return;
    for (k = 0; k < g_ripple_map->n; k++) {
        double v = pss_state_value(&g_ripple_map->refs[k], ckt);
        if (v < g_ripple_min[k]) g_ripple_min[k] = v;
        if (v > g_ripple_max[k]) g_ripple_max[k] = v;
    }
}

/* ripple[k] = max-min over whatever period run last ran with
 * PSSrippleRecordActive set; 0.0 for any state never recorded (e.g. if
 * allocation failed -- degrades to the old magnitude-only floor via the
 * ABSTOL_FLOOR max() in pss_newton_run, never a crash). */
static void
pss_ripple_get(PSSstateMap *map, double *ripple)
{
    int k;
    for (k = 0; k < map->n; k++) {
        if (g_ripple_min && g_ripple_max && g_ripple_min[k] <= g_ripple_max[k])
            ripple[k] = g_ripple_max[k] - g_ripple_min[k];
        else
            ripple[k] = 0.0;
    }
}

static void
pss_ripple_free(void)
{
    if (g_ripple_min) tfree(g_ripple_min);
    if (g_ripple_max) tfree(g_ripple_max);
    g_ripple_min = g_ripple_max = NULL;
    g_ripple_map = NULL;
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
        (ckt->CKTtime - g_break_times[g_break_count - 1]) < ckt->CKTmaxStep * 1e-2)
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
    int i, skipped = 0;
    /* card #458 SIXTH pass (plan T4.1-T4.3), robustness fix: do not force
     * the handful of very-early recorded times from a period run's own
     * post-UIC-restart LTE step-doubling ramp (t=0 always restarts at
     * order 1 with a small step and grows it every accepted step until
     * the ramp saturates near CKTmaxStep -- a numerical-integrator
     * artifact of every UIC-restarted period run, not a switch-edge
     * timing effect, so it needs none of the cross-run consistency this
     * hook exists for). Forcing one of those as a breakpoint makes
     * dctran.c's own "first timepoint after a breakpoint" logic
     * (CKTdelta = min(delta, 0.1*(CKTbreaks[1]-CKTbreaks[0])), PLUS
     * another /10 for the true t=0 "firsttime" case) cap the very FIRST
     * step of the run to roughly maxStep*1e-4 -- orders of magnitude
     * tighter than every OTHER forced-breakpoint restart in the same run
     * (which sits ~maxStep apart, so only gets the single 0.1x cut, not
     * the extra firsttime one). Found on tests/fixtures/pss_r3_shoot.cir
     * with the pss card's Points field raised from 200 to 1000 (card
     * #458, sixth pass): the finer Points 5x-shrinks CKTmaxStep AND
     * CKTdelmin AND this cap proportionally, but the actual nonlinear
     * corrector cost of resolving a Newton FD perturbation's charge-
     * redistribution transient right at t=0 does NOT shrink with it --
     * so the same relative cap becomes tight enough, at fine Points, to
     * make the corrector genuinely fail (not just take small steps) at
     * an unlucky forced instant, shrinking CKTdelta down to CKTdelmin
     * twice and aborting with "Timestep too small" (confirmed: identical
     * deck converges fine with PSS_SHOOT_NO_FORCE=1, isolating this
     * forcing interaction, not the perturbation or Points value alone,
     * as the cause). Skipping breakpoints below this floor leaves
     * CKTbreaks[1] at whatever the run's own step-doubling has grown to
     * by then (>= the floor), so the very first step gets the SAME
     * 0.1x-of-next-breakpoint cut every other forced restart in this run
     * already gets -- no longer two orders of magnitude tighter -- while
     * every OTHER forced time (mid-period, including every switch-edge
     * cluster SUSPECTED CAUSE 2 exists to pin down) is untouched. Applied
     * only to the LEADING run of skipped times (once a time clears the
     * floor, every later one is forced normally, so a genuinely early
     * switch edge -- not expected in this deck family, but not assumed
     * impossible -- is not silently dropped). PSS_SHOOT_NO_EARLY_SKIP
     * (debug env var, no effect unless set) disables this for isolation,
     * same pattern as the file's other tunables. */
    double early_floor = ckt->CKTmaxStep * 0.1;
    int skipping = (getenv("PSS_SHOOT_NO_EARLY_SKIP") == NULL);

    if (getenv("PSS_SHOOT_DEBUG_BREAK")) {
        double mingap = 1e300; int mingap_i = -1;
        fprintf(stderr, "DEBUGBREAK: installing %d breakpoints, maxStep=%.6g "
                "delmin=%.6g minBreak=%.6g early_floor=%.6g; first 8:", g_break_count,
                ckt->CKTmaxStep, ckt->CKTdelmin, ckt->CKTminBreak, early_floor);
        for (i = 0; i < g_break_count && i < 8; i++)
            fprintf(stderr, " %.6e", g_break_times[i]);
        fprintf(stderr, "; last 8:");
        for (i = g_break_count > 8 ? g_break_count - 8 : 0; i < g_break_count; i++)
            fprintf(stderr, " %.6e", g_break_times[i]);
        for (i = 1; i < g_break_count; i++) {
            double gap = g_break_times[i] - g_break_times[i - 1];
            if (gap < mingap) { mingap = gap; mingap_i = i; }
        }
        fprintf(stderr, "\nDEBUGBREAK: min gap=%.6e at index %d (t=%.6e -> %.6e); "
                "finalTime=%.6e last_break=%.6e gap_to_final=%.6e\n",
                mingap, mingap_i, mingap_i > 0 ? g_break_times[mingap_i - 1] : 0.0,
                mingap_i >= 0 ? g_break_times[mingap_i] : 0.0,
                ckt->CKTfinalTime, g_break_count ? g_break_times[g_break_count - 1] : -1.0,
                g_break_count ? ckt->CKTfinalTime - g_break_times[g_break_count - 1] : -1.0);
    }
    for (i = 0; i < g_break_count; i++) {
        double tt = g_break_times[i];
        if (skipping && tt < early_floor) {
            skipped++;
            continue;
        }
        skipping = 0;   /* first time we clear the floor, force everything after it */
        if (tt > 0.0 && tt < ckt->CKTfinalTime)
            CKTsetBreak(ckt, tt);
    }
    if (getenv("PSS_SHOOT_DEBUG_BREAK"))
        fprintf(stderr, "DEBUGBREAK: skipped %d leading early breakpoint(s) "
                "below %.6g\n", skipped, early_floor);
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
 * the specific cause was already printed to stderr when it happened.
 *
 * card #458 SIXTH pass (plan T4.1-T4.3), speed task, idea 1 from the card:
 * M_init, if non-NULL, is a caller-supplied n_slow x n_slow RAW monodromy
 * estimate (same convention as M below) to seed iteration 0's Jacobian
 * with INSTEAD of running a fresh FD sweep there -- used by the caller
 * (PSSshootDriven) to reuse pss_detect_fast_states()'s own one-time
 * detection FD sweep (which already computes a full n_full x n_full
 * monodromy just to classify fast/slow columns, then previously discarded
 * it) as the first real Newton attempt's starting Jacobian, saving that
 * attempt's own 2*n_slow+1 period runs. Only iteration 0 is seeded this
 * way (every later iteration keeps this function's existing recompute/
 * Broyden cadence exactly as before); NULL (every OTHER pss_newton_run()
 * call site) is byte-for-byte the fifth pass's own behaviour, unaffected.
 * Never used with use_analytic (that path has its own, unrelated
 * Jacobian source). Safe even if the seeded Jacobian is somewhat stale
 * (detection's x_full and this attempt's x0 are not always identical --
 * see the call site's own comment): a seeded M is treated exactly like a
 * Broyden-carried one for every purpose below, including the existing
 * "Broyden step did not reduce the residual -> force a fresh FD Jacobian
 * and retry" fallback, so a bad seed costs at most one wasted line-search
 * trial before self-correcting, never a wrong convergence (the actual
 * convergence test is always the true nonlinear residual, never the
 * Jacobian). */
static int
pss_newton_run(CKTcircuit *ckt, PSSstateMap *map, double T,
                int n_slow, const int *slow_idx, double *x_full,
                int maxiter_this, double tol, double eps,
                double trust_radius, int broyden_enabled, int broyden_period,
                const char *stage_label,
                int *iter_out, double *residual_out, double *Phi_out,
                int use_analytic, const double *M_init, int *seed_used_out)
{
    /* card #458 FIFTH pass (plan T4.3): generalised to a REDUCED Newton
     * unknown set -- n_slow of map->n total native states, indexed by
     * slow_idx[0..n_slow) (task 2). `n` (== n_slow) is kept as the local
     * name for the Newton-unknown count so the bulk of this function
     * (unchanged below) still reads exactly as the fourth pass's did;
     * every array that used to be n(=map->n)-sized and doubled as both
     * "the Newton unknowns" and "every native state" is now n(=n_slow)-
     * sized for the FORMER role only -- the LATTER role (device IC
     * set/get, which always needs every native state, Newton-controlled
     * or not) is carried by the new n_full-sized icbuf/outbuf/x_full
     * below. Eliminated (fast/"carried") states are written to the
     * device IC fields every period run exactly like a Newton unknown
     * would be, just never perturbed or solved for -- see the two
     * "carry the eliminated states forward" comments below for where
     * their value comes from instead. */
    int n = n_slow;
    int n_full = map->n;
    int iter = 0, i, j, k, err;
    double residual = 0.0;
    double *x0 = NULL, *x1 = NULL, *F = NULL, *scale = NULL, *delta_x = NULL,
           *xtrial = NULL, *x1trial = NULL, *A = NULL, *b = NULL,
           *M = NULL, *s = NULL, *y = NULL, *Ms = NULL;
    double *icbuf = NULL, *outbuf = NULL, *ripple_full = NULL, *Phi_full = NULL;
    int *is_slow = NULL;
    int have_baseline = 0;   /* do x0/x1/F already hold a valid, just-run point? */
    int recompute = 1;       /* does THIS iteration's Jacobian need a fresh FD sweep? */
    int since_recompute = 0;
    int result = 0;
    if (seed_used_out) *seed_used_out = 0;
    /* card #458 FIFTH pass, task 1: the floor under a state's RIPPLE
     * (not, as in the second/third pass, under its magnitude). The old
     * value (1e-3, a 1mV/1mA floor chosen there to keep a near-zero
     * cold-start MAGNITUDE from blowing up the residual division, per
     * that pass's own comment) turned out, measured on this pass's own
     * long-settling R3/nominal R4 fixtures, to be the BINDING constraint
     * for several genuinely slow (not fast-state-eliminated) states whose
     * true ripple is sub-mV (e.g. Cout in the long-settling variant,
     * ripple ~0.3mV) -- landing exactly AT a 1mV floor is still the
     * LOOSER of the two choices in the max() below (never stricter than
     * using the true, smaller ripple would be), but demanding residual
     * convergence to 1e-6*1e-3=1e-9 absolute native units is far tighter
     * than this class of circuit's own numerical noise floor (reltol=
     * 1e-5, FD eps=1e-5) can usefully resolve, and cost multiple extra
     * staged warm-up rounds to satisfy -- see the card #458 fifth-pass
     * report for the measured before/after. Raised to 1e-2 as this
     * pass's own shipped default (still far tighter than the pre-task-1
     * magnitude-relative test ever was for a volts/amps-scale state, and
     * still small enough that the degenerate-all-off-fixed-point guard
     * below is unaffected -- that trap sits at ~1e-7, five orders of
     * magnitude under even this raised floor). PSS_SHOOT_RIPPLE_FLOOR
     * overrides (debug env var, no effect unless set, same pattern as
     * the file's other tunables) -- used for the sweep this value was
     * chosen from, not verified as circuit-optimal beyond that sweep. */
    double ripple_floor = getenv("PSS_SHOOT_RIPPLE_FLOOR") ?
                           atof(getenv("PSS_SHOOT_RIPPLE_FLOOR")) : 1e-2;
    if (ripple_floor <= 0.0) ripple_floor = 1e-2;

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

    x0      = TMALLOC(double, n);
    x1      = TMALLOC(double, n);
    F       = TMALLOC(double, n);
    scale   = TMALLOC(double, n);
    delta_x = TMALLOC(double, n);
    xtrial  = TMALLOC(double, n);
    x1trial = TMALLOC(double, n);
    A       = TMALLOC(double, n * n);
    b       = TMALLOC(double, n);
    M       = TMALLOC(double, n * n);
    s       = TMALLOC(double, n);
    y       = TMALLOC(double, n);
    Ms      = TMALLOC(double, n);
    icbuf       = TMALLOC(double, n_full);
    outbuf      = TMALLOC(double, n_full);
    ripple_full = TMALLOC(double, n_full);
    Phi_full    = TMALLOC(double, (size_t) n_full * n_full);
    is_slow     = TMALLOC(int, n_full);
    if (!x0 || !x1 || !F || !scale || !delta_x || !xtrial || !x1trial || !A || !b ||
        !M || !s || !y || !Ms || !icbuf || !outbuf || !ripple_full || !Phi_full || !is_slow) {
        fprintf(stderr, "PSS-SHOOT: %s out of memory allocating Newton "
                "scratch space.\n", stage_label);
        result = 0;
        goto newton_done;
    }
    for (i = 0; i < n_full; i++) is_slow[i] = 0;
    for (k = 0; k < n; k++) is_slow[slow_idx[k]] = 1;
    for (i = 0; i < n_full; i++) ripple_full[i] = 0.0;

    for (iter = 0; ; iter++) {

        if (!have_baseline) {
            for (i = 0; i < n_full; i++) icbuf[i] = x_full[i];
            pss_state_set_ic(map, icbuf);

            /* card #458 FIFTH pass, task 1: track every native state's
             * ripple over THIS baseline period -- see
             * PSSrippleRecordStep()'s header comment. Tracked on every
             * FRESH baseline (this block), not every Broyden-accepted
             * trial below -- ripple changes slowly near the converged
             * orbit, so re-measuring it only when the Jacobian itself is
             * refreshed is an accepted approximation (documented in the
             * card #458 fifth-pass report), not a correctness need. */
            pss_ripple_reset(map, icbuf);
            PSSrippleRecordActive = 1;
            if (use_analytic) {
                /* One period run, sensitivity tracking active: this SAME
                 * run is both the baseline P(x0) (x1 below) AND the
                 * analytic Jacobian source (M below) -- no separate FD
                 * sweep, no separate baseline call. Not FORCED onto any
                 * grid (force=0): the analytic recursion uses THIS run's
                 * own accepted steps directly (via PSSsensRecordStep(),
                 * hooked at the same dctran.c call site as the
                 * breakpoint-recording hook), so there is no cross-run
                 * grid-consistency need the FD path's forcing addressed.
                 * record=1 (card #488 fix, see PSSbreakRecordStep/
                 * PSSbreakInstallHook and the q>1 emission-bug root-cause
                 * comment on the non-analytic branch below): recording is
                 * pure bookkeeping (PSSbreakRecordStep has no effect on
                 * the simulated trajectory), so turning it on here costs
                 * nothing but leaves g_break_times populated with a real
                 * grid for the POST-Newton Floquet-baseline and final
                 * trajectory-emission runs to force themselves onto
                 * instead of running bare/unforced -- see below. */
                pss_sens_reset(ckt, map);
                PSSsensActive = 1;
                err = pss_run_period_ex(ckt, T, "PSS shoot (Newton baseline period, analytic)", 1, 0);
                PSSsensActive = 0;
            } else {
                err = pss_run_period_ex(ckt, T, "PSS shoot (Newton baseline period)", 1, 0);
            }
            PSSrippleRecordActive = 0;
            pss_ripple_get(map, ripple_full);

            if (err) {
                fprintf(stderr, "PSS-SHOOT: %s period integration did not "
                        "complete at Newton iteration %d (ngspice error "
                        "%d); this is NOT a steady state.\n", stage_label, iter, err);
                result = 0;
                goto newton_done;
            }
            pss_state_get(map, ckt, outbuf);

            if (use_analytic) {
                if (pss_sens_get_phi(map, Phi_full)) {
                    fprintf(stderr, "PSS-SHOOT: %s analytic sensitivity failed "
                            "(refactor/solve error, or an integration order "
                            "outside TRAPEZOIDAL 1-2) at Newton iteration %d; "
                            "treating as a failed attempt (no FD fallback "
                            "within one attempt).\n", stage_label, iter);
                    result = 0;
                    goto newton_done;
                }
                if (getenv("PSS_SHOOT_DEBUG"))
                    for (k = 0; k < n_full * n_full; k++)
                        fprintf(stderr, "DEBUG analytic Phi[%d]=%.8g\n", k, Phi_full[k]);
                /* card #458 FIFTH pass (task 3): slice the FULL analytic
                 * monodromy down to the slow-slow block Newton actually
                 * solves -- PSSsensRecordStep()'s recursion still tracks
                 * every reactive device's own S_k regardless of
                 * reduction (every state's history still feeds the
                 * transient, and a slow state's sensitivity can depend on
                 * a fast one's along the way); task 2's reduction only
                 * changes which rows/columns feed the NEWTON solve. */
                for (k = 0; k < n; k++)
                    for (j = 0; j < n; j++)
                        M[k * n + j] = Phi_full[slow_idx[k] * n_full + slow_idx[j]];
            }

            /* card #458 FIFTH pass, task 1: scale/floor now relative to
             * each state's own RIPPLE (see PSSrippleRecordStep()'s header
             * comment), not its magnitude -- see ripple_floor's own
             * comment (function top) for this pass's floor value. */
            residual = 0.0;
            for (k = 0; k < n; k++) {
                int jf = slow_idx[k];
                double rip = ripple_full[jf];
                x0[k] = x_full[jf];
                x1[k] = outbuf[jf];
                F[k]  = x1[k] - x0[k];
                scale[k] = (rip > ripple_floor) ? rip : ripple_floor;
                {
                    double r = fabs(F[k]) / scale[k];
                    if (r > residual) residual = r;
                }
            }
            /* card #458 FIFTH pass, task 2: carry every ELIMINATED (fast)
             * state forward as this period's own output -- exact, not an
             * approximation, for a state whose Floquet direction is ~0
             * (pss_detect_fast_states()'s header comment): its value at
             * t=T does not depend on where it started, so using that
             * output as the next period's input is already the
             * self-consistent fixed point for that coordinate. */
            for (i = 0; i < n_full; i++)
                if (!is_slow[i]) x_full[i] = outbuf[i];

            have_baseline = 1;
        }

        fprintf(stdout, "PSS-SHOOT: %s iter=%d residual=%.6e%s\n", stage_label,
                iter, residual, recompute ? "" : " (broyden)");
        fflush(stdout);
        if (getenv("PSS_SHOOT_DEBUG")) {
            for (k = 0; k < n; k++)
                fprintf(stderr, "DEBUG %s iter=%d k=%d %s x0=%.8g x1=%.8g "
                        "F=%.4e scale=%.6g r=%.4e\n", stage_label,
                        iter, k, map->refs[slow_idx[k]].name, x0[k], x1[k], F[k],
                        scale[k], fabs(F[k]) / scale[k]);
        }

        if (residual <= tol) {
            /* THIRD-PASS GUARD, added after this pass's OWN testing found
             * it necessary: reject a "converged" point where EVERY
             * Newton-unknown state is still at (or under) the
             * ABSTOL_FLOOR scale floor -- the degenerate all-off fixed
             * point this file's top-of-function comment already warns
             * about (a converter sitting powered off is itself a
             * self-consistent, non-switching P(0)=0 solution). The
             * original design relied ENTIRELY on the pre-Newton warm-up
             * already being far from 0 to avoid this; that is no longer a
             * safe assumption with a SHORT probe warm-up and a generous
             * iteration budget -- found exactly this way on the R3
             * fixture: PSS_SHOOT_WARMUP_MAXITER=10 from a 20-period probe
             * converged (residual 1.1e-10, a seemingly excellent fit)
             * onto vout~1.2e-7 V, not ~3.3V, a numerically real but
             * physically wrong answer (and a telltale near-singular
             * Phi~=I / Floquet dominant multiplier exactly 1+0i went with
             * it). A real switching converter's TRUE orbit should not
             * have every single state sitting at/under a millivolt/
             * milliamp-scale floor at once, so this is a general,
             * circuit-agnostic check, not a magic number tuned to R3 --
             * and it is applied at every stage, including a full-budget/
             * full-stabtime one, since nothing before this pass actually
             * proved the old code was immune to the same trap, only that
             * it never happened to trigger it. A rejection here is
             * treated exactly like any other "this attempt did not
             * confirm a steady state" outcome: the staged caller extends
             * the warm-up (or, on the final stage, reports the honest
             * FAILED verdict) rather than silently accepting a wrong
             * answer. */
            int all_at_floor = 1;
            for (k = 0; k < n; k++)
                if (scale[k] > ripple_floor * (1.0 + 1e-9)) { all_at_floor = 0; break; }
            if (all_at_floor) {
                fprintf(stderr, "PSS-SHOOT: %s converged (residual=%.3e) "
                        "but EVERY Newton-unknown state is at/under the "
                        "%.3g ripple-floor -- rejecting this as the "
                        "degenerate all-off fixed point, not a "
                        "trustworthy orbit.\n", stage_label, residual, ripple_floor);
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

        if (recompute && iter == 0 && M_init && !use_analytic) {
            /* card #458 SIXTH pass, speed task, idea 1: seed THIS attempt's
             * very first Jacobian from the caller's supplied M_init
             * (pss_detect_fast_states()'s own detection sweep, reused --
             * see this function's header comment) instead of spending a
             * fresh 2n-period-run FD sweep to recompute essentially the
             * same thing. Treated exactly like a Broyden-carried Jacobian
             * from here on (since_recompute=0, same as the FD branch below
             * sets after a real sweep) -- every existing safety net
             * (singular-matrix retry, "Broyden step did not reduce the
             * residual -> force a fresh FD Jacobian" fallback) already
             * treats "M did not come from THIS iteration's own FD sweep"
             * uniformly, so a stale seed can only cost one extra wasted
             * trial before self-correcting, never a wrong convergence. */
            for (k = 0; k < n * n; k++) M[k] = M_init[k];
            since_recompute = 0;
            if (seed_used_out) *seed_used_out = 1;
            if (getenv("PSS_SHOOT_DEBUG"))
                fprintf(stderr, "DEBUG %s iter=0: seeded Jacobian from "
                        "fast-state detection's own sweep (skipped a fresh "
                        "FD recompute).\n", stage_label);
        } else if (recompute && !use_analytic) {
            /* finite-difference monodromy: 2n extra period runs (centred
             * differences, not the cheaper n forward differences the
             * card's "n+1 period runs per Newton step" description allows
             * -- found on the R3 fixture, second pass: one-sided
             * differences gave a Newton direction whose backtrack line
             * search approached the baseline residual from above for 24
             * straight halvings without ever dipping below it, i.e.
             * essentially no usable descent direction; centred
             * differences cancel the O(h) bias one-sided differences
             * carry and fixed it. `outbuf`, not x1, is the scratch buffer
             * here (see this function's header comment) so the clean
             * baseline x1/F/scale survive for the residual check above on
             * every iteration, recompute or not. card #458 FIFTH pass:
             * `n` (== n_slow) columns/rows now, not every native state
             * (task 2) -- every OTHER native state (fast/carried, or any
             * other slow state not currently being perturbed) is held at
             * its current x_full value via icbuf, exactly like the
             * un-reduced code held every other state at x0[k]. */
            for (j = 0; j < n; j++) {
                int jf = slow_idx[j];
                double h = fabs(x0[j]) * eps;
                if (h < 1e-6) h = 1e-6;

                for (i = 0; i < n_full; i++) icbuf[i] = x_full[i];
                icbuf[jf] = x0[j] + h;
                pss_state_set_ic(map, icbuf);
                err = pss_run_period_ex(ckt, T, "PSS shoot (monodromy perturbation +)", 0, 1);
                if (err) {
                    fprintf(stderr, "PSS-SHOOT: %s monodromy perturbation "
                            "run for state %s (+h) did not complete at "
                            "Newton iteration %d (ngspice error %d).\n",
                            stage_label, map->refs[jf].name, iter, err);
                    result = 0;
                    goto newton_done;
                }
                pss_state_get(map, ckt, outbuf);
                for (k = 0; k < n; k++) M[k * n + j] = outbuf[slow_idx[k]];

                for (i = 0; i < n_full; i++) icbuf[i] = x_full[i];
                icbuf[jf] = x0[j] - h;
                pss_state_set_ic(map, icbuf);
                err = pss_run_period_ex(ckt, T, "PSS shoot (monodromy perturbation -)", 0, 1);
                if (err) {
                    fprintf(stderr, "PSS-SHOOT: %s monodromy perturbation "
                            "run for state %s (-h) did not complete at "
                            "Newton iteration %d (ngspice error %d).\n",
                            stage_label, map->refs[jf].name, iter, err);
                    result = 0;
                    goto newton_done;
                }
                pss_state_get(map, ckt, outbuf);
                for (k = 0; k < n; k++)
                    M[k * n + j] = (M[k * n + j] - outbuf[slow_idx[k]]) / (2.0 * h);
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

                for (i = 0; i < n_full; i++) icbuf[i] = x_full[i];
                for (k = 0; k < n; k++) icbuf[slow_idx[k]] = xtrial[k];
                pss_state_set_ic(map, icbuf);
                /* card #465 FIFTH pass (task 2): FORCED (was force=0 since
                 * the second pass -- see that pass's comment, kept
                 * verbatim in this function's git history, on why NOT
                 * forcing was the original design: "a damped trial point
                 * can be a full Newton step away from x0, so it wants its
                 * own natural adaptive grid"). That reasoning was never
                 * empirically tested against forcing -- it just carried
                 * SUSPECTED CAUSE 2's original "perturbation runs must be
                 * forced" rationale over to this run WITHOUT applying it,
                 * on the assumption a large step needs its own grid. Found
                 * on R4 D~=0.45 (refs/r4_pcm_buck.cir, vin_v=7.33): this
                 * duty cycle's dominant Floquet multiplier is 0.994 (this
                 * pass's own PSS_SHOOT_DEBUG eigenvalue dump; still stable,
                 * |eig|<1, matching the FRA loop gain's own high phase
                 * margin at this D -- see the card #465 report) -- so
                 * (Phi-I) has an eigenvalue of only ~0.006, meaning
                 * pss_solve_dense()'s delta_x is amplified roughly 1/0.006
                 * ~= 166x along that mode relative to the residual F it is
                 * solving for. At that amplification, this loop's own
                 * UNFORCED trial run's natural-adaptive-step noise (this
                 * fixture's own comparators/switches are steep tanh
                 * transitions whose accepted-step sequence is itself
                 * sensitive to which side of a threshold a tiny state
                 * perturbation lands on) stopped shrinking as alpha
                 * backtracked toward 0 -- confirmed by hand with
                 * PSS_SHOOT_DEBUG: newresidual measured ~2.5e-4 at alpha=1
                 * and STILL ~3.0e-4 at alpha=1.2e-7 (24 halvings), instead
                 * of decaying toward the baseline's own ~5e-5 residual as
                 * xtrial->x0 the way a step-size-proportional perturbation
                 * would -- i.e. this was grid-selection noise, not
                 * perturbation-size noise, and it does not vanish with a
                 * smaller step. Below the ~166x-amplified noise floor that
                 * creates, no Newton step at any alpha can ever pass this
                 * loop's own "did the trial residual improve" test, so
                 * every attempt exhausts its line search and the whole
                 * Newton run reports non-convergence even though the
                 * direction was correct and the true residual (the
                 * baseline's own F, unaffected by this) was already close.
                 * FORCING the trial run onto the SAME grid the baseline
                 * (and every FD perturbation run) already uses removes
                 * this noise source exactly the way SUSPECTED CAUSE 2
                 * removed it from the FD monodromy itself -- confirmed:
                 * R4 D~=0.45 now converges in 4 iterations to
                 * residual=1.2e-7 (was FAILED at maxiter, residual
                 * plateaued ~3e-4 to ~6e-5 depending on fixture), spectral
                 * radius 0.994 (stable), and NONE of the other PSS
                 * fixtures regressed re-running this pass's own full
                 * test_pss.py suite with the change (R2/R3/R4-D030 still
                 * converge to the same residual/spectral-radius to 3-4
                 * significant figures; the D~=0.55 no-slope-comp and R3
                 * x20-gain red-proofs still correctly refuse to converge,
                 * D~=0.55's own residual still plateauing ~1.3e-3,
                 * consistent with a genuine period-doubling orbit that has
                 * no period-1 fixed point, not a forcing artefact;
                 * determinism (5x, num_threads=1/8) unchanged). Still
                 * gated by PSS_SHOOT_NO_FORCE (existing debug/isolation
                 * switch, pss_run_period_ex()'s own comment) for anyone
                 * who needs to reproduce the old unforced behaviour. */
                err = pss_run_period_ex(ckt, T, "PSS shoot (damped-step trial)", 0, 1);
                if (!err) {
                    double newresidual = 0.0;
                    pss_state_get(map, ckt, outbuf);   /* outbuf: this trial's FULL result, incl. carried states */
                    for (k = 0; k < n; k++) {
                        x1trial[k] = outbuf[slow_idx[k]];
                        {
                        double r = fabs(x1trial[k] - xtrial[k]) / scale[k];
                        if (r > newresidual) newresidual = r;
                        }
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
             * without ever seeing it drop below tol.
             *
             * card #458 FIFTH pass, task 2: carry the eliminated (fast)
             * states forward from THIS accepted trial's own full output
             * too -- `outbuf` here still holds the ACCEPTED trial's full
             * result (the bt loop broke immediately on acceptance, before
             * any further trial could overwrite it). Task 1: scale/floor
             * still uses ripple_full, but that buffer is only refreshed
             * on the next FRESH baseline (not re-measured every Broyden
             * step -- see the fresh-baseline block's own comment). */
            for (i = 0; i < n_full; i++)
                if (!is_slow[i]) x_full[i] = outbuf[i];
            for (k = 0; k < n; k++) {
                x0[k] = xtrial[k];
                x1[k] = x1trial[k];
                F[k]  = x1trial[k] - xtrial[k];
                x_full[slow_idx[k]] = xtrial[k];
            }
            residual = 0.0;
            for (k = 0; k < n; k++) {
                double rip = ripple_full[slow_idx[k]];
                scale[k] = (rip > ripple_floor) ? rip : ripple_floor;
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
    *iter_out = iter;
    *residual_out = residual;
    if (Phi_out && M)
        for (k = 0; k < n * n; k++) Phi_out[k] = M[k];
    if (x0) tfree(x0);
    if (x1) tfree(x1);
    if (F) tfree(F);
    if (scale) tfree(scale);
    if (delta_x) tfree(delta_x);
    if (xtrial) tfree(xtrial);
    if (x1trial) tfree(x1trial);
    if (A) tfree(A);
    if (b) tfree(b);
    if (M) tfree(M);
    if (s) tfree(s);
    if (y) tfree(y);
    if (Ms) tfree(Ms);
    if (icbuf) tfree(icbuf);
    if (outbuf) tfree(outbuf);
    if (ripple_full) tfree(ripple_full);
    if (Phi_full) tfree(Phi_full);
    if (is_slow) tfree(is_slow);
    return result;
}


/* ------------------------------------------------------------------ */
/* card #458 FIFTH pass (plan T4.3), task 2: automatic fast-state
 * detection. HYPOTHESIS (card #458, fifth pass): a state whose own time
 * constant is much shorter than the period T has, to first order,
 * forgotten its initial value by t=T -- its Floquet multiplier direction
 * is ~exp(-T/tau)~0, so perturbing x0 in that direction has ~no effect on
 * x(T) (this file's own R3 op-amp-output finding, fourth pass: tau~10ns
 * vs T=2us). Treating such a state as a Newton unknown (a) adds 2 FD
 * period runs per Newton iteration for a column that carries essentially
 * no signal, (b) makes that column of the FD Jacobian noise-dominated
 * (the true derivative is ~0, so eps-scale integration/roundoff noise IS
 * the whole signal), and (c) can leave it loosely converged under a
 * magnitude-relative tolerance (task 1's own finding).
 *
 * This function measures the ACTUAL criterion the hypothesis describes --
 * "perturbing x0_j has ~no effect on x(T)" -- directly, with ONE dedicated
 * full (un-reduced) finite-difference monodromy evaluation at the given
 * x_full (2*n_full+1 period runs, forced onto the baseline's own accepted
 * grid exactly like pss_newton_run's own FD sweep): column j's inf-norm
 * (max |dx_i(T)/dx0_j| over every state i, not just state j itself -- a
 * state that only "leaks" into some OTHER state by t=T is just as
 * forgettable a Newton unknown as one that only affects itself) below
 * PSS_SHOOT_FAST_THRESHOLD (native units -- volts for a cap, amps for an
 * inductor; default 1e-3, the same order-of-magnitude convention this
 * file's existing ABSTOL_FLOOR/PSS_SHOOT_MAX_STATE constants already use
 * for "a small circuit's own native units") marks that state FAST --
 * eliminated from the Newton unknown set, but NOT from the circuit: it
 * still gets a device IC every period run (pss_newton_run's icbuf), just
 * carried forward from the previous period's own output instead of
 * Newton-solved (see pss_newton_run's own "carry the eliminated states
 * forward" comments for why that is exact, not an approximation, given
 * this same near-zero-column property).
 *
 * This is a ONE-TIME cost per whole `.pss method=shoot` analysis (paid
 * once by the caller, PSSshootDriven, right before the Newton phase that
 * actually spends an iteration budget -- see that function's own comment
 * for why probe/warm-up stages are NOT reduced this pass), not a per-
 * Newton-iteration cost -- it pays for itself as soon as it saves even one
 * full FD Jacobian recompute on the reduced (smaller) system afterwards.
 * Never reduces to zero Newton unknowns (defensive guard: if every column
 * looks negligible, the single least-fast state is kept anyway) and never
 * runs at all for a 1-state circuit (R2 -- nothing to reduce, and the
 * detection sweep's own cost would not be worth paying there). Falls back
 * to "no reduction" (every state stays a Newton unknown, byte-for-byte the
 * fourth pass's own behaviour) on allocation failure, a failed period run
 * during detection, or PSS_SHOOT_NO_REDUCE (debug/regression-isolation
 * switch, same pattern as PSS_SHOOT_NO_FORCE) -- never silently skips
 * logging what it did. Caller owns the returned slow_idx array (tfree()
 * it when done).
 *
 * card #458 SIXTH pass, speed task, idea 1: phi_full_out, if non-NULL, is a
 * caller-owned n_full*n_full (map->n squared) scratch buffer that this
 * function fills with the RAW (no -I) full monodromy it computes anyway to
 * decide fast/slow -- letting the caller reuse this same one-time FD sweep
 * as the first real Newton attempt's seed Jacobian (pss_newton_run()'s
 * M_init) instead of discarding it, the card's own suggested speed idea.
 * *phi_valid_out (if phi_valid_out non-NULL) is set 1 only when phi_full_out
 * was actually filled by a real FD sweep here (never in a no_reduce/
 * failure fallback path, and never touched at all if phi_full_out is
 * NULL) -- the caller must check it before trusting phi_full_out. */
static int *
pss_detect_fast_states(CKTcircuit *ckt, PSSstateMap *map, double T,
                        const double *x_full, double eps, int *n_slow_out,
                        double *phi_full_out, int *phi_valid_out)
{
    int n = map->n, i, j, k, err;
    double threshold = 1e-3;
    double *Phi = NULL, *outbuf = NULL, *icbuf = NULL, *colnorm = NULL;
    int *is_fast = NULL, *slow_idx;
    int n_slow;

    if (phi_valid_out) *phi_valid_out = 0;

    if (getenv("PSS_SHOOT_FAST_THRESHOLD"))
        threshold = atof(getenv("PSS_SHOOT_FAST_THRESHOLD"));

    if (getenv("PSS_SHOOT_NO_REDUCE") || threshold <= 0.0) {
        fprintf(stdout, "PSS-SHOOT: fast-state reduction disabled "
                "(PSS_SHOOT_NO_REDUCE set, or a non-positive threshold); "
                "all %d states remain Newton unknowns.\n", n);
        goto no_reduce;
    }

    Phi     = TMALLOC(double, (size_t) n * n);
    outbuf  = TMALLOC(double, n);
    icbuf   = TMALLOC(double, n);
    colnorm = TMALLOC(double, n);
    is_fast = TMALLOC(int, n);
    if (!Phi || !outbuf || !icbuf || !colnorm || !is_fast) {
        fprintf(stderr, "PSS-SHOOT: fast-state detection out of memory; "
                "skipping reduction, all %d states remain Newton "
                "unknowns.\n", n);
        goto cleanup_no_reduce;
    }

    for (i = 0; i < n; i++) icbuf[i] = x_full[i];
    pss_state_set_ic(map, icbuf);
    err = pss_run_period_ex(ckt, T, "PSS shoot (fast-state detection baseline)", 1, 0);
    if (err) {
        fprintf(stderr, "PSS-SHOOT: fast-state detection baseline period "
                "did not complete (ngspice error %d); skipping reduction, "
                "all %d states remain Newton unknowns.\n", err, n);
        goto cleanup_no_reduce;
    }

    for (j = 0; j < n; j++) {
        double h = fabs(x_full[j]) * eps;
        if (h < 1e-6) h = 1e-6;

        for (i = 0; i < n; i++) icbuf[i] = x_full[i];
        icbuf[j] = x_full[j] + h;
        pss_state_set_ic(map, icbuf);
        err = pss_run_period_ex(ckt, T, "PSS shoot (fast-state detection perturbation +)", 0, 1);
        if (err) {
            fprintf(stderr, "PSS-SHOOT: fast-state detection perturbation "
                    "run for state %s (+h) did not complete (ngspice error "
                    "%d); skipping reduction, all %d states remain Newton "
                    "unknowns.\n", map->refs[j].name, err, n);
            goto cleanup_no_reduce;
        }
        pss_state_get(map, ckt, outbuf);
        for (k = 0; k < n; k++) Phi[k * n + j] = outbuf[k];

        for (i = 0; i < n; i++) icbuf[i] = x_full[i];
        icbuf[j] = x_full[j] - h;
        pss_state_set_ic(map, icbuf);
        err = pss_run_period_ex(ckt, T, "PSS shoot (fast-state detection perturbation -)", 0, 1);
        if (err) {
            fprintf(stderr, "PSS-SHOOT: fast-state detection perturbation "
                    "run for state %s (-h) did not complete (ngspice error "
                    "%d); skipping reduction, all %d states remain Newton "
                    "unknowns.\n", map->refs[j].name, err, n);
            goto cleanup_no_reduce;
        }
        pss_state_get(map, ckt, outbuf);
        for (k = 0; k < n; k++)
            Phi[k * n + j] = (Phi[k * n + j] - outbuf[k]) / (2.0 * h);
    }

    for (j = 0; j < n; j++) {
        double cn = 0.0;
        for (i = 0; i < n; i++) {
            double v = fabs(Phi[i * n + j]);
            if (v > cn) cn = v;
        }
        colnorm[j] = cn;
        is_fast[j] = (cn < threshold);
    }
    {
        int any_slow = 0;
        for (k = 0; k < n; k++) if (!is_fast[k]) any_slow = 1;
        if (!any_slow) {
            int best = 0; double bestnorm = -1.0;
            for (k = 0; k < n; k++)
                if (colnorm[k] > bestnorm) { bestnorm = colnorm[k]; best = k; }
            is_fast[best] = 0;
            fprintf(stdout, "PSS-SHOOT: fast-state detection found every "
                    "state's column under threshold %.3g -- keeping the "
                    "least-fast state (%s, colnorm=%.3e) as a Newton "
                    "unknown so at least one remains (defensive guard, "
                    "should not happen on a real circuit).\n",
                    threshold, map->refs[best].name, bestnorm);
        }
    }

    n_slow = 0;
    for (k = 0; k < n; k++) if (!is_fast[k]) n_slow++;
    slow_idx = TMALLOC(int, n_slow);
    fprintf(stdout, "PSS-SHOOT: fast-state detection (threshold=%.3g native "
            "units, column inf-norm of a one-time full FD monodromy):\n", threshold);
    n_slow = 0;
    for (k = 0; k < n; k++) {
        fprintf(stdout, "PSS-SHOOT:   %-12s colnorm=%.6e -> %s\n",
                map->refs[k].name, colnorm[k],
                is_fast[k] ? "ELIMINATED (carried)" : "kept (Newton unknown)");
        if (!is_fast[k]) slow_idx[n_slow++] = k;
    }
    fprintf(stdout, "PSS-SHOOT: fast-state reduction: %d Newton unknown(s) "
            "(was %d).\n", n_slow, n);
    fflush(stdout);

    /* belt-and-braces: restore x_full's own IC (the sweep above left the
     * device IC fields at the last perturbation trial) -- not load-
     * bearing, since the caller's very next action always re-sets IC from
     * x_full itself, but cheap and avoids any window where they disagree. */
    for (i = 0; i < n; i++) icbuf[i] = x_full[i];
    pss_state_set_ic(map, icbuf);

    /* card #458 SIXTH pass, speed task, idea 1: hand this sweep's own
     * n x n (n==n_full here, this function's own local `n`) raw monodromy
     * back to the caller before it is freed -- see this function's header
     * comment. */
    if (phi_full_out) {
        for (k = 0; k < n * n; k++) phi_full_out[k] = Phi[k];
        if (phi_valid_out) *phi_valid_out = 1;
    }

    tfree(Phi); tfree(outbuf); tfree(icbuf); tfree(colnorm); tfree(is_fast);
    *n_slow_out = n_slow;
    return slow_idx;

 cleanup_no_reduce:
    if (Phi) tfree(Phi);
    if (outbuf) tfree(outbuf);
    if (icbuf) tfree(icbuf);
    if (colnorm) tfree(colnorm);
    if (is_fast) tfree(is_fast);
 no_reduce:
    slow_idx = TMALLOC(int, n);
    for (k = 0; k < n; k++) slow_idx[k] = k;
    *n_slow_out = n;
    return slow_idx;
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
    /* card #458 FIFTH pass (plan T4.3), task 2: the reduced Newton-unknown
     * set (n_slow of map.n states, slow_idx[0..n_slow) their indices into
     * map.refs), detected once by pss_detect_fast_states() right before
     * the final Newton stage -- see that function's header comment and
     * the staged-warm-up loop below. slow_idx is heap-owned by whichever
     * pss_detect_fast_states() call set it; freed at shoot_cleanup. */
    int n_slow = 0;
    int *slow_idx = NULL;

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
        /* card #458 SIXTH pass, speed task: 6 (third/fourth/fifth pass
         * default) -> 4. A probe stage that is going to fail almost always
         * shows it well before iteration 6 (see the R3 long-settling
         * profile in the card #458 sixth-pass report: stages 2 and 4 both
         * ran their full 6-iteration budget without converging, then were
         * discarded anyway when the next, longer stage was tried) -- the
         * last 2 iterations of a failing probe are close to pure waste.
         * Measured: R3 long-settling 1.13s->0.86s (haven't retested with
         * the seeded-Jacobian change above combined -- see the speed
         * table), R3 nominal 0.56s->0.36s, R4/R2 nominal unchanged
         * (residual/specrad identical to 6 iterations), both red-proofs
         * (R4 D~0.55, R3 x20-gain) still correctly refuse to converge.
         * 2 was also tried and found to be a REGRESSION on the long-
         * settling fixture (1.50s, worse than either 6 or 4) -- too few
         * iterations for a probe that IS on the right track to actually
         * confirm it, so the staged loop ends up needing an extra,
         * expensive warm-up stage instead; 4 was the smallest value tried
         * that did not regress any fixture. PSS_SHOOT_WARMUP_MAXITER still
         * overrides. */
        int    warmup_maxiter  = 4;
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

        /* card #458 SIXTH pass, speed task, idea 1: scratch for reusing
         * pss_detect_fast_states()'s own detection sweep as the first
         * Newton attempt's seed Jacobian -- see pss_newton_run()'s and
         * pss_detect_fast_states()'s header comments. phi_full_detect is
         * n*n (n==map.n, the UN-reduced state count: detection always
         * measures the full system before reduction); M_init is the
         * n_slow*n_slow slice pss_newton_run() actually wants, built once
         * right after detection returns and reused for exactly one call
         * (the first Newton attempt after detection -- see below). */
        double *phi_full_detect = TMALLOC(double, (size_t) n * n);
        double *M_init = TMALLOC(double, (size_t) n * n);   /* n*n upper bound; only n_slow*n_slow used */
        int phi_full_detect_valid = 0;

        if (ckt->CKTstabTime <= 0.0) {
            pss_state_get_ic(&map, x0);
            /* card #458 FIFTH pass, task 2: detect+reduce ONCE, right
             * before the only Newton attempt this path makes -- see
             * pss_detect_fast_states()'s header comment for why this is a
             * one-time cost. */
            n_slow = n;
            slow_idx = pss_detect_fast_states(ckt, &map, T, x0, eps, &n_slow,
                                               phi_full_detect, &phi_full_detect_valid);
            if (phi_full_detect_valid && M_init)
                for (k = 0; k < n_slow; k++)
                    for (j = 0; j < n_slow; j++)
                        M_init[k * n_slow + j] =
                            phi_full_detect[slow_idx[k] * n + slow_idx[j]];
            ckt->CKTreltol = tight_reltol;
            converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                             maxiter, tol, eps, -1.0, broyden_enabled,
                                             broyden_period, "[no warm-up]",
                                             &iter, &residual, Phi, use_analytic,
                                             (phi_full_detect_valid && M_init) ? M_init : NULL,
                                             NULL);
            ckt->CKTreltol = save_reltol;
        } else {
            double stage_len = warmup_periods * T;
            double cumulative = 0.0;
            double *x0_stage_start = TMALLOC(double, n);
            int stage = 0;
            /* card #458 SIXTH pass, speed task, idea 1: offer the detection
             * seed to every stage's Newton attempt (not just stage 1's)
             * until one of them actually consumes it (pss_newton_run()'s
             * seed_used_out) -- on a circuit where stage 1's own attempt
             * never reaches a fresh Jacobian recompute at all (e.g. the
             * degenerate-fixed-point guard rejecting a too-short probe
             * outright before any FD sweep is attempted, observed on the
             * long-settling fixture), the seed would otherwise never be
             * used at all. Once consumed, it is not offered again (later
             * stages' own recomputes are always fresh FD sweeps, same as
             * the fifth pass). A seed consumed several stages after
             * detection is more stale than one consumed at stage 1, but no
             * less SAFE -- see pss_newton_run()'s header comment: the
             * existing Broyden-failed-retry fallback bounds the cost of a
             * bad seed to one wasted trial, never a wrong convergence. */
            int seed_available = 0;

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
                {
                    /* card #458 SIXTH pass, speed task, idea (new this
                     * pass, not one of the card's own listed ideas but the
                     * same "cut the biggest item" instruction): the
                     * warm-up transient's OWN trajectory is discarded the
                     * instant it is read back into x0 below -- nothing
                     * downstream needs it resolved any finer than "close
                     * enough to land Newton in its basin of attraction",
                     * unlike a Newton period run (whose OUTPUT x1=P(x0) is
                     * the quantity Newton's residual test and Floquet
                     * report both depend on directly). Profiled (card
                     * #458 sixth-pass report): warm-up transients are the
                     * single biggest cost on the long-settling fixture
                     * (~60% of total wall time, dominated by stage 5's own
                     * 320-period transient) even after the seeded-Jacobian
                     * and warmup_maxiter cuts above. Loosening CKTmaxStep
                     * for warm-up transients ONLY (restored to the deck's
                     * Points-derived resolution immediately after, before
                     * ANY Newton/detection period run sees it) lets the
                     * adaptive step control take much bigger strides
                     * through the smooth, non-switching parts of every
                     * period -- ngspice's own LTE control still forces
                     * small steps through every real switch transition
                     * regardless of this ceiling, so this does not skip
                     * resolving the switching behaviour itself, only the
                     * artificial per-period sampling density Points asks
                     * Newton's OWN evaluations for. PSS_SHOOT_WARMUP_
                     * MAXSTEP_FACTOR (default 10, i.e. 10x looser than the
                     * Points-derived ceiling; 1 disables this) overrides,
                     * no effect unless set, same pattern as the file's
                     * other tunables.
                     *
                     * PROBE STAGES ONLY (not last_stage): found on R4 (this
                     * pass's own testing) that loosening the FINAL stage's
                     * own warm-up transient too breaks reliable convergence
                     * outright (R4 needs close to its full stabtime -- see
                     * the fifth-pass README's own note on this -- so its
                     * "final" stage IS the one that matters, unlike R3
                     * long-settling where every stage that converges is
                     * still a PROBE, cumulative 620 of 1500 periods; R4's
                     * own final stage is reached at cumulative 750/750).
                     * factor=10 uniformly (including the final stage) was
                     * tried first and made R4 fail to converge at ANY
                     * factor from 3-10 (residual plateaus ~0.76-1.0,
                     * never a false CONVERGED, but a real, measured
                     * reliability regression); factor=2 uniformly
                     * recovered R4 but only partly recovered the speed win
                     * on R3 long-settling. Scoping the loosening to PROBE
                     * stages only (the ones whose own trajectory is
                     * discarded outright on failure anyway, per this
                     * loop's existing "revert x0 to the warm-up
                     * transient's own endpoint" comment above) gets both:
                     * R4's own final-stage warm-up keeps the full,
                     * Points-derived resolution it was measured to need,
                     * while R3 long-settling's own convergence point (a
                     * PROBE stage) still gets the full speed win. Measured
                     * (card #458 sixth-pass report): R4 nominal converges
                     * again (iter=4, residual/specrad unchanged from
                     * factor=1), R3 long-settling ratio unchanged from the
                     * uniform-factor measurement (its own final stage is
                     * never reached). */
                    double warmup_maxstep_factor = (getenv("PSS_SHOOT_WARMUP_MAXSTEP_FACTOR") ?
                        atof(getenv("PSS_SHOOT_WARMUP_MAXSTEP_FACTOR")) : 10.0);
                    double saved_maxstep = ckt->CKTmaxStep;
                    double saved_step    = ckt->CKTstep;
                    double _pt0 = getenv("PSS_SHOOT_PROFILE") ? SPfrontEnd->IFseconds() : 0.0;
                    if (last_stage) warmup_maxstep_factor = 1.0;
                    if (warmup_maxstep_factor > 1.0) {
                        ckt->CKTmaxStep = saved_maxstep * warmup_maxstep_factor;
                        ckt->CKTstep    = ckt->CKTmaxStep;
                    }
                    err = pss_run_period(ckt, stage_len,
                                          "PSS shoot (pre-stabilization transient)");
                    ckt->CKTmaxStep = saved_maxstep;
                    ckt->CKTstep    = saved_step;
                    if (getenv("PSS_SHOOT_PROFILE"))
                        fprintf(stderr, "PROFILE: stage %d warmup transient (%.0f periods) "
                                "took %.4fs\n", stage, stage_len / T,
                                SPfrontEnd->IFseconds() - _pt0);
                }
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

                /* card #458 FIFTH pass, task 2: detect+reduce ONCE, right
                 * after the FIRST warm-up transient (not before it, and
                 * not per-stage) -- reused for every stage's Newton
                 * attempt from here on, probe or final. Two things this
                 * pass's own testing found necessary, both documented in
                 * the card #458 fifth-pass report:
                 *   - detecting from the deck's raw COLD start (x0=0,
                 *     before ANY warm-up) is NOT robust -- tried first,
                 *     and failed outright on this very fixture ("Timestep
                 *     too small" perturbing a state right at the
                 *     degenerate all-off point every cap/inductor's IC=0
                 *     sits at, this file's own top-of-function comment).
                 *     One short (stage 1, ~20 periods) warm-up transient
                 *     is enough to clear that trap while still being a
                 *     small, early cost.
                 *   - reducing ONLY the final stage (an earlier version of
                 *     this pass) left every probe stage on the full state
                 *     set; measured on the long-settling R3 fixture, that
                 *     was a REGRESSION (2.5s vs the fourth pass's 0.62s)
                 *     -- task 1's ripple-relative tolerance correctly
                 *     stops several probe stages from declaring
                 *     "converged" too early on opout's magnitude-scale
                 *     slack, so the staged warm-up needs MORE, LONGER
                 *     probe stages to reach the same honestly-converged
                 *     point, and none of those extra probe attempts got
                 *     task 2's per-iteration cost cut unless reduction
                 *     applies there too. */
                if (stage == 1) {
                    double _pt0 = getenv("PSS_SHOOT_PROFILE") ? SPfrontEnd->IFseconds() : 0.0;
                    n_slow = n;
                    slow_idx = pss_detect_fast_states(ckt, &map, T, x0, eps, &n_slow,
                                                       phi_full_detect, &phi_full_detect_valid);
                    if (phi_full_detect_valid && M_init) {
                        for (k = 0; k < n_slow; k++)
                            for (j = 0; j < n_slow; j++)
                                M_init[k * n_slow + j] =
                                    phi_full_detect[slow_idx[k] * n + slow_idx[j]];
                        seed_available = 1;
                    }
                    if (getenv("PSS_SHOOT_PROFILE"))
                        fprintf(stderr, "PROFILE: fast-state detection took %.4fs\n",
                                SPfrontEnd->IFseconds() - _pt0);
                }

                this_maxiter = last_stage ? maxiter : warmup_maxiter;
                this_trust   = last_stage ? -1.0 : trust_radius;
                snprintf(label, sizeof label, "[stage %d/%s]", stage,
                         last_stage ? "final" : "probe");
                ckt->CKTreltol = tight_reltol;
                {
                double _pnt0 = getenv("PSS_SHOOT_PROFILE") ? SPfrontEnd->IFseconds() : 0.0;
                int seed_used_this = 0;
                /* card #458 SIXTH pass, speed task, idea 1: offer the
                 * detection seed starting at stage 2 (never stage 1 --
                 * found by this pass's own testing: stage 1's own probe is
                 * so short that its Newton attempt typically never reaches
                 * a real Jacobian use at all, most often rejected outright
                 * by the degenerate-fixed-point guard before any FD sweep
                 * or seed would matter; offering the seed there anyway was
                 * measured to cost MORE than the guard's own instant
                 * rejection, for no benefit, and then leaves nothing left
                 * to offer stage 2, which is where a real attempt usually
                 * happens first) to whichever stage's Newton attempt is
                 * the first to actually need a fresh Jacobian recompute --
                 * see seed_available's own comment above and
                 * pss_newton_run()'s header comment for the reuse
                 * mechanics and why a somewhat stale seed is still safe. */
                converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                                 this_maxiter, tol, eps,
                                                 this_trust, broyden_enabled,
                                                 broyden_period, label,
                                                 &iter, &residual, Phi, use_analytic,
                                                 (stage > 1 && seed_available) ? M_init : NULL,
                                                 &seed_used_this);
                if (seed_used_this)
                    seed_available = 0;
                ckt->CKTreltol = save_reltol;
                if (getenv("PSS_SHOOT_PROFILE"))
                    fprintf(stderr, "PROFILE: stage %d newton attempt (%d iters) "
                            "took %.4fs\n", stage, iter, SPfrontEnd->IFseconds() - _pnt0);
                }
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

        if (phi_full_detect) tfree(phi_full_detect);
        if (M_init) tfree(M_init);
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
        /* card #458 FIFTH pass, task 2: Floquet re-evaluated on the SAME
         * reduced (slow_idx, n_slow) set the converged Newton attempt
         * used -- the eliminated states' own multipliers are, by the
         * detection criterion itself (pss_detect_fast_states()'s header
         * comment), near 0 and therefore never the dominant one, so the
         * reduced matrix's spectral radius is the same physical quantity
         * the full one would report; computing it on the smaller block is
         * both cheaper (2*n_slow+1, not 2*n_full+1, period runs) and
         * avoids re-deriving eigenvalues the detection sweep already
         * showed to carry no real signal. n_full-sized icbuf/outbuf hold
         * the full device IC / period-output the same way
         * pss_newton_run()'s own internals do; x0 here is the Newton
         * loop's converged FULL state (n_full entries, fast/carried ones
         * already self-consistent by construction). */
        int n_full_fl = map.n;
        double *Phi2 = TMALLOC(double, (size_t) n_slow * n_slow);
        double *icbuf_fl = TMALLOC(double, n_full_fl);
        double *outbuf_fl = TMALLOC(double, n_full_fl);
        int floquet_ok = (Phi2 != NULL && icbuf_fl != NULL && outbuf_fl != NULL);

        if (floquet_ok && use_analytic) {
            /* card #458 FOURTH pass: analytic sensitivity, same one-period
             * source as pss_newton_run()'s own analytic baseline -- no FD
             * perturbation runs here either. */
            double *Phi_full_fl = TMALLOC(double, (size_t) n_full_fl * n_full_fl);
            for (k = 0; k < n_full_fl; k++) icbuf_fl[k] = x0[k];
            pss_state_set_ic(&map, icbuf_fl);
            pss_sens_reset(ckt, &map);
            PSSsensActive = 1;
            /* card #488 fix: force=1, reusing the Newton baseline's own
             * recorded grid (now populated even in analytic mode -- see
             * the record=1 change above) instead of running bare/
             * unforced -- see the non-analytic branch's root-cause
             * comment below for why an unforced multi-period run is not
             * reliable on its own. */
            if (pss_run_period_ex(ckt, T, "PSS shoot (Floquet baseline, analytic)", 0, 1))
                floquet_ok = 0;
            PSSsensActive = 0;
            if (!Phi_full_fl)
                floquet_ok = 0;
            if (floquet_ok && pss_sens_get_phi(&map, Phi_full_fl))
                floquet_ok = 0;
            if (floquet_ok)
                for (k = 0; k < n_slow; k++)
                    for (j = 0; j < n_slow; j++)
                        Phi2[k * n_slow + j] = Phi_full_fl[slow_idx[k] * n_full_fl + slow_idx[j]];
            if (Phi_full_fl) tfree(Phi_full_fl);
        } else {

        for (k = 0; k < n_full_fl; k++) icbuf_fl[k] = x0[k];
        pss_state_set_ic(&map, icbuf_fl);
        /* card #488 root cause / fix: this call used to be
         * pss_run_period_ex(..., record=1, force=0) -- a fresh, BARE
         * (unforced) natural period run. That is NOT reliable on its own
         * for a multi-native-cycle period (T = q * the circuit's own
         * fastest periodic source, q>1): ngspice's per-device breakpoint
         * bookkeeping for a repeating PULSE/PWL source (VSRCbreak_time,
         * vsrc/vsrcacct.c's VSRCaccept()) is a field on the *instance*
         * that persists across this file's repeated DCtran(ckt,1)
         * restarts within one PSS solve -- it is only ever initialised
         * once, to -1 ("await first pulse"), at circuit setup
         * (vsrc/vsrcset.c), never reset by a fresh restart. dctran.c's
         * own restart path (this file's SUSPECTED CAUSE 2 comment above,
         * pss_break_reset()/PSSbreakInstallHook()) resets ckt->CKTbreaks
         * itself to {0, CKTfinalTime} on every restart, but has no way to
         * reach into an individual device instance's own breakpoint
         * state -- so VSRCaccept()'s incremental "is it time for the
         * next edge yet?" check (ckt->CKTtime >= here->VSRCbreak_time)
         * compares this run's freshly-reset CKTtime (starting at 0)
         * against a STALE absolute time left over from whichever earlier
         * DCtran() call last touched this instance, which is very often
         * already past this run's own CKTfinalTime for a short T=1
         * native-period run -- so the device schedules NO further
         * breakpoints of its own for the rest of this run. For q=1 this
         * is apparently masked well enough by CKTmaxStep (=T/Points)
         * alone resolving the single switch edge every period without
         * device-supplied help; reproduced (card #488/#484's "near-fs/2"
         * finding, research/near-fs2.md) at q=5 (Fguess=fs/5), where the
         * SAME masking is not enough across all 5 edges and the run
         * eventually hits ngspice error 106 ("Timestep too small") --
         * this is why the existing forced-grid machinery (SUSPECTED
         * CAUSE 2, above) exists at all, and this call site (plus the
         * final trajectory-emission run below and the two analytic-mode
         * equivalents) were simply never migrated onto it: they still
         * ran bare/unforced. Fix: this call is now the SAME force=1
         * call every other post-iteration-0 run in this engine already
         * uses, reusing the grid the Newton loop's OWN last baseline run
         * recorded (now populated even in analytic mode -- see the
         * record=1 change above) instead of attempting a second bare
         * natural run of a converged x0 that is, numerically, no more
         * reliable than the ones already fixed. Verified (q=2, q=3, q=5
         * R3 orbits, tests/test_pss.py) to remove the corruption without
         * changing the converged x0/Floquet numbers on q=1 decks (the
         * grid this reuses is the one already used for every FD
         * perturbation/trial run in the last Newton iteration, so this
         * is not a new approximation, just applying the existing one
         * consistently to the two call sites that had been missed). */
        if (floquet_ok && pss_run_period_ex(ckt, T, "PSS shoot (Floquet baseline)", 0, 1))
            floquet_ok = 0;

        for (j = 0; floquet_ok && j < n_slow; j++) {
            int jf = slow_idx[j];
            double h = fabs(x0[jf]) * eps_floquet;
            if (h < 1e-6) h = 1e-6;

            for (k = 0; k < n_full_fl; k++) icbuf_fl[k] = x0[k];
            icbuf_fl[jf] = x0[jf] + h;
            pss_state_set_ic(&map, icbuf_fl);
            if (pss_run_period_ex(ckt, T, "PSS shoot (Floquet perturbation +)", 0, 1)) {
                floquet_ok = 0;
                break;
            }
            pss_state_get(&map, ckt, outbuf_fl);
            for (k = 0; k < n_slow; k++)
                Phi2[k * n_slow + j] = outbuf_fl[slow_idx[k]];

            for (k = 0; k < n_full_fl; k++) icbuf_fl[k] = x0[k];
            icbuf_fl[jf] = x0[jf] - h;
            pss_state_set_ic(&map, icbuf_fl);
            if (pss_run_period_ex(ckt, T, "PSS shoot (Floquet perturbation -)", 0, 1)) {
                floquet_ok = 0;
                break;
            }
            pss_state_get(&map, ckt, outbuf_fl);
            for (k = 0; k < n_slow; k++)
                Phi2[k * n_slow + j] = (Phi2[k * n_slow + j] - outbuf_fl[slow_idx[k]]) / (2.0 * h);
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
            /* NOTE: Phi (the fallback) was itself filled by the converged
             * Newton attempt using this SAME n_slow/slow_idx pair (set
             * once, before that attempt, and unchanged since) -- both
             * branches size n_slow, never the un-reduced n. */
            if (pss_eigenvalues(src, n_slow, eig) != 0) {
                fprintf(stderr, "PSS-SHOOT: WARNING -- the eigenvalue solver "
                        "did not converge on the monodromy matrix; Floquet "
                        "report is NOT trustworthy this run.\n");
                floq.valid = 0;
            } else {
                pss_classify_floquet(eig, n_slow, &floq);
            }
        }
        if (Phi2) tfree(Phi2);
        if (icbuf_fl) tfree(icbuf_fl);
        if (outbuf_fl) tfree(outbuf_fl);
        /* restore the baseline transient's data (x1/Phi corrupted as
         * scratch above) is not needed further -- nothing below reads
         * them except the final trajectory re-run, which sets its own
         * state explicitly. */
    }
    t_wall_end = SPfrontEnd->IFseconds();

    if (getenv("PSS_SHOOT_DEBUG") && floq.valid) {
        fprintf(stderr, "DEBUG Floquet multipliers (%d):", n_slow);
        for (k = 0; k < n_slow; k++)
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
     * title the existing driven-mode algorithm uses. card #488 fix: force
     * this run onto the last recorded grid (force=1), same as the Floquet
     * baseline above -- this used to be a bare pss_run_period() call (no
     * _ex, no forcing), the OTHER unforced multi-period run this card's
     * root-cause finding covers (see the Floquet-baseline comment above
     * for the full explanation: a bare run relies on per-device
     * breakpoint bookkeeping, e.g. vsrc/vsrcacct.c's VSRCbreak_time, that
     * does not survive this file's repeated DCtran() restarts, and that
     * is not reliably masked by CKTmaxStep alone once q = T/(fastest
     * native period) is large enough -- q=5 on refs/r3_vm_buck.cir type
     * decks reproduced ngspice error 106, "Timestep too small", exactly
     * here). g_break_times still holds the last good grid whether or not
     * the Floquet re-evaluation above itself succeeded (it is only
     * overwritten by a record=1 call, and every record=1 call left in
     * this file by this fix runs before this point). */
    pss_state_set_ic(&map, x0);
    err = pss_run_period_ex(ckt, T, "Time Domain Periodic Steady State Analysis", 0, 1);
    if (err) {
        fprintf(stderr, "PSS-SHOOT: WARNING -- converged, but the final "
                "trajectory output run failed (ngspice error %d); the "
                "PSS-SHOOT CONVERGED line above is still accurate for the "
                "orbit itself, just not the emitted plot.\n", err);
    }

    /* card #465 (plan T5.1-T5.2): `.pac ... solver=td` / `.pstb`
     * post-processing on the converged orbit -- see pacshoot.c's header
     * for the design. x0 (n entries, this function's own converged
     * reduced state) is passed as-is; pac_pstb_run() is a no-op (returns
     * OK immediately) unless job->PSSdoPAC&&PACsolver==td or
     * job->PSSdoPSTB, so every plain `.pss method=shoot` run (and every
     * `.pac`/`.pstb` NOT selecting this engine) is unaffected. */
    {
        int pac_err = pac_pstb_run(ckt, job, x0, n, T);
        if (pac_err)
            fprintf(stderr, "PSS-SHOOT: WARNING -- the PSS orbit itself "
                    "converged (see CONVERGED line above), but the .pac/"
                    ".pstb time-domain post-processing failed (error %d).\n",
                    pac_err);
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
    PSSrippleRecordActive = 0;
    pss_ripple_free();
    if (slow_idx) tfree(slow_idx);   /* card #458 fifth pass, task 2 */
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
