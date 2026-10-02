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
#include "sw/swdefs.h"   /* card #494 (plan T6.1): saltation for ideal `sw` events */
#include "res/resdefs.h"   /* card #494 fifth pass: pss_salt_local_cond() below */
#include "dio/diodefs.h"   /* card #510: pss_salt_ind_jump() diode term below */
#include "vsrc/vsrcdefs.h"   /* card #510: pss_salt_node_is_stable() below */

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
/* card #498 2nd pass: exact event location/lagged sw is used by the FD-Jacobian
 * path only; the analytic (saltation) Jacobian keeps the stock sw decision so
 * baseline and trial runs stay mutually consistent. */
static int g_pss_lag_enabled = 1;
int PSSshootLagEnabled(void) { return g_pss_lag_enabled && getenv("PSS_NO_EVENT_LOCATE") == NULL; }

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

    PSSeventLocateActive = (getenv("PSS_NO_EVENT_LOCATE") == NULL) && g_pss_lag_enabled && !PSSsensActive;
    PSSswLagActive = PSSeventLocateActive;
    err = DCtran(ckt, 1);   /* restart = 1: fresh run from t=0 */
    PSSeventLocateActive = 0;
    PSSswLagActive = 0;

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
static unsigned char *g_break_is_event = NULL;   /* card #498 SECOND pass, see below */
static int g_break_count = 0, g_break_cap = 0;

/* card #505: full, UN-thinned record of the baseline run's accepted steps
 * (time and the integration order used for the step) for exact replay --
 * see PSSreplayStep(). */
int PSSreplayActive = 0;
int PSSlastStepOrder = 1;
static double *g_rp_t = NULL;
static unsigned char *g_rp_ord = NULL;
static int g_rp_n = 0, g_rp_cap = 0, g_rp_idx = 0;

static void
pss_break_reset(void)
{
    g_break_count = 0;
    g_rp_n = 0;
}

static void
pss_break_free(void)
{
    if (g_break_times) tfree(g_break_times);
    if (g_break_is_event) tfree(g_break_is_event);
    g_break_times = NULL;
    g_break_is_event = NULL;
    g_break_count = g_break_cap = 0;
    if (g_rp_t) tfree(g_rp_t);
    if (g_rp_ord) tfree(g_rp_ord);
    g_rp_t = NULL; g_rp_ord = NULL; g_rp_n = g_rp_cap = g_rp_idx = 0;
}

/* card #498 SECOND pass (root cause of the plateau the first pass's
 * report left unresolved -- see docs/pss-pac.md's "card #498" section for
 * the full derivation). Does the accepted step just recorded (CKTstate0,
 * this step, vs CKTstate1, the previous accepted step) correspond to a
 * REAL, conduction-changing state flip on any native `sw` instance -- the
 * SAME on/off classification PSSsaltFindEvents()'s own pss_salt_is_on()
 * uses (a REALLY_ON<->HYST_ON or REALLY_OFF<->HYST_OFF move, staying on
 * the same side of the threshold, is NOT a real event; only an ON<->OFF
 * crossing is). Duplicated here rather than calling pss_salt_is_on()/
 * PSSsaltFindEvents() directly: pss_salt_is_on() is defined later in this
 * file (function-pointer-free forward use would need a prototype either
 * way) and PSSsaltFindEvents() itself carries PSSsaltActive/restart-
 * filter/h_dot-history state that has nothing to do with this,
 * unconditional, grid-consistency use -- this function only asks "did a
 * switch instance's own conduction state change on this accepted step",
 * nothing else. See PSSbreakInstallHook()'s own comment below for WHY
 * this accepted time must then NOT be forced onto every other period run
 * in the same Newton iteration. Self-gating: a deck with no native `sw`
 * instance (CKTtypelook("Switch")<0, or no instance ever flips) always
 * returns 0 here, making this an exact no-op for every existing switch-
 * free fixture -- no separate has-switch gate needed for this part of
 * the card #498 fix (unlike the first pass's retry/shrink fallbacks,
 * which needed pss_ckt_has_switch() because THEY could rescue a state
 * that was never going to converge; this one only withholds a forced
 * breakpoint, which is a strictly narrower change). */
static int
pss_break_step_is_switch_event(CKTcircuit *ckt)
{
    int swtype;
    SWmodel *model;
    SWinstance *here;

    swtype = CKTtypelook("Switch");
    if (swtype < 0)
        return 0;
    for (model = (SWmodel *) ckt->CKThead[swtype]; model; model = SWnextModel(model))
        for (here = SWinstances(model); here; here = SWnextInstance(here)) {
            double now  = ckt->CKTstate0[here->SWswitchstate];
            double prev = ckt->CKTstate1[here->SWswitchstate];
            int now_on  = (now  == 1.0 /* REALLY_ON */ || now  == 3.0 /* HYST_ON */);
            int prev_on = (prev == 1.0 /* REALLY_ON */ || prev == 3.0 /* HYST_ON */);
            if (now_on != prev_on)
                return 1;
        }
    return 0;
}

/* Called from dctran.c right after CKTaccept() when PSSbreakRecordActive. */
void
PSSbreakRecordStep(CKTcircuit *ckt)
{
    int is_event;

    if (ckt->CKTtime <= 0.0)
        return;   /* t=0 is already a breakpoint; nothing to add */

    /* card #505: full accepted-step record for exact replay */
    if (g_rp_n >= g_rp_cap) {
        int newcap = g_rp_cap ? g_rp_cap * 2 : 1024;
        double *nt = TREALLOC(double, g_rp_t, newcap);
        unsigned char *no = TREALLOC(unsigned char, g_rp_ord, newcap);
        if (nt) g_rp_t = nt;
        if (no) g_rp_ord = no;
        if (nt && no) g_rp_cap = newcap;
    }
    if (g_rp_n < g_rp_cap) {
        g_rp_t[g_rp_n] = ckt->CKTtime;
        g_rp_ord[g_rp_n++] = (unsigned char) PSSlastStepOrder;
    }

    is_event = PSSswLagActive ? PSSeventLocateCrossed : pss_break_step_is_switch_event(ckt);

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
     * Newton monodromy estimate anyway. card #498 SECOND pass: a real
     * switch-event point is never thinned away by this clustering check
     * even if it would otherwise qualify -- it is recorded (flagged, see
     * below) unconditionally, so pss_break_step_is_switch_event() always
     * gets a chance to see the exact accepted point that flipped. */
    if (!is_event && g_break_count > 0 &&
        (ckt->CKTtime - g_break_times[g_break_count - 1]) < ckt->CKTmaxStep * 1e-2)
        return;

    if (g_break_count >= g_break_cap) {
        int newcap = g_break_cap ? g_break_cap * 2 : 256;
        double *nt = TREALLOC(double, g_break_times, newcap);
        unsigned char *ne = TREALLOC(unsigned char, g_break_is_event, newcap);
        if (!nt || !ne)
            return;   /* out of memory: degrade to unforced steps past here,
                          rather than crash -- the perturbed run may then
                          reintroduce some noise, but never a wrong answer */
        g_break_times = nt;
        g_break_is_event = ne;
        g_break_cap = newcap;
    }
    g_break_is_event[g_break_count] = (unsigned char) is_event;
    g_break_times[g_break_count++] = ckt->CKTtime;
}

/* card #505: exact step-sequence replay. Called from dctran.c once per FRESH
 * time step (not on a corrector-failure retry) while PSSreplayActive.
 * Overrides the step the breakpoint/LTE logic just chose with the next
 * recorded accepted time of the baseline run and the integration order that
 * step used; dctran.c also skips CKTtrunc() (no LTE rejection, no order
 * change) while replaying. Rationale (see docs/pss-pac.md "card #505"):
 * the older "force every recorded time as a breakpoint" scheme restarts
 * backward Euler at EVERY forced point and re-derives the sub-steps by LTE,
 * so a forced run at the very same x0 does NOT reproduce the free baseline
 * run (measured on R3d: 1e-4 relative in i(L1) at alpha->0), the
 * baseline-vs-trial mismatch floors the Newton residual (3e-4) and the
 * clustered/thinned sub-100 ps points are re-derived per run (map noise).
 * Replay makes P(x0) on the recorded grid a smooth, exactly reproducible
 * function of x0 (to the corrector tolerance). A corrector failure at a
 * replayed step still falls back to dctran.c's own divide-by-8 sub-stepping;
 * the next fresh step then again targets the recorded time. */
void
PSSreplayStep(CKTcircuit *ckt)
{
    double tgt;

    while (g_rp_idx < g_rp_n &&
           g_rp_t[g_rp_idx] - ckt->CKTtime <= ckt->CKTdelmin * 2.0)
        g_rp_idx++;
    if (g_rp_idx >= g_rp_n)
        return;   /* past the record: leave the step as the stock logic chose it */
    tgt = g_rp_t[g_rp_idx];
    ckt->CKTdelta = tgt - ckt->CKTtime;
    ckt->CKTorder = g_rp_ord[g_rp_idx];
    ckt->CKTbreak = 0;
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
     * every OTHER forced time is untouched. Applied only to the LEADING
     * run of skipped times (once a time clears the floor, every later one
     * is forced normally, so a genuinely early switch edge -- not
     * expected in this deck family, but not assumed impossible -- is not
     * silently dropped). PSS_SHOOT_NO_EARLY_SKIP (debug env var, no
     * effect unless set) disables this for isolation, same pattern as
     * the file's other tunables.
     *
     * card #498 SECOND pass: a genuine switch-edge crossing (g_break_
     * is_event[i], see pss_break_step_is_switch_event()'s own header
     * comment) is now NEVER forced here, on ANY run, regardless of the
     * early-skip floor above -- root cause of the plateau the first
     * pass's fallback could only paper over (see docs/pss-pac.md's
     * "card #498" section for the derivation): forcing a perturbed FD
     * column's period run to land at exactly the BASELINE's own switch-
     * crossing time silently discards d(t_event)/d(x0) -- precisely the
     * state-dependent-event sensitivity a comparator-driven switch's
     * monodromy needs (the classical "saltation" contribution) -- and
     * for a large-enough perturbation can additionally force the
     * corrector to try to land at a time the perturbed trajectory's OWN
     * switch decision has already made physically impossible (the
     * switch already flipped, or hasn't yet), producing exactly the
     * "Timestep too small ... trouble with swmod-instance" collapse the
     * first pass's report characterises. Leaving this one instant
     * unforced lets each period run (baseline AND every +h/-h column)
     * locate ITS OWN crossing via its own LTE/hysteresis logic --
     * exactly what a plain, unforced `.tran` of the same circuit already
     * does cleanly (see the first pass's report) -- while every
     * SURROUNDING smooth-region time is still forced as before, so the
     * cross-run grid-consistency benefit SUSPECTED CAUSE 2 was built for
     * is kept everywhere except the one instant it was never valid for.
     * Self-gating (see pss_break_step_is_switch_event()): a switch-free
     * deck never sets g_break_is_event[i], so this is an exact no-op for
     * every existing no-`sw` fixture. */

    double early_floor = ckt->CKTmaxStep * 0.1;
    int skipping = (getenv("PSS_SHOOT_NO_EARLY_SKIP") == NULL);
    int skipped_event = 0;

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
        if (g_break_is_event && g_break_is_event[i]) {
            skipped_event++;
            continue;   /* card #498 SECOND pass: never force a real switch
                           crossing -- see this function's header comment */
        }
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
                "below %.6g, %d switch-event breakpoint(s) (card #498 "
                "SECOND pass, never forced)\n", skipped, early_floor,
                skipped_event);
}

/* ------------------------------------------------------------------ */
/* card #498 SECOND pass: exact switch-event location inside a period run.
 *
 * ngspice's native `sw` flips its conduction state at whatever accepted
 * time point first sees the control voltage past the threshold, so the
 * switching instant is quantised to the adaptive time grid (error up to
 * one step, ~1 % of the period at 200 points). For a comparator that is
 * driven by a state (a closed PWM loop), the period map x0 -> x(T) is then
 * piecewise constant in the event time: a finite-difference monodromy
 * column sees either 0 or a whole-step jump of d(t_event)/d(x0), so the
 * Newton residual plateaus at a Points-dependent floor (1.8e-3 at 200,
 * 9e-4 at 2000). Fix: when a converged trial step flips a switch, reject
 * it and re-run the step to the interpolated crossing time (regula falsi
 * on the control voltage, bracket kept across trials); the flip is
 * accepted only when the control overshoot corresponds to <= ttol seconds.
 * Then the event time is a smooth function of x0 to ~ttol and FD columns
 * with h >> ttol*slope are consistent.
 * Returns 1 = reject this step and retry with *newdelta; 0 = carry on. */
int PSSeventLocateActive = 0;
int PSSeventLocateRestart = 0;   /* set when the step just accepted ended at a located crossing */
int PSSeventLocateCrossed = 0;   /* the accepted step ended at a located crossing (event breakpoint) */
long PSSeventLocateTrials = 0, PSSeventLocateEvents = 0;

static int g_ev_n = 0;
static double *g_ev_clo = NULL, *g_ev_chi = NULL;
static double g_ev_tacc = -1.0, g_ev_tlo = 0.0, g_ev_thi = 0.0;
static double g_ev_flo = 1.0, g_ev_fhi = 1.0;   /* Illinois weights */
static int g_ev_have_hi = 0, g_ev_iter = 0, g_ev_last_side = 0;
static double g_ev_req = -1.0;

/* With PSSswLagActive the switch states are constant over a step, so the
 * control value c_i(t) is smooth. A step "crosses" for instance i when, with
 * respect to its (frozen) state, the control has moved beyond the flip
 * threshold at the trial time but was not beyond it at the last accepted
 * point. Return 1 = reject, retry with *newdelta ending at the estimated
 * crossing (just past it); 0 = accept normally. */
int
PSSeventLocate(CKTcircuit *ckt, double *newdelta)
{
    int swtype = CKTtypelook("Switch");
    SWmodel *model;
    SWinstance *here;
    int n = 0, i, any_cross = 0, accept = 1;
    double t = ckt->CKTtime, tacc = ckt->CKTtime - ckt->CKTdelta;
    double ttol, tbest = 1e300;
    static double ttol_env = -1.0;

    PSSeventLocateRestart = 0;
    PSSeventLocateCrossed = 0;
    if (swtype < 0 || !PSSswLagActive)
        return 0;
    if (ttol_env < 0.0)
        ttol_env = getenv("PSS_EVENT_TTOL") ? atof(getenv("PSS_EVENT_TTOL")) : 1e-16;
    ttol = ttol_env;

    for (model = (SWmodel *) ckt->CKThead[swtype]; model; model = SWnextModel(model))
        for (here = SWinstances(model); here; here = SWnextInstance(here))
            n++;
    if (n == 0)
        return 0;
    if (n != g_ev_n) {
        g_ev_clo = TREALLOC(double, g_ev_clo, n);
        g_ev_chi = TREALLOC(double, g_ev_chi, n);
        g_ev_n = n;
        g_ev_tacc = -1.0;
    }
    if (tacc != g_ev_tacc) {   /* new accepted point: fresh bracket */
        g_ev_tacc = tacc;
        g_ev_tlo = tacc;
        g_ev_have_hi = 0;
        g_ev_iter = 0;
        g_ev_flo = g_ev_fhi = 1.0;
        g_ev_last_side = 0;
        i = 0;
        for (model = (SWmodel *) ckt->CKThead[swtype]; model; model = SWnextModel(model))
            for (here = SWinstances(model); here; here = SWnextInstance(here))
                g_ev_clo[i++] = ckt->CKTstate1[here->SWctrlvalue];
    }
    PSSeventLocateTrials++;
    {   /* did a switch actually change conduction state in this step? then
         * the next step restarts with a flat history (see dctran.c) */
        for (model = (SWmodel *) ckt->CKThead[swtype]; model; model = SWnextModel(model))
            for (here = SWinstances(model); here; here = SWnextInstance(here)) {
                double a = ckt->CKTstate0[here->SWswitchstate], b = ckt->CKTstate1[here->SWswitchstate];
                if (((a == 1.0 || a == 3.0)) != ((b == 1.0 || b == 3.0)))
                    PSSeventLocateRestart = 1;
            }
    }

    /* classify this trial */
    {
        double w = t - g_ev_tlo;
        i = 0;
        for (model = (SWmodel *) ckt->CKThead[swtype]; model; model = SWnextModel(model))
            for (here = SWinstances(model); here; here = SWnextInstance(here)) {
                double st = ckt->CKTstate0[here->SWswitchstate];
                int on = (st == 1.0 || st == 3.0);
                double c = ckt->CKTstate0[here->SWctrlvalue];
                double cacc = ckt->CKTstate1[here->SWctrlvalue];
                double vt = model->SWvThreshold, vh = model->SWvHysteresis;
                double thr = on ? vt - vh : vt + vh;
                double dir = on ? -1.0 : 1.0;
                g_ev_chi[i] = c;
                if (vh > 0.0 && dir * (cacc - thr) <= 0.0 && dir * (c - thr) > 0.0) {
                    double slope = (w > 0.0) ? fabs((c - g_ev_clo[i]) / w) : 0.0;
                    double over = dir * (c - thr);
                    any_cross = 1;
                    if (over > ttol * slope + 1e-15 * (fabs(thr) + 1e-3))
                        accept = 0;
                }
                i++;
            }
    }
    if (!any_cross) {
        if (t > g_ev_tlo) {
            i = 0;
            for (model = (SWmodel *) ckt->CKThead[swtype]; model; model = SWnextModel(model))
                for (here = SWinstances(model); here; here = SWnextInstance(here))
                    g_ev_clo[i++] = ckt->CKTstate0[here->SWctrlvalue];
            if (g_ev_last_side == 0 && g_ev_have_hi) g_ev_fhi *= 0.5; else g_ev_flo = 1.0;
            g_ev_last_side = 0;
            g_ev_tlo = t;
        }
        if (!g_ev_have_hi)
            return 0;          /* ordinary smooth step */
        if (fabs(t - g_ev_req) > 1e-12 * fabs(t)) {
            g_ev_have_hi = 0;  /* LTE shortened the step on its own: ordinary step */
            return 0;
        }
        goto estimate;         /* bracketed and still short of the crossing: refine */
    }

    if (g_ev_last_side == 1 && g_ev_have_hi) g_ev_flo *= 0.5; else g_ev_fhi = 1.0;
    g_ev_last_side = 1;
    g_ev_thi = t;
    g_ev_have_hi = 1;
    if (t - g_ev_tlo < 4.0 * ttol)
        accept = 1;
    if (accept || g_ev_iter >= 80) {
        if (accept) PSSeventLocateEvents++;
        PSSeventLocateCrossed = 1;
        return 0;
    }

estimate:
    {
        double w = g_ev_thi - g_ev_tlo;
        i = 0;
        for (model = (SWmodel *) ckt->CKThead[swtype]; model; model = SWnextModel(model))
            for (here = SWinstances(model); here; here = SWnextInstance(here)) {
                double chi = g_ev_chi[i], clo = g_ev_clo[i];
                double vt = model->SWvThreshold, vh = model->SWvHysteresis;
                double st = ckt->CKTstate0[here->SWswitchstate];
                int on = (st == 1.0 || st == 3.0);
                double thr = on ? vt - vh : vt + vh;
                double dir = on ? -1.0 : 1.0;
                if (vh > 0.0 && dir * (clo - thr) <= 0.0 && dir * (chi - thr) > 0.0) {
                    double flo = (clo - thr) * g_ev_flo, fhi = (chi - thr) * g_ev_fhi;
                    double aim = 0.5 * ttol * fabs((chi - clo) / w) * dir;
                    if (fhi != flo) {
                        double ti = g_ev_tlo + (w * (-flo + aim)) / (fhi - flo);
                        if (ti < tbest && ti > g_ev_tlo && ti < g_ev_thi)
                            tbest = ti;
                    }
                }
                i++;
            }
    }
    g_ev_iter++;
    if (tbest > 1e299)
        tbest = 0.5 * (g_ev_tlo + g_ev_thi);   /* bisection fallback */
    if (tbest - tacc < 2.0 * ckt->CKTdelmin || g_ev_thi - tbest < 1e-15 * g_ev_thi) {
        PSSeventLocateCrossed = 1;
        return 0;     /* cannot refine further: accept */
    }
    *newdelta = tbest - tacc;
    g_ev_req = tbest;
    return 1;
}

/* Wraps pss_run_period(): record=1 captures this run's own accepted grid
 * (resetting whatever was recorded before); force=1 installs the
 * currently-recorded grid as breakpoints on this run instead. Never both:
 * a run either DEFINES the grid for this Newton iteration (the baseline
 * P(x0) evaluation) or FOLLOWS a previously-defined one (every perturbed
 * and trial run in the same iteration). */
/* card #505: replay mode. PSS_REPLAY=1 forces on, 0 off; default (unset) off
 * until validated -- see pss_replay_enabled() callers. */
static int g_replay_mode = -1;   /* -1 = read env */
static int g_pss_has_switch;      /* defined (with its comment) further down */
static int g_replay_auto_used;    /* card #505: the automatic replay rung ran in this analysis */
static int
pss_replay_enabled(void)
{
    if (g_replay_mode < 0)
        g_replay_mode = getenv("PSS_REPLAY") ? atoi(getenv("PSS_REPLAY")) : 0;
    return g_replay_mode > 0;
}

static long g_pss_nruns = 0;
static double g_pss_warm_periods = 0.0;   /* card #537: warm-up transient periods in this analysis */   /* card #537: period runs made through pss_run_period_ex() in this analysis */
static int
pss_run_period_ex(CKTcircuit *ckt, double T, const char *jobname,
                   int record, int force)
{
    int err;

    g_pss_nruns++;

    if (getenv("PSS_SHOOT_NO_FORCE") != NULL)
        force = 0;   /* debug/isolation switch, card #458 second pass */

    if (record)
        pss_break_reset();
    PSSbreakRecordActive  = record;
    if (force && pss_replay_enabled() && g_rp_n > 0 && !g_pss_has_switch) {
        PSSreplayActive = 1;
        g_rp_idx = 0;
        PSSbreakInstallActive = 0;
    } else {
        PSSbreakInstallActive = force;
    }

    err = pss_run_period(ckt, T, jobname);

    PSSbreakRecordActive  = 0;
    PSSbreakInstallActive = 0;
    PSSreplayActive = 0;
    return err;
}


/* card #498: does this circuit contain a native ideal Switch (SW) instance
 * anywhere? Both of this card's fallbacks (pss_run_period_ex_safe()'s
 * looser-reltol retry, and pss_fd_perturb_side()'s shrinking-h retry) exist
 * ONLY to rescue the specific numerical pathology an ideal `sw` device's
 * hard threshold creates right at a saltation event (see those two
 * functions' own header comments, and patches-wip/salt/README.md, for the
 * full characterisation) -- found, on this card's own repro
 * (tests/fixtures/pss_r6s_shoot.cir), to NOT be a safe blanket change for
 * every deck: an early version of this card's fix (gating only on "did the
 * period run fail", with no device-type check) rescued OTHERWISE-honest
 * convergence failures on the pss_r4_d55_unstable.cir red-proof fixture
 * (peak-current-mode instability, no `sw`/`csw` device anywhere in that
 * deck) into a spurious CONVERGED verdict -- a looser reltol or a
 * different perturbation size is not a side-effect-free "try harder", it
 * is a genuinely different (less accurate) trajectory, and for a marginal/
 * oscillating operating point that is enough to fool the residual test the
 * red-proof exists to exercise. Gating both fallbacks behind "this
 * circuit actually has a switch" makes them a provable, checked no-op for
 * every fixture in this project's regression suite that doesn't use one
 * (confirmed: test_pss.py back to its pre-card-#498 40/44+4info baseline
 * with this gate in place -- see the card #498 report) while still
 * reaching the deck family (R2, R6/R6o/R6s, pac_r6o_ctrl2out, ...) the
 * fallback was built for. Scanned once per PSSshootDriven() call (see its
 * own call site below) -- CSwitch (current-controlled) is intentionally
 * NOT included: this project's own switch-event work is SW-only so far
 * (patches-wip/salt/README.md's own "NOT reached" list), so treating a
 * CSW-only deck as switch-free is the conservative (no fallback) choice,
 * not an assumption that CSW is safe either way. */
static int
pss_ckt_has_switch(CKTcircuit *ckt)
{
    int swtype = CKTtypelook("Switch");
    SWmodel *sm;

    if (swtype < 0)
        return 0;
    for (sm = (SWmodel *) ckt->CKThead[swtype]; sm; sm = SWnextModel(sm))
        if (SWinstances(sm))
            return 1;
    return 0;
}

/* card #498 (plan: pss_r6s_shoot.cir closed-loop convergence): the deck's
 * own reltol (save_reltol in PSSshootDriven, e.g. the 1e-3 default), set
 * here once per PSSshootDriven() call, right after that value is known,
 * and ONLY if pss_ckt_has_switch() above says this deck can even reach the
 * pathology this fallback targets -- see pss_run_period_ex_safe() below.
 * 0.0 (the initial value, and the value shoot_cleanup restores it to, and
 * what a switch-free deck always leaves it at) means "no fallback
 * available", making pss_run_period_ex_safe() an exact pass-through -- so
 * a caller that never sets this (every OTHER analysis path,
 * PSSshootFreerun(), and every switch-free driven deck) is completely
 * unaffected. */
static double g_pss_loose_reltol = 0.0;
static int g_dbg_jf = -1; static double g_dbg_h = 0.0;   /* debug context for the PSS_SHOOT_DEBUG retry message */
/* card #498 2nd pass: set by pss_run_period_ex_safe() when the LAST run it
 * made needed the looser-reltol retry; a converged verdict must never rest
 * on such a run (the deck's own tolerance is the only honest one). */
static int g_pss_last_used_loose = 0;

/* card #498: companion gate for pss_fd_perturb_side()'s own shrinking-h
 * retry -- see pss_ckt_has_switch()'s header comment for why this needs
 * the SAME gate as g_pss_loose_reltol. Set alongside it. */
static int g_pss_has_switch = 0;

/* card #498: every Newton-loop period run (baseline, FD monodromy
 * perturbation columns, damped-step trial, fast-state detection, Floquet
 * re-evaluation, final trajectory emission) runs at CKTreltol tightened to
 * PSS_SHOOT_RELTOL/1e-5 (see PSSshootDriven's own comment on tight_reltol
 * for why: it is what lets the shooting engine resolve a switching edge
 * tightly enough for the monodromy to mean anything). On
 * tests/fixtures/pss_r6s_shoot.cir (card #498's own repro) this tight
 * value made the corrector need a sub-attosecond step exactly at the
 * "second", fast-ramp-fall saltation event patches-wip/salt/README.md's
 * FIFTH/SIXTH passes already characterise for this same deck family --
 * "Timestep too small ... trouble with swmod-instance shs" -- even though
 * a plain `.tran` of the identical circuit, and this same shooter's own
 * (always-loose) pre-stabilization warm-up transient, sail through the
 * identical instant every period. Confirmed by hand (card #498): rerunning
 * ONLY the failing period run at the deck's own looser reltol (no other
 * change) removes the failure outright. This wrapper automates exactly
 * that recovery, and ONLY that: the first attempt is always at whatever
 * CKTreltol the caller already set (byte-for-byte the pre-card-#498
 * behaviour), so a deck whose period runs already succeed there -- every
 * existing fixture -- never reaches the retry branch at all. Only a run
 * that already FAILED gets a second, looser-tolerance attempt before the
 * failure is reported upward; CKTreltol is restored to the caller's value
 * either way (this must never leak a looser tolerance into a later,
 * unrelated period run). */
static int
pss_run_period_ex_safe(CKTcircuit *ckt, double T, const char *jobname,
                        int record, int force)
{
    int err = pss_run_period_ex(ckt, T, jobname, record, force);

    g_pss_last_used_loose = 0;
    if (err && g_pss_loose_reltol > 0.0 && ckt->CKTreltol < g_pss_loose_reltol) {
        double saved_reltol = ckt->CKTreltol;
        if (getenv("PSS_SHOOT_DEBUG"))
            fprintf(stderr, "DEBUG %s failed at reltol=%.3g (ngspice error "
                    "%d); retrying once at the deck's own reltol=%.3g "
                    "(card #498 guarded fallback) before giving up. [col jf=%d h=%.6g]\n",
                    jobname, saved_reltol, err, g_pss_loose_reltol, g_dbg_jf, g_dbg_h);
        ckt->CKTreltol = g_pss_loose_reltol;
        err = pss_run_period_ex(ckt, T, jobname, record, force);
        ckt->CKTreltol = saved_reltol;
        g_pss_last_used_loose = (err == 0);
    }
    return err;
}


/* card #498: run ONE side (sign_h = +1.0 or -1.0) of a monodromy FD
 * perturbation column, retrying at a shrinking step if the period run
 * fails outright (ngspice error, e.g. the same swmod-instance "Timestep
 * too small" failure pss_run_period_ex_safe()'s own header comment
 * describes, but surfacing INSIDE a perturbed trajectory instead of the
 * baseline -- found on pss_r6s_shoot.cir, card #498, once the baseline-run
 * failure above was fixed: a large-enough perturbation can add or remove
 * an entire switching event relative to the baseline's forced grid,
 * corrupting or outright breaking that one column). Up to
 * PSS_FD_RETRY_MAX halvings of h (default 5; PSS_FD_RETRY_MAX overrides,
 * no effect unless set, same pattern as this file's other tunables)
 * before giving up on this side. Exact no-op (one call, no retry) for
 * every column that already succeeds at the first, unshrunk h -- every
 * existing fixture's own FD columns, confirmed by the regression suite
 * (card #498 report). Returns 0 and fills outbuf/*h_used on success (the
 * LAST h tried, which is h_start on the very first attempt); returns the
 * failing ngspice error code, unchanged outbuf, otherwise. */
static int
pss_fd_perturb_side(CKTcircuit *ckt, PSSstateMap *map, double T,
                     const double *x_full, int n_full, int jf, double x0j,
                     double sign_h, double h_start, const char *jobname,
                     double *outbuf, double *h_used)
{
    static int retry_max_base = -1;
    double h = h_start;
    double *icbuf;
    int err = 1, tries;   /* retry_max >= 1 guarantees this is overwritten below */
    int retry_max;

    if (retry_max_base < 0) {
        retry_max_base = getenv("PSS_FD_RETRY_MAX") ? atoi(getenv("PSS_FD_RETRY_MAX")) : 5;
        if (retry_max_base < 1) retry_max_base = 1;
    }
    /* card #498: see pss_ckt_has_switch()'s header comment -- a switch-
     * free deck gets exactly one attempt, byte-for-byte the pre-card-#498
     * behaviour (no shrinking, no chance to mask a genuine non-convergence
     * behind a smaller perturbation). */
    retry_max = g_pss_has_switch ? retry_max_base : 1;

    icbuf = TMALLOC(double, n_full);
    if (!icbuf)
        return E_NOMEM;

    for (tries = 0; tries < retry_max; tries++) {
        int i;
        for (i = 0; i < n_full; i++) icbuf[i] = x_full[i];
        icbuf[jf] = x0j + sign_h * h;
        g_dbg_jf = jf; g_dbg_h = sign_h * h;
        pss_state_set_ic(map, icbuf);
        err = pss_run_period_ex_safe(ckt, T, jobname, 0, 1);
        if (!err) {
            pss_state_get(map, ckt, outbuf);
            *h_used = h;
            tfree(icbuf);
            return 0;
        }
        if (getenv("PSS_SHOOT_DEBUG"))
            fprintf(stderr, "DEBUG %s h=%.6g failed (ngspice error %d); "
                    "shrinking (try %d/%d).\n", jobname, h, err, tries + 1,
                    retry_max);
        h *= 0.25;
    }
    tfree(icbuf);
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

/* card #494 (plan T6.1): saltation (event) correction for ideal `sw`
 * voltage-controlled switches. See research/simplis-pop-ac.md S4/S6 and
 * research/ngspice-internals-pss-pac.md S4.5 for the theory this
 * implements; pssshoot.h's own comment above the PSSsaltEvent struct has
 * the one-paragraph summary.
 *
 * THIRD pass (card #494) rewrite -- see docs/pss-pac.md's "THIRD pass"
 * section for the full writeup. The first two passes routed the
 * correction through an ag0-scaled companion-model current injected into
 * the event step's own G+ag0*C solve (rhs[node] -= ag0*dG*v_switch*coeff,
 * solved through SMPsolve()) -- entangling the correction with whichever
 * tiny, LTE-limited step size the event happened to land on (observed:
 * ~1.1ps, ag0~9e11), and never resolved past a reproducible ~2.2x error.
 * This pass computes the STATE-SPACE jump directly instead and applies it
 * as a value jump to the propagated column BEFORE the normal per-step
 * recursion runs, per the parent's suggested approach (card #494 third-
 * pass brief) and Fang's formula written in state-jump form:
 *
 *       delta_x(t_e+) = delta_x(t_e-) - (f+ - f-) * dt_e
 *       dt_e = -(grad(h) . delta_x(t_e-)) / h_dot
 *
 *   (grad(h).f- folded into h_dot exactly as the first two passes already
 *   argued: on ngspice's own accepted-step grid, h_dot -- the finite
 *   difference of v_ctrl across the step, from SWctrlvalue's own
 *   CKTstate0/1 history SWload.c already maintains for its OWN hysteresis
 *   logic, reused not recomputed -- IS the total derivative along the
 *   trajectory, so there is no separate grad(h).f- term to add).
 *
 *   grad(h) is nonzero only at the switch's own two control rows
 *   (posCntrl, negCntrl) -- an ideal switch's control port has no
 *   electrical coupling to its power port -- so
 *   grad(h).delta_x(t_e-) = y_prev[posCntrl] - y_prev[negCntrl], with
 *   y_prev = g_sens_Sprev / g_pac_Cprev, this column's own propagated
 *   sensitivity from the END of the previous accepted step (exactly
 *   delta_x(t_e-)).
 *
 *   f+ - f-: KEY CHANGE this pass. f+/f- are the reduced-state (capacitor
 *   voltage / inductor current) vector fields dx/dt just after/before the
 *   event, evaluated AT THE SAME state x(t_e) with the two different
 *   (pre-/post-event) switch topologies. Only the switch's own current
 *   injection changes between the two topologies (v_switch is a state,
 *   continuous across the event by construction -- only its time
 *   derivative jumps); every OTHER branch current touching the switch's
 *   power nodes (any resistor to a node whose own voltage is also
 *   continuous, any inductor -- inductor current cannot jump, it would
 *   need an impulsive voltage) is therefore ALSO continuous across the
 *   event and carries NO jump. By KCL at each power node, the entire
 *   current-injection jump dI = dG*v_switch (dG = g_after-g_before, the
 *   instance's actual conductance change) must then be absorbed by that
 *   node's OWN directly-attached capacitance alone:
 *
 *       f+_node - f-_node = +-dI / C_local(node)
 *
 *   C_local(node) = the sum of CAPcapac over every native capacitor
 *   instance with a terminal AT that node (pss_salt_local_cap() below;
 *   0 if none is directly there, e.g. an ideal voltage source's own
 *   node -- correctly no jump: an exogenous source's own voltage is not
 *   a state that can jump). This is a genuine Farad capacitance, with NO
 *   ag0 and NO dependence on the accepted step size anywhere in it --
 *   the step-size-independence the card's validation step (1) asks for.
 *   (Caveat, noted honestly: this local-KCL shortcut assumes the switch's
 *   OTHER-end current paths are not themselves jumping at the same
 *   instant -- true for every R1-R6 fixture checked, since only one
 *   node's current is discontinuous per event; a circuit where a
 *   capacitor bridges two simultaneously-switching nodes would need the
 *   fuller reduced-Jacobian solve the card's brief also offered as an
 *   option, not implemented this pass.)
 *
 *   So the whole correction, computed from y_prev = g_sens_Sprev /
 *   g_pac_Cprev (the UNMODIFIED value from the end of the PREVIOUS
 *   accepted step, delta_x(t_e-)) but ADDED DIRECTLY onto y_new =
 *   g_sens_Snew / g_pac_Cnew (this step's already-solved result), AFTER
 *   the caller's existing dceq/dveq history loop + SMPsolve for column
 *   j/m -- NOT routed through that solve:
 *
 *       coeff = (y_prev[posCntrl] - y_prev[negCntrl]) / h_dot   (= -dt_e)
 *       y_new[posNode] += dfdt_pos * coeff
 *       y_new[negNode] += dfdt_neg * coeff
 *
 *   with dfdt_pos = -dI/C_local(posNode), dfdt_neg = +dI/C_local(negNode)
 *   (same sign convention the first two passes' RHS injection used, dI
 *   leaving posNode/entering negNode). An EARLIER version of this pass
 *   applied the jump to y_prev BEFORE the solve instead (letting the
 *   existing per-step recursion carry it through this step's own G+ag0*C
 *   matrix) -- found, by direct measurement on pac_r6o_ctrl2out.cir, to
 *   still be wrong (7.73 V/V vs the ~6 V/V references, down from 13.19
 *   but not fixed): that matrix divides the jump by (ag0*C_local +
 *   G_new_switch + ...) instead of leaving it as the pure algebraic state
 *   discontinuity Fang's formula calls for, whenever the just-flipped
 *   switch's own new conductance is comparable to ag0*C_local at that
 *   node -- exactly the "not negligible" case the SECOND pass's report
 *   already flagged, now root-caused. Adding the jump directly to the
 *   ALREADY-SOLVED y_new avoids this: the smooth part of the step
 *   propagates the OLD (pre-jump) history through the NEW topology
 *   exactly as any other accepted step would (negligible, since the step
 *   itself is tiny), and the pure jump is superposed on top, undiluted.
 *   ag0 appears only in that surrounding, already-correct smooth
 *   propagation, never in the jump term itself. This is real-valued in
 *   pssshoot.c (one coeff per monodromy column j) and complex-valued in
 *   pacshoot.c (one coeff per PAC column m, y_prev/y_new complex,
 *   h_dot/dfdt still real) -- the CALLER does the scalar multiply/state-
 *   jump in its own (real or complex) arithmetic; this file only LOCATES
 *   the event and hands back (posNode, negNode, posCntrl, negCntrl,
 *   v_switch, dG, h_dot, dfdt_pos, dfdt_neg).
 *
 *   Location precision: unchanged from the first two passes -- this event
 *   is detected AFTER the fact, from whichever accepted step SWload.c's
 *   own state machine flipped SWswitchstate on (no breakpoint-based exact
 *   crossing time; see the card #494 report's "what remains"). The error
 *   this introduces is O(h) in the step size, same order the surrounding
 *   discretise-then-linearise recursion already carries for a smooth
 *   device's own G(t)/C(t) -- NOT the same thing as the ag0/step-size
 *   CONTAMINATION this pass removed from the jump formula itself.
 *
 *   Restart filtering (new this pass): PSSsaltReset(), called by both
 *   pss_sens_reset() and pacshoot.c's pac_reset() at the start of every
 *   fresh period/recording run, makes the NEXT PSSsaltFindEvents() call
 *   unconditionally return 0 -- filters the spurious ~1e-10s event both
 *   prior passes found (CKTstate1 not carrying the converged orbit's own
 *   wrap-around history on the first accepted step after a restart) and
 *   left unfixed.
 *
 *   Scope: SW (voltage-controlled switch) only. CSW (current-controlled)
 *   is NOT implemented (see the report) -- its control "voltage" is a
 *   branch current through another device, a different CKTstate slot
 *   convention this pass did not reach. */

/* card #494 FIFTH pass (plan T6.1): two fixes, per the card's own
 * diagnosis (docs/pss-pac.md's "FOURTH pass" section -- the jump formula
 * above is correct; two SEPARATE bugs remained in how it was located and
 * propagated).
 *
 * (1) EXACT h_dot (fixes the event-LOCATION-precision gap the third/
 * fourth passes flagged but did not root-cause). h_dot's finite
 * difference (vctrl_now-vctrl_prev) is exact -- CKTstate0/CKTstate1 are
 * genuine, exact samples of v_ctrl at the current and previous accepted
 * points -- the bug was in the DENOMINATOR: this function used
 * ckt->CKTdelta, which by the time dctran.c calls this hook (right after
 * CKTaccept(), dctran.c line ~502/546) has ALREADY been overwritten with
 * CKTtrunc()'s suggestion for the NEXT trial step's size (dctran.c line
 * ~993, `ckt->CKTdelta = newdelta;`, executed after the just-accepted
 * step's own NIiter()/CKTtrunc() pass, before control returns to
 * `nextTime:`/CKTaccept()) -- NOT the size of the step that actually
 * produced the CKTstate0 vs CKTstate1 transition being differenced. This
 * is the EXACT SAME class of bug as the already-documented CKTorder trap
 * a few dozen lines above in PSSsensRecordStep() (card #458 fourth pass)
 * -- a "read after accept but already advanced for the next attempt"
 * variable -- just never previously traced through for CKTdelta. Measured
 * directly (PSS_SALT_TRACE): at the real duty-cycle event, CKTdelta's
 * post-CKTtrunc value is ~1.25-2x the step that was actually taken,
 * matching the third pass's own reltol table almost exactly (825443 V/s
 * measured vs 1030928 V/s exact = 1.249x; 521443 V/s vs exact = 1.977x --
 * both consistent with CKTtrunc's own step-growth heuristic, which
 * commonly doubles the step after a tiny LTE-forced one). Fix: track this
 * function's own history of ckt->CKTtime across calls (g_salt_prev_time,
 * reset by PSSsaltReset()) and use the ACTUAL elapsed time between
 * consecutive accepted points as h_dot's denominator instead of
 * ckt->CKTdelta -- self-contained (no dependence on dctran.c's internal
 * CKTdeltaOld[]/CKTorder bookkeeping, whose exact semantics vary across
 * the XSPICE/non-XSPICE #ifdef paths), and exactly what pacshoot.c's
 * PSSpacRecordStep() already computes for its OWN, unrelated purposes
 * (`h = t_new - g_pac_tprev`) -- corroborating evidence this is the right
 * quantity. On `pac_r6o_ctrl2out.cir`'s exactly-piecewise-linear ramp
 * segment this makes h_dot match the analytic slope (1.030928e6 V/s) to
 * well under the card's 0.1% target (limited only by double-precision
 * roundoff on values ~0.5V differing by ~1e-6V over a ~1e-12s bracket --
 * no cancellation problem at this magnitude), independent of reltol/step
 * size -- see the card #494 fifth-pass report for the measured numbers.
 * f+/f- themselves were NOT the problem (dfdt_pos/neg depend only on Vin,
 * dG, C_local -- all step/time-independent constants; v_switch is read at
 * the already very precisely located accepted point) so no interpolation
 * or CKTsetBreak-based bisection of the crossing itself was needed once
 * this denominator bug was fixed -- see the report for why this was
 * judged sufficient rather than also adding breakpoint machinery.
 *
 * (2) L-stable propagation after an event (fixes the STIFFNESS bug the
 * fourth pass's independent arithmetic check found is a SECOND, separate
 * defect on top of (1) -- fixing h_dot alone left the code's own
 * (analytically re-checked) result at 7.5-11.9 V/V, not the ~6 V/V
 * references). Node `sw`'s local capacitance against the switch's own
 * ON/OFF conductance (and any other resistor directly on that node) gives
 * it a local time constant tau = C_local/G_local many orders of magnitude
 * below the accepted step sizes the surrounding recursion's LTE control
 * (tuned to the SLOW L/C states) actually picks. The existing recursion's
 * TRAPEZOIDAL formula (order 2, `-ag1*Tp` term) is A-stable but not
 * L-stable: its amplification factor -> -1 (not 0) as |lambda*h| -> Inf,
 * so it does not damp the huge post-jump spike this file injects at that
 * node -- it rings, and the ringing leaks into the slow states through
 * the A matrix's own v1<->iL coupling (see `patches-wip/salt/
 * check_r6o_arithmetic.py`'s A matrix). Per the card's option (a) (a
 * fuller reduced-Jacobian quasi-static elimination, option (b), was
 * considered and rejected here -- see the report: it would replace Fang's
 * EXACT raw jump, which check_r6o_arithmetic.py already proved gives the
 * right answer (6.017 V/V) under EXACT propagation, with a different,
 * approximate model; the bug is in this file's DISCRETE propagation of an
 * already-correct jump, not in the jump target itself), this pass forces
 * the SAME "order 1 (backward Euler) right after a breakpoint" idiom
 * dctran.c's own XSPICE path already uses for the real trajectory
 * (`ckt->CKTorder = 1` after a breakpoint, dctran.c) onto THIS recursion
 * for a settle window sized from the event's own node: N_TAU time
 * constants of real elapsed time after the event (PSS_SALT_SETTLE_TAUS
 * below), independent of step count -- self-terminating and adaptive to
 * whatever tau the triggering node actually has, not a fixture-tuned
 * magic step count. Backward Euler's amplification factor 1/(1-lambda*h)
 * -> 0 as |lambda*h| -> Inf regardless of h's size, so even the first
 * forced step crushes most of the spike; the multi-tau window covers the
 * handful of still-small (still not fully L-stable-safe under order 2)
 * steps the real LTE controller takes while growing back toward the slow
 * states' own step scale. This changes ONLY this file's own S/T (pssshoot.c)
 * and C/T (pacshoot.c) recursion -- never ckt->CKTorder itself, so the
 * real x(t) trajectory these recursions read their step/order history
 * FROM is completely unaffected (an existing-fixture, no-sw-event deck
 * sees PSSsaltForceOrder1() return 0 on every call -- exact no-op,
 * confirmed by regression). */

#define PSS_SALT_SETTLE_TAUS 30.0   /* generous safety margin -- see report */

int PSSsaltActive = 1;   /* PSS_SALTATION=0 env var disables (red-proof toggle) */

static int g_salt_first_step_pending = 0;   /* card #494 third pass: restart filter */
static double g_salt_prev_time = -1.0;      /* card #494 fifth pass: h_dot's true-elapsed-time tracker */
static double g_salt_settle_until = -1.0;   /* card #494 fifth pass: force-order-1 window deadline */
static int g_salt_events_seen = 0;          /* card #494 sixth pass: see PSSsaltEventsSeen() */
static int g_salt_ind_warned = 0;           /* card #510 step 1: see pss_salt_ind_jump() */
/* card #535: terminal voltage of every sw instance at the PREVIOUS accepted
 * point (index = instance order). The saltation current jump dI = dG*v_switch
 * needs the switch voltage at the crossing, i.e. the continuous state just
 * BEFORE the flip. With exact event location (patch 0065) the flip is
 * applied in the step that FOLLOWS the located crossing, and that step is
 * deliberately short but not infinitesimal (30 ps after patch 0070), so the
 * end-of-step CKTrhsOld voltage has already slewed through the switch node's
 * ron*Csw transient (R6o: 4.8 V instead of 11.5 V at the first event, 7.2 V
 * instead of -0.007 V at the second). */
static double *g_salt_vsw_prev = NULL, *g_salt_vsw_cur = NULL;
static int g_salt_vsw_n = 0, g_salt_vsw_have_prev = 0;

void
PSSsaltReset(void)
{
    g_salt_first_step_pending = 1;
    g_salt_prev_time = -1.0;
    g_salt_vsw_have_prev = 0;
    g_salt_settle_until = -1.0;
    g_salt_events_seen = 0;
    g_salt_ind_warned = 0;   /* card #510 step 1: once per analysis (fresh period/recording run) */
}

int
PSSsaltEventsSeen(void)
{
    return g_salt_events_seen;
}

/* card #494 fifth pass: is ckt->CKTtime still inside the post-event
 * L-stable settle window opened by the most recent PSSsaltFindEvents()
 * call that found a real event? See the derivation comment above (task
 * 2). Called by PSSsensRecordStep()/PSSpacRecordStep() right after their
 * own ag1-based order inference to override it to 1 for the duration. */
int
PSSsaltForceOrder1(CKTcircuit *ckt)
{
    return (g_salt_settle_until > 0.0 && ckt->CKTtime < g_salt_settle_until);
}

static int
pss_salt_is_on(double state)
{
    return state == 1.0 /* REALLY_ON */ || state == 3.0 /* HYST_ON */;
}

/* card #494 third pass: total capacitance (Farads) directly attached to
 * `node` (either terminal of any native capacitor instance) -- see the
 * derivation comment above for why this is exactly the right denominator
 * for that node's own (f+-f-) state jump. 0 if node<=0 (ground/unused) or
 * no capacitor is directly there. */
static double
pss_salt_local_cap(CKTcircuit *ckt, int node)
{
    int captype;
    CAPmodel *cm;
    CAPinstance *ci;
    double c = 0.0;

    if (node <= 0)
        return 0.0;

    captype = CKTtypelook("Capacitor");
    if (captype < 0)
        return 0.0;

    for (cm = (CAPmodel *) ckt->CKThead[captype]; cm; cm = CAPnextModel(cm))
        for (ci = CAPinstances(cm); ci; ci = CAPnextInstance(ci))
            if (ci->CAPposNode == node || ci->CAPnegNode == node)
                c += ci->CAPcapac;
    return c;
}

/* card #494 fifth pass: total conductance (Siemens) directly attached to
 * `node` -- every native Resistor instance with a terminal there, plus
 * every Switch instance's OWN current conductance (SWcond -- SWload()
 * already reloaded it for this accepted point before this hook runs,
 * same call site) with a terminal there. Paired with pss_salt_local_cap()
 * to estimate that node's own local RC time constant (tau=C_local/G_local)
 * for PSSsaltForceOrder1()'s settle window (task 2's L-stable-propagation
 * fix) -- NOT used anywhere in the jump formula itself, which stays
 * exactly the pure Farads-based KCL quantity the third/fourth passes
 * derived. Scope: R and SW only (this project's own R1-R6 fixture
 * family) -- a node whose dominant local conductance comes from another
 * device type (a diode, a MOSFET's rds, ...) would need extending this,
 * not attempted -- see the report. 0 if node<=0 or nothing conductive is
 * directly there (a pure capacitor-to-nothing node -- correctly means no
 * settle window is needed/possible from this estimate; PSSsaltForceOrder1()
 * just stays inactive for that event, unchanged from before this pass). */
static double
pss_salt_local_cond(CKTcircuit *ckt, int node)
{
    int restype, swtype;
    RESmodel *rm;
    RESinstance *ri;
    SWmodel *sm;
    SWinstance *si;
    double g = 0.0;

    if (node <= 0)
        return 0.0;

    restype = CKTtypelook("Resistor");
    if (restype >= 0)
        for (rm = (RESmodel *) ckt->CKThead[restype]; rm; rm = RESnextModel(rm))
            for (ri = RESinstances(rm); ri; ri = RESnextInstance(ri))
                if (ri->RESposNode == node || ri->RESnegNode == node)
                    g += ri->RESconduct;

    swtype = CKTtypelook("Switch");
    if (swtype >= 0)
        for (sm = (SWmodel *) ckt->CKThead[swtype]; sm; sm = SWnextModel(sm))
            for (si = SWinstances(sm); si; si = SWnextInstance(si))
                if (si->SWposNode == node || si->SWnegNode == node)
                    g += si->SWcond;

    return g;
}

/* card #510 (silently-wrong-result bug, high priority): a switch power
 * node whose AC/sensitivity perturbation reaches it, with NO local
 * capacitor there (pss_salt_local_cap()==0, e.g. a hard-switched buck's
 * own switching node feeding straight into an inductor, no snubber cap)
 * was silently giving dfdt_pos/dfdt_neg=0 -- the THIRD pass's comment
 * calls this "correctly no jump: an exogenous source's own voltage is not
 * a state that can jump", true ONLY for a node that is itself an
 * independent source's own terminal. For a node that is instead the near
 * terminal of an INDUCTOR (tests/fixtures/pac_sw_ctrl_direct.cir: node
 * `sw`, S1's negNode, feeding L1 with no capacitor anywhere on `sw`), it
 * is WRONG: `.pac`/`.pxf` measured exactly 0+0j control-to-output at
 * every frequency with no warning, while a `.tran` shows ~13 V/V DC gain
 * and the saltation machinery's own fs/2-band warning fires (proving it
 * DOES see the event) -- confirmed via PSS_SALT_TRACE=1
 * (dfdt_pos=0 dfdt_neg=0 at the real duty-cycle event) and via two
 * bracketing experiments (docs/pss-pac.md's card #510 section): (a) an
 * RC softening the CONTROL path (Vctrl -> switch control) does NOT change
 * the zero result at all, at any RC -- ruling out the control-path
 * hypothesis the card's brief floated ("the AC source enters the event
 * function directly, no state in between" -- refuted: the state-driven
 * coeff/dtau term, y_prev[posCntrl]-y_prev[negCntrl], already reaches a
 * directly-forced control node correctly, via the SAME full-circuit
 * per-step SMPsolve() that handles it for R6o); (b) adding a capacitor
 * directly AT the switch's POWER node (`sw`) DOES make the result
 * nonzero, confirming the actual missing piece is the (f+-f-) JUMP
 * TARGET, not the event-time-shift coeff.
 *
 * Root cause: PSSsensRecordStep()/PSSpacRecordStep()'s smooth per-step
 * recursion only ever READS a state's history (g_sens_Sprev/g_pac_Cprev)
 * at STATE rows -- capacitor node pairs and inductor branch-current rows,
 * the `for (d=0; d<n; d++)` loop in both files. A jump written onto a
 * NON-state node's own row (e.g. `sw` with no local cap) is therefore
 * INERT for propagation: it is overwritten from scratch, correctly, by
 * the very next accepted step's own SMPsolve() (which recomputes every
 * algebraic node fresh from the, by then correctly-jumped, states) --
 * this is exactly WHY pss_salt_local_cap()'s "0 at an exogenous node"
 * case is harmless there (nothing forward-propagating was ever lost) but
 * is a genuine bug at an inductor-adjacent node with no cap: the ONE
 * state that DOES persist there (the inductor's own current) was never
 * given its jump at all.
 *
 * Fix (pss_salt_ind_jump() below): when posNode/negNode has no local
 * capacitor, look for a UNIQUE adjacent inductor and route the jump onto
 * ITS branch-current row instead, in inductor-current-derivative units
 * (f_L=(V(pos)-V(neg))/L). f_L+-f_L- = (V(node)+ - V(node)-)/L, and
 * V(node)+ is simply the ACTUAL, already-solved accepted-point value
 * (ckt->CKTrhsOld[node]) -- no approximation needed there at all. Only
 * V(node)- (the COUNTERFACTUAL "what would this node's voltage be at
 * this same instant, holding the ADJACENT INDUCTOR's current I_L at its
 * (continuous, unable to jump) actual value, if the triggering switch
 * had NOT flipped") needs deriving, from a LOCAL one-node KCL solve at
 * `node`: the triggering switch (conductance g_prev, to the verified-
 * stable `otherNode`) plus every OTHER device actually touching `node`
 * (pss_salt_ind_jump() below walks R/other-SW/diode, same scope as
 * pss_salt_local_cond()) plus the fixed I_L, using PRE-EVENT (CKTstate1,
 * NOT CKTstate0) small-signal Norton equivalents for every nonlinear one
 * of those -- full derivation (including the two diode sign cases,
 * cathode vs RS=0 anode) in docs/pss-pac.md's card #510 section.
 *
 *   V(node)- = [g_prev*V(otherNode) + K_other_prev - s*I_L] /
 *              (g_prev + G_other_prev)
 *
 * (s = +1 if `node` is the found inductor's OWN posNode, else -1;
 * K_other_prev/G_other_prev = the sum of every other device's own
 * PRE-event Norton current/conductance at `node`, from its CKTstate1
 * history). No ag0, no step size, anywhere -- same step-size-
 * independence property the capacitor case has.
 *
 * EARLIER VERSION OF THIS FIX (this pass, superseded, kept here for the
 * record -- see docs/pss-pac.md's card #510 section, "a second bug
 * found"): a simpler `dV=dI/(g_prev+G_other)` shortcut, using ONLY the
 * POST-event G_other, gave a WILDLY wrong (many orders of magnitude)
 * answer on this very fixture's turn-ON edge, where the freewheeling
 * diode at the same node is ITSELF effectively co-switching (conducting
 * heavily just before the edge, off just after) -- that shortcut's
 * implicit assumption ("every other device at the node has a FIXED
 * conductance across the event") silently fails whenever the "other"
 * device is a diode whose own bias flips in sympathy with the switch, a
 * completely ordinary situation in a synchronous/freewheeling topology.
 * Caught by re-validating against tools/pac_perturbed.py rather than
 * trusting the first formula that removed the original all-zero bug --
 * see the report. The Norton/prev-state version above does not have this
 * failure mode: it solves the counterfactual PRE-event topology's own
 * node equation properly (using that topology's OWN diode operating
 * point, not a value borrowed from the post-event one).
 *
 * When neither a capacitor nor a usable unique adjacent inductor is
 * found, and the node is not itself an exogenous (source/ground/
 * capacitor-anchored) node either, this is refused outright (a stderr
 * message, PSSsaltFindEvents() returns -1) rather than silently
 * returning 0 -- per the card's own brief. Scope limits, stated
 * honestly: only ONE inductor per node; the switch's OTHER terminal, and
 * every "other" device's OWN far terminal, must be independently
 * verified stable (pss_salt_node_is_stable()) -- refused otherwise, a
 * genuinely coupled multi-node event this local solve cannot resolve;
 * the diode term only covers the cathode side and an RS=0 anode; a
 * NON-triggering switch's own PRE-event conductance is approximated by
 * its CURRENT SWcond (no separate CKTstate history kept for that field,
 * unlike the diode's) -- fine unless that other switch is ALSO flipping
 * at this exact same instant, not exercised by any fixture checked;
 * CSW, XSPICE, MOSFET/BJT conductances at such a node are NOT modelled
 * and fall through to the same refusal. */

/* card #510: is `node` guaranteed NOT to be jumping at this same event --
 * safe to hold fixed as a reference in pss_salt_ind_jump()'s local KCL?
 * True for ground/unused (node<=0), a node with its own local
 * capacitance (a real state, by construction slowly-varying across one
 * instantaneous event), or a node that is directly an independent
 * voltage source's own terminal (algebraically fixed by the source
 * regardless of current -- the same "exogenous node" pss_salt_local_
 * cap()'s own comment already argues can never jump). */
static int
pss_salt_node_is_stable(CKTcircuit *ckt, int node)
{
    int vtype;
    VSRCmodel *vm;
    VSRCinstance *vi;

    if (node <= 0)
        return 1;
    if (pss_salt_local_cap(ckt, node) > 0.0)
        return 1;

    vtype = CKTtypelook("Vsource");
    if (vtype >= 0)
        for (vm = (VSRCmodel *) ckt->CKThead[vtype]; vm; vm = VSRCnextModel(vm))
            for (vi = VSRCinstances(vm); vi; vi = VSRCnextInstance(vi))
                if (vi->VSRCposNode == node || vi->VSRCnegNode == node)
                    return 1;

    return 0;
}

/* card #510: the inductor-routed forcing term for `node` (one of a
 * switch's two power terminals) -- see the derivation comment above.
 * `otherNode` is the switch's OTHER power terminal (held fixed, must be
 * pss_salt_node_is_stable()). `g_prev` is the caller's already-computed
 * pre-event conductance of the TRIGGERING switch itself. On success
 * (return 1) fills *out_row (the target inductor's INDbrEq) and
 * *out_djump (the (f+-f-) value, in Amps/second, to be scaled by coeff
 * exactly like the capacitor case). Returns 0 if there is genuinely
 * nothing to do (no adjacent inductor -- a purely algebraic node, no
 * persistent state, no jump needed, same as an exogenous node). Returns
 * -1 (REFUSE) and fills *out_why if this event is outside what this
 * function can safely resolve. */
static int
pss_salt_ind_jump(CKTcircuit *ckt, SWinstance *sw, int node, int otherNode,
                   double g_prev, double dt_actual, int *out_row, double *out_djump,
                   double *out_tau, const char **out_why)
{
    int indtype, n_found = 0;
    INDmodel *im;
    INDinstance *ii, *found = NULL;
    int found_sign = 0;
    int restype, swtype, diotype;
    RESmodel *rm; RESinstance *ri;
    SWmodel *sm; SWinstance *si;
    DIOmodel *dm; DIOinstance *di;
    double Gsum, Ksum, Vfar;
    double V_now, V_other, I_L, V_prev, dV;

    if (node <= 0)
        return 0;

    indtype = CKTtypelook("Inductor");
    if (indtype < 0)
        return 0;

    for (im = (INDmodel *) ckt->CKThead[indtype]; im; im = INDnextModel(im))
        for (ii = INDinstances(im); ii; ii = INDnextInstance(ii)) {
            if (ii->INDposNode == node) {
                found = ii; found_sign = +1; n_found++;
            } else if (ii->INDnegNode == node) {
                found = ii; found_sign = -1; n_found++;
            }
        }

    if (n_found == 0)
        return 0;   /* purely algebraic node -- correctly nothing to do */
    if (n_found > 1) {
        *out_why = "more than one inductor terminal at the switching node "
                   "(ambiguous which one absorbs the jump)";
        return -1;
    }
    if (!pss_salt_node_is_stable(ckt, otherNode)) {
        *out_why = "the switch's other terminal is not itself a verified-"
                   "stable node (no local capacitor, not a source terminal) "
                   "-- a simultaneous double jump this engine cannot resolve";
        return -1;
    }

    Gsum = 0.0; Ksum = 0.0;

    restype = CKTtypelook("Resistor");
    if (restype >= 0)
        for (rm = (RESmodel *) ckt->CKThead[restype]; rm; rm = RESnextModel(rm))
            for (ri = RESinstances(rm); ri; ri = RESnextInstance(ri))
                if (ri->RESposNode == node || ri->RESnegNode == node) {
                    int far = (ri->RESposNode == node) ? ri->RESnegNode : ri->RESposNode;
                    if (!pss_salt_node_is_stable(ckt, far)) {
                        *out_why = "a resistor's own far terminal at the switching "
                                   "node is not itself a verified-stable node";
                        return -1;
                    }
                    Vfar = (far <= 0) ? 0.0 : ckt->CKTrhsOld[far];
                    Gsum += ri->RESconduct;
                    Ksum += ri->RESconduct * Vfar;
                }

    swtype = CKTtypelook("Switch");
    if (swtype >= 0)
        for (sm = (SWmodel *) ckt->CKThead[swtype]; sm; sm = SWnextModel(sm))
            for (si = SWinstances(sm); si; si = SWnextInstance(si))
                if (si != sw && (si->SWposNode == node || si->SWnegNode == node)) {
                    int far = (si->SWposNode == node) ? si->SWnegNode : si->SWposNode;
                    if (!pss_salt_node_is_stable(ckt, far)) {
                        *out_why = "another switch's own far terminal at the switching "
                                   "node is not itself a verified-stable node";
                        return -1;
                    }
                    Vfar = (far <= 0) ? 0.0 : ckt->CKTrhsOld[far];
                    /* approximation: this OTHER switch's OWN pre-event
                     * conductance is not separately tracked (unlike the
                     * diode's CKTstate1 below) -- its CURRENT SWcond
                     * stands in for both; see the header comment. */
                    Gsum += si->SWcond;
                    Ksum += si->SWcond * Vfar;
                }

    diotype = CKTtypelook("Diode");
    if (diotype >= 0)
        for (dm = (DIOmodel *) ckt->CKThead[diotype]; dm; dm = DIOnextModel(dm))
            for (di = DIOinstances(dm); di; di = DIOnextInstance(di)) {
                int is_cathode = (di->DIOnegNode == node);
                int is_anode = (di->DIOposNode == node && di->DIOposPrimeNode == di->DIOposNode);
                int far;
                double Gd, Vd, Id;

                if (!is_cathode && !is_anode)
                    continue;

                far = is_cathode ? di->DIOposNode : di->DIOnegNode;
                if (!pss_salt_node_is_stable(ckt, far)) {
                    *out_why = "a diode's own far terminal at the switching node "
                               "is not itself a verified-stable node";
                    return -1;
                }
                Vfar = (far <= 0) ? 0.0 : ckt->CKTrhsOld[far];
                Gd = ckt->CKTstate1[di->DIOconduct];
                Vd = ckt->CKTstate1[di->DIOvoltage];   /* anode-cathode, PREV point */
                Id = ckt->CKTstate1[di->DIOcurrent];   /* anode->cathode current, PREV point */

                Gsum += Gd;
                if (is_cathode)
                    /* current injected into the cathode = +I_D(V_D),
                     * V_D = Vfar(anode) - V(node); see the header comment. */
                    Ksum += Id + Gd * (Vfar - Vd);
                else
                    /* node is the anode; current injected = -I_D(V_D),
                     * V_D = V(node) - Vfar(cathode). */
                    Ksum += -Id + Gd * (Vfar + Vd);
            }

    if (g_prev + Gsum <= 0.0) {
        *out_why = "degenerate switch model (g_prev+G_other_prev <= 0)";
        return -1;
    }

    V_now = ckt->CKTrhsOld[node];
    V_other = (otherNode <= 0) ? 0.0 : ckt->CKTrhsOld[otherNode];
    I_L = ckt->CKTrhsOld[found->INDbrEq];

    V_prev = (g_prev * V_other + Ksum - (double) found_sign * I_L) / (g_prev + Gsum);

    dV = V_now - V_prev;
    *out_djump = (double) found_sign * dV / found->INDinduct;
    *out_row = found->INDbrEq;

    /* card #510 step 2 (accuracy pass, superseding the step-1 L/R-based
     * *out_tau this comment replaces): settle-window sizing for
     * PSSsaltForceOrder1() below (the FIFTH pass's L-stable-propagation
     * fix). Step 1 used the SAME "N_TAU time constants of real elapsed
     * time" idiom the capacitor case uses, with tau = L/(g_prev+Gsum) (an
     * L/R time constant) in place of C/G -- dimensionally the obvious
     * analogue, but WRONG in practice: measured worst-case error against
     * tools/pac_perturbed.py on pac_sw_ctrl_direct.cir was ~3.6dB/3.1deg
     * with that tau (PSS_SALT_SETTLE_TAUS=30 * ~0.5us = a ~15us window,
     * 7.5 SWITCHING PERIODS at this fixture's 2us period -- forcing
     * backward Euler across most of the recorded trajectory, not a brief
     * post-event settle). Re-derived: the FIFTH pass's own stiffness
     * argument (trapezoidal's amplification factor -> -1, not 0, as
     * |lambda*h|->Inf, so it rings on a mode STIFF relative to the
     * accepted step h) needs a mode with |lambda|*h large -- an
     * INDUCTOR's own L/R tau being LARGE (as it almost always is by
     * converter design -- that is WHY an inductor filters) means
     * |lambda|=1/tau is SMALL, so |lambda*h|<<1 for this state's OWN
     * eigenvalue: no genuine L-stability issue on the jumped state
     * itself, unlike R6o's tiny-Csw node (tau~50ps, genuinely stiff
     * against a ~10ns step). What the window still measurably helps with
     * here (scale=0 -- no window at all -- measured WORSE, ~6dB) is
     * something else: not decaying a stiff ring, but not reusing the
     * PRE-event order-2 history term (`Tp`) across the jump for the
     * first step or two -- the same "order 1 for exactly the step(s)
     * right after a breakpoint" idiom dctran.c's own real-trajectory path
     * already uses, sized in STEPS, not tau. Fix: size the window
     * directly off THIS event's own accepted-step size (`dt_actual`,
     * already computed by the caller -- self-scaling, no L/R/fixture-
     * specific constant to tune) instead of tau: a fixed small number of
     * step-widths, found by sweeping PSS_SALT_IND_STEPS against
     * tools/pac_perturbed.py (1k/10k/50k/100k/200kHz) -- 2 steps met the
     * 0.5dB/3deg target at every point (worst 0.48dB/2.37deg); see
     * docs/pss-pac.md's card #510 section for the full sweep table. Env
     * var kept as a documented, default-off-behaviour-preserving (default
     * matches the value validated below) diagnostic for any future
     * re-tuning, same convention as PSS_SALT_TRACE/PSS_SALT_NOFILTER. */
    {
        static int n_checked = 0;
        static double n_steps = 2.0;
        if (!n_checked) {
            const char *e = getenv("PSS_SALT_IND_STEPS");
            if (e) n_steps = atof(e);
            n_checked = 1;
        }
        *out_tau = (n_steps * dt_actual) / PSS_SALT_SETTLE_TAUS;
    }

    /* card #510 step 1 (maintainer decision: MERGE WITH A WARNING, not hold --
     * see the report): this path is genuinely less validated than the
     * THIRD pass's capacitor case (0.001dB on R6o) -- measured worst case
     * ~3.6dB/3.1deg vs tools/pac_perturbed.py on pac_sw_ctrl_direct.cir,
     * see docs/pss-pac.md's card #510 section -- so every analysis that
     * actually takes it must say so, once, unmissably, on stderr. Gated
     * by g_salt_ind_warned (cleared by PSSsaltReset(), so once per fresh
     * period/recording run == once per .pac/.pstb/.pss/.pxf analysis, the
     * same granularity patch 0904's VALIDATED-BAND GUARD warning uses) --
     * never fires for a deck whose events all resolve via the capacitor
     * case (R6o and every other existing fixture), confirmed by
     * tests/test_pac.py's new check_pac_ind_jump_warning(). */
    if (!g_salt_ind_warned) {
        g_salt_ind_warned = 1;
        fprintf(stderr, "PAC-TD: WARNING -- saltation via inductor-routed jump at node %s "
                "(no capacitor at the switch node): validated to 0.5 dB/3 deg vs an "
                "independent reference on ONE fixture only (card #510); cross-check with tools/pac_perturbed.py "
                "-- see docs/pss-pac.md's card #510 section.\n",
                NODENAME(ckt, node));
    }
    return 1;
}

int
PSSsaltFindEvents(CKTcircuit *ckt, PSSsaltEvent *ev, int max)
{
    int swtype;
    SWmodel *model;
    SWinstance *here;
    int n = 0;
    double dt_actual, tau_max;
    static int env_checked = 0;

    if (!env_checked) {
        const char *e = getenv("PSS_SALTATION");
        if (e && strcmp(e, "0") == 0)
            PSSsaltActive = 0;
        env_checked = 1;
    }
    if (!PSSsaltActive || ckt->CKTtime <= 0.0 || ckt->CKTdelta <= 0.0)
        return 0;

    /* card #494 fifth pass: this function's own true-elapsed-time
     * tracker for h_dot's denominator (task 1 -- see the derivation
     * comment above PSSsaltActive's definition for why ckt->CKTdelta
     * itself is NOT safe to use here). Updated UNCONDITIONALLY, before
     * the restart-filter's early return below, so g_salt_prev_time is
     * always seeded from the previous call's ckt->CKTtime regardless of
     * whether that call reported an event. g_salt_prev_time<0 (only true
     * on the very first call after PSSsaltReset()) falls back to
     * ckt->CKTdelta -- harmless, since the restart filter unconditionally
     * discards that first call's own event anyway. */
    dt_actual = (g_salt_prev_time >= 0.0) ? (ckt->CKTtime - g_salt_prev_time) : ckt->CKTdelta;
    g_salt_prev_time = ckt->CKTtime;

    /* card #535: roll the per-switch terminal-voltage snapshot (see above). */
    {
        int st_ = CKTtypelook("Switch"), cnt = 0;
        SWmodel *m_;
        SWinstance *h_;
        if (st_ >= 0) {
            for (m_ = (SWmodel *) ckt->CKThead[st_]; m_; m_ = SWnextModel(m_))
                for (h_ = SWinstances(m_); h_; h_ = SWnextInstance(h_))
                    cnt++;
            if (cnt != g_salt_vsw_n) {
                g_salt_vsw_prev = TREALLOC(double, g_salt_vsw_prev, cnt > 0 ? cnt : 1);
                g_salt_vsw_cur  = TREALLOC(double, g_salt_vsw_cur,  cnt > 0 ? cnt : 1);
                g_salt_vsw_n = cnt;
                g_salt_vsw_have_prev = 0;
            }
            if (g_salt_vsw_have_prev)
                memcpy(g_salt_vsw_prev, g_salt_vsw_cur, (size_t) cnt * sizeof(double));
            cnt = 0;
            for (m_ = (SWmodel *) ckt->CKThead[st_]; m_; m_ = SWnextModel(m_))
                for (h_ = SWinstances(m_); h_; h_ = SWnextInstance(h_)) {
                    g_salt_vsw_cur[cnt] = ckt->CKTrhsOld[h_->SWposNode] - ckt->CKTrhsOld[h_->SWnegNode];
                    if (!g_salt_vsw_have_prev)
                        g_salt_vsw_prev[cnt] = g_salt_vsw_cur[cnt];
                    cnt++;
                }
            g_salt_vsw_have_prev = 1;
        }
    }

    /* card #494 third pass: restart filter -- see PSSsaltReset()'s own
     * header comment in pssshoot.h. Consumes the latch unconditionally,
     * even if this particular step has no real sw event, so it only ever
     * suppresses the FIRST call after a reset. card #510 DIAGNOSTIC ONLY
     * (not the fix, opt-in, default-off, zero behaviour change unless set):
     * PSS_SALT_NOFILTER=1 skips this suppression, to test the hypothesis
     * that a REAL switching edge landing at/near a period boundary (e.g.
     * tests/fixtures/pac_sw_ctrl_direct.cir's Vramp, whose fast fall-back
     * edge lands essentially exactly at t=T) is being discarded by this
     * filter along with the stale-history artifact it was designed for --
     * see docs/pss-pac.md's card #510 section. NOT part of the shipped
     * fix (not validated, no regression run under it) -- reverting this
     * one-line change is `git revert`/dropping just this hunk. */
    if (g_salt_first_step_pending) {
        g_salt_first_step_pending = 0;
        if (!(getenv("PSS_SALT_NOFILTER") && strcmp(getenv("PSS_SALT_NOFILTER"), "1") == 0))
            return 0;
    }

    tau_max = 0.0;
    swtype = CKTtypelook("Switch");
    if (swtype < 0)
        return 0;

    int sw_idx = -1;
    for (model = (SWmodel *) ckt->CKThead[swtype]; model && n < max;
         model = SWnextModel(model)) {
        for (here = SWinstances(model); here && n < max;
             here = SWnextInstance(here)) {
            double now = ckt->CKTstate0[here->SWswitchstate];
            double prev = ckt->CKTstate1[here->SWswitchstate];
            double vctrl_now, vctrl_prev, h_dot, g_now, g_prev, dG, v_switch, dI;
            double Cpos, Cneg, Gpos, Gneg, tau_pos, tau_neg;
            int posRow, negRow, indrc;
            double dfdt_pos, dfdt_neg;
            double tau_ind_pos = 0.0, tau_ind_neg = 0.0;
            const char *why = NULL;

            sw_idx++;
            if (pss_salt_is_on(now) == pss_salt_is_on(prev))
                continue;   /* hysteresis-band-only or no-change: no real event */

            vctrl_now = ckt->CKTstate0[here->SWctrlvalue];
            vctrl_prev = ckt->CKTstate1[here->SWctrlvalue];
            h_dot = (vctrl_now - vctrl_prev) / dt_actual;
            if (fabs(h_dot) < 1e-30)
                continue;   /* degenerate: avoid a divide-by-zero downstream */

            /* card #535: a genuine event has the flip threshold between the
             * control values at the two ends of this step (plain mode: the
             * crossing is inside the step; exact-location mode: the located
             * crossing is at its start, to ~1e-16 s). A flip with the control
             * far from the threshold is not a crossing but a switch whose
             * restarted state disagrees with its control -- R6o's hi switch
             * restarts OFF at t=0 with its control 0.545 V past the ON
             * threshold. The first-step filter (g_salt_first_step_pending)
             * catches that in plain mode, but with the lagged decision
             * (patch 0065) the corrective flip lands one step later, got
             * through the filter and injected a spurious jump at t=0 (6.157
             * instead of 6.2036 V/V). Positive hysteresis only (the
             * negative-hysteresis semantics differ). */
            if (model->SWvHysteresis > 0.0 && !getenv("PSS_SALT_NOCROSSCHK")) {
                double thr = model->SWvThreshold +
                             (pss_salt_is_on(now) ? model->SWvHysteresis : -model->SWvHysteresis);
                double clo = (vctrl_now < vctrl_prev) ? vctrl_now : vctrl_prev;
                double chi = (vctrl_now < vctrl_prev) ? vctrl_prev : vctrl_now;
                double slack = 1e-6 * (fabs(thr) + 1e-3) + 0.01 * (chi - clo);
                if (thr < clo - slack || thr > chi + slack) {
                    if (getenv("PSS_SALT_TRACE"))
                        fprintf(stderr, "SALT-TRACE t=%.9g flip of %s ignored: control %.9g -> %.9g does not cross threshold %.9g\n",
                                ckt->CKTtime, here->SWname, vctrl_prev, vctrl_now, thr);
                    continue;
                }
            }

            g_now  = pss_salt_is_on(now)  ? model->SWonConduct : model->SWoffConduct;
            g_prev = pss_salt_is_on(prev) ? model->SWonConduct : model->SWoffConduct;
            dG = g_now - g_prev;
            /* card #535: the switch voltage BEFORE the flip (previous accepted point), not the
             * end-of-event-step value; PSS_SALT_VEND=1 restores the old end-of-step read. */
            v_switch = (sw_idx >= 0 && sw_idx < g_salt_vsw_n && !getenv("PSS_SALT_VEND"))
                       ? g_salt_vsw_prev[sw_idx]
                       : ckt->CKTrhsOld[here->SWposNode] - ckt->CKTrhsOld[here->SWnegNode];
            dI = dG * v_switch;

            Cpos = pss_salt_local_cap(ckt, here->SWposNode);
            Cneg = pss_salt_local_cap(ckt, here->SWnegNode);

            /* card #510: capacitor case unchanged (THIRD pass, validated);
             * otherwise try the inductor-routed forcing term; otherwise
             * (no cap, no unique adjacent inductor) posRow/negRow stay 0
             * -- correctly no jump, e.g. an exogenous source node or a
             * purely algebraic/resistive one, see the derivation comment
             * above pss_salt_ind_jump(). A REFUSE (-1) from
             * pss_salt_ind_jump() aborts the whole call immediately. */
            if (Cpos > 0.0) {
                posRow = here->SWposNode;
                dfdt_pos = -dI / Cpos;
            } else {
                posRow = 0;
                dfdt_pos = 0.0;
                indrc = pss_salt_ind_jump(ckt, here, here->SWposNode, here->SWnegNode,
                                           g_prev, dt_actual, &posRow, &dfdt_pos, &tau_ind_pos, &why);
                if (indrc < 0) {
                    fprintf(stderr, "PSS-SALT: REFUSING saltation correction for "
                            "switch %s at t=%.9g (posNode=%d): %s -- this .pss/.pac/"
                            ".pxf run cannot be trusted for this deck; aborting rather "
                            "than return a silently-wrong 0 (card #510).\n",
                            here->SWname, ckt->CKTtime, here->SWposNode, why);
                    return -1;
                }
            }

            if (Cneg > 0.0) {
                negRow = here->SWnegNode;
                dfdt_neg = dI / Cneg;
            } else {
                negRow = 0;
                dfdt_neg = 0.0;
                indrc = pss_salt_ind_jump(ckt, here, here->SWnegNode, here->SWposNode,
                                           g_prev, dt_actual, &negRow, &dfdt_neg, &tau_ind_neg, &why);
                if (indrc < 0) {
                    fprintf(stderr, "PSS-SALT: REFUSING saltation correction for "
                            "switch %s at t=%.9g (negNode=%d): %s -- this .pss/.pac/"
                            ".pxf run cannot be trusted for this deck; aborting rather "
                            "than return a silently-wrong 0 (card #510).\n",
                            here->SWname, ckt->CKTtime, here->SWnegNode, why);
                    return -1;
                }
            }

            ev[n].posNode = here->SWposNode;
            ev[n].negNode = here->SWnegNode;
            ev[n].posCntrl = here->SWposCntrlNode;
            ev[n].negCntrl = here->SWnegCntrlNode;
            ev[n].v_switch = v_switch;
            ev[n].dG = dG;
            ev[n].h_dot = h_dot;
            ev[n].posRow = posRow;
            ev[n].negRow = negRow;
            ev[n].dfdt_pos = dfdt_pos;
            ev[n].dfdt_neg = dfdt_neg;
            ev[n].name = here->SWname;

            /* card #494 fifth pass (task 2): estimate this event's own
             * settle-window requirement from its power nodes' local RC --
             * see pss_salt_local_cond()'s header and PSSsaltForceOrder1(). */
            Gpos = pss_salt_local_cond(ckt, here->SWposNode);
            Gneg = pss_salt_local_cond(ckt, here->SWnegNode);
            tau_pos = (Cpos > 0.0 && Gpos > 0.0) ? Cpos / Gpos : 0.0;
            tau_neg = (Cneg > 0.0 && Gneg > 0.0) ? Cneg / Gneg : 0.0;
            if (tau_pos > tau_max) tau_max = tau_pos;
            if (tau_neg > tau_max) tau_max = tau_neg;
            /* card #510: same settle-window estimate, inductor-routed case
             * (tau_ind_pos/neg, filled by pss_salt_ind_jump() above -- 0
             * when that path did not apply). */
            if (tau_ind_pos > tau_max) tau_max = tau_ind_pos;
            if (tau_ind_neg > tau_max) tau_max = tau_ind_neg;

            n++;
        }
    }
    if (n > 0)
        g_salt_events_seen += n;   /* card #494 sixth pass: see PSSsaltEventsSeen() */
    if (n > 0 && tau_max > 0.0) {
        double deadline = ckt->CKTtime + PSS_SALT_SETTLE_TAUS * tau_max;
        if (deadline > g_salt_settle_until)
            g_salt_settle_until = deadline;
    }
    if (n > 0 && getenv("PSS_SALT_TRACE"))
        fprintf(stderr, "SALT-TRACE t=%.9g dt_actual=%.6g (CKTdelta=%.4g) n_events=%d ev0=%s dG=%.6g "
                "v_switch=%.6g h_dot=%.6g posCntrl=%d negCntrl=%d dfdt_pos=%.6g dfdt_neg=%.6g "
                "tau_max=%.4g settle_until=%.6g\n",
                ckt->CKTtime, dt_actual, ckt->CKTdelta, n, ev[0].name, ev[0].dG, ev[0].v_switch, ev[0].h_dot,
                ev[0].posCntrl, ev[0].negCntrl, ev[0].dfdt_pos, ev[0].dfdt_neg, tau_max, g_salt_settle_until);
    return n;
}

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
    PSSsaltReset();   /* card #494 third pass: fresh run -- filter the first step's event */
}

/* dctran.c hook, called right after CKTaccept() (same site as
 * PSSbreakRecordStep) whenever PSSsensActive. See the long design comment
 * above for the recursion this implements. */
void
PSSsensRecordStep(CKTcircuit *ckt)
{
    int n, neq, d, j, nd, order, i_salt;
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

    {
        /* card #494: locate any sw event at THIS accepted step once, reuse
         * for every column j below (see PSSsaltFindEvents()'s header
         * comment in this file for the derivation). */
        PSSsaltEvent salt_ev[PSS_SALT_MAX_EVENTS];
        int n_salt = PSSsaltFindEvents(ckt, salt_ev, PSS_SALT_MAX_EVENTS);

        /* card #510: PSSsaltFindEvents() already printed why; a REFUSE
         * here must fail this whole (PSS-sens/monodromy) pass cleanly,
         * not silently continue with an unhandled event's jump missing. */
        if (n_salt < 0) {
            g_sens_failed = 1;
            return;
        }

        /* card #494 fifth pass (task 2): force backward Euler (L-stable)
         * for the whole step, all columns, while a post-event settle
         * window from a PREVIOUS event is still open -- see
         * PSSsaltForceOrder1()'s header comment in pssshoot.h. Must run
         * AFTER PSSsaltFindEvents() above (which may just have opened or
         * extended the window from an event at THIS very step) and BEFORE
         * `order` is used below. No-op (order unchanged) for every step
         * of every deck with no sw event ever seen this run. */
        if (PSSsaltForceOrder1(ckt))
            order = 1;

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

    /* card #494 THIRD pass (rewrite -- see PSSsaltFindEvents()'s derivation
     * comment above): apply the saltation correction as a DIRECT VALUE JUMP
     * added onto g_sens_Snew, AFTER this step's own smooth dceq/dveq solve
     * above, NOT routed through it. Routing the jump through this step's
     * own SMPsolve (tried earlier this pass, before this rewrite)
     * measurably attenuated it whenever the just-flipped switch's own new
     * conductance is comparable to ag0*C_local at that node (exactly the
     * ~9x-not-negligible case the SECOND pass's report already flagged) --
     * the matrix divides the jump by (ag0*C + G_new + ...) instead of
     * leaving it as the pure state discontinuity Fang's formula calls for.
     * Adding it directly to the already-solved g_sens_Snew avoids that
     * entirely. coeff uses g_sens_Sprev (the UNMODIFIED value from the END
     * of the PREVIOUS accepted step -- delta_x(t_e-), exactly what the
     * formula calls for) for every event, so simultaneous events (e.g. two
     * complementary switches) never see a partially-applied correction. */
    if (n_salt > 0) {
        for (j = 0; j < n; j++) {
            double djump[PSS_SALT_MAX_EVENTS][2];   /* [i_salt][0]=posRow, [1]=negRow */
            for (i_salt = 0; i_salt < n_salt; i_salt++) {
                PSSsaltEvent *e = &salt_ev[i_salt];
                double Sp = pss_sens_get(g_sens_Sprev, e->posCntrl, n, j);
                double Sn = pss_sens_get(g_sens_Sprev, e->negCntrl, n, j);
                double coeff = (Sp - Sn) / e->h_dot;
                djump[i_salt][0] = e->dfdt_pos * coeff;
                djump[i_salt][1] = e->dfdt_neg * coeff;
            }
            for (i_salt = 0; i_salt < n_salt; i_salt++) {
                PSSsaltEvent *e = &salt_ev[i_salt];
                /* card #510: posRow/negRow (posNode/negNode for the
                 * THIRD pass's original capacitor case, an adjacent
                 * inductor's own INDbrEq row for the new inductor-routed
                 * case, or 0/skip) -- never posNode/negNode directly, see
                 * pssshoot.h's PSSsaltEvent comment for why a non-state
                 * node's own row is inert. */
                if (e->posRow > 0)
                    pss_sens_set(g_sens_Snew, e->posRow, n, j,
                                 pss_sens_get(g_sens_Snew, e->posRow, n, j) + djump[i_salt][0]);
                if (e->negRow > 0)
                    pss_sens_set(g_sens_Snew, e->negRow, n, j,
                                 pss_sens_get(g_sens_Snew, e->negRow, n, j) + djump[i_salt][1]);
            }
        }
    }
    }   /* card #494: closes the salt_ev/n_salt scope opened above */

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


/* ==================================================================== */
/* card #537: matrix-free Newton-GMRES shooting (Telichevesky, Kundert &
 * White, DAC 1995) -- opt-in `krylov=1` on the .pss/.pac/.pstb/.pxf/.pnoise
 * card or NGSPICE_PSS_KRYLOV=1 in the environment.
 *
 * Why: the dense engine builds the monodromy M one column per retained
 * state, 2 perturbed period runs per column (centred difference), every
 * Newton iteration -- and, on a converter with extracted parasitics
 * (R3 buck + 2418 R/C = ~1150 capacitor states), pss_detect_fast_states()'s
 * one-time full sweep of 2*1150+ period runs is what the fast-state
 * reduction costs to learn that only ~5 states matter. Parasitic states
 * decay inside the period (multiplier ~0), so M has only a handful of
 * non-negligible eigenvalues and a Krylov method needs only about
 * (slow modes + a few) products M*v, each ONE period run pair from the
 * current iterate. Nothing of M is assembled and no fast-state detection
 * runs; every state is a Newton unknown.
 *
 * What is shared with the dense engine (so F and the products are one map):
 * the baseline P(x0) block of pss_newton_run (recorded grid, forced rerun
 * for switch decks), pss_run_period_ex_safe(..., force=1) for every product
 * (same recorded grid / replay as the FD columns and the damped-step
 * trials), the ripple-relative scale[] that defines the residual norm, the
 * trust-radius cap, the damped line search, the degenerate-fixed-point
 * guard, the staged warm-up and the replay ladder (all unchanged: only the
 * linear solve for the Newton direction is replaced).
 *
 * Linear solve: inexact Newton, (M - I) dx = -F by GMRES in the
 * scale[]-weighted variables y = x/scale (the same weights the residual
 * test max_k |F_k|/scale_k uses), x0 = 0, modified Gram-Schmidt with one
 * re-orthogonalisation, Givens least squares. Stop at relative residual eta
 * (PSS_KRYLOV_ETA, default 1e-3) or at the Krylov dimension limit
 * (PSS_KRYLOV_MAXDIM, default 40, no restart); if the limit is reached with
 * a relative residual still above 0.5, that is stagnation and the whole
 * analysis falls back to the dense path (fast-state detection + dense
 * Newton, see PSSshootDriven), as does a product whose period runs fail.
 *
 * Product M*v: v (unit 2-norm in the weighted variables) maps to the
 * physical direction d_k = scale_k v_k; the perturbation is p = sigma d with
 * sigma such that max_k |p_k| / h_k = 1, h_k being the engine's own per-state
 * finite-difference step (|x|*eps, floor 1e-6, 1 % of the ripple for a state
 * smaller than its ripple -- the exact logic of the dense FD columns), i.e.
 * no state is ever moved by more than its dense FD step. Centred difference
 * (P(x0+p) - P(x0-p))/2 (the dense engine found one-sided differences give
 * no usable direction on R3); PSS_KRYLOV_ONESIDED=1 uses P(x0+p) - P(x0)
 * with the baseline's own P(x0). A failing product retries at sigma/4
 * (every deck: the product only sets the direction, never the verdict).
 *
 * Stability: no M exists, so the Floquet report comes from Arnoldi: after
 * convergence a dedicated Arnoldi process (products M*v from a fixed
 * pseudo-random start vector, grown until the spectral radius settles, at
 * most PSS_KRYLOV_ARNOLDI = 40 products, see pss_kry_arnoldi()) gives the
 * Hessenberg matrix whose eigenvalues (Ritz values) approximate the leading
 * Floquet multipliers. It is reported as an "Arnoldi estimate", with the
 * change of the spectral radius between the last and the last-but-4
 * dimension as the convergence indicator; if that change exceeds 1 %, or
 * the eigen-solve fails, stability is reported "not computed" (the radius is
 * never invented). When the Krylov dimension equals the number of states the
 * Hessenberg spectrum is the whole spectrum. */

static int g_pss_krylov = 0;        /* this analysis runs matrix-free */
static int g_pss_kry_failed = 0;    /* set when the last pss_newton_run failed in the Krylov solve itself */

typedef struct {
    CKTcircuit *ckt;
    PSSstateMap *map;
    double T;
    int n;
    const double *x0;     /* the point the products are taken at (n_full states) */
    const double *w;      /* weights (ripple-relative scale) */
    const double *hs;     /* per-state FD step ceiling */
    const double *P0;     /* P(x0) for the one-sided variant, else NULL */
    int onesided;
    double *icbuf, *outp, *outm;
    long nprod, nshrink;
} PSSkry;

static double
pss_kry_dot(const double *a, const double *b, int n)
{
    double s = 0.0;
    int i;
    for (i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* the dense engine's own per-state FD step (pss_newton_run FD columns) */
static double
pss_kry_hstep(double x0j, double ripple_j, double eps)
{
    double h = fabs(x0j) * eps;
    if (h < 1e-6) h = 1e-6;
    if (ripple_j > 0.0 && fabs(x0j) < ripple_j && 0.01 * ripple_j > h &&
        getenv("PSS_SHOOT_NO_BIGH") == NULL)
        h = 0.01 * ripple_j;
    return h;
}

/* Mv = M~ v in weighted variables. Returns 0, or the failing ngspice error. */
static int
pss_kry_matvec(PSSkry *c, const double *v, double *Mv)
{
    int n = c->n, k, tries, retry_max, err = 0;
    double nv = sqrt(pss_kry_dot(v, v, n)), mx = 0.0, sigma;
    static int retry_base = -1;

    if (retry_base < 0) {
        retry_base = getenv("PSS_FD_RETRY_MAX") ? atoi(getenv("PSS_FD_RETRY_MAX")) : 5;
        if (retry_base < 1) retry_base = 1;
    }
    /* a failing product retries at sigma/4 on every deck (unlike the dense FD
     * columns' switch-only rule): the product only sets the Newton DIRECTION,
     * the convergence test is always the true residual, so a smaller step
     * cannot turn a non-converging orbit into a converged one. */
    retry_max = retry_base;
    for (k = 0; k < n; k++) Mv[k] = 0.0;
    if (!(nv > 0.0)) return 0;
    for (k = 0; k < n; k++) {
        double r = fabs(c->w[k] * v[k]) / nv / c->hs[k];
        if (r > mx) mx = r;
    }
    if (!(mx > 0.0)) return 0;
    sigma = 1.0 / mx;

    for (tries = 0; tries < retry_max; tries++) {
        for (k = 0; k < n; k++) c->icbuf[k] = c->x0[k] + sigma * c->w[k] * v[k] / nv;
        pss_state_set_ic(c->map, c->icbuf);
        err = pss_run_period_ex_safe(c->ckt, c->T, "PSS shoot (Krylov product +)", 0, 1);
        if (!err) {
            pss_state_get(c->map, c->ckt, c->outp);
            if (!c->onesided) {
                for (k = 0; k < n; k++) c->icbuf[k] = c->x0[k] - sigma * c->w[k] * v[k] / nv;
                pss_state_set_ic(c->map, c->icbuf);
                err = pss_run_period_ex_safe(c->ckt, c->T, "PSS shoot (Krylov product -)", 0, 1);
                if (!err) pss_state_get(c->map, c->ckt, c->outm);
            }
        }
        if (!err) break;
        if (getenv("PSS_SHOOT_DEBUG"))
            fprintf(stderr, "DEBUG Krylov product sigma=%.6g failed (ngspice error %d); shrinking "
                    "(try %d/%d).\n", sigma, err, tries + 1, retry_max);
        sigma *= 0.25;
        c->nshrink++;
    }
    if (err) return err;
    for (k = 0; k < n; k++) {
        double q = c->onesided ? (c->outp[k] - c->P0[k]) : 0.5 * (c->outp[k] - c->outm[k]);
        Mv[k] = q * nv / (sigma * c->w[k]);
    }
    c->nprod++;
    return 0;
}

/* GMRES for (M~ - I) y = b. y gets the solution. *relres_out: final relative
 * residual estimate, *dim_out: Krylov dimension used. Returns 0, or nonzero if a
 * product failed. */
static int
pss_kry_gmres(PSSkry *c, const double *b, double *y, double eta, double abs_floor,
              int maxdim, double *relres_out, int *dim_out)
{
    int n = c->n, m = maxdim < n ? maxdim : n, i, j, k, jdim = 0, err = 0;
    double beta = sqrt(pss_kry_dot(b, b, n)), res;
    double *V = NULL, *H = NULL, *cs = NULL, *sn = NULL, *g = NULL, *Av = NULL, *cc = NULL;

    for (k = 0; k < n; k++) y[k] = 0.0;
    *relres_out = 0.0;
    *dim_out = 0;
    if (!(beta > 0.0) || m < 1) return 0;
    V  = TMALLOC(double, (size_t) (m + 1) * n);
    H  = TMALLOC(double, (size_t) (m + 1) * m);
    cs = TMALLOC(double, m);
    sn = TMALLOC(double, m);
    g  = TMALLOC(double, m + 1);
    Av = TMALLOC(double, n);
    cc = TMALLOC(double, m);
    if (!V || !H || !cs || !sn || !g || !Av || !cc) { err = E_NOMEM; goto done; }
    for (k = 0; k < (m + 1) * m; k++) H[k] = 0.0;
    for (k = 0; k < n; k++) V[k] = b[k] / beta;
    g[0] = beta;
    res = beta;
    for (j = 0; j < m; j++) {
        double hn;
        err = pss_kry_matvec(c, V + (size_t) j * n, Av);
        if (err) goto done;
        for (k = 0; k < n; k++) Av[k] -= V[(size_t) j * n + k];   /* (M~ - I) v_j */
        {
            int pass;
            for (pass = 0; pass < 2; pass++)
                for (i = 0; i <= j; i++) {
                    double h = pss_kry_dot(Av, V + (size_t) i * n, n);
                    H[i * m + j] += h;
                    for (k = 0; k < n; k++) Av[k] -= h * V[(size_t) i * n + k];
                }
        }
        hn = sqrt(pss_kry_dot(Av, Av, n));
        H[(j + 1) * m + j] = hn;
        for (i = 0; i < j; i++) {   /* earlier rotations */
            double t = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
            H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
            H[i * m + j] = t;
        }
        {
            double a = H[j * m + j], bb = H[(j + 1) * m + j], r = hypot(a, bb);
            if (r > 0.0) { cs[j] = a / r; sn[j] = bb / r; } else { cs[j] = 1.0; sn[j] = 0.0; }
            H[j * m + j] = r;
            H[(j + 1) * m + j] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] = cs[j] * g[j];
        }
        jdim = j + 1;
        res = fabs(g[j + 1]);
        if (getenv("PSS_KRYLOV_DEBUG"))
            fprintf(stderr, "DEBUG gmres j=%d relres=%.3e (hn=%.3e)\n", j, res / beta, hn);
        if (res <= eta * beta || res <= abs_floor || hn <= 1e-13 * beta)
            break;
        for (k = 0; k < n; k++) V[(size_t) (j + 1) * n + k] = Av[k] / hn;
    }
    for (i = jdim - 1; i >= 0; i--) {
        double s = g[i];
        for (k = i + 1; k < jdim; k++) s -= H[i * m + k] * cc[k];
        cc[i] = (H[i * m + i] != 0.0) ? s / H[i * m + i] : 0.0;
    }
    for (i = 0; i < jdim; i++)
        for (k = 0; k < n; k++) y[k] += cc[i] * V[(size_t) i * n + k];
    *relres_out = res / beta;
    *dim_out = jdim;
 done:
    if (V) tfree(V);
    if (H) tfree(H);
    if (cs) tfree(cs);
    if (sn) tfree(sn);
    if (g) tfree(g);
    if (Av) tfree(Av);
    if (cc) tfree(cc);
    return err;
}

/* Arnoldi estimate of the leading Floquet multipliers of M at c->x0, from a
 * fixed pseudo-random start vector, grown in steps: after dimension 8, 12, 16,
 * 20, 24, 32, 40, ... (up to mmax, and n) the spectral radius of the Hessenberg
 * block is compared with that of its leading block of dimension-4; the process
 * stops at the first checkpoint where they agree to 2e-3 (relative) on the same
 * side of 1, at dimension n (the whole spectrum), or at an invariant subspace
 * (breakdown, also exact). Fills ritz[0..*m_eff), *rad (radius at m_eff),
 * *rad_prev (radius of the leading (m_eff-4)-block, -1 if m_eff < 8), *exact.
 * Returns 0, or the failing error / 1000 for an eigen-solve failure. */
static int
pss_kry_arnoldi(PSSkry *c, int mmax, PSScomplex *ritz, int *m_eff, double *rad,
                double *rad_prev, int *exact)
{
    int n = c->n, i, j, k, err = 0, mm = 0, m;
    double *V = NULL, *H = NULL, *Av = NULL, *Hm = NULL;
    PSScomplex *e2 = NULL;
    double nrm = 0.0;

    *m_eff = 0; *rad = -1.0; *rad_prev = -1.0; *exact = 0;
    m = mmax > n ? n : mmax;
    if (m < 1) return 1;
    V  = TMALLOC(double, (size_t) (m + 1) * n);
    H  = TMALLOC(double, (size_t) (m + 1) * m);
    Hm = TMALLOC(double, (size_t) m * m);
    Av = TMALLOC(double, n);
    e2 = TMALLOC(PSScomplex, m);
    if (!V || !H || !Hm || !Av || !e2) { err = E_NOMEM; goto done; }
    for (k = 0; k < (m + 1) * m; k++) H[k] = 0.0;
    {
        unsigned long s = 2463534242UL;
        for (k = 0; k < n; k++) {
            s ^= s << 13; s &= 0xffffffffUL; s ^= s >> 17; s ^= s << 5; s &= 0xffffffffUL;
            V[k] = ((double) (s % 20001UL) / 10000.0) - 1.0;
            nrm += V[k] * V[k];
        }
        nrm = sqrt(nrm);
        for (k = 0; k < n; k++) V[k] /= nrm;
    }
    for (j = 0; j < m; j++) {
        double hn;
        int pass, dim, checkpoint, r1 = 0, r2 = 0, stop = 0;
        err = pss_kry_matvec(c, V + (size_t) j * n, Av);
        if (err) goto done;
        for (pass = 0; pass < 2; pass++)
            for (i = 0; i <= j; i++) {
                double h = pss_kry_dot(Av, V + (size_t) i * n, n);
                H[i * m + j] += h;
                for (k = 0; k < n; k++) Av[k] -= h * V[(size_t) i * n + k];
            }
        hn = sqrt(pss_kry_dot(Av, Av, n));
        H[(j + 1) * m + j] = hn;
        dim = j + 1;
        mm = dim;
        if (hn < 1e-12) { *exact = 1; stop = 1; }
        if (dim == n) *exact = 1;
        checkpoint = (dim >= 8 && (dim <= 24 ? dim % 4 == 0 : dim % 8 == 0));
        if (*exact || dim == m || checkpoint) {
            for (i = 0; i < dim; i++)
                for (k = 0; k < dim; k++) Hm[i * dim + k] = H[i * m + k];
            if (pss_eigenvalues(Hm, dim, ritz) != 0) { err = 1000; goto done; }
            r1 = 1;
            {
                double r = 0.0;
                for (i = 0; i < dim; i++) { double a = hypot(ritz[i].re, ritz[i].im); if (a > r) r = a; }
                *rad = r;
            }
            *rad_prev = -1.0;
            if (dim >= 8) {
                int m2 = dim - 4;
                double *Hs = TMALLOC(double, (size_t) m2 * m2);
                if (Hs) {
                    double r = 0.0;
                    for (i = 0; i < m2; i++)
                        for (k = 0; k < m2; k++) Hs[i * m2 + k] = H[i * m + k];
                    if (pss_eigenvalues(Hs, m2, e2) == 0) {
                        for (i = 0; i < m2; i++) { double a = hypot(e2[i].re, e2[i].im); if (a > r) r = a; }
                        *rad_prev = r;
                        r2 = 1;
                    }
                    tfree(Hs);
                }
            }
            *m_eff = dim;
            if (*exact) stop = 1;
            else if (r1 && r2 && fabs(*rad - *rad_prev) <= 2e-3 * (*rad > 1e-3 ? *rad : 1e-3) &&
                     (*rad - 1.0) * (*rad_prev - 1.0) > 0.0)
                stop = 1;
        }
        if (stop || dim == m) break;
        for (k = 0; k < n; k++) V[(size_t) (j + 1) * n + k] = Av[k] / hn;
    }
    (void) mm;
 done:
    if (V) tfree(V);
    if (H) tfree(H);
    if (Hm) tfree(Hm);
    if (Av) tfree(Av);
    if (e2) tfree(e2);
    return err;
}

/* the identity "slow set" (every state a Newton unknown) the Krylov mode uses */
static int *
pss_kry_identity(int n, int *n_slow_out)
{
    int *idx = TMALLOC(int, n), k;
    for (k = 0; k < n; k++) idx[k] = k;
    *n_slow_out = n;
    fprintf(stdout, "PSS-SHOOT: krylov=1: matrix-free Newton-GMRES on all %d states -- no fast-state "
            "detection sweep, no assembled monodromy (card #537).\n", n);
    fflush(stdout);
    return idx;
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
    double *plus_raw = NULL;   /* card #498: see the FD perturbation loop below */
    int *is_slow = NULL;
    int baseline_loose = 0;
    int have_baseline = 0;   /* do x0/x1/F already hold a valid, just-run point? */
    int recompute = 1;       /* does THIS iteration's Jacobian need a fresh FD sweep? */
    int since_recompute = 0;
    int result = 0;
    /* card #498: does M currently hold a complete, valid monodromy estimate
     * (fresh FD, seeded, or Broyden-carried) that a LATER failing FD column
     * could fall back on? Set true only once a full sweep -- fresh or with
     * some columns themselves falling back -- has completed without an
     * unrecoverable failure, or once an M_init seed is consumed; see
     * pss_fd_perturb_side()'s call site below for the fallback itself. */
    int m_has_prior = 0;
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

    /* card #537: matrix-free Newton-GMRES for the Newton direction (see the
     * block comment above pss_kry_matvec()). Only on the full, un-reduced
     * unknown set (the caller passes the identity slow_idx in this mode) and
     * never with the analytic Jacobian. No Broyden: there is no M to update;
     * the Krylov solve is redone every iteration from fresh products. */
    int kry = g_pss_krylov && !use_analytic && (n == n_full);
    PSSkry kctx;
    double *khs = NULL, *kb = NULL, *ky = NULL, *kic = NULL, *kop = NULL, *kom = NULL;
    memset(&kctx, 0, sizeof kctx);
    g_pss_kry_failed = 0;
    if (kry) {
        broyden_enabled = 0;
        khs = TMALLOC(double, n); kb = TMALLOC(double, n); ky = TMALLOC(double, n);
        kic = TMALLOC(double, n); kop = TMALLOC(double, n); kom = TMALLOC(double, n);
    }

    x0      = TMALLOC(double, n);
    x1      = TMALLOC(double, n);
    F       = TMALLOC(double, n);
    scale   = TMALLOC(double, n);
    delta_x = TMALLOC(double, n);
    xtrial  = TMALLOC(double, n);
    x1trial = TMALLOC(double, n);
    A       = TMALLOC(double, kry ? 1 : n * n);
    b       = TMALLOC(double, n);
    M       = TMALLOC(double, kry ? 1 : n * n);
    s       = TMALLOC(double, n);
    y       = TMALLOC(double, n);
    Ms      = TMALLOC(double, n);
    icbuf       = TMALLOC(double, n_full);
    outbuf      = TMALLOC(double, n_full);
    ripple_full = TMALLOC(double, n_full);
    Phi_full    = TMALLOC(double, kry ? (size_t) 1 : (size_t) n_full * n_full);
    is_slow     = TMALLOC(int, n_full);
    plus_raw    = TMALLOC(double, n);
    if (!x0 || !x1 || !F || !scale || !delta_x || !xtrial || !x1trial || !A || !b ||
        !M || !s || !y || !Ms || !icbuf || !outbuf || !ripple_full || !Phi_full || !is_slow ||
        !plus_raw || (kry && (!khs || !kb || !ky || !kic || !kop || !kom))) {
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
                err = pss_run_period_ex_safe(ckt, T, "PSS shoot (Newton baseline period, analytic)", 1, 0);
                PSSsensActive = 0;
            } else {
                err = pss_run_period_ex_safe(ckt, T, "PSS shoot (Newton baseline period)", 1, 0);
            }
            PSSrippleRecordActive = 0;
            pss_ripple_get(map, ripple_full);
            baseline_loose = g_pss_last_used_loose;

            if (err) {
                fprintf(stderr, "PSS-SHOOT: %s period integration did not "
                        "complete at Newton iteration %d (ngspice error "
                        "%d); this is NOT a steady state.\n", stage_label, iter, err);
                result = 0;
                goto newton_done;
            }
            pss_state_get(map, ckt, outbuf);

            /* card #498: with a sw in the circuit the recording run (free
             * adaptive grid) and the forced runs (FD columns, trials) walk
             * different grids, and the trapezoid error at a conduction
             * kink makes their period maps differ by ~1e-7 V (R2 q=1:
             * orbit residual floor 1.4e-6 vs tol 1e-6, line search could
             * not reduce it). Newton solved the recording run's F but was
             * judged by the forced trial's F. Take F from a forced run at
             * the same x0 instead, so F, the columns and the trials all
             * belong to one map. The recording run only supplies the grid
             * (and the ripple). If the forced rerun fails, keep the
             * recording run's result. PSS_NO_FORCED_BASELINE=1 restores
             * the old behaviour. */
            if (!use_analytic && pss_ckt_has_switch(ckt) && !getenv("PSS_NO_FORCED_BASELINE")) {
                double *ob2 = TMALLOC(double, n_full);
                if (ob2) {
                    pss_state_set_ic(map, icbuf);
                    if (!pss_run_period_ex_safe(ckt, T, "PSS shoot (Newton baseline, forced rerun)", 0, 1)) {
                        pss_state_get(map, ckt, ob2);
                        memcpy(outbuf, ob2, (size_t) n_full * sizeof(double));
                    }
                    tfree(ob2);
                }
            }

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

        if (residual <= tol && baseline_loose) {
            fprintf(stderr, "PSS-SHOOT: %s iter=%d residual %.3e is below tol, but this baseline "
                    "period run only completed at the looser retry reltol; not declaring convergence.\n", stage_label, iter, residual);
        }
        if (residual <= tol && !baseline_loose) {
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
            m_has_prior = 1;   /* card #498: seeded M is a valid fallback source */
            if (seed_used_out) *seed_used_out = 1;
            if (getenv("PSS_SHOOT_DEBUG"))
                fprintf(stderr, "DEBUG %s iter=0: seeded Jacobian from "
                        "fast-state detection's own sweep (skipped a fresh "
                        "FD recompute).\n", stage_label);
        } else if (recompute && !use_analytic && !kry) {
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
                double h_plus, h_minus;
                int err_plus, err_minus;
                if (h < 1e-6) h = 1e-6;
                /* card #505: a state whose magnitude is below its own ripple
                 * (a comparator-reset latch node, the foot of a dead-time RC
                 * swing) has period-map noise (~1e-5..3e-5, exponential in
                 * the comparator on-time) far above the 1e-6 floor step, so
                 * its FD column is noise/h ~ 5-30 (garbage diagonal ->
                 * bogus Floquet multipliers). Difference it with 1% of its
                 * ripple instead. Every state with |x| >= ripple keeps the
                 * old step exactly. PSS_SHOOT_NO_BIGH disables. */
                if (ripple_full[jf] > 0.0 && fabs(x0[j]) < ripple_full[jf] &&
                    0.01 * ripple_full[jf] > h && getenv("PSS_SHOOT_NO_BIGH") == NULL)
                    h = 0.01 * ripple_full[jf];
                h_plus = h_minus = h;

                /* card #498: each side gets its OWN shrink-retry (see
                 * pss_fd_perturb_side()'s header comment) instead of the
                 * single fixed-h attempt this loop used before. plus_raw
                 * holds the raw f(x0+h) result (NOT yet a derivative) until
                 * -h is also known to have succeeded -- M's own column j is
                 * touched only once both sides are in hand, so a failing
                 * -h (after a successful +h) never leaves M's fallback
                 * column (the previous iteration's own valid estimate,
                 * m_has_prior below) partially overwritten. */
                for (k = 0; k < n; k++) plus_raw[k] = 0.0;
                err_plus = pss_fd_perturb_side(ckt, map, T, x_full, n_full, jf,
                                                x0[j], 1.0, h,
                                                "PSS shoot (monodromy perturbation +)",
                                                outbuf, &h_plus);
                if (!err_plus)
                    for (k = 0; k < n; k++) plus_raw[k] = outbuf[slow_idx[k]];

                if (!err_plus) {
                    err_minus = pss_fd_perturb_side(ckt, map, T, x_full, n_full, jf,
                                                     x0[j], -1.0, h,
                                                     "PSS shoot (monodromy perturbation -)",
                                                     outbuf, &h_minus);
                } else {
                    err_minus = 0;   /* -h not even attempted; err_plus alone triggers fallback below */
                }

                /* card #505: ONE-SIDED difference rung. A single failing side
                 * (a corrector collapse inside ONE perturbed period run --
                 * "Timestep too small", see docs/pss-pac.md "card #505",
                 * root cause (a)) no longer kills the whole Newton attempt
                 * for a switch-free deck: the surviving side is differenced
                 * against a FORCED run at x0 itself (same forced grid as the
                 * perturbed run, so the free-baseline-vs-forced mismatch that
                 * a naive (x1_free - forced) difference would divide by h
                 * cancels), giving an O(h) one-sided column. Only reached
                 * when a side FAILED, so every column that succeeds on both
                 * sides is computed exactly as before. */
                if ((err_plus || err_minus) && !g_pss_has_switch &&
                    getenv("PSS_NO_ONESIDED") == NULL) {
                    int have_plus = !err_plus, err_x0;
                    double *side = TMALLOC(double, n);
                    double h_x0 = 0.0;
                    if (side) {
                        if (have_plus) {
                            for (k = 0; k < n; k++) side[k] = plus_raw[k];
                        } else {
                            /* +h failed: try -h alone (it was not attempted above) */
                            err_minus = pss_fd_perturb_side(ckt, map, T, x_full, n_full, jf,
                                                             x0[j], -1.0, h,
                                                             "PSS shoot (monodromy perturbation -)",
                                                             outbuf, &h_minus);
                            if (!err_minus)
                                for (k = 0; k < n; k++) side[k] = outbuf[slow_idx[k]];
                        }
                        if (have_plus || !err_minus) {
                            err_x0 = pss_fd_perturb_side(ckt, map, T, x_full, n_full, jf,
                                                          x0[j], 0.0, h,
                                                          "PSS shoot (monodromy forced baseline at x0)",
                                                          outbuf, &h_x0);
                            if (!err_x0) {
                                double hh = have_plus ? h_plus : h_minus;
                                for (k = 0; k < n; k++)
                                    M[k * n + j] = (have_plus ? (side[k] - outbuf[slow_idx[k]])
                                                              : (outbuf[slow_idx[k]] - side[k])) / hh;
                                fprintf(stderr, "PSS-SHOOT: %s monodromy perturbation %s run for state %s "
                                        "did not complete; used a one-sided difference for this column "
                                        "at Newton iteration %d (card #505).\n", stage_label,
                                        have_plus ? "-h" : "+h", map->refs[jf].name, iter);
                                err_plus = err_minus = 0;
                                tfree(side);
                                continue;
                            }
                        }
                        tfree(side);
                    }
                    /* fall through: could not rescue this column */
                }

                if (err_plus || err_minus) {
                    /* card #498: the "reuse the previous column" rescue is
                     * gated the same way as the retry itself -- see
                     * pss_ckt_has_switch()'s header comment. A switch-free
                     * deck's column either succeeds on its one attempt
                     * (retry_max==1 above) or falls straight to the
                     * original "no previous column, abort" branch below,
                     * byte-for-byte the pre-card-#498 behaviour. */
                    if (m_has_prior && g_pss_has_switch) {
                        fprintf(stderr, "PSS-SHOOT: %s monodromy perturbation "
                                "column for state %s did not complete even "
                                "after shrinking h (ngspice error +h=%d -h=%d) "
                                "at Newton iteration %d -- reusing the "
                                "previous iteration's column for this state "
                                "only (card #498 guarded fallback).\n",
                                stage_label, map->refs[jf].name, err_plus,
                                err_minus, iter);
                        /* M[*, j] already holds the previous (Broyden-carried
                         * or seeded) estimate for this column -- leave it
                         * untouched rather than overwrite with a partial or
                         * failed computation. */
                    } else {
                        fprintf(stderr, "PSS-SHOOT: %s monodromy perturbation "
                                "run for state %s did not complete (ngspice "
                                "error +h=%d -h=%d) at Newton iteration %d%s.\n",
                                stage_label, map->refs[jf].name, err_plus,
                                err_minus, iter,
                                g_pss_has_switch ?
                                    (m_has_prior ? "" : ", and no previous "
                                     "column exists yet to fall back on") :
                                    " (no switch device in this circuit -- "
                                    "card #498's guarded fallback stays "
                                    "inactive)");
                        result = 0;
                        goto newton_done;
                    }
                } else {
                    for (k = 0; k < n; k++)
                        M[k * n + j] = (plus_raw[k] - outbuf[slow_idx[k]]) / (h_plus + h_minus);
                }
            }
            if (getenv("PSS_SHOOT_DEBUG"))
                for (k = 0; k < n * n; k++)
                    fprintf(stderr, "DEBUG FD Phi[%d]=%.8g\n", k, M[k]);
            since_recompute = 0;
            m_has_prior = 1;   /* card #498: this sweep (fresh, possibly with
                                 * some columns falling back) is now a valid
                                 * fallback source for a LATER iteration. */
            /* M is the RAW monodromy estimate here (no -I baked in) --
             * kept that way throughout (Broyden update below preserves
             * it) so Phi_out is always a true Phi, and the -I is
             * subtracted only into the disposable solve copy A just
             * below, same as the pre-third-pass code did. */
        }

        if (kry) {
            /* card #537: inexact Newton by GMRES, no assembled M */
            double relres = 0.0, eta = 1e-3, abs_floor = 0.01 * tol;
            int kdim = 0, maxdim = 40, gerr, bad = 0;
            long runs0 = g_pss_nruns;
            if (getenv("PSS_KRYLOV_ETA")) eta = atof(getenv("PSS_KRYLOV_ETA"));
            if (eta <= 0.0 || eta >= 1.0) eta = 1e-3;
            if (getenv("PSS_KRYLOV_MAXDIM")) maxdim = atoi(getenv("PSS_KRYLOV_MAXDIM"));
            if (maxdim < 2) maxdim = 2;
            for (k = 0; k < n; k++) {
                khs[k] = pss_kry_hstep(x0[k], ripple_full[k], eps);
                kb[k] = -F[k] / scale[k];
            }
            kctx.ckt = ckt; kctx.map = map; kctx.T = T; kctx.n = n;
            kctx.x0 = x0; kctx.w = scale; kctx.hs = khs; kctx.P0 = x1;
            kctx.onesided = getenv("PSS_KRYLOV_ONESIDED") && atoi(getenv("PSS_KRYLOV_ONESIDED")) != 0;
            kctx.icbuf = kic; kctx.outp = kop; kctx.outm = kom;
            gerr = pss_kry_gmres(&kctx, kb, ky, eta, abs_floor, maxdim, &relres, &kdim);
            if (gerr) {
                fprintf(stdout, "PSS-SHOOT: %s iter=%d krylov: a product period run did not complete "
                        "(ngspice error %d); falling back to the dense path.\n", stage_label, iter, gerr);
                fflush(stdout);
                g_pss_kry_failed = 1;
                result = 0;
                goto newton_done;
            }
            for (k = 0; k < n; k++) {
                delta_x[k] = scale[k] * ky[k];
                if (!(fabs(delta_x[k]) < 1e300)) bad = 1;
            }
            fprintf(stdout, "PSS-SHOOT: %s iter=%d krylov: GMRES dim %d relres %.2e, %ld period runs in the solve\n",
                    stage_label, iter, kdim, relres, g_pss_nruns - runs0);
            fflush(stdout);
            if (bad || relres > 0.5) {
                fprintf(stdout, "PSS-SHOOT: %s iter=%d krylov: GMRES stagnated (relres %.2e after dim %d); "
                        "falling back to the dense path.\n", stage_label, iter, relres, kdim);
                fflush(stdout);
                g_pss_kry_failed = 1;
                result = 0;
                goto newton_done;
            }
        } else {
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
        }   /* card #537: end of the dense (non-Krylov) solve */

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
            /* card #525: doomed-probe early exit. A probe attempt
             * (trust_radius>0) whose Newton step is clipped so hard that
             * even taking the clipped fraction alpha on EVERY remaining
             * iteration covers less than one full Newton step cannot
             * reach the orbit within its budget: the caller discards a
             * failed probe's wandering and restarts the next, longer
             * warm-up from the stage-start state, so finishing this
             * attempt is pure waste (measured: r3_longsettle stages 3-4,
             * alpha ~0.06/0.17 against 4 iterations; stages 1-2 alpha
             * ~1e-6). The final stage (trust_radius<=0) is never cut.
             * PSS_SHOOT_NO_DOOMED_EXIT=1 restores the old behaviour. */
            if (trust_radius > 0.0 && alpha < 1.0 &&
                !getenv("PSS_SHOOT_NO_DOOMED_EXIT") &&
                alpha * (double)(maxiter_this - iter) < 1.0) {
                fprintf(stdout, "PSS-SHOOT: %s iter=%d trust-clipped step "
                        "(alpha=%.3g) cannot reach the orbit in the %d "
                        "remaining iteration(s); abandoning this probe.\n",
                        stage_label, iter, alpha, maxiter_this - iter);
                fflush(stdout);
                result = 0;
                goto newton_done;
            }

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
                err = pss_run_period_ex_safe(ckt, T, "PSS shoot (damped-step trial)", 0, 1);
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
                    if (getenv("PSS_SHOOT_DEBUG_TRIAL"))
                        for (k = 0; k < n; k++)
                            fprintf(stderr, "DEBUGTRIAL bt=%d k=%d %s xtrial-x0=%.4e x1trial=%.10g x1base=%.10g Ftrial=%.4e Fbase=%.4e\n",
                                    bt, k, map->refs[slow_idx[k]].name, xtrial[k] - x0[k], x1trial[k], x1[k],
                                    x1trial[k] - xtrial[k], F[k]);
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
    if (Phi_out && M && !kry)
        for (k = 0; k < n * n; k++) Phi_out[k] = M[k];
    if (khs) tfree(khs);
    if (kb) tfree(kb);
    if (ky) tfree(ky);
    if (kic) tfree(kic);
    if (kop) tfree(kop);
    if (kom) tfree(kom);
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
    if (plus_raw) tfree(plus_raw);
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
    pss_ripple_reset(map, icbuf);   /* card #505: ripple per state, for the large-step recheck below */
    PSSrippleRecordActive = 1;
    err = pss_run_period_ex_safe(ckt, T, "PSS shoot (fast-state detection baseline)", 1, 0);
    PSSrippleRecordActive = 0;
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
        err = pss_run_period_ex_safe(ckt, T, "PSS shoot (fast-state detection perturbation +)", 0, 1);
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
        err = pss_run_period_ex_safe(ckt, T, "PSS shoot (fast-state detection perturbation -)", 0, 1);
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
    /* card #505: NOISE-ROBUST detection. A tiny-magnitude state (a latch node
     * that a comparator resets to ~1e-2, a dead-time RC node at the foot of
     * its swing) is differenced with h = max(eps*|x|, 1e-6) ~ 1e-6, but the
     * period map's own corrector/event noise on such a state is ~1e-5..3e-5
     * (its end value depends exponentially on the comparator's on-time), so
     * the "column" it produces is noise/h ~ 5-30 instead of the ~0 a
     * forgetful state really has. Such a state is kept as a Newton unknown,
     * and its garbage diagonal entry then poisons both the Newton Jacobian
     * and the Floquet verdict (spectral radius 2..14 "UNSTABLE" on an orbit
     * whose true multiplier is 0.994; measured on the R4d fixtures, card
     * #505). For every state whose small-step column is above the
     * threshold, re-measure it with a step of 5% of its own ripple
     * (>> the noise) and keep the SMALLER of the two: a genuinely sensitive
     * state is above the threshold at both step sizes, a noisy one only at
     * the small step. Never makes a column larger, so it can only ever
     * reclassify noisy columns as fast. PSS_SHOOT_NO_BIGH disables. */
    if (getenv("PSS_SHOOT_NO_BIGH") == NULL) {
        double *ripple_d = TMALLOC(double, n);
        double *phi_small_save = TMALLOC(double, (size_t) n * n);
        if (ripple_d && phi_small_save) {
            memcpy(phi_small_save, Phi, (size_t) n * n * sizeof(double));
            pss_ripple_get(map, ripple_d);
            for (j = 0; j < n; j++) {
                double h_small = fabs(x_full[j]) * eps, h2, cn2 = 0.0;
                int ok = 1;
                if (h_small < 1e-6) h_small = 1e-6;
                h2 = 0.05 * ripple_d[j];
                if (colnorm[j] < threshold || h2 < 5.0 * h_small)
                    continue;
                for (i = 0; i < n; i++) icbuf[i] = x_full[i];
                icbuf[j] = x_full[j] + h2;
                pss_state_set_ic(map, icbuf);
                if (pss_run_period_ex_safe(ckt, T, "PSS shoot (fast-state large-step recheck +)", 0, 1)) ok = 0;
                else {
                    pss_state_get(map, ckt, outbuf);
                    for (k = 0; k < n; k++) Phi[k * n + j] = outbuf[k];
                }
                if (ok) {
                    icbuf[j] = x_full[j] - h2;
                    pss_state_set_ic(map, icbuf);
                    if (pss_run_period_ex_safe(ckt, T, "PSS shoot (fast-state large-step recheck -)", 0, 1)) ok = 0;
                    else {
                        pss_state_get(map, ckt, outbuf);
                        for (k = 0; k < n; k++) {
                            double v = fabs((Phi[k * n + j] - outbuf[k]) / (2.0 * h2));
                            if (v > cn2) cn2 = v;
                        }
                    }
                }
                /* restore the small-step column into Phi (the seed Jacobian
                 * must stay exactly what it was): recomputed below from the
                 * saved copy */
                if (ok && cn2 < colnorm[j]) {
                    if (cn2 < threshold)
                        fprintf(stdout, "PSS-SHOOT:   %-12s small-step colnorm=%.3e is finite-difference noise "
                                "(large-step colnorm=%.3e at h=%.3g) -> fast\n", map->refs[j].name, colnorm[j], cn2, h2);
                    colnorm[j] = cn2;
                    is_fast[j] = (cn2 < threshold);
                }
                for (k = 0; k < n; k++) Phi[k * n + j] = phi_small_save[k * n + j];
            }
        }
        if (ripple_d) tfree(ripple_d);
        if (phi_small_save) tfree(phi_small_save);
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
    int kry_requested = 0, kry_fellback = 0;   /* card #537 */

    t_wall_start = SPfrontEnd->IFseconds();

    /* card #537: matrix-free Newton-GMRES requested? card key krylov=1, or
     * NGSPICE_PSS_KRYLOV=1 in the environment (=0 in the environment forces
     * it off even if the card asks). Not with the analytic Jacobian. */
    g_pss_nruns = 0;
    g_pss_warm_periods = 0.0;
    g_pss_kry_failed = 0;
    {
        int want = job->PSSkrylov != 0;
        if (getenv("NGSPICE_PSS_KRYLOV"))
            want = atoi(getenv("NGSPICE_PSS_KRYLOV")) != 0;
        if (want && getenv("PSS_SHOOT_JACOBIAN") && strcmp(getenv("PSS_SHOOT_JACOBIAN"), "analytic") == 0) {
            fprintf(stdout, "PSS-SHOOT: krylov=1 is not combined with PSS_SHOOT_JACOBIAN=analytic; "
                    "using the analytic Jacobian.\n");
            want = 0;
        }
        g_pss_krylov = want;
        kry_requested = want;
    }

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
    /* card #498: the deck's own reltol is the fallback pss_run_period_ex_safe()
     * retries a failed period run at -- see that function's header comment.
     * Set once here (this function's only entry point that captures
     * save_reltol), reset to 0.0/0 (disabled) at shoot_cleanup below.
     * BOTH card #498 fallbacks are gated on pss_ckt_has_switch() -- see its
     * own header comment for why (the pss_r4_d55_unstable.cir red-proof
     * regression an ungated version of this fix caused). */
    g_pss_has_switch = pss_ckt_has_switch(ckt);
    g_pss_loose_reltol = g_pss_has_switch ? save_reltol : 0.0;
    memset(&floq, 0, sizeof floq);

    x0      = TMALLOC(double, n);
    x1      = TMALLOC(double, n);
    F       = TMALLOC(double, n);
    Phi     = TMALLOC(double, g_pss_krylov ? 1 : n * n);   /* card #537: no dense matrices in Krylov mode */
    delta_x = TMALLOC(double, n);
    xtrial  = TMALLOC(double, n);
    x1trial = TMALLOC(double, n);
    scale   = TMALLOC(double, n);
    A       = TMALLOC(double, g_pss_krylov ? 1 : n * n);
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
        g_pss_lag_enabled = !use_analytic;

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
        /* card #505: tighter node-voltage tolerance for the Newton attempts only
         * (PSS_SHOOT_VNTOL; 0 = the deck's own). The period map's noise floor
         * (what two runs at x0 and x0+1e-12 disagree by) is ~ vntol: 7e-7 V on
         * R3d's c1f against a convergence tolerance of 1e-6*0.0113 = 1.1e-8. */
        double tight_vntol = getenv("PSS_SHOOT_VNTOL") ? atof(getenv("PSS_SHOOT_VNTOL")) : 0.0;
        if (tight_vntol <= 0.0 || tight_vntol > save_volttol) tight_vntol = save_volttol;

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
        double *phi_full_detect = TMALLOC(double, g_pss_krylov ? (size_t) 1 : (size_t) n * n);
        double *M_init = TMALLOC(double, g_pss_krylov ? (size_t) 1 : (size_t) n * n);   /* n*n upper bound; only n_slow*n_slow used */
        int phi_full_detect_valid = 0;

        if (ckt->CKTstabTime <= 0.0) {
            pss_state_get_ic(&map, x0);
            /* card #458 FIFTH pass, task 2: detect+reduce ONCE, right
             * before the only Newton attempt this path makes -- see
             * pss_detect_fast_states()'s header comment for why this is a
             * one-time cost. */
            n_slow = n;
            slow_idx = g_pss_krylov ? pss_kry_identity(n, &n_slow) :
                       pss_detect_fast_states(ckt, &map, T, x0, eps, &n_slow,
                                               phi_full_detect, &phi_full_detect_valid);
            if (phi_full_detect_valid && M_init)
                for (k = 0; k < n_slow; k++)
                    for (j = 0; j < n_slow; j++)
                        M_init[k * n_slow + j] =
                            phi_full_detect[slow_idx[k] * n + slow_idx[j]];
            ckt->CKTreltol = tight_reltol;
            ckt->CKTvoltTol = tight_vntol;   /* card #505 */
            converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                             maxiter, tol, eps, -1.0, broyden_enabled,
                                             broyden_period, "[no warm-up]",
                                             &iter, &residual, Phi, use_analytic,
                                             (phi_full_detect_valid && M_init) ? M_init : NULL,
                                             NULL);
            ckt->CKTreltol = save_reltol;
            ckt->CKTvoltTol = save_volttol;   /* card #505 */
            if (!converged_flag && g_pss_krylov) {
                /* card #537: dense fallback (see the stage loop's identical block) */
                fprintf(stdout, "PSS-SHOOT: [no warm-up] the matrix-free Newton attempt did not converge%s; "
                        "continuing on the dense path (fast-state detection + dense monodromy).\n",
                        g_pss_kry_failed ? " (Krylov solve failed)" : "");
                fflush(stdout);
                g_pss_krylov = 0; g_pss_kry_failed = 0; kry_fellback = 1;
                tfree(phi_full_detect); tfree(M_init); tfree(Phi);
                phi_full_detect = TMALLOC(double, (size_t) n * n);
                M_init = TMALLOC(double, (size_t) n * n);
                Phi = TMALLOC(double, (size_t) n * n);
                tfree(slow_idx);
                n_slow = n;
                slow_idx = pss_detect_fast_states(ckt, &map, T, x0, eps, &n_slow, phi_full_detect, &phi_full_detect_valid);
                if (phi_full_detect_valid && M_init) {
                    for (k = 0; k < n_slow; k++)
                        for (j = 0; j < n_slow; j++)
                            M_init[k * n_slow + j] = phi_full_detect[slow_idx[k] * n + slow_idx[j]];
                }
                ckt->CKTreltol = tight_reltol;
                ckt->CKTvoltTol = tight_vntol;
                converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                                 maxiter, tol, eps, -1.0, broyden_enabled,
                                                 broyden_period, "[no warm-up/dense]",
                                                 &iter, &residual, Phi, use_analytic,
                                                 (phi_full_detect_valid && M_init) ? M_init : NULL, NULL);
                ckt->CKTreltol = save_reltol;
                ckt->CKTvoltTol = save_volttol;
            }
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
                    g_pss_warm_periods += stage_len / T;
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
                    slow_idx = g_pss_krylov ? pss_kry_identity(n, &n_slow) :
                               pss_detect_fast_states(ckt, &map, T, x0, eps, &n_slow,
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
                ckt->CKTvoltTol = tight_vntol;   /* card #505 */
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
                ckt->CKTvoltTol = save_volttol;   /* card #505 */
                if (getenv("PSS_SHOOT_PROFILE"))
                    fprintf(stderr, "PROFILE: stage %d newton attempt (%d iters) "
                            "took %.4fs\n", stage, iter, SPfrontEnd->IFseconds() - _pnt0);
                }
                /* card #537: the matrix-free attempt failed -- either its Krylov
                 * solve did (a product run did not complete, GMRES stagnated),
                 * or this is the final stage and Newton did not converge --
                 * continue this stage and every later one on the dense path
                 * (fast-state detection from the current x0, dense FD
                 * monodromy), exactly as an analysis without krylov=1 would.
                 * A probe stage that merely ran out of its short iteration
                 * budget stays on the Krylov path (the normal escalation). */
                if (!converged_flag && g_pss_krylov && (g_pss_kry_failed || last_stage)) {
                    int seed_dummy3 = 0;
                    fprintf(stdout, "PSS-SHOOT: [stage %d] the matrix-free Newton attempt did not converge%s; "
                            "continuing on the dense path (fast-state detection + dense monodromy).\n",
                            stage, g_pss_kry_failed ? " (Krylov solve failed)" : "");
                    fflush(stdout);
                    g_pss_krylov = 0; g_pss_kry_failed = 0; kry_fellback = 1;
                    tfree(phi_full_detect); tfree(M_init); tfree(Phi);
                    phi_full_detect = TMALLOC(double, (size_t) n * n);
                    M_init = TMALLOC(double, (size_t) n * n);
                    Phi = TMALLOC(double, (size_t) n * n);
                    tfree(slow_idx);
                    n_slow = n;
                    slow_idx = pss_detect_fast_states(ckt, &map, T, x0, eps, &n_slow,
                                                       phi_full_detect, &phi_full_detect_valid);
                    if (phi_full_detect_valid && M_init)
                        for (k = 0; k < n_slow; k++)
                            for (j = 0; j < n_slow; j++)
                                M_init[k * n_slow + j] =
                                    phi_full_detect[slow_idx[k] * n + slow_idx[j]];
                    snprintf(label, sizeof label, "[stage %d/%s+dense]", stage,
                             last_stage ? "final" : "probe");
                    ckt->CKTreltol = tight_reltol;
                    ckt->CKTvoltTol = tight_vntol;
                    converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                                     this_maxiter, tol, eps,
                                                     this_trust, broyden_enabled,
                                                     broyden_period, label,
                                                     &iter, &residual, Phi, use_analytic,
                                                     (phi_full_detect_valid && M_init) ? M_init : NULL,
                                                     &seed_dummy3);
                    ckt->CKTreltol = save_reltol;
                    ckt->CKTvoltTol = save_volttol;
                }
                /* card #505 (composition with #525): the final stage's very
                 * first baseline period run can fail outright (residual 0,
                 * no Newton iteration done: "Timestep too small" at the
                 * period's first edge) when the warm-up happened to end on
                 * an unlucky phase state -- R4d after #525's early probe
                 * exits changed the warm-up trajectory. Nothing is wrong
                 * with the orbit, only with that x0, so before any rung
                 * gives up, move x0 on by a few more free periods at the
                 * deck's own tolerance and retry (up to 3 times). Only
                 * reached on that no-iteration failure, so every deck that
                 * ran its baseline is untouched. PSS_SHOOT_NO_X0_NUDGE=1
                 * disables. */
                if (!converged_flag && last_stage && iter == 0 && residual == 0.0 &&
                    !getenv("PSS_SHOOT_NO_X0_NUDGE")) {
                    int nudge;
                    for (nudge = 1; nudge <= 3 && !converged_flag && iter == 0 && residual == 0.0; nudge++) {
                        int seed_dummy2 = 0;
                        double nlen = (4.0 * nudge + 1.0) * T;   /* 5, 9, 13 periods */
                        fprintf(stdout, "PSS-SHOOT: [stage %d/final] the baseline period run did not complete; "
                                "advancing x0 by %.0f free periods and retrying (nudge %d/3, card #505).\n",
                                stage, nlen / T, nudge);
                        fflush(stdout);
                        pss_state_set_ic(&map, x0_stage_start);
                        if (pss_run_period(ckt, nlen, "PSS shoot (x0 nudge)"))
                            break;
                        pss_state_get(&map, ckt, x0);
                        for (k = 0; k < n; k++) x0_stage_start[k] = x0[k];
                        pss_state_set_ic(&map, x0);
                        ckt->CKTreltol = tight_reltol;
                        ckt->CKTvoltTol = tight_vntol;
                        converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                                         this_maxiter, tol, eps,
                                                         this_trust, broyden_enabled,
                                                         broyden_period, label,
                                                         &iter, &residual, Phi, use_analytic,
                                                         NULL, &seed_dummy2);
                        ckt->CKTreltol = save_reltol;
                        ckt->CKTvoltTol = save_volttol;
                    }
                }
                /* card #505: replay rung. The classic forced-grid map (every
                 * recorded time forced as a breakpoint => backward-Euler
                 * restart + LTE re-derivation per gap) is not the same
                 * function as the free baseline run, so on a sharp-edged
                 * deck the trial residual bottoms out at that mismatch
                 * (R3d: 3e-4 plateau). Retry the final stage once with the
                 * exact step/order replay, which makes P(x0) reproduce the
                 * baseline. Only reached when the classic attempt FAILED,
                 * so every deck that converged before is untouched. */
                if (!converged_flag && last_stage && !pss_replay_enabled() &&
                    getenv("PSS_REPLAY") == NULL && g_rp_n > 0 && !g_pss_has_switch) {
                    int seed_dummy = 0;
                    g_replay_auto_used = 1;
                    fprintf(stdout, "PSS-SHOOT: [stage %d/final] the forced-grid Newton attempt did not "
                            "converge (residual %.3e); retrying once with exact step-sequence replay "
                            "(card #505).\n", stage, residual);
                    fflush(stdout);
                    /* continue from the classic attempt's own last ACCEPTED
                     * iterate (x0 is updated in place on every accepted step,
                     * and Newton only accepts residual-reducing steps), not
                     * from the warm-up endpoint: the classic run usually got
                     * within 1e-5..1e-4 of the orbit before it plateaued. */
                    g_replay_mode = 1;
                    snprintf(label, sizeof label, "[stage %d/final+replay]", stage);
                    ckt->CKTreltol = tight_reltol;
                    ckt->CKTvoltTol = tight_vntol;   /* card #505 */
                    converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                                     this_maxiter, tol, eps,
                                                     this_trust, broyden_enabled,
                                                     broyden_period, label,
                                                     &iter, &residual, Phi, use_analytic,
                                                     NULL, &seed_dummy);
                    ckt->CKTreltol = save_reltol;
                    ckt->CKTvoltTol = save_volttol;   /* card #505 */
                    if (!converged_flag) {
                        /* third rung: replay again from the warm-up endpoint
                         * (a different basin/grid than the classic attempt's
                         * last iterate; measured: R3d's shoot deck converges
                         * from here but plateaus at 2e-5 from the classic
                         * iterate, while its .pac/.pstb decks show the reverse). */
                        fprintf(stdout, "PSS-SHOOT: [stage %d/final+replay] did not converge either "
                                "(residual %.3e); retrying replay from the warm-up endpoint "
                                "(card #505).\n", stage, residual);
                        fflush(stdout);
                        for (k = 0; k < n; k++) x0[k] = x0_stage_start[k];
                        snprintf(label, sizeof label, "[stage %d/final+replay2]", stage);
                        ckt->CKTreltol = tight_reltol;
                        ckt->CKTvoltTol = tight_vntol;   /* card #505 */
                        converged_flag = pss_newton_run(ckt, &map, T, n_slow, slow_idx, x0,
                                                         this_maxiter, tol, eps,
                                                         this_trust, broyden_enabled,
                                                         broyden_period, label,
                                                         &iter, &residual, Phi, use_analytic,
                                                         NULL, &seed_dummy);
                        ckt->CKTreltol = save_reltol;
                        ckt->CKTvoltTol = save_volttol;   /* card #505 */
                    }
                    if (!converged_flag)
                        g_replay_mode = -1;   /* back to the env/default choice */
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
        ckt->CKTvoltTol = tight_vntol;   /* card #505 */

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
    if (g_pss_krylov) {
        /* card #537: no monodromy exists -- Arnoldi estimate of the leading
         * Floquet multipliers from products at the converged orbit point
         * (see the block comment above pss_kry_matvec()). */
        int m_req = 40, m_eff = 0, exact = 0, aerr, k2;
        double rad = -1.0, radp = -1.0;
        long r0 = g_pss_nruns;
        double *rip = TMALLOC(double, n), *wk = TMALLOC(double, n), *hk = TMALLOC(double, n);
        double *ic = TMALLOC(double, n), *op = TMALLOC(double, n), *om = TMALLOC(double, n);
        double *P0 = TMALLOC(double, n);
        double rfloor = getenv("PSS_SHOOT_RIPPLE_FLOOR") ? atof(getenv("PSS_SHOOT_RIPPLE_FLOOR")) : 1e-2;
        PSSkry kc;
        memset(&kc, 0, sizeof kc);
        if (rfloor <= 0.0) rfloor = 1e-2;
        if (getenv("PSS_KRYLOV_ARNOLDI")) m_req = atoi(getenv("PSS_KRYLOV_ARNOLDI"));
        if (m_req < 2) m_req = 2;
        floq.valid = 0;
        if (rip && wk && hk && ic && op && om && P0) {
            pss_ripple_get(&map, rip);
            for (k2 = 0; k2 < n; k2++) {
                wk[k2] = (rip[k2] > rfloor) ? rip[k2] : rfloor;
                hk[k2] = pss_kry_hstep(x0[k2], rip[k2], eps);
            }
            kc.ckt = ckt; kc.map = &map; kc.T = T; kc.n = n; kc.x0 = x0; kc.w = wk; kc.hs = hk;
            kc.icbuf = ic; kc.outp = op; kc.outm = om; kc.P0 = P0;
            kc.onesided = getenv("PSS_KRYLOV_ONESIDED") && atoi(getenv("PSS_KRYLOV_ONESIDED")) != 0;
            aerr = 0;
            if (kc.onesided) {
                pss_state_set_ic(&map, x0);
                aerr = pss_run_period_ex_safe(ckt, T, "PSS shoot (Floquet baseline)", 0, 1);
                if (!aerr) pss_state_get(&map, ckt, P0);
            }
            if (!aerr)
                aerr = pss_kry_arnoldi(&kc, m_req, eig, &m_eff, &rad, &radp, &exact);
            if (aerr) {
                fprintf(stdout, "PSS-SHOOT: stability NOT computed (Arnoldi estimate failed, error %d).\n", aerr);
            } else {
                int ok = 1;
                pss_classify_floquet(eig, m_eff, &floq);
                if (!exact) {
                    double diff = (radp >= 0.0) ? fabs(rad - radp) : 1e30;
                    if (!(diff <= 1e-2 * (rad > 1e-3 ? rad : 1e-3)) || (rad - 1.0) * (radp - 1.0) <= 0.0)
                        ok = 0;   /* not converged by the maximum dimension, or ambiguous around 1 */
                }
                char cmp[64];
                if (radp >= 0.0) snprintf(cmp, sizeof cmp, " (%.6f at dimension %d)", radp, m_eff - 4);
                else cmp[0] = '\0';
                fprintf(stdout, "PSS-SHOOT: Arnoldi estimate of the Floquet multipliers: dimension %d of %d%s, "
                        "%ld period runs; spectral radius %.6f%s%s\n",
                        m_eff, n, exact ? " (whole spectrum)" : "", g_pss_nruns - r0, rad, cmp,
                        ok ? "" : " -- NOT converged / ambiguous around 1: stability not computed");
                if (!ok) floq.valid = 0;
            }
        } else {
            fprintf(stdout, "PSS-SHOOT: stability NOT computed (out of memory in the Arnoldi estimate).\n");
        }
        if (rip) tfree(rip);
        if (wk) tfree(wk);
        if (hk) tfree(hk);
        if (ic) tfree(ic);
        if (op) tfree(op);
        if (om) tfree(om);
        if (P0) tfree(P0);
    } else {
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
            if (pss_run_period_ex_safe(ckt, T, "PSS shoot (Floquet baseline, analytic)", 0, 1))
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
        if (floquet_ok && pss_run_period_ex_safe(ckt, T, "PSS shoot (Floquet baseline)", 0, 1))
            floquet_ok = 0;

        for (j = 0; floquet_ok && j < n_slow; j++) {
            int jf = slow_idx[j];
            double h = fabs(x0[jf]) * eps_floquet;
            if (h < 1e-6) h = 1e-6;

            for (k = 0; k < n_full_fl; k++) icbuf_fl[k] = x0[k];
            icbuf_fl[jf] = x0[jf] + h;
            pss_state_set_ic(&map, icbuf_fl);
            if (pss_run_period_ex_safe(ckt, T, "PSS shoot (Floquet perturbation +)", 0, 1)) {
                floquet_ok = 0;
                break;
            }
            pss_state_get(&map, ckt, outbuf_fl);
            for (k = 0; k < n_slow; k++)
                Phi2[k * n_slow + j] = outbuf_fl[slow_idx[k]];

            for (k = 0; k < n_full_fl; k++) icbuf_fl[k] = x0[k];
            icbuf_fl[jf] = x0[jf] - h;
            pss_state_set_ic(&map, icbuf_fl);
            if (pss_run_period_ex_safe(ckt, T, "PSS shoot (Floquet perturbation -)", 0, 1)) {
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

    /* card #505: the automatic replay rung only rescues a run onto a STABLE
     * orbit. Its smooth map can also find the (mathematically real but
     * unobservable) period-1 fixed point of a circuit whose true steady state
     * is period-2 -- refs/README.md's D~0.55 no-slope-compensation R4, the
     * test_pss red-proof, converges there with spectral radius 1.18 -- and
     * the FD Floquet of a noisy deck can also read spuriously > 1. Neither
     * is accepted as a steady state (the classic path never reported them
     * either: it plateaued). */
    if (g_replay_auto_used && floq.valid && floq.unstable) {
        fprintf(stderr, "PSS-SHOOT: the replay rung reached a fixed point with spectral radius %.4f >= 1 "
                "(residual %.3e); not accepted as a steady state (card #505).\n",
                floq.spectral_radius, residual);
        fprintf(stderr, "PSS-SHOOT: FAILED to converge within %d Newton iterations (final residual %.6e, "
                "tol %.3e); this is NOT a steady state.\n", maxiter, residual, tol);
        goto shoot_fail;
    }

    {
        /* card #537: in Krylov mode the verdict token carries the "arnoldi"
         * marker (the spectral radius is an estimate from the Hessenberg
         * matrix, not from an assembled monodromy) and an unknown stability
         * says "not-computed" -- never a made-up radius. */
        char vbuf[64];
        const char *verdict;
        if (!floq.valid)
            verdict = g_pss_krylov ? "UNKNOWN(not-computed)" : "UNKNOWN(eigensolve-failed)";
        else {
            const char *base = floq.unstable ? (floq.period_doubling ? "UNSTABLE(period-doubling)" :
                                                 floq.oscillation ? "UNSTABLE(oscillation)" : "UNSTABLE")
                                             : "stable";
            if (g_pss_krylov) {
                snprintf(vbuf, sizeof vbuf, "%s[arnoldi]", base);
                verdict = vbuf;
            } else
                verdict = base;
        }
        fprintf(stdout, "\nPSS-SHOOT: CONVERGED iter=%d residual=%.6e "
                "spectral_radius=%.6f %s wall_time=%.3fs\n",
                iter, residual, floq.valid ? floq.spectral_radius : -1.0, verdict,
                t_wall_end - t_wall_start);
    }
    if (floq.valid)
        fprintf(stdout, "PSS-SHOOT: FLOQUET dominant_multiplier=%.6g%+.6gi "
                "angle_deg=%.3f%s\n", floq.dominant.re, floq.dominant.im,
                atan2(floq.dominant.im, floq.dominant.re) * 180.0 / M_PI,
                g_pss_krylov ? " (Arnoldi estimate)" : "");
    if (kry_requested || getenv("PSS_SHOOT_COUNT_RUNS"))
        fprintf(stdout, "PSS-SHOOT: period runs=%ld (+ %.0f warm-up transient periods), mode=%s\n",
                g_pss_nruns, g_pss_warm_periods,
                !kry_requested ? "dense" : kry_fellback ? "krylov-fell-back-to-dense" : "krylov");

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
    /* card #505: the emission follows the Newton run's own map (replay if the
     * replay rung converged, classic forced otherwise). Emitting via replay
     * after a CLASSIC convergence was tried and made R4d's orbit 1-2 % of
     * ripple off the .tran reference (vs 0.01 % classic). */
    err = pss_run_period_ex_safe(ckt, T, "Time Domain Periodic Steady State Analysis", 0, 1);
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
    g_pss_loose_reltol = 0.0;   /* card #498: disable the fallback outside this call */
    g_pss_has_switch = 0;
    g_pss_lag_enabled = 1;
    g_pss_krylov = 0; g_pss_kry_failed = 0;   /* card #537 */
    PSSbreakRecordActive = PSSbreakInstallActive = 0;
    PSSreplayActive = 0;
    g_replay_mode = -1;   /* card #505: replay choice is per analysis */
    g_replay_auto_used = 0;
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


/* ==================================================================== */
/* card #489 (plan T6.2): FREE-RUNNING (autonomous) shooting.

Design (see docs/pss-pac.md's free-running section for the full writeup;
this comment is the short version kept next to the code):

  - `.pss Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff
    [uic] method=shoot trigger=<node> [trigval=<v>] [trigdir=rise|fall]`.
    OscNode is unused here (same as the driven engine's own comment),
    Fguess is a rough PERIOD ESTIMATE (not an exact clock -- unlike the
    driven engine, this circuit's period is an OUTPUT, not an input).

  - The Poincare section is g(x,t) := v(trigger) - trigval = 0, crossed in
    the given direction (trigdir). x0 is defined as "the shooting state
    (native capacitor voltages / inductor currents, same state vector the
    driven engine uses) at the instant of a trigger crossing" -- so x0
    itself always satisfies g(x0)=0 at the converged fixed point, though
    intermediate Newton iterates and FD perturbations do not (see below).

  - T is NOT a separate Newton unknown (the task's "or use a bordered
    Newton with a phase condition if cleaner; justify" alternative was
    considered and NOT used): instead, exactly as SIMPLIS's POP algorithm
    does (research/simplis-pop-ac.md Sec.2, Sec.6 Stage 1), T falls out of
    EVENT LOCATION -- pss_run_until_trigger() below integrates from a
    given x0 until the FIRST qualifying crossing, however long that takes
    (up to MAX_PERIOD), and returns both the endpoint state x(T) and the
    measured T itself. This is simpler to implement correctly on top of
    this file's existing machinery (SAME period-map/Newton/dense-solve/
    eigenvalue building blocks the driven engine already uses) than a
    bordered (n+1)x(n+1) system, and -- see the next paragraph -- it also
    sidesteps needing the explicit terminal-projection formula.

  - Monodromy: the Newton unknown is x0 (n-dimensional, unconstrained
    during the search -- see below), and the residual is the ordinary
    shooting residual F(x0) = x(T(x0)) - x0, exactly the driven engine's
    own F, just with T now itself a function of x0. The monodromy M used
    in Newton is the plain central finite difference of this SAME map:
    M[:,j] = (P(x0+h*e_j) - P(x0-h*e_j)) / (2h), where P(x0) :=
    x(T(x0)) is the FULL run_until_trigger() round trip (integrate, then
    locate and land exactly on the NEXT crossing) -- run independently for
    each perturbed x0, so each perturbed trajectory finds ITS OWN
    (slightly different) crossing time. This matters: research/simplis-
    pop-ac.md Sec.6 Stage 1 gives the analytic alternative, a terminal
    projection M = (I - f(T) grad(g)^T / (grad(g).f(T))) * Psi(T), needed
    because the RAW variational-equation monodromy Psi(T) (state-to-state
    sensitivity at a FIXED T) has a spurious unit ("trivial") eigenvalue
    in the flow direction and must be projected onto the section before
    Newton/Floquet make sense. This file's FD scheme above needs NO
    separate projection step: because each perturbed run relocates its
    OWN crossing (not a fixed T), its output x(T(x0+h*e_j)) is BY
    CONSTRUCTION always exactly ON the section (g=0) for every column,
    the same way the baseline's output is -- so the FD difference quotient
    directly estimates the derivative of the SECTION-TO-SECTION return
    map, which IS M, mechanically reproducing what the projection formula
    computes analytically. (The perturbed x0+h*e_j ITSELF is generally
    slightly OFF the section -- e.g. perturbing along a direction with a
    component normal to g=0 -- but pss_run_until_trigger() does not
    require its INPUT to be on the section, only integrates from
    wherever it is given until the first qualifying crossing, so this is
    not a problem.) Consequence: M as computed here is RANK n-1, not n --
    it has exactly one eigenvalue collapsed to (near) ZERO, not one, in
    the direction the projection removed (see "Floquet" below); this is
    the expected, correct signature of a working terminal-projected
    monodromy, not a bug.

  - Saltation: NONE needed (design ruling #447 A, smooth switch models
    only, same restriction as the driven engine) -- every device in a
    circuit this engine is asked to run is a continuous function of its
    inputs, so "discretise-then-linearise is consistent" (research/
    ngspice-internals-pss-pac.md Sec.4.5) applies here exactly as it does
    to the driven engine's own FD monodromy: the per-step Jacobians the
    transient already computes carry the switching instant's dependence
    on the state, provided the timestep resolves the transition (same
    Points-driven CKTmaxStep discipline as the driven engine).

  - Floquet: eigenvalues of M (psseigen.c, unchanged). Because M is the
    ALREADY-PROJECTED monodromy (see above), its trivial multiplier sits
    near 0 (the direction the projection annihilates), NOT near 1 (where
    it would sit in the UNprojected Psi(T) -- the textbook "flow
    direction always has multiplier 1" fact). This engine identifies the
    trivial direction as whichever eigenvalue has the SMALLEST modulus
    and reports it separately, labelled "trivial (projected flow
    direction)"; the other n-1 eigenvalues are the genuine Floquet
    multipliers, and stability is max|lambda| < 1 over THOSE n-1 only
    (excluding the trivial one, same convention the task asks for, just
    relocated from 1 to 0 by this engine's particular way of computing
    M).

  - Honest failure: pss_run_until_trigger() returning "not found" (no
    qualifying crossing within MAX_PERIOD) is reported as FAILED, never a
    steady state -- this is what makes a collapsed-hysteresis-band or
    stable-DC-point red-proof deck come back "no periodic orbit" instead
    of silently reporting the (non-oscillating) DC operating point as if
    it were one. A trigger= node that does not exist, or is not wired to
    anything live, degrades the same way (its voltage never reaches
    trigval, or the run itself fails to converge) -- refused the same
    honest way, with no special-case code needed for it.

  - What this does NOT do (see the card #489 report for the full list):
    no Broyden/staged-warm-up sophistication (a single StabTime warm-up
    transient, then plain damped-step FD Newton -- the driven engine's
    fuller machinery was not ported here given the time available this
    pass); no analytic (discretise-then-linearise) Jacobian (FD only,
    same as the driven engine's own default); no fast-state reduction;
    `.pac`/`.pstb` on the converged free-running orbit are explicitly
    refused below with a clear message (the existing pacshoot.c engine's
    LPTV linearisation assumes a FIXED period and has no period-variation
    term, which a free-running orbit's small-signal response needs -- see
    docs/pss-pac.md's free-running section for what that term is and why
    it was not implemented this pass).
**********/

/* card #489: trigger-crossing detector. Same call-site/no-op-unless-active
 * pattern as PSSbreakRecordStep/PSSsensRecordStep/PSSrippleRecordStep
 * above -- see pssshoot.h. Freezes on the FIRST qualifying crossing after
 * t=0 (subsequent accepted steps are ignored once found), so repeated
 * calls within one long probe run never overwrite an already-found
 * bracket. t<=0 is skipped for the same reason the three existing hooks
 * skip it (ckt->CKTrhsOld is not trustworthy at the exact UIC-restart
 * instant) -- the practical effect is just that the very first accepted
 * step establishes "prev" with no crossing test, and testing starts from
 * the second accepted step on; MAX_PERIOD is many periods long, so this
 * costs nothing worth avoiding. */

int PSStrigActive = 0;

static PSSstateMap *g_trig_map = NULL;
static int    g_trig_node_num = -1;
static double g_trig_val = 0.5;
static int    g_trig_dir = 1;
static int    g_trig_found = 0;
static int    g_trig_have_prev = 0;
static double g_trig_prev_t = 0.0, g_trig_prev_v = 0.0;
static double *g_trig_prev_x = NULL;   /* [n] */
static double g_trig_t_lo = 0.0, g_trig_v_lo = 0.0;
static double g_trig_t_hi = 0.0, g_trig_v_hi = 0.0;
static double *g_trig_x_lo = NULL;     /* [n], snapshot at t_lo */

static void
pss_trig_reset(PSSstateMap *map, int node_num, double val, int dir)
{
    g_trig_map = map;
    g_trig_node_num = node_num;
    g_trig_val = val;
    g_trig_dir = dir;
    g_trig_found = 0;
    g_trig_have_prev = 0;
    if (!g_trig_prev_x) g_trig_prev_x = TMALLOC(double, map->n);
    if (!g_trig_x_lo)   g_trig_x_lo   = TMALLOC(double, map->n);
}

static void
pss_trig_free(void)
{
    if (g_trig_prev_x) tfree(g_trig_prev_x);
    if (g_trig_x_lo)   tfree(g_trig_x_lo);
    g_trig_prev_x = g_trig_x_lo = NULL;
    g_trig_map = NULL;
    g_trig_node_num = -1;
    g_trig_found = g_trig_have_prev = 0;
}

void
PSStrigRecordStep(CKTcircuit *ckt)
{
    double v;
    int k;

    if (!g_trig_map || g_trig_node_num < 0 || g_trig_found)
        return;
    if (ckt->CKTtime <= 0.0)
        return;

    v = ckt->CKTrhsOld[g_trig_node_num];

    if (g_trig_have_prev) {
        int cross = 0;
        if (g_trig_dir > 0 && g_trig_prev_v < g_trig_val && v >= g_trig_val)
            cross = 1;
        else if (g_trig_dir < 0 && g_trig_prev_v > g_trig_val && v <= g_trig_val)
            cross = 1;
        if (cross) {
            g_trig_t_lo = g_trig_prev_t; g_trig_v_lo = g_trig_prev_v;
            g_trig_t_hi = ckt->CKTtime;  g_trig_v_hi = v;
            for (k = 0; k < g_trig_map->n; k++)
                g_trig_x_lo[k] = g_trig_prev_x[k];
            g_trig_found = 1;
            return;   /* freeze: first qualifying crossing only */
        }
    }

    g_trig_prev_t = ckt->CKTtime;
    g_trig_prev_v = v;
    for (k = 0; k < g_trig_map->n; k++)
        g_trig_prev_x[k] = pss_state_value(&g_trig_map->refs[k], ckt);
    g_trig_have_prev = 1;
}

/* The free-running period map: integrate from x0 (written into the device
 * IC fields, same restart mechanism as pss_run_period()) until the FIRST
 * qualifying trigger crossing, refine the exact crossing time by regula
 * falsi/bisection (restarting from the bracket's pre-crossing state each
 * trial -- cheap, each trial spans at most one of the outer probe run's
 * own accepted steps), then read back the full shooting-state vector at
 * that exact instant. Returns 0 and fills x_T/T_measured on success, 1 on
 * failure (no crossing within max_period, or an integration error) --
 * failure is the caller's cue to report FAILED/non-convergence, never a
 * false steady state (see this file's card #489 header comment).
 *
 * card #489 SECOND pass (plan T6.2, task 5): x_prev_out/t_prev_out are
 * optional (NULL-able, every existing call site below passes NULL/NULL and
 * is byte-for-byte unaffected) outputs of the LAST accepted-step state/time
 * BEFORE the trigger bracket (g_trig_prev_x/g_trig_prev_t, which
 * PSStrigRecordStep() already freezes at exactly that point and which this
 * function would otherwise discard). Together with x_T/T_measured, this
 * gives a cheap backward-difference estimate of f(T) = dx/dt at the
 * trigger crossing on the converged orbit -- f(T) e_trig^T / (e_trig .
 * f(T)) is the flow-direction/terminal-projection term PSSshootFreerun()
 * uses below to enable `.pac`/`.pstb` on a free-running orbit (see that
 * function's own comment and docs/pss-pac.md's free-running section for
 * the derivation). No new run is needed for this: g_trig_prev_x is a
 * sample the SAME probe run already took. */
static int
pss_run_until_trigger(CKTcircuit *ckt, PSSstateMap *map, const double *x0,
                       int node_num, double trigval, int dir,
                       double max_period, double *x_T, double *T_measured,
                       double *x_prev_out, double *t_prev_out)
{
    int err, iter;
    double dt_lo, dt_hi, v_lo, v_hi, dt_mid = 0.0, v_mid = 0.0;

    pss_state_set_ic(map, x0);
    pss_trig_reset(map, node_num, trigval, dir);
    PSStrigActive = 1;
    err = pss_run_period(ckt, max_period, "PSS freerun (trigger probe)");
    PSStrigActive = 0;

    if (err || !g_trig_found) {
        pss_trig_free();
        return 1;
    }

    dt_lo = 0.0;                          v_lo = g_trig_v_lo;
    dt_hi = g_trig_t_hi - g_trig_t_lo;    v_hi = g_trig_v_hi;
    dt_mid = dt_hi; v_mid = v_hi;   /* fallback if the loop body never runs */

    for (iter = 0; iter < 40; iter++) {
        double denom = v_hi - v_lo;
        if (fabs(denom) < 1e-300)
            break;
        /* regula falsi, with a bisection fallback if it ever steps outside
         * the current bracket (can happen on a strongly curved crossing) */
        dt_mid = dt_lo + (dt_hi - dt_lo) * (trigval - v_lo) / denom;
        if (!(dt_mid > dt_lo && dt_mid < dt_hi))
            dt_mid = 0.5 * (dt_lo + dt_hi);

        pss_state_set_ic(map, g_trig_x_lo);
        err = pss_run_period(ckt, dt_mid, "PSS freerun (trigger refine)");
        if (err) { pss_trig_free(); return 1; }
        v_mid = ckt->CKTrhsOld[node_num];

        if (fabs(v_mid - trigval) <= 1e-12 * (1.0 + fabs(trigval)))
            break;
        if ((dir > 0 && v_mid < trigval) || (dir < 0 && v_mid > trigval)) {
            dt_lo = dt_mid; v_lo = v_mid;
        } else {
            dt_hi = dt_mid; v_hi = v_mid;
        }
        if (dt_hi - dt_lo < 1e-15 * (1.0 + dt_hi))
            break;
    }

    if (getenv("PSS_TRIG_TRACE"))
        fprintf(stderr, "TRIG517 t_lo=%.17g t_hi=%.17g v_lo0=%.17g v_hi0=%.17g dt_mid=%.17g v_mid=%.17g T=%.17g iters=%d\n",
                g_trig_t_lo, g_trig_t_hi, g_trig_v_lo, g_trig_v_hi, dt_mid, v_mid, g_trig_t_lo + dt_mid, iter);
    pss_state_set_ic(map, g_trig_x_lo);
    err = pss_run_period(ckt, dt_mid, "PSS freerun (trigger endpoint)");
    if (!err)
        pss_state_get(map, ckt, x_T);
    *T_measured = g_trig_t_lo + dt_mid;
    if (!err && x_prev_out && t_prev_out && g_trig_have_prev) {
        int k;
        for (k = 0; k < map->n; k++)
            x_prev_out[k] = g_trig_prev_x[k];
        *t_prev_out = g_trig_prev_t;
    }
    pss_trig_free();
    return err;
}


int
PSSshootFreerun(CKTcircuit *ckt, PSSan *job)
{
    PSSstateMap map;
    int n = 0, j, k, iter, err = OK;
    double *x0 = NULL, *x1 = NULL, *F = NULL, *M = NULL, *delta = NULL;
    double *xtrial = NULL, *x1trial = NULL, *scale = NULL, *ripple = NULL;
    double *A = NULL, *b = NULL;
    PSScomplex *eig = NULL;
    double T = 0.0, Tguess, max_period, tol, eps;
    double save_reltol, save_abstol, save_volttol;
    int node_num, dir;
    double trigval;
    int converged = 0, maxiter;
    double residual = 0.0;
    double t_wall_start, t_wall_end;

    map.n = 0; map.refs = NULL;

    if (job->PSSkrylov || (getenv("NGSPICE_PSS_KRYLOV") && atoi(getenv("NGSPICE_PSS_KRYLOV")) != 0))
        fprintf(stderr, "PSS-FREERUN: krylov=1 is not supported for free-running (trigger=) shooting "
                "(the Poincare-section map carries the period as an extra unknown); this run uses "
                "the dense monodromy (card #537).\n");

    if (!job->PSStrigNode) {
        fprintf(stderr, "PSS-FREERUN: internal error -- dispatched with no "
                "trigger node (should have been refused before this "
                "point).\n");
        return E_UNSUPP;
    }
    node_num = job->PSStrigNode->number;
    trigval  = job->PSStrigVal;
    dir      = job->PSStrigDir ? job->PSStrigDir : 1;

    if (pss_build_state_map(ckt, &map) || map.n == 0) {
        fprintf(stderr, "PSS-FREERUN: no capacitor/inductor states found -- "
                "nothing to shoot.\n");
        if (map.refs) tfree(map.refs);
        return E_UNSUPP;
    }
    n = map.n;

    Tguess = (job->PSSguessedFreq > 0.0) ? 1.0 / job->PSSguessedFreq : 1e-6;
    /* MAX_PERIOD: 10x the guess -- SIMPLIS's own free-running guidance
     * (research/simplis-pop-ac.md Sec.2: "Self-oscillating systems: ...
     * MAX_PERIOD should be about 10x"), looser than the driven-mode
     * (clocked) recommendation there (>=3x) since there is no known clock
     * to bound the search here. */
    max_period = 10.0 * Tguess;

    ckt->CKTmaxStep = Tguess / (job->PSSpoints > 0 ? (double) job->PSSpoints : 200.0);
    ckt->CKTdelmin  = 1e-9 * ckt->CKTmaxStep;
    ckt->CKTmode    = job->PSSmode | MODEUIC;
    ckt->CKTinitTime = 0.0;

    save_reltol  = ckt->CKTreltol;
    save_abstol  = ckt->CKTabstol;
    save_volttol = ckt->CKTvoltTol;
    if (1e-5 < save_reltol)
        ckt->CKTreltol = 1e-5;   /* same tightening rationale as the driven engine */

    x0 = TMALLOC(double, n); x1 = TMALLOC(double, n);
    F  = TMALLOC(double, n); M  = TMALLOC(double, (size_t) n * n);
    delta = TMALLOC(double, n);
    xtrial = TMALLOC(double, n); x1trial = TMALLOC(double, n);
    scale = TMALLOC(double, n); ripple = TMALLOC(double, n);
    A = TMALLOC(double, (size_t) n * n); b = TMALLOC(double, n);
    if (!x0 || !x1 || !F || !M || !delta || !xtrial || !x1trial ||
        !scale || !ripple || !A || !b) {
        err = E_NOMEM;
        goto fr_cleanup;
    }

    /* Warm-up: StabTime seconds of plain transient from the deck's own
     * device IC=/default state, mirroring the driven engine's own
     * pre-stabilisation budget (same card field) -- needed here because a
     * hysteretic/free-running converter's cold start can be far from its
     * eventual limit cycle (refs/r5_hysteretic_buck.cir's own header: "q
     * stays pinned at 1 ... until Vout first reaches vref_hi", settling
     * to steady switching only after ~400us on that deck) and this
     * engine, unlike the driven one, has no staged-warm-up/Broyden
     * escalation to fall back on if Newton is started too far from the
     * orbit. */
    pss_state_get_ic(&map, x0);
    if (job->PSSstabTime > 0.0) {
        pss_state_set_ic(&map, x0);
        err = pss_run_period(ckt, job->PSSstabTime, "PSS freerun (warm-up)");
        if (err) {
            fprintf(stderr, "PSS-FREERUN: WARNING -- the StabTime=%.6g s "
                    "warm-up transient failed (ngspice error %d); starting "
                    "Newton from the deck's own IC instead (may not "
                    "converge).\n", job->PSSstabTime, err);
        } else {
            pss_state_get(&map, ckt, x0);
        }
        err = OK;
    }

    tol     = (job->steady_coeff > 0.0) ? job->steady_coeff : 1e-6;
    maxiter = (job->sc_iter > 0) ? job->sc_iter : 50;
    /* card #503 SECOND pass: PSS_FREERUN_FD_EPS overrides the relative FD
     * step used for the Newton monodromy below (default unchanged, 1e-6)
     * -- experiment only, gated behind PSS_FREERUN_PAC_EXPERIMENTAL's own
     * general "not for production" status via this file's existing env-var
     * convention. See docs/pss-pac.md's card #503 SECOND pass section:
     * this was added to test whether the trigger-crossing FD monodromy is
     * noise-dominated at the default eps (tightening steady_coeff alone
     * was already shown NOT to help, and in fact to destabilise the
     * Floquet report -- a LARGER eps trades truncation bias for less FD
     * noise, the opposite direction). */
    eps     = getenv("PSS_FREERUN_FD_EPS") ? atof(getenv("PSS_FREERUN_FD_EPS")) : 1e-6;

    t_wall_start = SPfrontEnd->IFseconds();

    for (iter = 0; iter < maxiter; iter++) {
        double Tmeas;

        pss_ripple_reset(&map, x0);
        PSSrippleRecordActive = 1;
        err = pss_run_until_trigger(ckt, &map, x0, node_num, trigval, dir,
                                     max_period, x1, &Tmeas, NULL, NULL);
        PSSrippleRecordActive = 0;
        if (err) {
            fprintf(stderr, "PSS-FREERUN: no trigger crossing found within "
                    "MAX_PERIOD=%.6g s (10x the fguess estimate) at Newton "
                    "iteration %d -- either trigger=/trigval=/trigdir= "
                    "names the wrong node/threshold/direction, or this "
                    "circuit genuinely has no periodic orbit here (e.g. a "
                    "collapsed hysteresis band or a stable DC point that "
                    "never re-crosses); reporting FAILED, not a false "
                    "steady state.\n", max_period, iter);
            goto fr_fail;
        }
        pss_ripple_get(&map, ripple);
        T = Tmeas;

        residual = 0.0;
        for (k = 0; k < n; k++) {
            F[k] = x1[k] - x0[k];
            scale[k] = (ripple[k] > 1e-9) ? ripple[k] : 1e-9;
            { double r = fabs(F[k]) / scale[k]; if (r > residual) residual = r; }
        }
        fprintf(stdout, "PSS-FREERUN: iter=%d residual=%.6e period=%.9g "
                "fsw=%.6g\n", iter, residual, T, 1.0 / T);
        fflush(stdout);

        if (residual <= tol) { converged = 1; break; }

        /* central-FD monodromy -- see this file's card #489 header
         * comment for why this is already the terminal-projected M. */
        for (j = 0; j < n; j++) {
            double h = fabs(x0[j]) * eps; double Tp, Tm;
            if (h < 1e-9) h = 1e-9;

            for (k = 0; k < n; k++) xtrial[k] = x0[k];
            xtrial[j] += h;
            if (pss_run_until_trigger(ckt, &map, xtrial, node_num, trigval,
                                       dir, max_period, x1trial, &Tp,
                                       NULL, NULL)) {
                fprintf(stderr, "PSS-FREERUN: monodromy perturbation (+h, "
                        "state %d) lost the trigger crossing at iter=%d; "
                        "treating this attempt as non-convergent.\n", j, iter);
                goto fr_fail;
            }
            for (k = 0; k < n; k++) M[k * n + j] = x1trial[k];

            for (k = 0; k < n; k++) xtrial[k] = x0[k];
            xtrial[j] -= h;
            if (pss_run_until_trigger(ckt, &map, xtrial, node_num, trigval,
                                       dir, max_period, x1trial, &Tm,
                                       NULL, NULL)) {
                fprintf(stderr, "PSS-FREERUN: monodromy perturbation (-h, "
                        "state %d) lost the trigger crossing at iter=%d; "
                        "treating this attempt as non-convergent.\n", j, iter);
                goto fr_fail;
            }
            for (k = 0; k < n; k++)
                M[k * n + j] = (M[k * n + j] - x1trial[k]) / (2.0 * h);
        }

        for (k = 0; k < n * n; k++) A[k] = M[k];
        for (k = 0; k < n; k++) A[k * n + k] -= 1.0;
        for (k = 0; k < n; k++) b[k] = -F[k];
        if (pss_solve_dense(A, b, n, delta)) {
            fprintf(stderr, "PSS-FREERUN: singular (M-I) at iter=%d -- "
                    "cannot take a Newton step; reporting FAILED.\n", iter);
            goto fr_fail;
        }

        {
            double damp = 1.0;
            int bt, accepted = 0;
            for (bt = 0; bt < 10 && !accepted; bt++) {
                double newres = 0.0, Tt; int ok2;
                for (k = 0; k < n; k++) xtrial[k] = x0[k] + damp * delta[k];
                ok2 = !pss_run_until_trigger(ckt, &map, xtrial, node_num,
                                              trigval, dir, max_period,
                                              x1trial, &Tt, NULL, NULL);
                if (ok2) {
                    for (k = 0; k < n; k++) {
                        double r = fabs(x1trial[k] - xtrial[k]) / scale[k];
                        if (r > newres) newres = r;
                    }
                }
                if ((ok2 && newres < residual * 1.5) || bt == 9) {
                    for (k = 0; k < n; k++) x0[k] = xtrial[k];
                    accepted = 1;
                }
                damp *= 0.5;
            }
        }
    }

    if (!converged) {
        fprintf(stderr, "PSS-FREERUN: FAILED to converge within %d Newton "
                "iterations (final residual %.6e, tol %.3e); this is NOT a "
                "periodic orbit.\n", maxiter, residual, tol);
        goto fr_fail;
    }
    t_wall_end = SPfrontEnd->IFseconds();

    /* Floquet: see this file's card #489 header comment -- M is already
     * the projected monodromy, so its trivial direction sits near 0, not
     * 1; identify it as the smallest-modulus eigenvalue and report the
     * other n-1 as the genuine multipliers. */
    eig = TMALLOC(PSScomplex, n);
    if (eig && pss_eigenvalues(M, n, eig) == 0) {
        int trivial_idx = 0, kk;
        double trivial_mag = 1e300, specrad = -1.0;
        for (kk = 0; kk < n; kk++) {
            double mag = sqrt(eig[kk].re * eig[kk].re + eig[kk].im * eig[kk].im);
            if (mag < trivial_mag) { trivial_mag = mag; trivial_idx = kk; }
        }
        for (kk = 0; kk < n; kk++) {
            double mag;
            if (kk == trivial_idx) continue;
            mag = sqrt(eig[kk].re * eig[kk].re + eig[kk].im * eig[kk].im);
            if (mag > specrad) specrad = mag;
        }
        /* card #503: NEAR-DEGENERATE TRIVIAL EIGENVALUE diagnostic. The
         * "smallest-modulus = trivial (projected flow direction)" pick
         * above assumes that eigenvalue is well SEPARATED from the rest --
         * true on every fixture measured before this card, but R5
         * (refs/r5_hysteretic_buck.cir, Csw~1ps RC and Cq~1ns RC, both
         * orders of magnitude faster than the ~2.4us switching period)
         * was found this card to report a genuine near-zero COMPLEX-
         * CONJUGATE PAIR (both |.|~9.8e-6, at whatever two indices the
         * eigensolver happens to order them) instead of a single real
         * near-zero trivial eigenvalue the theory predicts -- i.e. TWO
         * fast/stiff directions collapse close to zero, not one, and the
         * "smallest modulus" tie-break picks one member of that pair
         * arbitrarily. This does not by itself change pacshoot.c's
         * trigger-node-based projection (pac_find_trig_idx() resolves
         * trig_idx from the trigger= NODE NAME, never from this
         * eigenvalue pick), but it is corroborating evidence (card #503's
         * report, docs/pss-pac.md) that R5's own finite-difference
         * monodromy sits close to its own numerical noise floor -- print
         * a loud warning rather than let a mislabelled "trivial" direction
         * pass silently, same honesty convention as this file's other
         * WARNING lines. Threshold: the two smallest magnitudes within
         * 5x of each other (conservative -- R5's own two multipliers were
         * IDENTICAL to 6 significant figures, a much tighter case than
         * this threshold requires). */
        if (n >= 3) {
            double second_mag = 1e300;
            for (kk = 0; kk < n; kk++) {
                double mag;
                if (kk == trivial_idx) continue;
                mag = sqrt(eig[kk].re * eig[kk].re + eig[kk].im * eig[kk].im);
                if (mag < second_mag) second_mag = mag;
            }
            if (second_mag > 0.0 && second_mag < 5.0 * trivial_mag) {
                fprintf(stderr, "PSS-FREERUN: WARNING -- the \"trivial\" "
                        "(projected flow direction) eigenvalue |.|=%.6g is "
                        "NOT well separated from the next-smallest genuine "
                        "Floquet multiplier |.|=%.6g (ratio %.3g < 5x) -- "
                        "the smallest-modulus tie-break may have picked the "
                        "wrong direction, or this orbit has more than one "
                        "near-zero (very fast/heavily damped) mode; the "
                        "Floquet report below may mislabel which multiplier "
                        "is \"trivial\". Common on circuits with more than "
                        "one time constant far faster than the switching "
                        "period (e.g. R5's Csw~1ps and Cq~1ns) -- see "
                        "docs/pss-pac.md's card #503 section.\n",
                        trivial_mag, second_mag, second_mag / trivial_mag);
            }
        }
        fprintf(stdout, "\nPSS-FREERUN: CONVERGED iter=%d residual=%.6e "
                "period=%.9g fsw=%.6g spectral_radius=%.6f %s "
                "wall_time=%.3fs\n", iter, residual, T, 1.0 / T,
                specrad, (n <= 1) ? "UNKNOWN(n=1, no non-trivial multiplier)" :
                (specrad < 1.0 ? "stable" : "UNSTABLE"),
                t_wall_end - t_wall_start);
        fprintf(stdout, "PSS-FREERUN: FLOQUET trivial(projected flow "
                "direction, index %d)=%.6g%+.6gi (expect near 0)\n",
                trivial_idx, eig[trivial_idx].re, eig[trivial_idx].im);
        for (kk = 0; kk < n; kk++) {
            double mag;
            if (kk == trivial_idx) continue;
            mag = sqrt(eig[kk].re * eig[kk].re + eig[kk].im * eig[kk].im);
            fprintf(stdout, "PSS-FREERUN: FLOQUET multiplier[%d]=%.6g%+.6gi "
                    "|.|=%.6g\n", kk, eig[kk].re, eig[kk].im, mag);
        }
    } else {
        fprintf(stderr, "PSS-FREERUN: WARNING -- the eigenvalue solver did "
                "not converge on the monodromy matrix; Floquet report is "
                "NOT trustworthy this run.\n");
        fprintf(stdout, "\nPSS-FREERUN: CONVERGED iter=%d residual=%.6e "
                "period=%.9g fsw=%.6g spectral_radius=UNKNOWN "
                "wall_time=%.3fs\n", iter, residual, T, 1.0 / T,
                t_wall_end - t_wall_start);
    }

    if (getenv("PSS_FREERUN_DEBUG")) {
        int kk;
        for (kk = 0; kk < n; kk++)
            fprintf(stderr, "DEBUG x0[%d] (%s) = %.10g\n", kk,
                    map.refs[kk].isCap ? "cap" : "ind", x0[kk]);
    }

    /* card #503 THIRD pass: PSS_FREERUN_X0_OVERRIDE (n whitespace-separated
     * doubles, SAME order as pss_build_state_map()'s own canonical
     * native-capacitor/inductor scan order -- the order PSS_FREERUN_DEBUG
     * just printed above, "DEBUG x0[k]") / PSS_FREERUN_T_OVERRIDE (one
     * double, seconds) replace this Newton loop's own converged x0/T with
     * caller-supplied values for EVERYTHING downstream (the emitted
     * trajectory AND the .pac/.pstb block below) -- test-only, to decide
     * the "hand-off precision" hypothesis directly: feed in an
     * independently, much more precisely converged x0/T (e.g. from a
     * tight-reltol .tran started at this Newton loop's own x0, run past
     * cycle-to-cycle stability) and see whether .pac/.pstb's own
     * disagreement against tools/pac_perturbed.py --freerun closes. If it
     * does, the cause is this Newton loop's own x0/T precision (fixable
     * inside this function); if it does not, the cause is in how the
     * free-running small-signal problem itself is posed (the period-
     * variation/phase-condition boundary condition in pacshoot.c), not in
     * hand-off precision. Never silent -- loud stderr line, exactly like
     * this file's other test-only overrides (PSS_FREERUN_FT_BACKWARD/
     * PSS_FREERUN_FT_DT). Only takes effect at all when
     * PSS_FREERUN_PAC_EXPERIMENTAL=1 is ALSO set (checked again below,
     * same gate the whole .pac/.pstb block already requires) -- this does
     * not change default (gated) behaviour either way, and does not touch
     * the CONVERGED/FLOQUET report already printed above (which reflects
     * this Newton loop's own x0/T, unmodified, so the large-signal
     * orbit/stability report stays honest even when the small-signal
     * .pac/.pstb input below is overridden). */
    if (getenv("PSS_FREERUN_X0_OVERRIDE") && getenv("PSS_FREERUN_PAC_EXPERIMENTAL")) {
        char *buf = strdup(getenv("PSS_FREERUN_X0_OVERRIDE"));
        char *tok = strtok(buf, " \t,");
        int kk = 0;
        while (tok && kk < n) {
            x0[kk] = atof(tok);
            kk++;
            tok = strtok(NULL, " \t,");
        }
        free(buf);
        if (kk != n) {
            fprintf(stderr, "PSS-FREERUN: WARNING -- PSS_FREERUN_X0_OVERRIDE "
                    "gave %d value(s), need exactly n=%d -- override NOT "
                    "fully applied, remaining components keep the Newton "
                    "loop's own converged value.\n", kk, n);
        }
        fprintf(stderr, "PSS-FREERUN: PSS_FREERUN_X0_OVERRIDE=1 -- x0 "
                "REPLACED by caller-supplied values (card #503 THIRD pass "
                "hand-off-precision test; NOT for production use):");
        for (kk = 0; kk < n; kk++)
            fprintf(stderr, " %.10g", x0[kk]);
        fprintf(stderr, "\n");
    }
    if (getenv("PSS_FREERUN_T_OVERRIDE") && getenv("PSS_FREERUN_PAC_EXPERIMENTAL")) {
        double T_new = atof(getenv("PSS_FREERUN_T_OVERRIDE"));
        fprintf(stderr, "PSS-FREERUN: PSS_FREERUN_T_OVERRIDE=1 -- T REPLACED "
                "%.12g -> %.12g (card #503 THIRD pass hand-off-precision "
                "test; NOT for production use)\n", T, T_new);
        T = T_new;
    }

    /* Emit the converged trajectory as an ordinary time-domain plot -- now
     * that T is converged/known, this is a single plain period run (same
     * mechanism as the driven engine's own emission), not another
     * trigger-relocation round trip. */
    pss_state_set_ic(&map, x0);
    err = pss_run_period(ckt, T, "Time Domain Periodic Steady State Analysis (free-running)");
    if (err) {
        fprintf(stderr, "PSS-FREERUN: WARNING -- converged, but the final "
                "trajectory output run failed (ngspice error %d); the "
                "CONVERGED line above is still accurate for the orbit "
                "itself, just not the emitted plot.\n", err);
    }

    /* card #489 SECOND pass (plan T6.2, task 5): `.pac`/`.pstb` on a
     * free-running orbit. pacshoot.c's LPTV linearisation (pac_pstb_run())
     * assumes a FIXED period T; a free-running orbit's small-signal
     * response needs the period-VARIATION term (how T itself shifts under
     * a perturbation feeds back into the periodicity boundary condition).
     * A correction was DERIVED and IMPLEMENTED this pass (docs/pss-pac.md's
     * free-running section has the full derivation): the SAME terminal
     * projection this file's own monodromy already needs for the
     * LARGE-signal (Newton/Floquet) problem, `Proj = I - f(T) e_trig^T /
     * (e_trig . f(T))`, applied to pacshoot.c's complex periodicity
     * boundary condition (Phi, z) instead of to the real FD monodromy M.
     *
     * NOT trusted by default (still gated -- see the "else" branch below):
     * validated on R5 control-to-output against an independent, tightly
     * self-consistent (amplitude-halving <=0.01dB) transient-injection
     * Hann-windowed-DFT reference at 20kHz (0.07dB/0.34deg from the
     * reference) and 60kHz (0.30dB/0.59deg) -- both inside the project's
     * usual 0.5dB/3deg target -- but at 150kHz the SAME comparison shows a
     * 3.7dB/9.6deg disagreement, worse than leaving the correction OFF
     * (3.0dB/-7.2deg uncorrected) though neither is close to right. Root
     * cause NOT identified this pass: plausible candidates are (1) the
     * single-local-step backward-difference f(T) estimate below being too
     * imprecise at R5's OWN extremely stiff switch-node time constant
     * (Csw's RC ~= 1ps, two orders of magnitude faster than the ~1ns
     * comparator/latch constant every other timing estimate in this file
     * assumes is the fast scale) -- the correction's own row weights
     * fT[d]/fT[trig_idx] could be amplifying that imprecision at higher
     * frequency; or (2) pacshoot.c's own documented (docs/pss-pac.md
     * "Accuracy") Jacobian-reuse floor, never previously characterised on
     * a circuit this stiff (R3/R4 do not have anything like Csw's 1ps
     * constant), being exposed for the first time now that .pac finally
     * runs on R5 at all. NOT distinguished this pass -- see the card #489
     * report. Because of this, .pac/.pstb on a free-running orbit STAYS
     * REFUSED by default (same honest-failure principle as before this
     * pass); PSS_FREERUN_PAC_EXPERIMENTAL=1 opts into the implementation
     * above for further investigation, never for production use. */
    if (job->PSSdoPAC || job->PSSdoPSTB) {
        double *fT = TMALLOC(double, n);
        double *x_prev = TMALLOC(double, n);
        double t_prev = 0.0, T_meas2 = 0.0;
        /* card #503 THIRD pass: POLISH x0/T for .pac/.pstb ONLY (a local
         * x0_pac/T_pac -- the trajectory emission just above and the
         * CONVERGED/FLOQUET report already printed both keep using the
         * plain Newton x0/T, untouched, so nothing already validated
         * changes). Root cause this pass: an independent, much more
         * precisely converged reference orbit (a tight-reltol .tran seeded
         * at this Newton loop's own x0, run to cycle-to-cycle stability --
         * see docs/pss-pac.md's card #503 THIRD pass section) differs from
         * this Newton loop's own x0/T by ~0.05% (T) and up to ~0.02% (the
         * inductor-current state) -- small in absolute terms, but feeding
         * that independent, high-precision x0/T into pac_pstb_run_ex()
         * below (PSS_FREERUN_X0_OVERRIDE/PSS_FREERUN_T_OVERRIDE, this same
         * pass) cuts the .pac disagreement against tools/pac_perturbed.py
         * --freerun by roughly 5-8x at every frequency tested (20-200kHz)
         * -- the DOMINANT cause is this hand-off precision, not the
         * period-variation formula itself (unchanged, still pacshoot.c's
         * card #489 derivation). Tightening THIS loop's own Newton
         * tolerance was already shown (card #503 first/second pass) to be
         * counterproductive (noise-dominated, sometimes non-convergent) --
         * but the orbit is heavily damped (spectral_radius ~0.05, ALL
         * non-trivial Floquet multipliers small), so a few PLAIN period-
         * map (successive-substitution, NO finite-difference Jacobian at
         * all) iterations starting from the Newton-converged x0 converge
         * geometrically at roughly that same ~0.05x-per-iteration rate --
         * cheap (each iteration is one more pss_run_until_trigger() call,
         * the SAME machinery Newton's own probe/monodromy steps already
         * use) and needs no new numerics. PSS_FREERUN_POLISH_ITERS
         * (default 5, 0 disables) controls the iteration count;
         * PSS_FREERUN_X0_OVERRIDE (above) takes precedence when given
         * (explicit external reference for testing) -- polish is skipped
         * in that case, not stacked on top of it. */
        double *x0_pac = x0;
        double T_pac = T;
        double *x0_polished = NULL;
        double save_reltol_polish = ckt->CKTreltol;
        int reltol_tightened = 0;
        /* card #517: the whole .pac/.pstb solver profile (reltol, trtol,
         * max step) is tightened/matched for the polish AND the recording
         * passes, and restored after -- see the block comment below. */
        double save_trtol_polish = ckt->CKTtrtol;
        double save_maxstep_polish = ckt->CKTmaxStep;
        double save_step_polish = ckt->CKTstep;
        double save_delmin_polish = ckt->CKTdelmin;
        long rec_points = 0;
        if (!getenv("PSS_FREERUN_X0_OVERRIDE")) {
            int polish_iters = getenv("PSS_FREERUN_POLISH_ITERS") ?
                atoi(getenv("PSS_FREERUN_POLISH_ITERS")) : 5;
            if (polish_iters > 0) {
                int pk, polish_ok = 1;
                /* card #503 THIRD pass, refinement: the polish iterations
                 * below reuse pss_run_until_trigger(), which runs at
                 * WHATEVER ckt->CKTreltol this function's own top-of-
                 * function tightening left in place -- capped at 1e-5 even
                 * when the deck's own reltol was looser (see save_reltol
                 * above). On R5 that 1e-5 cap was found this pass to be
                 * the dominant source of x0/T imprecision (an independent
                 * tight-reltol reference orbit differs from the capped-at-
                 * 1e-5 Newton x0/T by ~0.05% in T) -- but simply running
                 * the WHOLE shooting (Newton AND the recording pass
                 * pac_engine_run() uses) at a tight reltol throughout was
                 * measured to help less, and sometimes hurt, at higher
                 * frequency (a confound: the recording pass's own
                 * per-step-reused-Jacobian floor, docs/pss-pac.md
                 * "Accuracy", behaves differently at a very different
                 * step-size regime). Tightening reltol ONLY for these
                 * polish period-map evaluations -- restored to whatever it
                 * was immediately after, before pac_pstb_run_ex() runs the
                 * (unchanged-environment) recording pass -- isolates just
                 * the x0/T precision improvement, matching what actually
                 * closed most of the gap in this pass's own experiments
                 * (see docs/pss-pac.md). PSS_FREERUN_POLISH_RELTOL
                 * (default 1e-9) sets the target; never loosens (same
                 * "only tighten" convention as the existing save_reltol
                 * cap just above). Left in effect (see reltol_tightened,
                 * restored just before pac_pstb_run_ex() below) through
                 * the fT/f(T) estimation that follows this block too --
                 * that re-locates the SAME trigger crossing from x0_pac
                 * via the SAME pss_run_until_trigger(), so it benefits
                 * from the SAME tightened precision the polish loop just
                 * established, rather than dropping back to the coarser
                 * cap for that one call. */
                double polish_reltol = getenv("PSS_FREERUN_POLISH_RELTOL") ?
                    atof(getenv("PSS_FREERUN_POLISH_RELTOL")) : 1e-9;
                if (polish_reltol > 0.0 && polish_reltol < ckt->CKTreltol) {
                    ckt->CKTreltol = polish_reltol;
                    reltol_tightened = 1;
                }
                /* card #517: SELF-CONSISTENT solver profile. Root cause of
                 * the machine-to-machine (libm last-bit) disagreement of the
                 * free-running .pac (0.435 vs 0.556 dB at 100 kHz), see
                 * docs/pss-pac.md's card #517 section: the period T this
                 * polish converges to is a property of the SOLVER that
                 * measured it, not of the circuit -- a plain .tran of the
                 * same orbit at reltol=1e-9 gives T=2.36700733e-6 s with
                 * the default trtol=7 / max-step T/200 but 2.36710784e-6 s
                 * (converged; +42 ppm = 100 ps, comparable to R5's switch
                 * edge) with trtol=1 and a fine step -- so x0/T did not
                 * close the recording pass's own orbit, the recording
                 * pass's Phi did not map f(T) onto itself, and every
                 * discrete step-size decision (which libm last bits flip)
                 * moved the answer by up to 0.12 dB. Fix: run the polish
                 * AND the two recording passes with ONE tight profile
                 * (reltol 1e-9, trtol 1, max step T/pacpoints) so the
                 * orbit is converged to well below the edge time and both
                 * agree with each other. PSS_FREERUN_POLISH_TRTOL
                 * (default 1; 0 disables) overrides the trtol part. */
                {
                    double polish_trtol = getenv("PSS_FREERUN_POLISH_TRTOL") ?
                        atof(getenv("PSS_FREERUN_POLISH_TRTOL")) : 1.0;
                    if (polish_trtol > 0.0 && polish_trtol < ckt->CKTtrtol)
                        ckt->CKTtrtol = polish_trtol;
                    rec_points = (job->PACrecordPoints > 0) ? job->PACrecordPoints :
                                 ((job->PSSpoints > 0) ? job->PSSpoints : 200);
                }
                x0_polished = TMALLOC(double, n);
                for (k = 0; k < n; k++) x0_polished[k] = x0[k];
                for (pk = 0; pk < polish_iters; pk++) {
                    double *xnext = TMALLOC(double, n);
                    double Tnext = 0.0;
                    if (rec_points > 0) {   /* same max step the recording pass will use */
                        /* the period a solver measures carries an O(h^2)
                         * max-step error that trtol does not see: 5.1e-11 s
                         * at T/200, 6e-13 s at T/2000 (docs/pss-pac.md card
                         * #517), and the answer moves ~1.5 dB per 1e-9 s --
                         * so the polish uses at least T/5000. */
                        long pp = getenv("PSS_FREERUN_POLISH_POINTS") ?
                            atol(getenv("PSS_FREERUN_POLISH_POINTS")) :
                            (rec_points > 5000 ? rec_points : 5000);
                        ckt->CKTmaxStep = T_pac / (double) pp;
                        ckt->CKTstep    = ckt->CKTmaxStep;
                        ckt->CKTdelmin  = 1e-9 * ckt->CKTmaxStep;
                    }
                    if (!xnext || pss_run_until_trigger(ckt, &map, x0_polished,
                                                         node_num, trigval, dir,
                                                         max_period, xnext,
                                                         &Tnext, NULL, NULL)) {
                        if (xnext) tfree(xnext);
                        polish_ok = 0;
                        break;
                    }
                    for (k = 0; k < n; k++) x0_polished[k] = xnext[k];
                    T_pac = Tnext;
                    tfree(xnext);
                }
                if (polish_ok) {
                    fprintf(stderr, "PSS-FREERUN: polished x0/T for .pac/"
                            ".pstb with %d extra plain period-map "
                            "(successive-substitution, no FD Jacobian) "
                            "iteration(s) (card #503 THIRD pass): "
                            "T %.12g -> %.12g\n", polish_iters, T, T_pac);
                    x0_pac = x0_polished;
                } else {
                    fprintf(stderr, "PSS-FREERUN: WARNING -- polish "
                            "iteration lost the trigger crossing; using "
                            "the plain Newton-converged x0/T for .pac/"
                            ".pstb instead (unpolished).\n");
                    T_pac = T;
                }
            }
        }
        /* card #503 FOURTH pass (maintainer decision: ungate INSIDE the
         * validated band, refuse outside it by default): the THIRD pass's
         * x0/T polish (above) was validated to meet this project's usual
         * 0.5dB/3deg target against tools/pac_perturbed.py --freerun at
         * 20/60/100kHz (<=23.7% of R5's fsw) but NOT yet at 150kHz+
         * (35%+) -- see docs/pss-pac.md's card #503 THIRD/FOURTH pass
         * sections. 0.25 x fsw (using the POLISHED T_pac, the actual
         * period this analysis uses) is the maintainer's chosen cutoff:
         * inside the measured-clean range (comfortably above the 23.7%
         * point actually validated, comfortably below the 35% point that
         * was not), consistent with this file's other band-guard
         * conventions (patch 0046's saltation guard, card #493's now-
         * superseded 10% guard) of drawing the line with margin rather
         * than exactly at the last point measured clean. Checked against
         * the SWEEP'S OWN TOP frequency (job->PACfStop) -- refusing the
         * WHOLE analysis when ANY requested point could be unvalidated,
         * simpler and no less honest than a per-point skip (a deck that
         * wants SOME points below the band and some above should split
         * into two .pac cards, same as any other sweep-range constraint
         * this project's tools already impose). Never silent either way:
         * a clear stderr line whichever branch is taken, same "no
         * silently-unvalidated number" convention as every other guard in
         * this file. */
        {
            double fsw_pac = (T_pac > 0.0) ? 1.0 / T_pac : 0.0;
            double band_limit = 0.25 * fsw_pac;
            int experimental = (getenv("PSS_FREERUN_PAC_EXPERIMENTAL") != NULL);
            int above_band = (fsw_pac <= 0.0) || (job->PACfStop > band_limit);

            if (above_band && !experimental) {
                fprintf(stderr, "PSS-FREERUN: REFUSED -- the requested .pac/"
                        ".pstb sweep's top frequency (%.6g Hz) exceeds the "
                        "VALIDATED BAND for a free-running (trigger=) orbit "
                        "(0.25 x fsw=%.6g Hz = %.6g Hz, fsw measured from "
                        "this orbit's own polished period) -- validated to "
                        "0.5dB/3deg against an independent reference "
                        "(tools/pac_perturbed.py --freerun) up to ~23.7%% "
                        "of fsw on R5, not yet above that (card #503 THIRD "
                        "pass; docs/pss-pac.md). Not running -- this is NOT "
                        "a silently-unvalidated number. Set "
                        "PSS_FREERUN_PAC_EXPERIMENTAL=1 to run anyway "
                        "(unvalidated above the band; not for production "
                        "use), or lower the sweep's own top frequency to "
                        "stay inside the band.\n",
                        job->PACfStop, band_limit, band_limit);
            } else {
            if (above_band && experimental)
                fprintf(stderr, "PSS-FREERUN: PSS_FREERUN_PAC_EXPERIMENTAL=1 -- "
                        "running .pac/.pstb on this free-running orbit ABOVE "
                        "the validated 0.25 x fsw=%.6g Hz band (sweep top "
                        "%.6g Hz) -- UNVALIDATED there (root cause of the "
                        "150kHz+/35%%+ disagreement not identified; card "
                        "#503 THIRD pass). Not for production use -- see "
                        "docs/pss-pac.md's free-running section.\n",
                        band_limit, job->PACfStop);
            else
                fprintf(stderr, "PSS-FREERUN: running .pac/.pstb on this "
                        "free-running orbit (validated to 0.5dB/3deg "
                        "against an independent reference up to 0.25 x "
                        "fsw=%.6g Hz on R5; card #503 THIRD/FOURTH pass, "
                        "docs/pss-pac.md).\n", band_limit);
        if (fT && x_prev &&
            !pss_run_until_trigger(ckt, &map, x0_pac, node_num, trigval, dir,
                                    max_period, x1trial, &T_meas2,
                                    x_prev, &t_prev)) {
            double dt = T_meas2 - t_prev;
            /* card #493: root-cause investigation of the 150kHz (35% of
             * R5's fsw) disagreement flagged in the comment above this
             * block. PSS_FREERUN_FT_DEBUG=1 on the original backward-
             * difference below showed dt itself is already fine (~15.8ps
             * on R5, not the outer Points-grid's ~11.8ns) -- but the Csw
             * (switch-node) row still swings ~0.7V over that window
             * (~4.5e10 V/s average slope): the conductance transition
             * this node rides on is driven by q's own tanh(ksw*(q-0.5))
             * argument crossing zero over roughly a ~130ps window (from
             * fT[q]'s own measured ~2.5e8 V/s and ksw=30's 1/ksw~=0.033
             * transition width), so a 15.8ps BACKWARD-difference window
             * (positioned by wherever the outer probe run's own adaptive
             * step happened to land, not centred on the crossing) is only
             * ~8x smaller than that -- resolving the local slope with a
             * plausible O(1) curvature error, not merely fine sampling
             * noise. PAC_EXACT_JACOBIAN=1 (pacshoot.c's OWN documented
             * accuracy floor, unrelated to this file) was cross-checked
             * on the SAME point and changed the result by <1e-4 relative
             * -- ruling out that floor as a meaningful contributor here.
             *
             * FIX: replace the backward difference (whatever the outer
             * probe run's own bracketing step happened to give) with a
             * FORWARD difference of a fixed, small, controlled width
             * PSS_FREERUN_FT_DT (default 1e-13s, ~130x smaller than the
             * conductance transition width above and ~1e5x smaller than
             * the old dt) run from the converged crossing point x_T
             * itself -- f+ = f- (design ruling #447 A, no true event) so
             * evaluating the "+" side is exactly as valid as the "-" side
             * this file's monodromy FD columns already use, and a forward
             * probe from x_T (already the exact converged crossing state)
             * needs no re-location of the crossing the way a backward
             * difference does. PSS_FREERUN_FT_BACKWARD=1 keeps the OLD
             * (pre-card-#493) backward-difference behaviour for direct
             * A/B comparison; unset (default whenever
             * PSS_FREERUN_PAC_EXPERIMENTAL=1 is set at all, itself still
             * required) uses this new forward difference. Both remain
             * behind PSS_FREERUN_PAC_EXPERIMENTAL -- this does not change
             * default (gated) behaviour either way. */
            if (dt > 0.0 && !getenv("PSS_FREERUN_FT_BACKWARD")) {
                double dt_fine_req = getenv("PSS_FREERUN_FT_DT") ?
                    atof(getenv("PSS_FREERUN_FT_DT")) : 1e-13;
                double *x_fine = TMALLOC(double, n);
                int fine_err;
                pss_state_set_ic(&map, x1trial);
                fine_err = pss_run_period(ckt, dt_fine_req,
                                          "PSS freerun (f(T) forward "
                                          "difference, card #493)");
                if (!fine_err && x_fine) {
                    pss_state_get(&map, ckt, x_fine);
                    for (k = 0; k < n; k++)
                        fT[k] = (x_fine[k] - x1trial[k]) / dt_fine_req;
                    if (getenv("PSS_FREERUN_FT_DEBUG"))
                        fprintf(stderr, "DEBUG (card #493 forward diff) "
                                "dt_fine=%.6g (old backward dt=%.6g)\n",
                                dt_fine_req, dt);
                    dt = dt_fine_req;   /* only used below as the >0 gate */
                } else {
                    fprintf(stderr, "PSS-FREERUN: WARNING -- card #493's "
                            "fine forward f(T) probe failed (dt=%.3g s); "
                            "falling back to the old backward-difference "
                            "estimate.\n", dt_fine_req);
                    for (k = 0; k < n; k++)
                        fT[k] = (x1trial[k] - x_prev[k]) / dt;
                }
                if (x_fine) tfree(x_fine);
            } else if (dt > 0.0) {
                for (k = 0; k < n; k++)
                    fT[k] = (x1trial[k] - x_prev[k]) / dt;
            }
            if (dt > 0.0) {
                int pac_err;
                fprintf(stdout, "PSS-FREERUN: .pac/.pstb period-variation "
                        "term: f(T) estimated from a %.6g s %s "
                        "difference at the trigger crossing (card #489 "
                        "second pass / card #493 -- see docs/pss-pac.md).\n",
                        dt, getenv("PSS_FREERUN_FT_BACKWARD") ? "backward" : "forward");
                if (getenv("PSS_FREERUN_FT_DEBUG")) {
                    int kk;
                    fprintf(stderr, "DEBUG T_meas2=%.12g t_prev=%.12g dt=%.6g\n",
                            T_meas2, t_prev, dt);
                    for (kk = 0; kk < n; kk++)
                        fprintf(stderr, "DEBUG fT[%d]=%.9g x_T=%.9g x_prev=%.9g\n",
                                kk, fT[kk], x1trial[kk], x_prev[kk]);
                }
                if (getenv("PSS_FREERUN_TRACE") || getenv("PSS_FREERUN_PERT_T") ||
                    getenv("PSS_FREERUN_PERT_X0") || getenv("PSS_FREERUN_PERT_FT")) {
                    /* card #517 diagnostics: cross-machine sensitivity study.
                     * PERT_X0="idx:rel", PERT_FT="idx:rel", PERT_T="rel". */
                    const char *e;
                    int kk, id;
                    double rl;
                    if ((e = getenv("PSS_FREERUN_PERT_X0")) && sscanf(e, "%d:%lf", &id, &rl) == 2 && id >= 0 && id < n)
                        x0_pac[id] *= (1.0 + rl);
                    if ((e = getenv("PSS_FREERUN_PERT_FT")) && sscanf(e, "%d:%lf", &id, &rl) == 2 && id >= 0 && id < n)
                        fT[id] *= (1.0 + rl);
                    if ((e = getenv("PSS_FREERUN_PERT_T")) && sscanf(e, "%lf", &rl) == 1)
                        T_pac *= (1.0 + rl);
                    fprintf(stderr, "TRACE517 T_pac=%.17g\n", T_pac);
                    for (kk = 0; kk < n; kk++)
                        fprintf(stderr, "TRACE517 k=%d x0=%.17g fT=%.17g\n", kk, x0_pac[kk], fT[kk]);
                }
                pac_err = pac_pstb_run_ex(ckt, job, x0_pac, n, T_pac, fT, node_num);
                if (reltol_tightened) {   /* card #517: profile held through the recording passes */
                    ckt->CKTreltol = save_reltol_polish;
                    reltol_tightened = 0;
                }
                if (pac_err)
                    fprintf(stderr, "PSS-FREERUN: .pac/.pstb on the "
                            "free-running orbit did not complete (see "
                            "pacshoot.c's own diagnostic above); the "
                            "PSS-FREERUN orbit/Floquet result above this "
                            "is still valid on its own.\n");
            } else {
                fprintf(stderr, "PSS-FREERUN: WARNING -- could not estimate "
                        "f(T) at the trigger crossing (degenerate dt=%.3g); "
                        ".pac/.pstb skipped.\n", dt);
            }
        } else {
            fprintf(stderr, "PSS-FREERUN: WARNING -- could not re-locate "
                    "the trigger crossing from the converged x0 to estimate "
                    "the period-variation term; .pac/.pstb skipped.\n");
        }
            }   /* end of "inside the band, or experimental" branch */
        }   /* end of band-check block */
        ckt->CKTtrtol   = save_trtol_polish;       /* card #517 */
        ckt->CKTmaxStep = save_maxstep_polish;
        ckt->CKTstep    = save_step_polish;
        ckt->CKTdelmin  = save_delmin_polish;
        if (reltol_tightened)
            ckt->CKTreltol = save_reltol_polish;   /* catch-all: any early-exit
                                                     * branch above that skipped
                                                     * the restore right before
                                                     * pac_pstb_run_ex() */
        if (fT) tfree(fT);
        if (x_prev) tfree(x_prev);
        if (x0_polished) tfree(x0_polished);
    }

    err = OK;

 fr_cleanup:
    ckt->CKTreltol  = save_reltol;
    ckt->CKTabstol  = save_abstol;
    ckt->CKTvoltTol = save_volttol;
    PSStrigActive = 0;
    pss_trig_free();
    PSSrippleRecordActive = 0;
    pss_ripple_free();
    if (x0) tfree(x0);
    if (x1) tfree(x1);
    if (F) tfree(F);
    if (M) tfree(M);
    if (delta) tfree(delta);
    if (xtrial) tfree(xtrial);
    if (x1trial) tfree(x1trial);
    if (scale) tfree(scale);
    if (ripple) tfree(ripple);
    if (A) tfree(A);
    if (b) tfree(b);
    if (eig) tfree(eig);
    if (map.refs) tfree(map.refs);
    return err;

 fr_fail:
    err = E_ITERLIM;
    goto fr_cleanup;
}
