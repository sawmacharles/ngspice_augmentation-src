/**********
card #513 (plan T9.11, round 11 of transient acceleration): IMPLICIT,
Newton-corrected envelope following inside the fork, for long switching-
converter transients (soft-start / slow load-line ramps spanning
thousands of switching periods).

    .tranenv Fguess Points SC_iter Steady_coeff Tstop Mmax
             [envtol=<tol>] [phase=<0..1>] [maptol=<reltol>] [uic]

Fguess = clock frequency (the period T), Points = integration points per
period (sets the max step T/Points of every period map), SC_iter = Newton
iteration budget per envelope step (capped at 8), Steady_coeff = Newton
tolerance (scaled step / residual, default 1e-4), Tstop = total simulated
time, Mmax = largest jump in periods, envtol = local envelope error
tolerance (relative to the state scale, default 5e-4: 0.1-0.4 % cycle-
average error on the validation decks; 2e-3 gives 0.2-1 %), phase = sampling phase (period maps run
[k*T+phase*T, (k+1)*T+phase*T]; default 0; pick a quiet part of the cycle
for decks with hard edges at the clock boundary), maptol = transient reltol
of every period map (default 1e-5, never looser than the deck's reltol;
0 = the deck's reltol; see "Map tolerance" below), uic = start from the
deck's IC=/.ic (the only supported start; assumed).  A separate directive,
not an extension of `.tran`, so no existing `.tran` deck can be affected.
Needs a driven (clocked) circuit and at least one native C or L.

Round 10 (card #512, docs/envelope.md) built an explicit projective jump
in Python over savestate/loadstate; it failed the accuracy bar (stale
CKTstates, 2.5-5.5 %) and gave 1.2x (per-hop process overhead).  This is
Kundert/White/Sangiovanni-Vincentelli (ICCAD 1988) done natively.

State: the reduced state of `.pss method=shoot` -- every native capacitor
voltage and inductor current (n values).  Every period map F(x; t0) is a
fresh DCtran() restart from those values (uic), at an ABSOLUTE start time
t0 (dctran.c's PSSenvActive/PSSenvPeriodT0 generalisation of "a restart
starts at t=0"), so slow sources (soft-start ramps) evaluate correctly and
the clock is in phase at every landing.  Nothing from a previous map (in
particular no CKTstates history) survives into the next, so there is no
stale state to leave behind at a jump.

Envelope step (M an integer number of periods): solve for y = x(t_n+M*T)

    y - x_n = M * ( F(y; t_n + M*T) - y )                        (1)

by Newton, Jm = (1+M)*I - M*Phi with Phi ~ dF/dy the monodromy, finite-
differenced from n free period runs (step 2e-3 of the state scale; larger
steps are essential -- the transient solver's own convergence noise in F
divided by a tiny step gave garbage Jacobians) and CACHED across steps and
across M (chord Newton, one period run per iteration).  The landing period
and one confirmation period are then run for real and emitted.

Map tolerance (card #719).  Newton, the curvature predictor and the error
estimate all difference period maps: the implicit step multiplies the noise
of F(y) - y by M, and the curvature estimate divides period-to-period
differences.  At the default transient reltol (1e-3) a closed-loop converter
with sub-nanosecond switch edges (refs r3d: tanh switches, 1-2 ns
transitions) has a period map whose step-control noise is far above that
level: on the soft-start deck the sampled inductor current jitters by
0.5-3 mA from period to period where a reltol=1e-7 run shows 1e-6 A, finite-
difference Jacobian columns change by 40x with the step (1e-8 .. 1e-2 of the
state), and no jump was ever accepted.  The same map at reltol 1e-5 is smooth
to 1e-7 A (columns agree from 1e-6 to 1e-3), jumps are accepted from the
first hundred periods on, and the run is both faster and closer to the
reltol=1e-7 reference than a plain .tran at the default reltol.  So every
period map (bootstrap, Newton, FD, landing, plain) runs at reltol =
min(deck reltol, maptol).  A fixed step grid (no LTE control, same steps in
every map) was tried as well and dropped: it makes the map independent of
the step decisions, but a 10 ns grid on a 1 ns switch edge turns the map
into a steep staircase and the chain into a quantisation limit cycle.
ENV_MAPTOL_OFF=1 restores the deck's reltol (red proof).

Error control and M selection, all on the SLOW states only (Phi_kk >= 0.9:
the integrator-like states; fast re-slaved states -- switch node, dead-time
network, LC ripple -- are handled exactly by the implicit step and their
period-to-period jitter would otherwise pin M at 2-3):
  * M is capped by the curvature predictor  M <= sqrt(2*0.3*envtol*scale/
    |c|), c = last change of the per-period slope of the slow states;
  * after Newton, |y - (x_n + M*(x_n - x_n-1))| over the slow states must be
    <= envtol*scale, else M shrinks by sqrt(envtol/err) and the step is
    refused (plain real period instead);
  * after an accepted step the next M is M * min(2, 0.9 sqrt(envtol/err))
    (card #719; a plain doubling alternated accept/refuse), capped by Mmax;
  * right after a landing (y and F(y)) the curvature is taken between the
    jump's secant slope (y - x_n)/M and the landing period's slope, so the
    next jump needs no extra real period (card #719; ENV_NCHECK, default 1,
    adds confirmation periods);
  * a jump never lands within 2 periods of a known slow-source corner (any
    native PWL V/I source: soft-start knees, load steps); the corner is
    crossed by real periods.  A hidden (non-PWL) persistent event inside a
    skipped span shows up in F(y) at the landing and is refused by the
    error estimate; a transient wholly inside a skipped span is invisible
    to every envelope method (docs/envelope.md);
  * Newton failure or a period-map error -> plain real period, then
    exponential back-off (4, 8 ... 128 real periods) before the next
    attempt, so a deck whose map is too rough for Newton (e.g. a closed-
    loop deck run with maptol=0 at reltol 1e-3: Phi entries change by
    10-100x with the FD step) degrades to a plain period-by-period chain
    instead of paying for a failed solve every period;
  * a history-relative smoothness gate (delta vs median of the last 8
    per-period deltas) forces M=1 on a sudden change.

Output: every REAL simulated period (bootstrap, landings, confirmations,
plain periods) is appended to ONE ordinary time-domain plot (title of the
first period), so wrdata/plot/meas work as after `.tran`; time simply has a
gap across every skipped span.  Newton/FD scratch runs create no plot
(PSSenvOutMode 1; an earlier version made one plot per scratch run and
reached 20 GB).  A second small plot "Envelope Samples" holds the n states
at each real period boundary (not created when a `save` list leaves none of
its vectors; a warning says so).  Nothing is interpolated or invented.

Env vars (debug/tuning, no effect unless set): ENV_DEBUG=1 (decisions on
stderr; 4 also the curvature limit per period and the state that sets it),
ENV_NCHECK (real periods after a landing, default 1), ENV_GROWTH_TOL
(smoothness gate factor, default 5), ENV_MAPTOL_OFF=1 (period maps at the
deck's reltol), ENV_SCAN=<time> (map smoothness diagnostic, see env_scan()).

Not done: analytic (discretise-then-linearise) monodromy (pssshoot.c has
one; it is documented inaccurate near high-gain comparators, and was tried
here and dropped), fast-state reduction of the Newton system, Broyden
updates, non-PWL slow-source corner detection, restart from a non-uic
operating point, autonomous circuits.
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/pssdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/fteext.h"

#include "cap/capdefs.h"
#include "vsrc/vsrcdefs.h"
#include "isrc/isrcdefs.h"
#include "ind/inddefs.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pssshoot.h"

/* ---------------------------------------------------------------- */
/* State map: deliberately duplicated from pssshoot.c's own (private,
 * static) pss_build_state_map()/pss_state_get_ic()/pss_state_set_ic()/
 * pss_state_get() rather than exported from there -- see this file's
 * header comment ("Why round 10's landing kick does NOT reappear here")
 * for the reduced-state design choice itself; this duplication is what
 * keeps pssshoot.c's own diff for this card down to the single
 * PSSbreakPeriodT0 generalisation every existing breakpoint-forcing
 * caller needs regardless (see pssshoot.h), rather than opening a new
 * exported API surface on top of that file's own internals. */

typedef struct {
    int isCap;
    union { CAPinstance *cap; INDinstance *ind; } dev;
    char name[64];
} EnvStateRef;

typedef struct {
    int n;
    EnvStateRef *refs;
} EnvStateMap;

static int
env_build_state_map(CKTcircuit *ckt, EnvStateMap *map)
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

    map->refs = TMALLOC(EnvStateRef, n);
    if (!map->refs)
        return E_NOMEM;

    n = 0;
    if (captype >= 0)
        for (cm = (CAPmodel *) ckt->CKThead[captype]; cm; cm = CAPnextModel(cm))
            for (ci = CAPinstances(cm); ci; ci = CAPnextInstance(ci)) {
                map->refs[n].isCap = 1;
                map->refs[n].dev.cap = ci;
                snprintf(map->refs[n].name, sizeof map->refs[n].name, "C:%s", ci->CAPname);
                n++;
            }
    if (indtype >= 0)
        for (im = (INDmodel *) ckt->CKThead[indtype]; im; im = INDnextModel(im))
            for (ii = INDinstances(im); ii; ii = INDnextInstance(ii)) {
                map->refs[n].isCap = 0;
                map->refs[n].dev.ind = ii;
                snprintf(map->refs[n].name, sizeof map->refs[n].name, "L:%s", ii->INDname);
                n++;
            }
    return OK;
}

static void
env_map_free(EnvStateMap *map)
{
    if (map->refs) tfree(map->refs);
    map->refs = NULL;
    map->n = 0;
}

/* Deck's own IC=/. ic default (typically 0) -- the t=0 bootstrap seed,
 * same quantity pss_state_get_ic() reads. */
static void
env_state_get_ic(EnvStateMap *map, double *x)
{
    int k;
    for (k = 0; k < map->n; k++)
        x[k] = map->refs[k].isCap ? map->refs[k].dev.cap->CAPinitCond
                                   : map->refs[k].dev.ind->INDinitCond;
}

static void
env_state_set_ic(EnvStateMap *map, const double *x)
{
    int k;
    for (k = 0; k < map->n; k++) {
        if (map->refs[k].isCap)
            map->refs[k].dev.cap->CAPinitCond = x[k];
        else
            map->refs[k].dev.ind->INDinitCond = x[k];
    }
}

static double
env_state_value(const EnvStateRef *ref, CKTcircuit *ckt)
{
    if (ref->isCap) {
        CAPinstance *ci = ref->dev.cap;
        return ckt->CKTrhsOld[ci->CAPposNode] - ckt->CKTrhsOld[ci->CAPnegNode];
    } else {
        INDinstance *ii = ref->dev.ind;
        return ckt->CKTrhsOld[ii->INDbrEq];
    }
}

static void
env_state_get(EnvStateMap *map, CKTcircuit *ckt, double *x)
{
    int k;
    for (k = 0; k < map->n; k++)
        x[k] = env_state_value(&map->refs[k], ckt);
}

/* ---------------------------------------------------------------- */
/* Small dense linear solve (Gauss elimination, partial pivoting) --
 * deliberately duplicated from pssshoot.c's own pss_solve_dense() for the
 * same "keep this file self-contained" reason as the state map above; n
 * is never more than a handful of states here. A is n*n row-major and is
 * destroyed; returns 0 on success, 1 if singular to working precision. */
static int
env_solve_dense(double *A, double *b, int n, double *x)
{
    int i, j, k, piv;
    for (i = 0; i < n; i++) {
        double best = fabs(A[i * n + i]);
        piv = i;
        for (k = i + 1; k < n; k++) {
            double v = fabs(A[k * n + i]);
            if (v > best) { best = v; piv = k; }
        }
        if (best < 1e-300)
            return 1;
        if (piv != i) {
            for (j = 0; j < n; j++) {
                double t = A[i * n + j]; A[i * n + j] = A[piv * n + j]; A[piv * n + j] = t;
            }
            { double t = b[i]; b[i] = b[piv]; b[piv] = t; }
        }
        for (k = i + 1; k < n; k++) {
            double f = A[k * n + i] / A[i * n + i];
            if (f == 0.0) continue;
            for (j = i; j < n; j++)
                A[k * n + j] -= f * A[i * n + j];
            b[k] -= f * b[i];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (j = i + 1; j < n; j++)
            s -= A[i * n + j] * x[j];
        if (fabs(A[i * n + i]) < 1e-300)
            return 1;
        x[i] = s / A[i * n + i];
    }
    return 0;
}

/* ---------------------------------------------------------------- */
/* One-period map, generalised to an arbitrary ABSOLUTE start time t0 --
 * see pssshoot.h's PSSenvActive/PSSenvPeriodT0/PSSbreakPeriodT0 and
 * dctran.c's own (narrowly gated) restart-branch generalisation. This is
 * pssshoot.c's pss_run_period()/pss_run_period_ex() pattern, duplicated
 * and extended here rather than exported from there (see this file's
 * header comment), plus the t0 handling those two functions never needed.
 *
 * len is the run's own LENGTH (T for a single period-map evaluation used
 * by Newton; K*T for a multi-period "kept output" segment -- this engine
 * only ever passes T, see the header comment on why batching was not
 * attempted this round). record=1 defines this run's own accepted-time
 * grid (pssshoot.c's g_break_times, via the exported PSSbreakReset()/
 * PSSbreakRecordActive); force=1 makes this run retrace a PREVIOUSLY
 * recorded grid (the SAME FD-Jacobian noise-removal reasoning
 * pssshoot.h's own header documents for `.pss method=shoot`). Every call
 * opens and closes its own ngspice plot -- see the header comment
 * "Output" section for why this engine does not try to avoid that. */
/* Reset every native V/I source's own per-instance "next breakpoint"
 * bookkeeping (VSRCbreak_time/ISRCbreak_time, set to -1 by vsrcset.c at
 * setup = "schedule the next corner from CKTtime") before every period
 * map. That field otherwise persists across DCtran() restarts, and this
 * engine restarts at arbitrary absolute times that are not monotonic
 * (Newton evaluates the map at t_land, then a refused jump re-runs at
 * t_n): a stale value either suppresses a PULSE/PWL source's own edges
 * or makes CKTsetBreak() panic ("breakpoint in the past"). */
static void
env_reset_source_breaks(CKTcircuit *ckt)
{
    int vt = CKTtypelook("Vsource");
    int it = CKTtypelook("Isource");
    if (vt >= 0) {
        VSRCmodel *m; VSRCinstance *h;
        for (m = (VSRCmodel *) ckt->CKThead[vt]; m; m = VSRCnextModel(m))
            for (h = VSRCinstances(m); h; h = VSRCnextInstance(h))
                h->VSRCbreak_time = -1.0;
    }
    if (it >= 0) {
        ISRCmodel *m; ISRCinstance *h;
        for (m = (ISRCmodel *) ckt->CKThead[it]; m; m = ISRCnextModel(m))
            for (h = ISRCinstances(m); h; h = ISRCnextInstance(h))
                h->ISRCbreak_time = -1.0;
    }
}

static double env_scale(double v, double floor_v);

/* card #719: count of the time steps the period maps rejected (diagnostic). */
static long g_env_rejected = 0;

static int
env_run_period_at(CKTcircuit *ckt, EnvStateMap *map, const double *x_in,
                   double t0, double len, const char *jobname,
                   int record, int force)
{
    JOB *saveJob = ckt->CKTcurJob;
    TRANan localjob;
    char namebuf[192];
    int err;

    env_state_set_ic(map, x_in);
    env_reset_source_breaks(ckt);

    if (getenv("PSS_SHOOT_NO_FORCE") != NULL)
        force = 0;   /* same debug/isolation switch pssshoot.c's own period runner honours */

    if (record)
        PSSbreakReset();
    PSSbreakRecordActive  = record;
    PSSbreakInstallActive = force;
    PSSbreakPeriodT0       = t0;
    PSSenvActive            = 1;
    PSSenvOutMode           = (strncmp(jobname, "ENV Newton", 10) == 0) ? 1 : 2;
    PSSenvPeriodT0           = t0;

    memset(&localjob, 0, sizeof localjob);
    strncpy(namebuf, jobname, sizeof namebuf - 1);
    namebuf[sizeof namebuf - 1] = '\0';
    localjob.JOBtype = saveJob->JOBtype;
    localjob.JOBname = namebuf;

    ckt->CKTcurJob    = (JOB *) &localjob;
    ckt->CKTfinalTime = t0 + len;
    ckt->CKTinitTime  = t0;
    ckt->CKTtime      = t0;
    ckt->CKTmode     |= MODEUIC;   /* envelope following always restarts each period from x */

    {
        long rej0 = ckt->CKTstat->STATrejected;
        err = DCtran(ckt, 1);
        g_env_rejected += ckt->CKTstat->STATrejected - rej0;
    }

    if (PSSenvOutMode == 2 && localjob.TRANplot)
        PSSenvKeptPlot = localjob.TRANplot;   /* stays open: one stitched plot */
    else if (localjob.TRANplot)
        SPfrontEnd->OUTendPlot(localjob.TRANplot);

    ckt->CKTcurJob = saveJob;
    PSSbreakRecordActive  = 0;
    PSSbreakInstallActive = 0;
    PSSenvActive            = 0;
    PSSenvOutMode           = 0;
    return err;
}

/* card #719, diagnostic (ENV_SCAN=<time>): at the first real period start at
 * or after <time>, print finite-difference columns dF/dx_j of the one-period
 * map for relative steps 1e-8 .. 1e-2 of each state's scale, and a line scan
 * F(x + s e_j) for 41 values of s in [-1e-4, 1e-4] scale (ENV_SCAN_LINE=j,
 * default all j), on stderr.  A smooth map gives columns that agree for small
 * steps; a rough one gives columns growing like noise/step. */
static void
env_scan(CKTcircuit *ckt, EnvStateMap *map, int n, const double *x, double t0,
         double T, double floor_v)
{
    double *F0 = TMALLOC(double, n), *Fp = TMALLOC(double, n), *xp = TMALLOC(double, n);
    static const double eps[] = { 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2 };
    int j, k, e, i;
    const char *ls = getenv("ENV_SCAN_LINE");
    int jline = ls ? atoi(ls) : -1;

    if (!F0 || !Fp || !xp) goto out;
    if (env_run_period_at(ckt, map, x, t0, T, "ENV Newton (scan)", 0, 0)) goto out;
    env_state_get(map, ckt, F0);
    fprintf(stderr, "ENVSCAN t0=%.12g n=%d reltol=%g\n", t0, n, ckt->CKTreltol);
    for (k = 0; k < n; k++)
        fprintf(stderr, "ENVSCAN state %d %s x=%.15g F=%.15g scale=%.6g\n", k,
                map->refs[k].name, x[k], F0[k], env_scale(x[k], floor_v));
    for (j = 0; j < n; j++) {
        double sc = env_scale(x[j], floor_v);
        for (e = 0; e < (int) (sizeof eps / sizeof eps[0]); e++) {
            double h = eps[e] * sc;
            for (k = 0; k < n; k++) xp[k] = x[k];
            xp[j] += h;
            if (env_run_period_at(ckt, map, xp, t0, T, "ENV Newton (scan)", 0, 0)) {
                fprintf(stderr, "ENVSCAN col %d eps=%g FAILED\n", j, eps[e]);
                continue;
            }
            env_state_get(map, ckt, Fp);
            fprintf(stderr, "ENVSCAN col %d eps=%g", j, eps[e]);
            for (k = 0; k < n; k++)
                fprintf(stderr, " %.9g", (Fp[k] - F0[k]) / h);
            fprintf(stderr, "\n");
        }
        if (jline >= 0 && jline != j) continue;
        for (i = -20; i <= 20; i++) {
            double s = 1e-4 * sc * (double) i / 20.0;
            for (k = 0; k < n; k++) xp[k] = x[k];
            xp[j] += s;
            if (env_run_period_at(ckt, map, xp, t0, T, "ENV Newton (scan)", 0, 0)) continue;
            env_state_get(map, ckt, Fp);
            fprintf(stderr, "ENVSCAN line %d s=%.6g", j, s);
            for (k = 0; k < n; k++)
                fprintf(stderr, " %.15g", Fp[k]);
            fprintf(stderr, "\n");
        }
    }
out:
    if (F0) tfree(F0);
    if (Fp) tfree(Fp);
    if (xp) tfree(xp);
}

/* ---------------------------------------------------------------- */
/* "Envelope Samples" plot: a small, separate (t, x[0..n)) trace recording
 * EVERY real sample this engine ever computes (bootstrap, landing,
 * confirmation, passthrough) -- see the header comment's "Output"
 * section. Built directly from SPfrontEnd's pointwise-plot primitives
 * (the SAME ones CKTdump()/DCtran() use internally), not through
 * DCtran(), since this plot's own width (n states) and cadence (one row
 * per REAL period, potentially spanning a skipped multi-period jump) do
 * not match any single ordinary transient run. */
static runDesc *
env_open_samples_plot(CKTcircuit *ckt, JOB *job, EnvStateMap *map)
{
    IFuid timeUid;
    IFuid *nameList;
    int i, err;
    runDesc *plot = NULL;

    nameList = TMALLOC(IFuid, map->n);
    if (!nameList)
        return NULL;
    for (i = 0; i < map->n; i++)
        SPfrontEnd->IFnewUid(ckt, &nameList[i], NULL, map->refs[i].name, UID_OTHER, NULL);
    SPfrontEnd->IFnewUid(ckt, &timeUid, NULL, "time", UID_OTHER, NULL);

    err = SPfrontEnd->OUTpBeginPlot(ckt, job, "Envelope Samples (card #513)",
                                     timeUid, IF_REAL,
                                     map->n, nameList, IF_REAL, &plot);
    tfree(nameList);
    if (err)
        return NULL;
    return plot;
}

static void
env_emit_sample(runDesc *plot, double t, const double *x, int n)
{
    IFvalue refData, valData;
    if (!plot)
        return;
    refData.rValue = t;
    valData.v.numValue = n;
    valData.v.vec.rVec = (double *) x;
    SPfrontEnd->OUTpData(plot, &refData, &valData);
}

/* ---------------------------------------------------------------- */
/* Smoothness gate / local-error-estimate helpers -- max-abs-relative
 * norm (relative to a floor, so a near-zero state doesn't blow up the
 * ratio), and an 8-entry rolling median of ACCEPTED real per-period
 * deltas (round 10's own history-relative smoothness test, see the
 * header comment; docs/envelope.md's own finding on why a fixed-
 * relative-to-magnitude test is unsafe applies identically here). */

static double
env_maxrel(const double *a, const double *b, int n, double floor_v)
{
    int k;
    double worst = 0.0;
    for (k = 0; k < n; k++) {
        double scale = fabs(b[k]);
        if (scale < floor_v) scale = floor_v;
        double e = fabs(a[k] - b[k]) / scale;
        if (e > worst) worst = e;
    }
    return worst;
}

/* Newton residual norm: max_k |R[k]| / max(|y[k]|, floor) -- R is not
 * itself "a value relative to some other value" the way env_maxrel's two
 * arguments are, so this is its own small helper rather than a call to
 * env_maxrel with a synthetic second argument. */
static double
env_residual_norm(const double *R, const double *y, int n)
{
    int k;
    double worst = 0.0;
    for (k = 0; k < n; k++) {
        double scale = fabs(y[k]);
        double e;
        if (scale < 1e-3) scale = 1e-3;
        e = fabs(R[k]) / scale;
        if (e > worst) worst = e;
    }
    return worst;
}

#define ENV_HIST_N 8

typedef struct {
    double buf[ENV_HIST_N];
    int n, next;
} EnvHist;

static void
env_hist_push(EnvHist *h, double v)
{
    h->buf[h->next] = v;
    h->next = (h->next + 1) % ENV_HIST_N;
    if (h->n < ENV_HIST_N) h->n++;
}

static double
env_hist_median(const EnvHist *h)
{
    double tmp[ENV_HIST_N];
    int i, j, n = h->n;
    if (n == 0) return 0.0;
    for (i = 0; i < n; i++) tmp[i] = h->buf[i];
    /* insertion sort -- n <= 8 */
    for (i = 1; i < n; i++) {
        double key = tmp[i];
        j = i - 1;
        while (j >= 0 && tmp[j] > key) { tmp[j + 1] = tmp[j]; j--; }
        tmp[j + 1] = key;
    }
    return (n % 2) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
}


/* ---------------------------------------------------------------- */
/* Implicit envelope-step Newton solve (see the header comment, "The
 * algorithm").  Solves (1+M)*y - M*F(y; t_land) - x_n = 0, F the ordinary
 * (free, adaptive-step) one-period map started at absolute time t_land.
 *
 * Chord Newton: the monodromy Phi ~ dF/dy is finite-differenced from n
 * free period runs and CACHED across envelope steps and across M
 * (Jm = (1+M)I - M*Phi is rebuilt from it), so an iteration costs ONE
 * period run.  Phi is refreshed when it is older than ENV_JAC_MAXAGE
 * solves, when an iteration fails to halve the residual, or when a step
 * is rejected.  The FD step is a large one (2e-3 of the state scale):
 * an earlier design with tiny steps (1e-5) on forced grids gave garbage
 * entries because the transient solver's own convergence-tolerance noise
 * in F is divided by the step; analytic sensitivities (pssshoot.c) were
 * tried and are documented as inaccurate near high-gain comparator
 * nonlinearities (docs/pss-pac.md).
 *
 * Converged when the scaled step or the scaled residual falls below tol.
 * Returns 0 on convergence; *nmaps counts period runs. */

#define ENV_JAC_MAXAGE 12

typedef struct {
    double *Phi;    /* n*n, row-major */
    int valid;
    int age;
} EnvJac;

static double
env_scale(double v, double floor_v)
{
    double a = fabs(v);
    return a < floor_v ? floor_v : a;
}

/* max_k |R_k| / ((1+M) * scale_k): an error-in-y-equivalent norm */
static double
env_res_norm(const double *R, const double *y, int n, double M, double floor_v)
{
    int k;
    double worst = 0.0;
    for (k = 0; k < n; k++) {
        double e = fabs(R[k]) / ((1.0 + M) * env_scale(y[k], floor_v));
        if (e > worst) worst = e;
    }
    return worst;
}

/* Phi by forward FD of the free map about (y, F0) at t0; n period runs */
static int
env_phi_fd(CKTcircuit *ckt, EnvStateMap *map, int n, const double *y,
           const double *F0, double t0, double T, double floor_v,
           EnvJac *jac, long *nmaps)
{
    double *yp = TMALLOC(double, n), *Fp = TMALLOC(double, n);
    int j, k, err = 0;
    if (!yp || !Fp) { err = 1; goto out; }
    for (j = 0; j < n; j++) {
        double h = 2e-3 * env_scale(y[j], floor_v);
        for (k = 0; k < n; k++) yp[k] = y[k];
        yp[j] += h;
        err = env_run_period_at(ckt, map, yp, t0, T, "ENV Newton (FD probe)", 0, 0);
        (*nmaps)++;
        if (err) goto out;
        env_state_get(map, ckt, Fp);
        for (k = 0; k < n; k++) jac->Phi[k * n + j] = (Fp[k] - F0[k]) / h;
    }
    jac->valid = 1;
    jac->age = 0;
out:
    if (yp) tfree(yp);
    if (Fp) tfree(Fp);
    return err;
}

static int
env_newton_step(CKTcircuit *ckt, EnvStateMap *map, int n,
                const double *x_n, double M, double t_land, double T,
                double tol, int maxiter, double floor_v, EnvJac *jac,
                double *y_io, long *nmaps)
{
    double *F0 = TMALLOC(double, n), *R0 = TMALLOC(double, n);
    double *Jm = TMALLOC(double, (size_t) n * n), *delta = TMALLOC(double, n);
    double *yt = TMALLOC(double, n), *Ft = TMALLOC(double, n);
    double *Rt = TMALLOC(double, n), *rhs = TMALLOC(double, n);
    int iter, j, k, err, converged = 0, refreshed = 0;
    double res0 = 1e300;
    int debug = getenv("ENV_DEBUG") != NULL;

    if (!F0 || !R0 || !Jm || !delta || !yt || !Ft || !Rt || !rhs) {
        fprintf(stderr, "ENV: out of memory in Newton scratch allocation.\n");
        goto done;
    }
    err = env_run_period_at(ckt, map, y_io, t_land, T, "ENV Newton (map)", 0, 0);
    (*nmaps)++;
    if (err) goto done;
    env_state_get(map, ckt, F0);
    for (k = 0; k < n; k++) R0[k] = (1.0 + M) * y_io[k] - M * F0[k] - x_n[k];
    res0 = env_res_norm(R0, y_io, n, M, floor_v);

    for (iter = 0; iter < maxiter; iter++) {
        double alpha = 1.0, stepn;
        int bt, accepted = 0;

        if (debug)
            fprintf(stderr, "ENV Newton: iter=%d M=%.6g res=%.3e jac_age=%d\n",
                    iter, M, res0, jac->valid ? jac->age : -1);
        if (res0 < tol) { converged = 1; break; }

        if (!jac->valid || jac->age >= ENV_JAC_MAXAGE) {
            if (env_phi_fd(ckt, map, n, y_io, F0, t_land, T, floor_v, jac, nmaps)) goto done;
            refreshed = 1;
        }
        for (k = 0; k < n; k++)
            for (j = 0; j < n; j++)
                Jm[k * n + j] = ((k == j) ? (1.0 + M) : 0.0) - M * jac->Phi[k * n + j];
        for (k = 0; k < n; k++) rhs[k] = -R0[k];
        if (env_solve_dense(Jm, rhs, n, delta) != 0) { jac->valid = 0; goto done; }
        stepn = 0.0;
        for (k = 0; k < n; k++) {
            double e = fabs(delta[k]) / env_scale(y_io[k], floor_v);
            if (e > stepn) stepn = e;
        }

        for (bt = 0; bt < 4; bt++) {
            for (k = 0; k < n; k++) yt[k] = y_io[k] + alpha * delta[k];
            err = env_run_period_at(ckt, map, yt, t_land, T, "ENV Newton (map)", 0, 0);
            (*nmaps)++;
            if (!err) {
                double rt;
                env_state_get(map, ckt, Ft);
                for (k = 0; k < n; k++)
                    Rt[k] = (1.0 + M) * yt[k] - M * Ft[k] - x_n[k];
                rt = env_res_norm(Rt, yt, n, M, floor_v);
                if (rt < res0 || rt < tol) {
                    accepted = 1;
                    for (k = 0; k < n; k++) { y_io[k] = yt[k]; F0[k] = Ft[k]; R0[k] = Rt[k]; }
                    if (rt > 0.5 * res0 && rt >= tol) jac->age = ENV_JAC_MAXAGE;   /* chord too slow: refresh next iteration */
                    res0 = rt;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            if (!refreshed) { jac->age = ENV_JAC_MAXAGE; continue; }   /* stale Phi: refresh once and retry */
            if (debug) fprintf(stderr, "ENV Newton: step rejected M=%.6g\n", M);
            goto done;
        }
        if (res0 < tol || stepn * alpha < tol) { converged = 1; break; }
    }
    jac->age++;

done:
    if (F0) tfree(F0);
    if (R0) tfree(R0);
    if (Jm) tfree(Jm);
    if (delta) tfree(delta);
    if (yt) tfree(yt);
    if (Ft) tfree(Ft);
    if (Rt) tfree(Rt);
    if (rhs) tfree(rhs);
    return converged ? 0 : 1;
}

/* Slow-state mask: out[k] = v[k] for the states whose one-period map is
 * (nearly) the identity on its own coordinate, Phi_kk >= ENV_SLOW_DIAG
 * (the integrator-like states: they ARE the envelope), 0 for the fast,
 * re-slaved-every-period states (dead-time network, switch node, LC
 * ripple), which an implicit envelope step handles exactly and whose
 * period-to-period jitter otherwise dominates any predictor-error
 * estimate (docs/envelope.md: raw all-state curvature held M at 2-3 over
 * the whole ramp).  Falls back to all states if none qualify. */
#define ENV_SLOW_DIAG 0.9
static void
env_slow_mask(const double *Phi, int n, const double *v, double *out)
{
    int k, any = 0;
    for (k = 0; k < n; k++)
        if (Phi[k * n + k] >= ENV_SLOW_DIAG) any = 1;
    for (k = 0; k < n; k++)
        out[k] = (!any || Phi[k * n + k] >= ENV_SLOW_DIAG) ? v[k] : 0.0;
}

/* Monodromy seed before the first jump: one base run + n FD runs. */
static int
env_phi_at(CKTcircuit *ckt, EnvStateMap *map, int n, const double *y,
           double t0, double T, double floor_v, EnvJac *jac, long *nmaps)
{
    double *F0 = TMALLOC(double, n);
    int err;
    if (!F0) return 1;
    err = env_run_period_at(ckt, map, y, t0, T, "ENV Newton (Phi map)", 0, 0);   /* scratch: no output (card #719: was kept, a duplicate period in the plot) */
    (*nmaps)++;
    if (!err) {
        env_state_get(map, ckt, F0);
        err = env_phi_fd(ckt, map, n, y, F0, t0, T, floor_v, jac, nmaps);
    }
    tfree(F0);
    return err;
}

/* Earliest corner strictly after t_from (with a small guard) of any
 * native PWL V/I source: those are the KNOWN slow events of a deck
 * (soft-start knees, load/line steps).  An envelope jump never straddles
 * one -- its landing is capped at that corner and the corner itself is
 * crossed by ordinary real periods -- because an event inside a skipped
 * span is invisible to every estimate made from the two sides (see
 * docs/envelope.md, "Refusal behaviour").  Repeating PWL (r=) is treated
 * as non-repeating.  Returns 1e300 if there is none. */
static double
env_next_pwl_corner(CKTcircuit *ckt, double t_from)
{
    double best = 1e300;
    int vt = CKTtypelook("Vsource");
    int it = CKTtypelook("Isource");
    int i;
    if (vt >= 0) {
        VSRCmodel *m; VSRCinstance *h;
        for (m = (VSRCmodel *) ckt->CKThead[vt]; m; m = VSRCnextModel(m))
            for (h = VSRCinstances(m); h; h = VSRCnextInstance(h))
                if (h->VSRCfunctionType == PWL && h->VSRCcoeffs)
                    for (i = 0; i < h->VSRCfunctionOrder; i += 2) {
                        double tc = h->VSRCcoeffs[i] + h->VSRCrdelay;
                        if (tc > t_from && tc < best) best = tc;
                    }
    }
    if (it >= 0) {
        ISRCmodel *m; ISRCinstance *h;
        for (m = (ISRCmodel *) ckt->CKThead[it]; m; m = ISRCnextModel(m))
            for (h = ISRCinstances(m); h; h = ISRCnextInstance(h))
                if (h->ISRCfunctionType == PWL && h->ISRCcoeffs)
                    for (i = 0; i < h->ISRCfunctionOrder; i += 2) {
                        double tc = h->ISRCcoeffs[i] + h->ISRCrdelay;
                        if (tc > t_from && tc < best) best = tc;
                    }
    }
    return best;
}

/* ---------------------------------------------------------------- */
/* Top-level driver: bootstrap two real periods, then per period either
 * an implicit envelope jump of M periods (predictor-limited M, Newton,
 * post-check) followed by real landing/confirmation periods, or one
 * plain real period (M=1).  See the header comment. */
int
ENVenvelopeRun(CKTcircuit *ckt, PSSan *job)
{
    EnvStateMap map;
    EnvJac jac;
    int n, err, k;
    double T = 1.0 / job->PSSguessedFreq;
    double Tstop = job->EnvTstop;
    int Mmax = job->EnvMmax > 1 ? job->EnvMmax : 1;
    double envtol = job->EnvTol > 0.0 ? job->EnvTol : 0.0005;
    double tol = job->steady_coeff > 1e-7 ? job->steady_coeff : 1e-4;
    int maxiter = job->sc_iter > 0 ? (job->sc_iter < 8 ? job->sc_iter : 8) : 8;
    double growth_tol = getenv("ENV_GROWTH_TOL") ? atof(getenv("ENV_GROWTH_TOL")) : 5.0;
    int n_check = getenv("ENV_NCHECK") ? atoi(getenv("ENV_NCHECK")) : 1;   /* card #719: was 2 */
    double phase = job->EnvPhase > 0.0 && job->EnvPhase < 1.0 ? job->EnvPhase : 0.0;
    double t_wall_start, t_wall_end;
    double *x_pp = NULL, *x_prev = NULL, *x_n = NULL, *x_next = NULL,
           *y = NULL, *y_pred = NULL, *cvec = NULL, *cfil = NULL;
    double t_n = 0.0;
    int consec = 0;              /* consecutive real states ending at x_n */
    int after_jump = 0;          /* card #719: the last step was an accepted jump */
    double M_last = 1.0, *x_j0 = NULL;
    EnvHist hist;
    runDesc *samples_plot = NULL;
    long periods_simulated = 0, n_jumps = 0, n_refused_growth = 0,
         n_refused_newton = 0, n_refused_error = 0, n_mid_warn = 0,
         n_maps = 0, n_skipped = 0;
    long total_periods_target;
    double M_try;
    int backoff_left = 0, backoff_len = 4, err_streak = 0;
    int debug = getenv("ENV_DEBUG") != NULL;

    double t_scan = getenv("ENV_SCAN") ? atof(getenv("ENV_SCAN")) : -1.0;
    /* card #719: period maps at a tight transient tolerance (see the header
     * comment, "Map tolerance"); ENV_MAPTOL_OFF=1 keeps the deck's reltol
     * (the card #513 behaviour, red proof). */
    double save_reltol = ckt->CKTreltol;
    double maptol = job->EnvMapTolGiven ? job->EnvMapTol : 1e-5;

    g_env_rejected = 0;
    if (n_check < 1) n_check = 1;
    memset(&hist, 0, sizeof hist);
    memset(&jac, 0, sizeof jac);
    memset(&map, 0, sizeof map);

    if (job->PSSguessedFreq <= 0.0 || Tstop <= 0.0) {
        fprintf(stderr, "ENV: .tranenv needs fguess>0 and Tstop>0.\n");
        return E_UNSUPP;
    }
    err = env_build_state_map(ckt, &map);
    if (err) return err;
    n = map.n;
    if (n == 0) {
        fprintf(stderr, "ENV: no native capacitor/inductor found -- .tranenv "
                "needs at least one (same reduced state as .pss method=shoot).\n");
        env_map_free(&map);
        return E_UNSUPP;
    }
    x_pp = TMALLOC(double, n); x_prev = TMALLOC(double, n);
    x_n = TMALLOC(double, n); x_next = TMALLOC(double, n);
    y = TMALLOC(double, n); y_pred = TMALLOC(double, n);
    jac.Phi = TMALLOC(double, (size_t) n * n);
    cvec = TMALLOC(double, n); cfil = TMALLOC(double, n);
    x_j0 = TMALLOC(double, n);
    if (!cvec || !cfil || !x_j0) { err = E_NOMEM; goto env_cleanup; }
    if (!x_pp || !x_prev || !x_n || !x_next || !y || !y_pred || !jac.Phi) {
        fprintf(stderr, "ENV: out of memory allocating state vectors.\n");
        err = E_NOMEM;
        goto env_cleanup;
    }

    if (maptol > 0.0 && maptol < ckt->CKTreltol && !getenv("ENV_MAPTOL_OFF"))
        ckt->CKTreltol = maptol;
    t_wall_start = SPfrontEnd->IFseconds();
    total_periods_target = (long) (Tstop / T + 0.5);
    fprintf(stdout, "\nENV: implicit envelope following started -- T=%.6g s "
            "(fs=%.6g Hz), Tstop=%.6g s (~%ld periods), n_states=%d, Mmax=%d, "
            "envtol=%.4g, newton_tol=%.4g, maxiter=%d, n_check=%d, map reltol=%g\n",
            T, job->PSSguessedFreq, Tstop, total_periods_target, n,
            Mmax, envtol, tol, maxiter, n_check, ckt->CKTreltol);

    PSSenvKeptPlot = NULL;
    samples_plot = env_open_samples_plot(ckt, (JOB *) job, &map);
    if (!samples_plot)
        fprintf(stderr, "ENV: WARNING -- could not open the 'Envelope Samples' plot.\n");

    /* Bootstrap: [0, phase*T] (if phase>0) and two real periods from
     * t0 = phase*T (uic, the deck's IC=/.ic). */
    env_state_get_ic(&map, x_next);
    env_emit_sample(samples_plot, 0.0, x_next, n);
    if (phase > 0.0) {
        err = env_run_period_at(ckt, &map, x_next, 0.0, phase * T, "Envelope Following (card #513)", 0, 0);
        if (err) { fprintf(stderr, "ENV: phase run failed (error %d).\n", err); goto env_cleanup; }
        env_state_get(&map, ckt, x_next);
        env_emit_sample(samples_plot, phase * T, x_next, n);
    }
    err = env_run_period_at(ckt, &map, x_next, phase * T, T, "Envelope Following (card #513)", 0, 0);
    if (err) { fprintf(stderr, "ENV: first period failed (error %d).\n", err); goto env_cleanup; }
    env_state_get(&map, ckt, x_prev);
    env_emit_sample(samples_plot, phase * T + T, x_prev, n);
    periods_simulated++; n_maps++;
    err = env_run_period_at(ckt, &map, x_prev, phase * T + T, T, "Envelope Following (card #513)", 0, 0);
    if (err) { fprintf(stderr, "ENV: second period failed (error %d).\n", err); goto env_cleanup; }
    env_state_get(&map, ckt, x_n);
    t_n = phase * T + 2.0 * T;
    env_emit_sample(samples_plot, t_n, x_n, n);
    periods_simulated++; n_maps++;
    consec = 2;
    env_hist_push(&hist, env_maxrel(x_n, x_prev, n, 1e-3));

    M_try = 2.0;   /* first jump attempts start small and double */

/* advance one real period: (x_pp,x_prev,x_n) <- (x_prev,x_n,x_next) */
#define ENV_ADVANCE(tnew)                                                    \
    do {                                                                     \
        for (k = 0; k < n; k++) { x_pp[k] = x_prev[k]; x_prev[k] = x_n[k];   \
                                  x_n[k] = x_next[k]; }                      \
        t_n = (tnew); consec++;                                              \
        env_emit_sample(samples_plot, t_n, x_n, n);                          \
    } while (0)

    while (t_n < Tstop - 1e-6 * T && (Tstop - t_n) >= T * (1.0 - 1e-6)) {
        long remaining = (long) ((Tstop - t_n) / T + 1e-6);
        double M = 1.0, floor_v = 1e-3, maxabs = 0.0;
        int jumped = 0;
        double grow = 2.0;

        for (k = 0; k < n; k++) if (fabs(x_n[k]) > maxabs) maxabs = fabs(x_n[k]);
        if (0.01 * maxabs > floor_v) floor_v = 0.01 * maxabs;
        if (t_scan >= 0.0 && t_n >= t_scan) {
            env_scan(ckt, &map, n, x_n, t_n, T, floor_v);
            t_scan = -1.0;
        }

        /* choose M: doubling ladder, capped by Mmax, the remaining time
         * (the landing period must fit) and the curvature predictor:
         * the explicit projection over M periods has error ~ M^2*|c|/2
         * with c the last per-period change of the per-period slope. */
        if (backoff_left > 0) backoff_left--;   /* after refusals: plain real periods only */
        if (Mmax > 1 && backoff_left == 0 && remaining >= 3 &&
            (consec >= 3 || (consec == 2 && after_jump))) {
            double Mc = M_try, Mpred = 1e30;
            if (!jac.valid) {
                /* seed the monodromy cache once, at the current state */
                if (env_phi_at(ckt, &map, n, x_n, t_n, T, floor_v, &jac, &n_maps) != 0)
                    Mc = 1.0;
            }
            if (Mc > Mmax) Mc = Mmax;
            if (Mc > (double) (remaining - 2)) Mc = (double) (remaining - 2);
            if (jac.valid) {
                if (consec >= 3)
                    for (k = 0; k < n; k++)
                        cvec[k] = (x_n[k] - x_prev[k]) - (x_prev[k] - x_pp[k]);
                else
                    /* card #719: right after a landing (y at t_j0 + M_last*T,
                     * F(y) one period later) the per-period slope change is
                     * taken between the jump's secant slope and the landing
                     * period's slope, (M_last+1)/2 periods apart -- no extra
                     * real period is needed before the next jump, and the
                     * estimate averages the map noise over the jump. */
                    for (k = 0; k < n; k++)
                        cvec[k] = ((x_n[k] - x_prev[k]) - (x_prev[k] - x_j0[k]) / M_last)
                                  / (0.5 * (M_last + 1.0));
                env_slow_mask(jac.Phi, n, cvec, cfil);
                int klim = -1;
                for (k = 0; k < n; k++) {
                    double c = fabs(cfil[k]);
                    double mk = c > 1e-300 ? sqrt(2.0 * envtol * 0.3 * env_scale(x_n[k], floor_v) / c) : 1e30;
                    if (mk < Mpred) { Mpred = mk; klim = k; }
                }
                M = floor(Mc < Mpred ? Mc : Mpred);
                if (debug && getenv("ENV_DEBUG")[0] == '4')
                    fprintf(stderr, "ENV: t=%.6g Mc=%.4g Mpred=%.4g (%s)\n", t_n, Mc, Mpred,
                            klim >= 0 ? map.refs[klim].name : "-");
            }
            if (M < 1.0) M = 1.0;
        }
        if (M > 1.0) {   /* never land past a known slow-source corner */
            double tc = env_next_pwl_corner(ckt, t_n + 0.5 * T);
            if (tc < 1e299) {
                double Mcorner = floor((tc - t_n) / T) - 2.0;   /* land >= 2 periods before it */
                if (Mcorner < M) {
                    if (debug) fprintf(stderr, "ENV: t=%.6g corner at %.6g caps M %.4g -> %.4g\n", t_n, tc, M, Mcorner);
                    M = Mcorner < 1.0 ? 1.0 : Mcorner;
                }
            }
        }
        if (M > 1.0 && hist.n >= 2) {
            double delta_now = env_maxrel(x_n, x_prev, n, 1e-3);
            double med = env_hist_median(&hist);
            if (med > 0.0 && delta_now > growth_tol * med) {
                M = 1.0;
                n_refused_growth++;
                if (debug)
                    fprintf(stderr, "ENV: refuse jump at t=%.6g (growth) delta=%.4g median=%.4g\n",
                            t_n, delta_now, med);
            }
        }

        if (M >= 2.0) {
            for (k = 0; k < n; k++)
                y[k] = y_pred[k] = x_n[k] + M * (x_n[k] - x_prev[k]);
            if (env_newton_step(ckt, &map, n, x_n, M, t_n + M * T, T, tol, maxiter,
                                floor_v, &jac, y, &n_maps) != 0) {
                n_refused_newton++;
                M_try = M / 2.0 < 2.0 ? 2.0 : M / 2.0;
                /* exponential back-off: a deck whose one-period map is not
                 * smooth enough for Newton must not pay for a failed
                 * solve every period (docs/envelope.md) */
                backoff_left = backoff_len;
                if (backoff_len < 128) backoff_len *= 2;
                if (debug) fprintf(stderr, "ENV: Newton failed M=%.6g t=%.6g -> back-off %d periods\n", M, t_n, backoff_left);
            } else {
                double e;
                for (k = 0; k < n; k++) cvec[k] = y[k] - y_pred[k];
                env_slow_mask(jac.Phi, n, cvec, cfil);
                e = 0.0;
                for (k = 0; k < n; k++) {
                    double ek = fabs(cfil[k]) / env_scale(y[k], floor_v);
                    if (ek > e) e = ek;
                }
                if (e > envtol) {
                    double f = sqrt(envtol / e) * 0.9;
                    n_refused_error++;
                    M_try = floor(M * f);
                    if (M_try < 2.0) M_try = 2.0;
                    if (++err_streak >= 3) {   /* repeated error refusals: back off too */
                        backoff_left = backoff_len;
                        if (backoff_len < 128) backoff_len *= 2;
                        err_streak = 0;
                    }
                    if (debug) fprintf(stderr, "ENV: err_est=%.4g > envtol=%.4g M=%.6g t=%.6g -> M_try=%.6g\n",
                                       e, envtol, M, t_n, M_try);
                } else {
                    jumped = 1;
                    /* card #719: next M from this step's error estimate (it
                     * grows like M^2), at most doubled -- a plain doubling
                     * after a step that used 60 % of envtol was refused next
                     * time and alternated accept/refuse on the PCM deck */
                    grow = e > 0.0 ? 0.9 * sqrt(envtol / e) : 2.0;
                    if (grow > 2.0) grow = 2.0;
                    if (grow < 1.0) grow = 1.0;
                    if (debug) fprintf(stderr, "ENV: accept M=%.6g t=%.6g err_est=%.4g\n", M, t_n, e);
                }
            }
        }

        after_jump = 0;
        if (jumped) {
            double t_land = t_n + M * T;
            int c;
            for (k = 0; k < n; k++) x_j0[k] = x_n[k];   /* card #719: jump start, for the curvature after the landing */
            M_last = M;
            err = env_run_period_at(ckt, &map, y, t_land, T, "Envelope Following (card #513)", 0, 0);
            n_maps++;
            if (err) {
                fprintf(stderr, "ENV: landing run failed (error %d) at t=%.6g; "
                        "falling back to a plain period.\n", err, t_land);
                n_refused_newton++;
                jumped = 0;
            } else {
                env_state_get(&map, ckt, x_next);
                n_jumps++; n_skipped += (long) M - 1;
                periods_simulated++;
                /* history restarts: states (y at t_land, F(y) at t_land+T) */
                for (k = 0; k < n; k++) x_n[k] = y[k];
                env_emit_sample(samples_plot, t_land, x_n, n);
                consec = 1;
                env_hist_push(&hist, env_maxrel(x_next, y, n, 1e-3));
                ENV_ADVANCE(t_land + T);
                for (c = 1; c < n_check && (Tstop - t_n) >= T * (1.0 - 1e-6); c++) {
                    double d_new, med;
                    err = env_run_period_at(ckt, &map, x_n, t_n, T, "Envelope Following (card #513)", 0, 0);
                    n_maps++;
                    if (err) { fprintf(stderr, "ENV: confirmation period failed (error %d) at t=%.6g.\n", err, t_n); break; }
                    env_state_get(&map, ckt, x_next);
                    periods_simulated++;
                    d_new = env_maxrel(x_next, x_n, n, 1e-3);
                    med = env_hist_median(&hist);
                    if (med > 0.0 && d_new > growth_tol * med) {
                        n_mid_warn++;
                        fprintf(stderr, "ENV: WARNING -- possible event after the landing at "
                                "t=%.6g (delta=%.4g, median=%.4g); the jump is not undone.\n",
                                t_n, d_new, med);
                    }
                    env_hist_push(&hist, d_new);
                    ENV_ADVANCE(t_n + T);
                }
                M_try = floor(grow * M);
                if (M_try < 2.0) M_try = 2.0;
                if (M_try > Mmax) M_try = Mmax;
                backoff_len = 4; err_streak = 0;
                after_jump = 1;
            }
        }

        if (!jumped) {
            err = env_run_period_at(ckt, &map, x_n, t_n, T, "Envelope Following (card #513)", 0, 0);
            n_maps++;
            if (err) {
                fprintf(stderr, "ENV: period failed (error %d) at t=%.6g -- stopping early.\n", err, t_n);
                goto env_cleanup;
            }
            env_state_get(&map, ckt, x_next);
            periods_simulated++;
            env_hist_push(&hist, env_maxrel(x_next, x_n, n, 1e-3));
            ENV_ADVANCE(t_n + T);
        }
    }

    if (Tstop - t_n > 1e-6 * T) {   /* final partial period (only when phase > 0 or Tstop is not a multiple of T) */
        err = env_run_period_at(ckt, &map, x_n, t_n, Tstop - t_n, "Envelope Following (card #513)", 0, 0);
        if (err) fprintf(stderr, "ENV: final partial period failed (error %d).\n", err);
        n_maps++;
    }
    t_wall_end = SPfrontEnd->IFseconds();
    fprintf(stdout, "ENV: DONE periods_simulated=%ld/%ld (%.1f%%) skipped=%ld jumps=%ld "
            "period_maps=%ld refused(newton)=%ld refused(error)=%ld refused(growth)=%ld "
            "post_landing_warnings=%ld rejected_steps=%ld map_reltol=%g wall_time=%.3fs\n",
            periods_simulated, total_periods_target,
            100.0 * (double) periods_simulated / (double) total_periods_target,
            n_skipped, n_jumps, n_maps, n_refused_newton, n_refused_error,
            n_refused_growth, n_mid_warn, g_env_rejected, ckt->CKTreltol, t_wall_end - t_wall_start);
    err = OK;

env_cleanup:
    ckt->CKTreltol = save_reltol;
    if (PSSenvKeptPlot) { SPfrontEnd->OUTendPlot((runDesc *) PSSenvKeptPlot); PSSenvKeptPlot = NULL; }
    if (samples_plot) SPfrontEnd->OUTendPlot(samples_plot);
    if (x_pp) tfree(x_pp);
    if (x_prev) tfree(x_prev);
    if (x_n) tfree(x_n);
    if (x_next) tfree(x_next);
    if (y) tfree(y);
    if (y_pred) tfree(y_pred);
    if (jac.Phi) tfree(jac.Phi);
    if (cvec) tfree(cvec);
    if (cfil) tfree(cfil);
    if (x_j0) tfree(x_j0);
    env_map_free(&map);
    return err;
}
