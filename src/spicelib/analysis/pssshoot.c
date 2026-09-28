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


/* Spectral radius of Phi (n*n row-major) by power iteration, plus a
 * cheap period-doubling flag: if the dominant eigenvalue is real and
 * negative, the power-iteration vector's sign flips every pass, which
 * shows up as a negative dot product between successive (normalised)
 * iterates once the direction has settled down. */
static double
pss_spectral_radius(double *Phi, int n, int *flip_flag)
{
    double *v, *w;
    int i, j, iter;
    double lambda = 0.0, prevlambda = -1.0;

    *flip_flag = 0;
    if (n == 1) {
        *flip_flag = (Phi[0] < 0.0);
        return fabs(Phi[0]);
    }

    v = TMALLOC(double, n);
    w = TMALLOC(double, n);
    if (!v || !w) {
        if (v) tfree(v);
        if (w) tfree(w);
        return 0.0;
    }
    for (i = 0; i < n; i++)
        v[i] = 1.0 / sqrt((double) n);

    for (iter = 0; iter < 300; iter++) {
        double norm = 0.0, dot;
        for (i = 0; i < n; i++) {
            double s = 0.0;
            for (j = 0; j < n; j++)
                s += Phi[i * n + j] * v[j];
            w[i] = s;
        }
        for (i = 0; i < n; i++)
            norm += w[i] * w[i];
        norm = sqrt(norm);
        if (norm < 1e-300) {
            lambda = 0.0;
            break;
        }
        dot = 0.0;
        for (i = 0; i < n; i++)
            dot += v[i] * (w[i] / norm);
        for (i = 0; i < n; i++)
            v[i] = w[i] / norm;
        lambda = norm;
        *flip_flag = (dot < 0.0);
        if (iter > 5 && fabs(lambda - prevlambda) < 1e-9 * (lambda > 1e-300 ? lambda : 1.0))
            break;
        prevlambda = lambda;
    }
    tfree(v);
    tfree(w);
    return lambda;
}


/* ------------------------------------------------------------------ */

int
PSSshootDriven(CKTcircuit *ckt, PSSan *job)
{
    double T = 1.0 / ckt->CKTguessedFreq;
    PSSstateMap map = { 0, NULL };
    int n, k, j, iter, err, maxiter, flip_flag = 0;
    double tol, eps, residual = 0.0, specrad = 0.0;
    double t_wall_start, t_wall_end;
    int converged_flag = 0;
    double *x0 = NULL, *x1 = NULL, *F = NULL, *Phi = NULL, *delta_x = NULL,
           *xtrial = NULL, *x1trial = NULL, *scale = NULL, *A = NULL, *b = NULL;

    t_wall_start = SPfrontEnd->IFseconds();

    /* Integration resolution within a period: match the sampling the
     * deck asked for (Points), or a fixed default. */
    ckt->CKTmaxStep  = (job->PSSpoints > 0) ? T / (double) job->PSSpoints : T / 200.0;
    ckt->CKTstep     = ckt->CKTmaxStep;
    ckt->CKTdelmin   = 1e-9 * ckt->CKTmaxStep;
    ckt->CKTminBreak = 0.0;
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

    if (!x0 || !x1 || !F || !Phi || !delta_x || !xtrial || !x1trial || !scale || !A || !b) {
        err = E_NOMEM;
        goto shoot_cleanup;
    }

    /* Newton's starting guess. The deck's raw IC=/.ic values (usually all
     * zero on the reactive states) are very often themselves an exact,
     * degenerate fixed point of the one-period map -- a converter sitting
     * powered-off with its compensator/comparator at whatever CKTop finds
     * from a cold start is self-consistent and non-switching, so P(0)=0
     * and Newton "converges" to that in a couple of iterations without
     * ever seeing the real orbit (found running the first R3 smoke test
     * this way: it reported CONVERGED with v(vout)~1e-7, not ~3.3V). A
     * plain `.tran` from the same cold start does NOT get stuck there,
     * because its comparator/ramp/op-amp dynamics carry continuous history
     * across many periods; a single one-period Newton evaluation has none
     * of that. So: reuse the existing `stabtime` PSS parameter (already
     * parsed from the .pss card, previously only meaningful to the old
     * algorithm's own STABILIZATION phase) to run one ordinary transient
     * first and seed x0 from wherever THAT lands -- close enough to the
     * true orbit for Newton to find it instead of the degenerate one. */
    if (ckt->CKTstabTime > 0.0) {
        fprintf(stdout, "\nPSS-SHOOT: pre-stabilizing with a %.6g s "
                "transient (stabtime, ~%.0f periods) before the first "
                "Newton guess...\n", ckt->CKTstabTime, ckt->CKTstabTime / T);
        fflush(stdout);
        err = pss_run_period(ckt, ckt->CKTstabTime,
                              "PSS shoot (pre-stabilization transient)");
        if (err) {
            fprintf(stderr, "PSS-SHOOT: FAILED -- the pre-stabilization "
                    "transient (stabtime=%g) did not complete (ngspice "
                    "error %d).\n", ckt->CKTstabTime, err);
            goto shoot_cleanup;
        }
        pss_state_get(&map, ckt, x0);
    } else {
        pss_state_get_ic(&map, x0);
    }

    maxiter = (job->sc_iter > 0) ? job->sc_iter : 20;
    tol     = (job->steady_coeff > 0.0) ? job->steady_coeff : 1e-6;
    /* Finite-difference perturbation size for the monodromy columns. Must
     * sit above the numerical noise floor of the transient itself (each
     * internal Newton solve only converges to CKTreltol, default 1e-3) or
     * the finite difference is dominated by noise, not signal -- found on
     * the R3 fixture: eps=1e-6 gave a Newton direction that could not
     * reduce the residual even at very small step sizes (the damped
     * backtracking exhausted all 8 halvings), consistent with a noisy
     * monodromy estimate rather than a genuinely bad direction. */
    eps     = 1e-2;

    fprintf(stdout, "\nPSS-SHOOT: shooting-Newton PSS, period T=%g s "
            "(f=%g Hz), %d state(s):", T, ckt->CKTguessedFreq, n);
    for (k = 0; k < n; k++)
        fprintf(stdout, " %s", map.refs[k].name);
    fprintf(stdout, "; max %d Newton iterations, tol %g\n", maxiter, tol);
    fflush(stdout);

    for (iter = 0; ; iter++) {

        pss_state_set_ic(&map, x0);
        err = pss_run_period(ckt, T, "PSS shoot (Newton baseline period)");
        if (err) {
            fprintf(stderr, "PSS-SHOOT: FAILED -- period integration did "
                    "not complete at Newton iteration %d (ngspice error "
                    "%d); this is NOT a steady state.\n", iter, err);
            goto shoot_fail;
        }
        pss_state_get(&map, ckt, x1);

        /* Scale: max(|x0|,|x1|), floored at ABSTOL_FLOOR -- not the
         * near-zero floor (1e-12) tried first. Found on the R3 fixture:
         * a compensator summing-node voltage (C3in) sits at ~0.9 mV in the
         * true orbit, so a purely relative scale amplifies its few-tens-
         * of-microvolt error into a ~0.057 "residual", swamping every
         * other (properly O(volts)) state and making the max-norm
         * residual both numerically noisy (small absolute noise -> huge
         * relative swings near zero) and a bad line-search metric --
         * backtracking could not find ANY step, at any alpha down to
         * 1e-7, that reduced this dominated-by-noise max. A 1 mV floor
         * (chosen against this circuit's ripple scales -- vout's ripple
         * is ~5 mV per refs/README.md, so this floor is smaller than a
         * physically meaningful ripple and only suppresses the near-zero
         * node's blow-up) fixed it outright. */
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
        fprintf(stdout, "PSS-SHOOT: iter=%d residual=%.6e\n", iter, residual);
        fflush(stdout);
        if (getenv("PSS_SHOOT_DEBUG")) {
            for (k = 0; k < n; k++)
                fprintf(stderr, "DEBUG iter=%d k=%d %s x0=%.8g x1=%.8g "
                        "F=%.4e scale=%.6g r=%.4e\n",
                        iter, k, map.refs[k].name, x0[k], x1[k], F[k],
                        scale[k], fabs(F[k]) / scale[k]);
        }

        if (residual <= tol) {
            converged_flag = 1;
            break;
        }
        if (iter >= maxiter)
            break;   /* out of iterations -- reported as FAILED below */

        /* finite-difference monodromy: 2n extra period runs (centred
         * differences, not the cheaper n forward differences the card's
         * "n+1 period runs per Newton step" description allows -- found on
         * the R3 fixture that one-sided differences gave a Newton
         * direction whose backtrack line search approached the baseline
         * residual from above for 24 straight halvings without ever
         * dipping below it, i.e. essentially no usable descent direction;
         * centred differences cancel the O(h) bias one-sided differences
         * carry and fixed it. Costs 2n period runs instead of n+1 per
         * Newton step -- still cheap at this n. */
        for (j = 0; j < n; j++) {
            double h = fabs(x0[j]) * eps;
            if (h < 1e-6)
                h = 1e-6;

            for (k = 0; k < n; k++)
                xtrial[k] = x0[k];
            xtrial[j] = x0[j] + h;
            pss_state_set_ic(&map, xtrial);
            err = pss_run_period(ckt, T, "PSS shoot (monodromy perturbation +)");
            if (err) {
                fprintf(stderr, "PSS-SHOOT: FAILED -- monodromy "
                        "perturbation run for state %s (+h) did not "
                        "complete at Newton iteration %d (ngspice error "
                        "%d).\n", map.refs[j].name, iter, err);
                goto shoot_fail;
            }
            pss_state_get(&map, ckt, x1trial);
            for (k = 0; k < n; k++)
                Phi[k * n + j] = x1trial[k];   /* stash P(x0+h e_j), col j */

            for (k = 0; k < n; k++)
                xtrial[k] = x0[k];
            xtrial[j] = x0[j] - h;
            pss_state_set_ic(&map, xtrial);
            err = pss_run_period(ckt, T, "PSS shoot (monodromy perturbation -)");
            if (err) {
                fprintf(stderr, "PSS-SHOOT: FAILED -- monodromy "
                        "perturbation run for state %s (-h) did not "
                        "complete at Newton iteration %d (ngspice error "
                        "%d).\n", map.refs[j].name, iter, err);
                goto shoot_fail;
            }
            pss_state_get(&map, ckt, x1);   /* reuse x1 as scratch for P(x0-h e_j) */
            for (k = 0; k < n; k++)
                Phi[k * n + j] = (Phi[k * n + j] - x1[k]) / (2.0 * h);
        }

        /* x1 was overwritten as perturbation scratch above -- recompute the
         * true baseline P(x0) once more for the F(x0) used below (cheap:
         * one extra period run, traded for centred-difference accuracy). */
        pss_state_set_ic(&map, x0);
        err = pss_run_period(ckt, T, "PSS shoot (Newton baseline, post-monodromy)");
        if (err) {
            fprintf(stderr, "PSS-SHOOT: FAILED -- baseline period "
                    "integration did not complete after the monodromy "
                    "sweep at Newton iteration %d (ngspice error %d).\n",
                    iter, err);
            goto shoot_fail;
        }
        pss_state_get(&map, ckt, x1);

        for (k = 0; k < n * n; k++) A[k] = Phi[k];
        for (k = 0; k < n; k++) A[k * n + k] -= 1.0;
        for (k = 0; k < n; k++) b[k] = -F[k];

        if (pss_solve_dense(A, b, n, delta_x)) {
            fprintf(stderr, "PSS-SHOOT: FAILED -- the monodromy matrix "
                    "(Phi - I) is singular at Newton iteration %d; Newton "
                    "cannot proceed.\n", iter);
            goto shoot_fail;
        }
        if (getenv("PSS_SHOOT_DEBUG")) {
            fprintf(stderr, "DEBUG iter=%d x0:", iter);
            for (k = 0; k < n; k++) fprintf(stderr, " %.6g", x0[k]);
            fprintf(stderr, "\nDEBUG iter=%d x1:", iter);
            for (k = 0; k < n; k++) fprintf(stderr, " %.6g", x1[k]);
            fprintf(stderr, "\nDEBUG iter=%d F:", iter);
            for (k = 0; k < n; k++) fprintf(stderr, " %.6g", F[k]);
            fprintf(stderr, "\nDEBUG iter=%d delta_x:", iter);
            for (k = 0; k < n; k++) fprintf(stderr, " %.6g", delta_x[k]);
            fprintf(stderr, "\nDEBUG iter=%d Phi diag:", iter);
            for (k = 0; k < n; k++) fprintf(stderr, " %.6g", Phi[k*n+k]);
            fprintf(stderr, "\n");
        }

        {
            double alpha = 1.0;
            int bt, accepted = 0;
            /* 24, not 8: found on the R3 fixture that a Newton direction
             * dominated by a fast/stiff coordinate (Phi diagonal ~1e-7 for
             * the switch-node capacitor, i.e. it forgets its own initial
             * value almost completely every period, but couples strongly
             * to the other columns) can need quite a few halvings before
             * the trial residual actually dips below the baseline -- the
             * backtrack sequence was still monotonically approaching (not
             * yet below) the baseline residual at the old limit of 8. */
            for (bt = 0; bt < 24; bt++) {
                for (k = 0; k < n; k++)
                    xtrial[k] = x0[k] + alpha * delta_x[k];

                pss_state_set_ic(&map, xtrial);
                err = pss_run_period(ckt, T, "PSS shoot (damped-step trial)");
                if (!err) {
                    double newresidual = 0.0;
                    pss_state_get(&map, ckt, x1trial);
                    /* Reuse the SAME scale[] the baseline residual above was
                     * measured with (not a fresh max(|xtrial|,|x1trial|) per
                     * trial) -- a per-trial rescale let the metric drift
                     * along with the trial point itself and was part of why
                     * the line search on the R3 fixture could not find a
                     * usable step (compounding the near-zero-state scaling
                     * bug fixed above). A fixed reference scale for the
                     * whole Newton iteration is the standard, well-behaved
                     * choice. */
                    for (k = 0; k < n; k++) {
                        double r = fabs(x1trial[k] - xtrial[k]) / scale[k];
                        if (r > newresidual) newresidual = r;
                    }
                    if (getenv("PSS_SHOOT_DEBUG"))
                        fprintf(stderr, "DEBUG iter=%d bt=%d alpha=%g "
                                "newresidual=%.6e\n", iter, bt, alpha, newresidual);
                    /* Strict decrease only -- an earlier version also
                     * force-accepted once alpha got small regardless of
                     * whether the step helped, which let Newton "wander"
                     * at a roughly constant residual for many iterations
                     * instead of either converging or failing cleanly
                     * (found on the R3 fixture: residual oscillated
                     * around 0.02-0.06 for 25 iterations without ever
                     * reaching tol). Failing here is the honest outcome
                     * when even a tiny step does not help. */
                    if (newresidual < residual) {
                        accepted = 1;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if (!accepted) {
                fprintf(stderr, "PSS-SHOOT: FAILED -- the damped Newton "
                        "step did not reduce the orbit residual at "
                        "iteration %d; no plain-transient fallback is "
                        "implemented (card #458, what was not done).\n",
                        iter);
                goto shoot_fail;
            }
            for (k = 0; k < n; k++)
                x0[k] = xtrial[k];
        }
    }

    if (!converged_flag) {
        fprintf(stderr, "PSS-SHOOT: FAILED to converge within %d Newton "
                "iterations (final residual %.6e, tol %.3e); this is NOT "
                "a steady state.\n", maxiter, residual, tol);
        goto shoot_fail;
    }

    /* Floquet report: recompute Phi once more at x0 (now converged) with a
     * smaller, steadier perturbation than the eps=1e-2 Newton ended up
     * needing. eps=1e-2 was needed during the iteration to get a
     * perturbation clearly above the adaptive-step transient's own
     * noise floor (found on the R3 fixture: eps=1e-4 gave a Newton
     * direction that could not improve the residual at all, plateauing
     * well short of tol -- see the eps history above); that size is more
     * than the local linearisation really wants for an eigenvalue
     * estimate on these steep tanh switches, since the perturbed run can
     * cross into a visibly different (if still locally linear-ish)
     * regime. At the converged x0 (not mid-iteration) a smaller,
     * more classically-scaled eps is worth trying for this one purpose;
     * fall back to the Newton-loop Phi if this second pass fails to run. */
    {
        double eps_floquet = 1e-3;
        double *Phi2 = TMALLOC(double, n * n);
        int floquet_ok = (Phi2 != NULL);

        pss_state_set_ic(&map, x0);
        if (floquet_ok && pss_run_period(ckt, T, "PSS shoot (Floquet baseline)"))
            floquet_ok = 0;
        if (floquet_ok)
            pss_state_get(&map, ckt, x1);

        for (j = 0; floquet_ok && j < n; j++) {
            double h = fabs(x0[j]) * eps_floquet;
            if (h < 1e-6) h = 1e-6;

            for (k = 0; k < n; k++) xtrial[k] = x0[k];
            xtrial[j] = x0[j] + h;
            pss_state_set_ic(&map, xtrial);
            if (pss_run_period(ckt, T, "PSS shoot (Floquet perturbation +)")) {
                floquet_ok = 0;
                break;
            }
            pss_state_get(&map, ckt, x1trial);
            for (k = 0; k < n; k++)
                Phi2[k * n + j] = x1trial[k];

            for (k = 0; k < n; k++) xtrial[k] = x0[k];
            xtrial[j] = x0[j] - h;
            pss_state_set_ic(&map, xtrial);
            if (pss_run_period(ckt, T, "PSS shoot (Floquet perturbation -)")) {
                floquet_ok = 0;
                break;
            }
            pss_state_get(&map, ckt, x1trial);
            for (k = 0; k < n; k++)
                Phi2[k * n + j] = (Phi2[k * n + j] - x1trial[k]) / (2.0 * h);
        }

        if (floquet_ok) {
            specrad = pss_spectral_radius(Phi2, n, &flip_flag);
        } else {
            fprintf(stderr, "PSS-SHOOT: WARNING -- the dedicated Floquet "
                    "re-evaluation (eps=%g) failed to run; reporting the "
                    "spectral radius from the last Newton-step monodromy "
                    "(eps=%g) instead, which may be less accurate.\n",
                    eps_floquet, eps);
            specrad = pss_spectral_radius(Phi, n, &flip_flag);
        }
        if (Phi2) tfree(Phi2);
        /* restore the baseline transient's data (x1/Phi corrupted as
         * scratch above) is not needed further -- nothing below reads
         * them except the final trajectory re-run, which sets its own
         * state explicitly. */
    }
    t_wall_end = SPfrontEnd->IFseconds();

    fprintf(stdout, "\nPSS-SHOOT: CONVERGED iter=%d residual=%.6e "
            "spectral_radius=%.6f %s wall_time=%.3fs\n",
            iter, residual, specrad,
            specrad >= 1.0 ? (flip_flag ? "UNSTABLE(period-doubling)" : "UNSTABLE")
                           : "stable",
            t_wall_end - t_wall_start);

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
    if (map.refs) tfree(map.refs);
    return err;

 shoot_fail:
    err = E_ITERLIM;
    goto shoot_cleanup;
}
