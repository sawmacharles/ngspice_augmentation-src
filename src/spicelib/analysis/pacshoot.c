/**********
card #465 (plan T5.1-T5.2): time-domain shooting-based periodic AC (PAC,
sideband 0) and periodic loop gain (.pstb), for switching converters, built
on the shooting-Newton PSS (pssshoot.c) in this fork. Design rulings
#446 B (fork kept separate from upstream) and #447 A (smooth switch models
only, no saltation) apply, same as pssshoot.c.

READ FIRST: research/ngspice-internals-pss-pac.md sec 3.3 (why the existing
frequency-domain .pac -- dcpss.c's dense harmonic conversion matrix -- is
NOT extended here) and sec 5 (the LPTV linearisation and the time-domain
shooting-based PAC algorithm this file implements); research/simplis-pop-ac.md
sec 6-7 (the same algorithm in state-transition form); research/cadence-stb-
pstb.md sec 2.4 (the Tian two-injection probe formula) and sec 4 (pstb =
sideband-0 loop gain on the LPTV orbit).

===========================================================================
DESIGN
===========================================================================

Setting: PSSshootDriven() (pssshoot.c) has already converged a periodic
orbit x0 (the reduced state: one entry per native capacitor's terminal
voltage, one per native inductor's branch current -- SAME convention as
pssshoot.c's own PSSstateMap, independently rebuilt here by pac_build_map()
walking the SAME device lists in the SAME order) with known period T. Around
that orbit, C(t)=dq/dx and G(t)=df/dx are T-periodic (LPTV). A small input
b(t) = b*e^{jwt} (w = 2*pi*f, f the PAC/PSTB sweep frequency) produces a
response y(t) = ytilde(t)*e^{jwt} with ytilde T-periodic, equivalently the
Floquet boundary condition y(T) = e^{jwT} y(0).

Discretise on the PSS's OWN accepted-step grid with the SAME integration
rule the transient used (trapezoidal orders 1-2, backward Euler = order 1's
special case): at each step k,

    J_k y_k = (companion-history term built from y_{k-1}, same recursion
               capload.c/indload.c use for the real transient's ceq/veq)
              + b * e^{j*w*t_k}

with J_k = G_k + ag0_k*C_k EXACTLY the matrix the real transient's own
Newton iteration factored at that accepted step (ckt->CKTmatrix, reused via
SMPsolve -- no refactor). This is the SAME per-step recursion pssshoot.c's
own PSSsensRecordStep() (card #458 fourth pass, "analytic sensitivity")
already implements and validated on R2/R3, generalised three ways:
  (1) complex arithmetic (y_k, the companion-history terms, and J_k's
      inverse action are all propagated as (re,im) pairs -- J_k itself
      stays REAL, so each complex solve is two real SMPsolve() calls,
      exactly matching this card's design brief: "split real/imag parts
      since the matrix is real");
  (2) an added forcing term b*e^{j*w*t_k} (PSSsensRecordStep has none: it
      is the homogeneous sensitivity of the SAME state-transition map);
  (3) many "columns" (independent right-hand sides sharing ONE per-step
      factorization) propagated SIMULTANEOUSLY in ONE recording period
      run, not just the n homogeneous unit columns PSSsensRecordStep uses
      for the monodromy: this file ALSO runs, in the SAME pass, one
      forced complex column per requested sweep frequency (b=0 columns for
      the homogeneous monodromy Phi, b!=0 columns for the forced response
      z(w) -- see PACbatch/PSSpacRecordStep below). This is the "refactor
      once per accepted step of ONE final period, reuse for every
      frequency via back-substitution only" design the card's plan
      describes, EXCEPT that "refactor" here still means "reuse the
      transient's own last-Newton-iterate factorization" (pssshoot.c's
      existing, tested, but documented-stale J_k source -- see
      PSSsensRecordStep's header comment in pssshoot.c for the staleness
      finding on R3) rather than a fresh CKTload()+SMPluFac() at the
      exact converged point: that fresh-refactor option was tried by the
      PSS work (pssshoot.c's own comment) and found to be EITHER actively
      wrong (a hysteretic switch's history-dependent load) or too slow
      (R3, refactoring every accepted step) -- this file inherits that
      same, already-adjudicated trade-off rather than re-opening it.

Algorithm per forcing pattern (PAC: the netlist's own AC-flagged source,
captured once via the E-120/E-123 CKTacLoad(w=1) trick, SAME idiom
dcpss.c's pac_extract_harmonics() uses; PSTB: two Tian injections, voltage
then current, at the probe, SAME RHS convention stban.c's loop-mode uses):

  1. ONE recording period run (pac_run_period, unperturbed -- x0 is the
     REAL, converged IC every time; this bookkeeping is a passive parallel
     computation that never touches ckt's real state, so the transient's
     own adaptive-step decisions are IDENTICAL, bit for bit, on every one
     of these re-runs -- no breakpoint forcing needed, unlike the Newton/
     Floquet code's finite-difference perturbation runs) with columns =
     [ n homogeneous columns (b=0, w=0, y0 = unit vector e_d) for Phi ]
     + [ one forced column per sweep frequency f (b = forcing pattern,
         w = 2*pi*f, y0 = 0) for z(w) ].
     At t=T, Phi[d][j] and z_f[d] are read off the reduced-state rows.
  2. Dense complex solve per frequency (n small: 1 on a switched-RC, ~7 on
     R3/R4): (e^{jwT} I - Phi) xi = z(w).
  3. A SECOND recording period run, columns = one forced column per sweep
     frequency (b = same forcing pattern, w = 2*pi*f, y0 = xi_f from step
     2) -- by linearity this equals the ACTUAL periodic small-signal
     response, y0 correctly satisfying the Floquet boundary condition.
     This run has demodulation ACTIVE: at every accepted step, for every
     circuit row, it accumulates the sideband-0 trapezoidal quadrature
     Y0[row] += (1/T) * 0.5*h_k*(y_{k-1}[row]*e^{-jwt_{k-1}} +
                                  y_k[row]*e^{-jwt_k})
     (non-uniform-grid trapezoidal rule, exactly as the design brief asks
     for; PACbatch's t=0 seed value is also folded into the FIRST
     trapezoid segment).

PAC (`.pac ... method=shoot solver=td`): forcing = the netlist AC-source
RHS (B0). Output = sideband-0 Y0 at every circuit node, vs frequency (same
node-name convention as the existing frequency-domain .pac, via CKTnames).

PSTB (`.pstb ... probe=<Vname>`): TWO forcing patterns (Phi is shared,
computed once): RHS#1 = unit voltage injection in the probe branch
equation, RHS#2 = unit current injection at the probe's node_b -- BOTH
copied verbatim (same sign conventions) from stban.c's native `.stb`
loop-mode so the w->0 (DC) limit is BY CONSTRUCTION the same Tian formula
`.stb` already implements (this file duplicates stban.c's small Tian/
margin helpers rather than exposing them, since stban.c is another
worker's file this card must not edit). v1/i1 (from RHS#1's sideband-0
V(node_b)/I(branch)) and v2/i2 (from RHS#2's) feed Tian eq. 30 exactly as
stban.c does; margins are stban.c's own log-f-interpolation, phase-
unwrapping, worst-of-several-crossings algorithm, duplicated verbatim.

===========================================================================
WHAT THIS DOES NOT DO (see the card #465 report for the full list)
===========================================================================
  - No factorization CACHING across the frequency sweep beyond what a
    single recording period run gives "for free" for every frequency in
    that SAME run (steps 1 and 3 above): step 3 is one fresh period
    re-run, so cost scales as O(2 * one-period-transient) total for the
    WHOLE sweep for step 1+2, plus O(1 * one-period-transient) for step 3
    -- NOT the fully lower-triangular "store K factorizations, no more
    transient re-runs at all" design an infinite-memory implementation
    could reach; see the timing table in the report.
  - No matrix-free/Krylov path, no Arnoldi-basis recycling across
    frequency (research doc sec 5.3's efficiency point) -- the dense
    n x n boundary solve is used throughout (fine for the state counts
    here: 1-2 on R2, ~7 on R3/R4).
  - No sidebands other than 0 (PACmaxSideband is a frequency-domain-.pac-
    only field; this engine reports the same-frequency response only).
  - No pnoise on this engine. PXF (card #499, `.pxf ... method="shoot"
    solver="td"`, `pac_td_pxf()` below) IS implemented, as a batched
    forward per-source sweep reusing pac_engine_run() below unchanged --
    NOT the from-scratch backward-adjoint engine that name might suggest;
    see pac_td_pxf()'s own "SHIPPED DESIGN" comment for the alternative
    considered and why.
  - No device-based (gain-nulling) PSTB probe. PSTB's probe= is an
    existing 0V V-source, <inst>:<terminal> (auto-inserted at parse time
    by the same splice as .stb), or an instance of the lgprobe subcircuit.
  - GEAR integration is refused (same TRAPEZOIDAL-orders-1-2-only scope as
    PSSsensRecordStep, and for the same reason -- see that function).
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/pssdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/fteext.h"
#include "ngspice/smpdefs.h"
#include "ngspice/devdefs.h"
#include "ngspice/noisedef.h"  /* card #509: NOISEAN/Ndata for pac_engine_run_adjoint_noise() */
#include "ngspice/cpextern.h"  /* card #509: cp_getvar() for the `sqrnoise` flag */

#include "cap/capdefs.h"
#include "ind/inddefs.h"
#include "vsrc/vsrcdefs.h"
#include "sw/swdefs.h"   /* card #559: switch power nodes with device-model charge */
#include "dio/diodefs.h"   /* card #581: diode commutation note */
#include "isrc/isrcdefs.h"
#ifdef XSPICE
#include "ngspice/miftypes.h"
#include "ngspice/mif.h"
#include "ngspice/mifdefs.h"
#include "ngspice/mifcmdat.h"
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pssshoot.h"
#include "pacshoot.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* small complex helpers */

typedef struct { double re, im; } PACc;

extern int PSSdiscardRunData;   /* card #553, pssshoot.c */
static PACc pc_add(PACc a, PACc b) { PACc r = { a.re + b.re, a.im + b.im }; return r; }
static PACc pc_sub(PACc a, PACc b) { PACc r = { a.re - b.re, a.im - b.im }; return r; }
static PACc pc_mul(PACc a, PACc b)
{
    PACc r;
    r.re = a.re * b.re - a.im * b.im;
    r.im = a.re * b.im + a.im * b.re;
    return r;
}
static PACc pc_scale(PACc a, double s) { PACc r = { a.re * s, a.im * s }; return r; }
static PACc pc_div(PACc a, PACc b)
{
    PACc r;
    double d = b.re * b.re + b.im * b.im;
    if (d == 0.0) d = 1e-300;
    r.re = (a.re * b.re + a.im * b.im) / d;
    r.im = (a.im * b.re - a.re * b.im) / d;
    return r;
}
static double pc_mag(PACc a) { return sqrt(a.re * a.re + a.im * a.im); }
static double pc_ph(PACc a) { return atan2(a.im, a.re); }
static PACc pc_conjexp(double phase) { PACc r; r.re = cos(phase); r.im = -sin(phase); return r; }

/* ------------------------------------------------------------------ */
/* Reduced state map -- independently rebuilt, SAME device-list order as
 * pssshoot.c's (file-static) PSSstateMap, so x0[] indices line up. */

typedef struct {
    int isCap;
    union {
        CAPinstance *cap;
        INDinstance *ind;
    } dev;
} PACstateRef;

typedef struct {
    int n;               /* native reduced states (x0 / fT layout, same as pssshoot.c) */
    PACstateRef *refs;
    /* card #542: GENERIC (charge-based) state set. When the circuit has device-model
     * charge storage (MOS/BSIM/diode/BJT/... or mutual inductors) the Phi/xi/z state
     * vector is the set of matrix rows (node voltages / branch currents) that appear
     * as a COLUMN of the charge matrix C = dq/dx anywhere on the orbit (found by
     * pac_generic_prepare()'s prepass); n above stays the native count only (it still
     * sizes x0 and sets the recording runs' initial conditions). */
    int generic;
    int ns;
    int *srows;          /* [ns] ascending 1-based rows */
} PACstateMap;

#define PAC_NST(map) ((map)->generic ? (map)->ns : (map)->n)

static void
pac_map_free(PACstateMap *map)
{
    if (map->refs) tfree(map->refs);
    if (map->srows) tfree(map->srows);
    map->refs = NULL; map->srows = NULL;
    map->generic = 0; map->ns = 0;
}

static int
pac_build_map(CKTcircuit *ckt, PACstateMap *map)
{
    int captype = CKTtypelook("Capacitor");
    int indtype = CKTtypelook("Inductor");
    int n = 0;
    CAPmodel *cm; CAPinstance *ci;
    INDmodel *im; INDinstance *ii;

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
    map->generic = 0; map->ns = 0; map->srows = NULL;
    if (n == 0)
        return OK;
    map->refs = TMALLOC(PACstateRef, n);
    if (!map->refs)
        return E_NOMEM;

    n = 0;
    if (captype >= 0)
        for (cm = (CAPmodel *) ckt->CKThead[captype]; cm; cm = CAPnextModel(cm))
            for (ci = CAPinstances(cm); ci; ci = CAPnextInstance(ci)) {
                map->refs[n].isCap = 1;
                map->refs[n].dev.cap = ci;
                n++;
            }
    if (indtype >= 0)
        for (im = (INDmodel *) ckt->CKThead[indtype]; im; im = INDnextModel(im))
            for (ii = INDinstances(im); ii; ii = INDnextInstance(ii)) {
                map->refs[n].isCap = 0;
                map->refs[n].dev.ind = ii;
                n++;
            }
    return OK;
}

static void
pac_state_set_ic(PACstateMap *map, const double *x)
{
    int k;
    for (k = 0; k < map->n; k++) {
        if (map->refs[k].isCap)
            map->refs[k].dev.cap->CAPinitCond = x[k];
        else
            map->refs[k].dev.ind->INDinitCond = x[k];
    }
}

static int
pac_state_row(const PACstateRef *ref)
{
    return ref->isCap ? ref->dev.cap->CAPposNode : ref->dev.ind->INDbrEq;
}

/* row of state d in the Phi/xi/z state vector (card #542: generic maps list rows directly) */
static int
pac_map_row(const PACstateMap *map, int d)
{
    return map->generic ? map->srows[d] : pac_state_row(&map->refs[d]);
}

/* card #489 SECOND pass (plan T6.2, task 5): resolve a `trigger=` node
 * (its raw circuit node number, job->PSStrigNode->number in pssshoot.c) to
 * an index into THIS map's reduced-state vector, for the free-running
 * period-variation projection (see pacshoot.h's pac_pstb_run_ex() comment
 * and docs/pss-pac.md's free-running section for the derivation). Only a
 * native capacitor's OWN positive terminal is handled -- the paradigm case
 * for every hysteretic/bang-bang/COT reference this project has (the
 * trigger is the latch/comparator-state capacitor itself, e.g. R5's `q`),
 * and the only case pss_run_until_trigger()'s own f(T) estimate lines up
 * with exactly (a reduced STATE, not a derived algebraic node). Returns
 * -1 (not -- e.g. -- 0, a valid index) if trig_node_num < 0 (ordinary
 * driven call) or no match is found (an inductor-current trigger, or a
 * trigger node not directly a reduced state -- e.g. separated from the
 * nearest cap by a series resistance) -- the caller must treat -1 as "no
 * correction available" and refuse rather than silently running the
 * uncorrected engine on what would then be a free-running orbit. */
static int
pac_find_trig_idx(const PACstateMap *map, int trig_node_num)
{
    int d;
    if (trig_node_num < 0)
        return -1;
    if (map->generic) {   /* card #544: the trigger node must be one of the state rows (a column of C somewhere on the orbit) */
        for (d = 0; d < map->ns; d++)
            if (map->srows[d] == trig_node_num)
                return d;
        return -1;
    }
    for (d = 0; d < map->n; d++)
        if (map->refs[d].isCap && map->refs[d].dev.cap->CAPposNode == trig_node_num)
            return d;
    return -1;
}

/* ------------------------------------------------------------------ */
/* card #465 SECOND pass (plan T5.3, task 3 "speed"): PAC_RECORD_POINTS
 * decouples the grid resolution of pacshoot.c's OWN two recording-period
 * reruns (pac_engine_run's "pass 1"/"pass 2", the ONLY place this file
 * calls pac_run_period) from job->PSSpoints, the grid pssshoot.c's
 * shooting-Newton convergence phase (fast-state detection, the Newton
 * iterations themselves, and the Floquet FD monodromy) ALSO uses via the
 * SAME ckt->CKTmaxStep/CKTstep.
 *
 * PROFILING FINDING (this pass): on pac_lti.cir (task 2's LTI accuracy
 * check, Points=200000 to reach the 1e-6 pac_sb0(a) target), of the
 * reported 113.4s total wall time, PSS-SHOOT's OWN "CONVERGED ...
 * wall_time=" line accounts for 112.4s (99%) and PAC-TD's OWN "done,
 * wall_time=" line for 1.2s (1%) -- see the card #465 report. The task's
 * "~110s for 31 frequency points" is therefore NOT this file's frequency
 * sweep (which is already O(2 period reruns) total, independent of the
 * number of frequencies, per this file's header design) -- it is
 * pssshoot.c's own ~10-12 period reruns (fast-state detection's FD
 * monodromy, the Newton iteration(s), the Floquet FD monodromy) each
 * paying the SAME Points=200000 resolution that only PAC's OWN
 * demodulation quadrature (task 2) actually needs. pssshoot.c is another
 * card's file (not edited here); it is not obviously safe or in scope to
 * make ITS OWN convergence phase run at a coarser grid than the deck's
 * Points= asks for.
 *
 * FIX (in scope: this file only): pac_run_period() OVERRIDES
 * ckt->CKTmaxStep/CKTstep/CKTdelmin to T/recPoints for JUST this one
 * recording-period call, restoring the PSS-converged values immediately
 * after. recPoints is resolved in this order (card #465 THIRD pass, task
 * 2 -- see pssdefs.h's PACrecordPoints comment):
 *   1. job->PACrecordPoints, if a `.pac`/`.pstb` card set `pacpoints=N`
 *      (N > 0) -- the primary way to ask for this now;
 *   2. else the PAC_RECORD_POINTS env var, if set (DEPRECATED: kept only
 *      as an override for scripts/decks written before pacpoints=
 *      existed -- see this project's docs for the deprecation note);
 *   3. else job->PSSpoints -- the SAME points= a deck already gave PSS's
 *      own convergence phase, so a deck that sets neither `pacpoints=`
 *      nor the env var gets an EXPLICIT override that is numerically
 *      identical to the second pass's "no override, inherit whatever
 *      ckt->CKTmaxStep PSS's Newton phase last left behind" default
 *      (PSS's own pssshoot.c also derives CKTmaxStep from T/PSSpoints),
 *      so old decks see byte-identical behaviour, now for a documented
 *      reason instead of by inheritance.
 * This is safe regardless of which source wins: (1) BOTH of
 * pac_engine_run's two recording passes for a given (br,bi,freqs) sweep
 * resolve recPoints the SAME way (same job, same call site pattern), so
 * their step grids still match each other exactly (the property
 * PSSpacRecordStep's per-step recursion needs -- see this file's header
 * -- is "pass 1 and pass 2 use the same grid as EACH OTHER", not "the
 * same grid PSS's Newton phase used"); (2) x0 (the converged orbit) and T
 * (the converged period) are UNCHANGED -- only the recording reruns'
 * internal step density changes. Lets a deck ask PSS for a COARSE, fast-
 * converging Points= (job->PSSpoints drives only the Newton/Floquet
 * reruns) while pacpoints= independently asks for the fine grid task 2's
 * demodulation accuracy needs -- see the card #465 report for the
 * measured before/after on pac_lti.cir. */
/* card #542: periodic-continuation start state of every recording run of the GENERIC path: the node/branch
 * vector at the END of the prepass period. DCtran's UIC start takes ckt->CKTrhs as the initial solution (NIiter
 * SWAPs it into CKTrhsOld) and the INITJCT load that follows initialises every device's charge history from the
 * DEVICE-level initial conditions (MOS icVDS/icVGS/icVBS, diode icVD, BJT icVBE/icVCE, ... derived by CKTic()
 * once, at job start, from zero node ICs). Left alone, nodes carried only by model charges start every recorded
 * period far from their periodic value: the base trajectory then has a startup transient the orbit does not
 * have (refs/README.md R7: "first ~15 ns of the orbit wrong") and the linearisation along it is off by about
 * (startup duration)/T (measured: MOS1 amplifier DC gain 2.7e-4, with the device ics re-derived 6e-8).
 * pac_run_period() therefore restores this vector and re-derives the device ics before each run. */
static double *g_pac_xinit = NULL;
static int g_pac_xinit_n = 0;

/* card #625: the recording passes used to inherit the PSS reltol (1e-5). Their LTE step grid then moves
 * discretely with any parameter (measured on rd_pstb_r3d: a 1e-10 relative change of a capacitor gave a new grid
 * every time, |T| scattered by 4e-3 rms and up to 0.15 dB, PM by 0.056 deg rms; the PSS orbit and the Newton
 * residual were not the cause: with the orbit fixed, the recording reltol alone moved |T(fc)| by 0.7 dB at 1e-4).
 * Rms of |T| against the recording reltol: 1e-5 4.0e-3, 1e-6 1.7e-4, 1e-7 3.8e-5, 1e-8 5.2e-7. With no pacreltol=
 * and no PAC_RECORD_RELTOL, the recording reltol is now min(inherited, rung), the rungs 1e-8, 1e-7, 1e-6 and
 * then the inherited value: a rung whose run stops with "timestep too small" (rd_pstb_r4d's dead-time edge at
 * <= 1e-7) is retried at the next one, with a note. PAC_RECORD_RELTOL_AUTO=0 restores the inherited value. */
static double g_pac_auto_reltol = 0.0;
static int g_pac_rec_tstep = 0;   /* card #625: a recording run of this attempt stopped with E_TIMESTEP */
/* card #625: the tighter recording reltol is for circuits without switch events. On ideal-switch decks the event
 * corrections (#559-#586: saltation, charge jumps, fast nodes, routed jumps) were derived and validated against the
 * perturbed-PSS references on the grid of the inherited reltol; on the 1e-8 grid two of them moved from 0.018 to
 * 0.12 dB and three red proofs stopped failing (s134 portfolio, 2026-10-14), so those decks keep the inherited
 * value: an ideal switch (SW, CSW) or any XSPICE code model in the circuit. */
static int
pac_ckt_has_events(CKTcircuit *ckt)
{
    int i;
    for (i = 0; i < DEVmaxnum; i++) {
        GENmodel *m;
        const char *dn;
        int ev;
        if (!DEVices[i] || !ckt->CKThead[i])
            continue;
        dn = DEVices[i]->DEVpublic.name;
        ev = dn && (!strcmp(dn, "Switch") || !strcmp(dn, "CSwitch"));
#ifdef XSPICE
        ev = ev || DEVices[i]->DEVpublic.cm_func != NULL;   /* a code model (DEVload is the .cm library's MIFload wrapper, not MIFload) */
#endif
        if (!ev)
            continue;
        for (m = ckt->CKThead[i]; m; m = m->GENnextModel)
            if (m->GENinstances)
                return 1;
    }
    return 0;
}

static int g_hid_nc, g_hid_nhy;   /* card #557, defined below with the hidden-state machinery */
static int pac_run_period_inner(CKTcircuit *ckt, PSSan *job, double T, const char *jobname);
static int
pac_run_period(CKTcircuit *ckt, PSSan *job, double T, const char *jobname)
{
    static int on = -1;
    double t0;
    int err;
    if (on < 0) on = getenv("PSS_RUN_STATS") != NULL;
    if (!on) return pac_run_period_inner(ckt, job, T, jobname);
    t0 = PSSrunStatsNow();
    err = pac_run_period_inner(ckt, job, T, jobname);
    PSSrunStatsAdd(jobname, PSSrunStatsNow() - t0);
    return err;
}

static int
pac_run_period_inner(CKTcircuit *ckt, PSSan *job, double T, const char *jobname)
{
    JOB *saveJob = ckt->CKTcurJob;
    TRANan localjob;
    char namebuf[128];
    int err;
    double saveMaxStep = ckt->CKTmaxStep, saveStep = ckt->CKTstep,
           saveDelmin = ckt->CKTdelmin;
    long recPoints;
    int overriding;

    if (job->PACrecordPoints > 0) {
        recPoints = job->PACrecordPoints;
    } else {
        const char *rpEnv = getenv("PAC_RECORD_POINTS");
        if (rpEnv && atol(rpEnv) > 0) {
            static int warned = 0;
            recPoints = atol(rpEnv);
            if (!warned) {
                fprintf(stderr, "PAC-TD: PAC_RECORD_POINTS is deprecated -- "
                        "use `pacpoints=%ld` on the .pac/.pstb card instead.\n",
                        recPoints);
                warned = 1;
            }
        } else {
            recPoints = job->PSSpoints;
        }
    }
    overriding = (recPoints > 0);

    memset(&localjob, 0, sizeof localjob);
    strncpy(namebuf, jobname, sizeof namebuf - 1);
    namebuf[sizeof namebuf - 1] = '\0';
    localjob.JOBtype = saveJob->JOBtype;
    localjob.JOBname = namebuf;

    if (overriding) {
        ckt->CKTmaxStep = T / (double) recPoints;
        ckt->CKTstep    = ckt->CKTmaxStep;
        ckt->CKTdelmin  = 1e-9 * ckt->CKTmaxStep;
    }

    ckt->CKTcurJob   = (JOB *) &localjob;
    ckt->CKTfinalTime = T;
    ckt->CKTinitTime  = 0.0;
    ckt->CKTtime      = 0.0;
    ckt->CKTmode     |= MODEUIC;

    if (g_pac_xinit && g_pac_xinit_n == SMPmatSize(ckt->CKTmatrix) + 1) {
        int di;
        memcpy(ckt->CKTrhs, g_pac_xinit, (size_t) g_pac_xinit_n * sizeof(double));
        /* device-level initial conditions (MOS icVDS/icVGS/icVBS, diode icVD, BJT icVBE/icVCE, ...) are what the
         * UIC start's INITJCT load uses for the history of every model charge; CKTic() derived them from the
         * (zero) node ICs once at job start. Re-derive them from the periodic state (devices whose ic the
         * deck gave explicitly keep it). Native C/L ICs are set by pac_state_set_ic() from the PSS orbit. */
        for (di = 0; di < DEVmaxnum; di++) {
            const char *dn;
            if (!DEVices[di] || !DEVices[di]->DEVsetic || !ckt->CKThead[di])
                continue;
            dn = DEVices[di]->DEVpublic.name;
            if (dn && (!strcmp(dn, "Capacitor") || !strcmp(dn, "Inductor")))
                continue;
            DEVices[di]->DEVsetic(ckt->CKThead[di], ckt);
        }
    }
    PSSeventLocateActive = PSSswLagActive = PSSshootLagEnabled();
    {
        double save_reltol = ckt->CKTreltol;
        if (job->PACrecordReltol > 0.0)
            ckt->CKTreltol = job->PACrecordReltol;
        else if (getenv("PAC_RECORD_RELTOL") && atof(getenv("PAC_RECORD_RELTOL")) > 0.0)
            ckt->CKTreltol = atof(getenv("PAC_RECORD_RELTOL"));   /* card #546: env form, tests */
        else if (g_pac_auto_reltol > 0.0 && g_pac_auto_reltol < ckt->CKTreltol)
            ckt->CKTreltol = g_pac_auto_reltol;                   /* card #625: default rung of the ladder */
#ifdef XSPICE
        /* card #551/#557: stateless models: no state to carry, model time 0.  With stateful analog models the carried hidden state
         * of the converged orbit (analog CKTstate values, model time) is put back, as for a restarted shooting run. */
        if (g_hid_nc > 0 || g_hid_nhy > 0) PSSxsPacBegin(ckt); else g_xs_toff = 0.0;
#endif
        /* card #553: the recording pass's plot is never read; keep one point per vector (PSS_KEEP_RUN_DATA=1: old behaviour) */
        PSSdiscardRunData = !(getenv("PSS_KEEP_RUN_DATA") && atoi(getenv("PSS_KEEP_RUN_DATA")) != 0);
        err = DCtran(ckt, 1);
        if (err == E_TIMESTEP) g_pac_rec_tstep = 1;
        PSSdiscardRunData = 0;
        ckt->CKTreltol = save_reltol;
#ifdef XSPICE
        if (g_hid_nc > 0 || g_hid_nhy > 0) PSSxsPacEnd();
#endif
    }
    PSSeventLocateActive = PSSswLagActive = 0;

    if (localjob.TRANplot)
        SPfrontEnd->OUTendPlot(localjob.TRANplot);

    ckt->CKTcurJob = saveJob;
    if (overriding) {
        ckt->CKTmaxStep = saveMaxStep;
        ckt->CKTstep    = saveStep;
        ckt->CKTdelmin  = saveDelmin;
    }
    return err;
}

/* ------------------------------------------------------------------ */
/* PACbatch: the per-step recursion state, generalising pssshoot.c's
 * analytic-sensitivity S_k/T_k recursion (see this file's header) to M
 * simultaneous complex, possibly-forced columns. */

int PSSpacActive = 0;

static PACstateMap *g_pac_map = NULL;
static int g_pac_n = 0, g_pac_neq = 0, g_pac_M = 0, g_pac_failed = 1;
static PACc *g_pac_Cprev = NULL, *g_pac_Cnew = NULL;    /* [(neq+1)*M] */
static PACc *g_pac_Tprev = NULL, *g_pac_Tnew = NULL;    /* [n*M] */
static double *g_pac_omega = NULL;                       /* [M] */
static int *g_pac_useforce = NULL;                       /* [M] */
static double *g_pac_br = NULL, *g_pac_bi = NULL;        /* [neq+1], shared forcing pattern */
static double g_pac_tprev = 0.0;
static int g_pac_demod_active = 0;
static PACc *g_pac_demod = NULL;                         /* [neq*M], sideband-0 accumulator */
static double g_pac_T = 0.0;                              /* period, for the 1/T normalisation */
static double *g_pac_rhs_re = NULL, *g_pac_rhs_im = NULL, *g_pac_scratch = NULL; /* [neq+1] */
static int g_pac_exact_jac = 0;   /* card #465 THIRD pass, task 1: see pac_refactor_at_accepted() */

/* ------------------------------------------------------------------ */
/* card #542: GENERIC CHARGE-BASED HISTORY.
 *
 * The linearised transient recursion is, for EVERY charge state (native C/L
 * AND device-model charges: MOS/BSIM4 intrinsic+overlap+junction, diode cj/tt,
 * BJT, JFET, mutual inductance, ...),
 *
 *     (G_k + ag0 C_k) y_k = b_k + ag0 Q_{k-1} + a1 T_{k-1}
 *     Q_k = C_k y_k            (dq, a full row vector: row r holds dq entering KCL row r)
 *     T_k = ag0 (Q_k - Q_{k-1}) - a1 T_{k-1}      (d qdot; a1 = ag1 for order 2, 0 for order 1)
 *
 * with J_k = G_k + ag0 C_k the loaded Jacobian and C_k = dq/dx the charge matrix of
 * the SAME load.  C_k is obtained generically, with no per-device code: every
 * charge-conserving device stamps its capacitive part through NIintegrate as
 * geq = ag0*cap, so the loaded matrix is exactly affine in ag0, and
 * C_k = (J(ag0) - J(ag0/2)) / (ag0/2) entry by entry.  A third load at 2*ag0
 * checks the affinity (device loads that are not affine in ag0, e.g. NQS models,
 * are REFUSED naming the device type).  The legacy native-only history (cap and
 * ind instances only: Cprev/Tprev per reduced state) is kept for circuits without
 * model charges (bit-identical to before) and under PAC_NATIVE_HISTORY=1 (red proof).
 * The state vector of Phi/xi/z is the set of rows that are a column of some C_k
 * (pac_generic_prepare()). */
typedef struct { int n, cap; int *ri, *ci; double *v; } PACtrip;

static int g_pac_gen = 0;                        /* generic recursion active for this pass */
static PACc *g_pac_Qprev = NULL, *g_pac_Qnew = NULL;   /* [(neq+1)*M] */
/* card #549: method=gear on the generic recursion.  The charge derivative of a BDF/Gear step of order p is
 * dq/dt|k = ag[0] Q_k + sum_{j=1..p} ag[j] Q_{k-j}, so  J(ag0) y_k = -sum_j ag[j] Q_{k-j}  (+ forcing): no derivative history T,
 * the charge history itself is kept (Q_{k-1} = g_pac_Qprev, Q_{k-2..k-6} = g_pac_Qh[2..6]).  The order of the step is the highest
 * nonzero ag[j] of the step just accepted (ckt->CKTorder is not trustworthy at this call site, as for trap).  Order 1 is
 * backward Euler, ag[1] = -ag[0]: the same formula.  Trap keeps its T recursion bit for bit. */
#define PAC_GEAR_MAXORD 6
static int g_pac_gear = 0;
static PACc *g_pac_Qh[PAC_GEAR_MAXORD + 1];   /* [2..6] allocated when gear */
static PACtrip g_pac_C = { 0, 0, NULL, NULL, NULL };   /* C_k of the current step */
static PACtrip g_pac_J1 = { 0, 0, NULL, NULL, NULL }, g_pac_J2 = { 0, 0, NULL, NULL, NULL },
               g_pac_J3 = { 0, 0, NULL, NULL, NULL };
static int g_pac_gen_warm = 0;                   /* pre-iterate loads to settle device limiting state (backward sweeps) */
static int g_pac_gen_check = 0;                  /* do the third-ag0 affinity check this pass */
static int g_pac_gen_bad = 0;                    /* nonaffine load seen */
static char g_pac_gen_badname[96] = "";
/* prepass (state-set discovery) */
static int g_pp_active = 0, g_pp_neq = 0, g_pp_steps = 0;
static unsigned char *g_pp_col = NULL;           /* [neq+1] 1 = column of some C_k */
static double g_pp_cmax = 0.0;
static void pac_gen_identify_bad(CKTcircuit *ckt, double ag0);
static int g_pp_susp[64], g_pp_nsusp = 0;       /* non-native device types with a DEVload */
static int g_pp_has_mutual = 0;                  /* a mutual inductor (K) is present */
static int g_pp_detect = 0, g_pp_detect_hit = 0; /* model-charge detection prepass state */
static void pac_detect_step(CKTcircuit *ckt);
static double *g_pp_xend = NULL;

/* card #557: HIDDEN-STATE td small-signal recursion for XSPICE analog code models with state (s_xfer, int, d_dt, slew, hyst).
 * The model's CKTstate slots (cm_analog_alloc'd: integrator outputs, old inputs, ...) are extra Floquet states h.  One accepted
 * step of the real transient is a map (x_k, h_k) = S(x_{k-1}..., h_{k-1}); linearised at the accepted point it reads
 *       J'_k y_k + B_k dh_{k-1} = (forcing + capacitive history),      dh_k = P_k y_k + R_k dh_{k-1}
 * where B = dF/dh1 (F = J x - rhs is the load residual, state1 slot perturbed), P = dh0/dx, R = dh0/dh1 and J' is J with the columns
 * of the model's input nodes replaced by the finite-difference dF/dx (the models' declared PARTIALs are not always the derivative:
 * s_xfer's PARTIAL(out,in) is the value of an integrator, not a partial).  Every matrix comes from central differences of
 * the model through full CKTload() calls, so no per-model code exists.  The slot set is found in the prepass: slots reachable
 * from the inputs (P or R chain) AND read by a later step (B or R chain) -- constants (coefficients copied state1 -> state0, time
 * stamps) and write-only slots (in/out copies) are dropped.  cm_analog_integrate uses g_xs_ag0_fix while ag0 is varied for C_k
 * so the models do not enter the charge matrix. */
static int g_hid_nc = 0;          /* candidate slots (0: no stateful XSPICE model) */
static int g_hid_nin = 0;         /* input rows of the stateful instances */
static int g_hid_nh = 0;          /* ACTIVE hidden states (virtual rows neq+1 .. neq+nh of Cprev/Cnew) */
static int g_hid_mode = 0;        /* 0 off, 1 prepass pattern, 2 forward (matrices per step) */
static int g_hid_slewact = 0;
static int g_hid_stride = 1;      /* prepass: pattern probed every stride-th step (all steps up to ~20000 per period) */
static char *g_hid_cslew = NULL;  /* candidate slot belongs to a 'slew' model: slew-rate limiting active on the orbit is refused (saltation) */
static int *g_hid_cslot = NULL, *g_hid_inrow = NULL, *g_hid_act = NULL;   /* candidate state indices; input rows; active -> candidate */
static double *g_hid_cscale = NULL;
static unsigned char *g_hid_patP = NULL, *g_hid_patR = NULL, *g_hid_patB = NULL;
static double *g_hid_Bh = NULL, *g_hid_P = NULL, *g_hid_R = NULL, *g_hid_Jc = NULL;   /* [(neq+1)*nh], [nh*nin], [nh*nh], [(neq+1)*nin] */
static int g_hid_nrow = 0;        /* neq + nh while the generic recursion runs */
static int g_hid_flips = 0;       /* discrete (Boolean) model state changes seen on the orbit */
#ifdef XSPICE
static int pac_hid_stage(CKTcircuit *ckt, double ag0, const double *xbase, const double *st0);
static PACtrip g_hid_J = { 0, 0, NULL, NULL, NULL };   /* scratch matrix of the finite-difference loads */
#endif



/* card #504: the true-adjoint engine's OWN forward Phi/G passes need the
 * SAME per-step J_k convention its backward sweep rebuilds (J_k at the
 * TRUE accepted x_k, not the legacy default's stale x_{k-1} -- see
 * pac_engine_run_adjoint()'s design comment for why this is NOT optional
 * the way PAC_EXACT_JACOBIAN is for v1). Independent of g_pac_exact_jac
 * (which pac_reset() unconditionally re-reads from the env var on every
 * call, so a plain assignment before pac_reset() would just be
 * overwritten) -- OR'd in at PSSpacRecordStep()'s own check, set/cleared
 * ONLY around the adjoint engine's own forward pac_run_period() calls. */
static int g_pac_force_exact_jac = 0;

/* card #465 FOURTH pass, task 1: the true y(0) at a NON-state row (any row
 * that is not a capacitor/inductor terminal the reduced state tracks, e.g.
 * a resistor node or a probe's branch equation) is NOT zero in general --
 * see pac_finish_demod_wraparound()'s header for the periodicity argument
 * that recovers it exactly, with no extra algebraic solve, from y(T) once
 * the period run has completed. g_pac_first_y/g_pac_first_t/g_pac_have_first
 * stash the FIRST accepted step's (t_1, y_1) so that closing computation
 * can be done after the loop. [(neq+1)*M], same layout as g_pac_Cprev. */
static PACc *g_pac_first_y = NULL;
static double g_pac_first_t = 0.0;
static int g_pac_have_first = 0;
static int g_pac_xls_have;   /* card #579, defined with pac_imp_roll() */

/* ------------------------------------------------------------------ */
/* card #504: TRAJECTORY RECORDER for the true-adjoint `.pxf ... adjoint=1`
 * engine (pac_engine_run_adjoint() below). Independent of g_pac_Cprev/etc
 * above (own reset/free, own lifetime -- spans MULTIPLE pac_reset()/
 * pac_free() cycles: the adjoint engine's own "Phi+trajectory" forward
 * pass uses the ordinary g_pac_* column machinery for Phi exactly like
 * pac_engine_run()'s pass 1 homogeneous columns, and ADDITIONALLY records
 * one entry here per accepted step -- see PSSpacRecordStep()'s own
 * g_pac_traj_active check, added right after order/t_new/h are computed).
 *
 * What is stored per accepted step k: the accepted point x_k itself
 * (ckt->CKTrhs, the TRUE just-accepted iterate -- SAME array
 * pac_refactor_at_accepted() already trusts is correct at this exact call
 * site, see that function's header) plus t_k/h_k/ag0_k/ag1_k/order_k --
 * everything pac_adj_refactor_at() (below, near the backward sweep) needs
 * to rebuild J_k = G(x_k)+ag0_k*C(x_k) LATER, in REVERSE order, via a
 * fresh CKTload()+SMPluFac() outside any live transient (no
 * save/restore-for-the-ongoing-transient's-sake concern the live-
 * interleaved pac_refactor_at_accepted() has -- there is no ongoing
 * transient by the time the backward sweep runs; pac_engine_run_adjoint()
 * itself saves/restores the handful of ckt-global scalars this touches,
 * see that function).
 *
 * Memory: O(N*(neq+2)) doubles, N = accepted steps in ONE period at
 * whatever recording resolution (PACrecordPoints/pacpoints=, same knob
 * pac_run_period() already uses) -- see the card #504 report's own
 * measured numbers on the three test_pxf.py fixtures. Grown by doubling;
 * never shrunk mid-job (freed only by pac_traj_free()). */
typedef struct {
    double t, h, ag0, ag1;
    int order;
    double *x;   /* [neq+1], 1-based, OWNED (malloc'd) copy of ckt->CKTrhs */
    /* card #542 (generic history only): the EXACT J_k (every stored matrix element) and charge
     * matrix C_k of the forward pass, kept so the backward sweep transposes precisely the matrices the
     * forward recursion used. Reloading them from the devices in reverse order is not exact for
     * devices whose load reads CKTstate1 (MOS1 Meyer capacitance averaging uses the previous step's
     * value) or limiting state. */
    int nJ, nC;
    int *Jri, *Jci; double *Jv;
    int *Cri, *Cci; double *Cv;
    /* card #544 (noise only): the converged CKTstate0 of the forward pass at this step.  DEVnoise of DIOnoise/BJTnoise/MOS3noise/
     * BSIM4noise reads the operating point (diode and BJT currents, vds/vgs) from CKTstate0, which the backward sweep must
     * therefore restore to the forward value at x_k (an iterated warm load from the stale end-of-period state does not always
     * converge through the junction limiting: a few steps then carried currents 10^6 too large). */
    double *st0; int nst0;
} PACtrajStep;
static int g_pac_traj_keep_st0 = 0;
static double g_pac_traj_bytes = 0.0;

static int g_pac_traj_active = 0;     /* PSSpacRecordStep() appends when set */
static PACtrajStep *g_pac_traj = NULL;
static int g_pac_traj_n = 0, g_pac_traj_cap = 0;
static int g_pac_traj_neq = 0;        /* neq at record time, for x[] sizing */
static int g_pac_traj_failed = 0;

/* card #544: one line on what the stored trajectory costs (PAC_GEN_TRACE=1 or PAC_TRAJ_REPORT=1). */
static void
pac_traj_report(const char *who)
{
    if (!getenv("PAC_GEN_TRACE") && !getenv("PAC_TRAJ_REPORT"))
        return;
    fprintf(stderr, "%s: stored trajectory %d steps, J_k/C_k %.1f MB + x_k %.1f MB (cap %s MB)\n", who, g_pac_traj_n,
            g_pac_traj_bytes / 1048576.0, (double) g_pac_traj_n * (double) (g_pac_traj_neq + 1) * 8.0 / 1048576.0,
            getenv("PAC_TRAJ_MAX_MB") ? getenv("PAC_TRAJ_MAX_MB") : "3000");
}

static void
pac_traj_free(void)
{
    int i;
    if (g_pac_traj) {
        for (i = 0; i < g_pac_traj_n; i++) {
            if (g_pac_traj[i].x) tfree(g_pac_traj[i].x);
            if (g_pac_traj[i].Jri) tfree(g_pac_traj[i].Jri);
            if (g_pac_traj[i].Jci) tfree(g_pac_traj[i].Jci);
            if (g_pac_traj[i].Jv) tfree(g_pac_traj[i].Jv);
            if (g_pac_traj[i].Cri) tfree(g_pac_traj[i].Cri);
            if (g_pac_traj[i].Cci) tfree(g_pac_traj[i].Cci);
            if (g_pac_traj[i].Cv) tfree(g_pac_traj[i].Cv);
            if (g_pac_traj[i].st0) tfree(g_pac_traj[i].st0);
        }
        tfree(g_pac_traj);
    }
    g_pac_traj = NULL;
    g_pac_traj_n = g_pac_traj_cap = 0;
    g_pac_traj_neq = 0;
    g_pac_traj_active = 0;
    g_pac_traj_failed = 0;
    g_pac_traj_bytes = 0.0;
}

static void
pac_traj_reset(int neq)
{
    pac_traj_free();
    g_pac_traj_neq = neq;
    g_pac_traj_active = 1;
}

/* Called from PSSpacRecordStep() at every accepted step (t>0) while
 * g_pac_traj_active -- appends (t,h,ag0,ag1,order, a COPY of ckt->CKTrhs)
 * to g_pac_traj[]. Sets g_pac_traj_failed on OOM (checked by the adjoint
 * engine before trusting the buffer). */
static void
pac_traj_append(CKTcircuit *ckt, double t, double h, double ag0, double ag1, int order)
{
    double *xcopy;
    if (g_pac_traj_failed)
        return;
    if (g_pac_traj_n >= g_pac_traj_cap) {
        int newcap = g_pac_traj_cap ? g_pac_traj_cap * 2 : 256;
        PACtrajStep *nu = TREALLOC(PACtrajStep, g_pac_traj, newcap);
        if (!nu) { g_pac_traj_failed = 1; return; }
        g_pac_traj = nu;
        g_pac_traj_cap = newcap;
    }
    xcopy = TMALLOC(double, g_pac_traj_neq + 1);
    if (!xcopy) { g_pac_traj_failed = 1; return; }
    memcpy(xcopy, ckt->CKTrhs, (size_t) (g_pac_traj_neq + 1) * sizeof(double));
    g_pac_traj[g_pac_traj_n].t = t;
    g_pac_traj[g_pac_traj_n].h = h;
    g_pac_traj[g_pac_traj_n].ag0 = ag0;
    g_pac_traj[g_pac_traj_n].ag1 = ag1;
    g_pac_traj[g_pac_traj_n].order = order;
    g_pac_traj[g_pac_traj_n].x = xcopy;
    g_pac_traj[g_pac_traj_n].nJ = g_pac_traj[g_pac_traj_n].nC = 0;
    g_pac_traj[g_pac_traj_n].Jri = g_pac_traj[g_pac_traj_n].Jci = NULL;
    g_pac_traj[g_pac_traj_n].Jv = NULL;
    g_pac_traj[g_pac_traj_n].Cri = g_pac_traj[g_pac_traj_n].Cci = NULL;
    g_pac_traj[g_pac_traj_n].Cv = NULL;
    g_pac_traj[g_pac_traj_n].st0 = NULL;
    g_pac_traj[g_pac_traj_n].nst0 = 0;
    if (g_pac_traj_keep_st0 && ckt->CKTnumStates > 0) {
        double *sc = TMALLOC(double, ckt->CKTnumStates);
        if (!sc) { g_pac_traj_failed = 1; return; }
        memcpy(sc, ckt->CKTstate0, (size_t) ckt->CKTnumStates * sizeof(double));
        g_pac_traj[g_pac_traj_n].st0 = sc;
        g_pac_traj[g_pac_traj_n].nst0 = ckt->CKTnumStates;
        g_pac_traj_bytes += (double) ckt->CKTnumStates * 8.0;
    }
    g_pac_traj_n++;
}

static PACc
pac_get(const PACc *C, int row, int M, int m)
{
    if (row <= 0 || !C) { PACc z = { 0.0, 0.0 }; return z; }
    return C[(size_t) row * M + m];
}
static void
pac_set(PACc *C, int row, int M, int m, PACc v)
{
    if (row <= 0 || !C) return;
    C[(size_t) row * M + m] = v;
}

static void
pac_free(void)
{
    if (g_pac_Cprev) tfree(g_pac_Cprev);
    if (g_pac_Cnew)  tfree(g_pac_Cnew);
    if (g_pac_Tprev) tfree(g_pac_Tprev);
    if (g_pac_Tnew)  tfree(g_pac_Tnew);
    if (g_pac_omega) tfree(g_pac_omega);
    if (g_pac_useforce) tfree(g_pac_useforce);
    if (g_pac_demod) tfree(g_pac_demod);
    if (g_pac_rhs_re) tfree(g_pac_rhs_re);
    if (g_pac_rhs_im) tfree(g_pac_rhs_im);
    if (g_pac_scratch) tfree(g_pac_scratch);
    if (g_pac_first_y) tfree(g_pac_first_y);
    if (g_pac_Qprev) tfree(g_pac_Qprev);
    if (g_pac_Qnew) tfree(g_pac_Qnew);
    g_pac_Qprev = g_pac_Qnew = NULL;
    { int jj; for (jj = 2; jj <= PAC_GEAR_MAXORD; jj++) { if (g_pac_Qh[jj]) tfree(g_pac_Qh[jj]); g_pac_Qh[jj] = NULL; } }
    g_pac_gear = 0;
    g_pac_gen = 0;
    g_pac_Cprev = g_pac_Cnew = NULL;
    g_pac_Tprev = g_pac_Tnew = NULL;
    g_pac_omega = NULL; g_pac_useforce = NULL; g_pac_demod = NULL;
    g_pac_rhs_re = g_pac_rhs_im = g_pac_scratch = NULL;
    g_pac_first_y = NULL;
    g_pac_first_t = 0.0;
    g_pac_have_first = 0;
    g_pac_map = NULL;
    g_pac_n = g_pac_neq = g_pac_M = 0;
    g_pac_failed = 1;
    g_pac_demod_active = 0;
}

/* br/bi: shared forcing pattern (size neq+1, 1-based), or NULL for no
 * forcing at all this pass. omega[M]/useforce[M]: per-column angular
 * frequency and whether the shared forcing applies to that column.
 * ic0[n*M]: initial (t=0) complex value at each reduced-state row d,
 * column m (row-major, index d*M+m). demod_active: accumulate sideband-0
 * quadrature at every row for every column (used for the "final" pass
 * only -- expensive-ish, skipped for the Phi/z pass). T: period, needed
 * for the demod 1/T normalisation. */
static int
pac_reset(CKTcircuit *ckt, PACstateMap *map, int M, const double *br, const double *bi,
          const double *omega, const int *useforce, const PACc *ic0, int demod_active, double T)
{
    int n = PAC_NST(map), neq, d, m;

    pac_free();

    neq = SMPmatSize(ckt->CKTmatrix);
    g_pac_map = map;
    g_pac_n = n;
    g_pac_neq = neq;
    g_pac_M = M;
    g_pac_T = T;
    g_pac_demod_active = demod_active;
    g_pac_tprev = 0.0;

    g_hid_nrow = neq + (map->generic ? g_hid_nh : 0);   /* card #557: virtual rows of the hidden states */
    g_pac_Cprev = TMALLOC(PACc, (size_t) (g_hid_nrow + 1) * M);
    g_pac_Cnew  = TMALLOC(PACc, (size_t) (g_hid_nrow + 1) * M);
    if (g_pac_Cnew) memset(g_pac_Cnew, 0, (size_t) (g_hid_nrow + 1) * M * sizeof(PACc));
    g_pac_gen = map->generic;
    {
        size_t tsz = g_pac_gen ? (size_t) (neq + 1) : (size_t) n;
        g_pac_Tprev = TMALLOC(PACc, tsz * M);
        g_pac_Tnew  = TMALLOC(PACc, tsz * M);
        if (g_pac_gen) {
            g_pac_Qprev = TMALLOC(PACc, (size_t) (neq + 1) * M);
            g_pac_Qnew  = TMALLOC(PACc, (size_t) (neq + 1) * M);
            if (g_pac_Qprev) memset(g_pac_Qprev, 0, (size_t) (neq + 1) * M * sizeof(PACc));
            if (g_pac_Qnew)  memset(g_pac_Qnew, 0, (size_t) (neq + 1) * M * sizeof(PACc));
            if (g_pac_Tprev) memset(g_pac_Tprev, 0, tsz * M * sizeof(PACc));
            if (g_pac_Tnew)  memset(g_pac_Tnew, 0, tsz * M * sizeof(PACc));
            g_pac_gear = (ckt->CKTintegrateMethod == GEAR);   /* card #549 */
            if (g_pac_gear) {
                int jj;
                for (jj = 2; jj <= PAC_GEAR_MAXORD; jj++) {
                    g_pac_Qh[jj] = TMALLOC(PACc, (size_t) (neq + 1) * M);
                    if (g_pac_Qh[jj]) memset(g_pac_Qh[jj], 0, (size_t) (neq + 1) * M * sizeof(PACc));
                }
            }
        }
    }
    g_pac_omega = TMALLOC(double, M);
    g_pac_useforce = TMALLOC(int, M);
    g_pac_rhs_re = TMALLOC(double, (size_t) (neq + 1));
    g_pac_rhs_im = TMALLOC(double, (size_t) (neq + 1));
    g_pac_scratch = TMALLOC(double, (size_t) (neq + 1));
    g_pac_br = (double *) br;   /* not owned, caller keeps alive for this pass */
    g_pac_bi = (double *) bi;
    if (demod_active) {
        g_pac_demod = TMALLOC(PACc, (size_t) neq * M);
        g_pac_first_y = TMALLOC(PACc, (size_t) (neq + 1) * M);
    }

    if (g_pac_gen && (!g_pac_Qprev || !g_pac_Qnew)) {
        pac_free();
        return 1;
    }
    if (g_pac_gear) {
        int jj;
        for (jj = 2; jj <= PAC_GEAR_MAXORD; jj++)
            if (!g_pac_Qh[jj]) { pac_free(); return 1; }
    }
    if (!g_pac_Cprev || !g_pac_Cnew || !g_pac_Tprev || !g_pac_Tnew || !g_pac_omega ||
        !g_pac_useforce || !g_pac_rhs_re || !g_pac_rhs_im || !g_pac_scratch ||
        (demod_active && (!g_pac_demod || !g_pac_first_y))) {
        pac_free();
        return 1;
    }

    memset(g_pac_Cprev, 0, (size_t) (g_hid_nrow + 1) * M * sizeof(PACc));
    memset(g_pac_Tprev, 0, (size_t) (g_pac_gen ? (neq + 1) : n) * M * sizeof(PACc));
    if (demod_active) {
        memset(g_pac_demod, 0, (size_t) neq * M * sizeof(PACc));
        memset(g_pac_first_y, 0, (size_t) (neq + 1) * M * sizeof(PACc));
    }
    g_pac_xls_have = 0;
    g_pac_have_first = 0;
    g_pac_first_t = 0.0;
    for (m = 0; m < M; m++) {
        g_pac_omega[m] = omega[m];
        g_pac_useforce[m] = useforce[m];
    }
    for (d = 0; d < n; d++) {
        int row = pac_map_row(map, d);
        for (m = 0; m < M; m++)
            pac_set(g_pac_Cprev, row, M, m, ic0[(size_t) d * M + m]);
    }
    g_pac_exact_jac = getenv("PAC_EXACT_JACOBIAN") != NULL;
    g_pac_gen_bad = 0;
    g_pac_failed = 0;
    PSSsaltReset();   /* card #494 third pass: fresh run -- filter the first step's event */
    return 0;
}

/* card #465 THIRD pass, task 1 (accuracy): .pstb vs native .stb on a
 * LINEAR circuit (pstb_linear.cir) was found to converge only O(1/Points),
 * not the O(h^2) a trapezoidal per-step recursion should give -- traced to
 * this file's own documented design trade-off (header, "PAC (`.pac ...
 * method=shoot solver=td`)" paragraph and "refactor" note): PSSpacRecordStep()
 * below solves against ckt->CKTmatrix AS LEFT by NIiter() (SMPsolve() only,
 * no refactor), which is G(x_{k-1})+ag0*C(x_{k-1}), NOT G(x_k)+ag0*C(x_k).
 *
 * Root cause, traced in src/maths/ni/niiter.c: on the converged
 * (MODEINITFLOAT, CKTnoncon==0) return path, NIiter() returns OK
 * immediately (line ~435) WITHOUT reaching its own
 * `SWAP(double*, ckt->CKTrhs, ckt->CKTrhsOld)` (line ~463, only on the
 * loop's continue path). So at this hook's call site (dctran.c, right
 * after CKTaccept()): ckt->CKTrhs holds x_k (the just-solved, TRUE
 * accepted point) but ckt->CKTrhsOld -- the array every device loader
 * reads node voltages from, e.g. capload.c's `*(ckt->CKTrhsOld+...)`,
 * resload.c likewise -- still holds x_{k-1}, the second-to-last Newton
 * iterate. ckt->CKTmatrix was last loaded+factored AT x_{k-1} (the top of
 * that final NIiter() loop pass), never refactored at x_k itself. This is
 * the SAME staleness pssshoot.c's PSSsensRecordStep() documents for the
 * PSS shooting-Newton monodromy and deliberately does NOT fix there (that
 * path pays this cost on every accepted step of EVERY period of EVERY
 * Newton iteration -- found too slow on R3, and a naive extra CKTload()
 * was found to measurably perturb R2's SW switch model's history-
 * dependent internal state -- see PSSsensRecordStep()'s own long comment).
 *
 * PAC's OWN recording pass is different in exactly the way that trade-off
 * turns on: pac_engine_run() calls pac_run_period() at most TWICE total
 * per `.pac`/`.pstb` invocation (pass 1: Phi + forced z; pass 2: re-
 * propagate + demodulate) -- PSS's Newton has already converged x0/T
 * before this file's code runs at all, so "refactor at every accepted
 * point" here costs at most 2x this file's own (already O(one period)
 * per the header) wall time, not O(Newton iterations x periods). The
 * card's plan calls this "affordable" -- validated below (accuracy
 * table, card #465 THIRD pass report).
 *
 * Guarded behind PAC_EXACT_JACOBIAN (env var, default unset = OFF, same
 * opt-in convention as pssshoot.c's PSS_SHOOT_JACOBIAN=analytic) rather
 * than made the unconditional default: validated to fix the O(1/Points)
 * convergence to true O(h^2) on the pure-LTI circuits this task's target
 * is stated against (pstb_linear.cir, pac_lti.cir -- no history-dependent
 * device state at all) but the save/restore in pac_refactor_at_accepted()
 * below (CKTrhs, CKTrhsOld, CKTstate0 -- everything an extra CKTload()
 * can touch beyond what this hook already owns) has only been exercised
 * end to end against R2/R3/R4's actual switching devices via this
 * project's existing pac_sb0(b)/pstb_loop regression checks, not proven
 * safe by the SAME kind of targeted, single-device validation
 * PSSsensRecordStep()'s own comment describes doing for the SW model --
 * see the card #465 THIRD pass report for exactly what was and was not
 * checked. Default OFF until that gap is closed. */
/* card #541: copy the currently loaded (unfactored) matrix into the last trajectory entry (native path). */
static void
pac_traj_store_J(CKTcircuit *ckt)
{
    PACtrajStep *st = &g_pac_traj[g_pac_traj_n - 1];
    int nnz = SMPnumElements(ckt->CKTmatrix);
    const char *mx = getenv("PAC_TRAJ_MAX_MB");
    double cap = (mx && atof(mx) > 0.0) ? atof(mx) * 1048576.0 : 3.0e9;
    if (g_pac_traj_bytes + (double) nnz * 16.0 > cap) {
        fprintf(stderr, "PAC-TD: adjoint trajectory would exceed %.0f MB (J_k of every step is kept for the "
                "exact transpose; set PAC_TRAJ_MAX_MB higher or use fewer pacpoints).\n", cap / 1048576.0);
        g_pac_traj_failed = 1;
        return;
    }
    st->Jri = TMALLOC(int, nnz + 1); st->Jci = TMALLOC(int, nnz + 1); st->Jv = TMALLOC(double, nnz + 1);
    if (!st->Jri || !st->Jci || !st->Jv) { g_pac_traj_failed = 1; return; }
    st->nJ = SMPcopyElements(ckt->CKTmatrix, st->Jri, st->Jci, st->Jv);
    g_pac_traj_bytes += (double) nnz * 16.0;
}

static int
pac_refactor_at_accepted(CKTcircuit *ckt)
{
    int neq = SMPmatSize(ckt->CKTmatrix);
    int nstates = ckt->CKTnumStates;
    double *save_rhs, *save_rhsold, *save_state0 = NULL;
    int err;
    long save_mode;
    int save_bypass, save_noncon;

    save_rhs = TMALLOC(double, neq + 1);
    save_rhsold = TMALLOC(double, neq + 1);
    if (nstates > 0)
        save_state0 = TMALLOC(double, nstates);
    if (!save_rhs || !save_rhsold || (nstates > 0 && !save_state0)) {
        if (save_rhs) tfree(save_rhs);
        if (save_rhsold) tfree(save_rhsold);
        if (save_state0) tfree(save_state0);
        return 1;
    }

    memcpy(save_rhs, ckt->CKTrhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(save_rhsold, ckt->CKTrhsOld, (size_t) (neq + 1) * sizeof(double));
    if (nstates > 0)
        memcpy(save_state0, ckt->CKTstate0, (size_t) nstates * sizeof(double));

    /* devices read node voltages from CKTrhsOld (see capload.c/resload.c) --
     * put x_k (the just-accepted CKTrhs) there so CKTload() evaluates
     * G(x_k)/C(x_k), not G(x_{k-1})/C(x_{k-1}). */
    memcpy(ckt->CKTrhsOld, ckt->CKTrhs, (size_t) (neq + 1) * sizeof(double));

    /* card #541: at this hook CKTmode is MODETRAN|MODEINITPRED (dctran has already switched for the NEXT
     * step), under which MOS1 and friends evaluate at the PREDICTED voltages extrapolated from CKTstate1/2,
     * not at CKTrhsOld = x_k: the "exact" J_k was J at a predicted point (0.4462 vs 0.4677 on the #539 bug 3
     * chain, 5 % on a threshold crossing). Load in MODEINITFLOAT, bypass off, like pac_gen_loads(). */
    save_mode = ckt->CKTmode; save_bypass = ckt->CKTbypass; save_noncon = ckt->CKTnoncon;
    ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODETRAN | MODEINITFLOAT;
    ckt->CKTbypass = 0;
    err = CKTload(ckt);
    ckt->CKTmode = save_mode;
    /* card #541: native adjoint -- keep the exact J_k the forward pass used (loaded at x_k from the
     * accepted device state) for the backward sweep, instead of reloading it in reverse order from device
     * state (limiting / previous-iterate history) that no longer belongs to step k. Without this the
     * adjoint of a MOS1 comparator chain was ~0 (#539 bug 3). */
    if (!err && g_pac_traj_active && g_pac_traj_n > 0 && !g_pac_traj_failed)
        pac_traj_store_J(ckt);
    if (!err)
        err = SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin);
    ckt->CKTbypass = save_bypass; ckt->CKTnoncon = save_noncon;

    /* restore CKTrhs/CKTrhsOld/CKTstate0 exactly as dctran.c's own
     * subsequent code (CKTdump right after this hook returns, and the
     * NEXT accepted step's own CKTload()) expects to find them. Only
     * ckt->CKTmatrix is left changed -- holding the fresh
     * G(x_k)+ag0*C(x_k) factorization PSSpacRecordStep()'s own
     * SMPsolve() calls (right after this function returns) need. */
    memcpy(ckt->CKTrhs, save_rhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
    if (nstates > 0)
        memcpy(ckt->CKTstate0, save_state0, (size_t) nstates * sizeof(double));

    tfree(save_rhs);
    tfree(save_rhsold);
    if (save_state0) tfree(save_state0);
    return err ? 1 : 0;
}


/* ------------------------------------------------------------------ */
/* card #542: generic charge matrix extraction and the generic recursion
 * (see the design note by g_pac_gen above). */

static int
pac_trip_fit(PACtrip *t, int n)
{
    if (n > t->cap) {
        int c = n + n / 4 + 64;
        if (t->ri) tfree(t->ri);
        if (t->ci) tfree(t->ci);
        if (t->v) tfree(t->v);
        t->ri = TMALLOC(int, c);
        t->ci = TMALLOC(int, c);
        t->v = TMALLOC(double, c);
        t->cap = c;
        if (!t->ri || !t->ci || !t->v) { t->cap = 0; t->n = 0; return 1; }
    }
    t->n = 0;
    return 0;
}

static void
pac_trip_free(PACtrip *t)
{
    if (t->ri) tfree(t->ri);
    if (t->ci) tfree(t->ci);
    if (t->v) tfree(t->v);
    t->ri = t->ci = NULL; t->v = NULL; t->n = t->cap = 0;
}

/* One load of the Jacobian with CKTag[0]=a, copied to t (every stored element). */
static int
pac_jload(CKTcircuit *ckt, double a, PACtrip *t)
{
    double sag0 = ckt->CKTag[0];
    int err, nnz;

    ckt->CKTag[0] = a;
#ifdef XSPICE
    g_xs_ag0_fix = sag0;   /* card #557: code-model integrals keep the true ag0 while the device charges are separated */
#endif
    err = CKTload(ckt);
#ifdef XSPICE
    g_xs_ag0_fix = 0.0;
#endif
    ckt->CKTag[0] = sag0;
    if (err)
        return 1;
    nnz = SMPnumElements(ckt->CKTmatrix);
    if (pac_trip_fit(t, nnz))
        return 1;
    t->n = SMPcopyElements(ckt->CKTmatrix, t->ri, t->ci, t->v);
    return 0;
}

/* Affinity test of one set of three loads (ag0/2, ag0, 2*ag0): J3 = 3 J1 - 2 J2. */
static int
pac_affine_ok(const PACtrip *J1, const PACtrip *J2, const PACtrip *J3, const char *tname)
{
    int i;
    /* test hook (red proof of the honesty guard): PAC_GEN_TEST_NONAFFINE=<device type name> makes that
     * type's load look non-affine in ag0 (the whole-circuit check fails if any is set) */
    {
        const char *tn = getenv("PAC_GEN_TEST_NONAFFINE");
        if (tn && (tname == NULL || !strcmp(tn, tname)))
            return 0;
    }
    if (J1->n != J2->n || J1->n != J3->n)
        return 0;
    for (i = 0; i < J1->n; i++) {
        double a = J1->v[i], b = J2->v[i], c = J3->v[i];
        double pred = 3.0 * a - 2.0 * b;
        double scale = fabs(a) + fabs(b) + fabs(c);
        if (fabs(c - pred) > 1e-7 * scale)
            return 0;
    }
    return 1;
}

/* A load that is not affine in ag0 was seen at ag0: find the device TYPE whose own load is the
 * culprit (loads each type alone into a cleared matrix).  Caller has saved ckt state. */
static void
pac_gen_identify_bad(CKTcircuit *ckt, double ag0)
{
    int i, neq = SMPmatSize(ckt->CKTmatrix), r;
    PACtrip A = { 0, 0, NULL, NULL, NULL }, B = A, C3 = A;

    g_pac_gen_badname[0] = '\0';
    for (i = 0; i < DEVmaxnum && !g_pac_gen_badname[0]; i++) {
        PACtrip *tt[3]; double ag[3]; int k, err = 0;
        if (!DEVices[i] || !DEVices[i]->DEVload || !ckt->CKThead[i])
            continue;
        tt[0] = &B; tt[1] = &A; tt[2] = &C3;
        ag[0] = 0.5 * ag0; ag[1] = ag0; ag[2] = 2.0 * ag0;
        for (k = 0; k < 3 && !err; k++) {
            double sag0 = ckt->CKTag[0];
            int nnz;
            SMPclear(ckt->CKTmatrix);
            for (r = 0; r <= neq; r++) ckt->CKTrhs[r] = 0.0;
            ckt->CKTag[0] = ag[k];
            err = DEVices[i]->DEVload(ckt->CKThead[i], ckt);
            ckt->CKTag[0] = sag0;
            if (err) break;
            nnz = SMPnumElements(ckt->CKTmatrix);
            if (pac_trip_fit(tt[k], nnz)) { err = 1; break; }
            tt[k]->n = SMPcopyElements(ckt->CKTmatrix, tt[k]->ri, tt[k]->ci, tt[k]->v);
        }
        if (!err && !pac_affine_ok(&A, &B, &C3, DEVices[i]->DEVpublic.name ? DEVices[i]->DEVpublic.name : "?"))
            snprintf(g_pac_gen_badname, sizeof g_pac_gen_badname, "%s",
                     DEVices[i]->DEVpublic.name ? DEVices[i]->DEVpublic.name : "?");
    }
    pac_trip_free(&A); pac_trip_free(&B); pac_trip_free(&C3);
    if (!g_pac_gen_badname[0])
        snprintf(g_pac_gen_badname, sizeof g_pac_gen_badname, "(not isolated to one device type)");
}

/* ---------------------------------------------------------------------------------------------------------------------------
 * card #557: hidden-state machinery (see the comment at g_hid_nc).  XSPICE only. */
#ifdef XSPICE
static int g_hid_nhy = 0;                 /* slots of hyst models (Boolean state; flips are events) */
static int *g_hid_hyslot = NULL;

static void
pac_hid_free(void)
{
    if (g_hid_cslot) tfree(g_hid_cslot);
    if (g_hid_inrow) tfree(g_hid_inrow);
    if (g_hid_act) tfree(g_hid_act);
    if (g_hid_cscale) tfree(g_hid_cscale);
    if (g_hid_cslew) tfree(g_hid_cslew);
    g_hid_cslew = NULL;
    if (g_hid_patP) tfree(g_hid_patP);
    if (g_hid_patR) tfree(g_hid_patR);
    if (g_hid_patB) tfree(g_hid_patB);
    if (g_hid_Bh) tfree(g_hid_Bh);
    if (g_hid_P) tfree(g_hid_P);
    if (g_hid_R) tfree(g_hid_R);
    if (g_hid_Jc) tfree(g_hid_Jc);
    if (g_hid_hyslot) tfree(g_hid_hyslot);
    g_hid_cslot = g_hid_inrow = g_hid_act = g_hid_hyslot = NULL;
    g_hid_cscale = g_hid_Bh = g_hid_P = g_hid_R = g_hid_Jc = NULL;
    g_hid_patP = g_hid_patR = g_hid_patB = NULL;
    g_hid_nc = g_hid_nin = g_hid_nh = g_hid_mode = g_hid_nrow = g_hid_flips = g_hid_nhy = 0;
}

static void
pac_hid_add_unique(int **list, int *n, int v)
{
    int i;
    for (i = 0; i < *n; i++)
        if ((*list)[i] == v) return;
    *list = TREALLOC(int, *list, *n + 1);
    (*list)[(*n)++] = v;
}

/* Candidate slots and input rows of every instance of a td-supported stateful analog code model.  Returns the number of candidate
 * slots (also g_hid_nc). */
static int
pac_hid_collect(CKTcircuit *ckt)
{
    int i, k, d, ci, j;
    pac_hid_free();
    if (getenv("PAC_HID_OFF"))
        return 0;
    for (i = 0; i < DEVmaxnum; i++) {
        GENmodel *m;
        const char *nm;
        int is_hyst;
        if (!DEVices[i] || !DEVices[i]->DEVpublic.cm_func || !ckt->CKThead[i] || !DEVices[i]->DEVpublic.name)
            continue;
        nm = DEVices[i]->DEVpublic.name;
        if (!PSSxsTdStateful(nm))
            continue;
        is_hyst = !strcmp(nm, "hyst");
        for (m = ckt->CKThead[i]; m; m = m->GENnextModel) {
            GENinstance *g;
            for (g = m->GENinstances; g; g = g->GENnextInstance) {
                MIFinstance *here = (MIFinstance *) g;
                if (!here->analog) continue;
                for (k = 0; k < here->num_state; k++)
                    for (d = 0; d < here->state[k].doubles; d++) {
                        int idx = here->state[k].index + d;
                        if (is_hyst) pac_hid_add_unique(&g_hid_hyslot, &g_hid_nhy, idx);
                        else {
                            int n0 = g_hid_nc;
                            pac_hid_add_unique(&g_hid_cslot, &g_hid_nc, idx);
                            if (g_hid_nc > n0) {
                                g_hid_cslew = TREALLOC(char, g_hid_cslew, g_hid_nc);
                                g_hid_cslew[g_hid_nc - 1] = !strcmp(nm, "slew");
                            }
                        }
                    }
                if (is_hyst) continue;   /* a discrete-state static nonlinearity: J exact, no hidden states */
                for (ci = 0; ci < here->num_conn; ci++) {
                    if (here->conn[ci]->is_null || !here->conn[ci]->is_input) continue;
                    for (j = 0; j < here->conn[ci]->size; j++) {
                        Mif_Port_Data_t *fast = here->conn[ci]->port[j];
                        if (fast->is_null) continue;
                        switch (fast->type) {
                        case MIF_VOLTAGE: case MIF_DIFF_VOLTAGE: case MIF_CONDUCTANCE: case MIF_DIFF_CONDUCTANCE:
                            if (fast->smp_data.pos_node > 0) pac_hid_add_unique(&g_hid_inrow, &g_hid_nin, fast->smp_data.pos_node);
                            if (fast->smp_data.neg_node > 0) pac_hid_add_unique(&g_hid_inrow, &g_hid_nin, fast->smp_data.neg_node);
                            break;
                        case MIF_CURRENT: case MIF_DIFF_CURRENT: case MIF_VSOURCE_CURRENT: case MIF_RESISTANCE: case MIF_DIFF_RESISTANCE:
                            if (fast->smp_data.ibranch > 0) pac_hid_add_unique(&g_hid_inrow, &g_hid_nin, fast->smp_data.ibranch);
                            break;
                        default:
                            break;
                        }
                    }
                }
            }
        }
    }
    if (g_hid_nc > 0) {
        g_hid_cscale = TMALLOC(double, g_hid_nc);
        g_hid_patP = TMALLOC(unsigned char, g_hid_nc);
        g_hid_patB = TMALLOC(unsigned char, g_hid_nc);
        g_hid_patR = TMALLOC(unsigned char, (size_t) g_hid_nc * g_hid_nc);
        memset(g_hid_cscale, 0, (size_t) g_hid_nc * sizeof(double));
    }
    return g_hid_nc;
}

/* One full load at (x, current states); F = J x - rhs over rows 1..neq, h0 = state0 of the candidate (or active) slots. */
static int
pac_hid_eval(CKTcircuit *ckt, double ag0, const double *x, double *F, double *h0, const int *slots, int ns)
{
    int neq = SMPmatSize(ckt->CKTmatrix), i, r;
    memcpy(ckt->CKTrhsOld, x, (size_t) (neq + 1) * sizeof(double));
    ckt->CKTnoncon = 0;
    if (pac_jload(ckt, ag0, &g_hid_J))
        return 1;
    for (r = 1; r <= neq; r++)
        F[r] = -ckt->CKTrhs[r];
    for (i = 0; i < g_hid_J.n; i++) {
        int rr = g_hid_J.ri[i], cc = g_hid_J.ci[i];
        if (rr > 0 && cc > 0) F[rr] += g_hid_J.v[i] * x[cc];
    }
    for (i = 0; i < ns; i++)
        h0[i] = ckt->CKTstate0[slots[i]];
    return 0;
}

/* Smallest relative response, over the entries that responded at all, of a central-difference pair: the step is enlarged until the
 * response is resolved (>= 1e-7 of its own entry) or three times.  Entries that did not respond (bit-identical) are no dependence. */
static double
pac_hid_resp(const double *Fp, const double *Fm, int n0, int n1)
{
    int r; double q = 0.0;
    for (r = n0; r <= n1; r++) {
        double d = fabs(Fp[r] - Fm[r]);
        if (d == 0.0) continue;
        d /= (fabs(Fp[r]) + fabs(Fm[r]) + 1e-300);
        if (d > q) q = d;
    }
    return q;
}

/* One central-difference probe: kind 0 = node row `idx` of x, kind 1 = CKTstate1[idx].  On return Fp/Fm/hp/hm hold the last pair and
 * *delta the step used. */
static int
pac_hid_probe(CKTcircuit *ckt, double ag0, double *x, int kind, int idx, double del0, const int *slots, int ns,
              double *Fp, double *Fm, double *hp, double *hm, double *delta)
{
    int neq = SMPmatSize(ckt->CKTmatrix), tries, err = 0;
    double base = (kind == 0) ? x[idx] : ckt->CKTstate1[idx];
    double del = del0;
    for (tries = 0; tries < 4; tries++) {
        double qF, qH;
        if (kind == 0) x[idx] = base + del; else ckt->CKTstate1[idx] = base + del;
        err = pac_hid_eval(ckt, ag0, x, Fp, hp, slots, ns);
        if (!err) {
            if (kind == 0) x[idx] = base - del; else ckt->CKTstate1[idx] = base - del;
            err = pac_hid_eval(ckt, ag0, x, Fm, hm, slots, ns);
        }
        if (kind == 0) x[idx] = base; else ckt->CKTstate1[idx] = base;
        if (err) return 1;
        qF = pac_hid_resp(Fp, Fm, 1, neq);
        qH = pac_hid_resp(hp - 1, hm - 1, 1, ns);
        if (tries == 3 || ((qF == 0.0 || qF >= 1e-7) && (qH == 0.0 || qH >= 1e-7)))
            break;
        { double q = (qF > 0.0 && qF < 1e-7) ? qF : 1e300; if (qH > 0.0 && qH < 1e-7 && qH < q) q = qH;
          del *= fmin(1e3, 1e-6 / q); }
    }
    *delta = del;
    return 0;
}

/* Central differences of the model.  mode 1 (prepass): pattern of dependences over all candidate slots.  mode 2: P, R, B and the
 * input columns of J for the ACTIVE slots.  Called from pac_gen_loads between its J2 and J1 loads; xbase = x_k, st0 = CKTstate0 of
 * the accepted step.  Returns 1 on a hard failure. */
static int
pac_hid_stage(CKTcircuit *ckt, double ag0, const double *xbase, const double *st0)
{
    int neq = SMPmatSize(ckt->CKTmatrix), nstates = ckt->CKTnumStates, ns, i, a, r, err = 0;
    const int *slots;
    double *Fp = NULL, *Fm = NULL, *hp = NULL, *hm = NULL, *x = NULL, *save_s0 = NULL, *save_s1 = NULL;
    int act = (g_hid_mode == 2);
    int nact = act ? g_hid_nh : g_hid_nc;
    int *aslot = NULL;
    double S = 0.0;

    if (g_hid_nc <= 0 || nstates <= 0)
        return 0;
    ns = nact;
    aslot = TMALLOC(int, ns > 0 ? ns : 1);
    for (a = 0; a < ns; a++)
        aslot[a] = act ? g_hid_cslot[g_hid_act[a]] : g_hid_cslot[a];
    slots = aslot;
    Fp = TMALLOC(double, neq + 1); Fm = TMALLOC(double, neq + 1);
    hp = TMALLOC(double, ns > 0 ? ns + 1 : 2); hm = TMALLOC(double, ns > 0 ? ns + 1 : 2);
    x = TMALLOC(double, neq + 1);
    save_s0 = TMALLOC(double, nstates); save_s1 = TMALLOC(double, nstates);
    if (!Fp || !Fm || !hp || !hm || !x || !save_s0 || !save_s1) { err = 1; goto out; }
    memcpy(x, xbase, (size_t) (neq + 1) * sizeof(double));
    memcpy(save_s0, ckt->CKTstate0, (size_t) nstates * sizeof(double));
    memcpy(save_s1, ckt->CKTstate1, (size_t) nstates * sizeof(double));

    if (!act)
        for (a = 0; a < ns; a++) {
            double v = fabs(st0[slots[a]]);
            if (v > g_hid_cscale[a]) g_hid_cscale[a] = v;
        }
    for (a = 0; a < g_hid_nc; a++)
        if (g_hid_cscale[a] > S) S = g_hid_cscale[a];
    /* inputs */
    for (i = 0; i < g_hid_nin && !err; i++) {
        int row = g_hid_inrow[i];
        double del = 0.0;
        err = pac_hid_probe(ckt, ag0, x, 0, row, 1e-6 * fmax(fabs(xbase[row]), 1e-3), slots, ns, Fp, Fm, hp, hm, &del);
        if (err) break;
        if (act) {
            for (r = 1; r <= neq; r++) g_hid_Jc[(size_t) r * g_hid_nin + i] = (Fp[r] - Fm[r]) / (2.0 * del);
            for (a = 0; a < ns; a++) g_hid_P[(size_t) a * g_hid_nin + i] = (hp[a] - hm[a]) / (2.0 * del);
        } else {
            for (a = 0; a < ns; a++)
                if (hp[a] != hm[a]) g_hid_patP[a] = 1;
        }
    }
    /* state slots */
    for (a = 0; a < ns && !err; a++) {
        int sl = slots[a], b;
        double h1 = save_s1[sl], del = 0.0;
        err = pac_hid_probe(ckt, ag0, x, 1, sl, 1e-6 * fmax(fmax(S, fabs(h1)), 1e-9), slots, ns, Fp, Fm, hp, hm, &del);
        if (err) break;
        if (act) {
            for (r = 1; r <= neq; r++) g_hid_Bh[(size_t) r * ns + a] = (Fp[r] - Fm[r]) / (2.0 * del);
            for (b = 0; b < ns; b++) g_hid_R[(size_t) b * ns + a] = (hp[b] - hm[b]) / (2.0 * del);
        } else {
            for (r = 1; r <= neq; r++)
                if (Fp[r] != Fm[r]) { g_hid_patB[a] = 1; break; }
            for (b = 0; b < ns; b++)
                if (hp[b] != hm[b]) g_hid_patR[(size_t) b * ns + a] = 1;
        }
    }
    memcpy(ckt->CKTstate0, save_s0, (size_t) nstates * sizeof(double));
    memcpy(ckt->CKTstate1, save_s1, (size_t) nstates * sizeof(double));
out:
    if (Fp) tfree(Fp); if (Fm) tfree(Fm); if (hp) tfree(hp); if (hm) tfree(hm); if (x) tfree(x);
    if (save_s0) tfree(save_s0); if (save_s1) tfree(save_s1); if (aslot) tfree(aslot);
    return err ? 1 : 0;
}

/* End of the prepass: reachable-from-inputs AND needed-by-a-later-step slots become the active set. */
static void
pac_hid_finalize(void)
{
    int nc = g_hid_nc, j, l, changed, na = 0;
    unsigned char *reach = NULL, *need = NULL;
    if (nc <= 0)
        return;
    reach = TMALLOC(unsigned char, nc); need = TMALLOC(unsigned char, nc);
    for (j = 0; j < nc; j++) { reach[j] = g_hid_patP[j]; need[j] = g_hid_patB[j]; }
    do {
        changed = 0;
        for (j = 0; j < nc; j++)
            for (l = 0; l < nc; l++) {
                if (j == l || !g_hid_patR[(size_t) j * nc + l]) continue;
                if (reach[l] && !reach[j]) { reach[j] = 1; changed = 1; }   /* j (new) depends on l */
                if (need[j] && !need[l]) { need[l] = 1; changed = 1; }       /* l feeds a needed slot */
            }
    } while (changed);
    g_hid_act = TMALLOC(int, nc);
    for (j = 0; j < nc; j++)
        if (reach[j] && need[j]) g_hid_act[na++] = j;
    g_hid_nh = na;
    g_hid_mode = 0;
    g_hid_slewact = 0;
    for (j = 0; j < na; j++)
        if (g_hid_cslew && g_hid_cslew[g_hid_act[j]]) g_hid_slewact = 1;
    if (getenv("PAC_HID_TRACE")) {
        int q;
        for (j = 0; j < nc; j++) {
            fprintf(stderr, "PAC-HID: cand %d state %d scale %.4g patP=%d patB=%d reach=%d need=%d R:", j, g_hid_cslot[j], g_hid_cscale[j], g_hid_patP[j], g_hid_patB[j], reach[j], need[j]);
            for (q = 0; q < nc; q++) fprintf(stderr, "%d", g_hid_patR[(size_t) j * nc + q]);
            fprintf(stderr, "\n");
        }
    }
    tfree(reach); tfree(need);
    if (getenv("PAC_HID_TRACE"))
        fprintf(stderr, "PAC-TD: hidden-state set: %d candidate slots -> %d active, %d input rows\n", nc, na, g_hid_nin);
}
#endif   /* XSPICE */

/* Load J(ag0) and C = dq/dx at the TRUE accepted point x_k = CKTrhs.  Leaves ckt->CKTmatrix holding
 * the factored J_k (if do_factor) -- exactly pac_refactor_at_accepted()'s contract -- and the
 * charge matrix in g_pac_C.  Everything the loads touch (CKTrhs, CKTrhsOld, CKTstate0, noncon,
 * mode, bypass) is saved and restored.  Returns 1 on a hard failure; a non-affine load sets
 * g_pac_gen_bad (+ g_pac_gen_badname). */
static int
pac_gen_loads(CKTcircuit *ckt, double a1, int check, int do_factor)
{
    int neq = SMPmatSize(ckt->CKTmatrix), nstates = ckt->CKTnumStates, i, err = 0, n;
    double *save_rhs, *save_rhsold, *save_state0 = NULL;
    long save_mode = ckt->CKTmode;
    int save_noncon = ckt->CKTnoncon, save_bypass = ckt->CKTbypass;
    PACtrip *C = &g_pac_C;
    double d12 = 0.5 * a1;
    int save_order = ckt->CKTorder;

    save_rhs = TMALLOC(double, neq + 1);
    save_rhsold = TMALLOC(double, neq + 1);
    if (nstates > 0)
        save_state0 = TMALLOC(double, nstates);
    if (!save_rhs || !save_rhsold || (nstates > 0 && !save_state0)) {
        if (save_rhs) tfree(save_rhs);
        if (save_rhsold) tfree(save_rhsold);
        if (save_state0) tfree(save_state0);
        return 1;
    }
    memcpy(save_rhs, ckt->CKTrhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(save_rhsold, ckt->CKTrhsOld, (size_t) (neq + 1) * sizeof(double));
    if (nstates > 0)
        memcpy(save_state0, ckt->CKTstate0, (size_t) nstates * sizeof(double));

    memcpy(ckt->CKTrhsOld, ckt->CKTrhs, (size_t) (neq + 1) * sizeof(double));
    ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODETRAN | MODEINITFLOAT;
    ckt->CKTbypass = 0;
#ifdef XSPICE
    /* card #557: cm_analog_integrate() switches on ckt->CKTorder, which is already the NEXT step's at this call site (see
     * PSSpacRecordStep); use the order of the step just accepted, recovered from ag as there. */
    if (g_hid_nc > 0 && ckt->CKTintegrateMethod == TRAPEZOIDAL)
        ckt->CKTorder = (fabs(ckt->CKTag[1] + ckt->CKTag[0]) < 1e-6 * fabs(ckt->CKTag[0])) ? 1 : 2;
#endif

    ckt->CKTnoncon = 0;
    if (g_pac_gen_warm) {
        /* backward sweep: CKTstate0 is stale (end of the recorded period). Device limiting reads it, so
         * iterate the load at fixed x_k until no device limits: state0 then equals the operating point
         * at x_k. (Forward, state0 IS the converged iterate and this is skipped.) */
        int it;
        for (it = 0; it < 80 && !err; it++) {
            ckt->CKTnoncon = 0;
            if (CKTload(ckt)) { err = 1; break; }
            if (!ckt->CKTnoncon) break;
        }
        ckt->CKTnoncon = 0;
    }
    if (!err)
        err = pac_jload(ckt, 0.5 * a1, &g_pac_J2);
    if (!err && ckt->CKTnoncon) {
        /* device limiting / convergence checks read CKTstate0 (the previous iterate's operating
         * point). Forward, state0 is the converged iterate and this never triggers; in the backward
         * sweep state0 is stale, so iterate the load at fixed x_k until no device limits (state0
         * then equals the operating point at x_k) -- otherwise the loaded G and C would be those of a
         * LIMITED (wrong) bias. */
        int it;
        for (it = 0; it < 60; it++) {
            ckt->CKTnoncon = 0;
            if (CKTload(ckt)) { err = 1; break; }
            if (!ckt->CKTnoncon) break;
        }
        if (!err)
            err = pac_jload(ckt, 0.5 * a1, &g_pac_J2);
    }
    if (!err && check)
        err = pac_jload(ckt, 2.0 * a1, &g_pac_J3);
#ifdef XSPICE
    if (!err && g_hid_nc > 0 && g_hid_mode) {      /* card #557: finite differences of the stateful code models */
        err = pac_hid_stage(ckt, a1, save_rhs, save_state0);
        memcpy(ckt->CKTrhsOld, save_rhs, (size_t) (neq + 1) * sizeof(double));
        ckt->CKTnoncon = 0;
    }
#endif
    if (!err)
        err = pac_jload(ckt, a1, &g_pac_J1);       /* last: the matrix now holds J_k */
    if (!err) {
        if (g_pac_J1.n != g_pac_J2.n || (check && g_pac_J1.n != g_pac_J3.n)) {
            err = 1;
        } else {
            n = g_pac_J1.n;
            if (pac_trip_fit(C, n)) {
                err = 1;
            } else {
                for (i = 0; i < n; i++) {
                    double a = g_pac_J1.v[i], b = g_pac_J2.v[i], diff = a - b;
                    if (g_pac_J1.ri[i] != g_pac_J2.ri[i] || g_pac_J1.ci[i] != g_pac_J2.ci[i]) { err = 1; break; }
                    if (fabs(diff) <= 1e-14 * (fabs(a) + fabs(b)))
                        continue;
                    C->ri[C->n] = g_pac_J1.ri[i];
                    C->ci[C->n] = g_pac_J1.ci[i];
                    C->v[C->n] = diff / d12;
                    C->n++;
                }
            }
        }
    }
    if (!err && check && !pac_affine_ok(&g_pac_J1, &g_pac_J2, &g_pac_J3, NULL)) {
        g_pac_gen_bad = 1;
        pac_gen_identify_bad(ckt, a1);
    }
#ifdef XSPICE
    if (!err && do_factor && g_hid_mode == 2 && g_hid_nh > 0) {
        /* card #557: the model's input columns of J_k are the finite-difference dF/dx (declared PARTIALs are not always the derivative) */
        int k, c, nin = g_hid_nin;
        unsigned char *seen = TMALLOC(unsigned char, (size_t) (neq + 1) * nin);
        int *cmap = TMALLOC(int, neq + 1);
        if (!seen || !cmap) err = 1;
        else {
            memset(seen, 0, (size_t) (neq + 1) * nin);
            for (c = 0; c <= neq; c++) cmap[c] = -1;
            for (k = 0; k < nin; k++) cmap[g_hid_inrow[k]] = k;
            for (i = 0; i < g_pac_J1.n; i++) {
                int r = g_pac_J1.ri[i], cc = g_pac_J1.ci[i];
                if (r > 0 && cc > 0 && cmap[cc] >= 0) {
                    g_pac_J1.v[i] = g_hid_Jc[(size_t) r * nin + cmap[cc]];
                    seen[(size_t) r * nin + cmap[cc]] = 1;
                }
            }
            for (k = 0; k < nin && !err; k++) {
                double mx = 0.0;
                int r;
                for (r = 1; r <= neq; r++) { double v = fabs(g_hid_Jc[(size_t) r * nin + k]); if (v > mx) mx = v; }
                for (r = 1; r <= neq; r++)
                    if (!seen[(size_t) r * nin + k] && fabs(g_hid_Jc[(size_t) r * nin + k]) > 1e-9 * mx) {
                        fprintf(stderr, "PAC-TD: stateful code model: the finite-difference Jacobian column of node row %d has an entry at row %d "
                                "outside the loaded matrix structure (%.3g); abandoning (card #557).\n", g_hid_inrow[k], r, g_hid_Jc[(size_t) r * nin + k]);
                        err = 1;
                    }
            }
            if (getenv("PAC_HID_NO_JCOL")) {   /* red proof: keep the declared PARTIALs (wrong for s_xfer) */
                SMPcopyElements(ckt->CKTmatrix, g_pac_J1.ri, g_pac_J1.ci, g_pac_J1.v);
            }
            else if (!err && SMPsetElements(ckt->CKTmatrix, g_pac_J1.n, g_pac_J1.ri, g_pac_J1.ci, g_pac_J1.v)) err = 1;
        }
        if (seen) tfree(seen);
        if (cmap) tfree(cmap);
    }
#endif
    if (!err && do_factor && !g_pac_gen_bad) {
        /* the last pac_jload() was the J1 load, but pac_gen_identify_bad() (bad path only) clobbers */
        err = SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin) ? 1 : 0;
    }

    ckt->CKTorder = save_order;
    ckt->CKTmode = save_mode;
    ckt->CKTbypass = save_bypass;
    ckt->CKTnoncon = save_noncon;
    memcpy(ckt->CKTrhs, save_rhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
    if (nstates > 0)
        memcpy(ckt->CKTstate0, save_state0, (size_t) nstates * sizeof(double));
    tfree(save_rhs); tfree(save_rhsold);
    if (save_state0) tfree(save_state0);
    return err ? 1 : 0;
}

/* Q = C y for the M complex columns of y (both [(neq+1)*M], row-major by row). */
static void
pac_C_mul(const PACtrip *C, const PACc *y, PACc *Q, int neq, int M)
{
    int i, m;
    memset(Q, 0, (size_t) (neq + 1) * M * sizeof(PACc));
    for (i = 0; i < C->n; i++) {
        int r = C->ri[i], c = C->ci[i];
        double v = C->v[i];
        if (r <= 0 || c <= 0) continue;
        for (m = 0; m < M; m++) {
            Q[(size_t) r * M + m].re += v * y[(size_t) c * M + m].re;
            Q[(size_t) r * M + m].im += v * y[(size_t) c * M + m].im;
        }
    }
}

static void
pac_gen_refuse_message(const char *who)
{
    fprintf(stderr, "%s: REFUSED -- the load of device type '%s' is not affine in the integration "
            "coefficient ag0 (3-point check at ag0/2, ag0, 2*ag0), so its charge matrix C = dq/dx cannot be "
            "separated from G generically. Not returning a curve that would silently omit that "
            "device's charge history (card #542).\n", who, g_pac_gen_badname);
}

/* card #544/#559: saltation with model charges.  A row r of the loaded charge matrix is "native" when the sums of |C| in its row and
 * in its column (diagonal and off-diagonal separately) are those of the native capacitors / inductors on it alone.  A power node of
 * a switch that is NOT native carries device-model charge (diode cj/tt, MOS junction/overlap/intrinsic, ...): since #559 its jump is
 * formed on the charge vector (the charge-rate jump of the switch current, no capacitance in it, then the state jump C dy = dQ over
 * the state rows), and the prepass records such nodes (pac_salq_note()).  Every other row the jump is applied to (a capacitor node
 * or an inductor branch, divided by the native local capacitance / inductance) must be native, else the analysis is REFUSED
 * (a mutual inductance on the branch, a model charge at a row the jump does not treat as a charge row). */
static int
pac_gen_row_native(CKTcircuit *ckt, int r, double *pad, double *pao, double *ped, double *peo)
{
    int captype = CKTtypelook("Capacitor"), indtype = CKTtypelook("Inductor");
    int j;
    double ed = 0.0, eo = 0.0, ad = 0.0, ao = 0.0, ad2 = 0.0, ao2 = 0.0;
    CAPmodel *cm; CAPinstance *ci;
    INDmodel *im; INDinstance *ii;
    if (captype >= 0)
        for (cm = (CAPmodel *) ckt->CKThead[captype]; cm; cm = CAPnextModel(cm))
            for (ci = CAPinstances(cm); ci; ci = CAPnextInstance(ci)) {
                int p = ci->CAPposNode, q = ci->CAPnegNode;
                if (p != r && q != r) continue;
                ed += fabs(ci->CAPcapac);
                if ((p == r ? q : p) > 0) eo += fabs(ci->CAPcapac);
            }
    if (indtype >= 0)
        for (im = (INDmodel *) ckt->CKThead[indtype]; im; im = INDnextModel(im))
            for (ii = INDinstances(im); ii; ii = INDnextInstance(ii))
                if (ii->INDbrEq == r)
                    ed += fabs(ii->INDinduct);
    for (j = 0; j < g_pac_C.n; j++) {
        double v = fabs(g_pac_C.v[j]);
        if (g_pac_C.ri[j] == r) { if (g_pac_C.ci[j] == r) ad += v; else ao += v; }
        if (g_pac_C.ci[j] == r && g_pac_C.ri[j] != r) ao2 += v;
        if (g_pac_C.ci[j] == r && g_pac_C.ri[j] == r) ad2 += v;
    }
    if (pad) *pad = ad;
    if (pao) *pao = ao;
    if (ped) *ped = ed;
    if (peo) *peo = eo;
    return !(fabs(ad - ed) > 1e-6 * (ad + ed) + 1e-30 || fabs(ao - eo) > 1e-6 * (ao + eo) + 1e-30 ||
             fabs(ao2 - eo) > 1e-6 * (ao2 + eo) + 1e-30 || fabs(ad2 - ed) > 1e-6 * (ad2 + ed) + 1e-30);
}

/* card #559: power nodes of ideal switches that carry device-model charge anywhere on the orbit (recorded by the prepass), with the
 * largest total diagonal |C| seen there.  PSSsaltModelCapHook reads them while a generic pass runs. */
static unsigned char *g_salq = NULL;
static double *g_salqC = NULL;
static double *g_salqCmin = NULL;   /* card #575: smallest device-model part (loaded diag |C| - native) of the node on the orbit; -1 = not seen */
static int g_salq_n = 0, g_salq_any = 0;

static int
pac_salq_off(void)
{
    static int off = -1;
    if (off < 0) off = getenv("PAC_SALT_NOQ") != NULL || getenv("PAC_SALT_NATIVEJUMP") != NULL;   /* red proofs: refusal / old native-only jump */
    return off;
}

static double
pac_salq_min_hook(CKTcircuit *ckt, int node)
{
    (void) ckt;
    if (!g_salqCmin || node <= 0 || node > g_salq_n || g_salqCmin[node] < 0.0) return 0.0;
    return g_salqCmin[node];
}

static double
pac_salq_hook(CKTcircuit *ckt, int node)
{
    (void) ckt;
    if (!g_pac_gen || !g_salq || node <= 0 || node > g_salq_n || pac_salq_off())
        return 0.0;
    return g_salq[node] ? g_salqC[node] : 0.0;
}

/* prepass step: note the switch power nodes whose loaded charge row is not native */
static void
pac_salq_note(CKTcircuit *ckt)
{
    int swtype = CKTtypelook("Switch");
    SWmodel *sm; SWinstance *si;
    if (!g_salq || swtype < 0) return;
    for (sm = (SWmodel *) ckt->CKThead[swtype]; sm; sm = SWnextModel(sm))
        for (si = SWinstances(sm); si; si = SWnextInstance(si)) {
            int k, nd[2];
            nd[0] = si->SWposNode; nd[1] = si->SWnegNode;
            for (k = 0; k < 2; k++) {
                double ad, ed, mp;
                int nat;
                if (nd[k] <= 0 || nd[k] > g_salq_n) continue;
                nat = pac_gen_row_native(ckt, nd[k], &ad, NULL, &ed, NULL);
                mp = ad - ed > 0.0 ? ad - ed : 0.0;
                if (g_salqCmin && (g_salqCmin[nd[k]] < 0.0 || mp < g_salqCmin[nd[k]])) g_salqCmin[nd[k]] = mp;
                if (!nat) {
                    g_salq[nd[k]] = 1; g_salq_any = 1;
                    if (ad > g_salqC[nd[k]]) g_salqC[nd[k]] = ad;
                }
            }
        }
}

static int
pac_gen_salt_guard(CKTcircuit *ckt, const PSSsaltEvent *ev, int n_ev)
{
    int i, k;
    if (getenv("PAC_SALT_NATIVEJUMP")) return 0;   /* red proof (card #559): the pre-#559 jump (divided by the native capacitance) without the refusal */
    for (i = 0; i < n_ev; i++) {
        int rows[4];
        /* posRow/negRow are where the jump is applied (capacitor node or inductor branch); the switch's own power nodes carry the
         * charge the jump would have to move when the jump is routed through an inductor (no capacitor at the switched node): a
         * model charge there makes the node a state the native jump does not describe, unless the event treats it as a charge
         * row (posQ/negQ, #559). */
        rows[0] = ev[i].posRow; rows[1] = ev[i].negRow; rows[2] = ev[i].posNode; rows[3] = ev[i].negNode;
        for (k = 0; k < 4; k++) {
            int r = rows[k];
            double ad, ao, ed, eo;
            if (r <= 0) continue;
            if ((r == ev[i].posNode && ev[i].posF) || (r == ev[i].negNode && ev[i].negF))
                continue;   /* card #575: fast node, the jump is on the inductor branch row */
            if ((r == ev[i].posNode && ev[i].posQ) || (r == ev[i].negNode && ev[i].negQ)) {
                /* card #559 validated envelope: the model charge sits on top of a NATIVE capacitor at the same node (a switch-node
                 * capacitor, Coss modelled as a C, ...), so the node is a charge state on the whole orbit and its fast transient is
                 * resolved by the steps.  A node whose only charge is device-model charge (reverse-biased diode: cj only, forward: tt*I)
                 * and that is small (tens of pF) turns the event into a ps-scale unresolved transient: pac_r6od with the node capacitor
                 * removed 0.06-0.10 dB off and drifting with reltol, a diode-only inductor-adjacent node 0.30 dB off and 3 dB off at reltol
                 * 1e-6, while the native-capacitor decks hold 0.013-0.02 dB at every reltol: refused (PAC_SALT_QONLY=1 overrides). */
                double ad, ao, ed, eo;
                (void) pac_gen_row_native(ckt, r, &ad, &ao, &ed, &eo);
                if ((ed <= 0.0 || ad <= 0.0) && !getenv("PAC_SALT_QONLY")) {
                    fprintf(stderr, "PAC-TD: REFUSED -- an ideal-switch (saltation) event (%s) at a node that carries device-model charge but "
                            "no native capacitor (loaded |C| diag %g, native %g): the state jump of such a node is validated only with a native "
                            "capacitor at the node (card #559: 0.06-0.3 dB and reltol-dependent without; PAC_SALT_QONLY=1 overrides, "
                            "unvalidated).\n", ev[i].name ? ev[i].name : "sw", ad, ed);
                    return 1;
                }
                continue;
            }
            if (!pac_gen_row_native(ckt, r, &ad, &ao, &ed, &eo)) {
                fprintf(stderr, "PAC-TD: REFUSED -- an ideal-switch (saltation) event (%s) at a row that carries device-model "
                        "charge or a mutual inductance (loaded |C| row %g/%g, native %g/%g): the state jump of the event is "
                        "derived for the native capacitor/inductor at the switched node only (cards #544, #559).\n",
                        ev[i].name ? ev[i].name : "sw", ad, ao, ed, eo);
                return 1;
            }
        }
    }
    return 0;
}

/* card #555: the matrix is real, so a complex right-hand side is solved as two real solves (re, im).  The homogeneous
 * (monodromy) columns of recording pass 1 -- n of the n + nfreq columns, 413..1500 on R7p -- have an exactly zero
 * imaginary part for the whole period (real unit-vector start, no forcing, real C_k): the second solve of a zero vector
 * returns zeros and costs as much as the first.  It is skipped when the imaginary right-hand side is identically zero
 * (same values; the sign of a zero may differ).  PAC_NO_ZERO_IM_SKIP=1 restores both solves. */
static void
pac_solve_cplx(CKTcircuit *ckt)
{
    static int skip = -1;
    int row, neq = g_pac_neq;
    if (skip < 0) skip = getenv("PAC_NO_ZERO_IM_SKIP") == NULL;
    SMPsolve(ckt->CKTmatrix, g_pac_rhs_re, g_pac_scratch);
    if (skip) {
        for (row = 0; row <= neq; row++)
            if (g_pac_rhs_im[row] != 0.0) break;
        if (row > neq) return;
    }
    SMPsolve(ckt->CKTmatrix, g_pac_rhs_im, g_pac_scratch);
}


/* card #558, round 17: the M columns of one step are independent solves against the same factored matrix.  They can be
 * collected (pac_mq_put) and solved with one multi-right-hand-side call (SMPsolveMulti: one klu_solve with M columns
 * under KLU, so L and U are read once per block of columns, not once per column).  The curves are byte-identical to the
 * column-by-column solves (tests/test_accel17.py), but the wall time is not measurably different on R7p .pstb (8
 * alternations on one machine: -3 %, inside the +-20 % run-to-run outliers of the machine): OPT-IN, PAC_MULTIRHS=1.
 * The imaginary parts that are identically zero are skipped exactly as in pac_solve_cplx. */
static double *g_pac_mq_re = NULL, *g_pac_mq_im = NULL, *g_pac_mq_tmp = NULL;
static size_t g_pac_mq_cap = 0;

static int
pac_mq_on(void)
{
    static int on = -1;
    if (on < 0) on = getenv("PAC_MULTIRHS") != NULL && atoi(getenv("PAC_MULTIRHS")) != 0;
    return on;
}

static int pac_hid_active(void);   /* card #557/#558: defined with the hidden-state code */

static void
pac_mq_put(int m)
{
    int ld = g_pac_neq + 1;
    size_t need = (size_t) ld * g_pac_M;
    if (g_pac_mq_cap < need) {
        g_pac_mq_re = TREALLOC(double, g_pac_mq_re, need);
        g_pac_mq_im = TREALLOC(double, g_pac_mq_im, need);
        g_pac_mq_tmp = TREALLOC(double, g_pac_mq_tmp, need);
        g_pac_mq_cap = need;
    }
    memcpy(g_pac_mq_re + (size_t) m * ld, g_pac_rhs_re, (size_t) ld * sizeof(double));
    memcpy(g_pac_mq_im + (size_t) m * ld, g_pac_rhs_im, (size_t) ld * sizeof(double));
}

static void
pac_mq_flush(CKTcircuit *ckt, int M)
{
    static int skip = -1;
    int ld = g_pac_neq + 1, neq = g_pac_neq, m, row, nz = 0;
    if (skip < 0) skip = getenv("PAC_NO_ZERO_IM_SKIP") == NULL;
    SMPsolveMulti(ckt->CKTmatrix, g_pac_mq_re, ld, M);
    if (skip) {
        for (m = 0; m < M; m++) {
            const double *c = g_pac_mq_im + (size_t) m * ld;
            for (row = 0; row <= neq; row++)
                if (c[row] != 0.0) break;
            if (row > neq) continue;
            memcpy(g_pac_mq_tmp + (size_t) nz * ld, c, (size_t) ld * sizeof(double));
            nz++;
        }
        if (nz) {
            SMPsolveMulti(ckt->CKTmatrix, g_pac_mq_tmp, ld, nz);
            nz = 0;
            for (m = 0; m < M; m++) {
                double *c = g_pac_mq_im + (size_t) m * ld;
                for (row = 0; row <= neq; row++)
                    if (c[row] != 0.0) break;
                if (row > neq) continue;
                memcpy(c, g_pac_mq_tmp + (size_t) nz * ld, (size_t) ld * sizeof(double));
                nz++;
            }
        }
    } else
        SMPsolveMulti(ckt->CKTmatrix, g_pac_mq_im, ld, M);
    for (m = 0; m < M; m++)
        for (row = 1; row <= neq; row++) {
            PACc v; v.re = g_pac_mq_re[(size_t) m * ld + row]; v.im = g_pac_mq_im[(size_t) m * ld + row];
            pac_set(g_pac_Cnew, row, M, m, v);
        }
}

#ifdef XSPICE
/* card #557: hidden-state terms of one step for column m: rhs -= B dh_{k-1} (before the solve), dh_k = P y_k + R dh_{k-1} (after). */
static int
pac_hid_active(void)
{
    return g_hid_nh > 0 && g_hid_mode == 2;
}

static void
pac_hid_rhs(int m, int M)
{
    int neq = g_pac_neq, a, row, nh = g_hid_nh;
    if (nh <= 0 || g_hid_mode != 2 || getenv("PAC_HID_RED")) return;   /* PAC_HID_RED=1: red proof, the hidden states do not act on the nodes */
    for (a = 0; a < nh; a++) {
        PACc hv = g_pac_Cprev[(size_t) (neq + 1 + a) * M + m];
        if (hv.re == 0.0 && hv.im == 0.0) continue;
        for (row = 1; row <= neq; row++) {
            double b = g_hid_Bh[(size_t) row * nh + a];
            if (b == 0.0) continue;
            g_pac_rhs_re[row] -= b * hv.re;
            g_pac_rhs_im[row] -= b * hv.im;
        }
    }
}

static void
pac_hid_update(int m, int M)
{
    int neq = g_pac_neq, a, b, i, nh = g_hid_nh, nin = g_hid_nin;
    if (nh <= 0 || g_hid_mode != 2) return;
    for (a = 0; a < nh; a++) {
        PACc v = { 0.0, 0.0 };
        for (i = 0; i < nin; i++) {
            double p = g_hid_P[(size_t) a * nin + i];
            PACc y;
            if (p == 0.0) continue;
            y = g_pac_Cnew[(size_t) g_hid_inrow[i] * M + m];
            v.re += p * y.re; v.im += p * y.im;
        }
        for (b = 0; b < nh; b++) {
            double r = g_hid_R[(size_t) a * nh + b];
            PACc y;
            if (r == 0.0) continue;
            y = g_pac_Cprev[(size_t) (neq + 1 + b) * M + m];
            v.re += r * y.re; v.im += r * y.im;
        }
        g_pac_Cnew[(size_t) (neq + 1 + a) * M + m] = v;
    }
}
#else
static int pac_hid_active(void) { return 0; }
static void pac_hid_rhs(int m, int M) { (void) m; (void) M; }
static void pac_hid_update(int m, int M) { (void) m; (void) M; }
#endif

/* card #559: ideal-switch event at a power node that carries device-model charge.  Fang's jump on the charge vector:
 *     dq+ = dq- - (F+ - F-) dt_e ,   dt_e = -(grad h . dx-) / h_dot ,   F+ - F- = (-dI at pos, +dI at neg)  [charge rate: the switch current]
 * so dQ(row) = coeff * dqdt with coeff = (y(posCntrl) - y(negCntrl)) / h_dot exactly as the native branch, but nothing is divided by a
 * capacitance.  The state jump that goes with it is C dy = dQ over the charge-coupled block of state rows (C_k at the event: a model
 * capacitance couples the node to its neighbours, e.g. Cgd to the gate: the charge of the neighbour does not jump, its voltage does,
 * a capacitive divider); y of those rows is moved by dy so the recorded perturbation and the Floquet closure see the jump.  Q_k = C_k y_k
 * is then Q_smooth + dQ.  Where the block is singular (a floating pair of nodes) y is moved by dQ / (diagonal C); Q_k is Q_smooth + dQ exactly
 * either way.  T_k = ag0 (Q_smooth - Q_{k-1}) - a1 T_{k-1} is formed from the SMOOTH part BEFORE the jump is added (the derivative of the flow
 * after the event, not dQ/h: the trapezoidal recursion would re-inject the spike as a1 T_k into the next step), so no forced order-1 window is
 * needed (PSS_SALT_QSTEPS=N forces N steps, experiments only: N = 0..2 give the same curve to 1e-3 dB on pac_r6od; the earlier tau-window
 * scheme moved the answer by 3 dB with the window length).
 * Returns 0 (no model-charged power node on these events), 1 (formed; pac_salt_qjump_finish() adds dQ and dy after T_k), -1 (no memory). */
static int
pac_salt_qjump_prepare(const PSSsaltEvent *ev, int n_ev, int neq, int M, PACc **pdq, PACc **pdy)
{
    int i, m, any = 0;
    size_t sz = (size_t) (neq + 1) * M;
    PACc *dq, *dy;
    unsigned char *inb, *isst;
    int *idx, *lst, nb = 0, ns, grew, j;
    double *A = NULL;
    for (i = 0; i < n_ev; i++)
        if (ev[i].posQ || ev[i].negQ) any = 1;
    if (!any) return 0;
    dq = TMALLOC(PACc, sz); dy = TMALLOC(PACc, sz);
    inb = TMALLOC(unsigned char, neq + 1); isst = TMALLOC(unsigned char, neq + 1);
    idx = TMALLOC(int, neq + 1); lst = TMALLOC(int, neq + 1);
    if (!dq || !dy || !inb || !isst || !idx || !lst) return -1;
    memset(dq, 0, sz * sizeof(PACc)); memset(dy, 0, sz * sizeof(PACc));
    memset(inb, 0, (size_t) neq + 1); memset(isst, 0, (size_t) neq + 1);
    ns = g_pac_map->ns;
    for (i = 0; i < ns; i++) { int r = g_pac_map->srows[i]; if (r > 0 && r <= neq) isst[r] = 1; }
    for (m = 0; m < M; m++)
        for (i = 0; i < n_ev; i++) {
            const PSSsaltEvent *e = &ev[i];
            PACc Sp = pac_get(g_pac_Cprev, e->posCntrl, M, m), Sn = pac_get(g_pac_Cprev, e->negCntrl, M, m);
            PACc coeff = pc_scale(pc_sub(Sp, Sn), 1.0 / e->h_dot);
            if (e->posQ && e->posRow > 0)
                dq[(size_t) e->posRow * M + m] = pc_add(dq[(size_t) e->posRow * M + m], pc_scale(coeff, e->dqdt_pos));
            if (e->negQ && e->negRow > 0)
                dq[(size_t) e->negRow * M + m] = pc_add(dq[(size_t) e->negRow * M + m], pc_scale(coeff, e->dqdt_neg));
        }
    /* charge-coupled block of state rows containing the jumped rows */
    for (i = 0; i < n_ev; i++) {
        if (ev[i].posQ && ev[i].posRow > 0 && isst[ev[i].posRow]) inb[ev[i].posRow] = 1;
        if (ev[i].negQ && ev[i].negRow > 0 && isst[ev[i].negRow]) inb[ev[i].negRow] = 1;
    }
    do {
        grew = 0;
        for (j = 0; j < g_pac_C.n; j++) {
            int r = g_pac_C.ri[j], cc = g_pac_C.ci[j];
            if (r <= 0 || cc <= 0 || g_pac_C.v[j] == 0.0 || !isst[r] || !isst[cc]) continue;
            if (inb[r] && !inb[cc]) { inb[cc] = 1; grew = 1; }
            if (inb[cc] && !inb[r]) { inb[r] = 1; grew = 1; }
        }
    } while (grew);
    for (i = 1; i <= neq; i++) { idx[i] = -1; if (inb[i]) { idx[i] = nb; lst[nb++] = i; } }
    {
        int ok = 0;
        if (nb > 0) A = TMALLOC(double, (size_t) nb * nb);
        if (A) {
            int a, b, k;
            memset(A, 0, (size_t) nb * nb * sizeof(double));
            for (j = 0; j < g_pac_C.n; j++) {
                int r = g_pac_C.ri[j], cc = g_pac_C.ci[j];
                if (r > 0 && cc > 0 && idx[r] >= 0 && idx[cc] >= 0) A[(size_t) idx[r] * nb + idx[cc]] += g_pac_C.v[j];
            }
            {   /* LU with partial pivoting kept in place, solved per column m */
                double amax = 0.0;
                int *piv = TMALLOC(int, nb);
                for (a = 0; a < nb * nb; a++) if (fabs(A[a]) > amax) amax = fabs(A[a]);
                ok = (piv != NULL && amax > 0.0);
                for (k = 0; ok && k < nb; k++) {
                    int p = k;
                    double pm = fabs(A[(size_t) k * nb + k]);
                    for (a = k + 1; a < nb; a++) if (fabs(A[(size_t) a * nb + k]) > pm) { pm = fabs(A[(size_t) a * nb + k]); p = a; }
                    if (pm <= 1e-10 * amax) { ok = 0; break; }
                    piv[k] = p;
                    if (p != k) for (b = 0; b < nb; b++) { double t = A[(size_t) k * nb + b]; A[(size_t) k * nb + b] = A[(size_t) p * nb + b]; A[(size_t) p * nb + b] = t; }
                    for (a = k + 1; a < nb; a++) {
                        double f = A[(size_t) a * nb + k] / A[(size_t) k * nb + k];
                        A[(size_t) a * nb + k] = f;
                        for (b = k + 1; b < nb; b++) A[(size_t) a * nb + b] -= f * A[(size_t) k * nb + b];
                    }
                }
                if (ok)
                    for (m = 0; m < M; m++) {
                        PACc *x = TMALLOC(PACc, nb);
                        if (!x) { ok = 0; break; }
                        for (a = 0; a < nb; a++) x[a] = dq[(size_t) lst[a] * M + m];
                        for (k = 0; k < nb; k++) {   /* apply the row permutation and forward-eliminate */
                            PACc t = x[k]; x[k] = x[piv[k]]; x[piv[k]] = t;
                            for (a = k + 1; a < nb; a++) {
                                x[a].re -= A[(size_t) a * nb + k] * x[k].re; x[a].im -= A[(size_t) a * nb + k] * x[k].im;
                            }
                        }
                        for (k = nb - 1; k >= 0; k--) {
                            for (b = k + 1; b < nb; b++) { x[k].re -= A[(size_t) k * nb + b] * x[b].re; x[k].im -= A[(size_t) k * nb + b] * x[b].im; }
                            x[k].re /= A[(size_t) k * nb + k]; x[k].im /= A[(size_t) k * nb + k];
                        }
                        for (a = 0; a < nb; a++) dy[(size_t) lst[a] * M + m] = x[a];
                        tfree(x);
                    }
                if (piv) tfree(piv);
            }
        }
        if (!ok) {   /* singular block: dy = dQ / diagonal C on the jumped rows (the residual is added to Q exactly) */
            memset(dy, 0, sz * sizeof(PACc));
            for (i = 0; i < n_ev; i++) {
                int sd, r = 0;
                for (sd = 0; sd < 2; sd++) {
                    double dg = 0.0;
                    r = sd ? (ev[i].negQ ? ev[i].negRow : 0) : (ev[i].posQ ? ev[i].posRow : 0);
                    if (r <= 0) continue;
                    for (j = 0; j < g_pac_C.n; j++) if (g_pac_C.ri[j] == r && g_pac_C.ci[j] == r) dg += g_pac_C.v[j];
                    if (dg == 0.0) continue;
                    for (m = 0; m < M; m++)
                        dy[(size_t) r * M + m] = pc_add(dy[(size_t) r * M + m], pc_scale(dq[(size_t) r * M + m], 1.0 / dg));
                }
            }
        }
    }
    if (A) tfree(A);
    tfree(inb); tfree(isst); tfree(idx); tfree(lst);
    *pdq = dq; *pdy = dy;
    return 1;
}

/* After T_k has been formed from the SMOOTH part of the step (Q_smooth = C y_smooth, so the derivative history carries no jump spike: the derivative
 * of the flow after the event, not dQ/h -- the trapezoidal recursion would otherwise re-inject the spike as a1 T_k into the next step), add the
 * jump itself: Q_k = Q_smooth + dQ (exactly, whatever the state jump dy could resolve) and y_k = y_smooth + dy. */
static void
pac_salt_qjump_finish(int neq, int M, PACc *dq, PACc *dy)
{
    size_t sz = (size_t) (neq + 1) * M, i;
    for (i = 0; i < sz; i++) {
        g_pac_Qnew[i].re += dq[i].re; g_pac_Qnew[i].im += dq[i].im;
        g_pac_Cnew[i].re += dy[i].re; g_pac_Cnew[i].im += dy[i].im;
    }
    tfree(dq); tfree(dy);
}

/* card #583: PAC_SALT_SMOOTHT=0 restores jumps-before-T with the order-1 window (red proof / comparison) */
static int
pac_smooth_t(void)
{
    static int v = -1;
    if (v < 0) { const char *e = getenv("PAC_SALT_SMOOTHT"); v = e ? atoi(e) : 1; }
    return v;
}

/* One accepted step of the generic recursion for all M columns. */
static int
pac_step_generic(CKTcircuit *ckt, int first, double ag0, double ag1, int order, double t_new,
                const PSSsaltEvent *salt_ev, int n_salt)
{
    int neq = g_pac_neq, M = g_pac_M, m, row, qjump = 0;
    PACc *dq = NULL, *dy = NULL;
    PACc *dy2 = NULL;
    int smooth_t = pac_smooth_t();
    double a1 = (order == 2) ? ag1 : 0.0;

#ifdef XSPICE
    if (g_hid_nh > 0) {
        g_hid_mode = 2;
        if (n_salt > 0) {
            fprintf(stderr, "PAC-TD: REFUSED -- an ideal-switch (saltation) event with a stateful XSPICE code model in the circuit: the "
                    "event jump of the hidden states is not derived (card #557).\n");
            return 1;
        }
    }
#endif
    if (pac_gen_loads(ckt, ag0, g_pac_gen_check, 1)) {
        fprintf(stderr, "PAC-TD: generic charge extraction/refactor failed at t=%g; abandoning this pass.\n", t_new);
        return 1;
    }
    if (g_pac_gen_bad) {
        pac_gen_refuse_message("PAC-TD");
        return 1;
    }
    if (n_salt > 0 && pac_gen_salt_guard(ckt, salt_ev, n_salt))
        return 1;
    if (first && getenv("PAC_GEN_TRACE"))
        fprintf(stderr, "PAC-TD: generic recursion, first step: nnz(C)=%d nnz(J)=%d M=%d ag0=%g\n",
                g_pac_C.n, g_pac_J1.n, M, ag0);
    if (first) {   /* Q_0 = C_1 y_0 (C_0 := C_1, see the design note), T_0 = 0 (first step is order 1) */
        pac_C_mul(&g_pac_C, g_pac_Cprev, g_pac_Qprev, neq, M);
        if (g_pac_gear) {
            int jj;
            for (jj = 2; jj <= PAC_GEAR_MAXORD; jj++)
                memcpy(g_pac_Qh[jj], g_pac_Qprev, (size_t) (neq + 1) * M * sizeof(PACc));
        }
    }
    if (g_pac_gear) {
        int pord = 0, jj;
        if (n_salt > 0) {
            fprintf(stderr, "PAC-TD: REFUSED -- an ideal-switch (saltation) event with method=gear: the event's state jump is "
                    "derived for the trapezoidal derivative history only (card #549); use method=trap (default).\n");
            return 1;
        }
        for (jj = 1; jj <= PAC_GEAR_MAXORD; jj++)
            if (ckt->CKTag[jj] != 0.0) pord = jj;
        if (pord < 1) pord = 1;
        {   /* red-proof hook (tests/test_pac_charge.py check_gear): PAC_GEAR_TEST_ORD1=1 drops every history term beyond Q_{k-1} */
            static int ord1 = -1;
            if (ord1 < 0) ord1 = getenv("PAC_GEAR_TEST_ORD1") ? 1 : 0;
            if (ord1) pord = 1;
        }
        for (m = 0; m < M; m++) {
            memset(g_pac_rhs_re, 0, (size_t) (neq + 1) * sizeof(double));
            memset(g_pac_rhs_im, 0, (size_t) (neq + 1) * sizeof(double));
            for (row = 1; row <= neq; row++) {
                double sr = 0.0, si = 0.0;
                for (jj = 1; jj <= pord; jj++) {
                    const PACc *Qj = (jj == 1 ? g_pac_Qprev : g_pac_Qh[jj]) + (size_t) row * M + m;
                    sr -= ckt->CKTag[jj] * Qj->re;
                    si -= ckt->CKTag[jj] * Qj->im;
                }
                g_pac_rhs_re[row] = sr;
                g_pac_rhs_im[row] = si;
            }
            if (g_pac_useforce[m] && g_pac_br && g_pac_bi) {
                double w = g_pac_omega[m];
                double ph = w * t_new;
                double c = cos(ph), s = sin(ph);
                for (row = 1; row <= neq; row++) {
                    double br = g_pac_br[row], bi = g_pac_bi[row];
                    if (br == 0.0 && bi == 0.0) continue;
                    g_pac_rhs_re[row] += br * c - bi * s;
                    g_pac_rhs_im[row] += br * s + bi * c;
                }
            }
            pac_hid_rhs(m, M);
            if (pac_mq_on() && !pac_hid_active())
                pac_mq_put(m);   /* card #558: solved for all columns after the loop (not with #557 hidden states: they update per column after the solve) */
            else {
                pac_solve_cplx(ckt);
                for (row = 1; row <= neq; row++) {
                    PACc v; v.re = g_pac_rhs_re[row]; v.im = g_pac_rhs_im[row];
                    pac_set(g_pac_Cnew, row, M, m, v);
                }
            }
            pac_hid_update(m, M);
        }
        if (pac_mq_on() && !pac_hid_active())
            pac_mq_flush(ckt, M);
        pac_C_mul(&g_pac_C, g_pac_Cnew, g_pac_Qnew, neq, M);
        {   /* shift the charge history: Q_{k-j} <- Q_{k-j+1}, Q_{k-1} <- Q_k */
            PACc *oldest = g_pac_Qh[PAC_GEAR_MAXORD];
            for (jj = PAC_GEAR_MAXORD; jj > 2; jj--) g_pac_Qh[jj] = g_pac_Qh[jj - 1];
            g_pac_Qh[2] = g_pac_Qprev;
            g_pac_Qprev = g_pac_Qnew;
            g_pac_Qnew = oldest;
        }
        return 0;
    }

    for (m = 0; m < M; m++) {
        memset(g_pac_rhs_re, 0, (size_t) (neq + 1) * sizeof(double));
        memset(g_pac_rhs_im, 0, (size_t) (neq + 1) * sizeof(double));
        for (row = 1; row <= neq; row++) {
            PACc Qp = g_pac_Qprev[(size_t) row * M + m];
            PACc Tp = g_pac_Tprev[(size_t) row * M + m];
            g_pac_rhs_re[row] = ag0 * Qp.re + a1 * Tp.re;
            g_pac_rhs_im[row] = ag0 * Qp.im + a1 * Tp.im;
        }
        if (g_pac_useforce[m] && g_pac_br && g_pac_bi) {
            double w = g_pac_omega[m];
            double ph = w * t_new;
            double c = cos(ph), s = sin(ph);
            for (row = 1; row <= neq; row++) {
                double br = g_pac_br[row], bi = g_pac_bi[row];
                if (br == 0.0 && bi == 0.0) continue;
                g_pac_rhs_re[row] += br * c - bi * s;
                g_pac_rhs_im[row] += br * s + bi * c;
            }
        }
        pac_hid_rhs(m, M);
        if (pac_mq_on() && !pac_hid_active())
            pac_mq_put(m);   /* card #558: solved for all columns after the loop (not with #557 hidden states: they update per column after the solve) */
        else {
            pac_solve_cplx(ckt);
            for (row = 1; row <= neq; row++) {
                PACc v; v.re = g_pac_rhs_re[row]; v.im = g_pac_rhs_im[row];
                pac_set(g_pac_Cnew, row, M, m, v);
            }
        }
        pac_hid_update(m, M);
    }
    if (pac_mq_on() && !pac_hid_active())
        pac_mq_flush(ckt, M);
    /* card #544: saltation (ideal sw/csw event) with model charges on the orbit.  The state jump of the native recursion
     * (Fang: delta_x+ = delta_x- - (f+ - f-) dt_e, dt_e = -(grad h . delta_x-)/h_dot, applied as a direct value jump at the
     * event's capacitor-node / inductor-branch rows) is exactly what the generic recursion needs too, once the charge vector
     * is formed from the JUMPED solution: Q_k = C_k y_k then carries dq = C_k (jump) and T_k = ag0 (Q_k - Q_{k-1}) - a1 T_{k-1}
     * the derivative history of it, as in the native branch.  It is only exact where the loaded charge matrix at the event rows
     * is that of the native capacitors/inductors alone (pac_gen_salt_guard() refuses otherwise), because the jump magnitude
     * (-dG v_sw / C_local) divides by the native local capacitance. */
    if (n_salt > 0) {
        int i_salt;
        for (m = 0; m < M; m++) {
            PACc djump[PSS_SALT_MAX_EVENTS][2];
            for (i_salt = 0; i_salt < n_salt; i_salt++) {
                const PSSsaltEvent *e = &salt_ev[i_salt];
                PACc Sp = pac_get(g_pac_Cprev, e->posCntrl, M, m);
                PACc Sn = pac_get(g_pac_Cprev, e->negCntrl, M, m);
                PACc coeff = pc_scale(pc_sub(Sp, Sn), 1.0 / e->h_dot);
                djump[i_salt][0] = pc_scale(coeff, e->dfdt_pos);
                djump[i_salt][1] = pc_scale(coeff, e->dfdt_neg);
            }
            for (i_salt = 0; i_salt < n_salt; i_salt++) {
                const PSSsaltEvent *e = &salt_ev[i_salt];
                /* card #583: a routed jump (inductor branch row, a capacitor fed through a snubber or a diode, a continuation) is
                 * added AFTER T has been formed from the smooth part, like the charge jump of #559: T then carries no dy / h spike
                 * and no order-1 step is needed after it.  A native capacitor jump at the switch's own node keeps the old order. */
                int rt_pos = smooth_t && e->posRow > 0 && !e->posQ && (e->dG == 0.0 || e->posRow != e->posNode);
                int rt_neg = smooth_t && e->negRow > 0 && !e->negQ && (e->dG == 0.0 || e->negRow != e->negNode);
                if (rt_pos || rt_neg) {
                    if (!dy2) { dy2 = TMALLOC(PACc, (size_t) (neq + 1) * M); if (!dy2) return 1; memset(dy2, 0, (size_t) (neq + 1) * M * sizeof(PACc)); }
                }
                if (rt_pos) dy2[(size_t) e->posRow * M + m] = pc_add(dy2[(size_t) e->posRow * M + m], djump[i_salt][0]);
                else if (e->posRow > 0 && !e->posQ)
                    pac_set(g_pac_Cnew, e->posRow, M, m, pc_add(pac_get(g_pac_Cnew, e->posRow, M, m), djump[i_salt][0]));
                if (rt_neg) dy2[(size_t) e->negRow * M + m] = pc_add(dy2[(size_t) e->negRow * M + m], djump[i_salt][1]);
                else if (e->negRow > 0 && !e->negQ)
                    pac_set(g_pac_Cnew, e->negRow, M, m, pc_add(pac_get(g_pac_Cnew, e->negRow, M, m), djump[i_salt][1]));
            }
        }
        qjump = pac_salt_qjump_prepare(salt_ev, n_salt, neq, M, &dq, &dy);
        if (qjump < 0) return 1;
    }
    pac_C_mul(&g_pac_C, g_pac_Cnew, g_pac_Qnew, neq, M);
    for (row = 1; row <= neq; row++)
        for (m = 0; m < M; m++) {
            size_t ix = (size_t) row * M + m;
            PACc Tn;
            Tn.re = ag0 * (g_pac_Qnew[ix].re - g_pac_Qprev[ix].re) - a1 * g_pac_Tprev[ix].re;
            Tn.im = ag0 * (g_pac_Qnew[ix].im - g_pac_Qprev[ix].im) - a1 * g_pac_Tprev[ix].im;
            g_pac_Tnew[ix] = Tn;
        }
    if (qjump > 0) pac_salt_qjump_finish(neq, M, dq, dy);   /* card #559: after T (see there) */
    if (dy2) {   /* card #583 */
        PACc *dq2 = TMALLOC(PACc, (size_t) (neq + 1) * M);
        size_t i2, sz2 = (size_t) (neq + 1) * M;
        if (!dq2) { tfree(dy2); return 1; }
        pac_C_mul(&g_pac_C, dy2, dq2, neq, M);
        for (i2 = 0; i2 < sz2; i2++) {
            g_pac_Qnew[i2].re += dq2[i2].re; g_pac_Qnew[i2].im += dq2[i2].im;
            g_pac_Cnew[i2].re += dy2[i2].re; g_pac_Cnew[i2].im += dy2[i2].im;
        }
        tfree(dq2); tfree(dy2);
    }
    { PACc *tmp = g_pac_Qprev; g_pac_Qprev = g_pac_Qnew; g_pac_Qnew = tmp; }
    return 0;
}

/* card #542: keep this step's exact J and C (the forward pass just used them) with the trajectory entry
 * that pac_traj_append() created for it. Refuses (sets g_pac_traj_failed) beyond PAC_TRAJ_MAX_MB. */
static void
pac_traj_attach_generic(void)
{
    PACtrajStep *st;
    double bytes;
    const char *mx = getenv("PAC_TRAJ_MAX_MB");
    double cap = (mx && atof(mx) > 0.0) ? atof(mx) * 1048576.0 : 3.0e9;
    if (g_pac_traj_failed || g_pac_traj_n < 1)
        return;
    st = &g_pac_traj[g_pac_traj_n - 1];
    bytes = (double) g_pac_J1.n * 16.0 + (double) g_pac_C.n * 16.0;
    if (g_pac_traj_bytes + bytes > cap) {
        fprintf(stderr, "PAC-TD: adjoint/noise trajectory would exceed %.0f MB (J_k and C_k of every step are kept "
                "for the exact transpose; set PAC_TRAJ_MAX_MB higher or use fewer pacpoints).\n", cap / 1048576.0);
        g_pac_traj_failed = 1;
        return;
    }
    st->nJ = g_pac_J1.n; st->nC = g_pac_C.n;
    st->Jri = TMALLOC(int, st->nJ + 1); st->Jci = TMALLOC(int, st->nJ + 1); st->Jv = TMALLOC(double, st->nJ + 1);
    st->Cri = TMALLOC(int, st->nC + 1); st->Cci = TMALLOC(int, st->nC + 1); st->Cv = TMALLOC(double, st->nC + 1);
    if (!st->Jri || !st->Jci || !st->Jv || !st->Cri || !st->Cci || !st->Cv) { g_pac_traj_failed = 1; return; }
    memcpy(st->Jri, g_pac_J1.ri, (size_t) st->nJ * sizeof(int));
    memcpy(st->Jci, g_pac_J1.ci, (size_t) st->nJ * sizeof(int));
    memcpy(st->Jv, g_pac_J1.v, (size_t) st->nJ * sizeof(double));
    memcpy(st->Cri, g_pac_C.ri, (size_t) st->nC * sizeof(int));
    memcpy(st->Cci, g_pac_C.ci, (size_t) st->nC * sizeof(int));
    memcpy(st->Cv, g_pac_C.v, (size_t) st->nC * sizeof(double));
    g_pac_traj_bytes += bytes;
}

/* card #579: the impulse an event leaves in an ALGEBRAIC output.  A quantity that is discontinuous at the event (the output of a boost
 * through its capacitor's ESR: i_D jumps by I_L; a switch current; the switch node itself when it has no capacitor) changes, when the
 * event moves by dt_e, by a pulse of area (y+ - y-) x (-dt_e): an impulse in the perturbation that no state carries and the recursion
 * never shows.  Its sideband-0 weight is J x coeff x exp(-j w t_e), J = the large-signal change of the row across the event step
 * (and across each continuation step of a fast node), coeff the event coefficient the state jumps use.  Rows that get a stiff
 * capacitor-node state jump at the switch itself are left out (the decaying spike of the recursion stands for the same thing there).
 * Buck outputs are continuous at the event (J ~ 0): unchanged.  PAC_NO_IMPULSE=1 switches it off (red proof). */
static double *g_pac_xls_prev = NULL;
static int g_pac_xls_n = 0, g_pac_xls_have = 0;
static int g_imp_on = 0, g_imp_pc = 0, g_imp_nc = 0, g_imp_nx = 0, g_imp_x[4 * PSS_SALT_MAX_EVENTS];
static double g_imp_hdot = 0.0;

/* card #581, #711: discontinuous conduction.  When an inductor's current runs out, the diode that carried it stops conducting somewhere
 * inside an accepted step; no switch flips and the step control does not treat it as an event.  Until card #711 the result of a DCM buck
 * was first order in the step (0.40 dB at 200 points per period): the cause was the pre-event node equation of the switch that turns on
 * (pss_salt_ind_jump()), which used the inductor current of the event step instead of the one before the flip.  It is fixed there; the
 * note stays as a pointer for decks with this behaviour, which were validated on a buck only.  Detected as an inductor whose current
 * stays below 1e-3 of its largest value for more than 1 % of the period; said once.  PAC_NO_DCM_NOTE=1 silences it. */
#define PAC_DCM_NIND 64
static double g_dcm_imax[PAC_DCM_NIND], g_dcm_tzero[PAC_DCM_NIND], g_dcm_tprev = -1.0;
static int g_dcm_said = 0;

static void
pac_dcm_note(CKTcircuit *ckt, int n_salt)
{
    int it = CKTtypelook("Inductor"), k = 0;
    INDmodel *im; INDinstance *ii;
    double h = (g_dcm_tprev >= 0.0 && ckt->CKTtime > g_dcm_tprev) ? ckt->CKTtime - g_dcm_tprev : 0.0;
    (void) n_salt;
    if (ckt->CKTtime < g_dcm_tprev) memset(g_dcm_tzero, 0, sizeof g_dcm_tzero);   /* a new period run */
    g_dcm_tprev = ckt->CKTtime;
    if (it < 0 || g_dcm_said || g_pac_T <= 0.0) return;
    for (im = (INDmodel *) ckt->CKThead[it]; im; im = INDnextModel(im))
        for (ii = INDinstances(im); ii && k < PAC_DCM_NIND; ii = INDnextInstance(ii), k++) {
            double a = fabs(ckt->CKTrhsOld[ii->INDbrEq]);
            if (a > g_dcm_imax[k]) g_dcm_imax[k] = a;
            if (g_dcm_imax[k] > 1e-9 && a < 1e-3 * g_dcm_imax[k]) g_dcm_tzero[k] += h;
            if (g_dcm_tzero[k] > 0.01 * g_pac_T && !getenv("PAC_NO_DCM_NOTE")) {
                fprintf(stderr, "PAC-TD: note -- the current of inductor %s is zero for part of the period (discontinuous conduction): the "
                        "instant its diode stops conducting is not located as an event; validated on a DCM buck "
                        "(0.02 dB at 200 points per period, card #711), not on other topologies. Raise the points of the .pss/.pac "
                        "line to check convergence.\n", ii->INDname);
                g_dcm_said = 1;
                return;
            }
        }
}

static void
pac_imp_note(const PSSsaltEvent *ev, int n_ev)
{
    int i, routed = 0;
    g_imp_on = 0; g_imp_nx = 0;
    if (n_ev <= 0 || getenv("PAC_NO_IMPULSE")) return;
    /* only a flip whose node is algebraic or fast (the jump went to an inductor branch row, or this is a continuation) leaves a
     * discontinuity: with a capacitor at the switch node everything but the switch current is continuous */
    for (i = 0; i < n_ev; i++)
        if (ev[i].dG == 0.0 || (ev[i].posRow > 0 && ev[i].posRow != ev[i].posNode) || (ev[i].negRow > 0 && ev[i].negRow != ev[i].negNode)) routed = 1;
    if (!routed) return;
    for (i = 0; i < n_ev; i++) {
        /* the same event function with the control terminals swapped (the two switches of a half bridge) is the same coefficient */
        if (i > 0 && ev[i].posCntrl == g_imp_nc && ev[i].negCntrl == g_imp_pc && ev[i].h_dot == -g_imp_hdot) {
            if (ev[i].dG != 0.0) {
                if (ev[i].posRow > 0 && ev[i].posRow == ev[i].posNode && !ev[i].posF) g_imp_x[g_imp_nx++] = ev[i].posRow;
                if (ev[i].negRow > 0 && ev[i].negRow == ev[i].negNode && !ev[i].negF) g_imp_x[g_imp_nx++] = ev[i].negRow;
            }
            continue;
        }
        if (i > 0 && (ev[i].posCntrl != g_imp_pc || ev[i].negCntrl != g_imp_nc || ev[i].h_dot != g_imp_hdot)) {
            static int warned = 0;
            if (!warned) {
                fprintf(stderr, "PAC-TD: note -- two different switches flip within one accepted step: the output impulse of that step "
                        "is not attributed (card #579).\n");
                warned = 1;
            }
            g_imp_on = 0;
            return;
        }
        g_imp_pc = ev[i].posCntrl; g_imp_nc = ev[i].negCntrl; g_imp_hdot = ev[i].h_dot; g_imp_on = 1;
        /* capacitor-node jump rows at the switch's own power nodes (native capacitor or charge jump): excluded */
        if (ev[i].dG == 0.0) continue;   /* a continuation or a far-capacitor row of the same flip, not the switch itself */
        if (ev[i].posRow > 0 && ev[i].posRow == ev[i].posNode && !ev[i].posF) g_imp_x[g_imp_nx++] = ev[i].posRow;
        if (ev[i].negRow > 0 && ev[i].negRow == ev[i].negNode && !ev[i].negF) g_imp_x[g_imp_nx++] = ev[i].negRow;
        /* a charge-jump event: the recursion carries the node's own transient from the jumped charge, every row follows it */
        if (ev[i].posQ || ev[i].negQ) { g_imp_on = 0; return; }
    }
}

static void
pac_imp_apply(CKTcircuit *ckt, int neq, int M, double t_e)
{
    int m, row, k;
    if (!g_imp_on || !g_pac_xls_have || g_pac_xls_n != neq || fabs(g_imp_hdot) < 1e-30) return;
    for (m = 0; m < M; m++) {
        PACc Sp = pac_get(g_pac_Cprev, g_imp_pc, M, m);
        PACc Sn = pac_get(g_pac_Cprev, g_imp_nc, M, m);
        PACc coeff = pc_mul(pc_scale(pc_sub(Sp, Sn), 1.0 / g_imp_hdot), pc_conjexp(g_pac_omega[m] * t_e));
        for (row = 1; row <= neq; row++) {
            double J = ckt->CKTrhsOld[row] - g_pac_xls_prev[row];
            size_t idx = (size_t) (row - 1) * M + m;
            int skip = 0;
            for (k = 0; k < g_imp_nx; k++) if (g_imp_x[k] == row) skip = 1;
            if (skip || J == 0.0) continue;
            g_pac_demod[idx] = pc_add(g_pac_demod[idx], pc_scale(coeff, J));
        }
    }
}

static void
pac_imp_roll(CKTcircuit *ckt, int neq)
{
    if (g_pac_xls_n != neq) {
        g_pac_xls_prev = TREALLOC(double, g_pac_xls_prev, neq + 1);
        g_pac_xls_n = g_pac_xls_prev ? neq : 0;
    }
    /* while a routed jump is being continued the reference stays where the last impulse was given (the steps in between are part
     * of the same event: their change is given with the next continuation event or its final flush) */
    if (g_pac_xls_prev && (!g_pac_xls_have || g_imp_on || !PSSsaltContActive())) {
        memcpy(g_pac_xls_prev, ckt->CKTrhsOld, (size_t) (neq + 1) * sizeof(double)); g_pac_xls_have = 1;
    }
}

/* prepass step: accumulate the set of columns of C (state-set discovery) and run the affinity check. */
static void
pac_prepass_step(CKTcircuit *ckt)
{
    int i;
    if (ckt->CKTtime <= 0.0 || g_pac_gen_bad || !g_pp_col)
        return;
#ifdef XSPICE
    if (g_hid_nhy > 0 && ckt->CKTstate0 && ckt->CKTstate1) {   /* hyst: a flip of the Boolean state is an event the recursion does not carry */
        int hi;
        for (hi = 0; hi < g_hid_nhy; hi++)
            if (*(int *) &ckt->CKTstate0[g_hid_hyslot[hi]] != *(int *) &ckt->CKTstate1[g_hid_hyslot[hi]]) { g_hid_flips++; break; }
    }
    if (g_hid_nc > 0)
        g_hid_mode = ((g_pp_steps % g_hid_stride) == 0) ? 1 : 0;
#endif
    if (pac_gen_loads(ckt, ckt->CKTag[0], 1, 0)) {
        g_pac_gen_bad = -1;   /* hard failure */
        return;
    }
    g_pp_steps++;
    if (g_salq) pac_salq_note(ckt);   /* card #559 */
#ifdef XSPICE
    if (g_hid_nc > 0 && getenv("PAC_HID_TRACE") && (g_pp_steps <= 2 || ckt->CKTtime >= 0.999999 * ckt->CKTfinalTime)) {
        int hi;
        fprintf(stderr, "PAC-HID: prepass step %d t=%.9g slots:", g_pp_steps, ckt->CKTtime);
        for (hi = 0; hi < g_hid_nc && hi < 12; hi++) fprintf(stderr, " %.6g", ckt->CKTstate0[g_hid_cslot[hi]]);
        fprintf(stderr, "\n");
    }
#endif
    if (g_pp_xend)
        memcpy(g_pp_xend, ckt->CKTrhs, (size_t) (g_pp_neq + 1) * sizeof(double));
    for (i = 0; i < g_pac_C.n; i++) {
        int c = g_pac_C.ci[i];
        if (c > 0 && c <= g_pp_neq) g_pp_col[c] = 1;
        if (fabs(g_pac_C.v[i]) > g_pp_cmax) g_pp_cmax = fabs(g_pac_C.v[i]);
    }
}

void
PSSpacRecordStep(CKTcircuit *ckt)
{
    int n, neq, M, d, m, row, order;
    double ag0, ag1, t_new, h;
    int is_first_call;

    if (g_pp_active) {   /* card #542: state-set discovery prepass */
        if (g_pp_detect)
            pac_detect_step(ckt);
        else
            pac_prepass_step(ckt);
        return;
    }
    if (g_pac_failed || !g_pac_map)
        return;
    if (ckt->CKTtime <= 0.0)
        return;   /* t=0 base case already set by pac_reset() */

    /* card #465 FOURTH pass, task 1: g_pac_tprev is set to exactly 0.0 by
     * pac_reset() and only ever reassigned to a real t_new at the tail of
     * this function, so this test identifies the very first accepted-step
     * call of this recording period run -- see pac_finish_demod_wraparound()
     * for why the demod quadrature's first segment is deferred rather than
     * accumulated here using g_pac_Cprev's (wrong, at non-state rows) t=0
     * seed value. */
    is_first_call = (g_pac_tprev == 0.0);

    if (ckt->CKTintegrateMethod == GEAR && g_pac_gear && g_pac_gen && !g_pac_traj_active) {
        /* card #549: handled by the generic recursion (charge history); the trap-order check below does not apply */
    } else if (ckt->CKTintegrateMethod == GEAR && g_pac_traj_active) {
        fprintf(stderr, "PAC-TD: method=gear is not supported by the adjoint/noise recursion (.pnoise, .pxf adjoint=1: the "
                "transposed sweep replays the trapezoidal derivative history); use method=trap (default). Abandoning this pass.\n");
        g_pac_failed = 1;
        return;
    } else if (ckt->CKTintegrateMethod != TRAPEZOIDAL || ckt->CKTorder < 1 || ckt->CKTorder > 2) {
        fprintf(stderr, "PAC-TD: analytic propagation only implements "
                "TRAPEZOIDAL orders 1-2 (got method=%d order=%d at t=%g); "
                "abandoning this pass.\n", ckt->CKTintegrateMethod, ckt->CKTorder, ckt->CKTtime);
        g_pac_failed = 1;
        return;
    }

    n = g_pac_n; neq = g_pac_neq; M = g_pac_M;
    ag0 = ckt->CKTag[0]; ag1 = ckt->CKTag[1];
    /* same order-detection-from-ag1 trick as PSSsensRecordStep (pssshoot.c) --
     * see that function's header comment for why ckt->CKTorder itself is
     * not trustworthy at this call site. */
    order = (fabs(ag1 + ag0) < 1e-6 * fabs(ag0)) ? 1 : 2;
    t_new = ckt->CKTtime;
    h = t_new - g_pac_tprev;
    if (getenv("PAC_STEP_TRACE"))   /* card #517 diagnostics */
        fprintf(stderr, "STEP517 t=%.17g h=%.17g\n", t_new, h);

    /* card #504: adjoint-engine trajectory recording, independent of the
     * M-column machinery below -- see g_pac_traj_active's own comment.
     * Uses ckt->CKTrhs (the true accepted x_k, untouched by anything in
     * this function), not g_pac_Cprev/Cnew. */
    if (g_pac_traj_active)
        pac_traj_append(ckt, t_new, h, ag0, ag1, order);
    /* card #465 THIRD pass, task 1 investigation: debug trace (PAC_ORDER_TRACE,
     * no-op unless set) used to rule out "order pinned at 1" as the cause of
     * pstb_linear.cir's O(1/Points) error -- it is NOT: on a 20000-point
     * recording pass only the first 6 of ~80000 calls (t < 2e-11 s, well
     * before the first real breakpoint) are order 1, the rest are trapezoidal
     * order 2 throughout. Left in since it is cheap and the question comes
     * up again with every accuracy pass on this file. */
    if (getenv("PAC_ORDER_TRACE")) {
        static int calls = 0, n1 = 0, n2 = 0;
        if (order == 1) n1++; else n2++;
        calls++;
        if (calls <= 5 || calls % 500 == 0)
            fprintf(stderr, "PACORDER call=%d t=%.9g h=%.6g ag0=%.6g ag1=%.6g order=%d n1=%d n2=%d\n",
                    calls, t_new, h, ag0, ag1, order, n1, n2);
    }


    /* card #465 THIRD pass, task 1: refactor ckt->CKTmatrix at the TRUE
     * accepted point x_k before using it below -- see
     * pac_refactor_at_accepted()'s header for why this is needed and why
     * it is safe to do here (unlike pssshoot.c's Newton-loop-scoped
     * case). No-op unless PAC_EXACT_JACOBIAN is set. */
    if (!g_pac_gen && (g_pac_exact_jac || g_pac_force_exact_jac)) {
        if (pac_refactor_at_accepted(ckt)) {
            fprintf(stderr, "PAC-TD: exact-Jacobian refactor failed at "
                    "t=%g; abandoning this pass.\n", ckt->CKTtime);
            g_pac_failed = 1;
            return;
        }
    }

    {
    /* card #494 (plan T6.1): locate any sw event at THIS accepted step
     * once, reuse for every column m below -- see pssshoot.c's
     * PSSsaltFindEvents()/PSSsensRecordStep() for the shared derivation
     * (real there, complex here: same formula, coeff is now a PACc). */
    PSSsaltEvent salt_ev[PSS_SALT_MAX_EVENTS];
    PSSsaltNoRoutedWindow = (g_pac_gen && pac_smooth_t()) ? 1 : 0;   /* card #583, before the events of this step are formed */
    int n_salt = PSSsaltFindEvents(ckt, salt_ev, PSS_SALT_MAX_EVENTS);
    int i_salt;

    /* card #510: PSSsaltFindEvents() already printed why; REFUSE this
     * whole PAC-TD pass rather than silently continue missing an event's
     * forcing term. */
    if (n_salt < 0) {
        g_pac_failed = 1;
        return;
    }
    pac_imp_note(salt_ev, n_salt);   /* card #579 */
    pac_dcm_note(ckt, n_salt);       /* card #581 */

    /* card #494 fifth pass (task 2): same L-stable-settle-window override
     * as PSSsensRecordStep() (pssshoot.c) -- see PSSsaltForceOrder1()'s
     * header comment in pssshoot.h. Must run after PSSsaltFindEvents()
     * above and before `order` is used below; no-op for any run that
     * never saw a real sw event. */
    if (PSSsaltForceOrder1(ckt))
        order = 1;

    if (g_pac_gen) {
        /* card #542/#544: generic charge-based recursion; saltation events are applied inside it (guarded). */
        if (n_salt > 0) {
            static int warned_salt = 0;
            if (!warned_salt && getenv("PAC_GEN_TRACE")) {
                fprintf(stderr, "PAC-TD: saltation event(s) on the generic path (card #544), first at t=%g\n", t_new);
                warned_salt = 1;
            }
        }
        if (pac_step_generic(ckt, is_first_call, ag0, ag1, order, t_new, salt_ev, n_salt)) {
            g_pac_failed = 1;
            return;
        }
        if (g_pac_traj_active)
            pac_traj_attach_generic();
    } else {
    for (m = 0; m < M; m++) {
        memset(g_pac_rhs_re, 0, (size_t) (neq + 1) * sizeof(double));
        memset(g_pac_rhs_im, 0, (size_t) (neq + 1) * sizeof(double));

        for (d = 0; d < n; d++) {
            PACstateRef *ref = &g_pac_map->refs[d];
            PACc Tp = g_pac_Tprev[(size_t) d * M + m];
            if (ref->isCap) {
                CAPinstance *ci = ref->dev.cap;
                int p = ci->CAPposNode, ngc = ci->CAPnegNode;
                double C = ci->CAPcapac;
                PACc Sp = pac_get(g_pac_Cprev, p, M, m);
                PACc Sn = pac_get(g_pac_Cprev, ngc, M, m);
                PACc diff = pc_sub(Sp, Sn);
                PACc dceq = pc_scale(diff, -ag0 * C);
                if (order == 2) dceq = pc_add(dceq, pc_scale(Tp, -ag1));
                if (p > 0)   { g_pac_rhs_re[p] -= dceq.re; g_pac_rhs_im[p] -= dceq.im; }
                if (ngc > 0) { g_pac_rhs_re[ngc] += dceq.re; g_pac_rhs_im[ngc] += dceq.im; }
            } else {
                INDinstance *ii = ref->dev.ind;
                int ibr = ii->INDbrEq;
                double L = ii->INDinduct;
                PACc Si = pac_get(g_pac_Cprev, ibr, M, m);
                PACc dveq = pc_scale(Si, -ag0 * L);
                if (order == 2) dveq = pc_add(dveq, pc_scale(Tp, -ag1));
                g_pac_rhs_re[ibr] += dveq.re; g_pac_rhs_im[ibr] += dveq.im;
            }
        }

        if (g_pac_useforce[m] && g_pac_br && g_pac_bi) {
            double w = g_pac_omega[m];
            double ph = w * t_new;
            double c = cos(ph), s = sin(ph);
            for (row = 1; row <= neq; row++) {
                double br = g_pac_br[row], bi = g_pac_bi[row];
                if (br == 0.0 && bi == 0.0) continue;
                g_pac_rhs_re[row] += br * c - bi * s;
                g_pac_rhs_im[row] += br * s + bi * c;
            }
        }

        if (pac_mq_on() && !pac_hid_active())
            pac_mq_put(m);   /* card #558: solved for all columns after the loop */
        else {
            pac_solve_cplx(ckt);
            for (row = 1; row <= neq; row++) {
                PACc v; v.re = g_pac_rhs_re[row]; v.im = g_pac_rhs_im[row];
                pac_set(g_pac_Cnew, row, M, m, v);
            }
        }
    }
    if (pac_mq_on() && !pac_hid_active())
        pac_mq_flush(ckt, M);

    /* card #494 THIRD pass (rewrite -- see pssshoot.c's PSSsaltFindEvents()
     * derivation comment): apply the saltation correction as a DIRECT VALUE
     * JUMP added onto g_pac_Cnew, AFTER this step's own smooth dceq/dveq
     * solve above, NOT routed through it. Routing the jump through this
     * step's own SMPsolve (tried earlier this pass, before this rewrite)
     * measurably attenuated it whenever the just-flipped switch's own new
     * conductance is comparable to ag0*C_local at that node (exactly the
     * ~9x-not-negligible case the SECOND pass's report already flagged) --
     * the matrix divides the jump by (ag0*C + G_new + ...) instead of
     * leaving it as the pure state discontinuity Fang's formula calls for.
     * Adding it directly to the already-solved g_pac_Cnew avoids that
     * entirely. coeff uses g_pac_Cprev (the UNMODIFIED value from the END
     * of the PREVIOUS accepted step -- delta_x(t_e-), exactly what the
     * formula calls for) for every event, so simultaneous events (e.g. two
     * complementary switches) never see a partially-applied correction. */
    if (n_salt > 0) {
        for (m = 0; m < M; m++) {
            PACc djump[PSS_SALT_MAX_EVENTS][2];   /* [i_salt][0]=posRow, [1]=negRow */
            for (i_salt = 0; i_salt < n_salt; i_salt++) {
                PSSsaltEvent *e = &salt_ev[i_salt];
                PACc Sp = pac_get(g_pac_Cprev, e->posCntrl, M, m);
                PACc Sn = pac_get(g_pac_Cprev, e->negCntrl, M, m);
                PACc coeff = pc_scale(pc_sub(Sp, Sn), 1.0 / e->h_dot);
                djump[i_salt][0] = pc_scale(coeff, e->dfdt_pos);
                djump[i_salt][1] = pc_scale(coeff, e->dfdt_neg);
            }
            for (i_salt = 0; i_salt < n_salt; i_salt++) {
                PSSsaltEvent *e = &salt_ev[i_salt];
                /* card #510: posRow/negRow -- posNode/negNode for the
                 * capacitor case, or an adjacent inductor's INDbrEq row
                 * for the new inductor-routed case, or 0/skip; see
                 * pssshoot.h's PSSsaltEvent comment. */
                if (e->posRow > 0)
                    pac_set(g_pac_Cnew, e->posRow, M, m,
                            pc_add(pac_get(g_pac_Cnew, e->posRow, M, m), djump[i_salt][0]));
                if (e->negRow > 0)
                    pac_set(g_pac_Cnew, e->negRow, M, m,
                            pc_add(pac_get(g_pac_Cnew, e->negRow, M, m), djump[i_salt][1]));
            }
        }
    }
    }   /* card #542: closes the legacy native branch */
    }   /* card #494: closes the salt_ev/n_salt scope opened above */

    for (d = 0; d < (g_pac_gen ? 0 : n); d++) {
        PACstateRef *ref = &g_pac_map->refs[d];
        for (m = 0; m < M; m++) {
            PACc val;
            if (ref->isCap) {
                CAPinstance *ci = ref->dev.cap;
                int p = ci->CAPposNode, ngc = ci->CAPnegNode;
                double C = ci->CAPcapac;
                PACc dn = pc_sub(pac_get(g_pac_Cnew, p, M, m), pac_get(g_pac_Cnew, ngc, M, m));
                PACc dop = pc_sub(pac_get(g_pac_Cprev, p, M, m), pac_get(g_pac_Cprev, ngc, M, m));
                val = pc_scale(pc_sub(dn, dop), ag0 * C);
            } else {
                INDinstance *ii = ref->dev.ind;
                int ibr = ii->INDbrEq;
                double L = ii->INDinduct;
                PACc dn = pac_get(g_pac_Cnew, ibr, M, m);
                PACc dop = pac_get(g_pac_Cprev, ibr, M, m);
                val = pc_scale(pc_sub(dn, dop), ag0 * L);
            }
            if (order == 2)
                val = pc_add(val, pc_scale(g_pac_Tprev[(size_t) d * M + m], -ag1));
            g_pac_Tnew[(size_t) d * M + m] = val;
        }
    }

    if (g_pac_demod_active && !is_first_call)
        pac_imp_apply(ckt, neq, M, t_new);   /* card #579: before Cprev is swapped (the event coefficient reads it) */
    pac_imp_roll(ckt, neq);

    if (g_pac_demod_active) {
        if (is_first_call) {
            /* card #465 FOURTH pass, task 1: y(0) at a non-state row is NOT
             * known yet (g_pac_Cprev's t=0 seed is only correct at the n
             * reduced-state rows -- see pac_reset()) -- stash (t_1, y_1) and
             * let pac_finish_demod_wraparound() supply the correct [0,t_1]
             * segment, exactly, once y(T) closes the loop (periodicity). */
            for (m = 0; m < M; m++)
                for (row = 1; row <= neq; row++)
                    pac_set(g_pac_first_y, row, M, m, pac_get(g_pac_Cnew, row, M, m));
            g_pac_first_t = t_new;
            g_pac_have_first = 1;
        } else {
            for (m = 0; m < M; m++) {
                double w = g_pac_omega[m];
                PACc ph_old = pc_conjexp(w * g_pac_tprev);
                PACc ph_new = pc_conjexp(w * t_new);
                for (row = 1; row <= neq; row++) {
                    PACc yold = pac_get(g_pac_Cprev, row, M, m);
                    PACc ynew = pac_get(g_pac_Cnew, row, M, m);
                    PACc term = pc_scale(pc_add(pc_mul(yold, ph_old), pc_mul(ynew, ph_new)), 0.5 * h);
                    size_t idx = (size_t) (row - 1) * M + m;
                    g_pac_demod[idx] = pc_add(g_pac_demod[idx], term);
                }
            }
        }
    }

    { PACc *tmp = g_pac_Cprev; g_pac_Cprev = g_pac_Cnew; g_pac_Cnew = tmp; }
    { PACc *tmp = g_pac_Tprev; g_pac_Tprev = g_pac_Tnew; g_pac_Tnew = tmp; }
    g_pac_tprev = t_new;
}

/* ------------------------------------------------------------------ */
/* card #465 FOURTH pass, task 1: DISCRETISATION FIX for the pac_sb0(d)
 * .pstb-vs-native-.stb O(1/Points) [NOT O(h^2)] error.
 *
 * ROOT CAUSE (traced this pass, replacing the THIRD pass's
 * PAC_EXACT_JACOBIAN Jacobian-staleness fix attempt, which did not touch
 * this error -- it is a linear circuit, PSSsensRecordStep's/this file's
 * own J_k staleness only matters for history-dependent nonlinear devices):
 * pac_reset() seeds g_pac_Cprev's t=0 row for the demodulation pass with
 * the TRUE y(0) only at the n reduced-state rows (the caller-supplied
 * ic0[], x0 or xi -- a capacitor/inductor terminal IS one of those rows,
 * pac_state_row()). Every OTHER row (a plain circuit node or a branch
 * equation, e.g. pac_td_pstb()'s branch_probe output, OR pac_td_pac()'s
 * many non-capacitor/inductor node outputs) is left at pac_reset()'s
 * memset() zero. That zero is WRONG in general -- a non-state row's true
 * value at t=0 is whatever the algebraic KCL constraint gives, given the
 * (correct) state ICs and the forcing b(0) -- and PSSpacRecordStep()'s
 * demodulation quadrature (this file's header, step 3) used it anyway for
 * the FIRST trapezoid segment [0, t_1], the ONLY segment that reads
 * g_pac_Cprev's un-propagated seed. Every OTHER segment [t_{k-1}, t_k]
 * (k=2..N) uses two values that both came out of a genuine SMPsolve() at
 * an accepted step -- correct for ALL rows, state and non-state alike --
 * so only this one segment was ever wrong. A wrong O(1) value multiplying
 * 0.5*h_1 in a sum normalised by 1/T contributes an O(h_1/T) = O(1/Points)
 * error to Y0 -- exactly the empirically observed scaling (test_pac.py's
 * "pstb_loop" header, 8.3e-4 at Points=2000 down to 5.77e-7 only by
 * Points=3,000,000, textbook 1/Points).
 *
 * DERIVATION of the exact fix (no extra algebraic solve needed): this
 * file's header already establishes that re-running the SAME period from
 * the SAME converged x0 reproduces the transient's adaptive-step decisions
 * bit for bit. Consider "period rep 2", i.e. relabelling t -> t+T and
 * continuing: because b(t+T) = b*e^{jw(t+T)} = e^{jwT} b(t) and the state
 * seed at the start of rep 2 is x(T) = e^{jwT} x0 (exactly, since xi was
 * solved for precisely that Floquet boundary condition), and because
 * PSSpacRecordStep()'s per-step map from (previous state row values,
 * forcing at t_k) to y_k (ALL rows, via ONE full SMPsolve()) is LINEAR and
 * uses the IDENTICAL companion coefficients/Jacobian/time grid at
 * "relative step k of rep 2" as at "relative step k of rep 1" (same
 * physical point in the periodic orbit) -- by linearity y_k^{rep2} =
 * e^{jwT} y_k^{rep1} for EVERY row, not just the n state rows. Taking
 * k = N (the last accepted step, t=T): y(T) computed by THIS run IS
 * "y_N^{rep1}", and by the SAME relation applied one step earlier,
 * "y_0^{rep2}" (= the true y(0) this file needs, continuing seamlessly
 * from where rep 1 left off) satisfies y(T) = y_N^{rep1} = e^{jwT} *
 * y_0^{rep2's own analogue, i.e. THIS run's true y(0)}. Hence:
 *
 *     y(0) = e^{-jwT} * y(T)                                   (EXACT)
 *
 * for every row, state or not -- sanity check: at the n state rows this
 * reproduces y(0)=xi exactly, since y(T)[state] = e^{jwT} xi BY
 * CONSTRUCTION (the boundary solve in pac_engine_run), so this fix leaves
 * the (already-correct) state-row t=0 values completely unchanged and only
 * corrects the non-state rows.
 *
 * y(T) is exactly g_pac_Cprev after PSSpacRecordStep()'s last call's final
 * swap (the period run has just returned) -- no extra CKTload()/SMPsolve()
 * needed. This function must run AFTER pac_run_period() returns success
 * and BEFORE pac_free()/the next pac_reset() overwrites g_pac_Cprev. */
static void
pac_finish_demod_wraparound(void)
{
    int m, row;

    if (!g_pac_demod_active || !g_pac_demod || !g_pac_have_first || !g_pac_first_y)
        return;

    for (m = 0; m < g_pac_M; m++) {
        double w = g_pac_omega[m];
        PACc ph1 = pc_conjexp(w * g_pac_first_t);       /* e^{-j*w*t_1} */
        PACc ephaseT = pc_conjexp(w * g_pac_T);          /* e^{-j*w*T} */
        for (row = 1; row <= g_pac_neq; row++) {
            PACc yT = pac_get(g_pac_Cprev, row, g_pac_M, m);      /* y(T), all rows exact */
            PACc y0 = pc_mul(yT, ephaseT);                         /* y(0) = e^{-jwT} y(T) */
            PACc y1 = pac_get(g_pac_first_y, row, g_pac_M, m);
            PACc term = pc_scale(pc_add(y0, pc_mul(y1, ph1)), 0.5 * g_pac_first_t);
            size_t idx = (size_t) (row - 1) * g_pac_M + m;
            g_pac_demod[idx] = pc_add(g_pac_demod[idx], term);
        }
    }
}

/* ------------------------------------------------------------------ */
/* Small dense complex Gauss elimination with partial pivoting -- n is
 * always small here (1-2 on a switched RC, ~7 on R3/R4). */
static int
pac_complex_solve(PACc *A, PACc *rhs, int n, PACc *x)
{
    int i, j, k, piv;
    double best;

    for (k = 0; k < n; k++) {
        piv = k; best = pc_mag(A[k * n + k]);
        for (i = k + 1; i < n; i++) {
            double m = pc_mag(A[i * n + k]);
            if (m > best) { best = m; piv = i; }
        }
        if (best < 1e-300)
            return 1;   /* singular */
        if (piv != k) {
            for (j = 0; j < n; j++) { PACc t = A[k * n + j]; A[k * n + j] = A[piv * n + j]; A[piv * n + j] = t; }
            { PACc t = rhs[k]; rhs[k] = rhs[piv]; rhs[piv] = t; }
        }
        for (i = k + 1; i < n; i++) {
            PACc f = pc_div(A[i * n + k], A[k * n + k]);
            for (j = k; j < n; j++)
                A[i * n + j] = pc_sub(A[i * n + j], pc_mul(f, A[k * n + j]));
            rhs[i] = pc_sub(rhs[i], pc_mul(f, rhs[k]));
        }
    }
    for (i = n - 1; i >= 0; i--) {
        PACc s = rhs[i];
        for (j = i + 1; j < n; j++)
            s = pc_sub(s, pc_mul(A[i * n + j], x[j]));
        x[i] = pc_div(s, A[i * n + i]);
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Capture the netlist AC-source stimulus B0 (SAME idiom as dcpss.c's
 * pac_extract_harmonics(), Enhancement-123: CKTload(DCOP|INITSMSIG) at
 * the CURRENT ckt state, then CKTacLoad(w=1) with the matrix in complex
 * mode -- the matrix content this clobbers is immediately superseded by
 * the fresh DCtran() the very next recording period run performs, so no
 * separate save/restore is needed). br/bi: [neq+1], 1-based. */
/* card #533: under `.options klu` the complex stamp mode this function switches
 * the matrix into (DEVbindCSCComplex: every device's stamp pointers re-bound
 * to the interleaved complex CSC values) must be undone before returning: the
 * recording/shooting runs that follow are REAL transients and would otherwise
 * stamp through complex pointers into a real Ax (singular matrix at once).
 * Same conversion acan.c does at the end of .ac (DEVbindCSCComplexToReal,
 * KLUmatrixIsComplex = KLUmatrixReal). SPARSE path unchanged (byte-identical:
 * nothing is compiled in or executed for it). */
#ifdef KLU
static void
pac_klu_back_to_real(CKTcircuit *ckt)
{
    int i;
    if (ckt->CKTmatrix->CKTkluMODE &&
        ckt->CKTmatrix->SMPkluMatrix->KLUmatrixIsComplex) {
        for (i = 0; i < DEVmaxnum; i++)
            if (DEVices[i] && DEVices[i]->DEVbindCSCComplexToReal && ckt->CKThead[i])
                DEVices[i]->DEVbindCSCComplexToReal(ckt->CKThead[i], ckt);
        ckt->CKTmatrix->SMPkluMatrix->KLUmatrixIsComplex = KLUmatrixReal;
    }
}
#endif

static int
pac_capture_b0(CKTcircuit *ckt, double *br, double *bi, int neq)
{
    int i;
    long savedMode = ckt->CKTmode;
    int err = 0;
#ifdef XSPICE
    Mif_Analysis_t saved_anal;
#endif

#ifdef KLU
    if (ckt->CKTmatrix->CKTkluMODE) {
        if (!ckt->CKTmatrix->SMPkluMatrix->KLUmatrixIsComplex) {
            for (i = 0; i < DEVmaxnum; i++)
                if (DEVices[i] && DEVices[i]->DEVbindCSCComplex && ckt->CKThead[i])
                    DEVices[i]->DEVbindCSCComplex(ckt->CKThead[i], ckt);
            ckt->CKTmatrix->SMPkluMatrix->KLUmatrixIsComplex = KLUMatrixComplex;
        }
    } else
#endif
        spSetComplex(ckt->CKTmatrix->SPmatrix);

    ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODEDCOP | MODEINITSMSIG;
    if (CKTload(ckt)) err = 1;
    if (!err) {
        ckt->CKTomega = 1.0;
        ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODEAC;
#ifdef XSPICE
        /* card #551: MIFload (the XSPICE DEVacLoad) takes its analysis type from g_mif_info, which only CKTload/acan set;
         * after the CKTload above it still says DC, so every analog code model stamped its DC RHS (offset, output -
         * partial*input) into the "small-signal" source vector and the td PAC response carried a constant error. */
        saved_anal = g_mif_info.circuit.anal_type;
        g_mif_info.circuit.anal_type = MIF_AC;
#endif
        if (CKTacLoad(ckt)) err = 1;
#ifdef XSPICE
        g_mif_info.circuit.anal_type = saved_anal;
#endif
    }
    if (!err)
        for (i = 1; i <= neq; i++) {
            br[i] = ckt->CKTrhs[i];
            bi[i] = ckt->CKTirhs[i];
        }
    ckt->CKTmode = savedMode;
#ifdef KLU
    pac_klu_back_to_real(ckt);
#endif
    return err;
}

/* ------------------------------------------------------------------ */
/* card #499 (plan T5.4): `.pxf ... method="shoot" solver="td"` -- one
 * time-domain shooting-based transfer function PER independent V/I source,
 * to a single fixed output, reusing pac_engine_run() (and therefore every
 * already-validated per-step recursion/boundary-solve/demodulation line in
 * this file) completely UNCHANGED -- see this file's own "SHIPPED DESIGN"
 * note above pac_td_pxf() below for why a batched-forward-column sweep was
 * chosen over a from-scratch backward/adjoint sweep of the per-step
 * recursion (the literal reading of "one adjoint solve").
 *
 * pac_build_src_list() walks every native Vsource/Isource instance in the
 * circuit (independent of whether the netlist itself gave it an AC value --
 * .pxf's whole point, confirmed against Spectre's own pxf semantics, is
 * NOT needing a deck's own AC tags) and returns a pointer to each one's OWN
 * VSRCacReal/VSRCacImag (or ISRCacReal/ISRCacImag) pair -- the field
 * VSRCacLoad()/ISRCacLoad() actually stamp into the RHS (computed from
 * VSRCacMag/VSRCacPhase at VSRCtemp()/ISRCtemp() time, card #499 found;
 * poking ACmag/ACphase directly would need a fresh temp pass to take
 * effect, poking *Real/*Imag directly does not). */
typedef struct {
    GENinstance *inst;
    double *preal, *pimag;
    IFuid name;
} PACsrcRef;

static int
pac_build_src_list(CKTcircuit *ckt, PACsrcRef **out)
{
    int vtype = CKTtypelook("Vsource");
    int itype = CKTtypelook("Isource");
    int n = 0, cap = 0;
    PACsrcRef *list = NULL;
    VSRCmodel *vm; VSRCinstance *vi;
    ISRCmodel *im; ISRCinstance *ii;

    if (vtype >= 0)
        for (vm = (VSRCmodel *) ckt->CKThead[vtype]; vm; vm = VSRCnextModel(vm))
            for (vi = VSRCinstances(vm); vi; vi = VSRCnextInstance(vi)) {
                if (n >= cap) { cap = cap ? cap * 2 : 8; list = TREALLOC(PACsrcRef, list, cap); }
                list[n].inst = (GENinstance *) vi;
                list[n].preal = &vi->VSRCacReal;
                list[n].pimag = &vi->VSRCacImag;
                list[n].name = vi->VSRCname;
                n++;
            }
    if (itype >= 0)
        for (im = (ISRCmodel *) ckt->CKThead[itype]; im; im = ISRCnextModel(im))
            for (ii = ISRCinstances(im); ii; ii = ISRCnextInstance(ii)) {
                if (n >= cap) { cap = cap ? cap * 2 : 8; list = TREALLOC(PACsrcRef, list, cap); }
                list[n].inst = (GENinstance *) ii;
                list[n].preal = &ii->ISRCacReal;
                list[n].pimag = &ii->ISRCacImag;
                list[n].name = ii->ISRCname;
                n++;
            }
    *out = list;
    return n;
}

/* Isolate ONE source (list[idx]) as if it alone had "AC 1" in the deck --
 * zero every OTHER source's *Real/*Imag first (so a deck with its own
 * unrelated AC-flagged source, or the "lgprobe hidden AC 1" trap
 * docs/pss-pac.md documents for .pac, cannot leak into a .pxf per-source
 * reading), restore all of them (regardless of what pac_capture_b0()
 * itself did/didn't touch) on the way out. Reuses pac_capture_b0()'s own
 * CKTload(DCOP|INITSMSIG)+CKTacLoad(w=1) idiom verbatim -- guarantees
 * byte-identical stamping to what an ordinary `.pac`/`AC 1` run against
 * the SAME single source would produce (the done-when this card's report
 * validates against). */
static int
pac_capture_source_b0(CKTcircuit *ckt, PACsrcRef *list, int nsrc, int idx,
                       double *br, double *bi, int neq)
{
    double *saveR, *saveI;
    int i, err;

    if (idx < 0 || idx >= nsrc)
        return 1;
    saveR = TMALLOC(double, nsrc);
    saveI = TMALLOC(double, nsrc);
    if (!saveR || !saveI) {
        if (saveR) tfree(saveR);
        if (saveI) tfree(saveI);
        return 1;
    }
    for (i = 0; i < nsrc; i++) {
        saveR[i] = *list[i].preal;
        saveI[i] = *list[i].pimag;
        *list[i].preal = (i == idx) ? 1.0 : 0.0;
        *list[i].pimag = 0.0;
    }

    err = pac_capture_b0(ckt, br, bi, neq);

    for (i = 0; i < nsrc; i++) {
        *list[i].preal = saveR[i];
        *list[i].pimag = saveI[i];
    }
    tfree(saveR);
    tfree(saveI);
    return err;
}

/* ------------------------------------------------------------------ */
/* Frequency-sweep enumeration (SAME stepType convention as the existing
 * .pac card: 0=LIN total points, 1=DEC, 2=OCT points-per-decade/octave). */
static int
pac_sweep_count(int stepType, int np, double fstart, double fstop)
{
    int cnt = 0;
    double freq, mult;
    if (fstart <= 0.0 || fstop < fstart || np < 1) return 0;
    mult = (stepType == 1) ? pow(10.0, 1.0 / np) : (stepType == 2) ? pow(2.0, 1.0 / np) : 0.0;
    for (freq = fstart; freq <= fstop * (1.0 + 1e-9); ) {
        cnt++;
        if (stepType == 0) { if (np <= 1) break; freq += (np > 1) ? (fstop - fstart) / (np - 1) : 0.0; }
        else freq *= mult;
        if (cnt > 100000) break;   /* runaway guard */
    }
    return cnt;
}
static void
pac_sweep_fill(int stepType, int np, double fstart, double fstop, double *freqs, int cnt)
{
    double freq, mult;
    int i = 0;
    mult = (stepType == 1) ? pow(10.0, 1.0 / np) : (stepType == 2) ? pow(2.0, 1.0 / np) : 0.0;
    for (freq = fstart; i < cnt; ) {
        freqs[i++] = freq;
        if (stepType == 0) { if (np <= 1) break; freq += (np > 1) ? (fstop - fstart) / (np - 1) : 0.0; }
        else freq *= mult;
    }
}


/* ------------------------------------------------------------------ */
/* card #542: decide native vs generic history and, for generic, find the state set.
 *
 * Native-only circuits (every storage element is a C or L instance and nothing else carries
 * charge) keep the legacy per-C/L recursion bit-for-bit.  Otherwise (any device type outside the
 * memoryless/native list, e.g. MOS, BJT, diode, JFET; or a mutual inductor) the generic recursion
 * is used.  PAC_FORCE_GENERIC=1 forces it on native circuits (cross-check), PAC_NATIVE_HISTORY=1
 * forces the legacy history on any circuit (red proof: it reproduces the card #539 bug).
 * Returns 0 (ok, map->generic set or not), or an SPICE error code to return (REFUSED). */
/* card #582: is there an ideal switch with a power node that has neither a native capacitor nor a voltage source on it? */
static int
pac_switch_node_without_cap(CKTcircuit *ckt)
{
    int st = CKTtypelook("Switch"), ct = CKTtypelook("Capacitor"), vt = CKTtypelook("Vsource"), k;
    SWmodel *sm; SWinstance *si;
    if (st < 0) return 0;
    /* card #585: every circuit with an ideal switch.  The legacy per-C/L recursion is 12.4 dB / 94 deg off on a half bridge whose
     * switch node has neither capacitor nor diode (same events and orbit as the deck with a chargeless diode there, which is right; not
     * a stale matrix: reloading and refactoring before the legacy solve changes nothing; cause not found), and nothing ties that
     * defect to the node being capacitor-less.  The generic recursion agrees with it to 1e-6 where the legacy one is validated and is
     * right where it is not.  PAC_LEGACY_SW=1 keeps the narrower rule of cards #582 / #583 below (comparison). */
    if (!getenv("PAC_LEGACY_SW") && ckt->CKThead[st]) return 1;
    for (sm = (SWmodel *) ckt->CKThead[st]; sm; sm = SWnextModel(sm))
        for (si = SWinstances(sm); si; si = SWnextInstance(si))
            for (k = 0; k < 2; k++) {
                int nd = k ? si->SWnegNode : si->SWposNode, held = 0;
                if (nd <= 0) continue;
                if (ct >= 0) {
                    CAPmodel *cm; CAPinstance *ci;
                    for (cm = (CAPmodel *) ckt->CKThead[ct]; cm && !held; cm = CAPnextModel(cm))
                        for (ci = CAPinstances(cm); ci && !held; ci = CAPnextInstance(ci))
                            if (ci->CAPposNode == nd || ci->CAPnegNode == nd) held = 1;
                }
                if (vt >= 0 && !held) {
                    VSRCmodel *vm; VSRCinstance *vi;
                    for (vm = (VSRCmodel *) ckt->CKThead[vt]; vm && !held; vm = VSRCnextModel(vm))
                        for (vi = VSRCinstances(vm); vi && !held; vi = VSRCnextInstance(vi))
                            if (vi->VSRCposNode == nd || vi->VSRCnegNode == nd) held = 1;
                }
                if (!held) return 1;
                {   /* card #583: a diode at the switch node: the node may be treated as a fast node (boost-type, #580), whose routed
                     * jumps are added after T in the generic recursion only */
                    int dt = CKTtypelook("Diode");
                    DIOmodel *dm; DIOinstance *di;
                    if (dt >= 0)
                        for (dm = (DIOmodel *) ckt->CKThead[dt]; dm; dm = DIOnextModel(dm))
                            for (di = DIOinstances(dm); di; di = DIOnextInstance(di))
                                if (di->DIOposNode == nd || di->DIOnegNode == nd) return 1;
                }
            }
    return 0;
}

static int
pac_circuit_has_model_charges(CKTcircuit *ckt, const char **delay_dev)
{
    static const char *const nocharge[] = { "Resistor", "Capacitor", "Inductor", "Vsource",
        "Isource", "VCVS", "VCCS", "CCVS", "CCCS", "Switch", "CSwitch", "ASRC", NULL };
    int i, k, has = 0;
    *delay_dev = NULL;
    g_pp_nsusp = 0; g_pp_has_mutual = 0;
    for (i = 0; i < DEVmaxnum; i++) {
        const char *nm;
        int known = 0;
        if (!DEVices[i] || !ckt->CKThead[i] || !DEVices[i]->DEVpublic.name)
            continue;
        nm = DEVices[i]->DEVpublic.name;
        if (!strcmp(nm, "Tranline") || !strcmp(nm, "TransLine") || !strcmp(nm, "CplLines") ||
            !strcmp(nm, "LTRA")) {
            *delay_dev = nm;
            return 1;
        }
        for (k = 0; nocharge[k]; k++)
            if (!strcmp(nm, nocharge[k])) { known = 1; break; }
        if (known)
            continue;
        has = 1;
        if (!strcmp(nm, "mutual"))
            g_pp_has_mutual = 1;
        else if (DEVices[i]->DEVload && g_pp_nsusp < 64)
            g_pp_susp[g_pp_nsusp++] = i;
    }
    return has;
}

/* Non-quasi-static MOS models.  BSIM4 (trnqsmod=1) and BSIM3v3.3 (nqsmod=1) carry the channel-charge relaxation in an extra circuit
 * node (qNode, "charge deficit"), stamped with charge terms proportional to ag0 (the affine-in-ag0 guard of the generic path passes) and
 * conductive terms, so its relaxation is an ordinary state of the generic charge recursion; since card #559 it is accepted:
 * td .pac vs a transient FRA of the same trnqsmod=1 model (LTI bias, 1..30 MHz, wt = 1 at 10 MHz): 1e-4 relative or better, 1e-5 typical
 * (BSIM4 and BSIM3); on a 1 MHz LO-switched BSIM4 stage 0.003 dB / 0.01 deg vs the transient FRA (docs/pss-pac.md).  `.ac` is NOT the
 * reference: its acnqsmod is a different (single-pole) approximation, 20-40 % away from the transient model at omega*tau ~ 1.
 * acnqsmod alone (trnqsmod=0) does not change the transient model, hence not a td analysis: accepted with a note.  An NQS-capable device
 * of any OTHER family is still refused (not validated), and so is every NQS model in `.pnoise` (strict: the noise of the relaxation node is
 * not validated).  PAC_NQS_REFUSE=1 restores the pre-#559 refusal of every NQS model (red proof).
 * Returns the refusal text (static) or NULL; *note is set to a text to print once when an NQS model is accepted. */
static const char *
pac_find_nqs_device(CKTcircuit *ckt, const char **note, int strict)
{
    int i, k;
    static char buf[96], nbuf[160];
    *note = NULL;
    for (i = 0; i < DEVmaxnum; i++) {
        SPICEdev *d = DEVices[i];
        GENmodel *m;
        const char *dn;
        if (!d || !ckt->CKThead[i] || !d->DEVmodAsk || !d->DEVpublic.modelParms || !d->DEVpublic.numModelParms)
            continue;
        dn = d->DEVpublic.name ? d->DEVpublic.name : "?";
        for (k = 0; k < *d->DEVpublic.numModelParms; k++) {
            const char *kw = d->DEVpublic.modelParms[k].keyword;
            int is_tr, is_ac;
            if (!kw)
                continue;
            is_tr = !strcasecmp(kw, "trnqsmod") || !strcasecmp(kw, "nqsmod");
            is_ac = !strcasecmp(kw, "acnqsmod");
            if (!is_tr && !is_ac)
                continue;
            for (m = ckt->CKThead[i]; m; m = m->GENnextModel) {
                IFvalue v;
                memset(&v, 0, sizeof v);
                if (d->DEVmodAsk(ckt, m, d->DEVpublic.modelParms[k].id, &v) != OK || v.iValue == 0)
                    continue;
                if (getenv("PAC_NQS_REFUSE") || (is_tr && (strict || (strcmp(dn, "BSIM4") && strcmp(dn, "BSIM3"))))) {
                    snprintf(buf, sizeof buf, "%s (%s=%d)", dn, kw, v.iValue);
                    return buf;
                }
                if (is_tr) {
                    snprintf(nbuf, sizeof nbuf, "%s %s=%d: the channel-charge relaxation node is an ordinary state of the generic charge recursion (card #559)", dn, kw, v.iValue);
                    *note = nbuf;
                } else if (!*note) {
                    snprintf(nbuf, sizeof nbuf, "%s %s=%d affects .ac only: this td analysis linearises the transient model (quasi-static unless trnqsmod/nqsmod is on; card #559)", dn, kw, v.iValue);
                    *note = nbuf;
                }
            }
        }
    }
    return NULL;
}

/* Model-charge detection (one prepass period, cheap): a circuit whose only non-native devices carry no charge
 * (e.g. a diode with cj=0, tt=0 next to ideal switches) keeps the legacy native recursion -- bit-identical to
 * before and still able to use saltation. Each suspicious device TYPE is loaded alone (matrix cleared) at ag0
 * and ag0/2; any element that differs is a charge (C = dJ/dag0 != 0). */
static void
pac_detect_step(CKTcircuit *ckt)
{
    int neq = SMPmatSize(ckt->CKTmatrix), nstates = ckt->CKTnumStates, ti, i, r;
    double *save_rhs, *save_rhsold, *save_state0 = NULL;
    long save_mode = ckt->CKTmode;
    int save_noncon = ckt->CKTnoncon, save_bypass = ckt->CKTbypass;
    double a1 = ckt->CKTag[0];
    PACtrip A = { 0, 0, NULL, NULL, NULL }, B = A;

    if (ckt->CKTtime <= 0.0 || g_pp_detect_hit)
        return;
    save_rhs = TMALLOC(double, neq + 1);
    save_rhsold = TMALLOC(double, neq + 1);
    if (nstates > 0) save_state0 = TMALLOC(double, nstates);
    if (!save_rhs || !save_rhsold || (nstates > 0 && !save_state0)) {
        g_pp_detect_hit = 1;   /* cannot tell: take the generic path (safe) */
        goto out;
    }
    memcpy(save_rhs, ckt->CKTrhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(save_rhsold, ckt->CKTrhsOld, (size_t) (neq + 1) * sizeof(double));
    if (nstates > 0) memcpy(save_state0, ckt->CKTstate0, (size_t) nstates * sizeof(double));
    memcpy(ckt->CKTrhsOld, ckt->CKTrhs, (size_t) (neq + 1) * sizeof(double));
    ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODETRAN | MODEINITFLOAT;
    ckt->CKTbypass = 0;
    for (ti = 0; ti < g_pp_nsusp && !g_pp_detect_hit; ti++) {
        int dev = g_pp_susp[ti], k, err = 0;
        PACtrip *tt[2]; double ag[2];
        tt[0] = &A; tt[1] = &B; ag[0] = a1; ag[1] = 0.5 * a1;
        for (k = 0; k < 2 && !err; k++) {
            int nnz;
            SMPclear(ckt->CKTmatrix);
            for (r = 0; r <= neq; r++) ckt->CKTrhs[r] = 0.0;
            ckt->CKTag[0] = ag[k];
            err = DEVices[dev]->DEVload(ckt->CKThead[dev], ckt);
            ckt->CKTag[0] = a1;
            if (err) break;
            nnz = SMPnumElements(ckt->CKTmatrix);
            if (pac_trip_fit(tt[k], nnz)) { err = 1; break; }
            tt[k]->n = SMPcopyElements(ckt->CKTmatrix, tt[k]->ri, tt[k]->ci, tt[k]->v);
        }
        if (err || A.n != B.n) { g_pp_detect_hit = 1; break; }
        for (i = 0; i < A.n; i++)
            if (fabs(A.v[i] - B.v[i]) > 1e-14 * (fabs(A.v[i]) + fabs(B.v[i]))) { g_pp_detect_hit = 1; break; }
    }
    ckt->CKTmode = save_mode;
    ckt->CKTbypass = save_bypass;
    ckt->CKTnoncon = save_noncon;
    ckt->CKTag[0] = a1;
    memcpy(ckt->CKTrhs, save_rhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
    if (nstates > 0) memcpy(ckt->CKTstate0, save_state0, (size_t) nstates * sizeof(double));
out:
    if (save_rhs) tfree(save_rhs);
    if (save_rhsold) tfree(save_rhsold);
    if (save_state0) tfree(save_state0);
    pac_trip_free(&A); pac_trip_free(&B);
}

/* card #542: periodic-continuation start of the recording runs (see g_pac_xinit). */
static void
pac_init_free(void)
{
    if (g_pac_xinit) tfree(g_pac_xinit);
    g_pac_xinit = NULL;
    g_pac_xinit_n = 0;
}

static void
pac_init_from_prepass(CKTcircuit *ckt)
{
    int neq = SMPmatSize(ckt->CKTmatrix);
    pac_init_free();
    if (getenv("PAC_NO_IC_INJECT") || !g_pp_xend)
        return;
    g_pac_xinit = TMALLOC(double, neq + 1);
    if (g_pac_xinit) {
        memcpy(g_pac_xinit, g_pp_xend, (size_t) (neq + 1) * sizeof(double));
        g_pac_xinit_n = neq + 1;
    }
}

/* state snapshot around the detection period */
typedef struct {
    double *rhs, *rhsold, *states[8];
    int nstates, nrhs, nst;
    double ag[7], deltaOld[7], delta, time;
    long mode;
    int order;
    int ok;
} PACsnap;

static void
pac_snap_take(CKTcircuit *ckt, PACsnap *sn)
{
    int i, neq = SMPmatSize(ckt->CKTmatrix);
    memset(sn, 0, sizeof *sn);
    sn->nrhs = neq + 1;
    sn->nstates = ckt->CKTnumStates;
    sn->nst = ckt->CKTmaxOrder + 2;
    if (sn->nst > 8) sn->nst = 8;
    sn->rhs = TMALLOC(double, sn->nrhs);
    sn->rhsold = TMALLOC(double, sn->nrhs);
    if (!sn->rhs || !sn->rhsold) return;
    memcpy(sn->rhs, ckt->CKTrhs, (size_t) sn->nrhs * sizeof(double));
    memcpy(sn->rhsold, ckt->CKTrhsOld, (size_t) sn->nrhs * sizeof(double));
    for (i = 0; i < sn->nst && sn->nstates > 0; i++) {
        if (!ckt->CKTstates[i]) continue;
        sn->states[i] = TMALLOC(double, sn->nstates);
        if (!sn->states[i]) return;
        memcpy(sn->states[i], ckt->CKTstates[i], (size_t) sn->nstates * sizeof(double));
    }
    for (i = 0; i < 7; i++) { sn->ag[i] = ckt->CKTag[i]; sn->deltaOld[i] = ckt->CKTdeltaOld[i]; }
    sn->delta = ckt->CKTdelta; sn->time = ckt->CKTtime; sn->mode = ckt->CKTmode; sn->order = ckt->CKTorder;
    sn->ok = 1;
}

static void
pac_snap_restore(CKTcircuit *ckt, PACsnap *sn, int restore)
{
    int i;
    if (restore && sn->ok) {
        memcpy(ckt->CKTrhs, sn->rhs, (size_t) sn->nrhs * sizeof(double));
        memcpy(ckt->CKTrhsOld, sn->rhsold, (size_t) sn->nrhs * sizeof(double));
        for (i = 0; i < sn->nst; i++)
            if (sn->states[i] && ckt->CKTstates[i])
                memcpy(ckt->CKTstates[i], sn->states[i], (size_t) sn->nstates * sizeof(double));
        for (i = 0; i < 7; i++) { ckt->CKTag[i] = sn->ag[i]; ckt->CKTdeltaOld[i] = sn->deltaOld[i]; }
        ckt->CKTdelta = sn->delta; ckt->CKTtime = sn->time; ckt->CKTmode = sn->mode; ckt->CKTorder = sn->order;
    }
    if (sn->rhs) tfree(sn->rhs);
    if (sn->rhsold) tfree(sn->rhsold);
    for (i = 0; i < 8; i++) if (sn->states[i]) tfree(sn->states[i]);
}

static int
pac_generic_prepare(CKTcircuit *ckt, PSSan *job, PACstateMap *map, const double *x0, double T,
                    int trig_node_num, const char *who)
{
    const char *delay = NULL;
    int has, err, neq, i, ns, npass, force_gen;
    static int warned_native = 0;
    double tprep0 = SPfrontEnd->IFseconds();

    /* card #551/#557: stateless analog XSPICE models have a td small-signal model (J exact); s_xfer/int/d_dt/slew/hyst through the hidden-state
     * recursion (.pac/.pstb/.pxf), all other stateful or digital models are refused; .pnoise (kind 3) refuses the stateful ones too. */
    if (PSSxspicePolicy(ckt, !strcmp(who, "PNOISE-TD") ? 3 : 2, who))
        return E_UNSUPP;
#ifdef XSPICE
    pac_hid_collect(ckt);
    if ((g_hid_nc > 0 || g_hid_nhy > 0) && ckt->CKTintegrateMethod == GEAR) {
        fprintf(stderr, "%s: REFUSED -- a stateful XSPICE code model (s_xfer int d_dt slew hyst) with method=gear: the hidden-state recursion "
                "is derived for the trapezoidal / backward-Euler integral of cm_analog_integrate only (card #557); use method=trap.\n", who);
        pac_hid_free();
        return E_UNSUPP;
    }
#endif

    if (getenv("PAC_NATIVE_HISTORY")) {
        if (!warned_native) {
            fprintf(stderr, "%s: PAC_NATIVE_HISTORY set -- using the LEGACY native C/L-only history "
                    "(device-model charges ignored: wrong on MOS/BJT/diode circuits, card #539 red proof)\n", who);
            warned_native = 1;
        }
#ifdef XSPICE
        pac_hid_free();
#endif
        return 0;
    }
    has = pac_circuit_has_model_charges(ckt, &delay);
    if (delay) {
        fprintf(stderr, "%s: REFUSED -- device type '%s' carries delay history that is not a charge state; "
                "td small-signal analysis of it is not implemented (card #542).\n", who, delay);
        return E_UNSUPP;
    }
    {
        const char *note = NULL;
        const char *nqs = pac_find_nqs_device(ckt, &note, !strcmp(who, "PNOISE-TD"));
        if (nqs && !getenv("PAC_ALLOW_NQS")) {
            fprintf(stderr, "%s: REFUSED -- %s: the non-quasi-static model carries an internal charge-relaxation state "
                    "that is not validated here (card #542; BSIM4 trnqsmod and BSIM3 nqsmod are accepted in .pac/.pstb/.pxf since "
                    "card #559, not in .pnoise). Use the quasi-static model (trnqsmod=0 / nqsmod=0). PAC_ALLOW_NQS=1 overrides "
                    "(result then unvalidated).\n", who, nqs);
            return E_UNSUPP;
        }
        if (note)
            fprintf(stderr, "%s: non-quasi-static MOS model: %s\n", who, note);
    }
    /* card #549: method=gear is carried by the generic (charge-history) recursion only, so it is taken for native decks too */
    force_gen = getenv("PAC_FORCE_GENERIC") != NULL || ckt->CKTintegrateMethod == GEAR;
    if (!force_gen && !getenv("PAC_NO582_GENERIC") && pac_switch_node_without_cap(ckt)) {
        /* card #582: an ideal switch at a node with no capacitor (and no voltage source): the legacy per-C/L recursion was 12 dB / 94 deg
         * off on a half bridge without dead time, node capacitor or diode (same events, same orbit; cause not found), the generic
         * recursion is right there and equal or better on every chargeless matrix deck: use it for such circuits */
        force_gen = 1;
    }
#ifdef XSPICE
    if (g_hid_nc > 0 || g_hid_nhy > 0) force_gen = 1;   /* card #557 */
#endif
    if (!has && !force_gen)
        return 0;
    neq = SMPmatSize(ckt->CKTmatrix);
    if (!force_gen && !g_pp_has_mutual) {
        /* the only non-native devices may carry no charge at all: look (one cheap period) before leaving the
         * legacy native path, which is bit-identical to before and supports saltation / free-running */
        /* The detection period must leave no trace on the legacy flow that follows (its first recorded period starts
         * from whatever vectors and device states are left in the circuit): snapshot and restore them. */
        PACsnap snap;
        pac_snap_take(ckt, &snap);
        pac_state_set_ic(map, x0);
        g_pp_detect = 1; g_pp_detect_hit = 0; g_pp_active = 1; PSSpacActive = 1;
        err = pac_run_period(ckt, job, T, "TD-shoot model-charge detection (card #542)");
        PSSpacActive = 0; g_pp_active = 0; g_pp_detect = 0;
        pac_snap_restore(ckt, &snap, g_pp_detect_hit ? 0 : 1);
        if (err) {
            fprintf(stderr, "%s: model-charge detection run failed (err=%d)\n", who, err);
            return E_INTERN;
        }
        if (!g_pp_detect_hit)
            return 0;   /* no device-model charge anywhere on the orbit: native C/L history is exact */
    }
    /* card #544: free-running (trigger=) orbits with device-model charges are accepted. Their .pac/.pstb is the fixed-period
     * engine on the polished orbit (the period-variation projection is opt-in, PAC_FREERUN_PROJ, and is not formed for a
     * generic state set: f(T) over the charge rows is not available). The trigger node must be a state row (below). */
    g_pp_col = TMALLOC(unsigned char, neq + 1);
    if (!g_pp_col) return E_NOMEM;
    if (g_pp_xend) tfree(g_pp_xend);
    g_pp_xend = TMALLOC(double, neq + 1);
    if (!g_pp_xend) { tfree(g_pp_col); g_pp_col = NULL; return E_NOMEM; }
    g_pp_neq = neq;
    /* card #559: switch power nodes with device-model charge (recorded by the prepass steps) */
    PSSsaltModelCapHook = NULL;
    if (g_salq) { tfree(g_salq); g_salq = NULL; }
    if (g_salqC) { tfree(g_salqC); g_salqC = NULL; }
    if (g_salqCmin) { tfree(g_salqCmin); g_salqCmin = NULL; }
    PSSsaltModelCapMinHook = NULL;
    g_salq_n = 0; g_salq_any = 0;
    if (CKTtypelook("Switch") >= 0 && ckt->CKThead[CKTtypelook("Switch")]) {
        g_salq = TMALLOC(unsigned char, neq + 1);
        g_salqC = TMALLOC(double, neq + 1);
        g_salqCmin = TMALLOC(double, neq + 1);
        if (g_salqCmin) { int i_; for (i_ = 0; i_ <= neq; i_++) g_salqCmin[i_] = -1.0; }
        if (g_salq && g_salqC) g_salq_n = neq;
        else { if (g_salq) tfree(g_salq); if (g_salqC) tfree(g_salqC); g_salq = NULL; g_salqC = NULL; }
    }
    /* two prepass periods: the first from whatever state the PSS left, the second from the first one's
     * end state (periodic continuation, see g_pac_xinit); the second one is what the recording passes
     * reproduce, so its steps/columns define the state set. */
    for (npass = 0; npass < 2; npass++) {
        memset(g_pp_col, 0, (size_t) neq + 1);
        if (g_salq) { memset(g_salq, 0, (size_t) neq + 1); memset(g_salqC, 0, (size_t) (neq + 1) * sizeof(double)); g_salq_any = 0; }
        memset(g_pp_xend, 0, (size_t) (neq + 1) * sizeof(double));
        g_pp_steps = 0; g_pp_cmax = 0.0;
        g_pac_gen_bad = 0; g_pac_gen_badname[0] = '\0';
#ifdef XSPICE
        g_hid_flips = 0;
        if (g_hid_nc > 0) {
            memset(g_hid_patP, 0, (size_t) g_hid_nc);
            memset(g_hid_patB, 0, (size_t) g_hid_nc);
            memset(g_hid_patR, 0, (size_t) g_hid_nc * g_hid_nc);
        }
#endif
        pac_state_set_ic(map, x0);
        g_pp_active = 1;
        PSSpacActive = 1;
        err = pac_run_period(ckt, job, T, "TD-shoot charge-state prepass (card #542)");
        PSSpacActive = 0;
        g_pp_active = 0;
        if (g_pac_gen_bad > 0 || err || g_pac_gen_bad < 0 || g_pp_steps < 1)
            break;
        pac_init_from_prepass(ckt);
#ifdef XSPICE
        g_hid_stride = (g_pp_steps > 20000) ? (g_pp_steps + 19999) / 20000 : 1;
#endif
    }
#ifdef XSPICE
    g_hid_stride = 1;
#endif
    if (g_pac_gen_bad > 0) {
        pac_gen_refuse_message(who);
        tfree(g_pp_col); g_pp_col = NULL;
        tfree(g_pp_xend); g_pp_xend = NULL;
        pac_init_free();
        return E_UNSUPP;
    }
    if (err || g_pac_gen_bad < 0 || g_pp_steps < 1) {
        fprintf(stderr, "%s: charge-state prepass failed (err=%d steps=%d)\n", who, err, g_pp_steps);
        tfree(g_pp_col); g_pp_col = NULL;
        tfree(g_pp_xend); g_pp_xend = NULL;
        pac_init_free();
        return E_INTERN;
    }
    ns = 0;
    for (i = 1; i <= neq; i++) if (g_pp_col[i]) ns++;
#ifdef XSPICE
    if (g_hid_nc > 0) pac_hid_finalize();
    if (g_hid_slewact && !getenv("PAC_HID_ALLOW_SLEW")) {   /* PAC_HID_ALLOW_SLEW=1: experiment (the error of the missing saltation) */
        fprintf(stderr, "%s: REFUSED -- slew-rate limiting of a 'slew' code model is active on the orbit: the change between the slewing and the "
                "tracking regime is an event whose small-signal map (saltation) is not derived (card #557); a slew model that never slews is "
                "accepted (its small-signal gain is the tracking one).\n", who);
        tfree(g_pp_col); g_pp_col = NULL; tfree(g_pp_xend); g_pp_xend = NULL;
        pac_init_free(); pac_hid_free();
        return E_UNSUPP;
    }
    if (g_hid_flips > 0 && !getenv("PAC_HID_ALLOW_FLIP")) {   /* PAC_HID_ALLOW_FLIP=1: experiment */
        fprintf(stderr, "%s: REFUSED -- the Boolean state of a 'hyst' code model changes on the orbit (%d steps with a state flip): the flip is an "
                "event whose small-signal map (saltation) is not derived (card #557); a hyst model that stays on one branch is accepted.\n",
                who, g_hid_flips);
        tfree(g_pp_col); g_pp_col = NULL; tfree(g_pp_xend); g_pp_xend = NULL;
        pac_init_free(); pac_hid_free();
        return E_UNSUPP;
    }
#endif
    map->srows = TMALLOC(int, ns + g_hid_nh + 1);
    if (!map->srows) { tfree(g_pp_col); g_pp_col = NULL; return E_NOMEM; }
    ns = 0;
    for (i = 1; i <= neq; i++) if (g_pp_col[i]) map->srows[ns++] = i;
#ifdef XSPICE
    if (g_hid_nh > 0) {
        int a;
        g_hid_Bh = TMALLOC(double, (size_t) (neq + 1) * g_hid_nh);
        g_hid_P = TMALLOC(double, (size_t) g_hid_nh * g_hid_nin);
        g_hid_R = TMALLOC(double, (size_t) g_hid_nh * g_hid_nh);
        g_hid_Jc = TMALLOC(double, (size_t) (neq + 1) * g_hid_nin);
        if (!g_hid_Bh || !g_hid_P || !g_hid_R || !g_hid_Jc) { tfree(g_pp_col); g_pp_col = NULL; return E_NOMEM; }
        for (a = 0; a < g_hid_nh; a++) map->srows[ns++] = neq + 1 + a;   /* virtual rows of the hidden states */
    }
#endif
    map->ns = ns;
    map->generic = 1;
    if (g_salq && g_salq_any) {
        int sn, nq = 0;
        for (sn = 1; sn <= g_salq_n; sn++) if (g_salq[sn]) nq++;
        PSSsaltModelCapHook = pac_salq_hook;
        PSSsaltModelCapMinHook = pac_salq_min_hook;
        fprintf(stderr, "%s: ideal-switch power nodes with device-model charge (card #559): %d node(s); their saltation jump is formed on the "
                "charge vector%s\n", who, nq, pac_salq_off() ? " [DISABLED by PAC_SALT_NOQ / PAC_SALT_NATIVEJUMP: native-only jump]" : "");
    }
    tfree(g_pp_col); g_pp_col = NULL;
    tfree(g_pp_xend); g_pp_xend = NULL;
    fprintf(stderr, "%s: generic charge-based history (card #542): %d native C/L + model-charge states -> "
            "%d state rows (columns of dq/dx on the orbit, %d steps, affine-in-ag0 check passed, max|C|=%.3g), "
            "prepass wall_time=%.3fs\n",
            who, map->n, ns - g_hid_nh, g_pp_steps, g_pp_cmax, SPfrontEnd->IFseconds() - tprep0);
    if (g_hid_nh > 0)
        fprintf(stderr, "%s: stateful XSPICE code models (card #557): %d hidden states (of %d state slots), %d input rows; "
                "finite-difference linearisation of the models, state carried between recording runs\n", who, g_hid_nh, g_hid_nc, g_hid_nin);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Core engine: given ONE forcing pattern (br/bi, shared across the whole
 * sweep) and NFREQ frequencies, returns Phi (n*n, real part only -- see
 * this file's header, the homogeneous columns never leave the real
 * subspace) and, for every requested output ROW (1-based, e.g. a node or
 * a branch equation) and frequency, the demodulated sideband-0 response
 * Y0. Runs exactly 2 period re-runs total (batched over ALL frequencies),
 * per this file's header design. Returns 0 on success. */
static int
pac_engine_run(CKTcircuit *ckt, PSSan *job, PACstateMap *map, const double *x0, double T,
                const double *br, const double *bi,
                const double *freqs, int nfreq,
                const int *outrows, int nout,
                PACc *Phi_out /* n*n, may be NULL */,
                const PACc *Phi_in /* card #555: n*n monodromy of an earlier run on the SAME orbit, else NULL: pass 1 then
                                      propagates only the forced columns */,
                PACc *Y0_out  /* nout*nfreq */,
                const double *fT /* card #489 second pass: n, NULL = ordinary
                                     driven (fixed-period) boundary condition */,
                int trig_idx /* index of the trigger state in fT/map's own
                                order, -1 = no correction (driven case, or
                                fT==NULL) */)
{
    int n = PAC_NST(map), neq = SMPmatSize(ckt->CKTmatrix);
    int Mst = Phi_in ? 0 : n;   /* homogeneous (state) columns of pass 1 */
    int M1 = Mst + nfreq, f, d, j, err;
    double *omega1, *ic_zero_placeholder = NULL;
    int *useforce1;
    PACc *ic0_1, *Phi, *z;
    PACc *xi;

    if (n <= 0 || nfreq <= 0)
        return 1;

    omega1 = TMALLOC(double, M1);
    useforce1 = TMALLOC(int, M1);
    ic0_1 = TMALLOC(PACc, (size_t) n * M1);
    Phi = TMALLOC(PACc, (size_t) n * n);
    z = TMALLOC(PACc, (size_t) n * nfreq);
    xi = TMALLOC(PACc, (size_t) n * nfreq);
    NG_IGNORE(ic_zero_placeholder);

    if (!omega1 || !useforce1 || !ic0_1 || !Phi || !z || !xi) {
        if (omega1) tfree(omega1); if (useforce1) tfree(useforce1);
        if (ic0_1) tfree(ic0_1); if (Phi) tfree(Phi); if (z) tfree(z); if (xi) tfree(xi);
        return 1;
    }

    memset(ic0_1, 0, (size_t) n * M1 * sizeof(PACc));
    for (d = 0; d < Mst; d++) {
        omega1[d] = 0.0;
        useforce1[d] = 0;
        for (j = 0; j < n; j++) {
            PACc v = { (d == j) ? 1.0 : 0.0, 0.0 };
            ic0_1[(size_t) j * M1 + d] = v;   /* column d: unit at reduced-state row d */
        }
    }
    for (f = 0; f < nfreq; f++) {
        omega1[Mst + f] = 2.0 * M_PI * freqs[f];
        useforce1[Mst + f] = 1;
    }

    pac_state_set_ic(map, x0);
    if (pac_reset(ckt, map, M1, br, bi, omega1, useforce1, ic0_1, 0, T)) {
        err = 1; goto cleanup;
    }
    PSSpacActive = 1;
    err = pac_run_period(ckt, job, T, "TD-shoot recording pass 1 (monodromy + forced)");
    PSSpacActive = 0;
    if (err || g_pac_failed) { err = 1; goto cleanup; }

    for (d = 0; d < n; d++) {
        int row_d = pac_map_row(map, d);
        if (!Phi_in)
            for (j = 0; j < n; j++)
                Phi[(size_t) d * n + j] = pac_get(g_pac_Cprev, row_d, M1, j);
        for (f = 0; f < nfreq; f++)
            z[(size_t) d * nfreq + f] = pac_get(g_pac_Cprev, row_d, M1, Mst + f);
    }
    if (Phi_in) memcpy(Phi, Phi_in, (size_t) n * n * sizeof(PACc));
    /* card #465 FOURTH pass, task 4 (red-proof, TEST-ONLY, default off):
     * PAC_REDPROOF_FREEZE_PHI substitutes the identity matrix for the real,
     * per-step-propagated LPTV monodromy Phi computed above -- a stand-in
     * for the WRONG "averaged/frozen" linearisation the card's plan calls
     * out (throwing away the periodic Jacobian's own time-modulation,
     * exactly what an averaged/frozen-operating-point small-signal model
     * does). This must make the boundary-value solve below wrong (it
     * changes the (e^{jwT}I-Phi) xi = z(w) equation actual LTI systems
     * satisfy only when Phi=I, i.e. w=0), and is asserted BY THE TEST
     * SUITE (tests/test_pac.py's check_redproof_frozen_phi()) to fail
     * pac_sb0(b) (R2 exact) and pstb_loop (c) (R3) -- proving those checks
     * would actually catch this class of regression, not just that they
     * currently pass. Never set in production; not documented in any
     * user-facing help text. */
    if (getenv("PAC_REDPROOF_FREEZE_PHI")) {
        for (d = 0; d < n; d++)
            for (j = 0; j < n; j++) {
                PACc v = { (d == j) ? 1.0 : 0.0, 0.0 };
                Phi[(size_t) d * n + j] = v;
            }
    }
    if (Phi_out) memcpy(Phi_out, Phi, (size_t) n * n * sizeof(PACc));
    if (getenv("PXF_ADJ_DEBUG")) {
        fprintf(stderr, "PAC-V1-DEBUG: n=%d exact_jac=%d\n", n, g_pac_exact_jac);
        for (d = 0; d < n; d++) {
            for (j = 0; j < n; j++)
                fprintf(stderr, "  Phi[%d][%d]=%.15g%+.15gi", d, j,
                        Phi[(size_t) d * n + j].re, Phi[(size_t) d * n + j].im);
            fprintf(stderr, "\n");
        }
    }
    pac_free();

    /* card #489 SECOND pass (plan T6.2, task 5): free-running period-
     * variation correction. Applied to LOCAL copies (PhiB/zB below) used
     * ONLY for the boundary-value solve immediately below -- the stored
     * Phi_out (if requested) stays the RAW monodromy, and pass 2's
     * re-propagation below is unaffected (see docs/pss-pac.md's
     * free-running section: the re-propagation/demodulation correction is
     * second order in the response and does not need this term).
     * Proj = I - f(T) e_trig^T / (e_trig . f(T)) is never built as a
     * dense matrix: (Proj . A)[d][*] = A[d][*] - (fT[d]/fT[trig_idx]) *
     * A[trig_idx][*] is a rank-1 update using ROW trig_idx of A, applied
     * to every row d (including d==trig_idx itself, which becomes
     * identically zero -- the expected rank-(n-1) signature, same as
     * pssshoot.c's own FD-projected monodromy M). */
    {
        PACc *PhiB = Phi, *zB = z;
        PACc *PhiCorr = NULL, *zCorr = NULL;
        /* card #489 SECOND pass, task 5 validation (TEST-ONLY, default off,
         * same convention as PAC_REDPROOF_FREEZE_PHI/PAC_REDPROOF_SIGNFLIP
         * above): PAC_FREERUN_NO_PROJ skips the period-variation projection
         * below even when fT/trig_idx were supplied -- used to confirm the
         * term actually changes the answer (not a silent no-op) by A/B
         * comparison against the same independent reference this pass's
         * own report validates against. Never set in production. */
        if (map->generic && getenv("PAC_FREERUN_PROJ"))
            fprintf(stderr, "PAC-TD: PAC_FREERUN_PROJ ignored: the period-variation projection is not formed on the generic (model-charge) state set (card #544).\n");
        if (!map->generic && fT && trig_idx >= 0 && fT[trig_idx] != 0.0 && getenv("PAC_FREERUN_PROJ") && !getenv("PAC_FREERUN_NO_PROJ")) {
            double denom = fT[trig_idx];
            PhiCorr = TMALLOC(PACc, (size_t) n * n);
            zCorr = TMALLOC(PACc, (size_t) n * nfreq);
            if (!PhiCorr || !zCorr) {
                if (PhiCorr) tfree(PhiCorr); if (zCorr) tfree(zCorr);
                err = 1; goto cleanup;
            }
            for (d = 0; d < n; d++) {
                double coeff = fT[d] / denom;
                for (j = 0; j < n; j++)
                    PhiCorr[(size_t) d * n + j] =
                        pc_sub(Phi[(size_t) d * n + j],
                               pc_scale(Phi[(size_t) trig_idx * n + j], coeff));
                for (f = 0; f < nfreq; f++)
                    zCorr[(size_t) d * nfreq + f] =
                        pc_sub(z[(size_t) d * nfreq + f],
                               pc_scale(z[(size_t) trig_idx * nfreq + f], coeff));
            }
            PhiB = PhiCorr; zB = zCorr;
        }

    /* dense boundary solve per frequency: (e^{jwT} I - Phi) xi = z(w)
     * (PhiB/zB: the free-running-corrected boundary condition above when
     * fT/trig_idx were given, otherwise identically Phi/z -- the ordinary
     * driven-engine equation, byte-for-byte unchanged). */
    for (f = 0; f < nfreq; f++) {
        PACc *A = TMALLOC(PACc, (size_t) n * n);
        PACc *rhs = TMALLOC(PACc, n);
        PACc *sol = TMALLOC(PACc, n);
        double w = 2.0 * M_PI * freqs[f];
        double wT = w * T;
        PACc ephase = { cos(wT), sin(wT) };
        /* card #494 sixth pass (plan T6.1): VALIDATED-BAND GUARD. A
         * saltation-corrected result (ideal sw/csw, PSSsaltActive with
         * real events found this run -- PSSsaltEventsSeen(), an exact
         * no-op/0 for every existing fixture with no such event) is
         * cross-validated against tools/pac_perturbed.py, an independent
         * nonlinear perturbed-periodic-steady-state reference sharing no
         * code with this engine, to 0.5dB/3deg from 100Hz up to at least
         * 99.5% of the switching Nyquist frequency fs/2 (T here is the
         * deck's own switching period, the driven .pss's native clock --
         * see docs/pss-pac.md's "SIXTH pass" section for the full
         * cross-check table). The disagreement grows sharply in the last
         * few percent below fs/2 and becomes structural (a genuine
         * sideband-0 vs sideband-(-1) aliasing degeneracy for a REAL
         * system probed at exactly half the switching frequency, common
         * to every sideband-0-only method -- this engine's own general
         * "Known limits" already flag it for smooth-switch decks too, not
         * new to saltation) AT fs/2 itself (measured: 3.2dB/41deg off at
         * f=fs/2 exactly on the reference fixture). 0.90 (not 0.995) is
         * used as the trigger ratio for margin: only two points near the
         * edge were actually measured this pass (0.9877 clean, 0.9950
         * already at 0.29dB of the 0.5dB budget), and this correction has
         * only been validated on ONE fixture (pac_r6o_ctrl2out.cir) --
         * warn well before the measured budget is actually exhausted
         * rather than assume the margin generalises. Never silent: one
         * print per affected frequency point, matching this file's other
         * per-point WARNING convention (boundary-solve-singular, just
         * below). */
        if (PSSsaltActive && PSSsaltEventsSeen() > 0 && T > 0.0) {
            double fs_nyquist = 0.5 / T;
            double ratio = freqs[f] / fs_nyquist;
            if (ratio >= 0.90)
                fprintf(stderr, "PAC-TD: WARNING -- saltation-corrected result at f=%.6g Hz "
                        "is inside the UNVALIDATED margin near the switching Nyquist "
                        "frequency (fs/2=%.6g Hz, f/(fs/2)=%.4f >= 0.90); cross-validated "
                        "to 0.5dB/3deg only up to ~99.5%% of fs/2 on the reference fixture, "
                        "a structural sideband-0 aliasing breakdown is expected approaching "
                        "and at fs/2 itself -- see docs/pss-pac.md's saltation section, "
                        "SIXTH pass, and patches-wip/salt/README.md.\n",
                        freqs[f], fs_nyquist, ratio);
        }
        /* card #493: FREE-RUNNING VALIDATED-BAND GUARD. .pac/.pstb on a
         * free-running (trigger=) orbit is EXPERIMENTAL and refused by
         * default (PSSshootFreerun()'s own gate, pssshoot.c); reaching
         * this code at all means the caller already opted in with
         * PSS_FREERUN_PAC_EXPERIMENTAL=1. This card cross-validated BOTH
         * the period-variation-corrected result (fT/trig_idx given, no
         * PAC_FREERUN_NO_PROJ) AND the uncorrected one (PAC_FREERUN_NO_PROJ
         * set, or fT missing) against tools/pac_perturbed.py's --freerun
         * backend (an independent, nonlinear Hann-windowed perturbed-
         * transient DFT sharing no code with this engine) on R5
         * (refs/r5_hysteretic_buck.cir, fsw~=422.68kHz): agreement is
         * inside the project's usual 0.5dB/3deg target at 20kHz (4.7% of
         * fsw, +0.07dB/+0.08deg) and 60kHz (14.2%, +0.34dB/+0.61deg,
         * corrected) but breaks down sharply by 100kHz (23.7%,
         * +1.1dB/+3.0deg) and is wildly wrong by 150-200kHz (35%/47%,
         * several dB and tens of degrees) -- BOTH corrected and
         * uncorrected, growing together, not diverging from each other.
         * Root cause NOT in the period-variation projection term itself
         * (ruled out this card: a from-scratch, controlled fine forward-
         * difference f(T) estimate, PSS_FREERUN_FT_BACKWARD unset -- see
         * pssshoot.c's PSSshootFreerun() -- barely changes the result;
         * PAC_EXACT_JACOBIAN=1, this file's own unrelated documented
         * accuracy-floor toggle, changes the 150kHz result by <1e-4
         * relative; pacpoints=4000 vs the fixture's default 200 changes
         * the wrong answer to a DIFFERENT wrong answer, not a smaller
         * one) -- consistent with a general per-step-linearisation
         * accuracy floor on R5's own extremely stiff switch node (Csw,
         * ~1ps local RC, far faster than any R3/R4 state) that is
         * independent of trigger=/free-running specifically, just never
         * stress-tested against dynamics this fast before this project's
         * R5 fixtures existed -- see docs/pss-pac.md's card #493 section
         * for the full numbers table and what was and was not tried.
         * 0.10 (10% of the measured free-running switching frequency, NOT
         * fs/2 -- there is no fixed fs here) is used as the trigger ratio:
         * conservative relative to the ~14% point actually measured
         * clean, since only ONE fixture/output was checked. */
        if (trig_idx >= 0 && T > 0.0) {
            double fsw = 1.0 / T;
            double ratio = freqs[f] / fsw;
            if (ratio >= 0.10)
                fprintf(stderr, "PAC-TD: WARNING -- free-running (trigger=) "
                        ".pac/.pstb result at f=%.6g Hz is inside the "
                        "UNVALIDATED margin relative to the measured "
                        "free-running switching frequency (fsw=%.6g Hz, "
                        "f/fsw=%.4f >= 0.10); cross-validated to 0.5dB/3deg "
                        "only up to ~14%% of fsw on R5, breaking down "
                        "sharply by ~24%% and badly wrong by 35%% -- both "
                        "the period-variation-corrected and the uncorrected "
                        "result are affected, root cause not isolated to "
                        "the correction term -- see docs/pss-pac.md's card "
                        "#493 section.\n", freqs[f], fsw, ratio);
        }
        if (!A || !rhs || !sol) {
            if (A) tfree(A); if (rhs) tfree(rhs); if (sol) tfree(sol);
            if (PhiCorr) tfree(PhiCorr); if (zCorr) tfree(zCorr);
            err = 1; goto cleanup;
        }
        for (d = 0; d < n; d++)
            for (j = 0; j < n; j++) {
                PACc v = pc_scale(PhiB[(size_t) d * n + j], -1.0);
                if (d == j) v = pc_add(v, ephase);
                A[(size_t) d * n + j] = v;
            }
        for (d = 0; d < n; d++) rhs[d] = zB[(size_t) d * nfreq + f];
        if (pac_complex_solve(A, rhs, n, sol)) {
            fprintf(stderr, "PAC-TD: WARNING -- boundary solve singular at f=%.6g Hz "
                    "(e^{jwT}=1 is a Floquet multiplier there: an undamped oscillation "
                    "or an unstable/marginal orbit); result at this frequency skipped.\n",
                    freqs[f]);
            for (d = 0; d < n; d++) { PACc nan = { strtod("NaN", NULL), strtod("NaN", NULL) }; sol[d] = nan; }
        }
        for (d = 0; d < n; d++)
            xi[(size_t) d * nfreq + f] = sol[d];
        tfree(A); tfree(rhs); tfree(sol);
    }
    if (PhiCorr) tfree(PhiCorr);
    if (zCorr) tfree(zCorr);
    }

    /* pass 2: re-propagate from y0=xi per frequency, demodulate at the
     * requested output rows. */
    {
        int M2 = nfreq;
        double *omega2 = TMALLOC(double, M2);
        int *useforce2 = TMALLOC(int, M2);
        PACc *ic0_2 = TMALLOC(PACc, (size_t) n * M2);
        if (!omega2 || !useforce2 || !ic0_2) {
            if (omega2) tfree(omega2); if (useforce2) tfree(useforce2); if (ic0_2) tfree(ic0_2);
            err = 1; goto cleanup;
        }
        for (f = 0; f < nfreq; f++) {
            omega2[f] = 2.0 * M_PI * freqs[f];
            useforce2[f] = 1;
            for (d = 0; d < n; d++)
                ic0_2[(size_t) d * M2 + f] = xi[(size_t) d * nfreq + f];
        }
        pac_state_set_ic(map, x0);
        if (pac_reset(ckt, map, M2, br, bi, omega2, useforce2, ic0_2, 1, T)) {
            tfree(omega2); tfree(useforce2); tfree(ic0_2);
            err = 1; goto cleanup;
        }
        PSSpacActive = 1;
        err = pac_run_period(ckt, job, T, "TD-shoot recording pass 2 (re-propagate)");
        PSSpacActive = 0;
        if (!err && !g_pac_failed) {
            /* card #465 FOURTH pass, task 1: supply the [0,t_1] demod
             * segment that PSSpacRecordStep() deferred -- see
             * pac_finish_demod_wraparound()'s header for the exact
             * (periodicity-based) closure. Must run before pac_free()
             * below (needs g_pac_Cprev, still holding y(T)). */
            pac_finish_demod_wraparound();
            for (f = 0; f < nfreq; f++) {
                double invT = 1.0 / T;
                int o;
                for (o = 0; o < nout; o++) {
                    PACc v = g_pac_demod[(size_t) (outrows[o] - 1) * M2 + f];
                    Y0_out[(size_t) o * nfreq + f] = pc_scale(v, invT);
                }
            }
            err = 0;
        } else {
            err = 1;
        }
        pac_free();
        tfree(omega2); tfree(useforce2); tfree(ic0_2);
    }

 cleanup:
    tfree(omega1); tfree(useforce1); tfree(ic0_1); tfree(Phi); tfree(z); tfree(xi);
    return err;
}

/* ------------------------------------------------------------------ */
/* card #504: TRUE ADJOINT `.pxf ... adjoint=1` engine.
 *
 * DERIVATION (full working kept in the card #504 report; summary here).
 * pac_engine_run()'s TWO forward period passes, for one source's forcing
 * pattern b (br/bi), compute:
 *   pass 1: Phi (n x n, homogeneous columns, b=0) and z(w) = E_read(u_N(b))
 *           (the b-forced, ZERO-IC trajectory's final state at the n
 *           reduced-state rows, E_read = pac_state_row()'s OWN
 *           posNode/branch-only convention -- NOT the Delta=p-n difference
 *           the INTERNAL per-step recursion below uses; the two coincide
 *           for every fixture here since every capacitor's negNode is
 *           ground, see the design note by pac_state_row() itself.)
 *   boundary solve: xi(w) = (e^{jwT}I - Phi)^{-1} z(w)
 *   pass 2: re-propagate from y0=xi WITH the SAME b forcing, demodulate at
 *           the output row -> Y0(w) = xi^T G(w) + h(w,b), where G(w)[d] is
 *           "the Y0 a homogeneous-IC=e_d, b=0 pass 2 would give" and h(w,b)
 *           is pass 2's OWN b-forced (zero-IC) demodulated contribution.
 *
 * Since pass 1's z(w) IS "the b-forced zero-IC trajectory's E_read(u_N)"
 * and pass 2's h(w,b) IS "the SAME b-forced zero-IC trajectory's
 * demodulated output" -- both read off the literal SAME trajectory u(b)
 * (b appears identically, at every step, in both -- this file's own
 * header notes re-running from the converged x0 reproduces the transient
 * bit for bit) -- superposition gives, for the SAME output row and SAME w:
 *
 *   Y0(w,b) = [E_read(u_N(b))]^T M(w)  +  (1/T) sum_k alpha_k[out] u_k(b)[out]
 *
 * where M(w) := ((e^{jwT}I - Phi)^T)^{-1} G(w) (the TASK'S "transpose of
 * the periodic boundary condition") and alpha_k[out] is the RAW (pre-1/T)
 * trapezoidal-plus-wraparound demodulation weight PSSpacRecordStep()/
 * pac_finish_demod_wraparound() already compute (see their own headers;
 * re-derived in closed periodic form below, coef_k = ph2_k * 0.5*(h_k +
 * h_{k+1 mod N}), h_{N+1}:=h_1 -- the expected periodic-trapezoidal
 * quadrature weight).
 *
 * u_k(b) is itself generated by ONE linear per-step recursion (b entering
 * identically, phase-modulated, at every step) -- Y0(w,.) as a function of
 * b is therefore a single COMPLEX-LINEAR functional, representable as
 * Y0(w,b) = b^T Lambda(w) for one Lambda(w) in C^neq that does NOT depend
 * on b at all. Lambda(w) is computed by a BACKWARD (reverse-time) sweep of
 * the per-step recursion's adjoint, folding in M(w) at the LAST step only
 * (k=N) and the alpha_k[out] demod weight at EVERY step -- the "transposed
 * per-step solve" (SMPtSolve, new this card) and "transpose of the
 * per-step history operator" (the trapezoidal order-2 two-term T_k
 * recursion, transposed below into beta/gamma "gather" formulas that are
 * the EXACT mirror of PSSpacRecordStep()'s own dceq/dveq/T_k "scatter"
 * formulas -- re-derived directly from that code, not from physics
 * intuition, to avoid a docs/sens.md-style transpose-or-sign bug) the
 * card's brief calls for. Once Lambda(w) is known, EVERY requested
 * source's PXF is a single dot product b_source^T Lambda(w) -- the "K
 * sources for the cost of ~1 sweep" payoff.
 *
 * J_k (needed at arbitrary, REVERSE order): RECOMPUTED on demand from a
 * recorded trajectory (g_pac_traj[], card #504 above -- x_k plus
 * t_k/h_k/ag0_k/ag1_k/order_k, captured during this function's OWN
 * "Phi + trajectory" forward pass) via a fresh CKTload()+SMPluFac() per
 * backward step, mirroring pac_refactor_at_accepted()'s existing,
 * validated "put x_k where devices read node voltages from, reload,
 * refactor" recipe -- NOT deep-cloned/stored factorizations (memory: see
 * the card #504 report's measured numbers; O(N*(neq+2)) doubles for the
 * trajectory alone was judged the simpler, safer choice for one session,
 * not proven the ONLY viable one -- storing factorizations remains valid
 * future work if per-step refactor cost ever dominates for a much larger
 * circuit than this project's own fixtures). Because this recompute-based
 * J_k is INTERNALLY CONSISTENT only with itself (x_k, not the live
 * transient's possibly-stale x_{k-1} the DEFAULT, non-PAC_EXACT_JACOBIAN
 * v1 path uses -- see pac_refactor_at_accepted()'s own header), this
 * engine's OWN forward Phi/G passes ALSO force the exact-Jacobian
 * convention (g_pac_force_exact_jac, independent of the PAC_EXACT_JACOBIAN
 * env var) so its own forward and backward halves use the IDENTICAL J_k
 * sequence -- a genuine, well-defined adjoint of ONE map, not an
 * accidental mix of two. Consequence: comparing this engine against v1's
 * OWN DEFAULT (stale-Jacobian, O(1/Points)-accurate) output is only
 * accurate to O(1/Points), NOT floating-point noise; tests/test_pxf.py's
 * adjoint checks therefore run v1 ALSO under PAC_EXACT_JACOBIAN=1 for the
 * comparison (both sides then mathematically exact, agreement limited
 * only by roundoff) -- see that script and the card #504 report's
 * "validation" section for the numbers either way.
 *
 * SCOPE CUT (refused, not silently wrong): saltation (ideal-switch/event)
 * corrections and the free-running (trigger=) period-variation projection
 * are NOT transposed by this engine -- see "what was skipped" in the card
 * #504 report. Detected via PSSsaltEventsSeen()>0 (after the Phi pass) or
 * trig_idx>=0 (free-running caller) and returned as PAC_ADJOINT_REFUSED
 * (2) so pac_td_pxf() falls back to the v1 per-source engine transparently
 * -- correctness is preserved on every fixture (pxf_r2_switched_rc.cir,
 * the one test_pxf.py fixture with real saltation events, is answered by
 * v1 even when `adjoint=1`/auto picks this engine, so the done-when's "==
 * v1" is trivially true there by construction, not because saltation was
 * actually transposed). */

#define PAC_ADJOINT_OK       0
#define PAC_ADJOINT_FAIL     1
#define PAC_ADJOINT_REFUSED  2

/* card #504: `.pxf` without an explicit `adjoint=` picks the adjoint
 * engine automatically once the requested source count exceeds this --
 * "default when K is above a small number" per the card's own brief. A
 * small K is cheap enough either way that the v1 per-source engine (no
 * new machinery in the request path) stays the default. */
#define PXF_ADJOINT_AUTO_THRESHOLD 4

/* Rebuild J_k = G(x_k)+ag0_k*C(x_k) into ckt->CKTmatrix (factored, ready
 * for SMPtSolve) from one recorded trajectory step -- the backward-sweep
 * analogue of pac_refactor_at_accepted(), but standalone (no live
 * transient to protect -- caller saves/restores the handful of ckt-global
 * scalars this touches, see pac_engine_run_adjoint()). */
/* card #544: operating point of every device at the accepted point x_k of a stored step, for the noise evaluation
 * (DEVnoise reads instance values -- gm, gds, ... -- and CKTstate0: diode/BJT currents, vds/vgs).  The backward sweep never
 * advances the device state, so this restores the converged forward CKTstate0 (when it was stored) and loads once at x_k with
 * the integration coefficients ZERO: the capacitive (displacement) currents the transient load adds to the junction current
 * state (NIintegrate with stale CKTstate1 history gave 10^6-fold currents on the first steps and ~1 % elsewhere) vanish, which is
 * the DC junction current plain `.noise` uses for the shot/flicker density.  CKTrhs/CKTrhsOld/CKTag/mode are restored. */
static int
pac_noise_warm_load(CKTcircuit *ckt, const PACtrajStep *step)
{
    int neq = SMPmatSize(ckt->CKTmatrix), it, err = 0;
    double *save_rhs = TMALLOC(double, neq + 1), *save_rhsold = TMALLOC(double, neq + 1);
    double sag0 = ckt->CKTag[0], sag1 = ckt->CKTag[1];
    long save_mode = ckt->CKTmode;
    int save_bypass = ckt->CKTbypass, save_noncon = ckt->CKTnoncon;
    if (getenv("PAC_NOISE_STALE_STATE")) {   /* test only (red proof): the pre-#544 behaviour, device values of whatever load ran last */
        if (save_rhs) tfree(save_rhs);
        if (save_rhsold) tfree(save_rhsold);
        return 0;
    }
    if (!save_rhs || !save_rhsold) {
        if (save_rhs) tfree(save_rhs);
        if (save_rhsold) tfree(save_rhsold);
        return 1;
    }
    memcpy(save_rhs, ckt->CKTrhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(save_rhsold, ckt->CKTrhsOld, (size_t) (neq + 1) * sizeof(double));
    memcpy(ckt->CKTrhsOld, step->x, (size_t) (neq + 1) * sizeof(double));
    if (step->st0 && step->nst0 == ckt->CKTnumStates)
        memcpy(ckt->CKTstate0, step->st0, (size_t) step->nst0 * sizeof(double));
    ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODETRAN | MODEINITFLOAT;
    ckt->CKTbypass = 0;
    ckt->CKTag[0] = 0.0; ckt->CKTag[1] = 0.0;
    for (it = 0; it < 80 && !err; it++) {
        ckt->CKTnoncon = 0;
        err = CKTload(ckt);
        if (!ckt->CKTnoncon) break;
    }
    ckt->CKTag[0] = sag0; ckt->CKTag[1] = sag1;
    ckt->CKTmode = save_mode;
    ckt->CKTbypass = save_bypass;
    ckt->CKTnoncon = save_noncon;
    memcpy(ckt->CKTrhs, save_rhs, (size_t) (neq + 1) * sizeof(double));
    memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
    tfree(save_rhs); tfree(save_rhsold);
    return err ? 1 : 0;
}

static int
pac_adj_refactor_at(CKTcircuit *ckt, const PACtrajStep *step, int devparams)
{
    int neq = SMPmatSize(ckt->CKTmatrix);

    memcpy(ckt->CKTrhsOld, step->x, (size_t) (neq + 1) * sizeof(double));
    ckt->CKTtime  = step->t;
    ckt->CKTdelta = step->h;
    ckt->CKTag[0] = step->ag0;
    ckt->CKTag[1] = step->ag1;
    if (devparams && step->nJ > 0 && pac_noise_warm_load(ckt, step))
        return 1;
    if (step->nJ > 0) {      /* card #541: the forward pass's own J_k */
        if (SMPsetElements(ckt->CKTmatrix, step->nJ, step->Jri, step->Jci, step->Jv))
            return 1;
        return SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin) ? 1 : 0;
    }
    if (CKTload(ckt))
        return 1;
    return SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin) ? 1 : 0;
}

/* card #542: generic-history analogue: x_k from the trajectory, exact J_k factored and the charge matrix
 * C_k in g_pac_C (see pac_gen_loads). */
static int
pac_adj_refactor_generic(CKTcircuit *ckt, const PACtrajStep *step, int devparams)
{
    int neq = SMPmatSize(ckt->CKTmatrix);
    memcpy(ckt->CKTrhs, step->x, (size_t) (neq + 1) * sizeof(double));
    ckt->CKTtime  = step->t;
    ckt->CKTdelta = step->h;
    ckt->CKTag[0] = step->ag0;
    ckt->CKTag[1] = step->ag1;
    if (devparams && pac_noise_warm_load(ckt, step))
        return 1;
    if (pac_trip_fit(&g_pac_C, step->nC)) return 1;
    memcpy(g_pac_C.ri, step->Cri, (size_t) step->nC * sizeof(int));
    memcpy(g_pac_C.ci, step->Cci, (size_t) step->nC * sizeof(int));
    memcpy(g_pac_C.v, step->Cv, (size_t) step->nC * sizeof(double));
    g_pac_C.n = step->nC;
    if (SMPsetElements(ckt->CKTmatrix, step->nJ, step->Jri, step->Jci, step->Jv))
        return 1;
    return SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin) ? 1 : 0;
}

static int
pac_engine_run_adjoint(CKTcircuit *ckt, PSSan *job, PACstateMap *map, const double *x0, double T,
                        const double *freqs, int nfreq, int outrow,
                        PACsrcRef *srclist, int nsrc_total, const int *use_idx, int nsrc_use,
                        PACc *Y0_out /* [nsrc_use * nfreq] */,
                        int trig_idx)
{
    int n = PAC_NST(map), neq = SMPmatSize(ckt->CKTmatrix);
    int d, e, f, k, s, row, err;
    double *br_all = NULL, *bi_all = NULL;    /* [nsrc_use*(neq+1)] */
    PACc *Phi = NULL, *G = NULL, *Mw = NULL;
    PACc *Lambda = NULL;                       /* [nfreq*(neq+1)] */
    PACc *abarS = NULL, *abarT = NULL;          /* [n*nfreq] */
    double *tdu_re = NULL, *tdu_im = NULL, *scratch = NULL;   /* [neq+1] */
    double *Cdiag = NULL;
    double save_ag0 = 0.0, save_ag1 = 0.0, save_time = 0.0, save_delta = 0.0;
    double *save_rhsold = NULL;
    int have_ckt_saved = 0;

    /* card #557: a free-running (trigger=) orbit is transposed like a driven one: the forward operator on it is the fixed-period
     * LPTV engine on the polished orbit (the period-variation projection is opt-in, PAC_FREERUN_PROJ), so its exact transpose is
     * this sweep.  Only the opt-in projection (not transposed) keeps the refusal. */
#ifdef XSPICE
    if (g_hid_nh > 0) {   /* card #557: the transposed sweep does not carry the hidden states of stateful code models */
        fprintf(stderr, "PXF-TD-ADJ: stateful XSPICE code model in the circuit: the adjoint sweep does not transpose its hidden states (card #557).\n");
        return PAC_ADJOINT_REFUSED;
    }
#endif
    if (trig_idx >= 0 && getenv("PAC_FREERUN_PROJ") && !getenv("PAC_FREERUN_NO_PROJ"))
        return PAC_ADJOINT_REFUSED;
    if (n <= 0 || nfreq <= 0 || nsrc_use <= 0)
        return PAC_ADJOINT_FAIL;

    /* ---- capture every requested source's isolated unit stamp FIRST,
     * before this function's own pac_reset()/pac_run_period() calls touch
     * ckt any further -- SAME timing v1's per-source loop uses (capture,
     * then propagate), just batched over every source up front instead of
     * interleaved one source at a time. */
    br_all = TMALLOC(double, (size_t) nsrc_use * (neq + 1));
    bi_all = TMALLOC(double, (size_t) nsrc_use * (neq + 1));
    if (!br_all || !bi_all) { err = PAC_ADJOINT_FAIL; goto done; }
    memset(br_all, 0, (size_t) nsrc_use * (neq + 1) * sizeof(double));
    memset(bi_all, 0, (size_t) nsrc_use * (neq + 1) * sizeof(double));
    pac_state_set_ic(map, x0);
    for (s = 0; s < nsrc_use; s++) {
        if (pac_capture_source_b0(ckt, srclist, nsrc_total, use_idx[s],
                                   br_all + (size_t) s * (neq + 1),
                                   bi_all + (size_t) s * (neq + 1), neq)) {
            fprintf(stderr, "PXF-TD-ADJ: failed to isolate source '%s'\n",
                    (char *) srclist[use_idx[s]].name);
            err = PAC_ADJOINT_FAIL; goto done;
        }
    }

    /* ---- pass A: Phi (n x n) + trajectory recording, ONE period run,
     * homogeneous columns only (b=0), no demod. */
    {
        double *omega1 = TMALLOC(double, n);
        int *useforce1 = TMALLOC(int, n);
        PACc *ic0_1 = TMALLOC(PACc, (size_t) n * n);
        if (!omega1 || !useforce1 || !ic0_1) {
            if (omega1) tfree(omega1); if (useforce1) tfree(useforce1); if (ic0_1) tfree(ic0_1);
            err = PAC_ADJOINT_FAIL; goto done;
        }
        memset(ic0_1, 0, (size_t) n * n * sizeof(PACc));
        for (d = 0; d < n; d++) {
            omega1[d] = 0.0;
            useforce1[d] = 0;
            for (e = 0; e < n; e++) {
                PACc v = { (d == e) ? 1.0 : 0.0, 0.0 };
                ic0_1[(size_t) e * n + d] = v;
            }
        }
        pac_state_set_ic(map, x0);
        pac_traj_reset(neq);
        if (pac_reset(ckt, map, n, NULL, NULL, omega1, useforce1, ic0_1, 0, T)) {
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_traj_free();
            err = PAC_ADJOINT_FAIL; goto done;
        }
        PSSpacActive = 1;
        g_pac_force_exact_jac = 1;
        err = pac_run_period(ckt, job, T, "TD-shoot ADJOINT pass A (Phi + trajectory)");
        g_pac_force_exact_jac = 0;
        PSSpacActive = 0;
        g_pac_traj_active = 0;   /* stop recording -- buffer itself kept */
        pac_traj_report("PXF-TD");

        if (err || g_pac_failed || g_pac_traj_failed) {
            int memrefuse = (!err && !g_pac_failed && g_pac_traj_failed && map->generic);
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_free(); pac_traj_free();
            err = memrefuse ? PAC_ADJOINT_REFUSED : PAC_ADJOINT_FAIL; goto done;
        }

        /* saltation/free-running scope cut (see header) -- refuse
         * gracefully, caller falls back to v1. */
        if (PSSsaltActive && PSSsaltEventsSeen() > 0) {
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_free(); pac_traj_free();
            err = PAC_ADJOINT_REFUSED; goto done;
        }

        Phi = TMALLOC(PACc, (size_t) n * n);
        if (!Phi) {
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_free(); pac_traj_free();
            err = PAC_ADJOINT_FAIL; goto done;
        }
        for (d = 0; d < n; d++) {
            int row_d = pac_map_row(map, d);
            for (e = 0; e < n; e++)
                Phi[(size_t) d * n + e] = pac_get(g_pac_Cprev, row_d, n, e);
        }
        tfree(omega1); tfree(useforce1); tfree(ic0_1);
        pac_free();   /* g_pac_traj_* is independent, survives this */
    }

    if (g_pac_traj_n < 1) {
        tfree(Phi); pac_traj_free();
        err = PAC_ADJOINT_FAIL; goto done;
    }

    /* ---- pass B: G(w)[d][f] for every (homogeneous column d, frequency
     * f), batched M=n*nfreq, demod active. */
    G = TMALLOC(PACc, (size_t) n * nfreq);
    if (!G) { tfree(Phi); pac_traj_free(); err = PAC_ADJOINT_FAIL; goto done; }
    {
        int M2 = n * nfreq;
        double *omega2 = TMALLOC(double, M2);
        int *useforce2 = TMALLOC(int, M2);
        PACc *ic0_2 = TMALLOC(PACc, (size_t) n * M2);
        if (!omega2 || !useforce2 || !ic0_2) {
            if (omega2) tfree(omega2); if (useforce2) tfree(useforce2); if (ic0_2) tfree(ic0_2);
            tfree(Phi); tfree(G); pac_traj_free();
            err = PAC_ADJOINT_FAIL; goto done;
        }
        memset(ic0_2, 0, (size_t) n * M2 * sizeof(PACc));
        for (d = 0; d < n; d++)
            for (f = 0; f < nfreq; f++) {
                int col = d * nfreq + f;
                omega2[col] = 2.0 * M_PI * freqs[f];
                useforce2[col] = 0;
                { PACc v = { 1.0, 0.0 }; ic0_2[(size_t) d * M2 + col] = v; }
            }
        pac_state_set_ic(map, x0);
        if (pac_reset(ckt, map, M2, NULL, NULL, omega2, useforce2, ic0_2, 1, T)) {
            tfree(omega2); tfree(useforce2); tfree(ic0_2);
            tfree(Phi); tfree(G); pac_traj_free();
            err = PAC_ADJOINT_FAIL; goto done;
        }
        PSSpacActive = 1;
        g_pac_force_exact_jac = 1;
        err = pac_run_period(ckt, job, T, "TD-shoot ADJOINT pass B (G(w), homogeneous+demod)");
        g_pac_force_exact_jac = 0;
        PSSpacActive = 0;
        if (err || g_pac_failed) {
            tfree(omega2); tfree(useforce2); tfree(ic0_2);
            pac_free(); tfree(Phi); tfree(G); pac_traj_free();
            err = PAC_ADJOINT_FAIL; goto done;
        }
        pac_finish_demod_wraparound();
        {
            double invT = 1.0 / T;
            for (d = 0; d < n; d++)
                for (f = 0; f < nfreq; f++) {
                    PACc v = g_pac_demod[(size_t) (outrow - 1) * M2 + (d * nfreq + f)];
                    G[(size_t) d * nfreq + f] = pc_scale(v, invT);
                }
        }
        pac_free();
        tfree(omega2); tfree(useforce2); tfree(ic0_2);
    }

    /* ---- M(w) = ((e^{jwT}I - Phi)^T)^{-1} G(w), one small dense complex
     * solve per frequency using the EXPLICIT TRANSPOSE of the boundary
     * matrix pac_engine_run() itself builds ((e^{jwT}I-Phi), not
     * transposed there). */
    Mw = TMALLOC(PACc, (size_t) n * nfreq);
    if (!Mw) { tfree(Phi); tfree(G); pac_traj_free(); err = PAC_ADJOINT_FAIL; goto done; }
    {
        PACc *AT = TMALLOC(PACc, (size_t) n * n);
        PACc *rhs = TMALLOC(PACc, n);
        PACc *sol = TMALLOC(PACc, n);
        if (!AT || !rhs || !sol) {
            if (AT) tfree(AT); if (rhs) tfree(rhs); if (sol) tfree(sol);
            tfree(Phi); tfree(G); tfree(Mw); pac_traj_free();
            err = PAC_ADJOINT_FAIL; goto done;
        }
        for (f = 0; f < nfreq; f++) {
            double w = 2.0 * M_PI * freqs[f];
            double wT = w * T;
            PACc ephase = { cos(wT), sin(wT) };
            for (d = 0; d < n; d++)
                for (e = 0; e < n; e++) {
                    /* AT[d][e] = A[e][d] = (e==d?ephase:0) - Phi[e][d] */
                    PACc v = pc_scale(Phi[(size_t) e * n + d], -1.0);
                    if (e == d) v = pc_add(v, ephase);
                    AT[(size_t) d * n + e] = v;
                }
            for (d = 0; d < n; d++) rhs[d] = G[(size_t) d * nfreq + f];
            if (pac_complex_solve(AT, rhs, n, sol)) {
                fprintf(stderr, "PXF-TD-ADJ: WARNING -- transposed boundary solve "
                        "singular at f=%.6g Hz; result at this frequency skipped.\n",
                        freqs[f]);
                for (d = 0; d < n; d++) { PACc nan = { strtod("NaN", NULL), strtod("NaN", NULL) }; sol[d] = nan; }
            }
            for (d = 0; d < n; d++)
                Mw[(size_t) d * nfreq + f] = sol[d];
        }
        tfree(AT); tfree(rhs); tfree(sol);
    }
    if (getenv("PXF_ADJ_DEBUG")) {
        fprintf(stderr, "PXF-ADJ-DEBUG: n=%d nfreq=%d traj_n=%d\n", n, nfreq, g_pac_traj_n);
        for (d = 0; d < n; d++) {
            for (e = 0; e < n; e++)
                fprintf(stderr, "  Phi[%d][%d]=%.15g%+.15gi", d, e,
                        Phi[(size_t) d * n + e].re, Phi[(size_t) d * n + e].im);
            fprintf(stderr, "\n");
        }
        for (d = 0; d < n; d++)
            fprintf(stderr, "  G[%d][f0]=%.6g%+.6gi  Mw[%d][f0]=%.6g%+.6gi\n", d,
                    G[(size_t) d * nfreq + 0].re, G[(size_t) d * nfreq + 0].im,
                    d, Mw[(size_t) d * nfreq + 0].re, Mw[(size_t) d * nfreq + 0].im);
    }
    tfree(Phi); tfree(G);

    /* ---- backward adjoint sweep over the recorded trajectory, batched
     * over every frequency (ONE CKTload()+SMPluFac() refactor per step,
     * reused for all nfreq transposed solves at that step -- mirrors the
     * forward passes' own "one factorization per accepted step, many
     * columns" reuse). */
    Lambda = TMALLOC(PACc, (size_t) nfreq * (neq + 1));
    abarS  = TMALLOC(PACc, (size_t) (map->generic ? neq + 1 : n) * nfreq);
    abarT  = TMALLOC(PACc, (size_t) (map->generic ? neq + 1 : n) * nfreq);
    tdu_re = TMALLOC(double, (size_t) (neq + 1));
    tdu_im = TMALLOC(double, (size_t) (neq + 1));
    scratch = TMALLOC(double, (size_t) (neq + 1));
    Cdiag  = TMALLOC(double, n);
    save_rhsold = TMALLOC(double, (size_t) (neq + 1));
    if (!Lambda || !abarS || !abarT || !tdu_re || !tdu_im || !scratch || !Cdiag || !save_rhsold) {
        err = PAC_ADJOINT_FAIL; goto done;
    }
    memset(Lambda, 0, (size_t) nfreq * (neq + 1) * sizeof(PACc));
    memset(abarS, 0, (size_t) (map->generic ? neq + 1 : n) * nfreq * sizeof(PACc));
    memset(abarT, 0, (size_t) (map->generic ? neq + 1 : n) * nfreq * sizeof(PACc));

    memcpy(save_rhsold, ckt->CKTrhsOld, (size_t) (neq + 1) * sizeof(double));
    save_time = ckt->CKTtime; save_delta = ckt->CKTdelta;
    save_ag0 = ckt->CKTag[0]; save_ag1 = ckt->CKTag[1];
    have_ckt_saved = 1;

    for (k = g_pac_traj_n - 1; k >= 0; k--) {
        const PACtrajStep *step = &g_pac_traj[k];
        int is_last = (k == g_pac_traj_n - 1);
        double h_next = is_last ? g_pac_traj[0].h : g_pac_traj[k + 1].h;
        double ag1eff = (step->order == 2) ? step->ag1 : 0.0;

        if (map->generic ? pac_adj_refactor_generic(ckt, step, 0) : pac_adj_refactor_at(ckt, step, 0)) {
            fprintf(stderr, "PXF-TD-ADJ: refactor failed at t=%.9g (backward step %d)\n",
                    step->t, k);
            err = PAC_ADJOINT_FAIL; goto done;
        }

        for (e = 0; !map->generic && e < n; e++) {
            const PACstateRef *ref = &map->refs[e];
            Cdiag[e] = ref->isCap ? step->ag0 * ref->dev.cap->CAPcapac
                                   : step->ag0 * ref->dev.ind->INDinduct;
        }

        for (f = 0; f < nfreq; f++) {
            double w = 2.0 * M_PI * freqs[f];
            PACc ph1 = { cos(w * step->t), sin(w * step->t) };       /* e^{+jwt_k} */
            PACc ph2 = pc_conjexp(w * step->t);                       /* e^{-jwt_k} */
            /* the (1/T) normalisation pac_engine_run()'s own pass 2 applies
             * ONCE, at the very end, to the whole raw demod sum (Y0_out =
             * pc_scale(g_pac_demod[...], invT)) -- G(w) above already has
             * it baked in (it is read off exactly that way), so the SAME
             * raw-weight formula here needs it applied explicitly, or this
             * term is 1/T too large relative to the M(w) contribution. */
            PACc coef = pc_scale(ph2, 0.5 * (step->h + h_next) / T);

            memset(tdu_re, 0, (size_t) (neq + 1) * sizeof(double));
            memset(tdu_im, 0, (size_t) (neq + 1) * sizeof(double));
            tdu_re[outrow] += coef.re; tdu_im[outrow] += coef.im;

            if (is_last) {
                /* card #504 bug fix (found validating pac_r3_ctrl2out.cir,
                 * which -- unlike pxf_lti.cir/pxf_r2_switched_rc.cir -- has
                 * a capacitor whose negNode is NOT ground, C1f between
                 * compin and opout): M(w) is dual to z(w) = E_read(u_N(b))
                 * (pac_engine_run()'s own z[d][f] = pac_get(Cprev, row_d,
                 * ...), row_d = pac_state_row() = posNode/branch ONLY --
                 * the SAME convention Phi and G(w) above are read off with,
                 * NOT the posNode-minus-negNode difference the INTERNAL
                 * per-step recursion's rhs/T_k formulas use just below and
                 * in PSSpacRecordStep() itself. The two conventions
                 * coincide whenever a capacitor's negNode is ground (every
                 * cap in pxf_lti.cir/pxf_r2_switched_rc.cir -- why this bug
                 * passed those two fixtures) but NOT in general -- z(w)
                 * itself never reads the negNode, so M(w)'s own scatter
                 * here must not either. */
                for (e = 0; e < n; e++) {
                    PACc Mv = Mw[(size_t) e * nfreq + f];
                    int row = pac_map_row(map, e);
                    if (row > 0) { tdu_re[row] += Mv.re; tdu_im[row] += Mv.im; }
                }
            }

            if (map->generic) {
                /* card #542: z += C_k^T (alpha + ag0 beta) over every row (alpha/beta: adjoints of Q_k, T_k
                 * from the later steps, see the generic recursion in the design note) */
                int ti;
                for (ti = 0; ti < g_pac_C.n; ti++) {
                    int r = g_pac_C.ri[ti], c = g_pac_C.ci[ti];
                    PACc al, be;
                    double v = g_pac_C.v[ti];
                    if (r <= 0 || c <= 0) continue;
                    al = abarS[(size_t) r * nfreq + f];
                    be = abarT[(size_t) r * nfreq + f];
                    tdu_re[c] += v * (al.re + step->ag0 * be.re);
                    tdu_im[c] += v * (al.im + step->ag0 * be.im);
                }
            }
            for (e = 0; !map->generic && e < n; e++) {
                const PACstateRef *ref = &map->refs[e];
                PACc abT = abarT[(size_t) e * nfreq + f];
                PACc abS_ext = abarS[(size_t) e * nfreq + f];
                PACc abS_tot = pc_add(abS_ext, pc_scale(abT, Cdiag[e]));
                if (ref->isCap) {
                    int p = ref->dev.cap->CAPposNode, ngc = ref->dev.cap->CAPnegNode;
                    if (p > 0)   { tdu_re[p] += abS_tot.re; tdu_im[p] += abS_tot.im; }
                    if (ngc > 0) { tdu_re[ngc] -= abS_tot.re; tdu_im[ngc] -= abS_tot.im; }
                } else {
                    int ibr = ref->dev.ind->INDbrEq;
                    tdu_re[ibr] += abS_tot.re; tdu_im[ibr] += abS_tot.im;
                }
            }

            SMPtSolve(ckt->CKTmatrix, tdu_re, scratch);
            SMPtSolve(ckt->CKTmatrix, tdu_im, scratch);
            /* tdu_re/tdu_im now hold nu_k (this step's, this frequency's
             * transposed-solve result), for EVERY row 1..neq. */

            for (row = 1; row <= neq; row++) {
                PACc nu = { tdu_re[row], tdu_im[row] };
                PACc contrib = pc_mul(nu, ph1);
                size_t idx = (size_t) f * (neq + 1) + row;
                Lambda[idx] = pc_add(Lambda[idx], contrib);
            }

            if (map->generic) {
                /* gamma = nu - beta_old; alpha_{k-1} = ag0 gamma; beta_{k-1} = a1 gamma */
                for (row = 1; row <= neq; row++) {
                    PACc be = abarT[(size_t) row * nfreq + f];
                    PACc gm = { tdu_re[row] - be.re, tdu_im[row] - be.im };
                    abarS[(size_t) row * nfreq + f] = pc_scale(gm, step->ag0);
                    abarT[(size_t) row * nfreq + f] = pc_scale(gm, ag1eff);
                }
            }
            for (e = 0; !map->generic && e < n; e++) {
                const PACstateRef *ref = &map->refs[e];
                PACc beta, gamma;
                if (ref->isCap) {
                    int p = ref->dev.cap->CAPposNode, ngc = ref->dev.cap->CAPnegNode;
                    PACc nup = (p > 0)   ? (PACc){ tdu_re[p], tdu_im[p] }     : (PACc){ 0.0, 0.0 };
                    PACc nun = (ngc > 0) ? (PACc){ tdu_re[ngc], tdu_im[ngc] } : (PACc){ 0.0, 0.0 };
                    PACc diffnu = pc_sub(nup, nun);
                    beta  = pc_scale(diffnu, step->ag0 * ref->dev.cap->CAPcapac);
                    gamma = pc_scale(diffnu, ag1eff);
                } else {
                    int ibr = ref->dev.ind->INDbrEq;
                    PACc nub = { tdu_re[ibr], tdu_im[ibr] };
                    beta  = pc_scale(nub, -step->ag0 * ref->dev.ind->INDinduct);
                    gamma = pc_scale(nub, -ag1eff);
                }
                {
                    PACc abT_old = abarT[(size_t) e * nfreq + f];
                    PACc new_abS = pc_sub(beta, pc_scale(abT_old, Cdiag[e]));
                    PACc new_abT = pc_add(gamma, pc_scale(abT_old, -ag1eff));
                    abarS[(size_t) e * nfreq + f] = new_abS;
                    abarT[(size_t) e * nfreq + f] = new_abT;
                }
            }
        }
    }

    memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
    ckt->CKTtime = save_time; ckt->CKTdelta = save_delta;
    ckt->CKTag[0] = save_ag0; ckt->CKTag[1] = save_ag1;
    have_ckt_saved = 0;

    /* ---- every requested source's PXF: a plain dot product against
     * Lambda(w), b^T Lambda (bilinear, NOT conjugated -- see docs/sens.md's
     * own "Y instead of Y^T" red-proof note; the SAME convention). */
    for (s = 0; s < nsrc_use; s++) {
        const double *br = br_all + (size_t) s * (neq + 1);
        const double *bi = bi_all + (size_t) s * (neq + 1);
        for (f = 0; f < nfreq; f++) {
            double sre = 0.0, sim = 0.0;
            for (row = 1; row <= neq; row++) {
                if (br[row] == 0.0 && bi[row] == 0.0) continue;
                PACc L = Lambda[(size_t) f * (neq + 1) + row];
                sre += br[row] * L.re - bi[row] * L.im;
                sim += br[row] * L.im + bi[row] * L.re;
            }
            Y0_out[(size_t) s * nfreq + f].re = sre;
            Y0_out[(size_t) s * nfreq + f].im = sim;
        }
    }

    tfree(Mw); pac_traj_free();
    err = PAC_ADJOINT_OK;

 done:
    if (have_ckt_saved) {
        memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
        ckt->CKTtime = save_time; ckt->CKTdelta = save_delta;
        ckt->CKTag[0] = save_ag0; ckt->CKTag[1] = save_ag1;
    }
    if (br_all) tfree(br_all);
    if (bi_all) tfree(bi_all);
    if (Lambda) tfree(Lambda);
    if (abarS) tfree(abarS);
    if (abarT) tfree(abarT);
    if (tdu_re) tfree(tdu_re);
    if (tdu_im) tfree(tdu_im);
    if (scratch) tfree(scratch);
    if (Cdiag) tfree(Cdiag);
    if (save_rhsold) tfree(save_rhsold);
    return err;
}

/* ------------------------------------------------------------------ */
/* card #509: `.pnoise ... method=shoot solver=td [dev=1]` -- periodic
 * (cyclostationary) noise on the shooting-PSS orbit, output = the
 * time-averaged PSD at sideband 0 of one output node (Spectre's own
 * `pnoise ... noisetype=timeaverage`).
 *
 * DERIVATION. pac_engine_run_adjoint() above computes, for a fixed output
 * row and frequency w, a single Lambda(w) = Sum_k nu_k(w) * e^{+j*w*t_k},
 * where nu_k(w) is that function's own per-step SMPtSolve() result BEFORE
 * the e^{+jwt_k} phase multiply -- by its own derivation (see its header),
 * Y0(w,b) = b^T Lambda(w) for a source forced IDENTICALLY (the same
 * spatial pattern b, phase-modulated e^{jwt_k}) at every step. Substituting
 * Lambda(w) = Sum_k nu_k*e^{jwt_k} and regrouping,
 *
 *   Y0(w,b) = b^T Sum_k (nu_k * e^{jwt_k}) = Sum_k (b*e^{jwt_k})^T nu_k
 *
 * and b*e^{jwt_k} IS this file's own per-step forced-column convention
 * (the SAME discretisation pac_capture_source_b0's br/bi get phase-
 * modulated through at every step -- see this file's header, "b(t) =
 * b*e^{jwt}") -- so nu_k itself is exactly dY0(w)/d(beta_k), the
 * sensitivity of Y0(w) to a UNIT, UNMODULATED forcing injected into the
 * per-step RHS at step k ONLY: the discrete-time periodic transfer
 * function H(t_k,f) this card's own brief calls for ("a noise current
 * injected over step k has the same stamp as a small-signal source at
 * step k") -- the SAME per-step RHS convention as a stationary AC source,
 * just not phase-modulated across steps (a genuinely LOCALISED injection
 * at ONE step, not a continuous sinusoid).
 *
 * For a noise generator g with a FIXED (bias-independent) terminal
 * incidence pattern b_g (e.g. drain-minus-source for a MOS channel
 * current -- EXACTLY NevalSrc()'s own posNode/negNode convention) and a
 * cyclostationary density S_g(t_k,f) (the device's own DEVnoise routine,
 * evaluated at the step-k operating point), the Telichevesky/Kundert
 * (1996)/Demir-Roychowdhury cyclostationary-noise integral
 *
 *   onoise(f) = (1/T) INT_0^T |H(t,f)|^2 S(t,f) dt
 *
 * discretises, on this engine's OWN periodic-trapezoidal grid, to
 *
 *   onoise(f) = T * Sum_k (1/w_k) * Sum_g |b_g . nu_k(f)|^2 S_g(t_k,f)
 *
 * w_k = 0.5*(h_k+h_next) is the SAME periodic-trapezoidal weight the
 * boundary term's own "coef" below uses for the OUTPUT side -- used here,
 * DIVIDED rather than multiplied in, for the INJECTION side (see the
 * SMPtSolve() call site below, right before dens_k is computed, for the
 * full white-noise-sample-variance derivation of why: in short, nu_k(f) is
 * the response to a DETERMINISTIC unit beta_k, but a noise sample n_k
 * consistent with continuous PSD S_k needs Var(n_k)=S_k/w_k, not S_k, so
 * the w_k lands in the DENOMINATOR of dens_k = S_k*|nu_k|^2, not as an
 * extra multiplicative quadrature weight the way it is for the output
 * tap). Two other candidates were tried and found WRONG by this card's own
 * pxf_lti.cir-vs-plain-.noise check before this one: "no extra factor at
 * all" (under-reports by ~N, the accepted-step count) and "multiply by
 * w_k/T" (the output tap's OWN weight, applied here too by a naive direct
 * analogy -- under-reports by ~N^2). "T * dens_k/w_k" matches plain .noise
 * to 6e-6 relative at Points=2000 -- see "Validation" in docs/pss-pac.md's
 * `.pnoise` section for the full numbers (all three candidates, several
 * grid resolutions) and the convergence-with-Points table.
 *
 * Sum_g |b_g.nu_k(f)|^2 S_g(t_k,f) is EXACTLY what NevalSrc()/DEVnoise()
 * already compute when ckt->CKTrhs/CKTirhs are loaded with nu_k(f) at the
 * step-k bias (NevalSrc: gain = |CKTrhs[node1]-CKTrhs[node2]|^2, *noise =
 * gain*S(bias)) -- so this engine needs NO new device-noise physics: it
 * reuses pac_engine_run_adjoint()'s OWN backward sweep (Phi/G(w)/M(w)
 * passes, and the per-step abarS/abarT/Cdiag/beta/gamma recursion --
 * copied verbatim below, NOT shared by a refactor, so this card touches
 * ZERO lines of the already-validated .pxf adjoint engine) and, at each
 * (k,f), loads nu_k into ckt->CKTrhs/CKTirhs with the bias ALREADY at x_k
 * (pac_adj_refactor_at() -- the SAME CKTload() that rebuilds J_k also
 * refreshes every device's own bias-dependent noise parameters, e.g.
 * RESconduct/REScurrent/gm) and calls the ordinary
 * DEVnoise(N_DENS,N_CALC,...) loop -- EXACTLY the idiom dcpss.c's own
 * frequency-domain pnoise_sweep()/QPnoiseAnalyze() already use for their
 * own (different) adjoint vector.
 *
 * SCOPE (v1, card #509): restricted to resistor thermal, diode shot, and
 * MOS1/BSIM4 channel thermal+flicker noise -- the only DEVnoise routines
 * this derivation has been checked against (pac_noise_check_devices()
 * below refuses, loudly, rather than silently omitting or mis-scaling,
 * any OTHER device type with noise-producing instances in the circuit).
 * Saltation (ideal-switch) events and free-running orbits are refused
 * (PAC_ADJOINT_REFUSED) exactly like `.pxf ... adjoint=1` -- see that
 * function's own "SCOPE CUT" -- with NO v1/non-adjoint fallback (unlike
 * .pxf: there is no meaningful "one pass per source" analogue for noise,
 * which comes continuously from every device, not a short list of named
 * sources). Input-referred (inoise) is NOT computed by this engine (v1
 * scope; the existing frequency-domain `.pnoise` still reports it via
 * PnInSrc). See docs/pss-pac.md's `.pnoise` section for the full
 * derivation, validation, and what was skipped. */


/* card #556: folded flicker (1/f) noise on a time-varying orbit.
 *
 * A device noise source n(t) = g(t) w(t), w stationary with a power-law density S_w(f') = A / |f'|^ef, enters a periodic
 * system.  The time-averaged output density at f is  Sum_q |G_q(f)|^2 S_w(f - q fs)  where G_q are the Fourier coefficients
 * of the periodic function a(t) = g(t) H(t,f) (H = response at f to a unit source at t).  The recursion above evaluates the
 * device density at the OUTPUT frequency f for every sideband, i.e. Sum_q |G_q|^2 S_w(f): exact for white noise (Parseval) and
 * for flicker when a(t) has no harmonics (LTI bias, slow modulation), but it overstates a fast-modulated flicker source by
 * S_w(f)/S_w(f - q fs).  Per flicker generator c and output frequency f this stores z_k = u_k sqrt(S_k) e^{+j w t_k} (u_k =
 * complex gain of the generator's two nodes in nu_k, S_k = its density at f; |z_k|^2 = the density the sweep accumulated), the
 * exact Parseval pairing  Sum_q |Sum_k z_k sinc(q ws w_k/2) e^{-j q ws tc_k}|^2 = T Sum_k |z_k|^2 / w_k  holds for the
 * piecewise-constant (zero-order-hold, interval w_k centred at tc_k) a(t) the engine's own discretisation represents, and
 * the correction  Delta = Sum_q |G_q|^2 (rho_q - 1),  rho_q = (f/|f - q fs|)^ef, is added to the accumulated density.  The
 * sum runs over |q| <= Qmax (env PNOISE_FOLD_QMAX, default 2048) or until the unsummed part of the Parseval energy is
 * below 1e-9; the rest is added with the density of a q^-2 tail, (f/(Q fs))^ef / (1+ef).  White and thermal sources are not
 * touched (no flicker slot, or no flicker current: the correction is skipped and the result is bit-identical). */
extern int g_nev_gain_cap, g_nev_gain_n, g_nev_gain_max;
extern double *g_nev_gain_re, *g_nev_gain_im;

typedef struct {
    int nfl, nfreq, N;
    int *slot, *tot;           /* [nfl] output-vector index of the flicker slot and of its instance total (-1: none) */
    double *ef;                /* [nfl] exponent: S ~ f^-ef (from the density at f and 2f) */
    int *ef_set;
    int nonpow;                /* the density is not a power law in f for some generator: folding not applied */
    int warned_sw;
    double *zr, *zi;           /* [(c*nfreq + f)*N + k] */
    double *E;                 /* [c*nfreq + f] energy accumulated by the sweep: Sum_k (T/w_k) |z_k|^2 */
    double *tc, *wk;           /* [N] */
    double *capr, *capi;       /* [nfl] N_GAIN capture buffer */
} PACfold;

static void
pac_fold_free(PACfold *F)
{
    if (!F) return;
    if (F->slot) tfree(F->slot);
    if (F->tot) tfree(F->tot);
    if (F->ef) tfree(F->ef);
    if (F->ef_set) tfree(F->ef_set);
    if (F->zr) tfree(F->zr);
    if (F->zi) tfree(F->zi);
    if (F->E) tfree(F->E);
    if (F->tc) tfree(F->tc);
    if (F->wk) tfree(F->wk);
    if (F->capr) tfree(F->capr);
    if (F->capi) tfree(F->capi);
    memset(F, 0, sizeof(*F));
}

/* Apply the correction.  onoise_out[f] and (if not NULL) gen_dens[slot*nfreq + f] are corrected in place. */
static void
pac_fold_apply(PACfold *F, const double *freqs, double T, double *onoise_out, double *gen_dens)
{
    int N = F->N, nfreq = F->nfreq, nfl = F->nfl;
    int np = nfl * nfreq, p, k, q, qmax = 2048;
    double ws = 2.0 * M_PI / T, fs = 1.0 / T;
    double *pr, *pi, *rr, *ri, *cr, *ci, *sr, *si, *bw, *sinc;
    double *head, *sum, *clampc;
    char *active;
    int nactive = 0, qlast = 0, warned_clamp = 0;
    const double fl_floor = 1e-3 * fs;      /* the 1/f density of a sideband at |f'| < 1e-3 fs is that of 1e-3 fs (S(0) is infinite) */
    const char *qe = getenv("PNOISE_FOLD_QMAX");
    int dbg = getenv("PNOISE_FOLD_DEBUG") != NULL;

    if (qe && atoi(qe) > 0) qmax = atoi(qe);
    pr = TMALLOC(double, N); pi = TMALLOC(double, N);
    rr = TMALLOC(double, N); ri = TMALLOC(double, N);
    cr = TMALLOC(double, N); ci = TMALLOC(double, N);
    sr = TMALLOC(double, N); si = TMALLOC(double, N);
    bw = TMALLOC(double, N); sinc = TMALLOC(double, N);
    head = TMALLOC(double, np); sum = TMALLOC(double, np); clampc = TMALLOC(double, np);
    active = TMALLOC(char, np);
    if (!pr || !pi || !rr || !ri || !cr || !ci || !sr || !si || !bw || !sinc || !head || !sum || !clampc || !active) goto out;

    for (k = 0; k < N; k++) {
        pr[k] = 1.0; pi[k] = 0.0;
        rr[k] = cos(ws * F->tc[k]); ri[k] = -sin(ws * F->tc[k]);     /* e^{-j ws tc_k} */
        bw[k] = 0.5 * ws * F->wk[k];
        cr[k] = 1.0; ci[k] = 0.0;                                     /* e^{+j q bw_k} */
        sr[k] = cos(bw[k]); si[k] = sin(bw[k]);
    }
    for (p = 0; p < np; p++) {
        int c = p / nfreq;
        const double *zr = F->zr + (size_t) p * N, *zi = F->zi + (size_t) p * N;
        double Gr = 0.0, Gi = 0.0;
        active[p] = 0; head[p] = sum[p] = clampc[p] = 0.0;
        if (!(F->E[p] > 0.0) || !(F->ef[c] > 1e-12)) continue;
        for (k = 0; k < N; k++) { Gr += zr[k]; Gi += zi[k]; }
        head[p] = sum[p] = Gr * Gr + Gi * Gi;
        active[p] = 1; nactive++;
    }
    for (q = 1; q <= qmax && nactive > 0; q++) {
        for (k = 0; k < N; k++) {
            double t;
            t = pr[k] * rr[k] - pi[k] * ri[k]; pi[k] = pr[k] * ri[k] + pi[k] * rr[k]; pr[k] = t;
            t = cr[k] * sr[k] - ci[k] * si[k]; ci[k] = cr[k] * si[k] + ci[k] * sr[k]; cr[k] = t;
            sinc[k] = (bw[k] > 0.0) ? ci[k] / ((double) q * bw[k]) : 1.0;
            if (q % 256 == 0) {                       /* renormalise the rotating phasors */
                double m = 1.0 / sqrt(pr[k] * pr[k] + pi[k] * pi[k]);
                pr[k] *= m; pi[k] *= m;
                m = 1.0 / sqrt(cr[k] * cr[k] + ci[k] * ci[k]);
                cr[k] *= m; ci[k] *= m;
            }
        }
        for (p = 0; p < np; p++) {
            int c = p / nfreq, f = p % nfreq;
            const double *zr, *zi;
            double Pr = 0.0, Pi = 0.0, Mr = 0.0, Mi = 0.0, gp, gm, fa, ef, rp, rm;
            if (!active[p]) continue;
            zr = F->zr + (size_t) p * N; zi = F->zi + (size_t) p * N;
            for (k = 0; k < N; k++) {
                double ar = zr[k] * sinc[k], ai = zi[k] * sinc[k];
                Pr += ar * pr[k] - ai * pi[k];  Pi += ar * pi[k] + ai * pr[k];      /* a e^{-j q ws tc} */
                Mr += ar * pr[k] + ai * pi[k];  Mi += ai * pr[k] - ar * pi[k];      /* a e^{+j q ws tc} */
            }
            gp = Pr * Pr + Pi * Pi;  gm = Mr * Mr + Mi * Mi;
            ef = F->ef[c];  fa = freqs[f];
            {
                double dp = fabs(fa - q * fs), dm = fabs(fa + q * fs);
                if (dp < fl_floor) { dp = fl_floor; clampc[p] += gp * pow(fa / dp, ef); }
                if (dm < fl_floor) { dm = fl_floor; clampc[p] += gm * pow(fa / dm, ef); }
                rp = pow(fa / dp, ef);  rm = pow(fa / dm, ef);
            }
            sum[p] += gp * rp + gm * rm;
            head[p] += gp + gm;
            if (F->E[p] - head[p] <= 1e-9 * F->E[p]) { active[p] = 0; nactive--; }
        }
        qlast = q;
    }
    for (p = 0; p < np; p++) {
        int c = p / nfreq, f = p % nfreq;
        double tail, nw, delta;
        if (!(F->E[p] > 0.0) || !(F->ef[c] > 1e-12)) continue;
        tail = F->E[p] - head[p];
        if (tail < 0.0) tail = 0.0;
        nw = sum[p] + (active[p] ? tail * pow(freqs[f] / ((qlast + 1) * fs), F->ef[c]) / (1.0 + F->ef[c]) : 0.0);
        delta = nw - F->E[p];
        if (clampc[p] > 0.01 * nw && !warned_clamp) {
            fprintf(stderr, "PNOISE-TD: WARNING output frequency %.6g Hz is within 1e-3 fs of a harmonic of fs = %.6g Hz and the "
                    "orbit modulates a flicker source at that harmonic: the 1/f density at zero frequency is infinite, evaluated "
                    "at 1e-3 fs instead (%.0f %% of this value)\n", freqs[f], fs, 100.0 * clampc[p] / nw);
            warned_clamp = 1;
        }
        if (dbg)
            fprintf(stderr, "PNOISE-FOLD: gen %d f=%.6g E=%.6e folded=%.6e ratio=%.6g ef=%.4g head/E=%.9f Q=%d\n",
                    c, freqs[f], F->E[p], nw, nw / F->E[p], F->ef[c], head[p] / F->E[p], qlast);
        onoise_out[f] += delta;
        if (gen_dens) {
            gen_dens[(size_t) F->slot[c] * nfreq + f] += delta;
            if (F->tot[c] >= 0) gen_dens[(size_t) F->tot[c] * nfreq + f] += delta;
        }
    }
out:
    if (pr) tfree(pr); if (pi) tfree(pi); if (rr) tfree(rr); if (ri) tfree(ri);
    if (cr) tfree(cr); if (ci) tfree(ci); if (sr) tfree(sr); if (si) tfree(si);
    if (bw) tfree(bw); if (sinc) tfree(sinc); if (head) tfree(head); if (sum) tfree(sum); if (clampc) tfree(clampc); if (active) tfree(active);
}

/* card #509: restrict this engine to DEVnoise routines this derivation has
 * actually been checked against -- refuse (loudly) rather than silently
 * mis-account any other noisy device type. Matched by FUNCTION POINTER
 * (not device-name string), so an aliased/renamed model table still
 * refuses correctly. */
static int
pac_noise_check_devices(CKTcircuit *ckt)
{
    extern int RESnoise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    extern int DIOnoise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    extern int MOS1noise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    extern int BSIM4noise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    /* card #544: BJT, JFET, MOS3 and BSIM3 use the same NevalSrc()/gain idiom; validated against plain .noise
     * (tests/test_pac_charge.py check_pnoise_dev) before being admitted. */
    extern int BJTnoise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    extern int JFETnoise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    extern int MOS3noise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    extern int BSIM3noise(int, int, GENmodel *, CKTcircuit *, Ndata *, double *);
    int i, bad = 0;

    for (i = 0; i < DEVmaxnum; i++) {
        if (!DEVices[i] || !DEVices[i]->DEVnoise || !ckt->CKThead[i])
            continue;
        if (DEVices[i]->DEVnoise == RESnoise || DEVices[i]->DEVnoise == DIOnoise ||
            DEVices[i]->DEVnoise == MOS1noise || DEVices[i]->DEVnoise == BSIM4noise ||
            DEVices[i]->DEVnoise == BJTnoise || DEVices[i]->DEVnoise == JFETnoise ||
            DEVices[i]->DEVnoise == MOS3noise || DEVices[i]->DEVnoise == BSIM3noise)
            continue;
        fprintf(stderr, "PNOISE-TD: REFUSED -- device type '%s' has noise-producing "
                "instances in this circuit but is outside this engine's v1 scope "
                "(card #509/#544: resistor, diode, BJT, JFET, MOS1/MOS3/BSIM3/BSIM4 only) -- refusing rather than silently omitting or "
                "mis-scaling its noise contribution. See docs/pss-pac.md's .pnoise "
                "section.\n",
                DEVices[i]->DEVpublic.name ? DEVices[i]->DEVpublic.name : "?");
        bad = 1;
    }
    return bad;
}

static int
pac_engine_run_adjoint_noise(CKTcircuit *ckt, PSSan *job, PACstateMap *map,
                              const double *x0, double T,
                              const double *freqs, int nfreq, int outrow,
                              int want_table, int trig_idx,
                              double *onoise_out /* [nfreq], caller-zeroed */,
                              double **gen_dens_out /* OUT [ngen*nfreq], malloc'd, NULL if !want_table/ngen==0 */,
                              IFuid **gen_names_out /* OUT [ngen], malloc'd, NULL if !want_table/ngen==0 */,
                              int *ngen_out)
{
    int n = PAC_NST(map), neq = SMPmatSize(ckt->CKTmatrix);
    int d, e, f, k, row, err;
    PACc *Phi = NULL, *G = NULL, *Mw = NULL;
    double *tdu_re = NULL, *tdu_im = NULL, *scratch = NULL;
    double *Cdiag = NULL;
    PACc *abarS = NULL, *abarT = NULL;
    double save_ag0 = 0.0, save_ag1 = 0.0, save_time = 0.0, save_delta = 0.0;
    double *save_rhsold = NULL;
    int have_ckt_saved = 0, have_job_saved = 0;
    NOISEAN nj;
    Ndata   nd;
    JOB    *oldJob = NULL;
    int     ngen = 0;
    double *gen_dens = NULL;
    IFuid  *gen_names = NULL;
    IFuid  *names_local = NULL;            /* card #556: the slot names (also kept when no table was asked for) */
    PACfold fold;
    int     fold_on = getenv("PNOISE_FLICKER_OLD") == NULL;   /* env: restore the pre-#556 flicker evaluation (red proof) */

    memset(&fold, 0, sizeof(fold));
    memset(&nd, 0, sizeof(nd));
    *gen_dens_out = NULL; *gen_names_out = NULL; *ngen_out = 0;

    if (trig_idx >= 0 && getenv("PAC_FREERUN_PROJ") && !getenv("PAC_FREERUN_NO_PROJ"))
        return PAC_ADJOINT_REFUSED;   /* card #557: only the opt-in period-variation projection is not transposed */
    if (n <= 0 || nfreq <= 0)
        return PAC_ADJOINT_FAIL;

    /* ---- pass A: Phi (n x n) + trajectory recording -- IDENTICAL design
     * to pac_engine_run_adjoint()'s own pass A above (duplicated, not
     * shared -- see this function's header). */
    {
        double *omega1 = TMALLOC(double, n);
        int *useforce1 = TMALLOC(int, n);
        PACc *ic0_1 = TMALLOC(PACc, (size_t) n * n);
        if (!omega1 || !useforce1 || !ic0_1) {
            if (omega1) tfree(omega1); if (useforce1) tfree(useforce1); if (ic0_1) tfree(ic0_1);
            return PAC_ADJOINT_FAIL;
        }
        memset(ic0_1, 0, (size_t) n * n * sizeof(PACc));
        for (d = 0; d < n; d++) {
            omega1[d] = 0.0;
            useforce1[d] = 0;
            for (e = 0; e < n; e++) {
                PACc v = { (d == e) ? 1.0 : 0.0, 0.0 };
                ic0_1[(size_t) e * n + d] = v;
            }
        }
        pac_state_set_ic(map, x0);
        g_pac_traj_keep_st0 = 1;
        pac_traj_reset(neq);
        if (pac_reset(ckt, map, n, NULL, NULL, omega1, useforce1, ic0_1, 0, T)) {
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_traj_free();
            return PAC_ADJOINT_FAIL;
        }
        PSSpacActive = 1;
        g_pac_force_exact_jac = 1;
        err = pac_run_period(ckt, job, T, "TD-shoot NOISE pass A (Phi + trajectory)");
        g_pac_force_exact_jac = 0;
        PSSpacActive = 0;
        g_pac_traj_active = 0;
        g_pac_traj_keep_st0 = 0;
        pac_traj_report("PNOISE-TD");

        if (err || g_pac_failed || g_pac_traj_failed) {
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_free(); pac_traj_free();
            return PAC_ADJOINT_FAIL;
        }
        if (PSSsaltActive && PSSsaltEventsSeen() > 0) {
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_free(); pac_traj_free();
            return PAC_ADJOINT_REFUSED;
        }

        Phi = TMALLOC(PACc, (size_t) n * n);
        if (!Phi) {
            tfree(omega1); tfree(useforce1); tfree(ic0_1);
            pac_free(); pac_traj_free();
            return PAC_ADJOINT_FAIL;
        }
        for (d = 0; d < n; d++) {
            int row_d = pac_map_row(map, d);
            for (e = 0; e < n; e++)
                Phi[(size_t) d * n + e] = pac_get(g_pac_Cprev, row_d, n, e);
        }
        tfree(omega1); tfree(useforce1); tfree(ic0_1);
        pac_free();
    }

    if (getenv("PNOISE_DEBUG"))
        fprintf(stderr, "PNOISE-DEBUG: g_pac_traj_n=%d T=%.9g avg_h=%.9g\n",
                g_pac_traj_n, T, g_pac_traj_n > 0 ? T / g_pac_traj_n : 0.0);

    if (g_pac_traj_n < 1) {
        tfree(Phi); pac_traj_free();
        return PAC_ADJOINT_FAIL;
    }

    /* ---- pass B: G(w)[d][f] -- IDENTICAL design to
     * pac_engine_run_adjoint()'s own pass B above. */
    G = TMALLOC(PACc, (size_t) n * nfreq);
    if (!G) { tfree(Phi); pac_traj_free(); return PAC_ADJOINT_FAIL; }
    {
        int M2 = n * nfreq;
        double *omega2 = TMALLOC(double, M2);
        int *useforce2 = TMALLOC(int, M2);
        PACc *ic0_2 = TMALLOC(PACc, (size_t) n * M2);
        if (!omega2 || !useforce2 || !ic0_2) {
            if (omega2) tfree(omega2); if (useforce2) tfree(useforce2); if (ic0_2) tfree(ic0_2);
            tfree(Phi); tfree(G); pac_traj_free();
            return PAC_ADJOINT_FAIL;
        }
        memset(ic0_2, 0, (size_t) n * M2 * sizeof(PACc));
        for (d = 0; d < n; d++)
            for (f = 0; f < nfreq; f++) {
                int col = d * nfreq + f;
                omega2[col] = 2.0 * M_PI * freqs[f];
                useforce2[col] = 0;
                { PACc v = { 1.0, 0.0 }; ic0_2[(size_t) d * M2 + col] = v; }
            }
        pac_state_set_ic(map, x0);
        if (pac_reset(ckt, map, M2, NULL, NULL, omega2, useforce2, ic0_2, 1, T)) {
            tfree(omega2); tfree(useforce2); tfree(ic0_2);
            tfree(Phi); tfree(G); pac_traj_free();
            return PAC_ADJOINT_FAIL;
        }
        PSSpacActive = 1;
        g_pac_force_exact_jac = 1;
        err = pac_run_period(ckt, job, T, "TD-shoot NOISE pass B (G(w))");
        g_pac_force_exact_jac = 0;
        PSSpacActive = 0;
        if (err || g_pac_failed) {
            tfree(omega2); tfree(useforce2); tfree(ic0_2);
            pac_free(); tfree(Phi); tfree(G); pac_traj_free();
            return PAC_ADJOINT_FAIL;
        }
        pac_finish_demod_wraparound();
        {
            double invT = 1.0 / T;
            for (d = 0; d < n; d++)
                for (f = 0; f < nfreq; f++) {
                    PACc v = g_pac_demod[(size_t) (outrow - 1) * M2 + (d * nfreq + f)];
                    G[(size_t) d * nfreq + f] = pc_scale(v, invT);
                }
        }
        pac_free();
        tfree(omega2); tfree(useforce2); tfree(ic0_2);
    }

    /* ---- M(w) -- IDENTICAL design to pac_engine_run_adjoint()'s own
     * boundary solve above. */
    Mw = TMALLOC(PACc, (size_t) n * nfreq);
    if (!Mw) { tfree(Phi); tfree(G); pac_traj_free(); return PAC_ADJOINT_FAIL; }
    {
        PACc *AT = TMALLOC(PACc, (size_t) n * n);
        PACc *rhs = TMALLOC(PACc, n);
        PACc *sol = TMALLOC(PACc, n);
        if (!AT || !rhs || !sol) {
            if (AT) tfree(AT); if (rhs) tfree(rhs); if (sol) tfree(sol);
            tfree(Phi); tfree(G); tfree(Mw); pac_traj_free();
            return PAC_ADJOINT_FAIL;
        }
        for (f = 0; f < nfreq; f++) {
            double w = 2.0 * M_PI * freqs[f];
            double wT = w * T;
            PACc ephase = { cos(wT), sin(wT) };
            for (d = 0; d < n; d++)
                for (e = 0; e < n; e++) {
                    PACc v = pc_scale(Phi[(size_t) e * n + d], -1.0);
                    if (e == d) v = pc_add(v, ephase);
                    AT[(size_t) d * n + e] = v;
                }
            for (d = 0; d < n; d++) rhs[d] = G[(size_t) d * nfreq + f];
            if (pac_complex_solve(AT, rhs, n, sol)) {
                fprintf(stderr, "PNOISE-TD: WARNING -- transposed boundary solve "
                        "singular at f=%.6g Hz; result at this frequency skipped.\n",
                        freqs[f]);
                for (d = 0; d < n; d++) { PACc nan = { strtod("NaN", NULL), strtod("NaN", NULL) }; sol[d] = nan; }
            }
            for (d = 0; d < n; d++)
                Mw[(size_t) d * nfreq + f] = sol[d];
        }
        tfree(AT); tfree(rhs); tfree(sol);
    }
    tfree(Phi); tfree(G);

    /* ---- minimal NOISEAN context for the device noise routines (they
     * cast ckt->CKTcurJob to NOISEAN* and read NStpsSm) -- same idiom
     * dcpss.c's pnoise_sweep() uses. */
    memset(&nj, 0, sizeof(nj));
    nj.output = job->PnOutNode;
    nj.outputRef = job->PnOutNode;
    nj.NstartFreq = freqs[0];
    nj.NstopFreq = freqs[nfreq - 1];
    nj.NnumSteps = nfreq;
    nj.NstpType = job->PACstepType;
    nj.NStpsSm = (want_table || fold_on) ? 1 : 0;
    nj.JOBname = "pnoise-td";
    nd.prtSummary = 0;
    oldJob = ckt->CKTcurJob;
    ckt->CKTcurJob = (JOB *) &nj;
    have_job_saved = 1;

    /* ---- backward-sweep buffers (same sizes pac_engine_run_adjoint()
     * allocates for its own backward sweep). */
    Cdiag  = TMALLOC(double, n);
    abarS  = TMALLOC(PACc, (size_t) (map->generic ? neq + 1 : n) * nfreq);
    abarT  = TMALLOC(PACc, (size_t) (map->generic ? neq + 1 : n) * nfreq);
    tdu_re = TMALLOC(double, (size_t) (neq + 1));
    tdu_im = TMALLOC(double, (size_t) (neq + 1));
    scratch = TMALLOC(double, (size_t) (neq + 1));
    save_rhsold = TMALLOC(double, (size_t) (neq + 1));
    if (!Cdiag || !abarS || !abarT || !tdu_re || !tdu_im || !scratch || !save_rhsold) {
        err = PAC_ADJOINT_FAIL; goto done;
    }
    memset(abarS, 0, (size_t) (map->generic ? neq + 1 : n) * nfreq * sizeof(PACc));
    memset(abarT, 0, (size_t) (map->generic ? neq + 1 : n) * nfreq * sizeof(PACc));

    memcpy(save_rhsold, ckt->CKTrhsOld, (size_t) (neq + 1) * sizeof(double));
    save_time = ckt->CKTtime; save_delta = ckt->CKTdelta;
    save_ag0 = ckt->CKTag[0]; save_ag1 = ckt->CKTag[1];
    have_ckt_saved = 1;

    /* ---- naming pass: one summary slot per noise generator, only if a
     * per-device table was requested (same as dcpss.c's own
     * pnoise_sweep()/NStpsSm=0 default when it isn't). Bias: the trajectory's
     * own first recorded step -- a naming pass names generators by
     * instance/model (RESnoise/DIOnoise/MOS1noise/BSIM4noise's own N_OPEN
     * cases all confirmed bias-independent), so any valid in-circuit
     * operating point works; step 0 is the cheapest one already at hand. */
    if (want_table || fold_on) {
        double dummy = 0.0;
        if (pac_adj_refactor_at(ckt, &g_pac_traj[0], 1)) {
            fprintf(stderr, "PNOISE-TD: naming-pass refactor failed\n");
            err = PAC_ADJOINT_FAIL; goto done;
        }
        for (d = 0; d < DEVmaxnum; d++)
            if (DEVices[d] && DEVices[d]->DEVnoise && ckt->CKThead[d])
                DEVices[d]->DEVnoise(N_DENS, N_OPEN, ckt->CKThead[d], ckt, &nd, &dummy);
        ngen = nd.numPlots;
        if (ngen > 0) {
            if (want_table) {
                gen_dens = TMALLOC(double, (size_t) ngen * nfreq);
                memset(gen_dens, 0, (size_t) ngen * nfreq * sizeof(double));
            }
            nd.outpVector = TMALLOC(double, ngen);
            nd.prtSummary = 1;
            names_local = nd.namelist;
            if (want_table) {
                gen_names = nd.namelist;   /* transfer ownership from nd */
                nd.namelist = NULL;
            }
        }
    }

    /* card #556: folded flicker.  Flicker slots = names with "1overf" (not the sidewall one); each owns one N_GAIN call of
     * the device noise routines, in the same instance order. */
    if (fold_on && ngen > 0) {
        int nfl = 0, s2, c;
        int Ns = g_pac_traj_n;
        for (s2 = 0; s2 < ngen; s2++) {
            const char *nm = (const char *) names_local[s2];
            if (nm && strstr(nm, "1overf") && !strstr(nm, "1overfsw")) nfl++;
        }
        if (nfl > 0 && (double) nfl * nfreq * Ns * 16.0 > 3.0e9) {
            fprintf(stderr, "PNOISE-TD: WARNING folded-flicker table would need %.1f GB -- flicker evaluated at the output "
                    "frequency only (the pre-card-#556 treatment)\n", (double) nfl * nfreq * Ns * 16.0 / 1e9);
            nfl = 0;
        }
        if (nfl > 0) {
            fold.nfl = nfl; fold.nfreq = nfreq; fold.N = Ns;
            fold.slot = TMALLOC(int, nfl); fold.tot = TMALLOC(int, nfl);
            fold.ef = TMALLOC(double, nfl); fold.ef_set = TMALLOC(int, nfl);
            fold.zr = TMALLOC(double, (size_t) nfl * nfreq * Ns); fold.zi = TMALLOC(double, (size_t) nfl * nfreq * Ns);
            fold.E = TMALLOC(double, (size_t) nfl * nfreq);
            fold.tc = TMALLOC(double, Ns); fold.wk = TMALLOC(double, Ns);
            fold.capr = TMALLOC(double, nfl); fold.capi = TMALLOC(double, nfl);
            if (!fold.slot || !fold.tot || !fold.ef || !fold.ef_set || !fold.zr || !fold.zi || !fold.E || !fold.tc ||
                !fold.wk || !fold.capr || !fold.capi) {
                err = PAC_ADJOINT_FAIL; goto done;
            }
            memset(fold.zr, 0, (size_t) nfl * nfreq * Ns * sizeof(double));
            memset(fold.zi, 0, (size_t) nfl * nfreq * Ns * sizeof(double));
            memset(fold.E, 0, (size_t) nfl * nfreq * sizeof(double));
            memset(fold.ef_set, 0, nfl * sizeof(int));
            memset(fold.ef, 0, nfl * sizeof(double));
            c = 0;
            for (s2 = 0; s2 < ngen; s2++) {
                const char *nm = (const char *) names_local[s2];
                if (nm && strstr(nm, "1overf") && !strstr(nm, "1overfsw")) {
                    int t2, tot = -1;
                    size_t ln = strlen(nm);
                    fold.slot[c] = s2;
                    for (t2 = s2 + 1; t2 < ngen; t2++) {
                        const char *tn = (const char *) names_local[t2];
                        size_t lt = tn ? strlen(tn) : 0;
                        if (tn && lt < ln && strncmp(tn, nm, lt) == 0) { tot = t2; break; }
                    }
                    fold.tot[c] = tot;
                    c++;
                }
            }
        }
    }

    /* ---- backward adjoint sweep -- Phi/G/M(w) reuse aside, the per-step
     * recursion below (Cdiag/abarS/abarT/beta/gamma, the SMPtSolve() pair,
     * the M(w)/output-tap scatter into tdu_re/tdu_im) is copied VERBATIM
     * from pac_engine_run_adjoint()'s own backward sweep above -- see that
     * function's header for the derivation of every term; this function
     * does not re-derive them, it only reuses nu_k = (tdu_re,tdu_im) (the
     * per-(k,f) SMPtSolve() result, BEFORE any e^{+jwt_k} phase multiply)
     * as this card's own H(t_k,f) -- see THIS function's header above for
     * why that (not Lambda(w)) is the quantity card #509 needs. */
    for (k = g_pac_traj_n - 1; k >= 0; k--) {
        const PACtrajStep *step = &g_pac_traj[k];
        int is_last = (k == g_pac_traj_n - 1);
        double h_next = is_last ? g_pac_traj[0].h : g_pac_traj[k + 1].h;
        double ag1eff = (step->order == 2) ? step->ag1 : 0.0;

        if (map->generic ? pac_adj_refactor_generic(ckt, step, 1) : pac_adj_refactor_at(ckt, step, 1)) {
            fprintf(stderr, "PNOISE-TD: refactor failed at t=%.9g (backward step %d)\n",
                    step->t, k);
            err = PAC_ADJOINT_FAIL; goto done;
        }

        for (e = 0; !map->generic && e < n; e++) {
            const PACstateRef *ref = &map->refs[e];
            Cdiag[e] = ref->isCap ? step->ag0 * ref->dev.cap->CAPcapac
                                   : step->ag0 * ref->dev.ind->INDinduct;
        }

        for (f = 0; f < nfreq; f++) {
            double w = 2.0 * M_PI * freqs[f];
            PACc ph2 = pc_conjexp(w * step->t);
            PACc coef = pc_scale(ph2, 0.5 * (step->h + h_next) / T);

            memset(tdu_re, 0, (size_t) (neq + 1) * sizeof(double));
            memset(tdu_im, 0, (size_t) (neq + 1) * sizeof(double));
            tdu_re[outrow] += coef.re; tdu_im[outrow] += coef.im;

            if (is_last) {
                for (e = 0; e < n; e++) {
                    PACc Mv = Mw[(size_t) e * nfreq + f];
                    int row2 = pac_map_row(map, e);
                    if (row2 > 0) { tdu_re[row2] += Mv.re; tdu_im[row2] += Mv.im; }
                }
            }

            if (map->generic) {
                /* card #542: z += C_k^T (alpha + ag0 beta) over every row (alpha/beta: adjoints of Q_k, T_k
                 * from the later steps, see the generic recursion in the design note) */
                int ti;
                for (ti = 0; ti < g_pac_C.n; ti++) {
                    int r = g_pac_C.ri[ti], c = g_pac_C.ci[ti];
                    PACc al, be;
                    double v = g_pac_C.v[ti];
                    if (r <= 0 || c <= 0) continue;
                    al = abarS[(size_t) r * nfreq + f];
                    be = abarT[(size_t) r * nfreq + f];
                    tdu_re[c] += v * (al.re + step->ag0 * be.re);
                    tdu_im[c] += v * (al.im + step->ag0 * be.im);
                }
            }
            for (e = 0; !map->generic && e < n; e++) {
                const PACstateRef *ref = &map->refs[e];
                PACc abT = abarT[(size_t) e * nfreq + f];
                PACc abS_ext = abarS[(size_t) e * nfreq + f];
                PACc abS_tot = pc_add(abS_ext, pc_scale(abT, Cdiag[e]));
                if (ref->isCap) {
                    int p = ref->dev.cap->CAPposNode, ngc = ref->dev.cap->CAPnegNode;
                    if (p > 0)   { tdu_re[p] += abS_tot.re; tdu_im[p] += abS_tot.im; }
                    if (ngc > 0) { tdu_re[ngc] -= abS_tot.re; tdu_im[ngc] -= abS_tot.im; }
                } else {
                    int ibr = ref->dev.ind->INDbrEq;
                    tdu_re[ibr] += abS_tot.re; tdu_im[ibr] += abS_tot.im;
                }
            }

            SMPtSolve(ckt->CKTmatrix, tdu_re, scratch);
            SMPtSolve(ckt->CKTmatrix, tdu_im, scratch);
            /* tdu_re/tdu_im now hold nu_k(f) = dY0(w)/d(beta_k) -- card
             * #509's own per-step DETERMINISTIC transfer, exactly as
             * pac_engine_run_adjoint() computes for Lambda(w) (before that
             * function's own e^{+jwt_k} phase multiply, which a stationary
             * source's per-step forcing needs and an isolated noise sample
             * does not -- see this function's header).
             *
             * White noise, unlike a deterministic AC source, cannot be
             * "sampled" as beta_k = i_n(t_k) directly (continuous white
             * noise has a Dirac-delta autocorrelation -- no well-defined
             * pointwise sample). The physically consistent discretisation
             * (same one ngspice's OWN trnoise white generator and the
             * TK96/Demir-Roychowdhury discrete-noise literature use):
             * represent i_n(t) over step k as a ZERO-ORDER-HOLD value n_k
             * (entering the per-step recursion exactly like a deterministic
             * beta_k -- "same stamp as a small-signal source at step k",
             * this card's own brief), with Var(n_k) = S_k / w_k so that the
             * CHARGE n_k*w_k over the step reproduces the true continuous
             * process's own charge variance S_k*w_k (a white process
             * integrated over an interval of length w_k has variance
             * S_k*w_k) -- w_k = 0.5*(h_k+h_next) is the SAME periodic-
             * trapezoidal weight "coef" above already uses for the OUTPUT
             * side, used here symmetrically for the INJECTION side. Since
             * the n_k are independent across k (white/uncorrelated), the
             * output sideband-0 response Y0_random(f) = Sum_k nu_k(f)*n_k
             * has variance Sum_k |nu_k(f)|^2 * Var(n_k) = Sum_k
             * S_k*|nu_k(f)|^2 / w_k -- this is dens_k/w_k below (dens_k =
             * S_k*|nu_k|^2, via the ORDINARY NevalSrc "gain*S" formula with
             * gain from nu_k -- w_k has NOT been folded into nu_k itself,
             * see the SMPtSolve() calls above). The remaining T is this
             * card's own empirically-confirmed normalisation, found by the
             * LTI-vs-plain-.noise check below (an EARLIER version of this
             * loop tried BOTH "no extra weight at all" and "multiply by
             * w_k/T (mirroring the output tap's own coef)" -- both wrong,
             * off by ~N and ~N^2 respectively, N = accepted-step count;
             * "T * dens_k/w_k", below, matches plain .noise to 6e-6
             * relative on pxf_lti.cir at Points=2000, see docs/pss-pac.md's
             * `.pnoise` section for the full derivation and the numbers all
             * three candidates gave). */

            {
                double dens = 0.0;
                double wk = 0.5 * (step->h + h_next);
                for (row = 0; row <= neq; row++) {
                    ckt->CKTrhs[row]  = tdu_re[row];
                    ckt->CKTirhs[row] = tdu_im[row];
                }
                nd.freq = freqs[f];
                nd.delFreq = 0.0;
                if (ngen > 0) {
                    nd.outNumber = 0;
                    if (fold.nfl > 0 && !fold.nonpow) {
                        g_nev_gain_cap = 1; g_nev_gain_n = 0; g_nev_gain_max = fold.nfl;
                        g_nev_gain_re = fold.capr; g_nev_gain_im = fold.capi;
                    }
                    for (d = 0; d < DEVmaxnum; d++)
                        if (DEVices[d] && DEVices[d]->DEVnoise && ckt->CKThead[d])
                            DEVices[d]->DEVnoise(N_DENS, N_CALC, ckt->CKThead[d], ckt, &nd, &dens);
                    g_nev_gain_cap = 0;
                    if (gen_dens)
                        for (d = 0; d < ngen; d++)
                            gen_dens[(size_t) d * nfreq + f] += nd.outpVector[d] * T / wk;
                    if (fold.nfl > 0 && !fold.nonpow) {
                        if (g_nev_gain_n != fold.nfl) {
                            fprintf(stderr, "PNOISE-TD: WARNING %d flicker gain evaluations for %d flicker slots -- flicker "
                                    "evaluated at the output frequency only (pre-card-#556 treatment)\n",
                                    g_nev_gain_n, fold.nfl);
                            fold.nonpow = 1;
                        } else {
                            double dummy2 = 0.0;
                            int cc;
                            double *F1v = TMALLOC(double, fold.nfl);
                            int need2 = 0;
                            for (cc = 0; cc < fold.nfl; cc++) {
                                F1v[cc] = nd.outpVector[fold.slot[cc]];
                                if (F1v[cc] > 0.0) need2 = 1;
                            }
                            if (need2) {
                                /* density at 2f: the exponent of a power law (and a test that it is one) */
                                nd.freq = 2.0 * freqs[f];
                                nd.outNumber = 0;
                                for (d = 0; d < DEVmaxnum; d++)
                                    if (DEVices[d] && DEVices[d]->DEVnoise && ckt->CKThead[d])
                                        DEVices[d]->DEVnoise(N_DENS, N_CALC, ckt->CKThead[d], ckt, &nd, &dummy2);
                                nd.freq = freqs[f];
                            }
                            for (cc = 0; cc < fold.nfl; cc++) {
                                double u_re = fold.capr[cc], u_im = fold.capi[cc];
                                double u2 = u_re * u_re + u_im * u_im, Fv = F1v[cc];
                                size_t idx = ((size_t) cc * nfreq + f) * fold.N + k;
                                if (Fv > 0.0 && u2 > 0.0) {
                                    double F2 = nd.outpVector[fold.slot[cc]];
                                    double sc = sqrt(Fv / u2);
                                    double ph_re = cos(w * step->t), ph_im = sin(w * step->t);
                                    double zre = u_re * sc, zim = u_im * sc;
                                    fold.zr[idx] = zre * ph_re - zim * ph_im;
                                    fold.zi[idx] = zre * ph_im + zim * ph_re;
                                    fold.E[(size_t) cc * nfreq + f] += Fv * T / wk;
                                    if (F2 > 0.0) {
                                        double e = log(Fv / F2) / log(2.0);
                                        if (!fold.ef_set[cc]) { fold.ef[cc] = e; fold.ef_set[cc] = 1; }
                                        else if (fabs(e - fold.ef[cc]) > 1e-6 * (1.0 + fabs(e)))
                                            fold.nonpow = 1;
                                    } else
                                        fold.nonpow = 1;
                                }
                            }
                            tfree(F1v);
                            fold.tc[k] = step->t + 0.25 * (h_next - step->h);
                            fold.wk[k] = wk;
                        }
                    }
                } else {
                    for (d = 0; d < DEVmaxnum; d++)
                        if (DEVices[d] && DEVices[d]->DEVnoise && ckt->CKThead[d])
                            DEVices[d]->DEVnoise(N_DENS, N_CALC, ckt->CKThead[d], ckt, &nd, &dens);
                }
                onoise_out[f] += dens * T / wk;
                if (f == 0 && getenv("PNOISE_DEBUG_STEPS"))
                    fprintf(stderr, "PNOISE-STEP: k=%d t=%.9g h=%.9g hnext=%.9g dens=%.17g nu0_re=%.9g nu0_im=%.9g running=%.17g\n",
                            k, step->t, step->h, h_next, dens, tdu_re[outrow], tdu_im[outrow], onoise_out[f]);
            }

            if (map->generic) {
                /* gamma = nu - beta_old; alpha_{k-1} = ag0 gamma; beta_{k-1} = a1 gamma */
                for (row = 1; row <= neq; row++) {
                    PACc be = abarT[(size_t) row * nfreq + f];
                    PACc gm = { tdu_re[row] - be.re, tdu_im[row] - be.im };
                    abarS[(size_t) row * nfreq + f] = pc_scale(gm, step->ag0);
                    abarT[(size_t) row * nfreq + f] = pc_scale(gm, ag1eff);
                }
            }
            for (e = 0; !map->generic && e < n; e++) {
                const PACstateRef *ref = &map->refs[e];
                PACc beta, gamma;
                if (ref->isCap) {
                    int p = ref->dev.cap->CAPposNode, ngc = ref->dev.cap->CAPnegNode;
                    PACc nup = (p > 0)   ? (PACc){ tdu_re[p], tdu_im[p] }     : (PACc){ 0.0, 0.0 };
                    PACc nun = (ngc > 0) ? (PACc){ tdu_re[ngc], tdu_im[ngc] } : (PACc){ 0.0, 0.0 };
                    PACc diffnu = pc_sub(nup, nun);
                    beta  = pc_scale(diffnu, step->ag0 * ref->dev.cap->CAPcapac);
                    gamma = pc_scale(diffnu, ag1eff);
                } else {
                    int ibr = ref->dev.ind->INDbrEq;
                    PACc nub = { tdu_re[ibr], tdu_im[ibr] };
                    beta  = pc_scale(nub, -step->ag0 * ref->dev.ind->INDinduct);
                    gamma = pc_scale(nub, -ag1eff);
                }
                {
                    PACc abT_old = abarT[(size_t) e * nfreq + f];
                    PACc new_abS = pc_sub(beta, pc_scale(abT_old, Cdiag[e]));
                    PACc new_abT = pc_add(gamma, pc_scale(abT_old, -ag1eff));
                    abarS[(size_t) e * nfreq + f] = new_abS;
                    abarT[(size_t) e * nfreq + f] = new_abT;
                }
            }
        }
    }

    memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
    ckt->CKTtime = save_time; ckt->CKTdelta = save_delta;
    ckt->CKTag[0] = save_ag0; ckt->CKTag[1] = save_ag1;
    have_ckt_saved = 0;

    if (fold.nfl > 0) {
        if (fold.nonpow) {
            fprintf(stderr, "PNOISE-TD: WARNING a flicker density is not a power law in frequency (BSIM4 fnoimod=1?) -- flicker "
                    "evaluated at the output frequency only (the pre-card-#556 treatment, exact for a bias that does not "
                    "vary fast)\n");
        } else {
            pac_fold_apply(&fold, freqs, T, onoise_out, gen_dens);
        }
    }

    err = PAC_ADJOINT_OK;
    *ngen_out = ngen;
    *gen_dens_out = gen_dens;  gen_dens = NULL;
    *gen_names_out = gen_names; gen_names = NULL;

 done:
    if (have_ckt_saved) {
        memcpy(ckt->CKTrhsOld, save_rhsold, (size_t) (neq + 1) * sizeof(double));
        ckt->CKTtime = save_time; ckt->CKTdelta = save_delta;
        ckt->CKTag[0] = save_ag0; ckt->CKTag[1] = save_ag1;
    }
    if (have_job_saved) ckt->CKTcurJob = oldJob;
    g_nev_gain_cap = 0;
    pac_fold_free(&fold);
    if (Mw) tfree(Mw);
    pac_traj_free();
    if (Cdiag) tfree(Cdiag);
    if (abarS) tfree(abarS);
    if (abarT) tfree(abarT);
    if (tdu_re) tfree(tdu_re);
    if (tdu_im) tfree(tdu_im);
    if (scratch) tfree(scratch);
    if (save_rhsold) tfree(save_rhsold);
    if (nd.outpVector) tfree(nd.outpVector);
    if (nd.namelist) tfree(nd.namelist);
    if (gen_dens) tfree(gen_dens);
    if (gen_names) tfree(gen_names);
    return err;
}

/* card #509: `.pnoise ... method=shoot solver=td` entry point, called from
 * pac_pstb_run_ex() below -- always the ordinary driven case (trig_idx
 * resolved from trig_node_num, exactly like pac_td_pxf()); free-running is
 * refused by pac_engine_run_adjoint_noise() itself if ever reached (see
 * that function's header -- "SCOPE"). */
static int
pac_td_pnoise(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T,
              const double *fT, int trig_node_num)
{
    PACstateMap map;
    int nfreq, outrow, err, trig_idx, want_table = job->PSSpnDevTable;
    double *freqs, *onoise;
    double *gen_dens = NULL; IFuid *gen_names = NULL; int ngen = 0;
    IFuid freqUid, *nameList = NULL;
    runDesc *plot = NULL;
    double t0, t1;
    int pn_sqr;

    NG_IGNORE(fT);

    if (!job->PnOutNode) {
        fprintf(stderr, "PNOISE-TD: no output node given -- skipped\n");
        return E_PARMVAL;
    }
    outrow = job->PnOutNode->number;
    if (outrow <= 0) {
        fprintf(stderr, "PNOISE-TD: output node '%s' is not in the matrix "
                "(literal ground?) -- skipped\n", job->PnOutNode->name);
        return E_PARMVAL;
    }
    if (job->PACfStart <= 0.0 || job->PACfStop < job->PACfStart || job->PACpoints < 1) {
        fprintf(stderr, "PNOISE-TD: bad frequency sweep -- skipped\n");
        return E_PARMVAL;
    }
    if (pac_noise_check_devices(ckt))
        return E_UNSUPP;

    map.n = n; map.refs = NULL; map.generic = 0; map.ns = 0; map.srows = NULL;
    if (pac_build_map(ckt, &map) || map.n != n) {
        fprintf(stderr, "PNOISE-TD: internal state-map mismatch -- skipped\n");
        pac_map_free(&map);
        return E_INTERN;
    }
    {
        int perr = pac_generic_prepare(ckt, job, &map, x0, T, trig_node_num, "PNOISE-TD");
        if (perr) { pac_map_free(&map); return perr; }
    }
    trig_idx = pac_find_trig_idx(&map, trig_node_num);
    if (trig_node_num >= 0 && trig_idx < 0) {
        fprintf(stderr, "PNOISE-TD: REFUSED -- the trigger= node (free-running orbit) is not one of the circuit's state "
                "rows (a capacitor / inductor state or a column of the charge matrix); pick a state node as trigger (card #557) -- "
                "see docs/pss-pac.md.\n");
        pac_map_free(&map);
        return E_UNSUPP;
    }

    nfreq = pac_sweep_count(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop);
    if (nfreq < 1) { pac_map_free(&map); return E_PARMVAL; }
    freqs = TMALLOC(double, nfreq);
    pac_sweep_fill(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop, freqs, nfreq);
    onoise = TMALLOC(double, nfreq);
    memset(onoise, 0, (size_t) nfreq * sizeof(double));

    fprintf(stderr, "PNOISE-TD: time-domain shooting-based periodic noise, %d "
                    "frequency point%s from %.6g to %.6g Hz, n=%d reduced state%s, "
                    "output node '%s' (row %d)%s\n",
            nfreq, nfreq == 1 ? "" : "s", job->PACfStart, job->PACfStop, n,
            n == 1 ? "" : "s", job->PnOutNode->name, outrow,
            want_table ? ", per-device table" : "");

    t0 = SPfrontEnd->IFseconds();
    err = pac_engine_run_adjoint_noise(ckt, job, &map, x0, T, freqs, nfreq, outrow,
                                        want_table, trig_idx, onoise,
                                        &gen_dens, &gen_names, &ngen);
    g_pac_traj_keep_st0 = 0;
    t1 = SPfrontEnd->IFseconds();

    if (err == PAC_ADJOINT_REFUSED) {
        fprintf(stderr, "PNOISE-TD: REFUSED -- saltation (ideal-switch) events seen "
                "on this orbit (or the opt-in PAC_FREERUN_PROJ projection on a free-running one); this engine needs the smooth-switch-only adjoint "
                "(card #504's own scope cut) -- no fallback engine exists for "
                "pnoise (unlike .pxf's adjoint=0). See docs/pss-pac.md.\n");
        pac_map_free(&map); tfree(freqs); tfree(onoise);
        if (gen_dens) tfree(gen_dens);
        if (gen_names) tfree(gen_names);
        return E_UNSUPP;
    }
    if (err != PAC_ADJOINT_OK) {
        fprintf(stderr, "PNOISE-TD: engine failed\n");
        pac_map_free(&map); tfree(freqs); tfree(onoise);
        if (gen_dens) tfree(gen_dens);
        if (gen_names) tfree(gen_names);
        return E_INTERN;
    }
    fprintf(stderr, "PNOISE-TD: done, wall_time=%.3fs\n", t1 - t0);

    /* card #193 convention (same as .noise / the frequency-domain
     * .pnoise): default reports sqrt(density) (V or A per sqrt(Hz)); `set
     * sqrnoise` reports the squared density (V^2/Hz or A^2/Hz) instead. */
    pn_sqr = cp_getvar("sqrnoise", CP_BOOL, NULL, 0);

    {
        int nvec = 1 + (want_table ? ngen : 0);
        int g, e3;

        nameList = TMALLOC(IFuid, nvec);
        SPfrontEnd->IFnewUid(ckt, &nameList[0], NULL, "onoise_spectrum", UID_OTHER, NULL);
        for (g = 0; g < (want_table ? ngen : 0); g++)
            nameList[1 + g] = gen_names[g];

        SPfrontEnd->IFnewUid(ckt, &freqUid, NULL, "frequency", UID_OTHER, NULL);
        e3 = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob,
                        "PNoise Analysis (time domain, sideband 0)",
                        freqUid, IF_REAL, nvec, nameList, IF_REAL, &plot);
        if (!e3) {
            int i;
            if (job->PACstepType != 0)
                SPfrontEnd->OUTattributes(plot, NULL, OUT_SCALE_LOG, NULL);
            for (i = 0; i < nfreq; i++) {
                IFvalue freqData, valueData;
                double *rowv = TMALLOC(double, nvec);
                int j;
                double v0 = onoise[i] > 0.0 ? onoise[i] : 0.0;
                rowv[0] = pn_sqr ? v0 : sqrt(v0);
                for (j = 0; j < (want_table ? ngen : 0); j++) {
                    double v = gen_dens[(size_t) j * nfreq + i];
                    if (v < 0.0) v = 0.0;
                    rowv[1 + j] = pn_sqr ? v : sqrt(v);
                }
                freqData.rValue = freqs[i];
                valueData.v.numValue = nvec;
                valueData.v.vec.rVec = rowv;
                SPfrontEnd->OUTpData(plot, &freqData, &valueData);
                FREE(rowv);
            }
            SPfrontEnd->OUTendPlot(plot);
        } else {
            err = 1;
        }
        tfree(nameList);
    }

    pac_map_free(&map); tfree(freqs); tfree(onoise);
    if (gen_dens) tfree(gen_dens);
    if (gen_names) tfree(gen_names);
    return err ? E_INTERN : OK;
}

/* ------------------------------------------------------------------ */
/* `.pac ... method=shoot solver=td`: sideband-0 response at every node,
 * forcing = the netlist AC-source stimulus B0. fT/trig_node_num: card #489
 * second pass, see pacshoot.h's pac_pstb_run_ex() comment -- NULL/-1 for
 * the ordinary driven case. */
static int
pac_td_pac(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T,
           const double *fT, int trig_node_num)
{
    PACstateMap map;
    int neq = SMPmatSize(ckt->CKTmatrix);
    double *br, *bi;
    int nfreq, *outrows, i, numNames, err, trig_idx;
    double *freqs;
    PACc *Y0;
    IFuid freqUid, *nameList = NULL;
    runDesc *pacPlot = NULL;
    double t0, t1;

    if (job->PACfStart <= 0.0 || job->PACfStop < job->PACfStart || job->PACpoints < 1) {
        fprintf(stderr, "PAC-TD: bad frequency sweep -- skipped\n");
        return E_PARMVAL;
    }

    map.n = n; map.refs = NULL; map.generic = 0; map.ns = 0; map.srows = NULL;
    if (pac_build_map(ckt, &map) || map.n != n) {
        fprintf(stderr, "PAC-TD: internal state-map mismatch -- skipped\n");
        pac_map_free(&map);
        return E_INTERN;
    }

    {
        int perr = pac_generic_prepare(ckt, job, &map, x0, T, trig_node_num, "PAC-TD");
        if (perr) { pac_map_free(&map); return perr; }
    }
    trig_idx = pac_find_trig_idx(&map, trig_node_num);
    if (trig_node_num >= 0 && trig_idx < 0) {
        fprintf(stderr, "PAC-TD: REFUSED -- the trigger= node (free-running "
                "orbit, card #489) is not a native capacitor's own positive "
                "terminal (or, on the generic model-charge path, a charge-carrying node row) among this circuit's states, so the "
                "period-variation terminal projection this engine needs "
                "(docs/pss-pac.md's free-running section) cannot be formed. "
                "Not running the fixed-period engine uncorrected on a "
                "moving-period orbit.\n");
        pac_map_free(&map);
        return E_UNSUPP;
    }

    br = TMALLOC(double, neq + 1);
    bi = TMALLOC(double, neq + 1);
    if (!br || !bi) { if (br) tfree(br); if (bi) tfree(bi); pac_map_free(&map); return E_NOMEM; }
    memset(br, 0, (size_t) (neq + 1) * sizeof(double));
    memset(bi, 0, (size_t) (neq + 1) * sizeof(double));
    if (pac_capture_b0(ckt, br, bi, neq)) {
        fprintf(stderr, "PAC-TD: failed to capture the netlist AC-source stimulus\n");
        tfree(br); tfree(bi); pac_map_free(&map);
        return E_INTERN;
    }

    nfreq = pac_sweep_count(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop);
    if (nfreq < 1) { tfree(br); tfree(bi); pac_map_free(&map); return E_PARMVAL; }
    freqs = TMALLOC(double, nfreq);
    pac_sweep_fill(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop, freqs, nfreq);

    if (CKTnames(ckt, &numNames, &nameList) || numNames != neq) {
        fprintf(stderr, "PAC-TD: CKTnames() mismatch -- skipped\n");
        tfree(br); tfree(bi); pac_map_free(&map); tfree(freqs);
        return E_INTERN;
    }
    outrows = TMALLOC(int, numNames);
    for (i = 0; i < numNames; i++) outrows[i] = i + 1;
    Y0 = TMALLOC(PACc, (size_t) numNames * nfreq);

    fprintf(stderr, "PAC-TD: time-domain shooting-based PAC, %d frequency point%s from "
                    "%.6g to %.6g Hz, n=%d reduced state%s, %d node%s; forcing: %s\n",
            nfreq, nfreq == 1 ? "" : "s", job->PACfStart, job->PACfStop, n, n == 1 ? "" : "s",
            numNames, numNames == 1 ? "" : "s", "netlist AC-flagged source(s)");

    t0 = SPfrontEnd->IFseconds();
    err = pac_engine_run(ckt, job, &map, x0, T, br, bi, freqs, nfreq, outrows, numNames, NULL, NULL, Y0,
                          fT, trig_idx);
    t1 = SPfrontEnd->IFseconds();
    fprintf(stderr, "PAC-TD: %s, wall_time=%.3fs (%.4fs/point)\n",
            err ? "FAILED" : "done", t1 - t0, nfreq > 0 ? (t1 - t0) / nfreq : 0.0);

    if (!err) {
        SPfrontEnd->IFnewUid(ckt, &freqUid, NULL, "frequency", UID_OTHER, NULL);
        err = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob, "PAC Analysis (time domain, sideband 0)",
                                        freqUid, IF_REAL, numNames, nameList, IF_COMPLEX, &pacPlot);
        if (!err) {
            if (job->PACstepType != 0)
                SPfrontEnd->OUTattributes(pacPlot, NULL, OUT_SCALE_LOG, NULL);
            for (i = 0; i < nfreq; i++) {
                IFvalue freqData, valueData;
                IFcomplex *data = TMALLOC(IFcomplex, numNames);
                int j;
                for (j = 0; j < numNames; j++) {
                    data[j].real = Y0[(size_t) j * nfreq + i].re;
                    data[j].imag = Y0[(size_t) j * nfreq + i].im;
                }
                freqData.rValue = freqs[i];
                valueData.v.numValue = numNames;
                valueData.v.vec.cVec = data;
                SPfrontEnd->OUTpData(pacPlot, &freqData, &valueData);
                FREE(data);
            }
            SPfrontEnd->OUTendPlot(pacPlot);
        }
    }

    tfree(br); tfree(bi); pac_map_free(&map); tfree(freqs); tfree(outrows); tfree(Y0); tfree(nameList);
    return err ? E_INTERN : OK;
}

/* ------------------------------------------------------------------ */
/* `.pstb probe=<Vname>`: Tian two-injection sideband-0 loop gain. Probe
 * resolution and RHS conventions copied verbatim from stban.c's native
 * loop-mode `.stb` (see that file for the derivation/validation this
 * mirrors) -- duplicated rather than shared because stban.c is another
 * worker's file this card must not edit (see pacshoot.c's header). */

typedef struct { double freq, value, aux; } PACcrossing;

static PACcrossing *
pac_push_crossing(PACcrossing *arr, int *n, int *cap, double freq, double value, double aux)
{
    if (*n >= *cap) {
        *cap = *cap ? *cap * 2 : 8;
        arr = TREALLOC(PACcrossing, arr, *cap);
    }
    arr[*n].freq = freq; arr[*n].value = value; arr[*n].aux = aux;
    (*n)++;
    return arr;
}

static double
pac_interp_logf(double f0, double f1, double y0, double y1, double target)
{
    double t;
    if (y1 == y0) return f0;
    t = (target - y0) / (y1 - y0);
    return f0 * pow(f1 / f0, t);
}

static int
pac_td_pstb(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T,
            const double *fT, int trig_node_num)
{
    PACstateMap map;
    int neq = SMPmatSize(ckt->CKTmatrix);
    GENinstance *inst;
    int branch_probe, node_a, node_b;
    int inj_node = 0;  /* node of the unit current injection / voltage read-back */
    int lgmode = 0;    /* probe is an lgprobe subcircuit instance */
    double *br1, *bi1, *br2, *bi2;
    int nfreq, i, err, trig_idx;
    double *freqs;
    int outrows[2];   /* inj_node, branch_probe */
    PACc *Phi, *Y1, *Y2;
    PACc *Tarr, *LGarr;
    double t0, t1;

    IFuid probeUid = job->PSTBprobeName;
    if (!probeUid) {
        fprintf(stderr, "PSTB: no probe= given -- skipped\n");
        return E_PARMVAL;
    }
    inst = CKTfndDev(ckt, probeUid);
    if (inst && inst->GENmodPtr->GENmodType == CKTtypelook("Vsource")) {
        /* plain probe=<Vname> (also what probe=<inst>:<terminal> reduces to: the 0V source
         * the parser spliced in series with that terminal) */
        branch_probe = CKTfndBranch(ckt, probeUid);
        if (branch_probe == 0) {
            fprintf(stderr, "PSTB: probe '%s' has no branch equation\n", (char *) probeUid);
            return E_NOTFOUND;
        }
        node_a = GENnode(inst)[0];
        node_b = GENnode(inst)[1];
        inj_node = node_b;
    } else if (!inst) {
        /* probe=<Xname>: an instance of the lgprobe subcircuit (stb/lgprobe.lib:
         *   Ii 0 x DC 0 AC 0 ; Vi x a DC 0 AC 1 ; Vb b x DC 0 ),
         * flattened to v.<Xname>.vi, v.<Xname>.vb, i.<Xname>.ii. */
        GENinstance *vi, *vb, *ii;
        char nm[300];
        const char *lc = (const char *) probeUid;
        size_t len = strlen(lc), q;
        if (len > 200) len = 200;
        #define PSTB_LGNAME(pre, post, dst) do { \
            snprintf(nm, sizeof nm, "%s.", pre); \
            for (q = 0; q < len; q++) nm[2 + q] = (char) tolower((unsigned char) lc[q]); \
            snprintf(nm + 2 + len, sizeof nm - 2 - len, ".%s", post); \
            dst = CKTfndDev(ckt, (IFuid) nm); } while (0)
        PSTB_LGNAME("v", "vi", vi); PSTB_LGNAME("v", "vb", vb); PSTB_LGNAME("i", "ii", ii);
        #undef PSTB_LGNAME
        if (!vi || !vb || !ii ||
            vi->GENmodPtr->GENmodType != CKTtypelook("Vsource") ||
            vb->GENmodPtr->GENmodType != CKTtypelook("Vsource") ||
            ii->GENmodPtr->GENmodType != CKTtypelook("Isource")) {
            fprintf(stderr, "PSTB: probe '%s' is neither an existing 0V V-source, nor an instance of the "
                    "lgprobe subcircuit (stb/lgprobe.lib: elements v.%s.vi, v.%s.vb, i.%s.ii), nor (parsed "
                    "earlier) an <inst>:<terminal>. Device-based (gain-nulling) probes are not supported "
                    "by this engine\n", (char *) probeUid, (char *) probeUid, (char *) probeUid,
                    (char *) probeUid);
            return E_NOTFOUND;
        }
        /* topology check: Vi x a, Vb b x, Ii 0 x -- one shared node x, current source from ground */
        if (GENnode(vi)[0] == 0 || GENnode(vi)[0] != GENnode(vb)[1] ||
            GENnode(ii)[1] != GENnode(vi)[0] || GENnode(ii)[0] != 0) {
            fprintf(stderr, "PSTB: probe '%s': v.%s.vi/vb and i.%s.ii do not have the lgprobe topology "
                    "(Ii 0 x; Vi x a; Vb b x) -- refused\n", (char *) probeUid, (char *) probeUid, (char *) probeUid);
            return E_UNSUPP;
        }
        branch_probe = CKTfndBranch(ckt, vi->GENname);
        if (branch_probe == 0) {
            fprintf(stderr, "PSTB: probe '%s': vi has no branch equation\n", (char *) probeUid);
            return E_NOTFOUND;
        }
        node_a = GENnode(vi)[1];
        node_b = GENnode(vb)[0];
        inj_node = GENnode(vi)[0];
        lgmode = 1;
    } else {
        fprintf(stderr, "PSTB: probe '%s' is not a V-source (supported: an existing 0V V-source, "
                "<inst>:<terminal>, or an instance of the lgprobe subcircuit)\n", (char *) probeUid);
        return E_NOTFOUND;
    }
    if (node_a == 0 || node_b == 0)
        fprintf(stderr, "PSTB: WARNING -- probe '%s' has a terminal at literal ground; "
                "the two-injection method cannot measure a loop broken there (same "
                "degeneracy as native .stb -- see stban.c). Reporting T=0 at every "
                "frequency.\n", (char *) probeUid);

    if (job->PACfStart <= 0.0 || job->PACfStop < job->PACfStart || job->PACpoints < 1) {
        fprintf(stderr, "PSTB: bad frequency sweep -- skipped\n");
        return E_PARMVAL;
    }

    map.n = n; map.refs = NULL; map.generic = 0; map.ns = 0; map.srows = NULL;
    if (pac_build_map(ckt, &map) || map.n != n) {
        fprintf(stderr, "PSTB: internal state-map mismatch -- skipped\n");
        pac_map_free(&map);
        return E_INTERN;
    }

    {
        int perr = pac_generic_prepare(ckt, job, &map, x0, T, trig_node_num, "PSTB");
        if (perr) { pac_map_free(&map); return perr; }
    }
    trig_idx = pac_find_trig_idx(&map, trig_node_num);
    if (trig_node_num >= 0 && trig_idx < 0) {
        fprintf(stderr, "PSTB: REFUSED -- the trigger= node (free-running "
                "orbit, card #489) is not a native capacitor's own positive "
                "terminal (or, on the generic model-charge path, a charge-carrying node row) among this circuit's states, so the "
                "period-variation terminal projection this engine needs "
                "(docs/pss-pac.md's free-running section) cannot be formed. "
                "Not running the fixed-period engine uncorrected on a "
                "moving-period orbit.\n");
        pac_map_free(&map);
        return E_UNSUPP;
    }

    br1 = TMALLOC(double, neq + 1); bi1 = TMALLOC(double, neq + 1);
    br2 = TMALLOC(double, neq + 1); bi2 = TMALLOC(double, neq + 1);
    if (!br1 || !bi1 || !br2 || !bi2) {
        if (br1) tfree(br1); if (bi1) tfree(bi1); if (br2) tfree(br2); if (bi2) tfree(bi2);
        pac_map_free(&map);
        return E_NOMEM;
    }
    memset(br1, 0, (size_t) (neq + 1) * sizeof(double)); memset(bi1, 0, (size_t) (neq + 1) * sizeof(double));
    memset(br2, 0, (size_t) (neq + 1) * sizeof(double)); memset(bi2, 0, (size_t) (neq + 1) * sizeof(double));
    /* RHS#1: unit voltage injection in the probe branch equation (stban.c:
     * ckt->CKTrhs[branch_probe] = -1.0). RHS#2: unit current injection at
     * node_b (stban.c: ckt->CKTrhs[node_b] += 1.0). */
    br1[branch_probe] = -1.0;
    if (lgmode) {
        /* lgprobe's Vi runs x -> a, the opposite orientation of a plain Vprobe a b (b = x, the upstream
         * side), so its branch row V(x)-V(a) = +1 is the same perturbation as the plain row -1, and its
         * branch current (into x) is the plain branch current negated (see the read-back below). */
        br1[branch_probe] = 1.0;
    }
    if (inj_node > 0) br2[inj_node] = 1.0;
    /* card #465 FOURTH pass, task 4 (red-proof, TEST-ONLY, default off):
     * PAC_REDPROOF_SIGNFLIP flips RHS#1's sign relative to stban.c's own
     * verbatim-matched convention (this file's header/this function's own
     * comment -- "BOTH copied verbatim (same sign conventions) from
     * stban.c's native `.stb` loop-mode"). Tian's eq. 30 combination is
     * sign-sensitive (asymmetric in v1/i1 vs v2/i2), so this must make
     * pstb_loop disagree with native `.stb` -- asserted by
     * tests/test_pac.py's check_redproof_signflip(). Never set in
     * production. */
    if (getenv("PAC_REDPROOF_SIGNFLIP"))
        br1[branch_probe] = -br1[branch_probe];

    nfreq = pac_sweep_count(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop);
    if (nfreq < 1) {
        tfree(br1); tfree(bi1); tfree(br2); tfree(bi2); pac_map_free(&map);
        return E_PARMVAL;
    }
    freqs = TMALLOC(double, nfreq);
    pac_sweep_fill(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop, freqs, nfreq);

    outrows[0] = inj_node; outrows[1] = branch_probe;
    Phi = TMALLOC(PACc, (size_t) PAC_NST(&map) * PAC_NST(&map));   /* card #555: the engine's n is PAC_NST (generic path: more rows than n) */
    Y1 = TMALLOC(PACc, (size_t) 2 * nfreq);
    Y2 = TMALLOC(PACc, (size_t) 2 * nfreq);
    Tarr = TMALLOC(PACc, nfreq);
    LGarr = TMALLOC(PACc, nfreq);

    fprintf(stderr, "PSTB: time-domain shooting-based loop gain, probe=%s, %d frequency "
                    "point%s from %.6g to %.6g Hz, n=%d reduced state%s\n",
            (char *) probeUid, nfreq, nfreq == 1 ? "" : "s", job->PACfStart, job->PACfStop,
            n, n == 1 ? "" : "s");

    t0 = SPfrontEnd->IFseconds();
    {
        /* card #555: both injections run on the same orbit and the same recorded trajectory, so the n homogeneous
         * (monodromy) columns of pass 1 are identical: the second run takes Phi from the first and propagates only
         * its forced columns.  PAC_NO_PHI_REUSE=1 runs both in full. */
        int reuse = (getenv("PAC_NO_PHI_REUSE") == NULL) && Phi != NULL;
        err = pac_engine_run(ckt, job, &map, x0, T, br1, bi1, freqs, nfreq, outrows, 2, reuse ? Phi : NULL, NULL, Y1,
                              fT, trig_idx);
        if (!err)
            err = pac_engine_run(ckt, job, &map, x0, T, br2, bi2, freqs, nfreq, outrows, 2, NULL, reuse ? Phi : NULL, Y2,
                                  fT, trig_idx);
    }
    t1 = SPfrontEnd->IFseconds();
    fprintf(stderr, "PSTB: %s, wall_time=%.3fs (%.4fs/point, 2 forcing patterns)\n",
            err ? "FAILED" : "done", t1 - t0, nfreq > 0 ? (t1 - t0) / (2 * nfreq) : 0.0);

    if (!err) {
        runDesc *pstbPlot = NULL;
        IFuid freqUid, vecuids[2];
        for (i = 0; i < nfreq; i++) {
            /* v1/i1 from RHS#1 (voltage injection), v2/i2 from RHS#2
             * (current injection); the branch-current sign flip (raw MNA
             * unknown -> "current the source delivers out of its + terminal")
             * matches stban.c's i1.re = -ckt->CKTrhs[branch_probe] exactly:
             * outrows[1]=branch_probe already IS that raw unknown's
             * sideband-0 value, so negate it here the same way. */
            PACc v1 = Y1[(size_t) 0 * nfreq + i];
            PACc i1 = pc_scale(Y1[(size_t) 1 * nfreq + i], lgmode ? 1.0 : -1.0);
            PACc v2 = Y2[(size_t) 0 * nfreq + i];
            PACc i2 = pc_scale(Y2[(size_t) 1 * nfreq + i], lgmode ? 1.0 : -1.0);
            PACc one = { 1.0, 0.0 };
            PACc diff = pc_sub(pc_mul(i1, v2), pc_mul(v1, i2));
            PACc two_diff = pc_scale(diff, 2.0);
            PACc s = pc_add(pc_add(two_diff, v1), i2);
            PACc t = pc_div(s, pc_sub(one, s));
            PACc lg = pc_scale(t, -1.0);
            if (node_a == 0 || node_b == 0) { t.re = t.im = lg.re = lg.im = 0.0; }
            Tarr[i] = t; LGarr[i] = lg;
        }

        /* card #465: the MARGINS plot is emitted BEFORE the main sweep
         * plot (opposite of the natural computation order -- margins are
         * derived FROM Tarr/LGarr, which are already fully computed by
         * this point) so that the per-frequency sweep plot ("PSTB
         * Analysis", the one a test/control script actually wants to
         * `wrdata`) is left as ngspice's CURRENT plot when this analysis
         * returns, addressable by bare vector names (loopgain/t) with no
         * need to guess or discover its auto-generated short plot name
         * (found while validating this card: that name depends on how
         * many OTHER plots -- PSS's own internal period-run plots
         * included -- were opened earlier in the SAME ngspice process,
         * which is not practically predictable from outside). */

        /* margins: log-f interpolation, phase unwrapping, worst-of-
         * several-crossings -- duplicated from stban.c's STBan() tail
         * (see that function for the derivation this mirrors verbatim). */
        {
            double *phs = TMALLOC(double, nfreq);
            double prevp = 0.0;
            PACcrossing *pmx = NULL, *gmx = NULL;
            int npm = 0, ngm = 0, capm = 0, capg = 0;
            int pmi = -1, gmi = -1;
            double pmbest = 0.0, gmbest = 0.0;
            double nan = strtod("NaN", NULL);

            for (i = 0; i < nfreq; i++) {
                double p = 180.0 / M_PI * pc_ph(Tarr[i]);
                if (i > 0) {
                    while (p - prevp > 180.0) p -= 360.0;
                    while (p - prevp < -180.0) p += 360.0;
                }
                phs[i] = p; prevp = p;
            }
            for (i = 0; i + 1 < nfreq; i++) {
                double f0 = freqs[i], f1 = freqs[i + 1];
                double m0 = 20.0 * log10(pc_mag(Tarr[i])), m1 = 20.0 * log10(pc_mag(Tarr[i + 1]));
                double p0 = phs[i], p1 = phs[i + 1];
                double lo, hi, target; int nn;
                if (m0 * m1 < 0 || (m1 == 0.0 && m0 != 0.0)) {
                    double fc = pac_interp_logf(f0, f1, m0, m1, 0.0);
                    double t2 = (0.0 - m0) / (m1 - m0);
                    double ph = p0 + t2 * (p1 - p0);
                    pmx = pac_push_crossing(pmx, &npm, &capm, fc, 180.0 + ph, 0.0);
                }
                lo = p0 < p1 ? p0 : p1; hi = p0 < p1 ? p1 : p0;
                nn = (int) ceil((lo + 180.0) / 360.0);
                while (-180.0 + 360.0 * nn <= hi) {
                    target = -180.0 + 360.0 * nn;
                    if (target != p0) {
                        double fc = pac_interp_logf(f0, f1, p0, p1, target);
                        double t2 = (target - p0) / (p1 - p0);
                        double mag = m0 + t2 * (m1 - m0);
                        gmx = pac_push_crossing(gmx, &ngm, &capg, fc, -mag, target);
                    }
                    nn++;
                }
            }
            for (i = 0; i < npm; i++) if (pmi < 0 || pmx[i].value < pmbest) { pmi = i; pmbest = pmx[i].value; }
            for (i = 0; i < ngm; i++) if (gmi < 0 || gmx[i].value < gmbest) { gmi = i; gmbest = gmx[i].value; }
            /* card #465 rebase (2026-09-29): stban.c patch 0013 (card #476) narrowed
             * its low-frequency positive-feedback check from "Re(T)<0 at the lowest
             * swept frequency alone" to arg0/mag0 thresholds, after the old form
             * false-triggered on healthy integrator-type loops. This block is a
             * verbatim duplicate of stban.c's STBan() tail (see this file's header,
             * "PSTB" paragraph, and stban.c's own comment at the same check) so it
             * is kept consistent with that NEW semantics here: freqs[] is always
             * ascending (pac_sweep_fill, fstart<=fstop enforced), so index 0 is the
             * lowest swept frequency, same as stban.c's Tarr[0]. NOTE: unlike
             * stban.c, this file's pmx[].aux is always pushed as a literal 0.0
             * (PAC's crossing detector never implemented patch 0011's rising/
             * falling classification), so no phaseMarginFalling vector is added
             * here -- adding one would just be a constant, misleading 0. */
            {
                double arg0 = (nfreq > 0) ? 180.0 / M_PI * pc_ph(Tarr[0]) : 0.0;
                double mag0 = (nfreq > 0) ? pc_mag(Tarr[0]) : 0.0;
                int low_f_posfb = (nfreq > 0) && (fabs(arg0) > 150.0) && (mag0 > 1.0);
                int low_f_posfb_low_gain = (nfreq > 0) && (fabs(arg0) > 150.0) && (mag0 <= 1.0);
                int low_f_integrator = (nfreq > 0) && (arg0 >= -120.0) && (arg0 <= -60.0);
                IFuid uids[7]; double outv[7]; runDesc *mplot = NULL; int merr;
                SPfrontEnd->IFnewUid(ckt, &uids[0], NULL, "phaseMargin", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[1], NULL, "phaseMarginFreq", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[2], NULL, "gainMargin", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[3], NULL, "gainMarginFreq", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[4], NULL, "lowFPositiveFeedback", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[5], NULL, "lowFPositiveFeedbackLowGain", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[6], NULL, "lowFIntegrator", UID_OTHER, NULL);
                outv[0] = pmi >= 0 ? pmx[pmi].value : nan;
                outv[1] = pmi >= 0 ? pmx[pmi].freq : nan;
                outv[2] = gmi >= 0 ? gmx[gmi].value : nan;
                outv[3] = gmi >= 0 ? gmx[gmi].freq : nan;
                outv[4] = low_f_posfb ? 1.0 : 0.0;
                outv[5] = low_f_posfb_low_gain ? 1.0 : 0.0;
                outv[6] = low_f_integrator ? 1.0 : 0.0;
                merr = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob, "PSTB Margins", NULL, 0, 7, uids, IF_REAL, &mplot);
                if (!merr) {
                    IFvalue refval, outdata;
                    refval.rValue = 0;
                    outdata.v.numValue = 7;
                    outdata.v.vec.rVec = outv;
                    SPfrontEnd->OUTpData(mplot, &refval, &outdata);
                    SPfrontEnd->OUTendPlot(mplot);
                }
                printf("pstb: ");
                if (pmi >= 0) printf("phase margin %.4g deg at %.6g Hz", pmx[pmi].value, pmx[pmi].freq);
                else printf("phase margin: no 0 dB crossing of T");
                if (npm > 1) printf(" (%d crossings found; the lowest margin is reported)", npm);
                printf("\n");
                printf("pstb: ");
                if (gmi >= 0) printf("gain margin %.4g dB at %.6g Hz", gmx[gmi].value, gmx[gmi].freq);
                else printf("gain margin: no -180 deg crossing of T");
                if (ngm > 1) printf(" (%d crossings found; the lowest margin is reported)", ngm);
                printf("\n");
                /* card #476's disjoint arg0-range elif chain, mirrored verbatim. */
                if (low_f_posfb) {
                    printf("pstb: WARNING: positive feedback at low frequency: loop is unstable "
                        "if |T|>1 there; phase/gain margins are not meaningful\n");
                } else if (low_f_posfb_low_gain) {
                    printf("pstb: note: positive feedback at low frequency, loop gain < 1\n");
                } else if (low_f_integrator) {
                    printf("pstb: note: integrator-like low-frequency behaviour\n");
                }
            }
            if (phs) tfree(phs);
            if (pmx) tfree(pmx);
            if (gmx) tfree(gmx);
        }

        SPfrontEnd->IFnewUid(ckt, &freqUid, NULL, "frequency", UID_OTHER, NULL);
        SPfrontEnd->IFnewUid(ckt, &vecuids[0], NULL, "loopGain", UID_OTHER, NULL);
        SPfrontEnd->IFnewUid(ckt, &vecuids[1], NULL, "T", UID_OTHER, NULL);
        err = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob, "PSTB Analysis",
                                        freqUid, IF_REAL, 2, vecuids, IF_COMPLEX, &pstbPlot);
        if (!err) {
            if (job->PACstepType != 0)
                SPfrontEnd->OUTattributes(pstbPlot, NULL, OUT_SCALE_LOG, NULL);
            for (i = 0; i < nfreq; i++) {
                IFvalue freqData, valueData;
                IFcomplex data[2];
                data[0].real = LGarr[i].re; data[0].imag = LGarr[i].im;
                data[1].real = Tarr[i].re;  data[1].imag = Tarr[i].im;
                freqData.rValue = freqs[i];
                valueData.v.numValue = 2;
                valueData.v.vec.cVec = data;
                SPfrontEnd->OUTpData(pstbPlot, &freqData, &valueData);
            }
            SPfrontEnd->OUTendPlot(pstbPlot);
        }
    }

    tfree(br1); tfree(bi1); tfree(br2); tfree(bi2); pac_map_free(&map); tfree(freqs);
    tfree(Phi); tfree(Y1); tfree(Y2); tfree(Tarr); tfree(LGarr);
    return err ? E_INTERN : OK;
}

/* ------------------------------------------------------------------ */
/* `.pxf ... method="shoot" solver="td"` (card #499): sideband-0 transfer
 * function from EVERY requested independent V/I source to ONE fixed
 * output, on the shooting PSS orbit.
 *
 * SHIPPED DESIGN (option chosen, and why -- see the card #499 report for
 * the full derivation of the alternative): the card's own brief frames
 * this as "one adjoint (transposed) solve of the time-domain PAC system
 * instead of one .pac per source", by direct analogy with docs/sens.md's
 * (card #475) two-transposed-solves-per-frequency adjoint. A genuine
 * from-scratch analogue of that here would need a BACKWARD sweep of this
 * file's own per-step recursion (transposed per-step solves against the
 * SAME J_k this file's forward recursion already reuses, via a new
 * SMPtSolve-style wrapper around Sparse's spSolveTransposed -- confirmed
 * available, TRANSPOSE=YES in this build -- composed with the transpose
 * of the per-step "history" operator H_k, which the ORDER-2 trapezoidal
 * companion model complicates: it is a 2-term recursion in y AND the
 * companion "T" array (PSSpacRecordStep()'s own g_pac_Tprev/Tnew), not a
 * simple 1-step Markov map) plus a REBUILD of every J_k on demand (this
 * file's own G(x)/C(x) is state-dependent -- pac_refactor_at_accepted()
 * already does exactly this rebuild, but only ever forward, at the
 * point the forward recursion is AT; a backward sweep needs it at
 * arbitrary points, out of order, which needs the whole state
 * trajectory x_k stored first). Each of those pieces is individually
 * tractable (the full derivation is in the card #499 report) but
 * NONE of it has any precedent or validation anywhere in this codebase,
 * and this project's OWN saltation work (docs/pss-pac.md's saltation
 * section) is a direct, sobering demonstration of what an under-derived
 * transposed/adjoint correction costs here: several passes, a real,
 * reproducible ~2.2x scaling discrepancy never root-caused, "a guessed
 * edit... would not have been trustworthy". Given that precedent, this
 * pass ships the LOWER-RISK, EQUALLY CORRECT option instead: reuse
 * pac_engine_run() (this file's own already-validated forward boundary-
 * value engine, the SAME one `.pac solver=td` uses -- ZERO changes to its
 * per-step recursion, boundary solve, or demodulation) once PER REQUESTED
 * SOURCE, each call supplying that ONE source's own isolated unit
 * stamp (pac_capture_source_b0() above) as the forcing pattern and the
 * SAME single output row as its only requested output (nout=1). This
 * still eliminates the dominant real cost "one .pac per source" would
 * pay -- a full separate PSS shooting-Newton RECONVERGENCE per source,
 * confirmed by docs/pss-pac.md's own timing note ("PSS Newton
 * convergence... was the entire cost, not this engine's own sweep") --
 * since PSS converges exactly ONCE for the whole .pxf job (same as
 * .pac/.pstb today) and only the cheap per-source forced-column sweep
 * (O(2 period reruns) each, independent of sweep length, same big-O this
 * file's header already gives .pac) repeats, K times for K sources. What
 * this does NOT do: share the Phi/monodromy computation (source-
 * independent) across sources -- each call recomputes it, a real but
 * modest cost (K x 2 period reruns instead of the theoretical 1+K) not
 * worth the extra pac_engine_run() surgery this pass chose not to risk.
 * A true backward-adjoint engine (amortising ALL sources behind one
 * extra sweep, no per-source period rerun at all) remains a valid, well-
 * scoped FUTURE improvement if the source count ever makes this pass's
 * O(K) cost the bottleneck -- not attempted this pass; see "what was
 * skipped" in the card #499 report. */
static int
pac_td_pxf(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T,
           const double *fT, int trig_node_num)
{
    PACstateMap map;
    int neq = SMPmatSize(ckt->CKTmatrix);
    int nfreq, outrow, err, trig_idx;
    double *freqs;
    PACsrcRef *srclist = NULL;
    int nsrc_total, *use_idx = NULL, nsrc_use, s;
    PACc *Y0;
    IFuid freqUid, *nameList = NULL;
    runDesc *pxfPlot = NULL;
    double t0, t1;

    if (!job->PxOutNode) {
        fprintf(stderr, "PXF-TD: no output node given -- skipped\n");
        return E_PARMVAL;
    }
    outrow = job->PxOutNode->number;
    if (outrow <= 0) {
        fprintf(stderr, "PXF-TD: output node '%s' is not in the matrix "
                "(literal ground?) -- skipped\n", job->PxOutNode->name);
        return E_PARMVAL;
    }
    if (job->PACfStart <= 0.0 || job->PACfStop < job->PACfStart || job->PACpoints < 1) {
        fprintf(stderr, "PXF-TD: bad frequency sweep -- skipped\n");
        return E_PARMVAL;
    }

    map.n = n; map.refs = NULL; map.generic = 0; map.ns = 0; map.srows = NULL;
    if (pac_build_map(ckt, &map) || map.n != n) {
        fprintf(stderr, "PXF-TD: internal state-map mismatch -- skipped\n");
        pac_map_free(&map);
        return E_INTERN;
    }
    {
        int perr = pac_generic_prepare(ckt, job, &map, x0, T, trig_node_num, "PXF-TD");
        if (perr) { pac_map_free(&map); return perr; }
    }
    trig_idx = pac_find_trig_idx(&map, trig_node_num);
    if (trig_node_num >= 0 && trig_idx < 0) {
        fprintf(stderr, "PXF-TD: REFUSED -- the trigger= node (free-running "
                "orbit, card #489) is not a native capacitor's own positive "
                "terminal (or, on the generic model-charge path, a charge-carrying node row) among this circuit's states, so the "
                "period-variation terminal projection this engine needs "
                "(docs/pss-pac.md's free-running section) cannot be formed. "
                "Not running the fixed-period engine uncorrected on a "
                "moving-period orbit.\n");
        pac_map_free(&map);
        return E_UNSUPP;
    }

    nsrc_total = pac_build_src_list(ckt, &srclist);
    if (nsrc_total <= 0) {
        fprintf(stderr, "PXF-TD: no independent V/I source instances found in "
                "the circuit -- skipped\n");
        pac_map_free(&map);
        return E_PARMVAL;
    }

    use_idx = TMALLOC(int, nsrc_total);
    nsrc_use = 0;
    if (job->PxfSourcesAll || job->PxfSourceCount == 0) {
        for (s = 0; s < nsrc_total; s++) use_idx[nsrc_use++] = s;
    } else {
        int k;
        for (k = 0; k < job->PxfSourceCount; k++) {
            GENinstance *inst = CKTfndDev(ckt, job->PxfSourceNames[k]);
            int found = -1;
            if (inst)
                for (s = 0; s < nsrc_total; s++)
                    if (srclist[s].inst == inst) { found = s; break; }
            if (found < 0) {
                fprintf(stderr, "PXF-TD: WARNING -- sources=\"%s\" is not an "
                        "independent V/I source instance in this circuit; "
                        "ignored\n", (char *) job->PxfSourceNames[k]);
                continue;
            }
            use_idx[nsrc_use++] = found;
        }
        if (nsrc_use == 0) {
            fprintf(stderr, "PXF-TD: no requested sources= matched an "
                    "independent V/I source -- skipped\n");
            pac_map_free(&map); tfree(srclist); tfree(use_idx);
            return E_PARMVAL;
        }
    }

    nfreq = pac_sweep_count(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop);
    if (nfreq < 1) { pac_map_free(&map); tfree(srclist); tfree(use_idx); return E_PARMVAL; }
    freqs = TMALLOC(double, nfreq);
    pac_sweep_fill(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop, freqs, nfreq);

    Y0 = TMALLOC(PACc, (size_t) nsrc_use * nfreq);
    nameList = TMALLOC(IFuid, nsrc_use);

    fprintf(stderr, "PXF-TD: time-domain shooting-based PXF, %d frequency "
                    "point%s from %.6g to %.6g Hz, n=%d reduced state%s, "
                    "output node '%s' (row %d), %d source%s\n",
            nfreq, nfreq == 1 ? "" : "s", job->PACfStart, job->PACfStop, n,
            n == 1 ? "" : "s", job->PxOutNode->name, outrow,
            nsrc_use, nsrc_use == 1 ? "" : "s");

    /* card #504: adjoint=1 (explicit), adjoint=0 (explicit, force v1), or
     * unset (auto: adjoint once nsrc_use is "above a small number" -- see
     * PXF_ADJOINT_AUTO_THRESHOLD). Free-running (trig_idx>=0) always
     * forces v1 -- pac_engine_run_adjoint() would refuse it anyway (see
     * that function's header), checked here too just to skip the attempt
     * and its own log line. */
    {
        int want_adjoint = job->PxfAdjointSet ? (job->PxfAdjoint != 0)
                                               : (nsrc_use > PXF_ADJOINT_AUTO_THRESHOLD);
        int adj_err = PAC_ADJOINT_REFUSED;

        if (trig_idx >= 0 && getenv("PAC_FREERUN_PROJ") && !getenv("PAC_FREERUN_NO_PROJ")) {
            if (want_adjoint)
                fprintf(stderr, "PXF-TD: adjoint engine REFUSED (the opt-in free-running period-variation projection, PAC_FREERUN_PROJ, is not "
                        "transposed) -- using the per-source engine.\n");
            want_adjoint = 0;   /* card #557: the opt-in projection is not transposed (default free-running forward = fixed-period engine) */
        }

        t0 = SPfrontEnd->IFseconds();
        err = OK;

        if (want_adjoint) {
            fprintf(stderr, "PXF-TD: adjoint=1 (%s) -- true backward/transposed engine, "
                    "one sweep for all %d source%s\n",
                    job->PxfAdjointSet ? "explicit" : "auto: source count > threshold",
                    nsrc_use, nsrc_use == 1 ? "" : "s");
            adj_err = pac_engine_run_adjoint(ckt, job, &map, x0, T, freqs, nfreq, outrow,
                                              srclist, nsrc_total, use_idx, nsrc_use, Y0, trig_idx);
            if (adj_err == PAC_ADJOINT_OK) {
                for (s = 0; s < nsrc_use; s++) {
                    char nm[160];
                    (void) snprintf(nm, sizeof(nm), "pxf_%s", (char *) srclist[use_idx[s]].name);
                    SPfrontEnd->IFnewUid(ckt, &nameList[s], NULL, nm, UID_OTHER, NULL);
                }
                t1 = SPfrontEnd->IFseconds();
                fprintf(stderr, "PXF-TD: adjoint done, wall_time=%.3fs "
                        "(%d source%s, %d frequenc%s -- one backward sweep, not %d)\n",
                        t1 - t0, nsrc_use, nsrc_use == 1 ? "" : "s",
                        nfreq, nfreq == 1 ? "y" : "ies", nsrc_use);
            } else if (adj_err == PAC_ADJOINT_REFUSED) {
                fprintf(stderr, "PXF-TD: adjoint engine REFUSED this circuit (saltation "
                        "events seen, or a free-running orbit -- see pacshoot.c's "
                        "pac_engine_run_adjoint() header) -- falling back to the v1 "
                        "per-source engine.\n");
            } else {
                fprintf(stderr, "PXF-TD: adjoint engine FAILED -- falling back to the v1 "
                        "per-source engine.\n");
            }
        }

        if (!want_adjoint || adj_err != PAC_ADJOINT_OK) {
            for (s = 0; s < nsrc_use; s++) {
                double *br = TMALLOC(double, neq + 1);
                double *bi = TMALLOC(double, neq + 1);
                int one_out[1];
                char nm[160];

                one_out[0] = outrow;
                if (!br || !bi) {
                    if (br) tfree(br);
                    if (bi) tfree(bi);
                    err = 1;
                    break;
                }
                memset(br, 0, (size_t) (neq + 1) * sizeof(double));
                memset(bi, 0, (size_t) (neq + 1) * sizeof(double));

                if (pac_capture_source_b0(ckt, srclist, nsrc_total, use_idx[s], br, bi, neq)) {
                    fprintf(stderr, "PXF-TD: failed to isolate source '%s' -- "
                            "skipped\n", (char *) srclist[use_idx[s]].name);
                    err = 1;
                } else {
                    PACc *Ycol = TMALLOC(PACc, nfreq);
                    int e2 = pac_engine_run(ckt, job, &map, x0, T, br, bi, freqs, nfreq,
                                             one_out, 1, NULL, NULL, Ycol, fT, trig_idx);
                    if (e2) {
                        fprintf(stderr, "PXF-TD: engine failed for source '%s'\n",
                                (char *) srclist[use_idx[s]].name);
                        err = 1;
                    } else {
                        int f;
                        for (f = 0; f < nfreq; f++)
                            Y0[(size_t) s * nfreq + f] = Ycol[f];
                    }
                    tfree(Ycol);
                }

                (void) snprintf(nm, sizeof(nm), "pxf_%s", (char *) srclist[use_idx[s]].name);
                SPfrontEnd->IFnewUid(ckt, &nameList[s], NULL, nm, UID_OTHER, NULL);

                tfree(br); tfree(bi);
            }
            t1 = SPfrontEnd->IFseconds();
            fprintf(stderr, "PXF-TD: %s, wall_time=%.3fs (%.4fs/source)\n",
                    err ? "completed with errors" : "done", t1 - t0,
                    nsrc_use > 0 ? (t1 - t0) / nsrc_use : 0.0);
        }
    }

    SPfrontEnd->IFnewUid(ckt, &freqUid, NULL, "frequency", UID_OTHER, NULL);
    {
        int e3 = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob,
                        "PXF Analysis (time domain, sideband 0)",
                        freqUid, IF_REAL, nsrc_use, nameList, IF_COMPLEX, &pxfPlot);
        if (!e3) {
            int i;
            if (job->PACstepType != 0)
                SPfrontEnd->OUTattributes(pxfPlot, NULL, OUT_SCALE_LOG, NULL);
            for (i = 0; i < nfreq; i++) {
                IFvalue freqData, valueData;
                IFcomplex *data = TMALLOC(IFcomplex, nsrc_use);
                int j;
                for (j = 0; j < nsrc_use; j++) {
                    data[j].real = Y0[(size_t) j * nfreq + i].re;
                    data[j].imag = Y0[(size_t) j * nfreq + i].im;
                }
                freqData.rValue = freqs[i];
                valueData.v.numValue = nsrc_use;
                valueData.v.vec.cVec = data;
                SPfrontEnd->OUTpData(pxfPlot, &freqData, &valueData);
                FREE(data);
            }
            SPfrontEnd->OUTendPlot(pxfPlot);
        } else if (!err) {
            err = 1;
        }
    }

    pac_map_free(&map); tfree(srclist); tfree(use_idx); tfree(freqs); tfree(Y0); tfree(nameList);
    return err ? E_INTERN : OK;
}

/* ------------------------------------------------------------------ */

int
pac_pstb_run_ex(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T,
                 const double *fT, int trig_node_num)
{
    int err = OK, a;
    static const double rungs[] = { 1e-8, 1e-7, 1e-6, 0.0 };
    const char *ae = getenv("PAC_RECORD_RELTOL_AUTO");
    int nr = (job->PACrecordReltol > 0.0 || (getenv("PAC_RECORD_RELTOL") && atof(getenv("PAC_RECORD_RELTOL")) > 0.0)
              || (ae && atoi(ae) == 0) || pac_ckt_has_events(ckt)) ? 1 : (int) (sizeof rungs / sizeof rungs[0]);

    /* card #625: each analysis runs at the first rung that completes; plots are only begun after its engine
     * runs succeed, so a retry leaves no partial plot. The prepass state is redone at the new rung. */
    for (a = 0; a < 4; a++) {
        int on = (a == 0) ? (job->PSSdoPAC && job->PACsolver == PAC_SOLVER_TD)
               : (a == 1) ? job->PSSdoPSTB
               : (a == 2) ? (job->PSSdoPXF && job->PACsolver == PAC_SOLVER_TD)
               : (job->PSSdoPnoise && job->PACsolver == PAC_SOLVER_TD);
        static const char *nm[] = { "pac", "pstb", "pxf", "pnoise" };
        int r, e = OK;
        if (!on) continue;
        for (r = 0; r < nr; r++) {
            g_pac_auto_reltol = (nr > 1) ? rungs[r] : 0.0;
            g_pac_rec_tstep = 0;
            if (a == 0)      e = pac_td_pac(ckt, job, x0, n, T, fT, trig_node_num);
            else if (a == 1) e = pac_td_pstb(ckt, job, x0, n, T, fT, trig_node_num);
            else if (a == 2) e = pac_td_pxf(ckt, job, x0, n, T, fT, trig_node_num);
            else             e = pac_td_pnoise(ckt, job, x0, n, T, fT, trig_node_num);
            if (!e || !g_pac_rec_tstep || r + 1 >= nr || (rungs[r + 1] == 0.0 && ckt->CKTreltol <= rungs[r]))
                break;
            fprintf(stdout, "Note: %s: the recording pass at reltol %g stopped (timestep too small); "
                    "retrying at %s (card #625; set pacreltol= to choose).\n", nm[a], rungs[r],
                    rungs[r + 1] > 0.0 ? (rungs[r + 1] == 1e-7 ? "1e-7" : "1e-6") : "the simulator reltol");
            pac_init_free();
#ifdef XSPICE
            pac_hid_free();
#endif
        }
        g_pac_auto_reltol = 0.0;
        if (e) err = e;
    }
    pac_init_free();   /* card #542 */
#ifdef XSPICE
    pac_hid_free();    /* card #557 */
#endif
    return err;
}

int
pac_pstb_run(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T)
{
    return pac_pstb_run_ex(ckt, job, x0, n, T, NULL, -1);
}
