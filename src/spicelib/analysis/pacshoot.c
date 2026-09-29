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
  - No PXF/adjoint, no pnoise on this engine.
  - No device-based (gain-nulling) PSTB probe, no probe=inst:terminal
    auto-insert for PSTB (probe=<Vname>, an existing 0V V-source, only).
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

#include "cap/capdefs.h"
#include "ind/inddefs.h"

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
    int n;
    PACstateRef *refs;
} PACstateMap;

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

/* ------------------------------------------------------------------ */
/* One-period re-run, restricted stack-local TRANan trick -- duplicated
 * from pssshoot.c's (file-static) pss_run_period() rather than exposing
 * it, for the same "don't change that file's visibility for this card"
 * reason as PACstateMap above. Caller must have already written x0 via
 * pac_state_set_ic(). */
static int
pac_run_period(CKTcircuit *ckt, double T, const char *jobname)
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
    ckt->CKTmode     |= MODEUIC;

    err = DCtran(ckt, 1);

    if (localjob.TRANplot)
        SPfrontEnd->OUTendPlot(localjob.TRANplot);

    ckt->CKTcurJob = saveJob;
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
    g_pac_Cprev = g_pac_Cnew = NULL;
    g_pac_Tprev = g_pac_Tnew = NULL;
    g_pac_omega = NULL; g_pac_useforce = NULL; g_pac_demod = NULL;
    g_pac_rhs_re = g_pac_rhs_im = g_pac_scratch = NULL;
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
    int n = map->n, neq, d, m;

    pac_free();

    neq = SMPmatSize(ckt->CKTmatrix);
    g_pac_map = map;
    g_pac_n = n;
    g_pac_neq = neq;
    g_pac_M = M;
    g_pac_T = T;
    g_pac_demod_active = demod_active;
    g_pac_tprev = 0.0;

    g_pac_Cprev = TMALLOC(PACc, (size_t) (neq + 1) * M);
    g_pac_Cnew  = TMALLOC(PACc, (size_t) (neq + 1) * M);
    g_pac_Tprev = TMALLOC(PACc, (size_t) n * M);
    g_pac_Tnew  = TMALLOC(PACc, (size_t) n * M);
    g_pac_omega = TMALLOC(double, M);
    g_pac_useforce = TMALLOC(int, M);
    g_pac_rhs_re = TMALLOC(double, (size_t) (neq + 1));
    g_pac_rhs_im = TMALLOC(double, (size_t) (neq + 1));
    g_pac_scratch = TMALLOC(double, (size_t) (neq + 1));
    g_pac_br = (double *) br;   /* not owned, caller keeps alive for this pass */
    g_pac_bi = (double *) bi;
    if (demod_active)
        g_pac_demod = TMALLOC(PACc, (size_t) neq * M);

    if (!g_pac_Cprev || !g_pac_Cnew || !g_pac_Tprev || !g_pac_Tnew || !g_pac_omega ||
        !g_pac_useforce || !g_pac_rhs_re || !g_pac_rhs_im || !g_pac_scratch ||
        (demod_active && !g_pac_demod)) {
        pac_free();
        return 1;
    }

    memset(g_pac_Cprev, 0, (size_t) (neq + 1) * M * sizeof(PACc));
    memset(g_pac_Tprev, 0, (size_t) n * M * sizeof(PACc));
    if (demod_active)
        memset(g_pac_demod, 0, (size_t) neq * M * sizeof(PACc));
    for (m = 0; m < M; m++) {
        g_pac_omega[m] = omega[m];
        g_pac_useforce[m] = useforce[m];
    }
    for (d = 0; d < n; d++) {
        int row = pac_state_row(&map->refs[d]);
        for (m = 0; m < M; m++)
            pac_set(g_pac_Cprev, row, M, m, ic0[(size_t) d * M + m]);
    }
    g_pac_failed = 0;
    return 0;
}

void
PSSpacRecordStep(CKTcircuit *ckt)
{
    int n, neq, M, d, m, row, order;
    double ag0, ag1, t_new, h;

    if (g_pac_failed || !g_pac_map)
        return;
    if (ckt->CKTtime <= 0.0)
        return;   /* t=0 base case already set by pac_reset() */

    if (ckt->CKTintegrateMethod != TRAPEZOIDAL || ckt->CKTorder < 1 || ckt->CKTorder > 2) {
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

        SMPsolve(ckt->CKTmatrix, g_pac_rhs_re, g_pac_scratch);
        SMPsolve(ckt->CKTmatrix, g_pac_rhs_im, g_pac_scratch);
        for (row = 1; row <= neq; row++) {
            PACc v; v.re = g_pac_rhs_re[row]; v.im = g_pac_rhs_im[row];
            pac_set(g_pac_Cnew, row, M, m, v);
        }
    }

    for (d = 0; d < n; d++) {
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

    if (g_pac_demod_active) {
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

    { PACc *tmp = g_pac_Cprev; g_pac_Cprev = g_pac_Cnew; g_pac_Cnew = tmp; }
    { PACc *tmp = g_pac_Tprev; g_pac_Tprev = g_pac_Tnew; g_pac_Tnew = tmp; }
    g_pac_tprev = t_new;
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
static int
pac_capture_b0(CKTcircuit *ckt, double *br, double *bi, int neq)
{
    int i;
    long savedMode = ckt->CKTmode;

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
    if (CKTload(ckt)) { ckt->CKTmode = savedMode; return 1; }
    ckt->CKTomega = 1.0;
    ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODEAC;
    if (CKTacLoad(ckt)) { ckt->CKTmode = savedMode; return 1; }

    for (i = 1; i <= neq; i++) {
        br[i] = ckt->CKTrhs[i];
        bi[i] = ckt->CKTirhs[i];
    }
    ckt->CKTmode = savedMode;
    return 0;
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
/* Core engine: given ONE forcing pattern (br/bi, shared across the whole
 * sweep) and NFREQ frequencies, returns Phi (n*n, real part only -- see
 * this file's header, the homogeneous columns never leave the real
 * subspace) and, for every requested output ROW (1-based, e.g. a node or
 * a branch equation) and frequency, the demodulated sideband-0 response
 * Y0. Runs exactly 2 period re-runs total (batched over ALL frequencies),
 * per this file's header design. Returns 0 on success. */
static int
pac_engine_run(CKTcircuit *ckt, PACstateMap *map, const double *x0, double T,
                const double *br, const double *bi,
                const double *freqs, int nfreq,
                const int *outrows, int nout,
                PACc *Phi_out /* n*n, may be NULL */,
                PACc *Y0_out  /* nout*nfreq */)
{
    int n = map->n, neq = SMPmatSize(ckt->CKTmatrix);
    int M1 = n + nfreq, f, d, j, err;
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
    for (d = 0; d < n; d++) {
        omega1[d] = 0.0;
        useforce1[d] = 0;
        for (j = 0; j < n; j++) {
            PACc v = { (d == j) ? 1.0 : 0.0, 0.0 };
            ic0_1[(size_t) j * M1 + d] = v;   /* column d: unit at reduced-state row d */
        }
    }
    for (f = 0; f < nfreq; f++) {
        omega1[n + f] = 2.0 * M_PI * freqs[f];
        useforce1[n + f] = 1;
    }

    pac_state_set_ic(map, x0);
    if (pac_reset(ckt, map, M1, br, bi, omega1, useforce1, ic0_1, 0, T)) {
        err = 1; goto cleanup;
    }
    PSSpacActive = 1;
    err = pac_run_period(ckt, T, "TD-shoot recording pass 1 (monodromy + forced)");
    PSSpacActive = 0;
    if (err || g_pac_failed) { err = 1; goto cleanup; }

    for (d = 0; d < n; d++) {
        int row_d = pac_state_row(&map->refs[d]);
        for (j = 0; j < n; j++)
            Phi[(size_t) d * n + j] = pac_get(g_pac_Cprev, row_d, M1, j);
        for (f = 0; f < nfreq; f++)
            z[(size_t) d * nfreq + f] = pac_get(g_pac_Cprev, row_d, M1, n + f);
    }
    if (Phi_out) memcpy(Phi_out, Phi, (size_t) n * n * sizeof(PACc));
    pac_free();

    /* dense boundary solve per frequency: (e^{jwT} I - Phi) xi = z(w) */
    for (f = 0; f < nfreq; f++) {
        PACc *A = TMALLOC(PACc, (size_t) n * n);
        PACc *rhs = TMALLOC(PACc, n);
        PACc *sol = TMALLOC(PACc, n);
        double w = 2.0 * M_PI * freqs[f];
        double wT = w * T;
        PACc ephase = { cos(wT), sin(wT) };
        if (!A || !rhs || !sol) { if (A) tfree(A); if (rhs) tfree(rhs); if (sol) tfree(sol); err = 1; goto cleanup; }
        for (d = 0; d < n; d++)
            for (j = 0; j < n; j++) {
                PACc v = pc_scale(Phi[(size_t) d * n + j], -1.0);
                if (d == j) v = pc_add(v, ephase);
                A[(size_t) d * n + j] = v;
            }
        for (d = 0; d < n; d++) rhs[d] = z[(size_t) d * nfreq + f];
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
        err = pac_run_period(ckt, T, "TD-shoot recording pass 2 (re-propagate)");
        PSSpacActive = 0;
        if (!err && !g_pac_failed) {
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
/* `.pac ... method=shoot solver=td`: sideband-0 response at every node,
 * forcing = the netlist AC-source stimulus B0. */
static int
pac_td_pac(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T)
{
    PACstateMap map;
    int neq = SMPmatSize(ckt->CKTmatrix);
    double *br, *bi;
    int nfreq, *outrows, i, numNames, err;
    double *freqs;
    PACc *Y0;
    IFuid freqUid, *nameList = NULL;
    runDesc *pacPlot = NULL;
    double t0, t1;

    if (job->PACfStart <= 0.0 || job->PACfStop < job->PACfStart || job->PACpoints < 1) {
        fprintf(stderr, "PAC-TD: bad frequency sweep -- skipped\n");
        return E_PARMVAL;
    }

    map.n = n; map.refs = NULL;
    if (pac_build_map(ckt, &map) || map.n != n) {
        fprintf(stderr, "PAC-TD: internal state-map mismatch -- skipped\n");
        if (map.refs) tfree(map.refs);
        return E_INTERN;
    }

    br = TMALLOC(double, neq + 1);
    bi = TMALLOC(double, neq + 1);
    if (!br || !bi) { if (br) tfree(br); if (bi) tfree(bi); tfree(map.refs); return E_NOMEM; }
    memset(br, 0, (size_t) (neq + 1) * sizeof(double));
    memset(bi, 0, (size_t) (neq + 1) * sizeof(double));
    if (pac_capture_b0(ckt, br, bi, neq)) {
        fprintf(stderr, "PAC-TD: failed to capture the netlist AC-source stimulus\n");
        tfree(br); tfree(bi); tfree(map.refs);
        return E_INTERN;
    }

    nfreq = pac_sweep_count(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop);
    if (nfreq < 1) { tfree(br); tfree(bi); tfree(map.refs); return E_PARMVAL; }
    freqs = TMALLOC(double, nfreq);
    pac_sweep_fill(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop, freqs, nfreq);

    if (CKTnames(ckt, &numNames, &nameList) || numNames != neq) {
        fprintf(stderr, "PAC-TD: CKTnames() mismatch -- skipped\n");
        tfree(br); tfree(bi); tfree(map.refs); tfree(freqs);
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
    err = pac_engine_run(ckt, &map, x0, T, br, bi, freqs, nfreq, outrows, numNames, NULL, Y0);
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

    tfree(br); tfree(bi); tfree(map.refs); tfree(freqs); tfree(outrows); tfree(Y0); tfree(nameList);
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
pac_td_pstb(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T)
{
    PACstateMap map;
    int neq = SMPmatSize(ckt->CKTmatrix);
    GENinstance *inst;
    int branch_probe, node_a, node_b;
    double *br1, *bi1, *br2, *bi2;
    int nfreq, i, err;
    double *freqs;
    int outrows[2];   /* node_b, branch_probe */
    PACc *Phi, *Y1, *Y2;
    PACc *Tarr, *LGarr;
    double t0, t1;

    IFuid probeUid = job->PSTBprobeName;
    if (!probeUid) {
        fprintf(stderr, "PSTB: no probe= given -- skipped\n");
        return E_PARMVAL;
    }
    inst = CKTfndDev(ckt, probeUid);
    if (!inst || inst->GENmodPtr->GENmodType != CKTtypelook("Vsource")) {
        fprintf(stderr, "PSTB: probe '%s' must be an existing 0V V-source (device-based "
                "and inst:terminal probes are not supported by this engine -- card #465 "
                "scope cut)\n", (char *) probeUid);
        return E_NOTFOUND;
    }
    branch_probe = CKTfndBranch(ckt, probeUid);
    if (branch_probe == 0) {
        fprintf(stderr, "PSTB: probe '%s' has no branch equation\n", (char *) probeUid);
        return E_NOTFOUND;
    }
    node_a = GENnode(inst)[0];
    node_b = GENnode(inst)[1];
    if (node_a == 0 || node_b == 0)
        fprintf(stderr, "PSTB: WARNING -- probe '%s' has a terminal at literal ground; "
                "the two-injection method cannot measure a loop broken there (same "
                "degeneracy as native .stb -- see stban.c). Reporting T=0 at every "
                "frequency.\n", (char *) probeUid);

    if (job->PACfStart <= 0.0 || job->PACfStop < job->PACfStart || job->PACpoints < 1) {
        fprintf(stderr, "PSTB: bad frequency sweep -- skipped\n");
        return E_PARMVAL;
    }

    map.n = n; map.refs = NULL;
    if (pac_build_map(ckt, &map) || map.n != n) {
        fprintf(stderr, "PSTB: internal state-map mismatch -- skipped\n");
        if (map.refs) tfree(map.refs);
        return E_INTERN;
    }

    br1 = TMALLOC(double, neq + 1); bi1 = TMALLOC(double, neq + 1);
    br2 = TMALLOC(double, neq + 1); bi2 = TMALLOC(double, neq + 1);
    if (!br1 || !bi1 || !br2 || !bi2) {
        if (br1) tfree(br1); if (bi1) tfree(bi1); if (br2) tfree(br2); if (bi2) tfree(bi2);
        tfree(map.refs);
        return E_NOMEM;
    }
    memset(br1, 0, (size_t) (neq + 1) * sizeof(double)); memset(bi1, 0, (size_t) (neq + 1) * sizeof(double));
    memset(br2, 0, (size_t) (neq + 1) * sizeof(double)); memset(bi2, 0, (size_t) (neq + 1) * sizeof(double));
    /* RHS#1: unit voltage injection in the probe branch equation (stban.c:
     * ckt->CKTrhs[branch_probe] = -1.0). RHS#2: unit current injection at
     * node_b (stban.c: ckt->CKTrhs[node_b] += 1.0). */
    br1[branch_probe] = -1.0;
    if (node_b > 0) br2[node_b] = 1.0;

    nfreq = pac_sweep_count(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop);
    if (nfreq < 1) {
        tfree(br1); tfree(bi1); tfree(br2); tfree(bi2); tfree(map.refs);
        return E_PARMVAL;
    }
    freqs = TMALLOC(double, nfreq);
    pac_sweep_fill(job->PACstepType, job->PACpoints, job->PACfStart, job->PACfStop, freqs, nfreq);

    outrows[0] = node_b; outrows[1] = branch_probe;
    Phi = TMALLOC(PACc, (size_t) n * n);
    Y1 = TMALLOC(PACc, (size_t) 2 * nfreq);
    Y2 = TMALLOC(PACc, (size_t) 2 * nfreq);
    Tarr = TMALLOC(PACc, nfreq);
    LGarr = TMALLOC(PACc, nfreq);

    fprintf(stderr, "PSTB: time-domain shooting-based loop gain, probe=%s, %d frequency "
                    "point%s from %.6g to %.6g Hz, n=%d reduced state%s\n",
            (char *) probeUid, nfreq, nfreq == 1 ? "" : "s", job->PACfStart, job->PACfStop,
            n, n == 1 ? "" : "s");

    t0 = SPfrontEnd->IFseconds();
    err = pac_engine_run(ckt, &map, x0, T, br1, bi1, freqs, nfreq, outrows, 2, Phi, Y1);
    if (!err)
        err = pac_engine_run(ckt, &map, x0, T, br2, bi2, freqs, nfreq, outrows, 2, NULL, Y2);
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
            PACc i1 = pc_scale(Y1[(size_t) 1 * nfreq + i], -1.0);
            PACc v2 = Y2[(size_t) 0 * nfreq + i];
            PACc i2 = pc_scale(Y2[(size_t) 1 * nfreq + i], -1.0);
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
            {
                IFuid uids[4]; double outv[4]; runDesc *mplot = NULL; int merr;
                SPfrontEnd->IFnewUid(ckt, &uids[0], NULL, "phaseMargin", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[1], NULL, "phaseMarginFreq", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[2], NULL, "gainMargin", UID_OTHER, NULL);
                SPfrontEnd->IFnewUid(ckt, &uids[3], NULL, "gainMarginFreq", UID_OTHER, NULL);
                outv[0] = pmi >= 0 ? pmx[pmi].value : nan;
                outv[1] = pmi >= 0 ? pmx[pmi].freq : nan;
                outv[2] = gmi >= 0 ? gmx[gmi].value : nan;
                outv[3] = gmi >= 0 ? gmx[gmi].freq : nan;
                merr = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob, "PSTB Margins", NULL, 0, 4, uids, IF_REAL, &mplot);
                if (!merr) {
                    IFvalue refval, outdata;
                    refval.rValue = 0;
                    outdata.v.numValue = 4;
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

    tfree(br1); tfree(bi1); tfree(br2); tfree(bi2); tfree(map.refs); tfree(freqs);
    tfree(Phi); tfree(Y1); tfree(Y2); tfree(Tarr); tfree(LGarr);
    return err ? E_INTERN : OK;
}

/* ------------------------------------------------------------------ */

int
pac_pstb_run(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T)
{
    int err = OK;

    if (job->PSSdoPAC && job->PACsolver == PAC_SOLVER_TD) {
        int e = pac_td_pac(ckt, job, x0, n, T);
        if (e) err = e;
    }
    if (job->PSSdoPSTB) {
        int e = pac_td_pstb(ckt, job, x0, n, T);
        if (e) err = e;
    }
    return err;
}
