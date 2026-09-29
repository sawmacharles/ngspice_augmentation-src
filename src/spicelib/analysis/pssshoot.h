/**********
card #458 (plan T4.1-T4.3, second pass): shared declarations between
pssshoot.c, psseigen.c and the two small dctran.c hooks that let the
monodromy computation force an IDENTICAL accepted-step time grid across
the baseline and perturbed period runs of one Newton iteration.

SUSPECTED CAUSE 2 (second-pass card, ground truth from the parent,
2026-09-28): the first pass's finite-difference monodromy was noisy
because each perturbed period run re-derives its OWN adaptive step
sequence (LTE-driven), so a tiny x0 perturbation can land the switching
edges at slightly different times in each of the 2n+1 period runs a
Newton step needs -- "switching-edge timing noise ~ step size", the same
effect flagged for the transient FRA. Recomputing G(t)/C(t) analytically
at every accepted step (the fully "replace" option in the card, S_k
propagated with the transient's own Jacobians per research/...:Sec.4.6)
needs deep, version-sensitive access to ngspice's live NIiter/companion-
model internals (order changes, predictor/corrector mode bits, trapezoidal
vs first-step backward-Euler, breakpoint order resets...) with real risk of
a silently-wrong recursion that LOOKS validated by luck on a 1-state
circuit but is wrong on R3/R4. Given the time available for this card,
the "supplement" option was implemented instead, which removes the SAME
root cause with far less internals risk: record every accepted time from
the Newton iteration's own baseline P(x0) run, then FORCE those exact
times as breakpoints on every other period run in that iteration (the 2n
perturbed runs and any re-baseline/trial runs) via ngspice's own, public,
untouched CKTsetBreak(). Because ngspice's breakpoint logic always forces
a step to land exactly on the next breakpoint once the natural step would
cross it, this makes the perturbed runs retrace the IDENTICAL time grid
as the baseline for small perturbations of a smooth (tanh-switch, design
ruling #447 A) circuit -- eliminating the adaptive-step mismatch that was
the diagnosed noise source, without hand-deriving or reimplementing any
part of the trapezoidal/Gear companion recursion. What this does NOT do:
it is still a finite-difference monodromy (bias ~ eps, not eliminated),
not the analytic S_k/J_k/C_k propagation the card's "replace" option
describes; see the card #458 report for what was and was not attempted of
that fuller design, and psstest/tools/pss_ref.py for the independent
check this was validated against on R2.
**********/

#ifndef ngspice_PSSSHOOT_H
#define ngspice_PSSSHOOT_H

#include "ngspice/cktdefs.h"

/* ---- psseigen.c: full eigenvalue decomposition of the (small) monodromy
 * matrix, replacing power iteration (SUSPECTED CAUSE 1). ---- */

typedef struct {
    double re, im;
} PSScomplex;

typedef struct {
    int valid;
    double spectral_radius;   /* max_i |lambda_i| */
    PSScomplex dominant;      /* the multiplier achieving it */
    int unstable;             /* spectral_radius >= 1.0 */
    int period_doubling;      /* a near-unit-circle multiplier is real and negative */
    int oscillation;          /* a near-unit-circle multiplier is complex (non-real) */
} PSSfloquet;

int  pss_eigenvalues(const double *Ain, int n, PSScomplex *eig);
void pss_classify_floquet(const PSScomplex *eig, int n, PSSfloquet *out);

/* ---- breakpoint-forcing hooks, called from dctran.c (SUSPECTED CAUSE 2).
 * Both are cheap no-ops (a single integer test) unless a `.pss
 * method=shoot` analysis is actively using them; every other analysis
 * (including the existing driven-mode `.pss`, and all of `.tran`) is
 * unaffected -- these flags default to 0/NULL and nothing in dctran.c's
 * own control flow changes. ---- */

extern int PSSbreakRecordActive;   /* dctran.c: record every accepted time */
extern int PSSbreakInstallActive;  /* dctran.c: force the recorded times as breakpoints */

void PSSbreakRecordStep(CKTcircuit *ckt);   /* called right after CKTaccept() */
void PSSbreakInstallHook(CKTcircuit *ckt);  /* called right after the restart-path breaks reset */

/* ---- card #458 FOURTH pass (plan T4.2): analytic (discretise-then-
 * linearise) sensitivity monodromy, replacing finite differences. See
 * pssshoot.c's PSSsensRecordStep() header comment for the full S_k/T_k
 * recursion and its derivation from niinteg.c/capload.c/indload.c. Same
 * no-op-unless-active pattern as the breakpoint hooks above: PSSsensActive
 * defaults to 0, so every other analysis (including method=shoot's own
 * default finite-difference Jacobian) is completely unaffected. ---- */

extern int PSSsensActive;   /* dctran.c: run the per-step sensitivity recursion */

void PSSsensRecordStep(CKTcircuit *ckt);    /* called right after CKTaccept(), same site as PSSbreakRecordStep */

/* ---- card #458 FIFTH pass (plan T4.3): per-state ripple tracking, used to
 * make the Newton convergence test relative to each state's own ripple
 * (task 1) instead of only its magnitude. Same no-op-unless-active pattern
 * as the hooks above. See pssshoot.c's PSSrippleRecordStep() header
 * comment. ---- */

extern int PSSrippleRecordActive;

void PSSrippleRecordStep(CKTcircuit *ckt);   /* called right after CKTaccept(), same site as the hooks above */

/* ---- card #489 (plan T6.2): free-running (autonomous) shooting's
 * Poincare-section trigger detector. Same no-op-unless-active pattern as
 * the hooks above. See pssshoot.c's PSStrigRecordStep() header comment
 * for the crossing-detection/bracketing it does. ---- */

extern int PSStrigActive;   /* dctran.c: watch the trigger node for a crossing */

void PSStrigRecordStep(CKTcircuit *ckt);   /* called right after CKTaccept(), same site as the hooks above */

#endif
