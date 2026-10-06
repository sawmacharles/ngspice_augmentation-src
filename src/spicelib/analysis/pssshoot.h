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
/* card #505: exact step-sequence REPLAY (recorded time AND integration order
 * per accepted step, no LTE control) for forced period runs; see
 * pssshoot.c's PSSreplayStep() header comment. */
extern int PSSreplayActive;      /* dctran.c: follow the recorded step/order sequence */
extern int PSSlastStepOrder;     /* dctran.c: integration order of the step just accepted */
void PSSreplayStep(CKTcircuit *ckt);   /* called once per fresh step, before the corrector */

void PSSbreakInstallHook(CKTcircuit *ckt);
extern void (*PSSxsPostInitHook)(CKTcircuit *ckt);   /* card #551: dctran.c calls it after the t=0 pass of a restarted run (XSPICE state put back) */
extern void (*PSSxsPostInitHook)(CKTcircuit *ckt);   /* card #551: dctran.c calls it after the t=0 pass of a restarted run (XSPICE state put back) */
extern int PSSeventLocateActive, PSSeventLocateRestart, PSSswLagActive, PSSeventLocateCrossed;   /* card #498 2nd pass: exact sw event location */
int PSSeventLocate(CKTcircuit *ckt, double *newdelta);
int PSSshootLagEnabled(void);   /* card #498: PAC/PSTB recording passes follow the PSS orbit's sw treatment */  /* called right after the restart-path breaks reset */

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
extern int PSShidActive;           /* card #546 */
extern int PSSxsTdStateful(const char *nm);   /* card #557 */
extern int PSSxsPacBegin(CKTcircuit *ckt);   /* card #557: put the carried XSPICE state back before a td recording run */
extern void PSSxsPacEnd(void);
extern int PSSxspicePolicy(CKTcircuit *ckt, int kind, const char *who);   /* card #551: XSPICE model whitelist (0 ok, E_UNSUPP refused) */
void PSShidRecordStep(CKTcircuit *ckt);

void PSSrippleRecordStep(CKTcircuit *ckt);   /* called right after CKTaccept(), same site as the hooks above */

/* ---- card #489 (plan T6.2): free-running (autonomous) shooting's
 * Poincare-section trigger detector. Same no-op-unless-active pattern as
 * the hooks above. See pssshoot.c's PSStrigRecordStep() header comment
 * for the crossing-detection/bracketing it does. ---- */

extern int PSStrigActive;   /* dctran.c: watch the trigger node for a crossing */

void PSStrigRecordStep(CKTcircuit *ckt);   /* called right after CKTaccept(), same site as the hooks above */

/* ---- card #494 (plan T6.1): saltation (event) correction for ideal `sw`
 * voltage-controlled switches, shared by pssshoot.c's PSSsensRecordStep()
 * (real monodromy Phi) and pacshoot.c's PSSpacRecordStep() (complex forced
 * response). See pssshoot.c's PSSsaltFindEvents() header comment for the
 * derivation (Fang's saltation matrix, delta_x+ = delta_x- - (f+-f-)*dt_e,
 * dt_e = -(grad(h).delta_x-)/h_dot).
 *
 * THIRD pass (card #494) rewrite: the first two passes routed the
 * correction through an ag0-scaled companion-model CURRENT injected into
 * the event step's own linear solve -- step-size dependent by
 * construction (ag0 ~ 1/h_step), and the leading suspect for the 2.2x
 * error both prior passes measured and could not resolve. This pass
 * instead precomputes (f+-f-) at the switch's own power nodes directly
 * from the LOCAL capacitance physically tied to each node (dfdt_pos/neg
 * below, in state-units/second, no ag0, no dependence on the step size
 * the event happens to land on -- see PSSsaltFindEvents()'s own comment),
 * and the caller applies it as a direct VALUE jump to the propagated
 * sensitivity/phasor state at the event (delta_x(t_e-) -> delta_x(t_e+)),
 * not as an RHS current routed through that step's own G+ag0*C solve.
 *
 * One event per SW instance that changed state
 * (CKTstate0[SWswitchstate] != CKTstate1[SWswitchstate]) at the CURRENT
 * accepted step. ---- */

typedef struct {
    int posNode, negNode;       /* switch power terminals */
    int posCntrl, negCntrl;     /* switch control terminals (h = v(posCntrl)-v(negCntrl)-vt) */
    double v_switch;            /* V(posNode)-V(negNode) at the accepted point */
    double dG;                  /* conductance jump, g_after - g_before */
    double h_dot;                /* d(v_ctrl)/dt across this step (finite difference,
                                   * denominator = the step's TRUE elapsed accepted time,
                                   * card #494 fifth pass -- see PSSsaltFindEvents()) */
    /* card #510: posRow/negRow is the STATE ROW the dfdt_pos/dfdt_neg jump
     * is actually applied to -- posNode/negNode themselves (the THIRD
     * pass's original capacitor case, unchanged) when a capacitor sits
     * directly on that node, OR an adjacent inductor's own branch-current
     * row (INDbrEq, card #510's new "input-forcing"/inductor-routed case,
     * see pss_salt_ind_jump() in pssshoot.c) when there is no capacitor
     * there but there IS a unique adjacent inductor, OR 0 (skip, no jump
     * applied -- an exogenous/purely-algebraic node, correctly a no-op)
     * when neither applies. Always use posRow/negRow, never posNode/
     * negNode directly, when applying the correction -- see
     * PSSsaltFindEvents()'s header comment for why a non-state node's own
     * row is inert (overwritten from scratch next step) even when
     * nonzero. */
    int posRow, negRow;
    double dfdt_pos;             /* (f+-f-) for posRow: capacitor case
                                   * -dG*v_switch/C_local(posNode); inductor
                                   * case, see pss_salt_ind_jump(); 0 if
                                   * posRow==0 */
    double dfdt_neg;             /* (f+-f-) for negRow, mirror of dfdt_pos */
    const char *name;           /* instance name, for diagnostics */
    /* card #559: posQ/negQ = 1 when that power node carries DEVICE-MODEL charge (diode cj/tt, MOS junction/overlap, ...) in the
     * generic td-PAC charge matrix; the jump of that node is then formed on the charge vector (dqdt_pos/neg = the charge-rate
     * jump of the switch current, -dI / +dI, no capacitance in it) by pac_step_generic(), not divided by a native capacitance.
     * dfdt_pos/neg of such a node stays the rate over the node's total diagonal charge (used for the settle estimate only). */
    int posQ, negQ;
    double dqdt_pos, dqdt_neg;
    /* card #575: posF/negF = 1 when that power node carries device-model charge and NO native capacitor, and the jump is routed to
     * the one inductor at the node (posRow/negRow = its branch row): the node is treated as a fast (slaved) variable, see
     * PSSsaltCont in pssshoot.c.  pac_gen_salt_guard() accepts the model charge at such a node. */
    int posF, negF;
} PSSsaltEvent;

/* card #559: set by pacshoot.c while a generic-charge td-PAC pass runs: returns the total diagonal |C| (Farads) of the node when its
 * loaded charge row carries more than the native capacitors/inductors on it (device-model charge), else 0. */
extern double (*PSSsaltModelCapHook)(CKTcircuit *ckt, int node);
/* card #575: the smallest device-model part of that node's diagonal charge on the orbit (the reverse-biased junction, typically) */
extern double (*PSSsaltModelCapMinHook)(CKTcircuit *ckt, int node);

int PSSsaltContActive(void);   /* card #579 */
extern int PSSsaltNoRoutedWindow;   /* card #583 */
#define PSS_SALT_MAX_EVENTS 8

extern int PSSsaltActive;   /* env-var gated (PSS_SALTATION, default on); see pssshoot.c */

/* Fills ev[0..return-1] (return <= PSS_SALT_MAX_EVENTS) with every sw/csw
 * event found at the circuit's CURRENT accepted step (ckt->CKTtime,
 * ckt->CKTdelta as the just-taken step size). Returns 0 if none (the
 * overwhelmingly common case -- every step of every deck with no ideal
 * switch, or whose switches didn't change state this step). Cheap unless
 * PSSsaltActive. Returns -1 (card #510) if a real event was found whose
 * forcing term this engine cannot compute safely (see pss_salt_ind_jump()
 * in pssshoot.c) -- a clear stderr message is printed before returning;
 * the caller (PSSsensRecordStep()/PSSpacRecordStep()) must treat this as
 * a hard failure (set its own g_sens_failed/g_pac_failed and return) --
 * REFUSE, not a silent wrong answer, per card #510's brief. */
int PSSsaltFindEvents(CKTcircuit *ckt, PSSsaltEvent *ev, int max);

/* card #494 third pass: call once at the START of every fresh
 * period/recording run (pss_sens_reset()/pac_reset(), both files) so
 * PSSsaltFindEvents() can filter the spurious ~1e-10s "event" both prior
 * passes found and did not fix -- CKTstate1 on the FIRST accepted step
 * after a shooting-engine period restart does not carry the converged
 * orbit's own wrap-around history, which fabricated a huge, bogus h_dot
 * every restart. The next PSSsaltFindEvents() call after PSSsaltReset()
 * unconditionally returns 0 (real duty-cycle edges land later in the
 * period; this costs one accepted step's worth of non-detection per
 * restart, the "ignore crossings within the first step" filter the
 * card's brief asks for). Also clears the fifth pass's own h_dot
 * elapsed-time tracker and post-event L-stable settle window (see
 * PSSsaltFindEvents()'s header comment) so neither carries state across a
 * period restart. */
void PSSsaltReset(void);

/* card #494 FIFTH pass (plan T6.1): the stiffness half of the fix (task
 * 2 in the card's brief). Node `sw` (and any node like it -- a small
 * capacitance directly against a large switch/resistor conductance) has
 * a local RC time constant many orders of magnitude below the accepted
 * step sizes the surrounding TRAPEZOIDAL recursion normally uses; TRAP is
 * A-stable but not L-stable, so it rings (amplification -> -1, not 0)
 * instead of damping the huge post-jump sensitivity spike PSSsaltFindEvents()
 * injects there, leaking gain into the slow states (see docs/pss-pac.md's
 * "FOURTH pass" section). PSSsaltFindEvents() estimates that node's own
 * tau = C_local/G_local (pss_salt_local_cap()/pss_salt_local_cond()) and
 * opens a settle window (ckt->CKTtime .. event_time + PSS_SALT_SETTLE_TAUS
 * * tau) every time it reports a real event; PSSsaltForceOrder1() below
 * is the query PSSsensRecordStep()/PSSpacRecordStep() call, right after
 * their own ag1-based order inference, to OVERRIDE that order to 1
 * (backward Euler, L-stable) for every step inside the window -- exactly
 * the same "force order 1 after a breakpoint" idiom dctran.c's own XSPICE
 * path already uses for the real trajectory (see dctran.c's `Cut
 * integration order if first timepoint after breakpoint`), applied here
 * to this file's SEPARATE sensitivity/phasor recursion instead, since the
 * real trajectory's own breakpoint/order handling has no knowledge of the
 * (much larger, saltation-injected) spike living only in that recursion.
 * A single call, one integer compare; returns 0 (no override) whenever no
 * event has been seen this run or the window has already closed. */
int PSSsaltForceOrder1(CKTcircuit *ckt);

/* card #494 SIXTH pass (plan T6.1): VALIDATED-BAND GUARD. Returns the
 * number of real saltation events PSSsaltFindEvents() has reported (n>0
 * calls, summed) since the last PSSsaltReset() -- i.e. whether THIS run's
 * circuit actually has an ideal-switch event a perturbation reaches at
 * all (0 for every existing fixture with no `sw`/`csw` event on the
 * perturbed path -- exact no-op for those, same as PSSsaltForceOrder1()).
 * Used by pacshoot.c's per-frequency-point loop to gate the near-fs/2
 * unvalidated-margin warning to saltation-affected decks only -- see
 * docs/pss-pac.md's "SIXTH pass" section for the cross-validated band
 * this threshold is drawn from (0.5dB/3deg to >=99.5% of fs/2 on
 * pac_r6o_ctrl2out.cir against tools/pac_perturbed.py, an independent
 * nonlinear reference sharing no code with this engine). */
int PSSsaltEventsSeen(void);

/* card #555: PSS_RUN_STATS diagnostics */
void PSSrunStatsAdd(const char *name, double dt);
double PSSrunStatsNow(void);

#endif
