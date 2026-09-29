/**********
card #465 (plan T5.1-T5.2): time-domain shooting-based periodic AC (PAC,
sideband 0) and periodic loop gain (.pstb) for switching converters, built
on the shooting-Newton PSS in pssshoot.c/psseigen.c. See pacshoot.c's
header comment for the full design (equations, per-step recursion, the
boundary-value solve, demodulation) and research/ngspice-internals-pss-pac.md
sec 5.3 / research/simplis-pop-ac.md sec 6-7 / research/cadence-stb-pstb.md
sec 2.4/4 for the theory this implements.

Entry point: pac_pstb_run(), called from pssshoot.c's PSSshootDriven() right
after the shooting Newton has CONVERGED (x0 = the converged reduced state,
same convention as pssshoot.c's own PSSstateMap: one entry per native
capacitor terminal voltage, one per native inductor branch current, in
device-list order -- pacshoot.c rebuilds this SAME map independently by
walking the SAME circuit device lists in the SAME order, rather than
depending on pssshoot.c's internal (file-static) PSSstateMap type, to avoid
changing that file's own visibility/ABI for this card).

PSSpacActive / PSSpacRecordStep(): the dctran.c hook, exactly the same
no-op-unless-active pattern as pssshoot.h's existing PSSbreakRecordActive/
PSSsensActive/PSSrippleRecordActive hooks (all called from the SAME site,
right after CKTaccept()). Every other analysis is unaffected.
**********/

#ifndef ngspice_PACSHOOT_H
#define ngspice_PACSHOOT_H

#include "ngspice/cktdefs.h"
#include "ngspice/pssdefs.h"

extern int PSSpacActive;
void PSSpacRecordStep(CKTcircuit *ckt);

/* Called once, after PSSshootDriven()'s Newton has converged. x0[0..n) is
 * the converged reduced state (cap voltage / inductor current, device-list
 * order), T is the period. Runs `.pac ... solver=td` and/or `.pstb`
 * post-processing according to job's PSSdoPAC/PACsolver/PSSdoPSTB flags
 * (a no-op, returning OK immediately, if neither applies). */
int pac_pstb_run(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T);

/* card #489 SECOND pass (plan T6.2, task 5): the general form, also
 * usable on a free-running (PSSshootFreerun()) orbit. fT[0..n) is the
 * vector field f(T) = dx/dt at the trigger crossing on the converged
 * orbit (device-list order, SAME as x0/PACstateMap), and trig_node_num is
 * the raw circuit node number of the `trigger=` node (job->PSStrigNode's
 * own ->number) -- both NULL/-1 for the ordinary driven case, in which
 * case this is EXACTLY pac_pstb_run()'s old (unmodified) behaviour: no
 * projection correction applied. When fT/trig_node_num ARE given, the
 * periodicity boundary condition (Phi, z) pac_engine_run() solves is
 * corrected by the terminal projection `Proj = I - f(T) e_trig^T /
 * (e_trig . f(T))` before the dense (e^{jwT}I-Phi)xi=z(w) solve (see
 * docs/pss-pac.md's free-running section for the derivation); this needs
 * trig_node_num to match a native capacitor's positive terminal among the
 * reduced states pacshoot.c's own map tracks -- if it does not (trigger
 * on an inductor branch, or an algebraic node not directly a reduced
 * state), this refuses cleanly (returns nonzero, prints why) rather than
 * silently running the uncorrected (wrong, for a free-running orbit)
 * fixed-period engine. pac_pstb_run() is a thin wrapper calling this with
 * fT=NULL, trig_node_num=-1. */
int pac_pstb_run_ex(CKTcircuit *ckt, PSSan *job, const double *x0, int n, double T,
                     const double *fT, int trig_node_num);

#endif
