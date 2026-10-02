/**********
Author: 2010-05 Stefano Perticaroli ``spertica''
Review: 2012-10 Francesco Lannutti
**********/

#ifndef ngspice_PSSDEFS_H
#define ngspice_PSSDEFS_H

#include "ngspice/jobdefs.h"
#include "ngspice/tskdefs.h"
    /*
     * PSSdefs.h - defs for pss analyses
     */

typedef struct {
    int JOBtype;
    JOB *JOBnextJob;
    char *JOBname;
    double PSSguessedFreq;
    CKTnode *PSSoscNode;
    double PSSstabTime;
    long PSSmode;
    long int PSSpoints;
    int PSSharms;
    runDesc *PSSplot_td;
    runDesc *PSSplot_fd;
    int sc_iter;
    double steady_coeff;

    /* card #458, plan T4.1-T4.3: shooting-Newton PSS for driven (clocked)
     * switching converters, as an alternative to the brute-force
     * transient-to-steady-state loop above. 0 = existing algorithm
     * (default, unchanged); 1 = shooting-Newton (`.pss ... method=shoot`).
     * Only valid for driven circuits (Enhancement-176 detection); the
     * autonomous/free-running case is explicitly out of scope (design
     * ruling #447 A) and is rejected with an error at run time. */
    int PSSmethod;

    /* Enhancement-119: the converged periodic operating point, retained past the
     * analysis as the substrate the periodic small-signal analyses (PAC / pnoise
     * / PXF) linearize around. PSS already samples the node voltages over one
     * period for its DFT but frees them; here they -- and the device states,
     * which the reactive Jacobian C(t) needs -- are kept on the job. Row-major
     * per sample: [unknown + sample*PSSopMsize] and [state + sample*PSSopNumStates]. */
    long   PSSopPoints;      /* number of time samples over one period */
    int    PSSopMsize;       /* matrix size (nodes + branch currents) */
    int    PSSopNumStates;   /* CKTnumStates captured per sample */
    double PSSopFreq;        /* converged fundamental frequency (Hz) */
    double *PSSopTimes;      /* [PSSopPoints] sample times across the period */
    double *PSSopVoltages;   /* [PSSopMsize * PSSopPoints] node voltages per sample */
    double *PSSopStates;     /* [PSSopNumStates * PSSopPoints] device states per sample */

    /* Enhancement-122: PAC (.pac) frequency sweep. The .pac card runs PSS, then
     * sweeps a small-signal input frequency and, at each point, solves the
     * harmonic conversion matrix (E-121) built from the retained operating point,
     * emitting the sideband-0 node responses as a complex plot vs frequency. */
    int    PSSdoPAC;         /* 1 if this job also runs a PAC sweep */
    double PACfStart;        /* input-frequency sweep start (Hz) */
    double PACfStop;         /* input-frequency sweep stop (Hz) */
    int    PACpoints;        /* points per decade/octave, or total for linear */
    int    PACstepType;      /* 0 = LIN, 1 = DEC, 2 = OCT */
    int    PACmaxSideband;   /* number of output conversion sidebands each side (0 = sb0 only) */

    /* card #465 (plan T5.1-T5.2): `.pac ... method=shoot solver=td` selects
     * the time-domain shooting-based PAC engine (pacshoot.c) in place of
     * the frequency-domain harmonic-conversion-matrix engine above, which
     * stays completely unchanged and is still the default. Only usable
     * together with PSSmethod==PSS_METHOD_SHOOT (refused otherwise -- the
     * time-domain engine needs the shooting orbit's own per-step
     * companion/Jacobian structure, which the driven brute-force PSS
     * algorithm does not produce). See pacshoot.c's header comment for
     * the design. */
    int    PACsolver;        /* 0 = frequency-domain (default), 1 = time-domain (card #465) */

    /* card #465: `.pstb` periodic loop gain (Tian two-injection on the
     * LPTV orbit, sideband 0), built on the SAME time-domain PAC engine.
     * A `.pstb` card creates its own PSS-type job (like `.pac`) with
     * PSSdoPSTB=1 and a probe name; it always needs method=shoot (refused
     * otherwise). Reuses the PAC sweep fields above (PACfStart/Stop/
     * points/stepType) for its own frequency sweep -- a `.pstb` job never
     * also sets PSSdoPAC, so there is no field clash. */
    int    PSSdoPSTB;        /* 1 if this job runs a .pstb loop-gain sweep */
    IFuid  PSTBprobeName;    /* name of the probe V-source instance (probe=) */

    /* Enhancement-124: periodic noise (.pnoise). Runs PSS, then folds each device's
     * noise through the conversion-matrix adjoint over all sidebands to get the
     * output noise spectrum. Reuses the PAC sweep fields for the frequency sweep. */
    int      PSSdoPnoise;    /* 1 if this job runs a pnoise sweep */
    CKTnode *PnOutNode;      /* pnoise output node (reference = ground) */
    IFuid    PnInSrc;        /* input source name, for the input-referred spectrum */
    int      PSSpnCyclo;     /* Enhancement-126: 1 = cyclostationary noise (per-sample
                              * bias + time-domain transfer, averaged over the period);
                              * 0 = stationary (noise PSD at one operating point) */

    /* card #509: `.pnoise ... method="shoot" solver="td" [dev=1]` -- the
     * shooting-PSS-orbit true-adjoint noise engine (pacshoot.c's
     * pac_td_pnoise()/pac_engine_run_adjoint_noise()), in place of the
     * frequency-domain harmonic-conversion-matrix engine above (unchanged,
     * still the default -- same additive PACsolver-gated pattern as
     * .pac/.pstb/.pxf solver=td). `dev=1` additionally emits a per-device
     * noise-generator contribution table (onoise_<inst><suffix> vectors,
     * same naming as plain .noise's own summary report) alongside
     * onoise_spectrum; 0/unset (default) emits onoise_spectrum only. */
    int      PSSpnDevTable;

    /* Enhancement-125: periodic transfer function (.pxf). The ADJOINT counterpart of
     * PAC: solve Hᵀ Ψ = e_{out,0} and dot Ψ with the netlist AC-source pattern to get
     * the transfer from the input to a fixed output at each sideband. Reuses the PAC
     * sweep fields + PACmaxSideband. */
    int      PSSdoPXF;       /* 1 if this job runs a PXF sweep */
    CKTnode *PxOutNode;      /* PXF output node (reference = ground) */

    /* card #499 (plan T5.4): time-domain shooting-based `.pxf ... method=
     * "shoot" solver="td"` -- the shooting-PAC engine's own per-source
     * periodic transfer function, in place of the frequency-domain
     * harmonic-conversion-matrix engine above (which stays completely
     * unchanged and is still the default -- same additive pattern as
     * PACsolver/PSTB on top of the existing .pac/.pnoise). One transfer
     * function per requested independent V/I source, to the SAME
     * PxOutNode above, emitted as "pxf_<source>" vectors in one plot --
     * see pacshoot.c's pac_td_pxf(). `sources=` is repeatable (same
     * reason as STBan's `sens=` -- stbdefs.h/docs/sens.md -- a comma list
     * risks the tokenizer comma-splitting hazard card #457 already hit);
     * `sources="all"` (the default, also case-insensitive, ALSO what an
     * omitted sources= means) short-circuits any prior/later name and
     * sweeps every independent V/I source instance in the circuit. */
    int      PxfSourcesAll;   /* 1 = sources="all" (default when unset) */
    IFuid   *PxfSourceNames;  /* [PxfSourceCount], repeated sources="Vname" */
    int      PxfSourceCount;

    /* card #504: `.pxf ... adjoint=1|0` -- selects pac_td_pxf_adjoint()
     * (src/spicelib/analysis/pacshoot.c) in place of the card #499
     * batched-forward-per-source engine above (pac_td_pxf()'s own loop
     * over pac_engine_run(), still the default/fallback -- see that
     * function's "SHIPPED DESIGN" comment, now superseded for the common
     * case by pac_td_pxf_adjoint()'s true backward/transposed sweep,
     * amortising every requested source behind ONE extra period
     * traversal instead of K). PxfAdjointSet distinguishes "not given on
     * this card" (0, auto: pac_td_pxf() picks the engine itself from the
     * requested source count, PXF_ADJOINT_AUTO_THRESHOLD in pacshoot.c)
     * from an explicit `adjoint=0` (force the v1 forward engine, e.g. for
     * A/B comparison) -- PxfAdjoint alone can't carry that distinction
     * since 0 is also PxfAdjoint's own zero-initialised default. */
    int      PxfAdjoint;      /* explicit value, only meaningful if PxfAdjointSet */
    int      PxfAdjointSet;   /* 1 if `adjoint=` appeared on this card at all */

    /* Enhancement-132: periodic S-parameters (.psp). Runs PSS, then for each RF
     * port injects a unit current in sideband 0 through the conversion matrix and
     * reads the per-sideband port waves, building the periodic scattering matrix
     * S^(k) = B^(k)*A^-1 vs input frequency. Reuses the PAC sweep fields
     * (PACfStart/Stop/points/stepType/maxSideband). For a time-invariant circuit
     * the sideband-0 result reduces to the ordinary .sp S-matrix. */
    int      PSSdoPSP;       /* 1 if this job runs a periodic S-parameter sweep */

    /* card #465 THIRD pass, task 2: `pacpoints=N` on `.pac`/`.pstb` -- a card
     * parameter for pacshoot.c's own recording-pass time grid (pac_run_period()
     * OVERRIDES ckt->CKTmaxStep/CKTstep/CKTdelmin to T/pacpoints for just its
     * two recording-period reruns), replacing the PAC_RECORD_POINTS env var
     * (second pass) as the primary way to ask for it. 0 = unset (default):
     * pac_run_period() then falls back to the deprecated PAC_RECORD_POINTS
     * env var if set, else to job->PSSpoints (the same grid PSS's own
     * convergence phase already uses by default -- so an old deck with
     * neither the card parameter nor the env var sees byte-identical
     * behaviour to before this field existed). See pacshoot.c's
     * pac_run_period() for the resolution order. */
    long int PACrecordPoints;

    /* card #489 (plan T6.2): free-running (autonomous) shooting. A
     * `.pss ... trigger=<node> [trigval=<v>] [trigdir=rise|fall]
     * method=shoot` card selects the Poincare-section engine
     * (PSSshootFreerun(), pssshoot.c) instead of the fixed-period driven
     * engine (PSSshootDriven()) -- selected whenever PSStrigNode is
     * non-NULL (trigger= was given), independent of whether a driven
     * source is also detected (a trigger takes precedence: the deck asked
     * for event-located shooting explicitly). fguess is then used only as
     * a rough period ESTIMATE (to size MAX_PERIOD and StabTime's own
     * pre-run), not an exact clock -- unlike the driven engine, where it
     * is the exact clock. See docs/pss-pac.md's free-running section for
     * the algorithm (event-located Poincare map, no explicit saltation --
     * design ruling #447 A restricts this engine to smooth switch models,
     * same as the driven engine, so every "event" here is a smooth
     * threshold crossing, not a true discontinuity). */
    CKTnode *PSStrigNode;    /* NULL = driven/fixed-period (default, unchanged) */
    double   PSStrigVal;     /* threshold the trigger node must cross (default 0.5) */
    int      PSStrigDir;     /* +1 = rise (default), -1 = fall */

    /* card #494 (plan T6.1, saltation): `driven=1` on a `.pss ...
     * method=shoot` card. Enhancement-176's own driven-mode detection
     * (dcpss.c) only recognises a native V/I source with a time function
     * (VSRCfuncTGiven -- a PULSE/SIN/... on a plain V or I source). An
     * XSPICE code-model clock (e.g. `d_pwm`, refs/r6_vm_buck_ideal.cir) is
     * an opaque MIF instance outside that scan -- a genuinely driven,
     * fixed-frequency circuit is otherwise refused with "needs a driven
     * (clocked) circuit". This is an explicit, opt-in admission that the
     * deck IS driven at the given fguess (the same honesty stance
     * `trigger=` already uses for the free-running case): it does not try
     * to introspect XSPICE code models, and does nothing for a genuinely
     * autonomous circuit unless the user asserts driven=1 themselves. */
    int      PSSforceDriven;

    /* card #537: `krylov=1` on the .pss/.pac/.pstb/.pxf/.pnoise card (or the
     * environment variable NGSPICE_PSS_KRYLOV=1): matrix-free Newton-GMRES
     * shooting for driven method=shoot -- the monodromy is never assembled,
     * each Krylov vector costs one (centred) perturbed period run. Default 0
     * (dense finite-difference monodromy, byte-for-byte the previous
     * behaviour). See pssshoot.c "card #537". */
    int      PSSkrylov;
} PSSan;

enum {
    PSS_METHOD_DRIVEN = 0,   /* existing brute-force periodicity-test algorithm */
    PSS_METHOD_SHOOT  = 1    /* card #458: shooting-Newton */
};

enum {
    PAC_SOLVER_FD = 0,       /* existing frequency-domain harmonic conversion matrix */
    PAC_SOLVER_TD = 1        /* card #465: time-domain shooting-based PAC/PSTB */
};

enum {
    GUESSED_FREQ = 1,
    STAB_TIME,
    OSC_NODE,
    PSS_POINTS,
    PSS_HARMS,
    PSS_UIC,
    SC_ITER,
    STEADY_COEFF,
    PSS_METHOD,
    PAC_DOPAC,
    PAC_FSTART,
    PAC_FSTOP,
    PAC_POINTS,
    PAC_STEPTYPE,
    PAC_MAXSB,
    PNOISE_DO,
    PNOISE_OUT,
    PNOISE_INSRC,
    PXF_DO,
    PXF_OUT,
    PNOISE_CYCLO,
    PSP_DOPSP,
    PAC_SOLVER,
    PSTB_DO,
    PSTB_PROBE,
    PAC_RECPOINTS,
    TRIG_NODE,
    TRIG_VAL,
    TRIG_DIR,
    PSS_FORCE_DRIVEN,
    PXF_SOURCES,
    PXF_ADJOINT,
    PNOISE_DEVTBL,
    PSS_KRYLOV,
};

#endif
