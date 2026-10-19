# Changes in this fork

Every commit on `main` after the upstream base `aaba50ba3f49` (branch `upstream-base`), grouped by the first release tag that contains it, oldest first. Numbers such as `#458` or `card #465` refer to an internal tracker; each message describes its change without it. Commits marked *(upstream)* are ngspice development fixes picked unchanged, with their original authors.

For an overview, see [README.md](README.md).

## s125 — Native .stb, shooting .pss, periodic small-signal (.pac/.pstb/.pxf/.pnoise, solver="td"), ideal-switch saltation, XSPICE hidden state, checkpoint, speed rounds, stock fixes

126 commits.

- `931abf7a5` stb: native Tian two-injection loop-stability analysis, .stb card and command
- `cdb17aa51` stb: fix plot-type-abbreviation collision with "op"
- `d79fdfa9d` stb: document the probe= quoting requirement on the .stb netlist card
- `4744f9bfe` stb: wire localgnd= into the current-injection RHS
- `8d44ec9cf` stb: probe=<inst>:<terminal> auto-insert, and probe2=/mode=dm\|cm parsing
- `fade86856` stb: native probe=<Vp>,<Vn> mode=dm\|cm two-injection math (stban.c)
- `4489c43a2` stb: device-based (gain-nulling) mode, probe=<mos1 instance>
- `fc43a3aa0` stb: fix two device-mode bugs found in validation (ground trash slot, gm read before op point)
- `3ff3d0eb1` stb: fix loop-mode ground-terminal trash-slot bug, warn on degenerate probe
- `6bbb3d7e7` stb: device-based (gain-nulling) mode for BSIM4, probe=<BSIM4 instance>
- `a7334c36d` stb: warn on positive feedback at low frequency; label rising vs falling PM crossings
- `aad3906cd` stb: refuse invalid n/fstart/fstop instead of silently defaulting; netlist-card auto-quote for probe=/probe2=/mode=/localgnd=
- `2e4f0a534` stb: narrow the low_f_positive_feedback warning, add low_f_integrator and low_f_positive_feedback_low_gain labels
- `136d93344` frontend: inp_casefix must not step past the line terminator after a quote (card #477)
- `c01830518` accel: hard-coded schedule(guided) for the BSIM4 OMP load loop (card #478 T9.2, candidate 0103)
- `2153d8840` numparam: explicit draw-once form for random .param expressions (card #482)
- `6822ccd04` pss: shooting-Newton method=shoot for driven PSS (card #458 T4.1-T4.3, WIP)
- `f0bd3c7e5` pss: shooting-Newton second pass -- eigensolver + forced-grid monodromy (card #458 T4.1-T4.3)
- `0ac10b729` pss: shooting-Newton third pass -- staged warm-up, Broyden, degenerate-fixedpoint guard (card #458 T4.1-T4.3)
- `f2a7083cd` pss: analytic (discretise-then-linearise) sensitivity monodromy, opt-in PSS_SHOOT_JACOBIAN=analytic (card #458 fourth pass, plan T4.2)
- `818007cb8` pss: fifth pass -- ripple-relative Newton tolerance, automatic fast-state reduction, analytic Jacobian on the reduced set (card #458 T4.1-T4.3)
- `a446ece59` pss: sixth pass -- long-settling speed >=10x, Points-robustness fix (card #458 T4.1-T4.3)
- `8cf9e3775` pac: time-domain shooting-based PAC (.pac solver=td) and periodic loop gain (.pstb), built on the shooting-Newton PSS (card #465, plan T5.1-T5.2)
- `19a71fddc` pac: PAC_RECORD_POINTS decouples PAC's recording-pass grid from PSS's own convergence-phase grid (card #465 SECOND pass, task 3 "speed")
- `908804f32` pac: keep PSTB Margins' low-frequency classification consistent with stban.c patch 0013 (card #465 THIRD pass, rebase)
- `5dca0d2a0` pac: exact-Jacobian refactor for the recording pass, opt-in PAC_EXACT_JACOBIAN (card #465 THIRD pass, task 1)
- `b3f4df421` pac: pacpoints=N card parameter on .pac/.pstb, PAC_RECORD_POINTS deprecated (card #465 THIRD pass, task 2)
- `208893727` pac: fix pstb/pac demod O(1/Points) discretisation error; red-proof toggles (card #465 FOURTH pass, tasks 1 and 4)
- `1975d76b4` pss: force the damped-step Newton trial onto the accepted grid, fixing R4 D~=0.45 non-convergence (card #465 FIFTH pass, task 2)
- `e2e1c38a5` accel: B-source cache moved OUTSIDE INPparseNode (card #483 T9.3, candidate 0104)
- `354efd8ba` accel: hard-coded schedule(guided) for the other BSIM3/BSIM3v32/HiSIM2/ BSIM-SOI/BSIM4v5-7 OMP load loops (card #483 T9.3 task 2, candidate 0105)
- `2116c7eae` checkpoint: keep the current breakpoint across loadstate, not just future ones
- `53b7ce898` checkpoint: dedupe near-equal breakpoints after keeping the current one
- `5e30359c1` stb: adjoint sensitivity (sens=all\|<name>) for loop-mode dT/dp, dPM/dp, dfc/dp, dGM/dp
- `de40777dc` stb: fix sgen early-break corrupting instance chains across sens= frequency points
- `bfc2cbe44` stb: repair DEVsetup's persistent per-instance matrix-element pointers after each sens= parameter
- `549d57ecc` stb: fix dGM/dp missing the gain-crossing frequency's own implicit shift with p
- `cbaea0dca` pss: fix q>1 method=shoot trajectory-EMISSION corruption (card #488)
- `f90cb58cb` pss: free-running (autonomous) shooting -- .pss ... trigger=<node> (card #489)
- `ca6a9561f` pac: free-running (autonomous) .pac/.pstb -- period-variation projection (card #489 SECOND pass, plan T6.2, task 5)
- `48d6dea3f` frontend: com_rset must replay ft_curckt->ci_vars across reset (card #492)
- `a85c1ac22` pss: driven=1 override for XSPICE-clocked decks (card #494, plan T6.1)
- `532bf7a2b` pss/pac: saltation (event) correction for ideal `sw` switches (card #494, plan T6.1)
- `13495866b` pss/pac: saltation THIRD pass -- state-space jump, no ag0 in the jump (card #494, plan T6.1)
- `e03a51dd9` pss/pac: saltation FIFTH pass -- exact h_dot denominator + L-stable post-event settle window (card #494, plan T6.1)
- `4dc74b3e3` pss/pac: saltation SIXTH pass -- validated-band guard near fs/2 (card #494, plan T6.1)
- `a86a5f234` pxf: time-domain shooting-based periodic transfer function (card #499, plan T5.4)
- `c4ed719f6` dcpss: frequency-domain .pss/.pxf read past the harmonics buffer with harmonics=1 (card #502)
- `a679b6fb2` accel round 7: opt-in autoklu (equation-count-threshold automatic KLU default), card #496
- `fdc56cc31` accel round 7: fix autoklu -- decide before DEVsetup/NIinit, not after
- `bc32b34bf` accel round 7: fix autoklu build error -- CKTlastNode is a CKTnode*, use ->number
- `283fde5f7` cktsopt: klu_memgrow_factor stored a comparison (0/1) instead of the value (card #501)
- `68fdf22d4` pss/pac: free-running .pac/.pstb -- f(T) forward-diff + validated-band guard (card #493)
- `1fe3152f9` stb: `.stb ... instant` -- linearise at a transient time point, no CKTop() re-solve (card #500)
- `f19fc66ac` pxf: true adjoint/transposed engine, .pxf ... adjoint=1 (card #504)
- `9964bf430` pssshoot.c: PSSshootFreerun() -- warn on a near-degenerate trivial Floquet eigenvalue (card #503)
- `53674ba71` pssshoot.c: PSSshootFreerun() -- gated PSS_FREERUN_FD_EPS diagnostic (card #503 SECOND pass)
- `bc0c05664` pssshoot.c: PSSshootFreerun() -- polish x0/T for .pac/.pstb hand-off (card #503 THIRD pass, fixes 20-100kHz)
- `8ef483d85` pssshoot.c: PSSshootFreerun() -- ungate .pac/.pstb inside the validated band (card #503 FOURTH pass)
- `e1174bb20` pssshoot/pacshoot: saltation input-forcing term for capacitor-free switch power nodes (card #510)
- `90b1008bd` pssshoot: warn once per analysis when saltation routes through an inductor (card #510 step 1)
- `785bc972a` pssshoot: re-derive the inductor-routed settle window from the accepted step, not L/R tau (card #510 step 2)
- `2177060e9` pssshoot: inductor-routed saltation warning states the validated accuracy after 0062 (0.5 dB/3 deg, one fixture) (card #510)
- `08e0eaceb` pss: guarded retry/fallback for ideal-switch period-run failures in the shooting-Newton engine (card #498)
- `b2c0f8893` pss: exact ideal-switch event location in shooting period runs (lagged sw decision + regula-falsi crossing time) and honest convergence gating (card #498, second pass)
- `f6360fab5` pss: exact sw event location/lag only for the FD-Jacobian path (card #498)
- `91dff8eba` pac/pstb: recording passes use the same exact sw event location as the PSS orbit (card #498)
- `463e41734` pss: drop the PSS_EVENT_DEBUG trace from PSSeventLocate (card #498)
- `12c051408` pss: with the lagged sw path, force the flip step's end time as a breakpoint; only the located crossing instant stays unforced (fixes q>1 shooting on R2 switched-RC, card #498 regression)
- `f6ca1e99f` pss: keep the step that starts at a located switch crossing short (3e-3*maxStep): trapezoid mis-weights the pre-flip capacitor current by h*dI/2 over it (R2 q>1 orbit error 3.5e-4 -> 5e-6 V) (card #498)
- `642888973` pss: with a sw in the circuit take the Newton baseline F from a forced rerun at x0, so F, FD columns and trials share one period map (R2 q=1 plateau 1.4e-6 after 6129a5187/8e201dd5f; recording run and forced runs differ by ~1e-7 V at the conduction kink) (card #498)
- `fb9460a0b` pnoise: true-adjoint periodic (cyclostationary) noise, `.pnoise ... method=shoot solver=td [dev=1]` (card #509)
- `d5af133c5` pssshoot: abandon doomed trust-clipped probe Newton attempts early (card #525)
- `6125662fc` pac/pstb free-running: self-consistent tight solver profile for polish+recording, polish at >=T/5000, projection off by default (card #517)
- `003c3038b` klu: refactor unscaled, pivots still chosen on the max-row-scaled matrix (card #520, round 13)
- `b9c2d59e0` pss: WIP exact step/order replay for forced period runs (PSS_REPLAY=1, opt-in for now) (card #505)
- `fed6b82e4` pss: one-sided FD rung, noise-robust fast-state detection, larger FD step for sub-ripple states, classic->replay final-stage fallback (card #505 WIP)
- `cd787de04` pss: replay ladder (classic -> replay from classic iterate -> replay from warm-up endpoint), replay-based orbit emission, PSS_SHOOT_VNTOL knob (card #505 WIP)
- `aaf6d2afd` pss: emit the orbit with the Newton run's own map (drop replay-only emission) (card #505 WIP)
- `ac114e562` pss: replay rung never accepts an unstable fixed point (keeps the D~0.55 period-2 red-proof refusing) (card #505)
- `37d40c73d` pss: final stage retries a failed first baseline run from an x0 advanced by a few free periods (R4d after #525's early probe exits; composition fix) (card #505)
- `ce845b3b8` klu: complex refactor unscaled too (SMPcLUfac), pivots from klu_z_factor still scaled (card #527, round 14)
- `1b8c050f4` res/cap: summed resistor stamps (real KLU tran) and R+C stamps (complex KLU AC) cached as a base and memcpy'd instead of re-scattered per load; invalidated by RES/CAP generation counters (card #527, round 14, prototype; NGSPICE_RESBASE=0 disables)
- `8fac51ee5` cktload/capload: capacitor base cache (Ax = Rbase + ag0*Cbase in one pass) and flat capacitor state/rhs load for real KLU tran; NGSPICE_CAPBASE=0 / NGSPICE_CAPFLAT=0 disable (card #532, round 15)
- `08b3a4c39` ifeval: flat B-source evaluator (value + derivative trees compiled once to a linear program in tree-walk order; ternary trees keep the tree walk); NGSPICE_FLATB=0 disables (card #532, round 15)
- `8af02bf4c` pss/pac saltation: ignore a flip whose control is not at the threshold; take v_switch before the flip (card #535)
- `861b5c1ac` pss/pac/pstb/pxf/pnoise: shooting engines run under .options klu (card #533)
- `53bf670e0` pss shoot: matrix-free Newton-GMRES shooting, opt-in krylov=1 / NGSPICE_PSS_KRYLOV=1 (card #537)
- `6009591b4` pss shoot (#538): automatic matrix-free shooting from 64 states, no Newton attempt on the first 130 warm-up periods, finer Arnoldi checkpoint
- `34f8bce88` pss freerun (#538): matrix-free Newton-GMRES on the Poincare return map for trigger= shooting (krylov=1 / automatic from 64 states), Arnoldi stability estimate
- `b0f168116` pss freerun (#540): correct Floquet verdict -- tight-profile stability pass with step-halving guard (UNKNOWN(not-computed) when it disagrees), dense Newton step 1e-4, free-running matrix-free from the first state (also with .pac/.pstb)
- `a6a273e38` pac td (#542): generic charge-based history for device-model charges (MOS/BSIM4/diode/BJT/mutual): C_k=dq/dx from two loads at ag0 and ag0/2 (third load at 2*ag0 = affinity guard, refuse naming the device type), state set = columns of C, exact transposed recursion for .pxf adjoint / .pnoise using stored J_k and C_k, periodic-continuation device ICs for the recording runs, native C/L path unchanged for native-only circuits (PAC_NATIVE_HISTORY=1 = old behaviour, PAC_FORCE_GENERIC=1 = cross-check)
- `f9483c196` pac td (#542): model-charge detection period (a chargeless diode beside ideal switches keeps the native path and saltation), prepass wall time in the generic message
- `64403b18e` pac td (#542): refuse BSIM3/BSIM4/SOI non-quasi-static models (internal charge-relaxation state is not a node charge; PAC_ALLOW_NQS=1 overrides)
- `16f5dc5eb` pac td (#542): model-charge detection period leaves no trace on the legacy flow (snapshot/restore of rhs, states, integration scalars): native circuits with chargeless diodes (R4d) bit-identical to s90 again
- `12abec545` #541: native td adjoint stores J_k and loads at x_k (INITFLOAT, not INITPRED); .pss accepts a PWL clock with r=
- `ed2938dbb` pnoise td (#544): device operating point at x_k for DEVnoise (stale CKTstate0 under a time-varying bias), BJT/JFET/MOS3/BSIM3 admitted; PAC_TRAJ_REPORT
- `5c4050dea` pac td (#544): ideal-switch (saltation) events with model charges on the orbit
- `723cc30c9` pss/pac free-running (#544): trigger= orbits with device-model charges (hidden-state crossing by interpolation, 1e-6 polish profile); .pac refusal lifted
- `dcce4fa12` pss shoot (#546): hidden-state continuation of the period restart on transistor-level circuits
- `64dec6ff4` pac td (#546): pacreltol= on .pac/.pstb/.pxf/.pnoise, RELTOL of the td recording passes only
- `7a7a8da7a` pss shoot/free-running (#549): hidden-state continuation of the trigger runs (R8), and of driven runs on diodes WITH charge
- `a71f3940d` pac td (#549): method=gear for the td recording pass (generic charge recursion with the charge history)
- `bb39583e6` pss shoot (#550): matrix-free replay rung before the dense fallback when the Krylov final stage plateaus near tol
- `08c3a660f` pss shoot (#550): matrix-free replay rung only for a stall near the target (residual <= 100 tol); the D=0.55 period-doubling deck keeps the dense path
- `aac779066` xspice (#551): PSS shoot carries XSPICE hidden state (analog code-model state, event system, model time) across period restarts; whitelist/refusal policy; td PAC B0 capture in AC analysis mode
- `9aaa585b5` xspice (#551): hidden-state periodicity check iterates the latched state to its fixed point; PSS_XS_NO_STATE red proof; validated-model table (Spice_Model_Name)
- `79dca476f` xspice (#551): periodicity check reports the baseline count, break at 1e-4
- `7d7facbc8` pss (#553): period runs keep no output plot - peak RSS of a matrix-free PSS on R7p 1.45 GB -> 40 MB
- `26c87c40a` pnoise (#556): folded flicker - each device's 1/f density is evaluated at the baseband frequency of every sideband, S_w(\|f - q fs\|), not at the output frequency f
- `5fd964376` accel round 16 (#555): R9 breakpoint table + noisyxspice lookup, td PAC pass 1 zero-imag solves and Phi reuse across .pstb injections, Floquet baseline, early Arnoldi stop, INDload typelook cache; PSS_RUN_STATS / PSS_KRYLOV_ARN_TRACE diagnostics
- `28fdc7a61` accel round 17 (#558): automatic KLU for shooting .pss, fork-parallel dense monodromy columns (opt-in), multi-RHS KLU solves in td PAC (opt-in)
- `4e742f37a` pss fork columns (#558): the free-running dense columns are not forked (R5 dense iteration not robust under fork); comment
- `1ec0df1ef` pnoise/pxf (#557): .pxf adjoint=1 and .pnoise on FREE-RUNNING (trigger=) orbits
- `201777039` xspice (#557): td .pac/.pstb/.pxf through stateful analog code models (s_xfer, int, d_dt; slew, hyst on a single branch)
- `b457d2913` pxf (#557): say so when adjoint=1 is not used because of the opt-in free-running projection
- `749cf531c` pac td (#559): ideal-switch saltation with device-model charge at a switch node that also has a native capacitor (charge-vector jump, state jump C dy = dQ over the coupled block, T from the smooth part); BSIM4 trnqsmod / BSIM3 nqsmod accepted in .pac/.pstb/.pxf (relaxation node is an ordinary state)
- `521060103` saltation (#568): the order-1 settle window after an ideal-switch event also ends after 2 accepted steps
- `9823ad2a8` saltation (#575): ideal-switch events at an inductor node without a dominant capacitor
- `232aff8f0` saltation (#579): boost and half-bridge events - diode-fed capacitor jump, output impulse, pre-flip voltage, telescoping continuation with hand-over
- `356410112` pac td (#581): say so when an inductor current is zero for part of the period (discontinuous conduction)
- `0cfd16790` saltation (#580): a capacitor node whose diode feeds another capacitor (boost) is treated as a fast node
- `1efc4a25f` saltation (#582): two switches flipping at one node in the same step are one event; generic recursion for capacitor-less switch nodes
- `acdaf51ce` saltation (#583): routed jumps are added after T is formed from the smooth part; no order-1 window for them
- `73aca98a3` pac td (#585): the generic recursion for every circuit with an ideal switch
- `63b0e81ab` pss shoot (#586): refuse a source that is not periodic from t = 0

## s131 — Five upstream device fixes picked from ngspice development, plus the VDMOS transient Jacobian

6 commits.

- `b42e78e5b` diode: correct instance parameter info *(upstream)*
- `ca5b45269` Remove a crash if den_coeff length is 1. *(upstream)*
- `86f8fb597` VDMOS AC/PZ bug detected by VADistiller fixed. *(upstream)*
- `543edb5a2` MOS2 derivatives fix. *(upstream)*
- `d453ac4dd` Derivatives bug in BSIM3v3 and BSIM4v5. *(upstream)*
- `c06ad3343` vdmos (#590): quasi-saturation derivative of the drain resistor in the transient Jacobian

## s133 — .stb sens=: the bias-point term and nonlinear devices through the full linearisation

1 commits.

- `24d91e38c` stb sens= (#595): the bias-point term of dT/dp, nonlinear-device parameters through the full linearisation, db/dp = 0

## s136 — td recording passes at reltol 1e-8 on circuits without switch events

1 commits.

- `3e2c83056` pac/pstb/pxf/pnoise td (#625): on circuits without switch events the recording passes run at reltol 1e-8 by default (ladder 1e-8, 1e-7, 1e-6, inherited on timestep too small); with the inherited 1e-5 the LTE grid moved with every parameter and \|T\| scattered 4e-3 rms (up to 0.15 dB) on rd_pstb_r3d, now 7e-7; r4 loop gain vs the perturbed-PSS reference at 50-200 kHz 0.08-0.46 deg -> 0.001-0.11 deg; decks with SW/CSW or XSPICE code models keep the inherited value (their event corrections were validated on that grid); PAC_RECORD_RELTOL_AUTO=0 restores the old value

## s137 — Absurd operating-point guard

1 commits.

- `60468b9df` cktop (#642): an OP stage that converges to absurd node voltages is set aside and the ladder goes on

## s138 — A linear sweep of 2 points gives both points

1 commits.

- `41f6598f1` ac/sp/stb (#651): a linear sweep of 2 points gives both points

## s139 — Transient noise follows the seed

1 commits.

- `e3d2dfdde` trnoise (#668): transient noise follows the seed -- setseed / .option seed refill the Wallace Gaussian pools from a stream derived from the seed

## s140 — .sens ac over a sweep is right at every point

1 commits.

- `eef49b3d9` sens (#672): .sens ac over a sweep of more than one point is right at every point; lin steps linearly

## s141 — .pss Floquet verdict on subcircuit netlists, .pstb probe forms

2 commits.

- `e6961c5f6` pss (#682): the Floquet re-evaluation uses the Newton columns' step rule
- `6b6923052` pstb (#717): probe=<inst>:<terminal> auto-insert and the lgprobe subcircuit
