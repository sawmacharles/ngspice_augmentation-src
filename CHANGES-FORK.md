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

## s143 — .pss with a digital clock, td .pac in discontinuous conduction, half bridge with dead time

3 commits. (`s142` was built but not released: its test run found the `.pss` change refusing a purely analog XSPICE circuit; the fix is part of the first commit below.)

- `7aa01cacd` pss (#712): period from a digital clock (d_osc / d_pwm) without driven=1; clock edge on the period boundary
- `8784b9570` pac td (#711): DCM buck at the default 200 points is within 0.03 dB of the perturbed reference (was 0.40 dB)
- `e34b94b08` pac td (#727): half bridge with dead time -- body-diode series resistance in the pre-event node equation, longer fast-node continuation

## s144 — Switch node with only device capacitance and an RC snubber in td .pac

1 commit.

- `7da64e19c` pac td (#713): charge-only switch node with an RC snubber runs (was refused); continuation jumps take the coefficient of their event

Nine open-loop asynchronous-buck decks of this kind are within 0.001-0.024 dB / 0.04-0.13 deg of the perturbed periodic steady state at 20/50/100 kHz; discontinuous conduction at 2000 points 0.033 dB / 0.08 deg (0.47 deg before the coefficient fix); closed-loop .pstb 0.002 dB / 0.045 deg. `PSS_SALT_NO713=1` restores the refusal, `PAC_NO713_COEF=1` the old coefficient timing.

## s145 — BSIM4 tnoimod noise nodes under .pnoise

1 commit.

- `e61e2e7c4` bsim4 (#729): create the tnoimod noise nodes for .pnoise and .sp donoise too (CKTnoiseAnalGiven), not only for a job named "Noise Analysis"; PSS_SALT_NO729=1 restores the old test

BSIM4 (4.5 to 4.8) created the internal source/drain nodes that carry the rs/rd part of the `tnoimod=1` channel noise (4.5/4.6: any `tnoimod != 0`) only when a plain `.noise` analysis existed. Under `.pnoise` (both engines) and `.sp donoise=1` that part was 0: a clocked BSIM4 common-source stage read 4.43 dB low against a transient-noise ensemble; now within 0.114 dB. A static stage: `.pnoise` per generator equals `.noise` to 4.7e-7. Analyses without noise are unchanged (no extra nodes).

## s146 — Capacitor at the switch node against an exact reference; no silent event drop

2 commits.

- `6d6f1e542` pac td (#728): capacitor nodes treated as fast (boost / buck-boost with a capacitor at the switch node) within 0.007 dB / 0.005 deg of an exact reference
- `97227fc47` pss salt (#728): event buffer 8 -> 32; a step that fills it, or an event with more snubber jumps than it holds, is refused instead of dropping the rest silently

Six decks (boost with 100 pF / 1 nF / 3 nF at the switch node, inverting buck-boost with 1 nF, chargeless or charged diode) were compared with an exact ODE reference (switch instants found exactly, error bar 6e-4 deg). s145 read up to 0.156 deg / 0.021 dB off; s146 0.005 deg / 0.007 dB at 20/50/100 kHz, also at 1000 and 4000 points and reltol 1e-6. The six causes, all in the parts given while the switch node is still moving, are listed in the first commit's message. `PSS_SALT_NO728=1` restores s145 exactly.

## s147 — .pnoise on non-quasi-static BSIM4 / BSIM3 models

1 commit.

- `5f92ad383` pac td pnoise (#715): accept BSIM4 trnqsmod=1 / BSIM3 nqsmod=1 in td .pnoise; PAC_NQS_PNOISE_REFUSE=1 restores the refusal

td `.pnoise` refused these models up to s146 although `.pac`, `.pstb` and `.pxf` accepted them. Neither model has a noise source on its charge-relaxation node, so the non-quasi-static model changes only the transfer to the output, which the td recursion already carries. Checked: a BSIM3 nqsmod=1 stage against its own `.noise` (which models the transient non-quasi-static behaviour) to 2.3e-6 of the total; a BSIM4 trnqsmod=1 stage per generator to 6.4e-8; a 1 MHz LO-switched BSIM4 trnqsmod=1 stage at 4-43 MHz against transient-noise ensembles, worst band 0.15 dB (resistor noise) and 0.14 dB (channel noise), where the quasi-static model is up to 3.9 dB off. Not checked: BSIM3 in a clocked circuit, flicker noise or tnoimod 1/2 together with the non-quasi-static model, and the option given on the device line instead of the model. **Known gap (all td analyses), fixed in s148:** the non-quasi-static option was detected only on the `.model` line; given on the device line of the other BSIM3 / BSIM4 versions (3.0, 3.1, 3.2, 4.5, 4.6, 4.7), which are refused when it is on the model line, the analysis ran unvalidated without a message. Circuits without such a model run unchanged code.

## s148 — Non-quasi-static option on the device line

1 commit.

- `51a47e520` pac td (#730): NQS flags on the device line are found like those on the .model line; a refused family no longer runs unvalidated when nqsmod/trnqsmod is given per instance

Up to s147 a BSIM4 4.7 device with `trnqsmod=1` on its device line (or BSIM3 3.2.4 with `nqsmod=1`) ran td `.pac` to the end without a refusal or a note, although the same option on the `.model` line is refused (those versions are not validated). Every device line is now asked too; the message names the instance. A validated BSIM4 4.8 device with the option on its line is accepted, with the note. `PAC_NQS_INSTANCE_OFF=1` restores the s147 behaviour. `.pnoise` was not affected: it refuses those device versions by its noise scope.

## s149 — Lossless transmission lines in shooting .pss and td small-signal analyses

2 commits.

- `50cfea2cc` tra (#716): periodic line history for the shooting analyses; the interpolation of the delayed excitation recorded per step
- `012d0b3d4` pss shoot, pac td (#716): lossless transmission lines in the shooting PSS and in td .pac/.pstb/.pxf

Up to s148 a lossless line (`T` element) restarted every shooting period from its initial history: `.pss method="shoot"` reported CONVERGED with an orbit 23-100 % of the peak-peak off the settled transient on RC-terminated lines, and 3.9 V off a 56 mV ripple on a buck with a line at its output; LTRA, TXL and coupled lines were accepted too (their error was not measured), and td `.pac/.pstb/.pxf` refused the circuit. The samples of the line's delayed history over one delay are now unknowns of a joint Newton with the circuit states, and of the td small-signal boundary problem. Checked: five orbits against a settled transient within 0.03-1.2 % of peak-peak (the longest-delay line is the 1.2 %, at the default 400 points; it shrinks with more points); td `.pac` against `.ac` within 1.3e-4 dB, `.pstb` against `.stb` within 6.3e-5 relative; a switching buck with an output line against a perturbed periodic steady state within 0.001 dB / 0.016°. Refused with a message: LTRA (other than R-G), TXL and coupled lines in shooting; a line in `.pnoise`, in a free-running `.pss`, or together with stateful XSPICE models; a line-table time not found among the recorded steps. `.pxf adjoint=1` falls back to the per-source engine. `PSS_TL_OFF=1` restores the s148 behaviour. Not checked: several lines in one circuit, coupled lines, a closed-loop converter with a line. Circuits without a transmission line run unchanged code.

## s150 — Small signal through digital XSPICE blocks

1 commit.

- `a5537f717` pac td, pss shoot (#714): small signal through digital XSPICE blocks; exact adc_bridge crossings

Up to s149 td `.pac/.pstb/.pxf` refused every event-driven XSPICE model; with the refusal lifted the comparator path gave exactly zero gain and decks with a `d_pwm` crashed the recording pass. The digital path is now linearised as edge-time shifts: an edge set by an `adc_bridge` threshold crossing moves by -dv/v' at the crossing; gates, flip-flops, latches and dividers carry the shift of the edge that triggered them, while an edge triggered by a fixed clock does not move; a moved `dac_bridge` edge injects (s_before - s_after)·dtau on its ramp; a `d_pwm` with an analog duty moves each transition by the shift of its previous one ± T·d'(c)·dc, carried over the period as an extra Floquet state. The `adc_bridge` crossings are located exactly in the periodic runs (`PSS_ADC_LOCATE=0` turns it off; `XS_ADC_LOCATE=1` in a plain `.tran`): a closed comparator loop's `.pss` converges in 3 iterations where s149 failed after 150. Checked against a brute-force FRA: open-loop comparator, latch, latch with an ideal switch, `d_pwm` duty and clock-only decks within 0.0009 dB at 2/25/50 kHz; a closed comparator-PWM loop's `.pstb` within 0.0012 dB / 0.0002° from 2 to 50 kHz, phase margin 38.350° at 33547 Hz vs 38.345° at 33536 Hz; two other closed loops within 0.0042 dB. Refused, naming the instance: simultaneous input edges of different origin, a multi-bit ADC reaching the analog side, a VCO, an edge or ramp still moving at the period end, a `d_pwm` whose control moves along the orbit, `.pnoise`. `.pxf adjoint=1` falls back to the per-source engine. `PAC_DIG_REFUSE=1` restores the s149 refusal. A closed loop through a `d_pwm` stays refused in `.pss`. Not checked: BSIM4 transistors together with digital blocks, the Krylov periodic solver on these decks. Plain `.tran` runs unchanged code unless `XS_ADC_LOCATE=1` is set.

## s151 — Envelope following for closed-loop converter soft-starts (`.tranenv`)

1 commit.

- `93f5bac64` tranenv (#513, #719): implicit envelope following for long converter transients, closed-loop soft-starts

A new directive, `.tranenv Fclock Points SC_iter Steady_coeff Tstop Mmax [envtol=] [phase=] [maptol=] [uic]`, simulates some switching periods and jumps over the others along the slow envelope: an implicit step over M periods solved by Newton on the one-period map, with a local error estimate against `envtol` (default 5e-4, relative to each state's scale). An earlier engine of this kind never accepted a jump on a closed-loop soft-start: at the deck's reltol 1e-3 the step control places the 1–2 ns switch edges differently from one period to the next, so the one-period map is rough (finite-difference columns 40x too small for small steps) and Newton fails. Every period map now runs at reltol = min(deck reltol, `maptol`), default 1e-5; `maptol=0` keeps the deck's tolerance and `ENV_MAPTOL_OFF=1` restores the earlier engine (the red proof: 0 jumps). Also: the confirmation period after a landing is optional (`ENV_NCHECK`, default 1); the next M grows by min(2, 0.9·sqrt(envtol/err)) instead of doubling; the Jacobian seed run no longer writes a duplicated period into the output. Measured against a plain `.tran` at 1 ns, reltol 1e-7, on every simulated period: a voltage-mode buck soft-start with a dead-time network (2000 periods) simulates 171 periods in 1.09 s, within 0.014 % / 0.16 % / 0.14 % in cycle-average output voltage / inductor current / output ripple; the deck's own `.tran 2n` takes 6.81 s and is 0.57 % / 7.6 % / 13.7 % off; the cheapest plain `.tran` at least as accurate as `envtol=1e-4` takes 6.06 s against 1.60 s (3.8x). A peak-current-mode buck soft-start draws level (1.0–1.05x at matched accuracy): its slow compensation integrator keeps the jumps short. Jumped periods are not output and not interpolated; a disturbance that starts and ends inside a jumped span is not seen (source corners are never jumped over). Plain `.tran` runs unchanged code and gave byte-identical output to s150 on three decks. Not checked: line or load ramps, a closed loop without dead time, decks whose `save` list leaves no vector for the Envelope Samples plot, memory use.

## s152 — `.pss` with an XSPICE clock: restarted period runs, a converter that is off

1 commit.

- `019777697` pss shoot, xspice (#741): d_osc init_phase 270/306 latch buck converges; no crash on a failed restart; all-off orbit accepted when the transient ends there

The test deck is an open-loop buck whose latch is set by a `d_osc` clock (duty 0.1) and reset by an `adc_bridge` comparator. Up to s151, with `init_phase` 270 or 306, the shooting `.pss` failed in its third restarted period run with "Too many analog/event-driven solution alternations", then crashed. The XSPICE operating point compared its alternation limit (hybrid outputs + 1) with a count that is zeroed once per job, and a shooting `.pss` solves the t = 0 point again in every restarted run, so the count piled up across runs. With those phases the `d_osc` toggles at t = 0 and needs a second alternation there. The limit now applies to the alternations of one operating point (`PSS_XS_ALTCOUNT_OFF=1` restores the old comparison). Both phases now converge in 3 iterations. Against a settled transient (2 ms, edges located exactly with `XS_ADC_LOCATE=1`), the cycle-average output is within 0.014 % of the ripple and the ripple within 0.001 %. The crash: on a failed transient operating point, the transient ended its output plot but kept the pointer, and the periodic run ended it a second time (use after free); the pointer is now cleared. With `init_phase` 0 or 90 the clock edge always meets an active reset, the latch never sets and the converter is off; `.pss` reported "FAILED to converge" with a residual below its tolerance, because a guard rejects every point where all states sit under the ripple floor. That point is now accepted when it is reached before any Newton step, starting exactly at the end of the warm-up transient (the transient itself is periodic there); a point Newton stepped onto is still rejected (`PSS_ALLOFF_ACCEPT_OFF=1` restores the rejection). Also: the analog state array grown for code models is zeroed, and the XSPICE snapshot carried between period runs no longer copies unused padding slots. Plain `.tran` output is byte-identical to s151 on six decks, and a `.dc` sweep with six digital transitions on one. Not checked: the all-off acceptance without a warm-up and in free-running `.pss`; a phase-90 comparison with a transient.

## s153 — td small-signal on transistor circuits: no false error message

1 commit.

- `da6980b62` pacshoot: model-charge detection once per run, no error for a recovered ladder rung (#753)

Since s134, td `.pac/.pstb/.pxf/.pnoise` on a circuit without ideal switches or XSPICE code models record their periods on a tolerance ladder: reltol 1e-8, then 1e-7, then 1e-6, then the deck's own value. The first period of each attempt also detects whether the device models carry charge. On a transistor-level BSIM4 buck the 1e-8 attempt stops after three points ("timestep too small") and the run printed `PAC-TD: model-charge detection run failed (err=106)`, although the next rung succeeded and the result was correct; detection was then repeated at 1e-7. Detection now runs once per analysis call at the deck's tolerance (an explicit `pacreltol=` or `PAC_RECORD_RELTOL` still applies), and its answer is reused by the retries; a detection or charge-state run that stops with "timestep too small" while a later rung remains prints nothing beyond the ladder's own retry note. When no rung is left the failure is reported as before. `PAC_DETECT_EACH_RUNG=1` restores the old behaviour. On the BSIM4 buck: `.pac` output byte-identical to s152, wall 30.9 s against 31.6 s; 15 periodic small-signal and noise test suites give the same values as s152, and an AddressSanitizer build is clean on two decks. Measured alongside: the 1e-7 recording that s134 introduced, not the failed attempt, is what makes this deck about 13 % slower than s133 (28.0 s with the deck's reltol 1e-4); it brings the transfer function closer to an injected-tone transient reference (worst 0.062 dB / 0.45° against 0.075 dB / 0.84°). `pacreltol=1e-4` gives the old speed and values. Not checked: a circuit where the charge-state prepass fails on every rung.

## s154 — closed loop through an XSPICE `d_pwm`; inverted `dac_bridge` inputs; stale Newton baseline

2 commits.

- `fc14db7a3` pss shoot, pac td (#751): closed loop through an XSPICE d_pwm -- the period follows the self-timed PWM; inverted dac_bridge inputs in td small signal
- `93e6de281` pssshoot: re-evaluate a stale Newton baseline before giving up an attempt (#768)

A `d_pwm` code model times each output edge from its previous edge, so in a closed loop its period follows the ripple on its control input. On a voltage-mode buck switching nominally at 500 kHz the period is 2.00865 µs, not 2 µs; stock ngspice 45.2 and s153 transients agree, and the waveform repeats every switching period (no subharmonic). A shooting `.pss` at 2 µs therefore had no steady state to find and reported "does not settle". When the circuit holds exactly one `d_pwm`, `.pss` now solves for the PWM's own period alongside the states (period lock; `PSS_PWM_LOCK_OFF=1` restores the fixed period). The period agrees with a 10 ms settled transient to 1.7e-5 relative, and the ripple to 0.001 %. Time-domain `.pac` on that orbit includes the edge-time dependence on the control input (`PAC_PWM_NO_CDOT=1` removes it); against stock ngspice's injected-tone response it agrees within 0.016 dB up to 21.5 kHz and within 0.163 dB above. A second fault, present since s150: a `dac_bridge` with an inverted input (`[~dout dout]`) took its small-signal level from the uninverted node, so td small signal through such a stage came out 12.9 dB low even in open loop; fixed (`PAC_DIG_INVERT_OFF=1` restores the old reading). The stability verdict on a period-locked orbit reads `UNKNOWN(self-timed-pwm)`: the multipliers the solver computes hold the PWM edge time fixed, so they are not the orbit's own. Not done: circuits with several `d_pwm`; `.pstb`, `.pxf` and `.pnoise` on these circuits. Not explained: on the open-loop deck at 100 kHz the injected-tone response differs from `.pac` by 2.9 dB (suspected: ringing in the tone's short measurement window, not checked).

The second commit: when a damped Newton step fails to reduce a residual that was taken on an earlier iterate's time grid, `.pss` now re-evaluates the residual at the current iterate before giving up the attempt (`PSS_REBASE_OFF=1` restores the old path). On a voltage-mode buck with a perturbed compensator, s153 accepted an orbit whose two compensator states were 0.345 mV high (3.2 % of the ripple), and the result depended on the system's math library; the orbit is now found in 7 iterations with two different math libraries. Plain transient output is byte-identical to s153 on seven decks; an AddressSanitizer build is clean on the new paths. Not done: a new deck that triggers the refusal one sensitivity test relied on (that check now skips).

## s155 — upstream sync to pre-master-48 of 2026-11-21; two regressions it would have caused fixed

86 commits: 82 upstream, 4 of this fork.

The ngspice development branch moved on from the 2026-10-07 commit this fork had picked up to. Every new upstream commit up to `47eb8ad02` (2026-11-21) is now on `main`, unchanged and with its original author, plus the bug 873 seed fix from 2026-10-05 that the earlier pick had left out. One was not taken because the fork already had the same fix (the `klu_memgrow_factor` option stored a comparison, not the value). Upstream changes that alter results:

- Diode level 1: derivative of the recombination current, generation current in reverse and breakdown, `IBV` scaled with area and `m`. BJT: quasi-saturation epi charge continuous, no division by `VO`, quasi-saturation scaled with area; `ISE`/`ISC` kept as given. BSIM4: STI stress and WPE mobility restored. MOS binning limits as in other simulators, `wnflag` handling, optional size on the outermost `lmax`/`wmax`.
- A top-level `.param` with `agauss`, `gauss`, `aunif`, `unif` or `limit` is drawn once per circuit load, as the ngspice manual documents (it was redrawn at every reference). The fork's `once()` form still works and now gives the same result as without it.
- Seeding: `setseed`, `option seed` and `.option seed` refill the transient-noise Gaussian pools in upstream's own way, so the fork's s139 helper for this is removed; seeded `agauss`/`unif` values are now upstream's (they differ from s138..s154).
- `ITL1`, `ITL2` and `ITL4` are honoured below 100 (defaults 100); breakpoint fixes (bugs 887, 888, a breakpoint at the final time); memristor `rinit` in `.ac`/`.dc`; `.save none`; `.lib` sections and nested includes; the bundled KLU is SuiteSparse 7.14.1 (was 3.7.0), or an installed SuiteSparse's KLU when one is found; OSDI header v0.5 and `osdi_bind`; CMake build next to autotools; several memory leaks and crashes.
- After a rejected time step, a retry starts again from the last accepted solution (bug 879).

Two of these broke fork features, and are fixed here:

- The retry rule above. In R4d's (peak-current-mode buck with dead time) `.pss` warm-up, every retry restarted from the accepted point at a latch edge built from B sources and diverged the same way down to "timestep too small"; a perturbed R3d `.pss` stopped converging. Retries now alternate: the first from the failed iterate (the old rule), the second from the accepted solution, and so on. Both decks converge again; the perturbed R3d gives s154's result exactly. Also, at the first point of a `uic` run upstream stored the loaded right-hand side instead of the solution as the accepted one (up to 11 apart on R4d); corrected. `DCTRAN_RETRY_ALT_OFF=1` and `DCTRAN_UIC_ACCEPT_OFF=1` restore upstream's behaviour.
- Upstream now calls each XSPICE code model's destroy callback when an analysis is set up again. Shooting `.pss` with XSPICE blocks restores the event state before every restarted period through the same setup routine, so the `d_pwm` code model lost its duty-cycle table and the `d_pwm` period lock failed. Restoring the state no longer triggers the callback; `.pac` through `d_pwm` is again byte-identical to s154. `XS_SNAP_DESTROY_ON=1` shows the failure.

Results that moved, all explained: R4d's compensation-node orbit at the default tolerance is 0.31 % of its 1.3 mV ripple from a settled transient (0.06 % on s154; at reltol 1e-5 it is within 0.10 %); the radius of an unstable orbit (R4 at duty 0.55, about 1.3) moves within the spread its own finite-difference variants already had (2.6 %); plain transients move at tolerance level (up to 0.9 % in cycle averages of a converter soft-start); checkpoint chaining against the same build's plain transient: 81 % of periods within 0.5 % (100 % on s154).

Also explained since s154: on the open-loop `d_pwm` buck at 100 kHz the injected-tone response differed from `.pac` by 2.9 dB because the tone was measured over only 5 of its periods, while the output filter (5 kHz, Q about 10) still rang from the tone's start. The difference falls as the window grows (1.0 dB at 25 periods, 0.27 dB at 100, 0.08 dB / 0.45° at 400): the `.pac` value was right.

- `22682ca62` Fix bug 873 main.c: initialise rndseed to 1 and and get random numbers by initw(). inp.c: extract options seed=val|random in netlist and script files. wallace.c: if initw is called a second time, re-allocate memory, not just overwrite it (memory leak!). *(upstream)*
- `ed6ef2159` Fix previous commit: Re-enable options in comfiles. *(upstream)*
- `e961e4d1f` Enable setting the seed by command 'setseed' also for tranoise etc. *(upstream)*
- `a8a776367` Update long comment to sharedsync() to better reflect the dctran.c code. Add similar to sharedspice.h. *(upstream)*
- `5fc25bd4c` Add prompting for new delta and return value to sync callback and fix errors in shared library example. *(upstream)*
- `2a1045b20` Fix bug #875 - ".save none sticks ...".  Using the shared library, after "save none", no vectors could be saved. Also make "save none" work in the program and fix memory leaks, dead code and simplify memory allocation. *(upstream)*
- `68b3a7375` Fix bug #874 "OPtran never finishes when a GetSyncData callback               is registered (shared library)" by removing the check for over-run in sharedsync().  That has also been seen to cause failure in transient simulation.  Make stopping more robust in the analysis code. *(upstream)*
- `e83b83f1c` Initialise random generator as in main.c *(upstream)*
- `47a67a08c` Fixes bug 876 *(upstream)*
- `f1a096277` Fixing bugs no. 877 and 878, reported by Mattia Fiumara: Enable true reverse and forward bias calculations. *(upstream)*
- `eaac23e28` dense: solve conflict of norm() definition collision for MSVS, workaround because norm() function is not used in source tree. *(upstream)*
- `f5f50907a` Add a CMake build system alongside autotools *(upstream)*
- `4aee08ea9` E-573: BJT: keep ISE/ISC as given, resolve multiplier form in BJTtemp *(upstream)*
- `f12e03748` .gitignore: extend to cmake generated files *(upstream)*
- `ab96c5c59` cmake: assume ELF on Linux (fixes Ninja install RPATH error with icx) *(upstream)*
- `23e324255` inpcom.c: evaluate .param with agauss etc. once, as documented *(upstream)*
- `e2bd0b069` Remove bug: really use rinit for ac/dc *(upstream)*
- `999ec837d` Remove linker warnins *(upstream)*
- `c9044943f` Fix blank default output file name (reported by Diarmuid Collins) and scan for .sv files as well as .s. *(upstream)*
- `49c76373a` Fix handling of missing file names in ghnggen. *(upstream)*
- `d9d71fe47` Response to bug report 879: Upon a rejected time point resort to the last accepted CKTrhs *(upstream)*
- `fe4cfd229` Remove the floor of 100 for the iteration count, so the user may select ITL1, ITL2, and ITL4. The defaults are set to 100. Reported in bug 879 by d b *(upstream)*
- `babbedd6b` typo *(upstream)*
- `8cb1dcf21` Diode level 1: breakdown current IBV is now scaled with area and m, consistent with level 3 *(upstream)*
- `291c6181b` typo *(upstream)*
- `597f7f6c8` If .option rshunt is found in the netlist, it overrides the spinit, .spceinit given values. If another netlist is sourced without .option rshunt, the global variable is removed, or reset to the spinit, .spiceinit value, if available. Same for cshunt. *(upstream)*
- `f6f1aa17f` path for circuit netlist is found automatically. *(upstream)*
- `01440fe7c` Remove unphysical ac gain: Set derivatives to 0 if limits are exceeded. Normalize the derivative in the logarithmic transition region. Add ac tests for xspice to make check. *(upstream)*
- `b1aba521e` dioload.c: fix derivative of the recombination current *(upstream)*
- `dd5148ea7` dioload.c: generation current (isr) in reverse and breakdown region *(upstream)*
- `caf262e03` Missing files for make check *(upstream)*
- `384bcda5d` bjt: make the quasi-saturation epi charge continuous in vrci *(upstream)*
- `86b9180e2` BJT: remove the division by VO in the quasi-saturation epi current *(upstream)*
- `5cdfdb1a3` BJT: scale the quasi-saturation block with the instance area *(upstream)*
- `7d6e95198` This repairs commit 50bfc21 "bsim4: reimplement mobility and binning extension" by applying STI stress and WPE reduced mobility instead only size dependent mobility. Credit to Vilem - see Discussion https://sourceforge.net/p/ngspice/discussion/127605/thread/56ed12d6f3/ *(upstream)*
- `d9aaf40c8` typo *(upstream)*
- `fc7cd0185` Correct the limits for binning of MOS model Make results equal to other simulators. Patch provided by Vilem Kledrowetz. *(upstream)*
- `600fe0788` Fixes bug 885 on option handling, improves code for bug 873. and replaces patch 136 by - Separating option handling, command 'option' in com_option.c,   command '.option seed = val' in inp.c - Adding a warning when '.option seed=xx' and 'option seed=yy' are given together with 'reset' - Resetting seed to .spiceinit values or 1, if circuits without .option seed= are loaded. *(upstream)*
- `42e6908fc` Fixes bugs 887 and 888 by setting a default 1e-15 for minBreak, reported by Andrew Sterian *(upstream)*
- `c0d697d65` Don't modify breakpoint which is at CKTfinalTime Patch 128 *(upstream)*
- `659f0611f` Catches the core dump situation when the CPL model exceeds 8 conductors and prints an error message. Provided by Rick Walker Support request 89 *(upstream)*
- `1a5e00980` Remove generated tests/xspice/analog/Makefile.in from git, it is created by autogen.sh *(upstream)*
- `cc34597b3` Add missing headers for make dist *(upstream)*
- `75102712c` When two plots are employed with different x-coordinates, the prev_xval corresponds to the second line, while it is compared with an xval of the first line. So even with monotonic plots, the direction-reversal will be triggered and an extra \n will be inserted. Fixes bug 864 *(upstream)*
- `c0189133b` CMake: no ctest with the Windows GUI binary, document it *(upstream)*
- `0f44fbf03` cmake: build without X11 if Xaw, Xmu, Xt or Xext are missing *(upstream)*
- `e7938d16e` Streamline seeding of random number generator *(upstream)*
- `13409893b` cmake: separate TMP per winflexbison run (parallel builds corrupted the parsers) *(upstream)*
- `ff52eeeb3` subckt.c: bin in double, with the same test as INPgetModBin() *(upstream)*
- `2cbbc59a8` Remove memory leak *(upstream)*
- `8d2786eda` Set the default values to 100, as the lower limit in niiter.c has been removed by commit d7e64e28d ("Remove the floor of 100 for the iteration count, so the user may select ITL1, ITL2, and ITL4. The defaults are set to 100. Reported in bug 879 by d b", 2026-10-23) *(upstream)*
- `9d3e9fd4c` Partial integration of Patch #132 by Sourceforge user George -   "Shared library: memory leak in ngGet_Evt_NodeInfo result cleanup." Also test code and a little tidying. *(upstream)*
- `2747447c2` Fix bug reported with Patch #133 -  "XSPICE: reused plot names after destroy all can reconnect stale event data." Also disable the undocumented behaviour of saving XSPICE node history for multiple simulations, unless variable "keepxs" is set. Tidy a bit. *(upstream)*
- `908e16112` Remove traces of nutmeg to improve code flavour and make "plot alle" work on non-VisualC builds. *(upstream)*
- `9e57b341e` Fix warnings from gcc 14.2. *(upstream)*
- `479d13ae2` Improve a warning message in d_cosim. *(upstream)*
- `e6c9cd0e5` In a library file, ngspice expands .include lines only inside the requested section. Nested .lib references work as before. Within one load, ngspice reuses a library file it has already read. It now reuses the file only for the same section, because a read for one section skips the includes of the others. If the library file itself lacks the requested section, ngspice rereads the library with every .include expanded, because an included file may define the section. The manual (chapt. 2.10) places each section inside the library file. The hs compatibility mode exists to read PDK libraries written for commercial simulators (chapt. 12.11), and commercial simulators reject a section inside an included file. ngspice keeps accepting one, but now prints a warning: *(upstream)*
- `3b7d2e2a9` INPgetModBin: wnflag without a value is wnflag=1 *(upstream)*
- `85b664222` subckt.c: apply wnflag when dropping binned model cards *(upstream)*
- `20a6942e7` Binning: optionally accept a size on the outermost lmax, wmax *(upstream)*
- `27b9fb697` KLU: separate the bundled SuiteSparse headers from the ngspice extensions *(upstream)*
- `f955ba55a` KLU: use the KLU of an installed SuiteSparse if available *(upstream)*
- `763093769` KLU: update the bundled SuiteSparse from 3.7.0 to 7.14.1 *(upstream)*
- `ccd106d8b` CMake: SuiteSparse package next to ngspice on Windows, install its DLLs *(upstream)*
- `89b0c0173` Remove MSVC declarations of scalbn and logb *(upstream)*
- `546f5ef1b` cmake: use Xft if available, as configure does (rotated plot text) *(upstream)*
- `3f4d262df` inpcom.c: diode with an unlabeled area: count the nodes correctly *(upstream)*
- `bc0803104` inpfindv.c: read double-quoted versions without overflow *(upstream)*
- `2a68b36e8` inpcom.c: diode node count: a value after a value is a node *(upstream)*
- `37a5f44fa` dev.c: grow the device lists of existing circuits *(upstream)*
- `da059bcac` dev.c: warn about OSDI modules with the same name *(upstream)*
- `a9be3d471` dev.c: show the .osdi file of an OSDI device *(upstream)*
- `263172a15` Revert the commits 8e03f9570 and f14a632a5, reapplied in the following commits *(upstream)*
- `98089a649` dev.c: warn about OSDI modules with the same name *(upstream)*
- `cacadb730` dev.c: show the .osdi file of an OSDI device *(upstream)*
- `d135fdb49` dev.c: same-named OSDI modules: the warning says what is affected *(upstream)*
- `c47224988` Fix three XSPICE memory leaks: cc9be5b3f mistakenly activated dead code in the ngspice program.  The bug was alreay in the shared library. Code in runcoms2.c and evtsetup.c is now removed.  Other ancient bugs: on a second or later analysis run of a circuit, event values left in output queues were lost; and code models with event nodes were reinitialied, losing previously allocated memory.  Check on variable "keepxs" moved to a more logical location. *(upstream)*
- `8e8258038` Optimize code sequence of the fourier command. Add a warning if no vector/expression is given. *(upstream)*
- `7e194a5ad`   In logicexp.c, fix lexer_scan which missed tokens, and extract_delay   could be stuck in an endless while loop.   In udevices.c, diagnostic messages could be mixed confusingly with   error messages.   In cktdest.c, CKTdestroy was missing FREE(ckt->CKTrhsAccepted) which   caused a memory leak. *(upstream)*
- `1eeeb48ca` OSDI: update osdi.h to v0.5, keep loading v0.3/v0.4 models *(upstream)*
- `d5a331825` osdi_bind: bind model levels to Verilog-A OSDI modules *(upstream)*
- `b38611c94` tests/osdi: build the dummy modules as DLL on Cygwin/MinGW *(upstream)*
- `f8f5e3cbe` trnoise: drop the #668 Wallace refill helper, superseded by upstream (com_sseed now calls initw on every seed path)
- `c0d6cb32a` dctran (#1192): retry after a rejected step -- alternate the accepted solution and the failed iterate; correct accepted solution at the first point of a uic run
- `afb999c32` xspice, pss shoot (#1192): restoring the carried event state no longer calls the code models' destroy callback
- `3af96aaae` dctran (#1192): retry order -- the failed iterate first, then the accepted solution


## s156 — free-running `.pac`/`.pstb`: line-to-output was up to 10 dB off; now within 0.07 dB of a transient

2 commits.

On a free-running converter (a hysteretic buck, `.pac ... trigger=`), td small-signal results depend on the time grid of the recording passes far more than on a clocked one: the orbit's phase (Floquet multiplier 1) amplifies their discretisation error by about f_sw/(2πf). Recorded on the deck's own 200 points per period, line-to-output on a 422 kHz hysteretic buck was 10.3 dB / 126° off a brute-force transient at 1.3 kHz, 1.8 dB / 40° at 10 kHz and still 0.6 dB / 7° at 56 kHz. The transient reference: the same circuit with a 1 V and a 0.5 V sine on the 12 V input, reltol 1e-7, a Hann-windowed Fourier component of the output; the two amplitudes agree within 0.05 dB.

When the `.pac`/`.pstb`/`.pxf`/`.pnoise` card sets no `pacpoints=`, a free-running orbit's recording passes now use at least 4000 points per period (`PAC_FREERUN_POINTS`, a notice says so; the shooting PSS keeps its own `points=`). Line-to-output is then within 0.07 dB / 0.2° of the transient from 1.3 to 56 kHz, for about 0.02 s more per sweep; control-to-output, validated before at 20 and 60 kHz, moves slightly closer to its reference (0.04 dB / 0.47° and 0.002 dB / 0.05°). A free-running grid below 4000 points prints a warning, and so does a frequency below 0.3 % of f_sw (lower than anything checked). `PAC_FREERUN_POINTS=0` restores the old grid. New test suite: `test_freerun_pac` (59 suites).

- `318dc3d0d` pacshoot (#1757): free-running .pac/.pstb warns once per run that the phase mode is not treated, and per point below 4.7 % of fsw
- `4157ab650` pacshoot (#1758): free-running .pac/.pstb records on at least 4000 points per period by default


## s157 — MOS binning: bin edges are reported, `.option binlegacy` restores the older choice

1 commit.

Since s155 (upstream 03ef19244, 2026-10-30) a size is binned by the standard rule min <= size < max, as in other simulators. Before, a size within 1e-9 of a limit was in the bin and the first card found won, so a device exactly on a bin edge got whichever card came first (in hs/spe mode the subcircuit expander first kept only the cards with min <= size < max in float, on w/nf). The choice matters: a synchronous buck of sky130 5 V devices, all six exactly on a bin edge, finished on s154 and stopped on s155/s156 with "Timestep too small" at a gate edge. Both cards give the same answer when the run gets through (with `.option rshunt=1e12`), so the stop is a convergence weakness that the other card exposes, not a wrong device.

- A warning for the first 5 instances per circuit that sit on a bin edge where the two rules take different cards, naming the card used and the other one.
- `.option binlegacy` (or `set binlegacy`) bins by the older rule, including its card pruning in hs/spe mode; that buck then finishes with the s154 results to all printed digits.
- The default stays the standard rule. New test suite: `test_binedge` (60 suites).

Tried and left out: a transient rescue that retries a would-be "Timestep too small" with a temporary shunt on every node. On that buck it never got past the failing instant.

- `689a8c6d2` binning (#1766): .option binlegacy restores the bin choice before 2026-10-30; warn at bin edges where the rules differ


## s158 — `.noise` input-referred total; `trrandom` on current sources

2 commits.

`inoise_total` is the integral over frequency of the output noise divided by the gain. Every noise source integrated each frequency interval with the gain of the interval's upper end at both ends, so the total was wrong wherever the gain changes over an interval, and it depended on the number of points per decade. A resistive divider with a filter capacitor, whose input-referred density is flat, read 6.459960e-06 instead of the plain divider's 5.757760e-06 (examples/noise/resdiv.cir); examples/various/bjtnoise.cir read 4.701227e-06 against 4.384e-06 by hand. Each point now keeps its gain and the lower end of the interval uses it: 5.757760e-06 and 4.384172e-06, the same at 10 and at 1000 points per decade. `onoise` and the spectra do not change. This covers every device's noise code (34 files) and the noise of `.sp`.

A current source with `trrandom` drew a new value only when a time point fell within 3 ulps of a multiple of TS while a breakpoint was taken, and never set its first breakpoint; next to a voltage source with the same TS (whose breakpoints drift by the minimum break per period) it stopped after two values and stayed constant. It now keeps the time of the next value and sets a breakpoint at each one, like the voltage source; with a seed, a current source alone gives the voltage source's sequence.

New test suite: `test_inoise_trrandom` (61 suites).

- `ed57bcb2a` noise (#1830): integrate the input-referred noise with the gain of each end point
- `8077f46d2` isrc (#1830): trrandom draws a new sample every TS, as the voltage source does


## s159 — `ddt()`; TXL and CPL lossy lines

2 commits.

`ddt()` in B sources was evaluated once per new time point, at the first Newton iteration, as (y_n - y_n-1) / (t_n-1 - t_n-2), with a zero Jacobian: the values were not converged, the step was the previous one, a rejected point stayed in the history, and the result lagged one step (1.05 V at the 2 V edges of examples/ddt/ddt-test3.cir, where ddt(v(2))/1000 must equal v(1)-v(2)). A B-source capacitor `i = C*ddt(v)` behind a diode stopped with "Timestep too small". Now the history holds the converged argument at accepted points, the result is computed at every iteration by the formula capacitors use (trapezoidal order 1 and 2, gear up to 6), and its Jacobian goes into the matrix, so `C*ddt(v(a,b))` equals the current of a capacitor C: ddt-test3 to 5e-9 V with trapezoidal and gear; the B-source capacitor runs and matches a real one within 0.01 % of the peak with rejected steps, also in `.pss` shooting (the `ddt()` history is not a shooting state, so `.pss` takes one retry there). `.option ddtlegacy` (or `set ddtlegacy`) restores the old `ddt()`.

The TXL and CPL lines kept their history times in whole picoseconds: a step that did not advance a whole ps skipped the history update, slopes used the truncated step, and the time overflowed after 2.1 ms. The convolution update also multiplied each term's c/x into the slope, so later terms used the product of all earlier ones. A step into an open TXL line settled to 1.00036 instead of 1 and drifted with the step size; a `txl` model with R/L above 1.6e9 (which runs as a single CPL line) settled 10 % low at 1 ps steps and 2.5 % high at 10 ps; the late tail of examples/TransmissionLines/cpl3_4_line.sp was 18 times a 400-section coupled ladder's. Now: exactly 1, within 0.3 %, within 8 %. The three-pole fit of both lines loses accuracy as the loss R*length/(2*Z0) grows (txl1_1_line.sp, an exact line as reference: 1 % at 0.23, 4 % at 0.35, 19 % at 0.72); a `txl` model above 0.2 now gets a warning that names the O element (LTRA), exact for an RLC line. Still open: cpl3_4_line.sp at steps of 5 ps and below (about 10 % off, as before).

New test suites: `test_ddt`, `test_txl_cpl` (63 suites).

- `2b934a7eb` ddt (#1831): derivative of the converged waveform by the circuit's integration formula, with its Jacobian
- `ff152499d` txl, cpl (#1832): history times as doubles, per-term slope in the convolution update; warn for lossy TXL models


## s160 — XSPICE microstrip lines (`mlin`, `cpmlin`)

1 commit.

`src/xspice/tlines/tline_common.h` defined `sqr(x)` as `(x*x)` (also `cubic`, `quadr`) with no parentheses round x, so `sqr(er - 1)` evaluated to -1 whatever er was, and `sqr(1 + 2*log10(1 + W/h))` to 1 + 4*log10(1 + W/h). That hit the Hammerstad-Jensen coupled-line quasi-static model (the r term), the Kirschning-Jansen coupled-line impedance dispersion (q16, q25, q29) and the Yamashita single-line dispersion. The coupled-line Kirschning P1 also read 0.27488 * (...) where the single-line P1 in the same file and the published formula have 0.27488 + (...). On examples/xspice/transmission_lines/cpmline_ac.cir |S31| at 4.2 GHz moves from 0.254 to 0.290. Against the model equations written out with real squares (lossless lines, closed-form impedance matrix), s159 was 9.0 % off on `cpmlin` (model=0, disp=0) and 10.5 % on `mlin disp=2`; s160 agrees to 1e-12. The other `mlin` dispersion models did not change. Reported on the public example deck by a downstream user.

New test suite: `test_tline_macros` (64 suites).

- `8b6727717` tlines (#1885): parenthesise the sqr/cubic/quadr macros; coupled-line Kirschning P1 with +


## s161 — BSIM4 version 4.5: negative effective junction perimeter

1 commit.

With `pd = 0` (or `ps = 0`) given and `permod = 1`, the BSIM4 4.5 code (`bsim4v5`) took the effective perimeter as PD - Weff*nf, which is negative. The sidewall junction capacitance cjswd*Pdeff was then negative; the charge code leaves a negative term out in reverse bias but adds it in forward bias, so the junction capacitance jumped at zero bias (with a large cjswd it went negative). sky130's 5 V subcircuits pass `pd = 0` and `ps = 0` with `permod = 1`. On a 5 V sky130 synchronous buck a step and Newton trace showed Newton cycling on the low-side switch's internal drain-body node, the junction voltage alternating between -0.030 and +0.013 V at every step size down to 1e-16 s, and the transient stopped with "timestep too small" in the dead time. The bin-edge card choice taken in at s155 only changed the trajectory: with either card, W moved 1 ppm off the edge also stopped, and a 1e12 ohm resistor on a node that already had 5 ohm to ground let it finish. BSIM 4.6.2 and later (ngspice's `bsim4v6`, `bsim4v7`, `bsim4`) set a negative effective perimeter to zero; `bsim4v5` now does the same. Reverse-bias junction charge is unchanged. `BSIM4V5_PERIM_LEGACY=1` in the environment keeps the old values.

New test suite: `test_bsim4v5_perim` (65 suites).

- `10dc5fa98` bsim4v5 (#1793): a negative effective junction perimeter is set to zero, as in BSIM 4.6.2 and later


## s162 — coupled CPL lines (regression of s159-s161)

1 commit.

s159 changed the CPL real-pole history update to the exact piecewise-linear convolution with each term's own c/x. The line predicts the present step with the trapezoidal rule, so the two no longer matched, and coupled lines (two or more) drifted and went unstable as the step shrank: two coupled lines with the parameters of examples/TransmissionLines/cpl3_4_line.sp read 5.05 V at 30 ns at tmax 50 ps, 5.40 V at 5 ps and 12.1 V at 1 ps, from a 5 V source (s158 and stock ngspice 4.99 V, a 400-section coupled ladder 4.985 V). The history now uses the same trapezoidal rule as the present step: 4.9856 to 4.9930 V from 50 ps to 1 ps. Single CPL lines and TXL lines do not use this path and are unchanged; the s159 note that CPL settled values were "now right" held for single lines only. The 4-line deck itself still depends on the step below 10 ps, as in stock ngspice (a different path). `CPL_HIST_S159=1` in the environment restores the s159-s161 update. Coupled-line results from s159, s160 and s161 should be re-run.

Test suite `test_txl_cpl`: two new checks.

- `be580ed80` cpl (#1850): advance the history with the load's trapezoidal rule; coupled lines no longer go unstable at small steps


## s163 — CPL lines independent of the time step

1 commit.

The CPL line predicts each step with the trapezoidal rule; s162 made the real-pole entries' history follow that rule. The entries whose fit has complex poles kept one exception: their real pole was advanced with the exact piecewise-linear update (the complex pair already used the trapezoidal rule). Lines with three or more conductors, whose fits have complex poles, therefore settled to a wrong level that moved with the step: examples/TransmissionLines/cpl3_4_line.sp, v(6) at 30 ns, read 5.034 / 5.038 / 4.508 / 6.584 / 3.716 V at tmax 50 / 10 / 5 / 2 / 1 ps; stock ngspice gives the same. A 400-section coupled lumped ladder gives 4.985 V. That real pole now uses the trapezoidal rule as well: 4.9851 to 4.9855 V at all five steps, the late tail (49 ns) within 1.6 % of the ladder. Single CPL lines (and TXL models with R/L above 1.6e9, which run as one) now settle exactly. `CPL_IMG_EXACT=1` in the environment restores the old update.

Test suite `test_txl_cpl`: the coupled-line check now runs at 50, 5 and 1 ps with 0.1 % on the ladder; one red-proof check added.

- `588873f58` cpl (#1850): advance the real pole of complex-pole entries with the load's trapezoidal rule; 3- and 4-line CPL no longer depend on the step
