# ngspice-augmentation

A fork of the [ngspice](https://ngspice.sourceforge.io/) circuit simulator that adds what
analog and power-electronics designers otherwise have to do by hand or in a commercial tool:
loop-stability analysis, fast periodic steady state of switching circuits, and periodic
small-signal analyses (transfer function, loop gain, noise). It also fixes defects in stock
ngspice that give wrong answers without a message, and makes common simulations faster
without changing their results.

- **Base:** the ngspice development branch (`pre-master-48`) at commit `aaba50ba3f49`
  (2026-09-27), kept on the branch [`upstream-base`](../../tree/upstream-base). Five later
  upstream fixes were picked on top.
- **This fork:** 151 commits on top of that base, the branch `main`. Releases are the tags
  `s125` … `s149`; `s149` is the latest.
- **Size of the change:** 96 files, about 24 000 lines added, 231 removed.
- **License:** unchanged, the ngspice licenses (see `COPYING`). Upstream authorship is kept in the
  history.
- The full list of changes is in [CHANGES-FORK.md](CHANGES-FORK.md). The upstream ngspice
  README is in [`README`](README).

## Why this fork exists

Stock ngspice is a good transient and AC simulator, but three things a designer needs every
day are missing or slow:

1. **Is my loop stable?** Stock has no stability analysis. Phase margin is found by hand,
   by breaking the loop and running two AC analyses and combining them, or by guessing.
2. **What does my switching converter settle to, and what is its loop gain?** A DC-DC
   converter has no DC operating point to linearise around. Its state repeats every switching
   period. The only stock route is a long transient (you guess when it has settled) and then
   one more transient per frequency with an injected tone. That takes minutes to hours and
   loses accuracy toward half the switching frequency.
3. **Can I trust the answer?** While building this fork we found several stock defects that
   return a wrong number with no warning (listed below). A wrong margin or ripple value is
   worse than none.

The fork's approach: add native analyses where stock has none; refuse with a message when an
analysis cannot be done correctly, rather than answer silently; check every result against an
independent method before shipping it.

## What was added

### 1. Loop stability: `.stb`

A native loop-gain analysis using Tian's two-injection method. It is exact for a loop broken
at a single point, including loading in both directions.

```spice
Vprobe fb fbtop DC 0                      ; a 0 V source spliced into the loop
.stb dec 50 1 1e9 probe="Vprobe"
```

- Reports the loop gain, crossover frequency, phase margin and gain margin. It warns on
  positive feedback at low frequency, and labels rising and falling phase crossings.
- Supports differential and common-mode loops (`probe2=`, `mode=dm|cm`), a local ground
  (`localgnd=`), and automatic probe insertion (`probe=<instance>:<terminal>`).
- **Device mode:** `probe=<MOSFET instance>` gives the loop gain around a transistor's own
  transconductance (gain nulling), for loops that cannot be broken at a wire. Works with
  MOS1 and BSIM4.
- **Sensitivities:** `sens=all|<name>` gives the derivative of the loop gain, phase margin,
  crossover and gain margin with respect to every parameter, from one adjoint solve.
- **`instant`:** linearise at a time point of a transient instead of the DC operating point.
- Bad arguments (`n`, `fstart`, `fstop`) are refused instead of silently replaced by
  defaults.

### 2. Periodic steady state by shooting: `.pss ... method=shoot`

Finds the periodic orbit of a driven circuit (such as a converter clocked by a source)
directly, with a shooting-Newton solver, instead of running a transient until it looks
settled. It reports convergence, so there is no guessing.

```spice
.pss 500e3 300e-6 vout 200 10 150 1e-6 method="shoot"
```

- Staged warm-up, Broyden updates, and a guard against degenerate fixed points. Fast states
  are reduced automatically.
- **Large circuits:** matrix-free Newton-GMRES (`krylov=1`), automatic from 64 states. On a
  1150-state extracted buck the analysis went from 152 s to 8.5 s, and peak memory on a
  parasitic-extracted buck fell from 1.45 GB to 40 MB.
- **Free-running circuits** (oscillators, hysteretic converters): `trigger=<node>` finds the
  period as well. Floquet stability is reported, with an `UNKNOWN` when the check is not
  conclusive.
- **Ideal switches (`sw`):** the exact switching instant is located, and a saltation
  (jump) correction is applied at each switch event. These events are where small-signal
  results of switching converters go wrong if handled naively.
- Carries hidden state across the period restart: device charges in transistor-level
  circuits, and the state of XSPICE code models.
- **Digital clocks:** a converter whose only clock is an XSPICE `d_pwm` or `d_osc` block
  takes its period from that block (open loop: orbit within 0.05 % of the ripple of a settled
  transient). A loop where the analog control sets the `d_pwm` duty is refused with a
  message, because the edge time is state the solver does not carry.
- Refuses a source that is not periodic from t = 0, rather than converging to the orbit of a
  different circuit.

### 3. Periodic small-signal analyses: `.pac`, `.pstb`, `.pxf`, `.pnoise` with `solver="td"`

Small-signal analysis around the periodic orbit, the way commercial RF and power simulators
do it. This is what a converter designer needs: control-to-output transfer function, loop
gain and margins, and noise of a circuit that never sits at a DC point.

```spice
.pac  100e3 100e-6 vx 20000 1 25 1e-6 dec 10 1k 130k method="shoot" solver="td"
.pstb 500e3 300e-6 vout 200 10 150 1e-6 dec 5 500 200000 probe="Vprobe" method="shoot"
.pxf  ... method="shoot" solver="td" sources=... adjoint=1
.pnoise ... OutNode InSrc dec 10 1k 1meg method="shoot" solver="td"
```

- **`.pac`**: periodic AC (transfer function through the switching circuit).
- **`.pstb`**: periodic loop gain and margins (Tian's method on the periodic orbit). The probe
  takes the same forms as `.stb`: a 0 V source, `probe="<inst>:<terminal>"` (the source is
  inserted for you), or an instance of the `lgprobe` subcircuit; all three give the same loop gain.
- **`.pxf`**: periodic transfer function from every source to one output. `adjoint=1`
  computes all of them in one transposed solve.
- **`.pnoise`**: cyclostationary noise. The noise of each device is evaluated at the bias it
  has at each instant of the orbit. Flicker noise is folded correctly from every sideband.
- Handles ideal-switch converters (buck, boost, half bridge, with or without dead time,
  synchronous or diode rectified), device-model charges (MOS, BSIM3/4 including their
  non-quasi-static models, diode, BJT), lossless transmission lines (`T` element), stateful
  XSPICE analog models, and free-running orbits.
- Says so when it cannot answer correctly: for example discontinuous conduction at a coarse
  step, lossy transmission lines (LTRA, TXL, coupled lines) in shooting `.pss`, a lossless line in
  `.pnoise`, or small-signal through digital blocks.

### 4. Speed (results unchanged)

Every speed change was required to give the same transient results as the unpatched base
(identical printed measures; waveforms identical or within 8e-9 of range).

- Resistor and capacitor matrix stamps are cached and copied instead of re-scattered every
  iteration. A flat capacitor load is used under KLU.
- Behavioural (B-source) expressions are compiled once into a flat program.
- Better OpenMP scheduling for the BSIM3, BSIM4, HiSIM2 and BSIM-SOI load loops.
- KLU refactors without scaling (pivots still chosen on the scaled matrix). KLU is chosen
  automatically for shooting analyses.
- Each speedup can be switched off at run time, for example with `NGSPICE_FLATB=0`,
  `NGSPICE_CAPBASE=0` or `NGSPICE_CAPFLAT=0`.

### 5. Robustness

- **Transient checkpoint and restart:** a resumed run is bit-exact with an uninterrupted one.
- **Absurd operating points:** a DC operating-point stage that "converges" to absurd node
  voltages (or non-finite ones) is set aside, and the next strategy is tried. The absurd point
  is kept, with a warning, only if nothing else converges. `OP_ABSURD_GUARD=0` restores stock
  behaviour.
- **Random parameters:** `.param` expressions can draw a random value once, explicitly.

## Stock ngspice defects fixed

Each of these gives a wrong result in stock ngspice 45.2 without an error message.

| Defect | Effect in stock | Fixed in |
|---|---|---|
| `reset` drops circuit options | A temperature sweep in one session runs at the wrong temperature. In one test, stock ran all 50 points at 0 °C when 25 should have been at 85 °C (30 % error on a ripple value). | s125 |
| `.sens ac` over more than one frequency | Resistor sensitivities are exactly 0 after the first point. A `lin` sweep steps geometrically (1e4, 7.5e7, 5.6e11 …). | s140 |
| Transient noise (`trnoise`) ignores `.option seed` / `setseed` | Seeded runs never repeat. Now two seeded runs give identical noise, and other random draws are unchanged. | s139 |
| A linear `.ac`, `.sp` or `.stb` sweep of 2 points | Only the first point is computed. | s138 |
| VDMOS quasi-saturation | The drain-resistor derivative is missing in the transient Jacobian (upstream fixed the AC part; the fork also fixed the transient part). | s131 |
| `klu_memgrow_factor` option | Stored as 0/1 instead of its value. | s125 |
| Frequency-domain `.pss` / `.pxf` with `harmonics=1` | Reads past a buffer. | s125 |
| Netlist parser: a dot-card line ending in a quote | Writes past the end of the line and corrupts the heap (crash at quit: `free(): invalid next size`). | s125 |

## Measured results

Summary measured on 2026-10-08 on one 32-thread Linux workstation. The fork build was `s125`, compared
with stock ngspice 45.2; the transient table below was re-measured on `s140`. Two threads per run, median of 5 runs (3 for the periodic part).

| Task | Stock ngspice | Fork | Gain |
|---|---|---|---|
| Plain transient, 11 decks (see the table below, re-measured on s140) | — | same deck, same settings | 1.06x – 1.71x faster (s140: 1.02x – 1.80x) |
| Steady state of a switching converter (versus a transient long enough to settle) | 1.1 – 109 s | 0.3 – 37 s | 3x – 62x |
| Control-to-output transfer function, 8 frequencies (stock: injected-tone transient sweep) | 31 – 838 s | 0.1 – 28 s | 30x – 275x |
| Loop gain of a closed-loop buck, 14 frequencies | 144 s | 0.32 s | 450x. Same crossover (33.5 kHz) and phase margin (38.4° vs 38.45°) |

### Transient speed: stock ngspice vs this fork

Measured 2026-10-18 with the latest release `s140`. The same deck text was run on three
builds: **stock** ngspice 45.2 (the current release), **base** (the unpatched upstream commit
this fork starts from, which separates upstream changes from ours), and the **fork**. Runs
were interleaved and rotated between builds, 5 runs per deck per build, 2 threads, one
simulation at a time, on a 32-thread Linux workstation that had other jobs running (load 2–6).
Wall time in seconds: median, with [min – max] below it.

| Deck | What it stresses | Stock 45.2 | Base | Fork s140 | Fork vs stock | Fork vs base |
|---|---|---|---|---|---|---|
| e_pex_tree | 11 111-node RC tree + 60 BSIM4, KLU, 25 µs | 5.69<br>[4.17 – 14.01] | 5.69<br>[4.20 – 12.43] | **3.17**<br>[2.70 – 6.16] | **1.80x** | 1.80x |
| c_pvt_loop | 50 temperature/corner points in one session | 3.96<br>[3.86 – 4.03] | 4.05<br>[3.91 – 4.15] | **2.29**<br>[2.23 – 2.36] | **1.73x** | 1.77x |
| h_pex_tree_sparse | smaller RC tree + BSIM4, KLU, 20 µs | 0.70<br>[0.69 – 0.73] | 0.71<br>[0.68 – 0.76] | **0.47**<br>[0.45 – 0.49] | **1.49x** | 1.51x |
| b_long_switch | peak-current-mode buck, behavioural switches, 20 ms (10 000 periods) | 15.87<br>[15.52 – 16.49] | 17.09<br>[16.71 – 17.46] | **11.35**<br>[11.08 – 11.99] | **1.40x** | 1.51x |
| r3_vm_buck | voltage-mode buck, 1.5 ms | 1.45<br>[1.40 – 1.50] | 1.48<br>[1.44 – 1.52] | **1.17**<br>[1.13 – 1.22] | **1.24x** | 1.27x |
| r4_pcm_buck | peak-current-mode buck, 6 ms | 7.73<br>[7.38 – 7.78] | 8.15<br>[7.92 – 8.26] | **6.36**<br>[6.09 – 6.43] | **1.22x** | 1.28x |
| r3d_vm_buck_deadtime | buck with dead time, 1.5 ms, 2 ns max step | 8.55<br>[8.28 – 8.74] | 8.74<br>[8.37 – 8.92] | **7.35**<br>[7.02 – 7.51] | **1.16x** | 1.19x |
| a_pex_buck | 12 100-node RC grid + 60 BSIM4, KLU, 25 µs | 59.04<br>[52.36 – 69.89] | 55.77<br>[46.89 – 66.72] | **53.32**<br>[52.40 – 61.40] | **1.11x** | 1.05x |
| r7_bsim_buck | BSIM4 transistor-level buck, 400 µs | 14.59<br>[14.42 – 14.63] | 15.01<br>[14.90 – 15.22] | **13.18**<br>[12.96 – 13.54] | **1.11x** | 1.14x |
| f_analog_block_bsim3_long | 600 BSIM3 transistors, 50 µs | 31.49<br>[30.69 – 31.99] | 31.38<br>[30.59 – 32.21] | **29.69**<br>[29.10 – 30.33] | **1.06x** | 1.06x |
| d_analog_block_long | 300 BSIM4 transistors, 50 µs | 20.79<br>[20.31 – 21.32] | 21.89<br>[21.72 – 22.04] | **20.38**<br>[20.00 – 20.62] | **1.02x** | 1.07x |

**Same answers?**
- **Fork vs base: yes.** Every printed measurement is identical on all 11 decks (300 values),
  and the waveforms are identical or differ by at most 1.2e-8 of their range. The speedups
  do not change results.
- **Fork vs fork, run to run:** identical.
- **Fork vs stock 45.2:** printed values differ by 1.7e-6 to 1.2 % (the largest is output ripple
  on the switching bucks). That difference comes from upstream ngspice changing between 45.2
  and its development branch, not from this fork: base and fork are equal. Switching edges
  land on different time steps. This table does not say which of the two is closer to the truth.
- **c_pvt_loop is the exception.** Stock 45.2 ran all 50 points at 0 °C because of its `reset`
  bug, while 25 points should be at 85 °C (up to 31 % difference in one ripple value). The fork
  runs them at the right temperatures.

**Reading:** the gain is largest where time goes into resistor/capacitor loading, behavioural
sources and per-point overhead (RC trees from layout extraction, switching converters, many
short runs): 1.2x to 1.8x. It is smallest where the time goes into the transistor model itself
(large BSIM blocks): 2–7 %. Memory use is the same on all three builds (within 2 %). On the two
large extracted-RC decks the min–max spread is wide because of the other jobs on the machine;
compare their medians with that in mind.

**Accuracy** was checked against an independent reference: a perturbed periodic steady state,
which uses no small-signal code. At the frequencies checked (10, 21.4 and 100 kHz on three bucks), the fork's `.pac` is within 0.044 dB and 0.18° of it. The
stock injected-tone route drifts toward half the switching frequency. On an ideal-switch buck at
100 kHz it was 3.5 dB and 46° off, and its answer moved by 11 dB with the injection amplitude.

**Where the fork is not better:**
- Plain transient on device-dominated circuits gains only 6–11 %.
- Transient accuracy is not claimed to be better than stock. The fork reproduces its upstream
  base exactly, and both differ from 45.2 by about 1 % of ripple on switching decks.
- A closed loop in which an analog control signal sets an XSPICE `d_pwm` duty cannot run in
  `.pss` (refused); open-loop digital clocks work since `s143`.
- **Earlier silent errors in the fork itself.** Builds before `s123` gave wrong `.pac`/`.pstb`
  results on several ideal-switch converters without a message. They were found by checking
  52 decks against the perturbed reference and fixed in `s123`. Since then all 47 open-loop
  decks agree within 0.067 dB and 0.32°, except three documented cases (up to 0.40 dB);
  since `s143` the discontinuous-conduction case is within 0.03 dB at the default 200 points.
  Builds before `s143` read the control-to-output gain of a half bridge with dead time 0.04 dB
  high (the body diode's series resistance was left out of the switch-event jump).
  Builds before `s144` refused a switch node whose only charge is the devices' (diode or body-diode
  capacitance, no capacitor at the node) when it had an RC snubber; it runs since `s144`, within
  0.024 dB and 0.13° on nine such decks, and a slow turn-off edge no longer delays the phase
  (0.47° to 0.08° at 100 kHz in discontinuous conduction).
  Builds before `s146` read the phase of a boost or buck-boost with a capacitor at its switch
  node up to 0.16° off (gain up to 0.02 dB) against an exact reference; since `s146` within
  0.007 dB and 0.005°. A switching step with more switch-event parts than the engine holds is
  now refused; before `s146` the extra parts (beyond 8) were dropped without a message.
  Builds before `s148` looked for a non-quasi-static option (`nqsmod`, `trnqsmod`) only on the
  `.model` line: given on the device line of a BSIM3 3.0–3.2 or BSIM4 4.5–4.7 device (versions
  refused with the option on the model line), a td small-signal analysis ran unvalidated without
  a message. Since `s148` the device line is checked too and such a device is refused.
  Builds before `s149` reported `.pss` with a lossless transmission line as converged while each
  period restarted the line from its initial history: on test circuits the orbit was 23–100 % of
  its peak-peak off the settled transient (a buck with a line at its output: 3.9 V off a 56 mV
  ripple); LTRA, TXL and coupled lines were accepted too (their error not measured). Since `s149` the line's
  delayed history is part of the periodic state: within 0.03–1.2 % of peak-peak on five circuits,
  and td `.pac`/`.pstb`/`.pxf` run on it (`.pac` against `.ac` within 1.3e-4 dB); lossy lines are refused.
  Builds before `s145` left out the source/drain part of the BSIM4 `tnoimod=1` channel noise
  in `.pnoise` (about 4.4 dB low on a test circuit): the model created those internal noise
  nodes only when the deck also had a plain `.noise` analysis. Since `s145` `.pnoise` (and
  `.sp` with noise) creates them too and agrees with a transient-noise ensemble within 0.3 dB;
  on older builds an inert `.noise` card on a separate node is a workaround.
  Builds before `s141` could report a stable converter as `UNSTABLE(period-doubling)`
  (Floquet multiplier −1.27 instead of 0.994) when the same circuit was written with
  subcircuits; the stability check now uses the same finite-difference step rule as the
  Newton iteration, and both forms read 0.994.

## How it is tested

Each release was run through the full set of test suites (56 at `s149`: stability, periodic analyses, noise,
sensitivities, XSPICE, speed, checkpointing, Monte Carlo and yield tools), also on a build with
AddressSanitizer, plus about 210 checks of documented example results. A fix comes with a
"red proof": its test is run on the previous build, or with the fix switched off, and must
fail there. These suites and the analysis tools built around the fork live in a separate repository
that is not public yet.

Not tested: Windows and macOS builds (Linux only); comparison with a commercial simulator;
an analytic reference for plain transient.

## Building

A Linux x86-64 binary of `s149` is attached to the [s149 release](../../releases/tag/s149)
(extract anywhere, run `bin/ngspice-aug`). To build, Linux, same as stock ngspice:

```sh
./autogen.sh
mkdir release && cd release
../configure --with-x --enable-klu --enable-xspice --enable-openmp CFLAGS="-m64 -O2"
make -j8
sudo make install      # or --prefix=... to keep it beside a system ngspice
```

Syntax notes:
- On a netlist card, put `probe=`, `method=` and `solver=` values in quotes, as in the
  examples above. The interactive commands do not need the quotes.
- The fork can sit beside a stock ngspice: install it to its own prefix.

## Status and history

Development has been active since 2026-09-27. Commit messages refer to internal tracker numbers
(`#458`, `card #465` …); those trackers are not public, but each commit message describes its
change on its own. The tags are:

| Tag | Main change |
|---|---|
| s125 | `.pss` shooting, `.pac/.pstb/.pxf/.pnoise` td, `.stb`, ideal-switch saltation, the speed rounds |
| s131 | five upstream device fixes (VDMOS, MOS2, BSIM3v3/BSIM4v5 derivatives, diode info) plus the VDMOS transient Jacobian |
| s133 | `.stb sens=`: bias-point term and nonlinear devices through the full linearisation |
| s136 | td recording passes at reltol 1e-8 on switch-free circuits: loop gain no longer scatters |
| s137 | absurd operating-point guard |
| s138 | 2-point linear sweeps |
| s139 | `trnoise` follows the seed |
| s140 | `.sens ac` over a sweep |
| s141 | `.pss` Floquet verdict on subcircuit netlists; `.pstb` probe forms of `.stb` |
| s143 | `.pss` with a digital clock; td `.pac` in discontinuous conduction; half bridge with dead time (`s142` was not released: its test run found a regression, fixed in `s143`) |
| s144 | td `.pac` on a switch node with only device capacitance and an RC snubber (was refused); continuation jumps at the event's coefficient |
| s145 | BSIM4 `tnoimod=1` channel noise complete in `.pnoise` |
| s146 | boost / buck-boost with a capacitor at the switch node within 0.005° of an exact reference; full event buffer refused, never dropped |
| s147 | `.pnoise` on non-quasi-static BSIM4 / BSIM3 models (was refused), checked against transient-noise ensembles |
| s148 | a non-quasi-static option on the device line is checked like one on the `.model` line (unvalidated versions ran without a message) |
| s149 | lossless transmission lines in `.pss` shooting and td `.pac/.pstb/.pxf` (the orbit was wrong without a message); lossy lines refused |
