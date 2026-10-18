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
- **This fork:** 138 commits on top of that base, the branch `main`. Releases are the tags
  `s125` … `s140`; `s140` is the latest.
- **Size of the change:** 89 files, about 22 000 lines added, 143 removed.
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
- **`.pstb`**: periodic loop gain and margins (Tian's method on the periodic orbit).
- **`.pxf`**: periodic transfer function from every source to one output. `adjoint=1`
  computes all of them in one transposed solve.
- **`.pnoise`**: cyclostationary noise. The noise of each device is evaluated at the bias it
  has at each instant of the orbit. Flicker noise is folded correctly from every sideband.
- Handles ideal-switch converters (buck, boost, half bridge, with or without dead time,
  synchronous or diode rectified), device-model charges (MOS, BSIM3/4, diode, BJT), stateful
  XSPICE analog models, and free-running orbits.
- Says so when it cannot answer correctly: for example discontinuous conduction at a coarse
  step, non-quasi-static MOS models, or small-signal through digital blocks.

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

Measured on 2026-10-08 on one 32-thread Linux workstation. The fork build was `s125`, compared
with stock ngspice 45.2. Two threads per run, median of 5 runs (3 for the periodic part).

| Task | Stock ngspice | Fork | Gain |
|---|---|---|---|
| Plain transient, 11 decks (PEX RC trees, switching converters, BSIM blocks) | — | same deck, same settings | 1.06x – 1.71x faster |
| Steady state of a switching converter (versus a transient long enough to settle) | 1.1 – 109 s | 0.3 – 37 s | 3x – 62x |
| Control-to-output transfer function, 8 frequencies (stock: injected-tone transient sweep) | 31 – 838 s | 0.1 – 28 s | 30x – 275x |
| Loop gain of a closed-loop buck, 14 frequencies | 144 s | 0.32 s | 450x. Same crossover (33.5 kHz) and phase margin (38.4° vs 38.45°) |

**Accuracy** was checked against an independent reference: a perturbed periodic steady state,
which uses no small-signal code. At the frequencies checked (10, 21.4 and 100 kHz on three bucks), the fork's `.pac` is within 0.044 dB and 0.18° of it. The
stock injected-tone route drifts toward half the switching frequency. On an ideal-switch buck at
100 kHz it was 3.5 dB and 46° off, and its answer moved by 11 dB with the injection amplitude.

**Where the fork is not better:**
- Plain transient on device-dominated circuits gains only 6–11 %.
- Transient accuracy is not claimed to be better than stock. The fork reproduces its upstream
  base exactly, and both differ from 45.2 by about 1 % of ripple on switching decks.
- A circuit clocked only by an XSPICE digital PWM block cannot run in `.pss`: there is no
  source to take the period from.
- **Earlier silent errors in the fork itself.** Builds before `s123` gave wrong `.pac`/`.pstb`
  results on several ideal-switch converters without a message. They were found by checking
  52 decks against the perturbed reference and fixed in `s123`. Since then all 47 open-loop
  decks agree within 0.067 dB and 0.32°, except three documented cases (up to 0.40 dB).

## How it is tested

Each release was run through 46 test suites (stability, periodic analyses, noise,
sensitivities, XSPICE, speed, checkpointing, Monte Carlo and yield tools), also on a build with
AddressSanitizer, plus about 210 checks of documented example results. A fix comes with a
"red proof": its test is run on the previous build, or with the fix switched off, and must
fail there. These suites and the analysis tools built around the fork live in a separate repository
that is not public yet.

Not tested: Windows and macOS builds (Linux only); comparison with a commercial simulator;
an analytic reference for plain transient.

## Building

Linux, same as stock ngspice:

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
