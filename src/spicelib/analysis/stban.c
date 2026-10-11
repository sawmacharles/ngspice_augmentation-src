/**********
ngspice_augmentation private fork -- card #456 (plan T3.2/T3.3), 2026-09-28.

Native `.stb` loop-stability analysis: Tian two-injection return ratio, one
complex-matrix load and factor per frequency, two back-solves (a Bode-return-
ratio analogue of the RFSPICE NIspPreload/NIspSolve split; see acan.c and
niaciter.c). Formula and sign convention: research/cadence-stb-pstb.md sec 2.4
and tools/stbcore.tian() in the ngspice_augmentation repo (Tian, Visvanathan,
Hantgan, Kundert, IEEE Circuits & Devices 17(1), 2001, eq. 30). Probe wiring
(inject at one terminal of an existing 0V V-source, read the other terminal's
voltage and the source's own branch current) follows the same pattern used by
the unmerged 2020 loop_analysis branch's LOOPan() (src/spicelib/analysis/
loopan.c on the ngspice `loop_analysis` branch, never merged upstream) for its
two AC excitations -- this file is an independent, smaller implementation
built for Tian's formula and a single LU factorization per frequency, not a
port of that file's body.

Probe convention (documented per card #456's requirement):
  probe=<Vname>, where Vname is an existing 0V V-source `Vname node_a node_b
  ...` already in the netlist, in series with the wire that breaks the loop.
  node_a (the source's first/positive node) is the "downstream" side (voltage
  is read back here); node_b (the second/negative node) is the "upstream"
  side (unit current is injected here for the current-injection run). Tian's
  T is symmetric under swapping which physical node is called a/b (see
  research/cadence-stb-pstb.md sec 2.4), so this choice does not bias the
  result; it only fixes which physical net tools report the two runs at.
  probe=<inst>:<terminal> (auto-inserting the probe by rebinding a device
  terminal, as the loop_analysis branch's LOOPpreset did) is NOT implemented
  in this fork -- deferred, see the plan T3.2/T3.3 report.

card #457 follow-up (2026-09-28): probe=<inst>:<terminal> IS now implemented
(at netlist-parse time, in dot_stb()/stb_probe_splice(), src/spicelib/parser/
inp2dot.c -- by the time STBan() runs, the synthesized probe is an ordinary
probe=<Vname>, indistinguishable here from a hand-written one). Also added:
native probe=<Vp>,<Vn> mode=dm|cm -- Vp (STBprobeName, "branch p") is read
exactly as the single-ended case; Vn (STBprobeName2, "branch n") is driven
in lockstep with the SAME two RHS injections, sign-flipped for mode=dm
(antiphase) and same-signed for mode=cm (in phase), but never read back.
This mirrors stb/diffprobe.lib's two correlated branches (tools/stb.py
--mode dm|cm), which documents why reading branch p alone gives the whole
coupled network's DM/CM loop gain: for an EXACTLY symmetric linear network,
a purely antisymmetric (DM) or purely symmetric (CM) excitation makes the
whole network's response equally antisymmetric/symmetric, so the ratio T is
the same read at either symmetric point. This is the standard differential/
common-mode half-circuit theorem -- exact for an exactly symmetric linear
circuit, NOT exact in general for an asymmetric network or where Spectre's
own diffprobe would coupled the two branches by DEPENDENT sources (an E/F
pair whose gain is the actual mismatch), which this native version does not
implement (see the card #457 report's "what was not done").
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/stbdefs.h"
#include "ngspice/devdefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/sperror.h"
/* card #475 (plan T10.7): adjoint sensitivity reuses cktsgen.c's generic
 * instance/model-parameter walker (sgen_init/sgen_next/sens_getp), the
 * SAME machinery src/spicelib/analysis/cktsens.c's existing DC/AC `.sens`
 * (a finite-difference-of-the-device's-OWN-stamp, one-solve-per-parameter
 * "direct method") is built on -- see this file's stb_sens_* block below
 * for how the adjoint method replaces its per-parameter SOLVE with a
 * per-parameter cheap DOT PRODUCT against two precomputed adjoint vectors
 * (2 extra solves per frequency total, not per parameter). sens_setp/
 * sens_temp in cktsens.c are file-static there (not exported); this file
 * defines its own tiny equivalents (stb_sens_setp/stb_sens_temp) rather
 * than change cktsens.c's linkage for a private, non-upstream fork. */
#include "ngspice/sensgen.h"
#include "ngspice/stringutil.h"

#include <math.h>

#ifdef XSPICE
#include "ngspice/evt.h"
#include "ngspice/enh.h"
#include "ngspice/mif.h"
#include "ngspice/evtproto.h"
#include "ngspice/ipctiein.h"
#endif

/* card #459 (plan T3.4): device-based ("gain-nulling") .stb, probe=<MOSFET
 * instance>. mos1defs.h resolves via the -I$(top_srcdir)/src/spicelib/
 * devices already on this directory's AM_CPPFLAGS (see Makefile.am) -- no
 * build-system change needed for this one extra include. This is the "small
 * per-model hook" the plan calls for: a helper that reads the instance's own
 * model-specific fields directly, not a new SPICEdev-wide function pointer
 * (SPICEdev is shared by 100+ device tables in this fork; a field there for
 * one device family was not worth the footprint for a private, non-upstream
 * fork). BSIM4 is NOT implemented (see this file's stb_device_probe() and
 * the card #459 report's "what was not done"): its DEVacLoad (b4acld.c)
 * mixes the gm/gmbs stamp with NQS and forward/reverse transcapacitance sums
 * in a way that does not visibly reduce to the same two-vector rank-1
 * pattern without more study than this card's budget allowed. */
#include "mos1/mos1defs.h"
/* card #460: BSIM4 device mode, probe=<BSIM4 instance>. bsim4def.h
 * resolves the same way mos1defs.h does above (this directory's
 * AM_CPPFLAGS already has -I$(top_srcdir)/src/spicelib/devices). */
#include "bsim4/bsim4def.h"

extern SPICEdev **DEVices;

typedef struct { double re, im; } STBcplx;

static STBcplx
stb_cadd(STBcplx a, STBcplx b) { STBcplx r = { a.re + b.re, a.im + b.im }; return r; }

static STBcplx
stb_csub(STBcplx a, STBcplx b) { STBcplx r = { a.re - b.re, a.im - b.im }; return r; }

static STBcplx
stb_cmul(STBcplx a, STBcplx b)
{
    STBcplx r;
    r.re = a.re * b.re - a.im * b.im;
    r.im = a.re * b.im + a.im * b.re;
    return r;
}

static STBcplx
stb_cdiv(STBcplx a, STBcplx b)
{
    STBcplx r;
    double d = b.re * b.re + b.im * b.im;
    r.re = (a.re * b.re + a.im * b.im) / d;
    r.im = (a.im * b.re - a.re * b.im) / d;
    return r;
}

static double
stb_cmag(STBcplx a) { return sqrt(a.re * a.re + a.im * a.im); }

static double
stb_cph(STBcplx a) { return atan2(a.im, a.re); }

/* interpolate the frequency at which a quantity linear in log(f) crosses
 * 'target', between (f0,y0) and (f1,y1); mirrors tools/stbcore._interp_logf */
static double
interp_logf(double f0, double f1, double y0, double y1, double target)
{
    double t;
    if (y1 == y0)
        return f0;
    t = (target - y0) / (y1 - y0);
    return f0 * pow(f1 / f0, t);
}

/* One crossing record, kept in a growable array (mirrors stbcore.margins'
 * pm_crossings / gm_crossings lists). card #475: idx/t2 (the bracket's
 * lower freqs[]/Tarr[]/dTdp_re[k][] index, and the SAME log-f
 * interpolation fraction used to find fc itself) let the post-sweep
 * margin-sensitivity pass re-use the identical bracket/interpolation a
 * crossing was found with, instead of re-deriving it. */
typedef struct { double freq, value, aux, t2; int idx; } CROSSING;

static CROSSING *
push_crossing(CROSSING *arr, int *n, int *cap, double freq, double value, double aux, double t2, int idx)
{
    if (*n >= *cap) {
        *cap = *cap ? *cap * 2 : 8;
        arr = TREALLOC(CROSSING, arr, *cap);
    }
    arr[*n].freq = freq;
    arr[*n].value = value;
    arr[*n].aux = aux;
    arr[*n].t2 = t2;
    arr[*n].idx = idx;
    (*n)++;
    return arr;
}

/* card #459: everything STBan's device-mode branch needs to null a
 * MOSFET's dominant controlled source, in the rank-1 form
 *   T = 1/(1 - s) - 1,  s = gm * v^T Y^-1 u
 * with the ordinary (gm-included) complex AC matrix Y factored ONCE per
 * frequency (the same factorization CKTacLoad/SMPcLUfac already produce for
 * loop mode) and ONE extra back-solve, w = Y^-1 u (research/cadence-stb-
 * pstb.md sec 2.5's "cheaper alternative", and this card's math derivation
 * -- see the card #459 report):
 *   u    = +1 at node_d, -1 at node_s   (the gm/gmbs current: drain row +,
 *          source row -)
 *   v_gm = +1 at node_g, -1 at node_s   (the controlling voltage: gate col
 *          +, source col -)
 * so s = gm * (w[node_g] - w[node_s]), read off the SAME solve vector
 * SMPcSolve already overwrites CKTrhs/CKTirhs with (matching stban.c's
 * existing loop-mode v1/v2/i1/i2 reads, immediately below in STBan()).
 *
 * gm-only vs gm+gmbs (documented per this card's requirement): production
 * nulls gm ALONE -- x = the drain-source VCCS gm*vgs, node_b/gmbs below are
 * computed but NOT added into s. This fork's read of "the dominant
 * controlled source" (research/cadence-stb-pstb.md sec 2.5's [DOC-gap] on
 * what Spectre itself nulls) is the transconductance alone: gmbs is a
 * second, materially smaller VCCS on a typical device and folding it in
 * changes the applicability argument (nulling gm+gmbs together is a rank-1
 * update too, since both share the same output nodes node_d/node_s -- see
 * the report -- but was not made the default). tests/test_stb_device.py's
 * body-effect check shows this choice is NOT negligible: gm-only vs
 * gm+gmbs read measurably different T on a source-degenerated stage. */
typedef struct {
    int node_d, node_s, node_g, node_b;
    double gm, gmbs;
} STBdevProbe;

/* Node convention, worked out by hand against MOS1acLoad's own stamp
 * (src/spicelib/devices/mos1/mos1acld.c -- DPdpPtr/SPspPtr/DPgPtr/SPgPtr/
 * DPbPtr/SPbPtr/DPspPtr/SPdpPtr, the xnrm/xrev-gated terms), not merely
 * assumed: the gm+gmbs contribution to the AC matrix is exactly the outer
 * product u * (gm*v_gm + gmbs*v_gmbs)^T given above, in NORMAL mode
 * (MOS1mode >= 0: node_d/node_s = the internal drain/source-prime nodes,
 * i.e. dNodePrime/sNodePrime -- NOT the external dNode/sNode: when RD/RS
 * != 0 the gm stamp sits on the internal nodes, on the far side of the
 * drain/source series resistance). In REVERSE mode (MOS1mode < 0, Vds < 0:
 * the model's own notion of "drain" and "source" swaps), the identical
 * stamps reappear with dNodePrime and sNodePrime exchanged throughout --
 * this function applies that substitution so callers never need to know
 * the sign of MOS1mode themselves. */
static int
stb_mos1_probe(GENinstance *inst, STBdevProbe *p)
{
    MOS1instance *here = (MOS1instance *) inst;
    bool normal = (here->MOS1mode >= 0);

    p->node_d = normal ? here->MOS1dNodePrime : here->MOS1sNodePrime;
    p->node_s = normal ? here->MOS1sNodePrime : here->MOS1dNodePrime;
    p->node_g = here->MOS1gNode;
    p->node_b = here->MOS1bNode;
    p->gm     = here->MOS1gm;
    p->gmbs   = here->MOS1gmbs;
    return 0;
}

/* card #460 (plan follow-up to #459): BSIM4 device mode, probe=<BSIM4
 * instance>. Derivation (see the card #460 report for the full worked
 * algebra): b4acld.c's AC stamp for gm ALONE (ignoring gmbs, gds, the
 * GIDL/GISL leakage conductances, and the optional gate-leakage/rgate/
 * rbody extra rows -- none of which touch the same four corners, see
 * below) reduces to EXACTLY the same rank-1 pattern as MOS1's, in BOTH
 * BSIM4mode branches:
 *   u    = +1 at dNodePrime, -1 at sNodePrime   (forward, BSIM4mode>=0)
 *          +1 at sNodePrime, -1 at dNodePrime   (reverse, BSIM4mode<0)
 *   v_gm = +1 at gNodePrime, -1 at (whichever of the above is "node_s")
 *   scalar = +gm  (BSIM4gm, UNCHANGED sign in both modes -- verified by
 *          hand against b4acld.c's Gmr/FwdSumr/RevSumr bookkeeping: the
 *          "node_d/node_s swap" picture below and BSIM4's own reverse-
 *          mode Gmr=-gmr-plus-unswapped-DP/SP-rows bookkeeping are two
 *          different-looking but ALGEBRAICALLY IDENTICAL ways of writing
 *          the same six nonzero matrix entries -- checked corner by
 *          corner, not assumed)
 * i.e. this function is stb_mos1_probe() with the field names changed;
 * the node_d/node_s swap-on-reverse trick reused here for BSIM4 gives
 * bit-for-bit the same 6 corner entries as b4acld.c's own reverse-mode
 * Gmr=-gmr bookkeeping produces, by direct hand expansion of both.
 *
 * Why the OTHER stamps in b4acld.c don't corrupt this: gds (channel
 * conductance, incl. its own NQS real/imag split) is a plain 2-terminal
 * admittance between dNodePrime/sNodePrime -- passive, stays in Y, not
 * nulled (exactly as MOS1's Rds is never nulled either). GIDL/GISL
 * (the ggidl-prefixed and ggisl-prefixed fields) and gate-leakage (the
 * gIxtot-prefixed fields) are separate, genuinely
 * small leakage-current stamps with their OWN control pattern (not a
 * scalar multiple of v_gm) -- they are not "the dominant controlled
 * source" in Tian's sense (nulling them would not render the device
 * passive; they are leakage, not gain) and are left in Y, same
 * treatment as gds. rdsMod/rgateMod/rbodyMod add entries at OTHER matrix
 * corners entirely (Dgp/Sgp/GMgp/DBdb/etc, never DPgp/DPsp/SPgp/SPsp/
 * DPdp/SPdp), so they don't interact with this derivation regardless of
 * their setting.
 *
 * Scope, deliberately narrow (card #460's own ask -- "quasi-static, NQS
 * off, first"): QUASI-STATIC ONLY. Under here->BSIM4acnqsMod, b4acld.c
 * turns Gm into a COMPLEX, frequency-shaped scalar gm*(T2-jT3) (T2,T3
 * from omega*BSIM4taunet) -- itself STILL a rank-1 update with the SAME
 * u/v (a per-frequency complex scalar in place of the real gm), so NQS
 * support is a small, mechanical follow-up (recompute T2/T3 from
 * ckt->CKTomega and BSIM4taunet here, same formula as b4acld.c) -- NOT
 * implemented in this card; instances with acnqsMod on are refused
 * (return nonzero) rather than silently nulling the wrong (real-only)
 * scalar. gmbs is NOT folded in (same "gm-only" production default as
 * MOS1, for the same reason -- see stb_mos1_probe()'s doc comment
 * above); it reduces to rank-1 too (same u, v_gmbs = +1@bNodePrime,
 * -1@node_s) but was not made the default here either. */
static int
stb_bsim4_probe(GENinstance *inst, STBdevProbe *p)
{
    BSIM4instance *here = (BSIM4instance *) inst;
    bool normal = (here->BSIM4mode >= 0);

    if (here->BSIM4acnqsMod)
        return 1;

    p->node_d = normal ? here->BSIM4dNodePrime : here->BSIM4sNodePrime;
    p->node_s = normal ? here->BSIM4sNodePrime : here->BSIM4dNodePrime;
    p->node_g = here->BSIM4gNodePrime;
    p->node_b = here->BSIM4bNodePrime;
    p->gm     = here->BSIM4gm;
    p->gmbs   = here->BSIM4gmbs;
    return 0;
}

/* ============================================================
 * card #475 (plan T10.7): adjoint sensitivity of the Tian loop-gain
 * probe readings (v1,i1,v2,i2 -> T -> margins) wrt element/model
 * parameters. Design: research/feature-survey.md's adjoint-sensitivity
 * item, docs/stb.md, and the card's own DESIGN note (see the report).
 *
 * For a small-signal quantity y = c^T x at one frequency, with
 * Y(p) x = b (b independent of p here -- the probe's own unit
 * injections), the ordinary forward/direct sensitivity is
 *   dx/dp = Y^-1 (db/dp - (dY/dp) x)
 * (what cktsens.c's existing DC/AC `.sens` computes, ONE extra solve
 * per parameter, reusing the already-factored Y). The ADJOINT form
 * reuses the same identity by moving the solve onto c instead of onto
 * dY/dp x:
 *   dy/dp = lambda^T (db/dp - (dY/dp) x),   Y^T lambda = c
 * -- ONE extra (TRANSPOSED, same factorisation) solve per OUTPUT
 * reading, independent of how many parameters p are swept, then a
 * cheap sparse matrix-vector product + dot product per parameter (no
 * solve at all). Four probe readings (v1,i1,v2,i2) but only two
 * DISTINCT c vectors (v1/v2 share c_v = e_{node_b}; i1/i2 share
 * c_i = -e_{branch_probe}), so exactly 2 adjoint solves per frequency
 * cover all four.
 *
 * dY/dp (times a small delta) for one isolated device/parameter is
 * built with EXACTLY cktsens.c's own technique (DEVsetup into a fresh,
 * per-parameter-isolated delta_Y, then DEVacLoad at p and at p+delta,
 * subtracting) -- this is not a fresh design, it is the same proven
 * mechanism reused so every IF_REAL, AC-capable instance/model
 * parameter of every device type (R, C, L, VCCS/VCVS gains, and any
 * MOSFET model parameter, e.g. w/l/vto) is already covered generically,
 * with no per-device-type code. What this does NOT cover: MOSFET small-
 * signal "device" quantities (gm/gds/caps) treated as independent knobs
 * -- those are operating-point-DERIVED, not settable instance/model
 * parameters, so they are outside this generic mechanism; see the card
 * #475 report's "what was not done".
 *
 * Loop mode only (device mode / dm|cm are refused together with sens=,
 * see STBan() below) and SPARSE (non-KLU) only: SMPcaSolve (the complex
 * TRANSPOSED solve this needs) and SMPmultiply only implement the
 * Sparse1.3 matrix path in this fork (spsolve.c/sputils.c) -- refused,
 * not silently wrong, when CKTkluMODE is on (see below). */

static double STBsens_Delta = 1e-6;      /* same defaults as cktsens.c's Sens_Delta */
static double STBsens_AbsDelta = 1e-6;

/* Simple glob match ('*' / '?'), mirroring cktsens.c's static scan()/
 * check_filter() (duplicated here rather than exported from cktsens.c --
 * both are private-fork WIP files and this keeps the two independent). */
static int
stb_sens_glob(const char *filter, const char *name)
{
    while (*filter && *name) {
        if (*filter == '*') {
            if (!filter[1])
                return 1;
            while (*name && !stb_sens_glob(filter + 1, name))
                ++name;
            return *name != '\0';
        }
        if (*filter == *name || *filter == '?')
            ++name, ++filter;
        else
            return 0;
    }
    return !*name && (!*filter || (*filter == '*' && !filter[1]));
}

/* namebuf convention, byte-for-byte mirroring cktsens.c's sens_sens() so
 * sens=<name> selects the SAME parameter a .sens run would print (minus
 * cktsens.c's "M1_vto" underscore form for a non-principal instance
 * param, which is kept here too -- only the "M1:vto" model-param colon
 * form is remapped to underscore for the OUTPUT VECTOR name, since ':'
 * is not safe in a vector name; the sens=<name> MATCH itself still uses
 * the colon form to stay identical to cktsens.c's own naming). */
static void
stb_sens_name(sgen *sg, char *namebuf, size_t buflen)
{
    if (!sg->is_instparam) {
        snprintf(namebuf, buflen, "%s:%s", sg->instance->GENname, sg->ptable[sg->param].keyword);
    } else if ((sg->ptable[sg->param].dataType & IF_PRINCIPAL) && sg->is_principle == 1) {
        snprintf(namebuf, buflen, "%s", sg->instance->GENname);
    } else {
        snprintf(namebuf, buflen, "%s_%s", sg->instance->GENname, sg->ptable[sg->param].keyword);
    }
}

static int
stb_sens_selected(STBAN *job, const char *namebuf)
{
    int k;
    if (job->STBsensAll)
        return 1;
    for (k = 0; k < job->STBsensCount; k++)
        if (stb_sens_glob((char *) job->STBsensNames[k], namebuf))
            return 1;
    return 0;
}

/* AC-only equivalents of cktsens.c's (file-static, unexported) sens_load/
 * sens_temp -- same one-line bodies, just not reachable from here. */
static int
stb_sens_load_ac(sgen *sg, CKTcircuit *ckt)
{
    int (*fn)(GENmodel *, CKTcircuit *) = DEVices[sg->dev]->DEVacLoad;
    if (!fn) return 1;
    return fn(sg->model, ckt);
}

static int
stb_sens_temp(sgen *sg, CKTcircuit *ckt)
{
    int (*fn)(GENmodel *, CKTcircuit *) = DEVices[sg->dev]->DEVtemperature;
    if (!fn) return 1;
    return fn(sg->model, ckt);
}

/* cktsens.c's sens_setp is file-static (see this block's top comment) --
 * same body, new name. */
static int
stb_sens_setp(sgen *sg, CKTcircuit *ckt, IFvalue *val)
{
    int pid, error = 0;
    NG_IGNORE(ckt);
    if (sg->is_instparam) {
        int (*fn)(int, IFvalue *, GENinstance *, IFvalue *) = DEVices[sg->dev]->DEVparam;
        pid = DEVices[sg->dev]->DEVpublic.instanceParms[sg->param].id;
        if (fn) error = fn(pid, val, sg->instance, NULL);
        else return 1;
    } else {
        int (*fn)(int, IFvalue *, GENmodel *) = DEVices[sg->dev]->DEVmodParam;
        pid = DEVices[sg->dev]->DEVpublic.modelParms[sg->param].id;
        if (fn) error = fn(pid, val, sg->model);
        else return 1;
    }
    return error;
}

/* complex "bilinear" dot product lambda^T v (ordinary transpose, NOT a
 * conjugated inner product -- see this block's top comment: the adjoint
 * identity needs Y^T, not Y^H). Skips index 0 (the MNA ground "trash"
 * slot -- same guard stban.c's loop-mode code already applies to
 * CKTrhs[0]/CKTirhs[0], see STBan() below). */
static STBcplx
stb_sens_dot(const double *lam_re, const double *lam_im,
             const double *v_re, const double *v_im, int size)
{
    STBcplx r = { 0.0, 0.0 };
    int j;
    for (j = 1; j <= size; j++) {
        r.re += lam_re[j] * v_re[j] - lam_im[j] * v_im[j];
        r.im += lam_re[j] * v_im[j] + lam_im[j] * v_re[j];
    }
    return r;
}

/* card #595: the bias-point term of dT/dp.
 *
 * Card #475 held the DC operating point fixed: dY/dp was only the
 * device's own stamp derivative at constant node voltages. For a
 * parameter that moves the bias (a MOSFET's w, a bias resistor, a
 * model parameter) the small-signal matrix also changes through the
 * operating point, dY/dp|total = dY/dp|x + (dY/dx)(dx_op/dp), and the
 * second term was missing (test_sens.py's M2_w case was labelled
 * EXPECT-DIVERGE for that reason).
 *
 * dx_op/dp: one DC Newton step from the converged x_op with the
 * parameter at p+delta and at p-delta (central difference). At the
 * converged point the Newton map N(p) = J(x_op,p)^-1 rhs(x_op,p) has
 * dN/dp = -J^-1 dF/dp = dx_op/dp (the residual is zero to Newton's own
 * accuracy, and its contribution cancels in the central difference).
 * Two real factorisations per parameter, once, before the sweep.
 *
 * (dY/dx)(dx_op/dp) x, per frequency: the circuit is re-linearised at
 * x_op + h*dx_op/dp and at x_op - h*dx_op/dp (CKTload MODEINITFLOAT to
 * put the shifted voltages into the device states, MODEINITSMSIG to
 * recompute the small-signal quantities from them, CKTacLoad), the
 * matrix multiplied by the forward solutions, and the two products
 * differenced. h is chosen so the largest node-voltage shift is
 * STB_SENS_BIAS_DV (default 1e-5 V). Afterwards CKTstate0 is restored
 * from a copy and MODEINITSMSIG re-run, so the next frequency's
 * linearisation (and T itself) is the same as without sens=.
 *
 * The same pass also fixes the fixed-op term of nonlinear devices: the
 * card #475 stamp difference re-ran only DEVacLoad at p+delta, but a
 * nonlinear device computes gm, gds and its charges in DEVload, so their
 * dependence on p (a MOSFET's gm on w, a diode's capacitance on cjo) was
 * missing as well. For every parameter of a device that is not one of the
 * linear elements (R, C, L, K, E, F, G, H, V, I) the whole derivative is
 * therefore taken jointly, Y(p + h, x_op + h*dx_op/dp) against
 * Y(p - h, x_op - h*dx_op/dp), and replaces the stamp difference; linear
 * elements keep the exact stamp difference and add the bias term only if
 * dx_op/dp is not zero.
 *
 * Parameters for which dx_op/dp is exactly zero on a linear element (a
 * capacitor, an AC-only resistor) skip the per-frequency work.
 *
 * A MODEL parameter of a model with more than one instance keeps card
 * #475's stamp difference (one warning lists them; the result is
 * incomplete): that term is per instance (cktsgen.c isolates one
 * instance), but a model parameter set outside that isolation, as a full
 * CKTload needs, moves every instance.
 * Instant mode (no DC solve) keeps the fixed-op term with a warning.
 * STB_SENS_FIXED_OP=1 restores card #475's behaviour (red proof). */
typedef struct {
    int dev, param, is_instparam;
    GENmodel *model;
    GENinstance *instance;
    double value;
} STBsensRec;

/* the linear elements, whose AC stamp is their whole small-signal
 * dependence on their own parameters (see the block comment above) */
static int
stb_sens_linear_dev(int dev)
{
    static const char *lin[] = { "Resistor", "Capacitor", "Inductor", "mutual", "VCVS", "VCCS",
                                 "CCCS", "CCVS", "Vsource", "Isource", NULL };
    const char *n = DEVices[dev] ? DEVices[dev]->DEVpublic.name : NULL;
    int i;
    if (!n) return 0;
    for (i = 0; lin[i]; i++)
        if (strcmp(n, lin[i]) == 0) return 1;
    return 0;
}

static int
stb_rec_setp(STBsensRec *r, double v, CKTcircuit *ckt)
{
    IFvalue val;
    int pid, error;
    val.rValue = v;
    if (r->is_instparam) {
        int (*fn)(int, IFvalue *, GENinstance *, IFvalue *) = DEVices[r->dev]->DEVparam;
        if (!fn) return 1;
        pid = DEVices[r->dev]->DEVpublic.instanceParms[r->param].id;
        error = fn(pid, &val, r->instance, NULL);
    } else {
        int (*fn)(int, IFvalue *, GENmodel *) = DEVices[r->dev]->DEVmodParam;
        if (!fn) return 1;
        pid = DEVices[r->dev]->DEVpublic.modelParms[r->param].id;
        error = fn(pid, &val, r->model);
    }
    if (DEVices[r->dev]->DEVtemperature)
        DEVices[r->dev]->DEVtemperature(r->model, ckt);
    return error;
}

/* one DC Newton step from xop: out = J(xop)^-1 rhs(xop) at the current
 * parameter values. The caller saves and restores CKTrhs/CKTrhsOld/
 * CKTmode and the device states. */
static int
stb_sens_newton_step(CKTcircuit *ckt, const double *xop, double *bufOld, double *out, int size)
{
    int error;
    memcpy(bufOld, xop, (size_t) (size + 1) * sizeof(double));
    ckt->CKTrhsOld = bufOld;
    ckt->CKTrhs = out;
    ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODEDCOP | MODEINITFLOAT;
    error = CKTload(ckt);
    if (error) return error;
    error = SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin);
    if (error == E_SINGULAR)
        error = SMPreorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTpivotRelTol, ckt->CKTdiagGmin);
    if (error) return error;
    SMPsolve(ckt->CKTmatrix, out, ckt->CKTrhsSpare);
    out[0] = 0.0;
    return OK;
}

/* the complex matrix linearised at xop + h*dx (see the block comment
 * above); leaves it in ckt->CKTmatrix, unfactored */
static int
stb_sens_shift_load(CKTcircuit *ckt, const double *xop, const double *dx, double h,
                    double *bufOld, double *scr_re, double *scr_im, int size)
{
    int j, error;
    long savemode = ckt->CKTmode;
    for (j = 0; j <= size; j++)
        bufOld[j] = dx ? xop[j] + h * dx[j] : xop[j];
    bufOld[0] = 0.0;
    ckt->CKTrhsOld = bufOld;
    ckt->CKTrhs = scr_re;
    ckt->CKTirhs = scr_im;
    ckt->CKTmode = (savemode & MODEUIC) | MODEDCOP | MODEINITFLOAT;
    error = CKTload(ckt);
    if (error) return error;
    ckt->CKTmode = (savemode & MODEUIC) | MODEDCOP | MODEINITSMSIG;
    error = CKTload(ckt);
    if (error) return error;
    ckt->CKTmode = (savemode & MODEUIC) | MODEAC;
    error = CKTacLoad(ckt);
    ckt->CKTmode = savemode;
    return error;
}

int
STBan(CKTcircuit *ckt, int restart)
{
    STBAN *job = (STBAN *) ckt->CKTcurJob;

    double freq, freqTol;
    int error, size, i, ignore;
    GENinstance *inst;
    int branch_probe, node_a, node_b;
    /* card #457: probe=<Vp>,<Vn> mode=dm|cm -- branch n, driven in lockstep
     * with branch p (this file's existing branch_probe/node_a/node_b) but
     * never read back; see this file's top-of-file comment. */
    GENinstance *inst2 = NULL;
    int branch_probe2 = 0, node_b2 = 0;
    double dmcm_sign = 1.0;

    /* card #459: device-based (gain-nulling) mode, probe=<MOSFET instance>.
     * See stb_mos1_probe() above and the card #459 report for the math. */
    int device_mode = 0;
    STBdevProbe devp;

    STBcplx v1, i1, v2, i2, s, t, lg;

    static runDesc *stbPlot = NULL;
    IFuid freqUid;
    IFuid vecuids[2];

    double *freqs = NULL;
    STBcplx *Tarr = NULL, *LGarr = NULL;
    int npts = 0, cap = 0;
    /* card #475: kept alive (not freed right after the margin/crossing
     * pass, unlike upstream's original scoping) so the dPM/dp, dfc/dp
     * block below can reuse the same unwrapped-phase trace the crossing
     * search itself used, rather than recomputing it. */
    double *phs = NULL;

    CROSSING *pmx = NULL, *gmx = NULL;
    int npm = 0, ngm = 0, capm = 0, capg = 0;

    /* card #475: adjoint sensitivity state. sens_requested is fixed at
     * entry (job->STBsensAll/STBsensCount, set by the parser -- see
     * stbdefs.h/stbsetp.c); nsens (the actual matched-parameter count) is
     * only known after the first dry enumeration pass, below, once the DC
     * operating point exists. sensNames/sensVecUids/dTdp_re/dTdp_im are
     * parallel nsens-length arrays; dTdp_re[k]/dTdp_im[k] grow in lockstep
     * with freqs[]/Tarr[] (same cap-doubling pattern) so a per-parameter
     * dT/dp trace can be read back at the same bracket indices the
     * existing PM/GM crossing search already found. */
    int sens_requested = (job->STBsensAll || job->STBsensCount > 0);
    int nsens = 0;
    char **sensNames = NULL;
    IFuid *sensVecUids = NULL;
    double **dTdp_re = NULL, **dTdp_im = NULL;
    SMPmatrix *sensDeltaY = NULL;
    double *sens_lamv_re = NULL, *sens_lamv_im = NULL;
    double *sens_lami_re = NULL, *sens_lami_im = NULL;
    double *sens_x1_re = NULL, *sens_x1_im = NULL;
    double *sens_x2_re = NULL, *sens_x2_im = NULL;
    double *sens_dI = NULL, *sens_dII = NULL;
    double *sens_mv_re = NULL, *sens_mv_im = NULL;
    runDesc *sensPlot = NULL;
    /* card #595: bias-point term state (see stb_rec_setp() above) */
    STBsensRec *sensRec = NULL;
    double **sensDxdp = NULL, *sensH = NULL;
    double *sens_xop = NULL, *sens_bufOld = NULL, *sens_state0 = NULL;
    double *sens_yp1_re = NULL, *sens_yp1_im = NULL, *sens_yp2_re = NULL, *sens_yp2_im = NULL;
    STBcplx *sensD = NULL;      /* nsens x 4: dv1, di1, dv2, di2 at this frequency */
    int *sensNoAc = NULL, *sensJoint = NULL;
    int sens_nbias = 0;

    /* ---- resolve the probe device ---- */
    inst = CKTfndDev(ckt, job->STBprobeName);
    if (!inst || inst->GENmodPtr->GENmodType < 0) {
        SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe '%s' not in circuit", job->STBprobeName);
        return E_NOTFOUND;
    }
    if (inst->GENmodPtr->GENmodType == CKTtypelook("Mos1") ||
        inst->GENmodPtr->GENmodType == CKTtypelook("BSIM4")) {
        /* card #459/#460: device mode (MOS1 or, card #460, BSIM4).
         * probe2/mode=dm|cm (loop-mode DM/CM, card #457) is meaningless
         * here -- refuse rather than silently ignore it (rule: no silent
         * numbers). */
        device_mode = 1;
        if (job->STBprobeName2Given || job->STBmodeGiven) {
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb probe '%s' is a device (device-mode gain-nulling); "
                "probe2=/mode=dm|cm (loop-mode differential probing) cannot "
                "be combined with a device-mode probe", job->STBprobeName);
            return E_PARMVAL;
        }
        if (job->STBlocalGndGiven)
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb probe '%s': localgnd= has no effect in device mode "
                "(there is no current-injection node to reference) -- ignored",
                job->STBprobeName);
        /* NOTE: devp (node indices AND gm/gmbs) is filled in below, AFTER
         * the DC operating point is computed -- NOT here. MOS1gm/MOS1gmbs
         * are only valid once CKTop()+CKTload() have run for this job; a
         * fresh (never-op'd) circuit reads back gm=gmbs=0 here (found
         * empirically while validating this card: a device-mode run with
         * no PRECEDING top-level "op" command silently read a stale
         * gm=0, giving T=0 for every frequency -- see the card #459
         * report's "what was not done" / bugs found). */
        branch_probe = 0; node_a = 0; node_b = 0;  /* unused in device mode */
    } else if (inst->GENmodPtr->GENmodType == CKTtypelook("Vsource")) {
        branch_probe = CKTfndBranch(ckt, job->STBprobeName);
        if (branch_probe == 0) {
            SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe '%s': no branch equation found", job->STBprobeName);
            return E_NOTFOUND;
        }
        node_a = GENnode(inst)[0];
        node_b = GENnode(inst)[1];
        /* card #460: loop-mode probe with EITHER terminal at literal
         * ground (node 0) is degenerate for Tian's two-injection method,
         * for a reason that has nothing to do with this file's own
         * node-0-trash-slot hazard (guarded below): the probe's own 0V
         * branch equation is V(node_a) - V(node_b) = <this run's
         * injection>. If node_a==0 (true ground, a fixed 0, not a row of
         * its own), that equation ALONE algebraically forces node_b to a
         * fixed value equal to -<injection>, for EVERY circuit,
         * regardless of anything attached elsewhere -- so run 1's
         * "voltage read back at node_b" is not a network response at
         * all, it is the injection itself, unchanged. If node_b==0, the
         * symmetric thing happens to run 2: injecting current AT literal
         * ground has nowhere to go but through the probe's own branch
         * (ground has no KCL row to receive it), so node_a's response
         * -- which this file does not even read, see the loop below --
         * or a naive read of node_b (=ground, always exactly 0 by
         * definition) carries no information either. Proved three ways
         * (algebraically from the branch equation above; empirically
         * from this build with STB_DEBUG=1 on tests/fixtures/
         * stb3_gnd_loop_linear.cir and stb3_gnd_loop_mos.cir, both
         * orientations; and independently via tools/stb.py's real
         * 3-element lgprobe subckt, which has its own separate internal
         * "x" node and STILL comes out identically degenerate for the
         * same algebraic reason) -- see the card #460 report. Both
         * degenerate cases converge on T=0 exactly once the node-0
         * read/write hazard below is also guarded (before that guard,
         * node_b==0 additionally read back STALE data and gave a small
         * WRONG nonzero T, not the honest degenerate 0 -- see the FAIL
         * before that guard is added below). Warn rather than silently
         * hand back a T that looks like "provably stable, huge margin"
         * when it is actually "this measurement carries no information"
         * (rule: no silent numbers) -- but still compute and report the
         * (honest, degenerate) T=0 rather than refuse outright: unlike
         * the device-mode probe2/mode= conflict above, there is a
         * single well-defined right answer here (T=0), not an
         * ill-posed request. */
        if (node_a == 0 || node_b == 0)
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb probe '%s': one terminal is literal ground (node 0) "
                "-- Tian's two-injection loop-based method cannot measure "
                "a loop broken there (the probe's own branch equation "
                "pins the OTHER terminal to the injection value itself, "
                "or the injected current has nowhere to go but the "
                "probe's own branch): loopgain/T will read exactly 0 at "
                "every frequency, which means 'not measured', not 'stable "
                "with a huge margin'. Break the loop one component "
                "further from true ground.", job->STBprobeName);
    } else {
        SPfrontEnd->IFerrorf(ERR_WARNING,
            "stb probe '%s' must be a voltage source (an existing 0V V-source "
            "in series with the loop-break wire) or a supported active device "
            "instance (device mode: mos1 or bsim4, quasi-static only, in "
            "this fork)", job->STBprobeName);
        return E_NOTFOUND;
    }

    /* ---- card #457: resolve the optional branch-n probe (probe=Vp,Vn
     * mode=dm|cm) ---- */
    if (job->STBprobeName2Given) {
        if (!job->STBmodeGiven) {
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb probe2 '%s' given without mode=dm|cm", job->STBprobeName2);
            return E_PARMVAL;
        }
        inst2 = CKTfndDev(ckt, job->STBprobeName2);
        if (!inst2 || inst2->GENmodPtr->GENmodType < 0) {
            SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe2 '%s' not in circuit", job->STBprobeName2);
            return E_NOTFOUND;
        }
        if (inst2->GENmodPtr->GENmodType != CKTtypelook("Vsource")) {
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb probe2 '%s' must be a voltage source (an existing 0V V-source "
                "in series with the loop-break wire)", job->STBprobeName2);
            return E_NOTFOUND;
        }
        branch_probe2 = CKTfndBranch(ckt, job->STBprobeName2);
        if (branch_probe2 == 0) {
            SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe2 '%s': no branch equation found", job->STBprobeName2);
            return E_NOTFOUND;
        }
        node_b2 = GENnode(inst2)[1];
        /* dm = antiphase (branch n's excitation is the negative of branch
         * p's, same as stb/diffprobe.lib's diffprobe_dm baking acphase=180
         * on its second branch); cm = in phase (diffprobe_cm, acphase=0). */
        dmcm_sign = (job->STBmode == STB_MODE_DM) ? -1.0 : 1.0;
        /* card #460: branch n's own node_b2==0 hits the same current-
         * injection write hazard as branch p's node_b (guarded below) --
         * warn the same way (branch n is never read back, only driven in
         * lockstep, but a degenerate branch n still perturbs the shared
         * two-injection solve if left unguarded). */
        if (GENnode(inst2)[0] == 0 || node_b2 == 0)
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb probe2 '%s': one terminal is literal ground (node 0) "
                "-- see the probe= warning above; the same degeneracy "
                "applies to the correlated branch", job->STBprobeName2);
    } else if (job->STBmodeGiven) {
        SPfrontEnd->IFerrorf(ERR_WARNING, "stb mode=dm|cm given without probe=<Vp>,<Vn>");
        return E_PARMVAL;
    }

    /* card #475: sens= scope checks -- loop mode (single-ended), Sparse
     * only. Refused up front (E_PARMVAL), never silently ignored or
     * silently wrong (rule: no silent numbers). */
    if (sens_requested) {
        if (device_mode) {
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb sens=: device-mode probes are not supported (loop mode only)");
            return E_PARMVAL;
        }
        if (job->STBprobeName2Given || job->STBmodeGiven) {
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb sens=: probe2=/mode=dm|cm (differential probing) is not supported "
                "(single-ended loop mode only)");
            return E_PARMVAL;
        }
#ifdef KLU
        if (ckt->CKTmatrix->CKTkluMODE) {
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb sens=: needs the Sparse solver (SMPcaSolve/SMPmultiply are not "
                "implemented for KLU in this fork) -- add '.options sparse' or remove "
                "'.option klu' and re-run");
            return E_PARMVAL;
        }
#endif
    }

    /* ---- DC operating point (identical to the closed loop: the probe is
     * an ordinary 0V branch already in the netlist, not a mutation) ---- */
    if (job->STBsaveFreq == 0 || restart) {

        if (job->STBnumberSteps < 1)
            job->STBnumberSteps = 1;

        switch (job->STBstepType) {
        case DECADE:
            if (job->STBstartFreq <= 0) {
                fprintf(stderr, "ERROR: STB startfreq <= 0\n");
                return E_PARMVAL;
            }
            job->STBfreqDelta = exp(log(10.0) / job->STBnumberSteps);
            break;
        case OCTAVE:
            if (job->STBstartFreq <= 0) {
                fprintf(stderr, "ERROR: STB startfreq <= 0\n");
                return E_PARMVAL;
            }
            job->STBfreqDelta = exp(log(2.0) / job->STBnumberSteps);
            break;
        case LINEAR:
            if (job->STBnumberSteps - 1 > 1)
                job->STBfreqDelta = (job->STBstopFreq - job->STBstartFreq) / (job->STBnumberSteps - 1);
            else
                job->STBfreqDelta = 0;
            break;
        default:
            return E_BADPARM;
        }

        /* card #500 (instant mode): SKIP the DC operating-point solve
         * entirely and reuse whatever state is already loaded in `ckt`.
         * The normal (default) path below calls CKTop()/EVTop() -- full
         * Newton -- which ALWAYS converges to the circuit's true DC
         * equilibrium given the present independent-source values,
         * discarding any prior .tran/IC=/.nodeset/savestate state (this is
         * exactly the limitation docs/stb-tran.md documents in detail,
         * verified three independent ways). "instant" exists specifically
         * so a `.control` sequence can do
         *   tran <step> {t_k} 0 uic
         *   stb dec ... probe=Vprobe instant
         * and have the second command linearise about the FIRST command's
         * own accepted transient solution at t=t_k (already sitting in
         * CKTrhsOld and in every device's own persistent small-signal
         * fields, e.g. MOS1gm/BSIM4gm/a diode's geq, set by the .tran's
         * own last DEVload call at t_k) instead of re-solving to
         * equilibrium. This is architecturally sound because CKTop's
         * Newton loop is a fixed-point SEARCH (it always converges to the
         * SAME true equilibrium regardless of the starting guess, barring
         * genuine bistability) -- skipping the search is the only way to
         * linearise about a point Newton would otherwise walk straight
         * past. The CKTload(MODEINITSMSIG) call just below still runs
         * (unconditionally, both paths): it is a SINGLE evaluation of the
         * device load routines at whatever node voltages are already
         * present, not an iteration, so calling it at the frozen t_k
         * voltages is exactly as valid as calling it at a converged DC
         * equilibrium -- the per-device algebraic quantities it computes
         * (gm, gds, an incremental capacitance) are pure functions of
         * terminal voltages, not of how those voltages were reached. What
         * the load computes UNDER MODEDCOP (no transient charge/current
         * integration -- reactive elements are treated as open/short and
         * folded into the matrix only via CKTacLoad's separate jwC/jwL
         * stamp, below) is unaffected by which analysis last touched the
         * circuit. The one thing instant mode does NOT get right: at a
         * genuine transient (non-equilibrium) point, this MODEDCOP-style
         * reload's OWN scratch RHS residual (computed from the frozen
         * voltages, ignoring capacitor currents) is not zero -- but
         * that residual is simply overwritten by CKTacLoad/the
         * probe-injection RHS a few lines below and never used: AC
         * analysis is a pure LINEARISATION (a single complex-matrix
         * solve of Y*x=b for the probe excitation), not a Newton search,
         * so a nonzero "DC" residual at a non-equilibrium point is
         * irrelevant to it -- see docs/stb.md's "instant mode" section
         * for the full KCL-residual argument.
         *
         * XSPICE event-driven instances are refused outright (not
         * silently wrong): EVTop() (the event-system equivalent of
         * CKTop, above) also re-solves to a fresh digital/event
         * equilibrium, and this fork's instant-mode patch does not
         * attempt an event-side equivalent of the CKTop skip -- an event
         * netlist with instant would otherwise appear to work but
         * silently mix a frozen analog state with a re-initialised
         * event state. */
        if (job->STBinstant) {
#ifdef XSPICE
            if (ckt->evt && ckt->evt->counts.num_insts != 0) {
                SPfrontEnd->IFerrorf(ERR_WARNING,
                    "stb instant: XSPICE event-driven instances are present -- "
                    "instant mode (card #500) only skips the analog CKTop() "
                    "solve, not EVTop()'s event-side equilibrium re-solve, so "
                    "the result would silently mix a frozen analog state with "
                    "a re-initialised event state. Not supported; refusing.");
                return E_PARMVAL;
            }
#endif
            fprintf(stdout,
                "\n stb instant: skipping the op-point solve, linearising "
                "about the circuit's CURRENT state (e.g. a preceding "
                "'tran ... uic')\n");
            if (getenv("STB_DEBUG_INSTANT")) {
                int dbgsize = SMPmatSize(ckt->CKTmatrix);
                int dbgi;
                fprintf(stderr, "STB_DEBUG_INSTANT CKTtime=%.9g CKTmode=0x%x size=%d\n",
                        ckt->CKTtime, ckt->CKTmode, dbgsize);
                for (dbgi = 1; dbgi <= dbgsize; dbgi++)
                    fprintf(stderr, "  CKTrhsOld[%d]=%.9g CKTrhs[%d]=%.9g\n",
                            dbgi, ckt->CKTrhsOld[dbgi], dbgi, ckt->CKTrhs[dbgi]);
            }
        } else {
#ifdef XSPICE
        g_mif_info.circuit.anal_type = MIF_DC;
        g_mif_info.circuit.anal_init = MIF_TRUE;
        if (ckt->evt->counts.num_insts != 0) {
            error = EVTop(ckt,
                (ckt->CKTmode & MODEUIC) | MODEDCOP | MODEINITJCT,
                (ckt->CKTmode & MODEUIC) | MODEDCOP | MODEINITFLOAT,
                ckt->CKTdcMaxIter, MIF_TRUE);
            EVTdump(ckt, IPC_ANAL_DCOP, 0.0);
            EVTop_save(ckt, MIF_TRUE, 0.0);
        } else
#endif
        if (!ckt->CKTnoopac) {
            error = CKTop(ckt,
                (ckt->CKTmode & MODEUIC) | MODEDCOP | MODEINITJCT,
                (ckt->CKTmode & MODEUIC) | MODEDCOP | MODEINITFLOAT,
                ckt->CKTdcMaxIter);
            if (error) {
                fprintf(stdout, "\nSTB operating point failed -\n");
                CKTncDump(ckt);
                return error;
            }
        } else {
            fprintf(stdout, "\n Linear circuit, option noopac given: no OP analysis\n");
        }
        }

        ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODEDCOP | MODEINITSMSIG;
        error = CKTload(ckt);
        if (error) return error;

#ifdef KLU
        if (ckt->CKTmatrix->CKTkluMODE) {
            if (!ckt->CKTmatrix->SMPkluMatrix->KLUmatrixIsComplex) {
                for (i = 0; i < DEVmaxnum; i++)
                    if (DEVices[i] && DEVices[i]->DEVbindCSCComplex && ckt->CKThead[i])
                        DEVices[i]->DEVbindCSCComplex(ckt->CKThead[i], ckt);
                ckt->CKTmatrix->SMPkluMatrix->KLUmatrixIsComplex = KLUMatrixComplex;
            }
        }
#endif

        {
            IFuid stbuid_lg, stbuid_t;
            SPfrontEnd->IFnewUid(ckt, &stbuid_lg, NULL, "loopgain", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &stbuid_t, NULL, "T", UID_OTHER, NULL);
            vecuids[0] = stbuid_lg;
            vecuids[1] = stbuid_t;
        }
        SPfrontEnd->IFnewUid(ckt, &freqUid, NULL, "frequency", UID_OTHER, NULL);
        error = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob,
            ckt->CKTcurJob->JOBname,
            freqUid, IF_REAL,
            2, vecuids, IF_COMPLEX,
            &stbPlot);
        if (error) return error;
        if (job->STBstepType != LINEAR)
            SPfrontEnd->OUTattributes(stbPlot, NULL, OUT_SCALE_LOG, NULL);

        /* card #475: dry-run parameter enumeration -- walks EVERY IF_REAL,
         * AC-capable instance/model parameter of every device (same
         * generic walk cktsens.c's DC/AC `.sens` uses), keeping only the
         * ones sens=all|<name> selects. Counted/named ONCE here so the
         * per-frequency loop below can re-walk with the SAME filter and
         * index results 0..nsens-1 by POSITION (deterministic: the
         * circuit's device/model/instance/param structure does not change
         * between these two walks) instead of a name lookup every time. */
        if (sens_requested) {
            sgen *sg;
            int senscap = 0;
            for (sg = sgen_init(ckt, 0); sg; sgen_next(&sg)) {
                char namebuf[256];
                stb_sens_name(sg, namebuf, sizeof namebuf);
                if (!stb_sens_selected(job, namebuf))
                    continue;
                if (nsens >= senscap) {
                    senscap = senscap ? senscap * 2 : 8;
                    sensNames  = TREALLOC(char *, sensNames, senscap);
                    sensVecUids = TREALLOC(IFuid, sensVecUids, senscap);
                    dTdp_re = TREALLOC(double *, dTdp_re, senscap);
                    dTdp_im = TREALLOC(double *, dTdp_im, senscap);
                    sensRec = TREALLOC(STBsensRec, sensRec, senscap);
                }
                sensRec[nsens].dev = sg->dev;
                sensRec[nsens].param = sg->param;
                sensRec[nsens].is_instparam = sg->is_instparam;
                sensRec[nsens].model = sg->model;
                sensRec[nsens].instance = sg->instance;
                sensRec[nsens].value = sg->value;
                sensNames[nsens] = copy(namebuf);
                dTdp_re[nsens] = NULL;
                dTdp_im[nsens] = NULL;
                {
                    char vecname[300];
                    int ci;
                    snprintf(vecname, sizeof vecname, "dT_d%s", namebuf);
                    for (ci = 0; vecname[ci]; ci++)
                        if (vecname[ci] == ':') vecname[ci] = '_';
                    SPfrontEnd->IFnewUid(ckt, &sensVecUids[nsens], NULL, vecname, UID_OTHER, NULL);
                }
                nsens++;
            }
            if (nsens == 0) {
                SPfrontEnd->IFerrorf(ERR_WARNING,
                    "stb sens=: no matching IF_REAL instance/model parameter found in the "
                    "circuit -- T itself is still computed and printed normally, but no "
                    "dT/dp vectors or margin sensitivities are produced");
            } else {
                IFuid sensFreqUid;
                SPfrontEnd->IFnewUid(ckt, &sensFreqUid, NULL, "frequency", UID_OTHER, NULL);
                error = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob,
                    (IFuid) "STB dT/dp", sensFreqUid, IF_REAL,
                    nsens, sensVecUids, IF_COMPLEX, &sensPlot);
                if (error) return error;
                if (job->STBstepType != LINEAR)
                    SPfrontEnd->OUTattributes(sensPlot, NULL, OUT_SCALE_LOG, NULL);

                /* persistent per-frequency scratch, allocated once (size+1
                 * doubles: index 0 is the MNA ground "trash" slot, same
                 * convention as CKTrhs/CKTirhs -- see this file's existing
                 * node-0 guards). sensDeltaY mirrors cktsens.c's delta_Y:
                 * a throwaway matrix used only to isolate one device's own
                 * dY/dp stamp (SMPcClear'd fresh for every parameter). */
                {
                    int vsize = SMPmatSize(ckt->CKTmatrix) + 1;
                    sensDeltaY = TMALLOC(SMPmatrix, 1);
                    error = SMPnewMatrix(sensDeltaY, vsize - 1);
                    if (error) return error;
                    sens_lamv_re = TMALLOC(double, vsize); sens_lamv_im = TMALLOC(double, vsize);
                    sens_lami_re = TMALLOC(double, vsize); sens_lami_im = TMALLOC(double, vsize);
                    sens_x1_re = TMALLOC(double, vsize); sens_x1_im = TMALLOC(double, vsize);
                    sens_x2_re = TMALLOC(double, vsize); sens_x2_im = TMALLOC(double, vsize);
                    sens_dI = TMALLOC(double, vsize); sens_dII = TMALLOC(double, vsize);
                    sens_mv_re = TMALLOC(double, vsize); sens_mv_im = TMALLOC(double, vsize);
                    sensD = TMALLOC(STBcplx, 4 * nsens);
                    sensNoAc = TMALLOC(int, nsens);
                    sensJoint = TMALLOC(int, nsens);
                    { int kk; for (kk = 0; kk < nsens; kk++) sensJoint[kk] = 0; }
                }

                /* card #595: dx_op/dp for every selected parameter */
                if (getenv("STB_SENS_FIXED_OP") && atoi(getenv("STB_SENS_FIXED_OP")) != 0) {
                    fprintf(stdout, "stb sens=: STB_SENS_FIXED_OP set -- the operating point is held "
                            "fixed (bias-point term omitted, card #475 behaviour)\n");
                } else if (job->STBinstant) {
                    SPfrontEnd->IFerrorf(ERR_WARNING,
                        "stb sens= with instant: no DC operating point, so the bias-point term of "
                        "dT/dp is omitted (the operating point is held fixed)");
                } else {
                    int vsize = SMPmatSize(ckt->CKTmatrix) + 1;
                    double *save_rhs = ckt->CKTrhs, *save_rhsOld = ckt->CKTrhsOld;
                    long save_mode = ckt->CKTmode;
                    int save_bypass = ckt->CKTbypass;
                    double *np = TMALLOC(double, vsize), *nm = TMALLOC(double, vsize);
                    int k, j;
                    const char *e = getenv("STB_SENS_BIAS_DV");
                    char shared_names[512] = "";
                    int nshared = 0;
                    double dvt = (e && atof(e) > 0.0) ? atof(e) : 1e-5;

                    sens_xop = TMALLOC(double, vsize);
                    sens_bufOld = TMALLOC(double, vsize);
                    memcpy(sens_xop, ckt->CKTrhsOld, (size_t) vsize * sizeof(double));
                    sens_xop[0] = 0.0;
                    sens_state0 = TMALLOC(double, ckt->CKTnumStates + 1);
                    memcpy(sens_state0, ckt->CKTstate0, (size_t) ckt->CKTnumStates * sizeof(double));
                    sensDxdp = TMALLOC(double *, nsens);
                    sensH = TMALLOC(double, nsens);
                    ckt->CKTbypass = 0;

                    for (k = 0; k < nsens && !error; k++) {
                        STBsensRec *r = &sensRec[k];
                        double delta_var = (r->value != 0.0) ? r->value * STBsens_Delta : STBsens_AbsDelta;
                        double vmax = 0.0, amax = 0.0;
                        CKTnode *nd;

                        sensDxdp[k] = NULL;
                        sensH[k] = 0.0;
                        sensJoint[k] = 0;
                        if (!r->is_instparam) {
                            int ninst = 0;
                            GENinstance *gi;
                            for (gi = r->model->GENinstances; gi; gi = gi->GENnextInstance)
                                ninst++;
                            if (ninst > 1) {
                                /* collected into one warning after the loop */
                                size_t l = strlen(shared_names);
                                if (nshared < 8 && l + strlen(sensNames[k]) + 3 < sizeof shared_names)
                                    snprintf(shared_names + l, sizeof shared_names - l, "%s%s",
                                             l ? ", " : "", sensNames[k]);
                                nshared++;
                                continue;
                            }
                        }
                        if (stb_rec_setp(r, r->value + delta_var, ckt)) continue;
                        error = stb_sens_newton_step(ckt, sens_xop, sens_bufOld, np, vsize - 1);
                        if (!error) {
                            stb_rec_setp(r, r->value - delta_var, ckt);
                            error = stb_sens_newton_step(ckt, sens_xop, sens_bufOld, nm, vsize - 1);
                        }
                        stb_rec_setp(r, r->value, ckt);
                        if (error) break;
                        for (j = 1; j < vsize; j++) {
                            np[j] = (np[j] - nm[j]) / (2.0 * delta_var);
                            if (fabs(np[j]) > amax) amax = fabs(np[j]);
                        }
                        if (amax == 0.0) {
                            /* does not move the bias */
                            if (!stb_sens_linear_dev(r->dev)) {
                                sensJoint[k] = 1;
                                sensH[k] = (r->value != 0.0) ? fabs(r->value) * 1e-5 : STBsens_AbsDelta;
                                sens_nbias++;
                            }
                            continue;
                        }
                        for (nd = ckt->CKTnodes; nd; nd = nd->next)
                            if (nd->type == SP_VOLTAGE && nd->number > 0 && nd->number < vsize
                                    && fabs(np[nd->number]) > vmax)
                                vmax = fabs(np[nd->number]);
                        if (vmax == 0.0) vmax = amax;
                        sensDxdp[k] = TMALLOC(double, vsize);
                        memcpy(sensDxdp[k], np, (size_t) vsize * sizeof(double));
                        sensDxdp[k][0] = 0.0;
                        sensH[k] = dvt / vmax;
                        if (r->value != 0.0 && sensH[k] > fabs(r->value) * 1e-4)
                            sensH[k] = fabs(r->value) * 1e-4;
                        sensJoint[k] = !stb_sens_linear_dev(r->dev);
                        sens_nbias++;
                        if (getenv("STB_SENS_DEBUG")) {
                            for (nd = ckt->CKTnodes; nd; nd = nd->next)
                                if (nd->number > 0 && nd->number < vsize)
                                    fprintf(stderr, "STB_SENS_DEBUG %s dx/dp %s = %.9g (x_op %.9g)\n",
                                            sensNames[k], nd->name, np[nd->number], sens_xop[nd->number]);
                        }
                    }

                    if (nshared > 0)
                        SPfrontEnd->IFerrorf(ERR_WARNING,
                            "stb sens=: %d model parameter(s) of a model with several instances keep the "
                            "fixed operating point and the AC-stamp-only derivative (INCOMPLETE: the bias "
                            "shift and the device's DC-load quantities such as gm are not included; give the "
                            "device a model of its own for the full derivative): %s%s",
                            nshared, shared_names, nshared > 8 ? ", ..." : "");

                    /* back to the original linearisation, bit for bit */
                    memcpy(ckt->CKTstate0, sens_state0, (size_t) ckt->CKTnumStates * sizeof(double));
                    ckt->CKTrhsOld = save_rhsOld;
                    ckt->CKTrhs = np;
                    ckt->CKTmode = (save_mode & MODEUIC) | MODEDCOP | MODEINITSMSIG;
                    if (!error)
                        error = CKTload(ckt);
                    ckt->CKTrhs = save_rhs;
                    ckt->CKTmode = save_mode;
                    ckt->CKTbypass = save_bypass;
                    tfree(np);
                    tfree(nm);
                    if (error) return error;
                    if (sens_nbias > 0) {
                        sens_yp1_re = TMALLOC(double, vsize); sens_yp1_im = TMALLOC(double, vsize);
                        sens_yp2_re = TMALLOC(double, vsize); sens_yp2_im = TMALLOC(double, vsize);
                    }
                }
            }
        }

        freq = job->STBstartFreq;
    } else {
        freq = job->STBsaveFreq;
        job->STBsaveFreq = 0;
        error = SPfrontEnd->OUTpBeginPlot(NULL, NULL, NULL, NULL, 0, 666, NULL, 666, &stbPlot);
        if (sens_requested) {
            SPfrontEnd->IFerrorf(ERR_WARNING,
                "stb sens=: pause/resume is not supported for a sensitivity run in this "
                "fork -- re-run the whole sweep in one call");
            return E_PARMVAL;
        }
    }

    switch (job->STBstepType) {
    case DECADE:
    case OCTAVE:
        freqTol = job->STBfreqDelta * job->STBstopFreq * ckt->CKTreltol;
        break;
    case LINEAR:
        freqTol = job->STBfreqDelta * ckt->CKTreltol;
        break;
    default:
        return E_BADPARM;
    }

    /* card #459/#460: NOW read gm/gmbs and the node indices off the MOS1
     * or BSIM4 instance -- after CKTop()+CKTload() above (fresh start) or
     * a prior call's op point (resume/pause: *gm lives on the persistent
     * instance struct, not this function's stack, so it is still valid
     * here without recomputing). Doing this any earlier reads a
     * never-computed (zero) gm on a fresh run -- see this file's comment
     * where device_mode is set, above. */
    if (device_mode) {
        bool is_bsim4 = (inst->GENmodPtr->GENmodType == CKTtypelook("BSIM4"));
        int perr = is_bsim4 ? stb_bsim4_probe(inst, &devp) : stb_mos1_probe(inst, &devp);
        if (perr) {
            if (is_bsim4)
                SPfrontEnd->IFerrorf(ERR_WARNING,
                    "stb probe '%s': BSIM4 device mode is quasi-static only in "
                    "this fork (card #460) -- this instance (or its model) has "
                    "acnqsmod=1 (NQS on); refusing rather than nulling the "
                    "wrong (real-only) scalar", job->STBprobeName);
            else
                SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe '%s': could not read device fields", job->STBprobeName);
            ckt->CKTcurrentAnalysis = 0;
            return E_NOTFOUND;
        }
    }

    ckt->CKTcurrentAnalysis = DOING_AC;
    size = SMPmatSize(ckt->CKTmatrix);

    while (freq <= job->STBstopFreq + freqTol) {
        if (SPfrontEnd->IFpauseTest()) {
            job->STBsaveFreq = freq;
            ckt->CKTcurrentAnalysis = 0;
            return E_PAUSE;
        }

        ckt->CKTomega = 2.0 * M_PI * freq;
        ckt->CKTmode = (ckt->CKTmode & MODEUIC) | MODEAC;

        /* load the complex matrix once; drop every AC source's own
         * stimulus (this analysis supplies the only excitation), same
         * pattern as loopan.c's per-run RHS clearing. */
        error = CKTacLoad(ckt);
        if (error) { ckt->CKTcurrentAnalysis = 0; return error; }
        for (i = 0; i <= size; i++) {
            ckt->CKTrhs[i] = 0.0;
            ckt->CKTirhs[i] = 0.0;
        }

        /* factor once (retry with reordering on a singular pivot, as
         * NIacIter/NIspPreload do) */
    refactor:
        if (ckt->CKTniState & NIACSHOULDREORDER) {
            error = SMPcReorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTpivotRelTol, &ignore);
            ckt->CKTniState &= ~NIACSHOULDREORDER;
            if (error) { ckt->CKTcurrentAnalysis = 0; return error; }
        } else {
            error = SMPcLUfac(ckt->CKTmatrix, ckt->CKTpivotAbsTol);
            if (error == E_SINGULAR) {
                ckt->CKTniState |= NIACSHOULDREORDER;
                goto refactor;
            }
            if (error) { ckt->CKTcurrentAnalysis = 0; return error; }
        }

        if (device_mode) {
            /* card #459: device-based gain-nulling, ONE extra back-solve.
             * RHS = u (the gm/gmbs output current's incidence vector):
             * +1 at the internal drain, -1 at the internal source (see
             * stb_mos1_probe()'s header comment for the normal/reverse-mode
             * node choice). No probe branch equation exists in this mode.
             *
             * Ground (node index 0) is NOT a real unknown: it is excluded
             * from the factored system, and CKTrhs[0]/CKTirhs[0] is a
             * "trash" slot that SMPcSolve does not touch -- so after the
             * solve it still holds whatever this loop last wrote there,
             * NOT the (always-exactly-0) ground voltage. Found empirically
             * while smoke-testing this card (a grounded MOSFET source, the
             * ordinary case, gave a frequency-INDEPENDENT s exactly equal
             * to gm -- i.e. w[node_s] was silently reading back the -1.0
             * this code itself had just written to the RHS, not a solved
             * value). Guard both the write (skip it -- not needed for
             * correctness, ground contributes nothing to a real solve
             * either way) and every read (node 0 reads back as exactly
             * 0+0j, never the trash slot). This is a general MNA fact, not
             * specific to device mode; the existing loop-mode code above
             * does not hit it because none of this fork's reference decks
             * probe a wire with a literal-ground terminal -- see the card
             * #459 report's "what was not done" for this latent issue. */
            for (i = 0; i <= size; i++) {
                ckt->CKTrhs[i] = 0.0;
                ckt->CKTirhs[i] = 0.0;
            }
            if (devp.node_d != 0) ckt->CKTrhs[devp.node_d] += 1.0;
            if (devp.node_s != 0) ckt->CKTrhs[devp.node_s] -= 1.0;
            SMPcSolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTirhs, ckt->CKTrhsSpare, ckt->CKTirhsSpare);
            /* s = gm * v_gm^T w, w = Y^-1 u, v_gm = +1 at gate, -1 at the
             * internal source (gm-only production default; gmbs NOT
             * folded in here -- see stb_mos1_probe()'s doc comment on the
             * STBdevProbe struct, above). */
            {
                double wg_re = devp.node_g ? ckt->CKTrhs[devp.node_g]  : 0.0;
                double wg_im = devp.node_g ? ckt->CKTirhs[devp.node_g] : 0.0;
                double ws_re = devp.node_s ? ckt->CKTrhs[devp.node_s]  : 0.0;
                double ws_im = devp.node_s ? ckt->CKTirhs[devp.node_s] : 0.0;
                s.re = devp.gm * (wg_re - ws_re);
                s.im = devp.gm * (wg_im - ws_im);
            }

            if (getenv("STB_DEBUG"))
                fprintf(stderr, "STB_DEBUG f=%g device node_d=%d node_s=%d node_g=%d node_b=%d "
                    "gm=%.10g gmbs=%.10g s=%.10g%+.10gj\n",
                    freq, devp.node_d, devp.node_s, devp.node_g, devp.node_b,
                    devp.gm, devp.gmbs, s.re, s.im);

            /* T = 1/(1-s) - 1 (research/cadence-stb-pstb.md sec 2.5's
             * rank-1 form, derived from the matrix determinant lemma --
             * see the card #459 report). loopGain = -T, same convention
             * as loop mode. */
            {
                STBcplx one = { 1.0, 0.0 };
                t = stb_csub(stb_cdiv(one, stb_csub(one, s)), one);
                lg.re = -t.re; lg.im = -t.im;
            }
        } else {

        /* RHS #1: unit voltage injection in the probe branch equation */
        ckt->CKTrhs[branch_probe] = -1.0;
        ckt->CKTirhs[branch_probe] = 0.0;
        /* card #457 probe=Vp,Vn mode=dm|cm: drive branch n's voltage
         * injection in lockstep, sign-flipped for dm (see dmcm_sign above).
         * branch n's own v/i are never read -- only branch p's (below). */
        if (job->STBprobeName2Given) {
            ckt->CKTrhs[branch_probe2] = -1.0 * dmcm_sign;
            ckt->CKTirhs[branch_probe2] = 0.0;
        }
        SMPcSolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTirhs, ckt->CKTrhsSpare, ckt->CKTirhsSpare);
        /* card #460: node index 0 (true ground) is not a real unknown --
         * it is excluded from the factored system, and CKTrhs[0]/
         * CKTirhs[0] is a "trash" slot SMPcSolve never touches, so it
         * still holds whatever this loop last wrote there. Ground's
         * voltage is always exactly 0 by definition; read that directly
         * rather than the trash slot (same hazard and same guard pattern
         * as the device-mode fix in c931b1f4c, and as devp.node_d/
         * devp.node_s above). Only node_b can be 0 here -- node_a is
         * never dereferenced into CKTrhs/CKTirhs in this file, so it
         * needs no guard; see the warning above where node_a/node_b were
         * resolved for why a node_a==0 probe is ALSO degenerate anyway
         * (T=0), just not via this trash-slot mechanism. */
        v1.re = node_b ? ckt->CKTrhs[node_b]  : 0.0;
        v1.im = node_b ? ckt->CKTirhs[node_b] : 0.0;
        /* The raw MNA branch unknown for a V-source is the current flowing
         * INTO its positive terminal (node_a) from the external circuit;
         * ngspice's printed "i(Vxxx)" (what tools/stbcore.tian is written
         * against, and what Tian's own convention assumes) is the current
         * the source delivers OUT of its positive terminal -- the negative
         * of that raw unknown. Both runs need this sign flip; confirmed
         * empirically against tools/stb.py on the R1a/R1b/R1d reference
         * decks with this same build (see tests/test_stb_native.py and the
         * card #456 report). */
        i1.re = -ckt->CKTrhs[branch_probe]; i1.im = -ckt->CKTirhs[branch_probe];

        /* card #475: save the whole x1 solution vector (not just v1/i1)
         * BEFORE it is overwritten by run 2's solve -- the adjoint-
         * sensitivity contraction below needs the full vector for the
         * (dY/dp)*x1 matrix-vector product. */
        if (sens_requested && nsens > 0)
            for (i = 0; i <= size; i++) {
                sens_x1_re[i] = ckt->CKTrhs[i];
                sens_x1_im[i] = ckt->CKTirhs[i];
            }

        /* RHS #2: unit current injection into the probe's other terminal
         * (node_b); probe branch equation back to 0 (unperturbed) */
        for (i = 0; i <= size; i++) {
            ckt->CKTrhs[i] = 0.0;
            ckt->CKTirhs[i] = 0.0;
        }
        /* card #460: guard the write too -- ground has no KCL row of its
         * own to receive this injection (it is the reference, not an
         * unknown), so writing here is a no-op for the real solve either
         * way; skip it rather than leave a stray value sitting in the
         * trash slot for the read below (and for the STB_DEBUG dump) to
         * pick up. */
        if (node_b) ckt->CKTrhs[node_b] += 1.0;
        /* localgnd=: reference the injected current to the given local
         * ground node instead of true ground (same two-terminal-source
         * pattern as the unmerged loop_analysis branch's LOOPan() used for
         * its optional refnode=). NOT independently validated numerically
         * (no reference deck here uses localgnd=) -- see the card #456
         * report's "what was not done". Guarded the same way (card #460):
         * localgnd=0 would just mean "local ground IS true ground", i.e.
         * no local-reference correction at all -- a no-op, not a special
         * case, once guarded. */
        if (job->STBlocalGndGiven && job->STBlocalGnd && job->STBlocalGnd->number)
            ckt->CKTrhs[job->STBlocalGnd->number] -= 1.0;
        /* card #457 probe=Vp,Vn mode=dm|cm: drive branch n's current
         * injection in lockstep (same dmcm_sign as run 1, above). card
         * #460: guarded the same way as node_b just above. */
        if (job->STBprobeName2Given && node_b2)
            ckt->CKTrhs[node_b2] += dmcm_sign * 1.0;
        SMPcSolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTirhs, ckt->CKTrhsSpare, ckt->CKTirhsSpare);
        /* card #460: same ground-read guard as v1, above. */
        v2.re = node_b ? ckt->CKTrhs[node_b]  : 0.0;
        v2.im = node_b ? ckt->CKTirhs[node_b] : 0.0;
        i2.re = -ckt->CKTrhs[branch_probe]; i2.im = -ckt->CKTirhs[branch_probe];

        /* card #475: save x2 the same way x1 was saved above. */
        if (sens_requested && nsens > 0)
            for (i = 0; i <= size; i++) {
                sens_x2_re[i] = ckt->CKTrhs[i];
                sens_x2_im[i] = ckt->CKTirhs[i];
            }

        if (getenv("STB_DEBUG"))
            fprintf(stderr, "STB_DEBUG f=%g node_a=%d node_b=%d branch=%d "
                "v1=%.10g%+.10gj i1=%.10g%+.10gj v2=%.10g%+.10gj i2=%.10g%+.10gj\n",
                freq, node_a, node_b, branch_probe,
                v1.re, v1.im, i1.re, i1.im, v2.re, v2.im, i2.re, i2.im);

        /* Tian eq. 30 (tools/stbcore.tian): s = 2*(i1*v2 - v1*i2) + v1 + i2
         *                                    T = s/(1-s), loopGain = -T */
        {
            STBcplx one = { 1.0, 0.0 };
            STBcplx diff = stb_csub(stb_cmul(i1, v2), stb_cmul(v1, i2));
            STBcplx two_diff = { 2.0 * diff.re, 2.0 * diff.im };
            s = stb_cadd(stb_cadd(two_diff, v1), i2);
            t = stb_cdiv(s, stb_csub(one, s));
            lg.re = -t.re; lg.im = -t.im;
        }
        }

        if (npts >= cap) {
            cap = cap ? cap * 2 : 64;
            freqs = TREALLOC(double, freqs, cap);
            Tarr = TREALLOC(STBcplx, Tarr, cap);
            LGarr = TREALLOC(STBcplx, LGarr, cap);
            if (sens_requested && nsens > 0) {
                int kk;
                for (kk = 0; kk < nsens; kk++) {
                    dTdp_re[kk] = TREALLOC(double, dTdp_re[kk], cap);
                    dTdp_im[kk] = TREALLOC(double, dTdp_im[kk], cap);
                }
            }
        }
        freqs[npts] = freq;
        Tarr[npts] = t;
        LGarr[npts] = lg;

        /* card #475: adjoint sensitivity, this frequency. Two extra
         * (TRANSPOSED, same factorisation) solves for lambda_v/lambda_i,
         * THEN a cheap per-parameter contraction (no further solves) --
         * see this file's design comment above stb_sens_dot(). */
        if (sens_requested && nsens > 0) {
            double *save_rhs = ckt->CKTrhs, *save_irhs = ckt->CKTirhs;
            SMPmatrix *save_matrix = ckt->CKTmatrix;
            sgen *sg;
            int k;
            STBcplx one;
            one.re = 1.0; one.im = 0.0;

            for (i = 0; i <= size; i++) { sens_lamv_re[i] = 0.0; sens_lamv_im[i] = 0.0; }
            if (node_b) sens_lamv_re[node_b] = 1.0;
            SMPcaSolve(ckt->CKTmatrix, sens_lamv_re, sens_lamv_im, ckt->CKTrhsSpare, ckt->CKTirhsSpare);

            for (i = 0; i <= size; i++) { sens_lami_re[i] = 0.0; sens_lami_im[i] = 0.0; }
            sens_lami_re[branch_probe] = -1.0;
            SMPcaSolve(ckt->CKTmatrix, sens_lami_re, sens_lami_im, ckt->CKTrhsSpare, ckt->CKTirhsSpare);

            ckt->CKTrhs = sens_dI; ckt->CKTirhs = sens_dII; ckt->CKTmatrix = sensDeltaY;

            k = 0;
            for (sg = sgen_init(ckt, 0); sg; sgen_next(&sg)) {
                char namebuf[256];
                IFvalue nvalue, ovalue;
                double delta_var;
                STBcplx dv1, di1, dv2, di2, dT;
                int (*setupfn)(SMPmatrix *, GENmodel *, CKTcircuit *, int *);

                stb_sens_name(sg, namebuf, sizeof namebuf);
                if (!stb_sens_selected(job, namebuf))
                    continue;
                /* card #475 bugfix (found while smoke-testing this card):
                 * sgen_next()'s isolation trick temporarily unlinks each
                 * instance from its model's chain and only relinks it as
                 * the generator advances PAST it -- an early `break` out
                 * of this loop (as an earlier version of this code did,
                 * once k>=nsens) abandons the walk before sgen_next() ever
                 * reaches its own "done" state, so the LAST instance(s)
                 * touched are left permanently unlinked from CKThead's
                 * model->GENinstances chains: the very next CKTacLoad()
                 * (next frequency point) then silently skips them, giving
                 * WRONG (missing-stamp) T *and* garbage/zero sensitivities
                 * from that point on. Never break here -- `continue`, the
                 * same pattern cktsens.c's own per-frequency sgen walk
                 * uses for its filtered-out parameters, always lets
                 * sgen_next() run the walk to completion so it can restore
                 * every temporarily-isolated instance. k can safely exceed
                 * nsens-1 only if this walk finds MORE matches than the
                 * dry run did (should not happen -- circuit structure is
                 * unchanged in between -- but the store below is guarded
                 * defensively anyway rather than indexing out of bounds). */
                if (k >= nsens) {
                    k++;
                    continue;
                }

                SMPcClear(sensDeltaY);
                for (i = 0; i <= size; i++) { sens_dI[i] = 0.0; sens_dII[i] = 0.0; }

                ckt->CKTnumStates = sg->istate;
                setupfn = DEVices[sg->dev]->DEVsetup;
                if (setupfn)
                    setupfn(sensDeltaY, sg->model, ckt, &ckt->CKTnumStates);

                sensNoAc[k] = 0;
                if (stb_sens_load_ac(sg, ckt)) {
                    sensNoAc[k] = 1;
                    k++; continue;  /* no AC load fn on this device -- dT/dp undefined, report 0 rather than garbage */
                }

                delta_var = (sg->value != 0.0) ? sg->value * STBsens_Delta : STBsens_AbsDelta;
                nvalue.rValue = sg->value + delta_var;
                stb_sens_setp(sg, ckt, &nvalue);

                SMPconstMult(sensDeltaY, -1.0);
                for (i = 0; i <= size; i++) { sens_dI[i] = -sens_dI[i]; sens_dII[i] = -sens_dII[i]; }

                stb_sens_temp(sg, ckt);
                stb_sens_load_ac(sg, ckt);

                ovalue.rValue = sg->value;
                stb_sens_setp(sg, ckt, &ovalue);
                stb_sens_temp(sg, ckt);

                /* card #475 bugfix (found while smoke-testing this card,
                 * the SECOND real bug this pass found -- see git log for
                 * the first, the early-break above): DEVsetup(sensDeltaY,
                 * ...) just above does not merely allocate STRUCTURAL
                 * positions in sensDeltaY -- for devices that cache
                 * matrix-element pointers on the INSTANCE itself (e.g. a
                 * capacitor's own posPosPtr-style fields, written once by
                 * DEVsetup and read directly by DEVload/DEVacLoad ever
                 * after, with no further SMPmakeElt call), it REPOINTS
                 * those cached pointers at sensDeltaY's storage -- a
                 * PERSISTENT, PER-INSTANCE side effect that survives this
                 * function restoring ckt->CKTmatrix back to the real Y.
                 * The very next frequency's ordinary CKTacLoad() then
                 * silently loads this device's real contribution into the
                 * now-defunct sensDeltaY instead of the real matrix --
                 * found empirically: T itself (not just dT/dp) went wrong
                 * from the SECOND swept frequency once 2+ sens= devices
                 * with this property were requested (a lone resistor
                 * alone did not show it -- RESsetup apparently does not
                 * cache pointers the same way; a lone capacitor did).
                 * cktsens.c's own DC/AC `.sens` never hits this because it
                 * calls CKTunsetup()+CKTsetup() fresh before EVERY
                 * frequency point, which re-runs every device's REAL
                 * DEVsetup(ckt->CKTmatrix, ...) and so re-fixes the
                 * pointers before the next real load -- STBan()'s
                 * existing per-frequency loop does no such reset (by
                 * design, it is much cheaper: one factor + two solves).
                 * Fix: immediately re-run this SAME device's DEVsetup
                 * against the REAL matrix (save_matrix) to restore its
                 * pointers before moving to the next parameter -- a
                 * throwaway numStates counter is used (seeded with this
                 * instance's own already-correct sg->istate) so this
                 * repair call cannot perturb ckt->CKTnumStates either. */
                if (setupfn) {
                    int throwaway_states = sg->istate;
                    setupfn(save_matrix, sg->model, ckt, &throwaway_states);
                }

                /* sensDeltaY/sens_dI/sens_dII now hold delta*(dY/dp),
                 * delta*(db/dp) for this ISOLATED device (sgen's own
                 * temporary instance-unlinking keeps DEVsetup/DEVacLoad
                 * from touching any sibling instance -- cktsgen.c). */

                /* card #595: spMultiply writes only the rows present in
                 * sensDeltaY's structure; clear the rest first */
                for (i = 0; i <= size; i++) { sens_mv_re[i] = 0.0; sens_mv_im[i] = 0.0; }
                SMPmultiply(sensDeltaY, sens_mv_re, sens_x1_re, sens_mv_im, sens_x1_im);
                /* card #595: db/dp is zero -- .stb drops every AC source's own
                 * stimulus (the probe injection is the only excitation), so the
                 * device's AC-stimulus stamp in sens_dI/sens_dII is not used
                 * (it gave Vin_acmag and the like a nonzero dT/dp before) */
                for (i = 0; i <= size; i++) {
                    sens_mv_re[i] = -sens_mv_re[i];
                    sens_mv_im[i] = -sens_mv_im[i];
                }
                {
                    STBcplx dvv = stb_sens_dot(sens_lamv_re, sens_lamv_im, sens_mv_re, sens_mv_im, size);
                    STBcplx dii = stb_sens_dot(sens_lami_re, sens_lami_im, sens_mv_re, sens_mv_im, size);
                    dv1.re = dvv.re / delta_var; dv1.im = dvv.im / delta_var;
                    di1.re = dii.re / delta_var; di1.im = dii.im / delta_var;
                }

                for (i = 0; i <= size; i++) { sens_mv_re[i] = 0.0; sens_mv_im[i] = 0.0; }
                SMPmultiply(sensDeltaY, sens_mv_re, sens_x2_re, sens_mv_im, sens_x2_im);
                /* card #595: db/dp is zero -- .stb drops every AC source's own
                 * stimulus (the probe injection is the only excitation), so the
                 * device's AC-stimulus stamp in sens_dI/sens_dII is not used
                 * (it gave Vin_acmag and the like a nonzero dT/dp before) */
                for (i = 0; i <= size; i++) {
                    sens_mv_re[i] = -sens_mv_re[i];
                    sens_mv_im[i] = -sens_mv_im[i];
                }
                {
                    STBcplx dvv = stb_sens_dot(sens_lamv_re, sens_lamv_im, sens_mv_re, sens_mv_im, size);
                    STBcplx dii = stb_sens_dot(sens_lami_re, sens_lami_im, sens_mv_re, sens_mv_im, size);
                    dv2.re = dvv.re / delta_var; dv2.im = dvv.im / delta_var;
                    di2.re = dii.re / delta_var; di2.im = dii.im / delta_var;
                }

                sensD[4 * k + 0] = dv1;
                sensD[4 * k + 1] = di1;
                sensD[4 * k + 2] = dv2;
                sensD[4 * k + 3] = di2;
                k++;
            }

            ckt->CKTrhs = save_rhs; ckt->CKTirhs = save_irhs; ckt->CKTmatrix = save_matrix;

            /* card #595: the bias-point term, (dY/dx)(dx_op/dp) x, by a
             * central difference of the re-linearised matrix (see
             * stb_rec_setp()'s block comment) */
            if (sens_nbias > 0) {
                double *save_rhsOld = ckt->CKTrhsOld;
                int save_bypass = ckt->CKTbypass;
                ckt->CKTbypass = 0;
                for (k = 0; k < nsens; k++) {
                    int sgn;
                    double h = sensH[k];
                    if ((!sensDxdp[k] && !sensJoint[k]) || sensNoAc[k]) continue;
                    if (sensJoint[k]) {
                        /* the joint derivative replaces the stamp difference */
                        int q;
                        for (q = 0; q < 4; q++) { sensD[4 * k + q].re = 0.0; sensD[4 * k + q].im = 0.0; }
                    }
                    for (sgn = 0; sgn < 2; sgn++) {
                        if (sensJoint[k])
                            stb_rec_setp(&sensRec[k], sensRec[k].value + (sgn ? -h : h), ckt);
                        error = stb_sens_shift_load(ckt, sens_xop, sensDxdp[k], sgn ? -h : h,
                                                    sens_bufOld, sens_dI, sens_dII, size);
                        if (sensJoint[k])
                            stb_rec_setp(&sensRec[k], sensRec[k].value, ckt);
                        if (error) break;
                        if (!sgn)
                            for (i = 0; i <= size; i++)
                                sens_yp1_re[i] = sens_yp1_im[i] = sens_yp2_re[i] = sens_yp2_im[i] = 0.0;
                        if (!sgn) {
                            SMPmultiply(ckt->CKTmatrix, sens_yp1_re, sens_x1_re, sens_yp1_im, sens_x1_im);
                            SMPmultiply(ckt->CKTmatrix, sens_yp2_re, sens_x2_re, sens_yp2_im, sens_x2_im);
                        } else {
                            int pass;
                            for (pass = 0; pass < 2; pass++) {
                                double *yr = pass ? sens_yp2_re : sens_yp1_re;
                                double *yi = pass ? sens_yp2_im : sens_yp1_im;
                                STBcplx dvv, dii;
                                for (i = 0; i <= size; i++) { sens_mv_re[i] = 0.0; sens_mv_im[i] = 0.0; }
                                if (pass)
                                    SMPmultiply(ckt->CKTmatrix, sens_mv_re, sens_x2_re, sens_mv_im, sens_x2_im);
                                else
                                    SMPmultiply(ckt->CKTmatrix, sens_mv_re, sens_x1_re, sens_mv_im, sens_x1_im);
                                /* -(Y+ x - Y- x) / 2h */
                                for (i = 0; i <= size; i++) {
                                    sens_mv_re[i] = (sens_mv_re[i] - yr[i]) / (2.0 * h);
                                    sens_mv_im[i] = (sens_mv_im[i] - yi[i]) / (2.0 * h);
                                }
                                dvv = stb_sens_dot(sens_lamv_re, sens_lamv_im, sens_mv_re, sens_mv_im, size);
                                dii = stb_sens_dot(sens_lami_re, sens_lami_im, sens_mv_re, sens_mv_im, size);
                                sensD[4 * k + 2 * pass + 0] = stb_cadd(sensD[4 * k + 2 * pass + 0], dvv);
                                sensD[4 * k + 2 * pass + 1] = stb_cadd(sensD[4 * k + 2 * pass + 1], dii);
                            }
                        }
                    }
                    if (error) break;
                }
                /* back to the operating-point linearisation */
                memcpy(ckt->CKTstate0, sens_state0, (size_t) ckt->CKTnumStates * sizeof(double));
                ckt->CKTrhsOld = save_rhsOld;
                ckt->CKTrhs = sens_dI; ckt->CKTirhs = sens_dII;
                {
                    long save_mode = ckt->CKTmode;
                    ckt->CKTmode = (save_mode & MODEUIC) | MODEDCOP | MODEINITSMSIG;
                    if (!error)
                        error = CKTload(ckt);
                    ckt->CKTmode = save_mode;
                }
                ckt->CKTrhs = save_rhs; ckt->CKTirhs = save_irhs;
                ckt->CKTbypass = save_bypass;
                if (error) { ckt->CKTcurrentAnalysis = 0; return error; }
            }

            for (k = 0; k < nsens; k++) {
                STBcplx dv1, di1, dv2, di2, dT;
                if (sensNoAc[k]) {
                    dTdp_re[k][npts] = 0.0; dTdp_im[k][npts] = 0.0;
                    continue;
                }
                dv1 = sensD[4 * k + 0]; di1 = sensD[4 * k + 1];
                dv2 = sensD[4 * k + 2]; di2 = sensD[4 * k + 3];
                /* s = 2*(i1*v2 - v1*i2) + v1 + i2
                 * ds/dp = 2*(di1*v2 + i1*dv2 - dv1*i2 - v1*di2) + dv1 + di2
                 * T = s/(1-s)  =>  dT/dp = ds/dp / (1-s)^2 */
                {
                    STBcplx term1 = stb_cmul(di1, v2);
                    STBcplx term2 = stb_cmul(i1, dv2);
                    STBcplx term3 = stb_cmul(dv1, i2);
                    STBcplx term4 = stb_cmul(v1, di2);
                    STBcplx inner = stb_csub(stb_cadd(term1, term2), stb_cadd(term3, term4));
                    STBcplx two_inner; two_inner.re = 2.0 * inner.re; two_inner.im = 2.0 * inner.im;
                    STBcplx dsdp = stb_cadd(stb_cadd(two_inner, dv1), di2);
                    STBcplx one_minus_s = stb_csub(one, s);
                    STBcplx denom = stb_cmul(one_minus_s, one_minus_s);
                    dT = stb_cdiv(dsdp, denom);
                }
                dTdp_re[k][npts] = dT.re;
                dTdp_im[k][npts] = dT.im;
            }
        }

        npts++;

        {
            IFvalue freqData, valueData;
            IFcomplex data[2];
            data[0].real = lg.re; data[0].imag = lg.im;
            data[1].real = t.re;  data[1].imag = t.im;
            freqData.rValue = freq;
            valueData.v.numValue = 2;
            valueData.v.vec.cVec = data;
            SPfrontEnd->OUTpData(stbPlot, &freqData, &valueData);
        }

        if (sens_requested && nsens > 0) {
            IFvalue freqData, valueData;
            IFcomplex *data = TMALLOC(IFcomplex, nsens);
            int kk;
            for (kk = 0; kk < nsens; kk++) {
                data[kk].real = dTdp_re[kk][npts - 1];
                data[kk].imag = dTdp_im[kk][npts - 1];
            }
            freqData.rValue = freq;
            valueData.v.numValue = nsens;
            valueData.v.vec.cVec = data;
            SPfrontEnd->OUTpData(sensPlot, &freqData, &valueData);
            tfree(data);
        }

        switch (job->STBstepType) {
        case DECADE:
        case OCTAVE:
            freq *= job->STBfreqDelta;
            if (job->STBfreqDelta == 1) goto endsweep;
            break;
        case LINEAR:
            freq += job->STBfreqDelta;
            if (job->STBfreqDelta == 0) goto endsweep;
            break;
        default:
            ckt->CKTcurrentAnalysis = 0;
            return E_INTERN;
        }
    }
endsweep:
    SPfrontEnd->OUTendPlot(stbPlot);
    stbPlot = NULL;
    if (sensPlot) {
        SPfrontEnd->OUTendPlot(sensPlot);
        sensPlot = NULL;
    }
    ckt->CKTcurrentAnalysis = 0;

    /* ---- margins: post-process the whole T(f) trace, mirroring
     * tools/stbcore.margins() exactly (log-f interpolation, unwrapped
     * phase, worst-of-several-crossings, "no crossing" is never a number) */
    {
        double prevp = 0.0;
        phs = TMALLOC(double, npts);
        for (i = 0; i < npts; i++) {
            double p = 180.0 / M_PI * stb_cph(Tarr[i]);
            if (i > 0) {
                while (p - prevp > 180.0) p -= 360.0;
                while (p - prevp < -180.0) p += 360.0;
            }
            phs[i] = p;
            prevp = p;
        }
        for (i = 0; i + 1 < npts; i++) {
            double f0 = freqs[i], f1 = freqs[i + 1];
            double m0 = 20.0 * log10(stb_cmag(Tarr[i])), m1 = 20.0 * log10(stb_cmag(Tarr[i + 1]));
            double p0 = phs[i], p1 = phs[i + 1];
            double lo, hi, target;
            int n;

            if (m0 * m1 < 0 || (m1 == 0.0 && m0 != 0.0)) {
                double fc = interp_logf(f0, f1, m0, m1, 0.0);
                double t2 = (0.0 - m0) / (m1 - m0);
                double ph = p0 + t2 * (p1 - p0);
                /* card #464: aux carries the crossing direction, mirroring
                 * tools/stbcore.margins' "falling" flag -- 1.0 when |T| is
                 * decreasing through this crossing (the usual stability-margin
                 * crossing), 0.0 when it is rising (PM here, even < 180 deg, is
                 * not the textbook margin -- see the printed warning below). */
                pmx = push_crossing(pmx, &npm, &capm, fc, 180.0 + ph, m1 < m0 ? 1.0 : 0.0, t2, i);
            }
            lo = p0 < p1 ? p0 : p1;
            hi = p0 < p1 ? p1 : p0;
            n = (int) ceil((lo + 180.0) / 360.0);
            while (-180.0 + 360.0 * n <= hi) {
                target = -180.0 + 360.0 * n;
                if (target != p0) {
                    double fc = interp_logf(f0, f1, p0, p1, target);
                    double t2 = (target - p0) / (p1 - p0);
                    double mag = m0 + t2 * (m1 - m0);
                    gmx = push_crossing(gmx, &ngm, &capg, fc, -mag, target, t2, i);
                }
                n++;
            }
        }
    }

    {
        int pmi = -1, gmi = -1;
        double pmbest = 0.0, gmbest = 0.0;
        /* card #464, narrowed by card #476: Re(T) < 0 alone at the lowest
         * swept frequency false-triggered on healthy integrator-type loops
         * (a few degrees of extra lag past -90 deg on a near-integrator
         * compensator branch makes Re(T) slightly negative with no real
         * instability -- found on tests/fixtures/compopt_r3_averaged_detuned.cir
         * at fstart=10 Hz). DECISION (maintainer, option A, card #476): let
         * arg0 = arg(T(fmin)) in degrees, mag0 = |T(fmin)|.
         *   - low_f_posfb (the WARNING): |arg0| > 150 deg AND mag0 > 1 --
         *     positive feedback WITH loop gain, cannot be stable in the
         *     usual sense. Mirrors tools/stbcore.margins' "low_f_positive_
         *     feedback" / "note" exactly (found on the pre-fix
         *     refs/r1c_fully_diff.cir, 2026-09-29: this printed "phase
         *     margin 274/186 deg" for a loop that was actually unstable at
         *     DC).
         *   - low_f_posfb_low_gain (a note, not the warning): |arg0| > 150
         *     deg but mag0 <= 1.
         *   - low_f_integrator (a neutral label): arg0 in [-120, -60] deg
         *     (-90 +/- 30), ordinary integrator-like low-frequency behaviour.
         * freqs[] is ascending (DECADE/OCTAVE/LINEAR sweeps only), so index 0
         * is always the lowest frequency swept. */
        double arg0 = (npts > 0) ? 180.0 / M_PI * stb_cph(Tarr[0]) : 0.0;
        double mag0 = (npts > 0) ? stb_cmag(Tarr[0]) : 0.0;
        int low_f_posfb = (npts > 0) && (fabs(arg0) > 150.0) && (mag0 > 1.0);
        int low_f_posfb_low_gain = (npts > 0) && (fabs(arg0) > 150.0) && (mag0 <= 1.0);
        int low_f_integrator = (npts > 0) && (arg0 >= -120.0) && (arg0 <= -60.0);
        for (i = 0; i < npm; i++)
            if (pmi < 0 || pmx[i].value < pmbest) { pmi = i; pmbest = pmx[i].value; }
        for (i = 0; i < ngm; i++)
            if (gmi < 0 || gmx[i].value < gmbest) { gmi = i; gmbest = gmx[i].value; }

        {
            IFuid uids[8];
            double outv[8];
            runDesc *mplot = NULL;
            double nan = strtod("NaN", NULL);

            SPfrontEnd->IFnewUid(ckt, &uids[0], NULL, "phaseMargin", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[1], NULL, "phaseMarginFreq", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[2], NULL, "gainMargin", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[3], NULL, "gainMarginFreq", UID_OTHER, NULL);
            /* card #464: lowFPositiveFeedback (1/0) and phaseMarginFalling
             * (1.0 = usual falling |T| crossing, 0.0 = rising, NaN = no PM
             * crossing at all) -- the scalar/flag mechanism the warning uses.
             * card #476 (additive): lowFPositiveFeedbackLowGain (1/0, the
             * note case: |arg0|>150 deg but mag0<=1) and lowFIntegrator
             * (1/0, the neutral -90+/-30 deg label). */
            SPfrontEnd->IFnewUid(ckt, &uids[4], NULL, "lowFPositiveFeedback", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[5], NULL, "phaseMarginFalling", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[6], NULL, "lowFPositiveFeedbackLowGain", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[7], NULL, "lowFIntegrator", UID_OTHER, NULL);

            outv[0] = pmi >= 0 ? pmx[pmi].value : nan;
            outv[1] = pmi >= 0 ? pmx[pmi].freq : nan;
            outv[2] = gmi >= 0 ? gmx[gmi].value : nan;
            outv[3] = gmi >= 0 ? gmx[gmi].freq : nan;
            outv[4] = low_f_posfb ? 1.0 : 0.0;
            outv[5] = pmi >= 0 ? pmx[pmi].aux : nan;
            outv[6] = low_f_posfb_low_gain ? 1.0 : 0.0;
            outv[7] = low_f_integrator ? 1.0 : 0.0;

            error = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob,
                (IFuid) "STB Margins", NULL, 0, 8, uids, IF_REAL, &mplot);
            if (!error) {
                IFvalue refval, outdata;
                refval.rValue = 0;
                outdata.v.numValue = 8;
                outdata.v.vec.rVec = outv;
                SPfrontEnd->OUTpData(mplot, &refval, &outdata);
                SPfrontEnd->OUTendPlot(mplot);
            }
        }

        printf("stb%s%s: ", job->JOBname ? " " : "", job->JOBname ? (char *) job->JOBname : "");
        if (pmi >= 0) {
            printf("phase margin %.4g deg at %.6g Hz (%s 0 dB crossing)",
                pmx[pmi].value, pmx[pmi].freq, pmx[pmi].aux != 0.0 ? "falling" : "rising");
            if (pmx[pmi].aux == 0.0)
                printf(" -- RISING |T| crossing, not the usual falling-crossing stability margin%s",
                    pmx[pmi].value > 180.0 ? " (a value over 180 deg here is not a real margin)" : "");
        } else {
            printf("phase margin: no 0 dB crossing of T");
        }
        if (npm > 1)
            printf(" (%d crossings found; the lowest margin is reported)", npm);
        printf("\n");

        printf("stb%s%s: ", job->JOBname ? " " : "", job->JOBname ? (char *) job->JOBname : "");
        if (gmi >= 0)
            printf("gain margin %.4g dB at %.6g Hz", gmx[gmi].value, gmx[gmi].freq);
        else
            printf("gain margin: no -180 deg crossing of T");
        if (ngm > 1)
            printf(" (%d crossings found; the lowest margin is reported)", ngm);
        printf("\n");

        /* card #476: exactly one of these three fires (disjoint arg0 ranges) --
         * mirrors tools/stb.py's describe() elif chain exactly. */
        if (low_f_posfb) {
            printf("stb%s%s: WARNING: positive feedback at low frequency: loop is unstable "
                "if |T|>1 there; phase/gain margins are not meaningful\n",
                job->JOBname ? " " : "", job->JOBname ? (char *) job->JOBname : "");
        } else if (low_f_posfb_low_gain) {
            printf("stb%s%s: note: positive feedback at low frequency, loop gain < 1\n",
                job->JOBname ? " " : "", job->JOBname ? (char *) job->JOBname : "");
        } else if (low_f_integrator) {
            printf("stb%s%s: note: integrator-like low-frequency behaviour\n",
                job->JOBname ? " " : "", job->JOBname ? (char *) job->JOBname : "");
        }

        /* card #475: dPM/dp, dfc/dp (implicit differentiation of the
         * crossing condition |T(fc)|=1 through the SAME bracket/log-f
         * interpolation the crossing search above used) and dGM/dp (a
         * plain chain-rule readout at the gain crossing's own fixed
         * phase point -- see the design comment above stb_sens_dot()).
         * "no crossing" (pmi/gmi < 0) reports NaN here too, same rule as
         * phaseMargin/gainMargin themselves: never a number for a
         * quantity that was not actually measured. */
        if (sens_requested && nsens > 0) {
            IFuid *suids = TMALLOC(IFuid, 3 * nsens);
            double *sout = TMALLOC(double, 3 * nsens);
            runDesc *splot = NULL;
            double snan = strtod("NaN", NULL);
            int kk;

            for (kk = 0; kk < nsens; kk++) {
                char nbuf[340];
                int ci;

                snprintf(nbuf, sizeof nbuf, "dPM_d%s", sensNames[kk]);
                for (ci = 0; nbuf[ci]; ci++) if (nbuf[ci] == ':') nbuf[ci] = '_';
                SPfrontEnd->IFnewUid(ckt, &suids[3 * kk + 0], NULL, nbuf, UID_OTHER, NULL);
                snprintf(nbuf, sizeof nbuf, "dfc_d%s", sensNames[kk]);
                for (ci = 0; nbuf[ci]; ci++) if (nbuf[ci] == ':') nbuf[ci] = '_';
                SPfrontEnd->IFnewUid(ckt, &suids[3 * kk + 1], NULL, nbuf, UID_OTHER, NULL);
                snprintf(nbuf, sizeof nbuf, "dGM_d%s", sensNames[kk]);
                for (ci = 0; nbuf[ci]; ci++) if (nbuf[ci] == ':') nbuf[ci] = '_';
                SPfrontEnd->IFnewUid(ckt, &suids[3 * kk + 2], NULL, nbuf, UID_OTHER, NULL);

                sout[3 * kk + 0] = snan;
                sout[3 * kk + 1] = snan;
                sout[3 * kk + 2] = snan;

                if (pmi >= 0) {
                    int idx = pmx[pmi].idx;
                    double t2v = pmx[pmi].t2;
                    double f0 = freqs[idx], f1 = freqs[idx + 1];
                    double x0 = log(f0), x1 = log(f1);
                    STBcplx Tc, dTc;
                    double Tcmag2;

                    Tc.re = Tarr[idx].re + t2v * (Tarr[idx + 1].re - Tarr[idx].re);
                    Tc.im = Tarr[idx].im + t2v * (Tarr[idx + 1].im - Tarr[idx].im);
                    dTc.re = dTdp_re[kk][idx] + t2v * (dTdp_re[kk][idx + 1] - dTdp_re[kk][idx]);
                    dTc.im = dTdp_im[kk][idx] + t2v * (dTdp_im[kk][idx + 1] - dTdp_im[kk][idx]);
                    Tcmag2 = Tc.re * Tc.re + Tc.im * Tc.im;

                    if (Tcmag2 > 0.0 && x1 != x0) {
                        /* d ln|T|/dp = Re(dT/dp / T) = Re(dT/dp * conj(T)) / |T|^2 */
                        double dGdp = (dTc.re * Tc.re + dTc.im * Tc.im) / Tcmag2;
                        double m0log = log(stb_cmag(Tarr[idx]));
                        double m1log = log(stb_cmag(Tarr[idx + 1]));
                        double dGdx = (m1log - m0log) / (x1 - x0);
                        /* d(arg T)/dp [rad] = Im(dT/dp / T) = Im(dT/dp*conj(T))/|T|^2 */
                        double dphidp_deg = (180.0 / M_PI) * (dTc.im * Tc.re - dTc.re * Tc.im) / Tcmag2;
                        double dphidx_deg = (phs[idx + 1] - phs[idx]) / (x1 - x0);
                        double fc = pmx[pmi].freq;

                        if (dGdx != 0.0) {
                            double dxdp = -dGdp / dGdx;
                            sout[3 * kk + 0] = dphidp_deg + dphidx_deg * dxdp;  /* dPM/dp */
                            sout[3 * kk + 1] = fc * dxdp;                       /* dfc/dp */
                        }
                    }
                }

                if (gmi >= 0) {
                    int idx = gmx[gmi].idx;
                    double t2v = gmx[gmi].t2;
                    double f0 = freqs[idx], f1 = freqs[idx + 1];
                    double x0 = log(f0), x1 = log(f1);
                    STBcplx Tc, dTc;
                    double Tcmag2;

                    Tc.re = Tarr[idx].re + t2v * (Tarr[idx + 1].re - Tarr[idx].re);
                    Tc.im = Tarr[idx].im + t2v * (Tarr[idx + 1].im - Tarr[idx].im);
                    dTc.re = dTdp_re[kk][idx] + t2v * (dTdp_re[kk][idx + 1] - dTdp_re[kk][idx]);
                    dTc.im = dTdp_im[kk][idx] + t2v * (dTdp_im[kk][idx + 1] - dTdp_im[kk][idx]);
                    Tcmag2 = Tc.re * Tc.re + Tc.im * Tc.im;

                    /* card #475 bugfix (found validating vs central FD:
                     * dGM/dp was wrong by O(1), far worse than dPM/dp's
                     * plain bracket-discretisation error): GM is read at
                     * the -180 deg phase crossing, and that crossing's
                     * OWN frequency shifts with p too -- same implicit-
                     * differentiation structure as fc above, just through
                     * the phase=-180 constraint instead of |T|=1. Missing
                     * this (the original version below only had the
                     * "evaluated at a FIXED f" partial term) is exactly
                     * "dropping half the chain rule", not a rounding
                     * artifact -- this file's own red-proof recipe
                     * (docs/sens.md) uses a REVERT of this exact fix as
                     * one of its two required red-proofs.
                     *   d(phase)/dp + d(phase)/df * df_gm/dp = 0
                     *   => df_gm/dp = -(d(phase)/dp) / (d(phase)/df)
                     *   dGM/dp = dGM/dp|_fixed_f + dGM/df|_fixed_f * df_gm/dp */
                    if (Tcmag2 > 0.0 && x1 != x0) {
                        double dGdp = (dTc.re * Tc.re + dTc.im * Tc.im) / Tcmag2;       /* d ln|T|/dp */
                        double m0log = log(stb_cmag(Tarr[idx]));
                        double m1log = log(stb_cmag(Tarr[idx + 1]));
                        double dGdx = (m1log - m0log) / (x1 - x0);                      /* d ln|T|/dx */
                        double dphidp_deg = (180.0 / M_PI) * (dTc.im * Tc.re - dTc.re * Tc.im) / Tcmag2;
                        double dphidx_deg = (phs[idx + 1] - phs[idx]) / (x1 - x0);
                        double fgm = gmx[gmi].freq;
                        double dGM_dp_fixed_f = -20.0 / log(10.0) * dGdp;
                        double dGM_df_fixed_f = -20.0 / log(10.0) * dGdx / fgm;
                        if (dphidx_deg != 0.0) {
                            double dfgm_dp = -dphidp_deg * fgm / dphidx_deg;
                            sout[3 * kk + 2] = dGM_dp_fixed_f + dGM_df_fixed_f * dfgm_dp;
                        }
                    }
                }
            }

            error = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob,
                (IFuid) "STB dMargins/dp", NULL, 0, 3 * nsens, suids, IF_REAL, &splot);
            if (!error) {
                IFvalue refval, outdata;
                refval.rValue = 0;
                outdata.v.numValue = 3 * nsens;
                outdata.v.vec.rVec = sout;
                SPfrontEnd->OUTpData(splot, &refval, &outdata);
                SPfrontEnd->OUTendPlot(splot);
            }
            tfree(suids);
            tfree(sout);
        }
    }

    if (freqs) tfree(freqs);
    if (Tarr) tfree(Tarr);
    if (LGarr) tfree(LGarr);
    if (pmx) tfree(pmx);
    if (gmx) tfree(gmx);
    if (phs) tfree(phs);

    /* card #475: sensitivity scratch, freed unconditionally -- every
     * pointer here is either NULL (sens not requested / nsens==0) or was
     * allocated exactly once, above. */
    if (nsens > 0) {
        int kk;
        for (kk = 0; kk < nsens; kk++) {
            tfree(sensNames[kk]);
            if (dTdp_re) tfree(dTdp_re[kk]);
            if (dTdp_im) tfree(dTdp_im[kk]);
        }
    }
    if (sensNames) tfree(sensNames);
    if (sensVecUids) tfree(sensVecUids);
    if (dTdp_re) tfree(dTdp_re);
    if (dTdp_im) tfree(dTdp_im);
    /* SMPdestroy() frees the matrix's own internals; the small SMPmatrix
     * wrapper itself is left un-freed here, matching cktsens.c's existing
     * DC/AC `.sens' (its own delta_Y, same TMALLOC+SMPdestroy pattern,
     * never tfree()s the wrapper either) -- not changed for this fork. */
    if (sensDeltaY) SMPdestroy(sensDeltaY);
    if (sens_lamv_re) tfree(sens_lamv_re);
    if (sens_lamv_im) tfree(sens_lamv_im);
    if (sens_lami_re) tfree(sens_lami_re);
    if (sens_lami_im) tfree(sens_lami_im);
    if (sens_x1_re) tfree(sens_x1_re);
    if (sens_x1_im) tfree(sens_x1_im);
    if (sens_x2_re) tfree(sens_x2_re);
    if (sens_x2_im) tfree(sens_x2_im);
    if (sens_dI) tfree(sens_dI);
    if (sens_dII) tfree(sens_dII);
    if (sens_mv_re) tfree(sens_mv_re);
    if (sens_mv_im) tfree(sens_mv_im);
    if (sensDxdp) {
        int kk;
        for (kk = 0; kk < nsens; kk++)
            if (sensDxdp[kk]) tfree(sensDxdp[kk]);
        tfree(sensDxdp);
    }
    if (sensH) tfree(sensH);
    if (sensRec) tfree(sensRec);
    if (sensD) tfree(sensD);
    if (sensNoAc) tfree(sensNoAc);
    if (sensJoint) tfree(sensJoint);
    if (sens_xop) tfree(sens_xop);
    if (sens_bufOld) tfree(sens_bufOld);
    if (sens_state0) tfree(sens_state0);
    if (sens_yp1_re) tfree(sens_yp1_re);
    if (sens_yp1_im) tfree(sens_yp1_im);
    if (sens_yp2_re) tfree(sens_yp2_re);
    if (sens_yp2_im) tfree(sens_yp2_im);

    return OK;
}
