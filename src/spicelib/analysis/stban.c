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
 * pm_crossings / gm_crossings lists). */
typedef struct { double freq, value, aux; } CROSSING;

static CROSSING *
push_crossing(CROSSING *arr, int *n, int *cap, double freq, double value, double aux)
{
    if (*n >= *cap) {
        *cap = *cap ? *cap * 2 : 8;
        arr = TREALLOC(CROSSING, arr, *cap);
    }
    arr[*n].freq = freq;
    arr[*n].value = value;
    arr[*n].aux = aux;
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

    CROSSING *pmx = NULL, *gmx = NULL;
    int npm = 0, ngm = 0, capm = 0, capg = 0;

    /* ---- resolve the probe device ---- */
    inst = CKTfndDev(ckt, job->STBprobeName);
    if (!inst || inst->GENmodPtr->GENmodType < 0) {
        SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe '%s' not in circuit", job->STBprobeName);
        return E_NOTFOUND;
    }
    if (inst->GENmodPtr->GENmodType == CKTtypelook("Mos1")) {
        /* card #459: device mode. probe2/mode=dm|cm (loop-mode DM/CM,
         * card #457) is meaningless here -- refuse rather than silently
         * ignore it (rule: no silent numbers). */
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
    } else {
        SPfrontEnd->IFerrorf(ERR_WARNING,
            "stb probe '%s' must be a voltage source (an existing 0V V-source "
            "in series with the loop-break wire) or a supported active device "
            "instance (device mode: mos1 only in this fork)", job->STBprobeName);
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
    } else if (job->STBmodeGiven) {
        SPfrontEnd->IFerrorf(ERR_WARNING, "stb mode=dm|cm given without probe=<Vp>,<Vn>");
        return E_PARMVAL;
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

        freq = job->STBstartFreq;
    } else {
        freq = job->STBsaveFreq;
        job->STBsaveFreq = 0;
        error = SPfrontEnd->OUTpBeginPlot(NULL, NULL, NULL, NULL, 0, 666, NULL, 666, &stbPlot);
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

    /* card #459: NOW read gm/gmbs and the node indices off the MOS1
     * instance -- after CKTop()+CKTload() above (fresh start) or a prior
     * call's op point (resume/pause: MOS1gm lives on the persistent
     * instance struct, not this function's stack, so it is still valid
     * here without recomputing). Doing this any earlier reads a
     * never-computed (zero) gm on a fresh run -- see this file's comment
     * where device_mode is set, above. */
    if (device_mode && stb_mos1_probe(inst, &devp)) {
        SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe '%s': could not read device fields", job->STBprobeName);
        ckt->CKTcurrentAnalysis = 0;
        return E_NOTFOUND;
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
        v1.re = ckt->CKTrhs[node_b];  v1.im = ckt->CKTirhs[node_b];
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

        /* RHS #2: unit current injection into the probe's other terminal
         * (node_b); probe branch equation back to 0 (unperturbed) */
        for (i = 0; i <= size; i++) {
            ckt->CKTrhs[i] = 0.0;
            ckt->CKTirhs[i] = 0.0;
        }
        ckt->CKTrhs[node_b] += 1.0;
        /* localgnd=: reference the injected current to the given local
         * ground node instead of true ground (same two-terminal-source
         * pattern as the unmerged loop_analysis branch's LOOPan() used for
         * its optional refnode=). NOT independently validated numerically
         * (no reference deck here uses localgnd=) -- see the card #456
         * report's "what was not done". */
        if (job->STBlocalGndGiven && job->STBlocalGnd)
            ckt->CKTrhs[job->STBlocalGnd->number] -= 1.0;
        /* card #457 probe=Vp,Vn mode=dm|cm: drive branch n's current
         * injection in lockstep (same dmcm_sign as run 1, above). */
        if (job->STBprobeName2Given)
            ckt->CKTrhs[node_b2] += dmcm_sign * 1.0;
        SMPcSolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTirhs, ckt->CKTrhsSpare, ckt->CKTirhsSpare);
        v2.re = ckt->CKTrhs[node_b];  v2.im = ckt->CKTirhs[node_b];
        i2.re = -ckt->CKTrhs[branch_probe]; i2.im = -ckt->CKTirhs[branch_probe];

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
        }
        freqs[npts] = freq;
        Tarr[npts] = t;
        LGarr[npts] = lg;
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
    ckt->CKTcurrentAnalysis = 0;

    /* ---- margins: post-process the whole T(f) trace, mirroring
     * tools/stbcore.margins() exactly (log-f interpolation, unwrapped
     * phase, worst-of-several-crossings, "no crossing" is never a number) */
    {
        double *phs = TMALLOC(double, npts);
        double prevp = 0.0;
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
                pmx = push_crossing(pmx, &npm, &capm, fc, 180.0 + ph, 0.0);
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
                    gmx = push_crossing(gmx, &ngm, &capg, fc, -mag, target);
                }
                n++;
            }
        }
        tfree(phs);
    }

    {
        int pmi = -1, gmi = -1;
        double pmbest = 0.0, gmbest = 0.0;
        for (i = 0; i < npm; i++)
            if (pmi < 0 || pmx[i].value < pmbest) { pmi = i; pmbest = pmx[i].value; }
        for (i = 0; i < ngm; i++)
            if (gmi < 0 || gmx[i].value < gmbest) { gmi = i; gmbest = gmx[i].value; }

        {
            IFuid uids[4];
            double outv[4];
            runDesc *mplot = NULL;
            double nan = strtod("NaN", NULL);

            SPfrontEnd->IFnewUid(ckt, &uids[0], NULL, "phaseMargin", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[1], NULL, "phaseMarginFreq", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[2], NULL, "gainMargin", UID_OTHER, NULL);
            SPfrontEnd->IFnewUid(ckt, &uids[3], NULL, "gainMarginFreq", UID_OTHER, NULL);

            outv[0] = pmi >= 0 ? pmx[pmi].value : nan;
            outv[1] = pmi >= 0 ? pmx[pmi].freq : nan;
            outv[2] = gmi >= 0 ? gmx[gmi].value : nan;
            outv[3] = gmi >= 0 ? gmx[gmi].freq : nan;

            error = SPfrontEnd->OUTpBeginPlot(ckt, ckt->CKTcurJob,
                (IFuid) "STB Margins", NULL, 0, 4, uids, IF_REAL, &mplot);
            if (!error) {
                IFvalue refval, outdata;
                refval.rValue = 0;
                outdata.v.numValue = 4;
                outdata.v.vec.rVec = outv;
                SPfrontEnd->OUTpData(mplot, &refval, &outdata);
                SPfrontEnd->OUTendPlot(mplot);
            }
        }

        printf("stb%s%s: ", job->JOBname ? " " : "", job->JOBname ? (char *) job->JOBname : "");
        if (pmi >= 0)
            printf("phase margin %.4g deg at %.6g Hz", pmx[pmi].value, pmx[pmi].freq);
        else
            printf("phase margin: no 0 dB crossing of T");
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
    }

    if (freqs) tfree(freqs);
    if (Tarr) tfree(Tarr);
    if (LGarr) tfree(LGarr);
    if (pmx) tfree(pmx);
    if (gmx) tfree(gmx);

    return OK;
}
