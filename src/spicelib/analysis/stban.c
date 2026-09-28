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

int
STBan(CKTcircuit *ckt, int restart)
{
    STBAN *job = (STBAN *) ckt->CKTcurJob;

    double freq, freqTol;
    int error, size, i, ignore;
    GENinstance *inst;
    int branch_probe, node_a, node_b;

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
    if (inst->GENmodPtr->GENmodType != CKTtypelook("Vsource")) {
        SPfrontEnd->IFerrorf(ERR_WARNING,
            "stb probe '%s' must be a voltage source (an existing 0V V-source "
            "in series with the loop-break wire)", job->STBprobeName);
        return E_NOTFOUND;
    }
    branch_probe = CKTfndBranch(ckt, job->STBprobeName);
    if (branch_probe == 0) {
        SPfrontEnd->IFerrorf(ERR_WARNING, "stb probe '%s': no branch equation found", job->STBprobeName);
        return E_NOTFOUND;
    }
    node_a = GENnode(inst)[0];
    node_b = GENnode(inst)[1];

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

        /* RHS #1: unit voltage injection in the probe branch equation */
        ckt->CKTrhs[branch_probe] = -1.0;
        ckt->CKTirhs[branch_probe] = 0.0;
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
