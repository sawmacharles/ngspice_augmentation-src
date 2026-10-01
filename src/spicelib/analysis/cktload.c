/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/
/*
 */

/* CKTload(ckt)
 * this is a driver program to iterate through all the various
 * load functions provided for the circuit elements in the
 * given circuit
 */

#include "ngspice/ngspice.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"

#ifdef XSPICE
#include "ngspice/enh.h"
/* gtri - add - wbk - 11/26/90 - add include for MIF global data */
#include "ngspice/mif.h"
/* gtri - end - wbk - 11/26/90 */
#endif

#ifdef KLU
/* card #527 (round 14): resistor-base cache for the real KLU transient load.
 * A resistor's stamp is constant (RESconduct changes only in REStemp), so the
 * sum of all resistor stamps is snapshotted once and memcpy'd into the matrix
 * instead of re-running RESload (4 scattered read-modify-writes per resistor)
 * on every Newton iteration.  Keyed on RESgeneration (REStemp and the CSC
 * binding functions bump it), the matrix storage and the device list.  Only
 * MODETRAN (REScurrent, read only by flicker noise, is then not refreshed).
 * NGSPICE_RESBASE=0 disables it. */
extern unsigned long RESgeneration;
static int rb_enabled = -1;
static int rb_resIdx = -2;
static int rb_noinval = 0;   /* test-only red-proof: NGSPICE_RESBASE=2 ignores invalidation */
static struct { double *Ax; unsigned int nz; unsigned long gen; void *head; double *base; CKTcircuit *ckt; } rb;

/* card #532 (round 15): capacitor-base cache, on top of the resistor base.
 * A linear capacitor's matrix stamp in MODETRAN is m*ag0*C (NIintegrate:
 * geq = ag0*cap), so the summed unit stamp Cbase is cached and
 * Ax = Rbase + ag0*Cbase is formed in one vector pass; CAPload then skips its
 * matrix writes (state/rhs work unchanged).  Keyed on CAPgeneration (CAPtemp,
 * CAPparam, CSC binds), matrix storage, ckt and the device list head.  Falls back
 * if a capacitor stamp cannot be expressed as an offset into Ax.
 * NGSPICE_CAPBASE=0 disables (NGSPICE_RESBASE=0 disables both); =2 is the
 * test-only no-invalidation switch. */
extern unsigned long CAPgeneration;
extern int CAPskipMatrix;
#include "../devices/cap/capflat.h"
static int cb_enabled = -1;
static int cb_capIdx = -2;
static int cb_noinval = 0;
static struct { double *Ax; unsigned int nz; unsigned long gen; void *head; double *base; CKTcircuit *ckt; int ok; CAPflat fl; } cb;
static int cb_flat = -1;   /* NGSPICE_CAPFLAT=0: keep the instance-list CAPload (matrix writes skipped) */
#endif

static int ZeroNoncurRow(SMPmatrix *matrix, CKTnode *nodes, int rownum);

// #define PER_DEVICE_STATS
int
CKTload(CKTcircuit *ckt)
{
    int i;
    int size;
    double startTime, td0;
    CKTnode *node;
    int error;
    int skipRes;
    int skipCap = -1;
    int flatCap = 0;
#ifdef STEPDEBUG
    int noncon;
#endif /* STEPDEBUG */

#ifdef XSPICE
    /* gtri - begin - Put resistors to ground at all nodes */
    /*   SMPmatrix  *matrix; maschmann : deleted , because unused */

    double gshunt;
    int num_nodes;

    /* gtri - begin - Put resistors to ground at all nodes */
#endif

    td0 = SPfrontEnd->IFseconds();

#ifdef PER_DEVICE_STATS
    ckt->CKTstat->devTimes[DEVmaxnum] += SPfrontEnd->IFseconds() - td0;
    ckt->CKTstat->devCounts[DEVmaxnum]++;
#endif

    startTime = SPfrontEnd->IFseconds();
    size = SMPmatSize(ckt->CKTmatrix);
    for (i = 0; i <= size; i++) {
        ckt->CKTrhs[i] = 0;
    }
    skipRes = -1;
#ifdef KLU
    if (rb_enabled < 0) {
        const char *e = getenv("NGSPICE_RESBASE");
        rb_enabled = (e && atoi(e) == 0) ? 0 : 1;
        rb_noinval = (e && atoi(e) == 2);
    }
    if (rb_enabled && (ckt->CKTmode & MODETRAN) && ckt->CKTmatrix->CKTkluMODE &&
        ckt->CKTmatrix->SMPkluMatrix && !ckt->CKTmatrix->SMPkluMatrix->KLUmatrixIsComplex &&
        ckt->CKTmatrix->SMPkluMatrix->KLUmatrixAx && ckt->CKTmatrix->SMPkluMatrix->KLUmatrixNZ) {
        if (rb_resIdx == -2) {
            rb_resIdx = -1;
            for (i = 0; i < DEVmaxnum; i++)
                if (DEVices[i] && DEVices[i]->DEVpublic.name &&
                    strcmp(DEVices[i]->DEVpublic.name, "Resistor") == 0 && DEVices[i]->DEVload)
                    rb_resIdx = i;
        }
        if (rb_resIdx >= 0 && ckt->CKThead[rb_resIdx]) {
            double *Ax = ckt->CKTmatrix->SMPkluMatrix->KLUmatrixAx;
            unsigned int nz = ckt->CKTmatrix->SMPkluMatrix->KLUmatrixNZ;
            if (rb.base && rb.ckt == ckt && rb.Ax == Ax && rb.nz == nz &&
                (rb_noinval || rb.gen == RESgeneration) && rb.head == (void *) ckt->CKThead[rb_resIdx]) {
                memcpy(Ax, rb.base, (size_t) nz * sizeof(double));
                skipRes = rb_resIdx;
            } else {
                SMPclear(ckt->CKTmatrix);
                error = DEVices[rb_resIdx]->DEVload(ckt->CKThead[rb_resIdx], ckt);
                if (error) return(error);
                free(rb.base);
                rb.base = TMALLOC(double, nz);
                memcpy(rb.base, Ax, (size_t) nz * sizeof(double));
                rb.Ax = Ax; rb.nz = nz; rb.gen = RESgeneration; rb.ckt = ckt;
                rb.head = (void *) ckt->CKThead[rb_resIdx];
                skipRes = rb_resIdx;
            }
            if (cb_enabled < 0) {
                const char *e = getenv("NGSPICE_CAPBASE");
                cb_enabled = (e && atoi(e) == 0) ? 0 : 1;
                cb_noinval = (e && atoi(e) == 2);
            }
            if (cb_enabled && cb_capIdx == -2) {
                cb_capIdx = -1;
                for (i = 0; i < DEVmaxnum; i++)
                    if (DEVices[i] && DEVices[i]->DEVpublic.name &&
                        strcmp(DEVices[i]->DEVpublic.name, "Capacitor") == 0 && DEVices[i]->DEVload)
                        cb_capIdx = i;
            }
            if (cb_enabled && cb_capIdx >= 0 && ckt->CKThead[cb_capIdx]) {
                if (!(cb.base && cb.ckt == ckt && cb.Ax == Ax && cb.nz == nz &&
                      (cb_noinval || cb.gen == CAPgeneration) &&
                      cb.head == (void *) ckt->CKThead[cb_capIdx])) {
                    free(cb.base);
                    CAPflatFree(&cb.fl);
                    cb.base = TMALLOC(double, nz);
                    cb.ok = CAPbuildBase(ckt->CKThead[cb_capIdx], Ax, nz, cb.base, &cb.fl);
                    cb.Ax = Ax; cb.nz = nz; cb.gen = CAPgeneration; cb.ckt = ckt;
                    cb.head = (void *) ckt->CKThead[cb_capIdx];
                }
                if (cb.ok) {
                    const double *rbase = rb.base, *cbase = cb.base;
                    const double a0 = ckt->CKTag[0];
                    unsigned int k;
                    for (k = 0; k < nz; k++)
                        Ax[k] = rbase[k] + a0 * cbase[k];
                    skipCap = cb_capIdx;
                    CAPskipMatrix = 1;
                    if (cb_flat < 0) { const char *e = getenv("NGSPICE_CAPFLAT"); cb_flat = !(e && atoi(e) == 0); }
                    flatCap = cb_flat;
                }
            }
        }
    }
#endif
    if (skipRes < 0)
        SMPclear(ckt->CKTmatrix);
#ifdef STEPDEBUG
    noncon = ckt->CKTnoncon;
#endif /* STEPDEBUG */

    for (i = 0; i < DEVmaxnum; i++) {
        if (i == skipRes) continue;
        if (DEVices[i] && DEVices[i]->DEVload && ckt->CKThead[i]) {
#ifdef PER_DEVICE_STATS
            td0 = SPfrontEnd->IFseconds();
#endif
#ifdef KLU
            if (i == skipCap && flatCap && CAPloadFlat(ckt, &cb.fl))
                error = OK;
            else
#endif
            error = DEVices[i]->DEVload (ckt->CKThead[i], ckt);
#ifdef KLU
            if (i == skipCap) CAPskipMatrix = 0;
#endif
#ifdef PER_DEVICE_STATS
            ckt->CKTstat->devTimes[i] += SPfrontEnd->IFseconds() - td0;
            ckt->CKTstat->devCounts[i]++;
#endif
            if (ckt->CKTnoncon)
                ckt->CKTtroubleNode = 0;
#ifdef STEPDEBUG
            if (noncon != ckt->CKTnoncon) {
                printf("device type %s nonconvergence\n",
                       DEVices[i]->DEVpublic.name);
                noncon = ckt->CKTnoncon;
            }
#endif /* STEPDEBUG */
            if (error) {
#ifdef KLU
                CAPskipMatrix = 0;
#endif
                return(error);
            }
        }
    }
#ifdef KLU
    CAPskipMatrix = 0;
#endif


#ifdef XSPICE
    /* gtri - add - wbk - 11/26/90 - reset the MIF init flags */

    /* init is set by CKTinit and should be true only for first load call */
    g_mif_info.circuit.init = MIF_FALSE;

    /* anal_init is set by CKTdoJob and is true for first call */
    /* of a particular analysis type */
    g_mif_info.circuit.anal_init = MIF_FALSE;

    /* gtri - end - wbk - 11/26/90 */
#endif

    if (ckt->CKTmode & MODEDC) {
        /* consider doing nodeset & ic assignments */
        if (ckt->CKTmode & (MODEINITJCT | MODEINITFIX)) {
            /* do nodesets */
            for (node = ckt->CKTnodes; node; node = node->next) {
                if (node->nsGiven) {
                    if (ZeroNoncurRow(ckt->CKTmatrix, ckt->CKTnodes,
                                      node->number)) {
                        ckt->CKTrhs[node->number] = 1.0e10 * node->nodeset *
                                                      ckt->CKTsrcFact;
                        *(node->ptr) = 1e10;
                    } else {
                        ckt->CKTrhs[node->number] = node->nodeset *
                                                      ckt->CKTsrcFact;
                        *(node->ptr) = 1;
                    }
                    /* DAG: Original CIDER fix. If above fix doesn't work,
                     * revert to this.
                     */
                    /*
                     *  ckt->CKTrhs[node->number] += 1.0e10 * node->nodeset;
                     *  *(node->ptr) += 1.0e10;
                     */
                }
            }
        }
        if ((ckt->CKTmode & MODETRANOP) && (!(ckt->CKTmode & MODEUIC))) {
            for (node = ckt->CKTnodes; node; node = node->next) {
                if (node->icGiven) {
                    if (ZeroNoncurRow(ckt->CKTmatrix, ckt->CKTnodes,
                                      node->number)) {
                        /* Original code:
                         ckt->CKTrhs[node->number] += 1.0e10 * node->ic;
                        */
                        ckt->CKTrhs[node->number] = 1.0e10 * node->ic *
                                                      ckt->CKTsrcFact;
                        *(node->ptr) += 1.0e10;
                    } else {
                        /* Original code:
                          ckt->CKTrhs[node->number] = node->ic;
                        */
                        ckt->CKTrhs[node->number] = node->ic*ckt->CKTsrcFact; /* AlansFixes */
                        *(node->ptr) = 1;
                    }
                    /* DAG: Original CIDER fix. If above fix doesn't work,
                     * revert to this.
                     */
                    /*
                     *  ckt->CKTrhs[node->number] += 1.0e10 * node->ic;
                     *  *(node->ptr) += 1.0e10;
                     */
                }
            }
        }
    }
    /* SMPprint(ckt->CKTmatrix, stdout); if you want to debug, this is a
    good place to start ... */

    ckt->CKTstat->STATloadTime += SPfrontEnd->IFseconds()-startTime;
    return(OK);
}

static int
ZeroNoncurRow(SMPmatrix *matrix, CKTnode *nodes, int rownum)
{
    CKTnode     *n;
    double      *x;
    int         currents;

    currents = 0;
    for (n = nodes; n; n = n->next) {
        x = (double *) SMPfindElt(matrix, rownum, n->number, 0);
        if (x) {
            if (n->type == SP_CURRENT)
                currents = 1;
            else
                *x = 0.0;
        }
    }

    return currents;
}
