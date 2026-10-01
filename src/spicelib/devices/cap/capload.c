/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: September 2003 Paolo Nenzi
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "capdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "capflat.h"

/* card #532 (round 15): set by CKTload (real KLU, MODETRAN) when the matrix
 * part of every capacitor stamp (m*ag0*C) has already been folded into the
 * matrix as Ax = Rbase + ag0*Cbase; CAPload then does only its state/rhs work. */
int CAPskipMatrix = 0;

void
CAPflatFree(CAPflat *fl)
{
    free(fl->pos); free(fl->neg); free(fl->q); free(fl->m); free(fl->cap);
    memset(fl, 0, sizeof(*fl));
}

/* The state/rhs half of CAPload for MODETRAN (cond1 false, trapezoidal order
 * 1 or 2) over the flattened instance arrays; the operations per instance are
 * those of CAPload + NIintegrate in the same order, so the state and rhs
 * results are bit-identical to the instance-list walk.  Returns 0 if the
 * ordinary CAPload has to be used instead. */
int
CAPloadFlat(CKTcircuit *ckt, const CAPflat *fl)
{
    const int *pos = fl->pos, *neg = fl->neg, *qq = fl->q;
    const double *mm = fl->m, *cc = fl->cap;
    double *s0 = ckt->CKTstate0, *s1 = ckt->CKTstate1, *rhs = ckt->CKTrhs;
    const double *old = ckt->CKTrhsOld;
    const double ag0 = ckt->CKTag[0], ag1 = ckt->CKTag[1];
    const int pred = (ckt->CKTmode & MODEINITPRED) != 0;
    const int itran = (ckt->CKTmode & MODEINITTRAN) != 0;
    const int n = fl->n;
    int k;
    if (!(ckt->CKTmode & MODETRAN) || ckt->CKTintegrateMethod != TRAPEZOIDAL ||
        (ckt->CKTorder != 1 && ckt->CKTorder != 2))
        return 0;
    if ((ckt->CKTmode & MODEDC) && (ckt->CKTmode & MODEINITJCT)) return 0;
    if ((ckt->CKTmode & MODEUIC) && itran) return 0;
    for (k = 0; k < n; k++) {
        const int q = qq[k], c = q + 1;
        double ceq;
        if (pred) {
            s0[q] = s1[q];
        } else {
            const double vcap = old[pos[k]] - old[neg[k]];
            s0[q] = cc[k] * vcap;
            if (itran) s1[q] = s0[q];
        }
        if (ckt->CKTorder == 1)
            s0[c] = ag0 * s0[q] + ag1 * s1[q];
        else
            s0[c] = - s1[c] * ag1 + ag0 * (s0[q] - s1[q]);
        ceq = s0[c] - ag0 * s0[q];
        if (itran) s1[c] = s0[c];
        rhs[pos[k]] -= mm[k] * ceq;
        rhs[neg[k]] += mm[k] * ceq;
    }
    return 1;
}

/* Matrix stamp of all capacitors with unit coefficient (m*capac, signs as in
 * CAPload) into base[0..nz).  Returns 0 (and the caller falls back to the
 * ordinary load) if any matrix pointer that CAPload would write is not inside
 * Ax[0..nz): such stamps (ground rows) cannot be reproduced as a base. */
int
CAPbuildBase(GENmodel *inModel, const double *Ax, unsigned int nz, double *base, CAPflat *fl)
{
    CAPmodel *model = (CAPmodel*)inModel;
    CAPinstance *here;
    int n = 0, k0 = 0;
    memset(base, 0, (size_t) nz * sizeof(double));
    memset(fl, 0, sizeof(*fl));
    for ( ; model != NULL; model = CAPnextModel(model))
        for (here = CAPinstances(model); here != NULL; here = CAPnextInstance(here))
            n++;
    fl->n = n;
    fl->pos = TMALLOC(int, n + 1); fl->neg = TMALLOC(int, n + 1); fl->q = TMALLOC(int, n + 1);
    fl->m = TMALLOC(double, n + 1); fl->cap = TMALLOC(double, n + 1);
    model = (CAPmodel*)inModel;
    for( ; model != NULL; model = CAPnextModel(model)) {
        for (here = CAPinstances(model); here != NULL; here = CAPnextInstance(here)) {
            double v = here->CAPm * here->CAPcapac;
            double *p[4];
            fl->pos[k0] = here->CAPposNode; fl->neg[k0] = here->CAPnegNode; fl->q[k0] = here->CAPqcap;
            fl->m[k0] = here->CAPm; fl->cap[k0] = here->CAPcapac; k0++;
            double sg[4] = { 1.0, 1.0, -1.0, -1.0 };
            int k;
            p[0] = here->CAPposPosPtr; p[1] = here->CAPnegNegPtr;
            p[2] = here->CAPposNegPtr; p[3] = here->CAPnegPosPtr;
            for (k = 0; k < 4; k++) {
                int a = (k == 0 || k == 2) ? here->CAPposNode : here->CAPnegNode;
                int b = (k == 0 || k == 3) ? here->CAPposNode : here->CAPnegNode;
                if (a <= 0 || b <= 0) continue;   /* trash-can stamp, never read */
                if (!p[k] || p[k] < Ax || p[k] >= Ax + nz) return 0;
                base[p[k] - Ax] += sg[k] * v;
            }
        }
    }
    return 1;
}

int
CAPload(GENmodel *inModel, CKTcircuit *ckt)
/* actually load the current capacitance value into the
 * sparse matrix previously provided
 */
{
    CAPmodel *model = (CAPmodel*)inModel;
    CAPinstance *here;
    int cond1;
    double vcap;
    double geq;
    double ceq;
    int error;
    double m;

    /* check if capacitors are in the circuit or are open circuited */
    if(ckt->CKTmode & (MODETRAN|MODEAC|MODETRANOP) ) {
        /* evaluate device independent analysis conditions */
        cond1=
            ( ( (ckt->CKTmode & MODEDC) &&
                (ckt->CKTmode & MODEINITJCT) )
              || ( ( ckt->CKTmode & MODEUIC) &&
                   ( ckt->CKTmode & MODEINITTRAN) ) ) ;
        /*  loop through all the capacitor models */
        for( ; model != NULL; model = CAPnextModel(model)) {

            /* loop through all the instances of the model */
            for (here = CAPinstances(model); here != NULL ;
                    here=CAPnextInstance(here)) {

                m = here->CAPm;

                if(cond1) {
                    vcap = here->CAPinitCond;
                } else {
                    vcap = *(ckt->CKTrhsOld+here->CAPposNode) -
                           *(ckt->CKTrhsOld+here->CAPnegNode) ;
                }
                if(ckt->CKTmode & (MODETRAN | MODEAC)) {
#ifndef PREDICTOR
                    if(ckt->CKTmode & MODEINITPRED) {
                        *(ckt->CKTstate0+here->CAPqcap) =
                            *(ckt->CKTstate1+here->CAPqcap);
                    } else { /* only const caps - no poly's */
#endif /* PREDICTOR */
                        *(ckt->CKTstate0+here->CAPqcap) = here->CAPcapac * vcap;
                        if((ckt->CKTmode & MODEINITTRAN)) {
                            *(ckt->CKTstate1+here->CAPqcap) =
                                *(ckt->CKTstate0+here->CAPqcap);
                        }
#ifndef PREDICTOR
                    }
#endif /* PREDICTOR */
                    error = NIintegrate(ckt,&geq,&ceq,here->CAPcapac,
                                        here->CAPqcap);
                    if(error) return(error);
                    if(ckt->CKTmode & MODEINITTRAN) {
                        *(ckt->CKTstate1+here->CAPccap) =
                            *(ckt->CKTstate0+here->CAPccap);
                    }
                    if (!CAPskipMatrix) {
                    *(here->CAPposPosPtr) += m * geq;
                    *(here->CAPnegNegPtr) += m * geq;
                    *(here->CAPposNegPtr) -= m * geq;
                    *(here->CAPnegPosPtr) -= m * geq;
                    }
                    *(ckt->CKTrhs+here->CAPposNode) -= m * ceq;
                    *(ckt->CKTrhs+here->CAPnegNode) += m * ceq;
                } else
                    *(ckt->CKTstate0+here->CAPqcap) = here->CAPcapac * vcap;
            }
        }
    }
    return(OK);
}

