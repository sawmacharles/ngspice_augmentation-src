/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/*
 */

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "tradefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"



int TRApssRestore = 0;   /* card #716: set by the shooting analyses while a run must start from the latched history */
int TRApssBreaks = 1;    /* card #716: replay the restored history's breakpoints (off in runs forced onto a recorded grid) */

/* card #716: the excitations of both ports at time t, by the three-point (quadratic Lagrange) interpolation of the delay table
 * at t - TD that the time-point prediction has always used (the code is the former MODEINITPRED branch of TRAload, unchanged);
 * the points and weights are kept in TRAip* for the time-domain small-signal recursion.  Returns 0 on a degenerate table (two
 * equal times: the excitations are left unchanged, as before). */
int
TRAinterp(TRAinstance *here, double t)
{
    double t1,t2,t3;
    double f1,f2,f3;
    int i;

    for(i=2;(i<here->TRAsizeDelay) &&
        (*(here->TRAdelays +3*i) <=
        (t-here->TRAtd));i++) {;/*loop does it*/}
    t1 = *(here->TRAdelays + (3*(i-2)));
    t2 = *(here->TRAdelays + (3*(i-1)));
    t3 = *(here->TRAdelays + (3*(i  )));
    if( (t2-t1)==0  || (t3-t2) == 0) return 0;
    f1 = (t - here->TRAtd - t2) *
         (t - here->TRAtd - t3) ;
    f2 = (t - here->TRAtd - t1) *
         (t - here->TRAtd - t3) ;
    f3 = (t - here->TRAtd - t1) *
         (t - here->TRAtd - t2) ;
    if((t2-t1)==0) { /* should never happen, but don't want
                      * to divide by zero, EVER... */
        f1=0;
        f2=0;
    } else {
        f1 /= (t1-t2);
        f2 /= (t2-t1);
    }
    if((t3-t2)==0) { /* should never happen, but don't want
                      * to divide by zero, EVER... */
        f2=0;
        f3=0;
    } else {
        f2 /= (t2-t3);
        f3 /= (t2-t3);
    }
    if((t3-t1)==0) { /* should never happen, but don't want
                      * to divide by zero, EVER... */
        f1=0;
        f2=0;
    } else {
        f1 /= (t1-t3);
        f3 /= (t1-t3);
    }
    here->TRAinput1 = f1 * *(here->TRAdelays + (3*(i-2))+1)
                    + f2 * *(here->TRAdelays + (3*(i-1))+1)
                    + f3 * *(here->TRAdelays + (3*(i  ))+1);
    here->TRAinput2 = f1 * *(here->TRAdelays + (3*(i-2))+2)
                    + f2 * *(here->TRAdelays + (3*(i-1))+2)
                    + f3 * *(here->TRAdelays + (3*(i  ))+2);
    here->TRAipT[0] = t1; here->TRAipT[1] = t2; here->TRAipT[2] = t3;
    here->TRAipF[0] = f1; here->TRAipF[1] = f2; here->TRAipF[2] = f3;
    here->TRAipOK = 1;
    return 1;
}

/*ARGSUSED*/
int
TRAload(GENmodel *inModel, CKTcircuit *ckt)
        /* actually load the current values into the 
         * sparse matrix previously provided 
         */
{
    TRAmodel *model = (TRAmodel *)inModel;
    TRAinstance *here;

    /*  loop through all the transmission line models */
    for( ; model != NULL; model = TRAnextModel(model)) {

        /* loop through all the instances of the model */
        for (here = TRAinstances(model); here != NULL ;
                here=TRAnextInstance(here)) {
            
            *(here->TRApos1Pos1Ptr) += here->TRAconduct;
            *(here->TRApos1Int1Ptr) -= here->TRAconduct;
            *(here->TRAneg1Ibr1Ptr) -= 1;
            *(here->TRApos2Pos2Ptr) += here->TRAconduct;
            *(here->TRAneg2Ibr2Ptr) -= 1;
            *(here->TRAint1Pos1Ptr) -= here->TRAconduct;
            *(here->TRAint1Int1Ptr) += here->TRAconduct;
            *(here->TRAint1Ibr1Ptr) += 1;
            *(here->TRAint2Int2Ptr) += here->TRAconduct;
            *(here->TRAint2Ibr2Ptr) += 1;
            *(here->TRAibr1Neg1Ptr) -= 1;
            *(here->TRAibr1Int1Ptr) += 1;
            *(here->TRAibr2Neg2Ptr) -= 1;
            *(here->TRAibr2Int2Ptr) += 1;
            *(here->TRApos2Int2Ptr) -= here->TRAconduct;
            *(here->TRAint2Pos2Ptr) -= here->TRAconduct;

            if(ckt->CKTmode & MODEDC) {
                *(here->TRAibr1Pos2Ptr) -= 1;
                *(here->TRAibr1Neg2Ptr) += 1;
                *(here->TRAibr1Ibr2Ptr) -= (1-ckt->CKTgmin)*here->TRAimped;
                *(here->TRAibr2Pos1Ptr) -= 1;
                *(here->TRAibr2Neg1Ptr) += 1;
                *(here->TRAibr2Ibr1Ptr) -= (1-ckt->CKTgmin)*here->TRAimped;
            } else {
                if (ckt->CKTmode & MODEINITTRAN) {
                    if(ckt->CKTmode & MODEUIC) {
                        here->TRAinput1 = here->TRAinitVolt2 + here->TRAinitCur2
                            * here->TRAimped;
                        here->TRAinput2 = here->TRAinitVolt1 + here->TRAinitCur1
                            * here->TRAimped;
                    } else {
                        here->TRAinput1 = 
                            ( *(ckt->CKTrhsOld+here->TRAposNode2)
                            - *(ckt->CKTrhsOld+here->TRAnegNode2) )
                            + ( *(ckt->CKTrhsOld+here->TRAbrEq2) 
                                *here->TRAimped);
                        here->TRAinput2 = 
                            ( *(ckt->CKTrhsOld+here->TRAposNode1)
                            - *(ckt->CKTrhsOld+here->TRAnegNode1) )
                            + ( *(ckt->CKTrhsOld+here->TRAbrEq1) 
                                *here->TRAimped);
                    }
                    *(here->TRAdelays  ) = -2*here->TRAtd;
                    *(here->TRAdelays+3) = -here->TRAtd;
                    *(here->TRAdelays+6) = 0;
                    *(here->TRAdelays+1) = *(here->TRAdelays +4) = 
                            *(here->TRAdelays+7) = here->TRAinput1;
                    *(here->TRAdelays+2) = *(here->TRAdelays +5) = 
                            *(here->TRAdelays+8) = here->TRAinput2;
                    here->TRAsizeDelay = 2;
                    here->TRAipOK = 0;
                    if (TRApssRestore && here->TRApssHistN >= 3) {
                        /* card #716: start from the latched periodic history */
                        int nh = here->TRApssHistN;
                        if (here->TRAallocDelay + 1 < nh) {
                            here->TRAallocDelay = nh + 4;
                            here->TRAdelays = TREALLOC(double, here->TRAdelays, (here->TRAallocDelay + 1) * 3);
                        }
                        memcpy(here->TRAdelays, here->TRApssHist, (size_t) nh * 3 * sizeof(double));
                        here->TRAsizeDelay = nh - 1;
                        (void) TRAinterp(here, ckt->CKTtime);
                        if (TRApssBreaks && !getenv("PSS_TL_NOBREAK")) {   /* the breakpoints TRAaccept() set when these points were accepted (a wave whose slope changed arrives
                             * at the other port TD later) lie in this run: set them again */
                            const double *hh = here->TRApssHist;
                            int j, q;
                            for (j = 2; j < nh; j++) {
                                double tj = hh[3 * j], tj1 = hh[3 * (j - 1)], tj2 = hh[3 * (j - 2)], when = tj1 + here->TRAtd;
                                int chg = 0;
                                if (tj <= tj1 || tj1 <= tj2 || when <= ckt->CKTtime || when >= ckt->CKTfinalTime)
                                    continue;
                                for (q = 1; q <= 2; q++) {
                                    double d1 = (hh[3 * j + q] - hh[3 * (j - 1) + q]) / (tj - tj1);
                                    double d2 = (hh[3 * (j - 1) + q] - hh[3 * (j - 2) + q]) / (tj1 - tj2);
                                    if (fabs(d1 - d2) >= here->TRAreltol * MAX(fabs(d1), fabs(d2)) + here->TRAabstol)
                                        chg = 1;
                                }
                                if (chg && CKTsetBreak(ckt, when))
                                    return E_INTERN;
                            }
                        }
                    }
                } else {
                    if(ckt->CKTmode & MODEINITPRED) {
                        if (!TRAinterp(here, ckt->CKTtime))
                            continue;
                    }
                }
            *(ckt->CKTrhs + here->TRAbrEq1) += here->TRAinput1;
            *(ckt->CKTrhs + here->TRAbrEq2) += here->TRAinput2;
            }
        }
    }
    return(OK);
}
