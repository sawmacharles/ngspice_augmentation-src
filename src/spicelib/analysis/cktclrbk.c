/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/*
 */

    /* CKTclrBreak(ckt)
     *   delete the first time from the breakpoint table for the given circuit
     */

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/sperror.h"
#include <string.h>
#include <stdlib.h>



int
CKTclrBreak(CKTcircuit *ckt)
{
    double *tmp;
    int j;
    static int inplace = -1;   /* card #555: NGSPICE_CLRBREAK_LEGACY=1 -> always the copying form */

    if (inplace < 0)
        inplace = !(getenv("NGSPICE_CLRBREAK_LEGACY") && atoi(getenv("NGSPICE_CLRBREAK_LEGACY")) != 0);

    if (inplace && ckt->CKTbreakSize > 256) {
        /* long table (a PSS period run with its recorded times forced: ~40000 entries): shift in place instead of
         * calloc + copy + free of the whole table at every accepted step (same contents, same size bookkeeping) */
        memmove(ckt->CKTbreaks, ckt->CKTbreaks + 1, (size_t) (ckt->CKTbreakSize - 1) * sizeof(double));
        ckt->CKTbreakSize--;
    } else
    if(ckt->CKTbreakSize >2) {
        tmp = TMALLOC(double, ckt->CKTbreakSize - 1);
        if(tmp == NULL) return(E_NOMEM);
        for(j=1;j<ckt->CKTbreakSize;j++) {
            tmp[j-1] = ckt->CKTbreaks[j];
        }
        FREE(ckt->CKTbreaks);
        ckt->CKTbreakSize--;
        ckt->CKTbreaks=tmp;
    } else {
        ckt->CKTbreaks[0] = ckt->CKTbreaks[1];
        ckt->CKTbreaks[1] = ckt->CKTfinalTime;
    }
    return(OK);
}
