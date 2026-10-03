/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

    /* CKTsetup(ckt)
     * this is a driver program to iterate through all the various
     * setup functions provided for the circuit elements in the
     * given circuit
     */

#include "ngspice/ngspice.h"
#define PSS_AUTOKLU_MIN_DEFAULT 128
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/fteext.h"

#ifdef XSPICE
#include "ngspice/enh.h"
#endif

#ifdef USE_OMP
#include <omp.h>
#include "ngspice/cpextern.h"
#endif

#define CKALLOC(var,size,type) \
    if(size && ((var = TMALLOC(type, size)) == NULL)){\
            return(E_NOMEM);\
}

int
CKTsetup(CKTcircuit *ckt)
{
    int i;
    int error;
#ifdef USE_OMP
    int nthreads = 2;
#endif

    SMPmatrix *matrix;

    if (!ckt->CKThead) {
        fprintf(stderr, "Error: No model list found, device setup not possible!\n");
        if (ft_stricterror)
            controlled_exit(EXIT_BAD);
        return E_PANIC;
    }
    if (!DEVices) {
        fprintf(stderr, "Error: No device list found, device setup not possible!\n");
        if (ft_stricterror)
            controlled_exit(EXIT_BAD);
        return E_PANIC;
    }

    ckt->CKTnumStates=0;

#ifdef WANT_SENSE2
    if(ckt->CKTsenInfo){
        error = CKTsenSetup(ckt);
        if (error)
            return(error);
    }
#endif

    if (ckt->CKTisSetup)
        return E_NOCHANGE;

#ifdef KLU
    /* card #496, transient-ACCELERATION round 7: opt-in automatic
     * default-solver selection ('.options autoklu'), round 5/6's ranked
     * #1 candidate (research/transient-acceleration.md round 7).
     * Round 5 measured that ngspice's own default (always SPARSE unless
     * a deck explicitly says '.options klu') hits a >>100x pathology on
     * large, badly-ordered matrices -- this defends a deck that never
     * thought to ask for KLU, without changing anything for a deck that
     * already states its solver explicitly (CKTkluExplicit) or is small
     * (below CKTautokluThreshold equations, where round 5 measured KLU's
     * own fixed per-call overhead can make it 6-34% SLOWER on tiny
     * matrices).
     *
     * MUST run here, before NIinit()/DEVsetup(), not after (a first
     * implementation that decided after DEVsetup, right before the
     * existing COO->CSC conversion below, SEGFAULTED -- found by running
     * it, not by inspection: klusmp.c's own SMPmakeElt() branches on
     * Matrix->CKTkluMODE at EVERY CALL, so every device's matrix-element
     * request during DEVsetup already commits to the SPARSE (spGetElement)
     * or KLU (COO linked-list) storage path for that element; deciding
     * afterwards leaves the KLU COO list empty and SMPconvertCOOtoCSC()
     * converts nothing, binding garbage. So the decision has to be made
     * before matrix stamping starts, when the only fixed-cost signal
     * available is ckt->CKTlastNode -- the running node-number counter
     * CKTmkNode() maintains as the circuit is PARSED, already valid here,
     * long before CKTsetup() is ever called (grepped: cktnewn.c). This is
     * an approximation of the final equation count (device-internal
     * nodes added during DEVsetup are not yet counted), not the fill-in
     * ratio round 5/6 also discuss (which needs a symbolic ordering pass
     * that itself needs the matrix built first -- circular, and the
     * reason round 5/6 judged that fuller version too large for one
     * round); close enough for a coarse, defensive, opt-in threshold. */
    if (ckt->CKTautoklu && !ckt->CKTkluExplicit && !ckt->CKTkluMODE
            && ckt->CKTlastNode) {
        double CKTautokluThresh =
            (ckt->CKTautokluThreshold > 0) ? ckt->CKTautokluThreshold : 2000.0;
        int CKTautokluNodes = ckt->CKTlastNode->number;
        if ((double)CKTautokluNodes >= CKTautokluThresh) {
            fprintf(stdout,
                "Note: .options autoklu switching to KLU (%d nodes >= "
                "autoklu_threshold=%g)\n", CKTautokluNodes, CKTautokluThresh);
            ckt->CKTkluMODE = 1;
        }
    }

    /* card #558, round 17: periodic (shooting) analyses pick KLU on their
     * own above a node-count rule (research/transient-acceleration.md sec. 22).
     * Same constraints as autoklu: before NIinit/DEVsetup, never over an
     * explicit '.options klu'/'.options sparse'.  Off: NGSPICE_PSS_AUTOKLU=0;
     * threshold: NGSPICE_PSS_AUTOKLU_MIN (nodes, default PSS_AUTOKLU_MIN_DEFAULT). */
    if (ckt->CKTperiodicShoot && !ckt->CKTkluExplicit && !ckt->CKTkluMODE
            && ckt->CKTlastNode) {
        const char *e = getenv("NGSPICE_PSS_AUTOKLU");
        if (!(e && e[0] == '0')) {
            const char *m = getenv("NGSPICE_PSS_AUTOKLU_MIN");
            double thr = (m && m[0]) ? atof(m) : (double) PSS_AUTOKLU_MIN_DEFAULT;
            int nn = ckt->CKTlastNode->number;
            if (thr > 0 && (double) nn >= thr) {
                fprintf(stdout,
                    "Note: periodic analysis switching to KLU (%d nodes >= %g; "
                    ".options sparse or NGSPICE_PSS_AUTOKLU=0 keeps SPARSE)\n", nn, thr);
                ckt->CKTkluMODE = 1;
            }
        }
    }
#endif

    error = NIinit(ckt);
    if (error)
        return(error);

    ckt->CKTisSetup = 1;

    matrix = ckt->CKTmatrix;

#ifdef USE_OMP
    if (!cp_getvar("num_threads", CP_NUM, &nthreads, 0))
        nthreads = 2;

    omp_set_num_threads(nthreads);
/*    if (nthreads == 1)
      printf("OpenMP: %d thread is requested in ngspice\n", nthreads);
    else
      printf("OpenMP: %d threads are requested in ngspice\n", nthreads);*/
#endif

#ifdef HAS_PROGREP
    SetAnalyse("Device Setup", 0);
#endif

    /* preserve CKTlastNode before invoking DEVsetup()
     * so we can check for incomplete CKTdltNNum() invocations
     * during DEVunsetup() causing an erronous circuit matrix
     *   when reinvoking CKTsetup()
     */
    ckt->prev_CKTlastNode = ckt->CKTlastNode;

    for (i=0;i<DEVmaxnum;i++) {
        if ( DEVices[i] && DEVices[i]->DEVsetup && ckt->CKThead[i] ) {
            error = DEVices[i]->DEVsetup (matrix, ckt->CKThead[i], ckt,
                    &ckt->CKTnumStates);
            if(error) return(error);
        }
    }

#ifdef KLU
    if (ckt->CKTmatrix->CKTkluMODE)
    {
        fprintf (stdout, "Using KLU as Direct Linear Solver\n") ;

        /* Convert the COO Storage to CSC for KLU and Fill the Binding Table */
        SMPconvertCOOtoCSC (matrix) ;

        /* Assign the KLU Pointers */
        for (i = 0 ; i < DEVmaxnum ; i++)
            if (DEVices [i] && DEVices [i]->DEVbindCSC && ckt->CKThead [i])
                DEVices [i]->DEVbindCSC (ckt->CKThead [i], ckt) ;

    } else {
        fprintf (stdout, "Using SPARSE 1.3 as Direct Linear Solver\n") ;
    }
#endif

    for(i=0;i<=MAX(2,ckt->CKTmaxOrder)+1;i++) { /* dctran needs 3 states as minimum */
        CKALLOC(ckt->CKTstates[i],ckt->CKTnumStates,double);
    }
#ifdef WANT_SENSE2
    if(ckt->CKTsenInfo){
        /* to allocate memory to sensitivity structures if
         * it is not done before */

        error = NIsenReinit(ckt);
        if(error) return(error);
    }
#endif
    if(ckt->CKTniState & NIUNINITIALIZED) {
        error = NIreinit(ckt);
        if(error) return(error);
    }

    return(OK);
}

int
CKTunsetup(CKTcircuit *ckt)
{
    int i, error, e2;
    CKTnode *node;

    error = OK;
    if (!ckt->CKTisSetup)
        return OK;

    for(i=0;i<=ckt->CKTmaxOrder+1;i++) {
        tfree(ckt->CKTstates[i]);
    }

    /* added by HT 050802*/
    for(node=ckt->CKTnodes;node;node=node->next){
        if(node->icGiven || node->nsGiven) {
            node->ptr=NULL;
        }
    }

    for (i=0;i<DEVmaxnum;i++) {
        if ( DEVices[i] && DEVices[i]->DEVunsetup && ckt->CKThead[i] ) {
            e2 = DEVices[i]->DEVunsetup (ckt->CKThead[i], ckt);
            if (!error && e2)
                error = e2;
        }
    }

    if (ckt->prev_CKTlastNode != ckt->CKTlastNode) {
        fprintf(stderr, "Internal Error: incomplete CKTunsetup(), this will cause serious problems, please report this issue !\n");
        controlled_exit(EXIT_FAILURE);
    }
    ckt->prev_CKTlastNode = NULL;

    ckt->CKTisSetup = 0;
    if(error) return(error);

    NIdestroy(ckt);
    /*
    if (ckt->CKTmatrix)
        SMPdestroy(ckt->CKTmatrix);
    ckt->CKTmatrix = NULL;
    */

    return OK;
}
