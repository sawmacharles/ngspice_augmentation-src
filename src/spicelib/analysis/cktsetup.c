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
    /* card #496, transient-ACCELERATION round 7: opt-in automatic
     * default-solver selection ('.options autoklu'), round 5/6's ranked
     * #1 candidate (research/transient-acceleration.md round 7).
     * Round 5 measured that ngspice's own default (always SPARSE unless
     * a deck explicitly says '.options klu') hits a >>100x pathology on
     * large, badly-ordered matrices (a_pex_buck: SPARSE did not finish
     * in 180s where KLU takes 44s) -- this defends a deck that never
     * thought to ask for KLU, without changing anything for a deck that
     * already states its solver explicitly (CKTkluExplicit) or is small
     * (below CKTautokluThreshold equations, where round 5 measured KLU's
     * own fixed per-call overhead can make it 6-34% SLOWER on tiny
     * matrices). Equation count is used as the threshold signal, not the
     * fill-in ratio round 5 also discusses: the fill-in-ratio version
     * needs a solver-specific symbolic ordering pass before the KLU/
     * SPARSE choice is even made, which touches the matrix-setup/CSC-
     * conversion path structurally and was explicitly judged too large
     * to implement and verify safely in one round (round 5 section 10.1.3,
     * round 6 section 11) -- this equation-count-only heuristic is the
     * deliberately smaller, safer subset of that candidate: it is
     * evaluated here, AFTER DEVsetup has fully stamped the matrix (so
     * SMPmatSize() below is exact, solver-agnostic -- see spsmp.c), and
     * BEFORE the COO->CSC/DEVbindCSC conversion just below, i.e. entirely
     * within this existing branch point and structurally identical
     * either way it resolves. */
    if (ckt->CKTautoklu && !ckt->CKTkluExplicit && !ckt->CKTmatrix->CKTkluMODE) {
        int CKTautokluNeqn = SMPmatSize(matrix);
        double CKTautokluThresh =
            (ckt->CKTautokluThreshold > 0) ? ckt->CKTautokluThreshold : 2000.0;
        if (CKTautokluNeqn >= CKTautokluThresh) {
            fprintf(stdout,
                "Note: .options autoklu switching to KLU (%d equations >= "
                "autoklu_threshold=%g)\n", CKTautokluNeqn, CKTautokluThresh);
            /* Keep ckt->CKTkluMODE and ckt->CKTmatrix->CKTkluMODE in sync:
             * several call sites elsewhere (cktic.c, cktacct.c, cktsens.c)
             * read ckt->CKTkluMODE directly rather than through the matrix,
             * on the assumption (true everywhere else in the tree) that
             * cktdojob.c/niinit.c already made the two agree -- checked by
             * grepping every CKTkluMODE site before adding this override,
             * not assumed. */
            ckt->CKTkluMODE = 1;
            ckt->CKTmatrix->CKTkluMODE = 1;
        }
    }

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
