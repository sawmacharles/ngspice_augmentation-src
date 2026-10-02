#ifndef ngspice_SMPDEFS_H
#define ngspice_SMPDEFS_H

/* Typedef removed by Francesco Lannutti (2012-02) to create the new SMPmatrix structure */
/*
typedef  struct MatrixFrame     SMPmatrix;
*/
typedef  struct MatrixFrame     MatrixFrame;
typedef  struct MatrixElement  *SMPelement;

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000  AlansFixes
**********/

#include <stdio.h>
#include <math.h>
#include "ngspice/complex.h"

#ifdef KLU
#include "ngspice/klu.h"
#include "ngspice/spmatrix.h"
#endif

/* SMPmatrix structure - Francesco Lannutti (2012-02) */
typedef struct sSMPmatrix {
    MatrixFrame *SPmatrix ;                /* pointer to sparse matrix */

#ifdef KLU
    KLUmatrix *SMPkluMatrix ;              /* KLU Pointer to the KLU Matrix Data Structure (only for CIDER, for the moment) */
    unsigned int CKTkluMODE:1 ;            /* KLU MODE parameter to enable KLU or not from the heuristic */
    #define CKTkluON 1                     /* KLU MODE ON definition */
    #define CKTkluOFF 0                    /* KLU MODE OFF definition */
    double CKTkluMemGrowFactor ;           /* KLU Memory Grow Factor - default = 1.2 */
#endif

} SMPmatrix ;


#ifdef KLU
void spDeterminant_KLU (SMPmatrix *, int *, double *, double *) ;
extern long KLUrgrowthFallbacks ;
extern int KLUrgrowthGuard ;    /* card #533: nonzero = refactor with tiny pivot growth is re-pivoted */
extern int KLUrefScaleForce ;   /* card #533: nonzero = klu_refactor keeps row scaling */
void SMPconvertCOOtoCSC (SMPmatrix *) ;

#ifdef CIDER
void SMPsolveKLUforCIDER (SMPmatrix *, double [], double [], double [], double []) ;
int SMPreorderKLUforCIDER (SMPmatrix *) ;
double *SMPmakeEltKLUforCIDER (SMPmatrix *, int, int) ;
void SMPclearKLUforCIDER (SMPmatrix *) ;
void SMPconvertCOOtoCSCKLUforCIDER (SMPmatrix *) ;
void SMPdestroyKLUforCIDER (SMPmatrix *) ;
int SMPnewMatrixKLUforCIDER (SMPmatrix *, int, unsigned int) ;
int SMPluFacKLUforCIDER (SMPmatrix *) ;
void SMPprintKLUforCIDER (SMPmatrix *, char *) ;
#endif

#else
int SMPaddElt (SMPmatrix *, int, int, double) ;
#endif

double * SMPmakeElt( SMPmatrix * , int , int );
void SMPcClear( SMPmatrix *);
void SMPclear( SMPmatrix *);
int SMPcLUfac( SMPmatrix *, double );
int SMPluFac( SMPmatrix *, double , double );
int SMPcReorder( SMPmatrix * , double , double , int *);
int SMPreorder( SMPmatrix * , double , double , double );
void SMPcaSolve(SMPmatrix *Matrix, double RHS[], double iRHS[],
		double Spare[], double iSpare[]);
void SMPcSolve( SMPmatrix *, double [], double [], double [], double []);
void SMPsolve( SMPmatrix *, double [], double []);
/* card #504: REAL transposed solve (A^T x = b, A the currently-factored
 * REAL matrix) -- the real-matrix analogue of SMPcaSolve (which is the
 * COMPLEX adjoint/transposed solve already used by docs/sens.md's card
 * #475 adjoint sensitivity). Dispatches Sparse (spSolveTransposed with
 * iRHS/iSolution NULL -- safe because Matrix->Complex is off whenever this
 * project's own transient/PAC code calls SMPsolve, so spSolveTransposed
 * takes its real-only branch, see spsolve.c) vs KLU (klu_tsolve, the real
 * counterpart of SMPcaSolve's klu_z_tsolve) exactly like SMPsolve/
 * SMPcaSolve already dispatch on Matrix->CKTkluMODE. */
void SMPtSolve( SMPmatrix *, double [], double []);
int SMPmatSize( SMPmatrix *);
int SMPnewMatrix( SMPmatrix *, int );
void SMPdestroy( SMPmatrix *);
int SMPpreOrder( SMPmatrix *);
void SMPprint( SMPmatrix * , char *);
void SMPprintRHS( SMPmatrix * , char *, double*, double*);
void SMPgetError( SMPmatrix *, int *, int *);
int SMPcProdDiag( SMPmatrix *, SPcomplex *, int *);
int SMPcDProd(SMPmatrix *Matrix, SPcomplex *pMantissa, int *pExponent);
SMPelement * SMPfindElt( SMPmatrix *, int , int , int );
int SMPcZeroCol(SMPmatrix *Matrix, int Col);
int SMPcAddCol(SMPmatrix *Matrix, int Accum_Col, int Addend_Col);
int SMPzeroRow(SMPmatrix *Matrix, int Row);
void SMPconstMult(SMPmatrix *, double);
void SMPmultiply(SMPmatrix *, double *, double *, double *, double *);
int SMPnumElements(SMPmatrix *);                         /* card #542 */
int SMPcopyElements(SMPmatrix *, int *, int *, double *); /* card #542: (ext row, ext col, value), unfactored real matrix */
int SMPsetElements(SMPmatrix *, int, const int *, const int *, const double *); /* card #542: clear + set the listed elements */

#ifdef CIDER
void SMPcSolveForCIDER (SMPmatrix *, double [], double [], double [], double []) ;
int SMPluFacForCIDER (SMPmatrix *) ;
int SMPnewMatrixForCIDER (SMPmatrix *, int, int) ;
void SMPsolveForCIDER (SMPmatrix *, double [], double []) ;
#endif

#endif

