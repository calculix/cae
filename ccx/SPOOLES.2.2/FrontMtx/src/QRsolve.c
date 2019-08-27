/*  QRsolve.c  */

#include "../FrontMtx.h"
#include "../../timings.h"

#define MYDEBUG  0
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   minimize ||b - Ax||_2 by
   solving (U^T + I) * D * (I + U) X = A^T * B
      where A = QR = QD(I + U)
   by calling FrontMtx_solve()

   note: if A is square, mtxX and mtxB must not be the same

   mtxmanager -- object that manages working storage
   cpus -- vector of cpu time breakdowns
      cpus[0] -- set up solves
      cpus[1] -- fetch rhs and store solution
      cpus[2] -- forward solve
      cpus[3] -- diagonal solve
      cpus[4] -- backward solve
      cpus[5] -- total solve time
      cpus[6] -- time to compute A^T * B
      cpus[7] -- total time

   created  -- 97may27, dkw
   modified -- 98may28, cca
   --------------------------------------------------------
*/
void
FrontMtx_QR_solve (
   FrontMtx        *frontmtx,
   InpMtx          *mtxA,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   SubMtxManager   *mtxmanager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) {
double   t0, t1, t2 ;
double   alpha[2] ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( frontmtx == NULL || mtxA == NULL || mtxX == NULL || mtxB == NULL 
   || cpus == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_QR_solve(%p,%p,%p,%p,%p,%d,%p)"
           "\n bad input\n", 
           frontmtx, mtxA, mtxX, mtxB, cpus, msglvl, msgFile) ;
   exit(-1) ;
}
/*
   --------------------
   compute X := A^T * B
   --------------------
*/
MARKTIME(t1) ;
DenseMtx_zero(mtxX) ;
alpha[0] = 1.0 ; alpha[1] = 0.0 ;
if ( FRONTMTX_IS_REAL(frontmtx) ) {
   InpMtx_nonsym_mmm_T(mtxA, mtxX, alpha, mtxB) ;
} else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
   InpMtx_nonsym_mmm_H(mtxA, mtxX, alpha, mtxB) ;
}
MARKTIME(t2) ;
cpus[6] = t2 - t1 ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n B") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile)  ;
   fprintf(msgFile, "\n A^T * B") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile)  ;
   fflush(msgFile) ;
}
/*
   -----------------------------------
   solve (U^T + I) * D * (I + U) Z = X
   where Z overwrites X
   -----------------------------------
*/
MARKTIME(t1) ;
FrontMtx_solve(frontmtx, mtxX, mtxX, mtxmanager, 
               cpus, msglvl, msgFile) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n computed mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile)  ;
   fflush(msgFile) ;
}
MARKTIME(t2) ;
cpus[7] = t2 - t0 ;

return ; }

/*--------------------------------------------------------------------*/
