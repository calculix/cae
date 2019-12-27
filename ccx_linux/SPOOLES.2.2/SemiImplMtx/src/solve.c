/*  solve.c  */

#include "../SemiImplMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to solve the linear system A X = B 
              using a semi-implicit factorization

   on return ---
     cpus[0] -- time to initialize working matrices
     cpus[1] -- time to load rhs 
     cpus[2] -- time for first solve with domains
     cpus[3] -- time to compute schur rhs
     cpus[4] -- time for schur solve
     cpus[5] -- time to compute domains' rhs
     cpus[6] -- time for second solve with domains
     cpus[7] -- time to store solution
     cpus[8] -- miscellaneous time 
     cpus[9] -- total time 

   return value --
      1 -- normal return
     -1 -- semimtx is NULL
     -2 -- mtxX is NULL
     -3 -- mtxB is NULL
     -4 -- mtxmanager is NULL
     -5 -- cpus is NULL

   created -- 98oct17, cca
   ------------------------------------------------
*/
int
SemiImplMtx_solve (
   SemiImplMtx     *semimtx,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   SubMtxManager   *mtxmanager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) {
DenseMtx   *T1, *T2 ;
double     t0, t1, t2 ;
double     alpha[2] = {-1.0, 0.0} ;
double     localcpus[6] ;
int        irow, ndomeqns, nrhs, nschureqns ;
int        *domCols, *domRows, *schurCols, *schurRows ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( semimtx == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_solve()"
           "\n semimtx is NULL\n") ;
   return(-1) ;
}
if ( mtxX == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_solve()"
           "\n mtxX is NULL\n") ;
   return(-2) ;
}
if ( mtxB == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_solve()"
           "\n mtxB is NULL\n") ;
   return(-3) ;
}
if ( mtxmanager == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_solve()"
           "\n mtxmanager is NULL\n") ;
   return(-4) ;
}
if ( cpus == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_solve()"
           "\n cpus is NULL\n") ;
   return(-5) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_solve()"
           "\n msglvl = %d, msgFile is NULL\n", msglvl) ;
   return(-6) ;
}
DVzero(10, cpus) ;
nrhs = mtxX->ncol ;
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(mtxB, "B", msgFile) ;
   DenseMtx_writeForMatlab(mtxX, "X", msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   create two DenseMtx objects to store the domain's 
   and schur complement's rhs and solution
   -------------------------------------------------
*/
MARKTIME(t1) ;
IV_sizeAndEntries(semimtx->domRowsIV, &ndomeqns, &domRows) ;
IV_sizeAndEntries(semimtx->schurRowsIV, &nschureqns, &schurRows) ;
T1 = DenseMtx_new() ;
DenseMtx_init(T1, mtxX->type, -1, -1, ndomeqns, nrhs, 1, ndomeqns) ;
T2 = DenseMtx_new() ;
DenseMtx_init(T2, mtxX->type, -1, -1, nschureqns, nrhs, 1, nschureqns) ;
MARKTIME(t2) ;
cpus[0] = t2 - t1 ;
/*
   ---------------------------------------
   load T1 with the domains' rhs and
   load T2 with the schur complement's rhs
   ---------------------------------------
*/
MARKTIME(t1) ;
for ( irow = 0 ; irow < ndomeqns ; irow++ ) {
   DenseMtx_copyRow(T1, irow, mtxB, domRows[irow]) ;
}
for ( irow = 0 ; irow < nschureqns ; irow++ ) {
   DenseMtx_copyRow(T2, irow, mtxB, schurRows[irow]) ;
}
MARKTIME(t2) ;
cpus[1] = t2 - t1 ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n T1 loaded with B1") ;
   DenseMtx_writeForHumanEye(T1, msgFile) ;
   fprintf(msgFile, "\n\n T2 loaded with B2") ;
   DenseMtx_writeForHumanEye(T2, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(T1, "B1", msgFile) ;
   DenseMtx_writeForMatlab(T2, "B2", msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   if symmetric then
      solve (U11^T+I)D11(I+U11) T1 = B1
   else if hermitian then
      solve (U11^H+I)D11(I+U11) T1 = B1
   else if nonsymmetric then
      solve (L11+I)D11(I+U11) T1 = B1
   end
   ------------------------------------
*/
MARKTIME(t1) ;
DVzero(6, localcpus) ;
FrontMtx_solve(semimtx->domainMtx, T1, T1, mtxmanager,
               localcpus, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[2] = t2 - t1 ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n T1 after first solve") ;
   DenseMtx_writeForHumanEye(T1, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(T1, "T1", msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------
   if symmetric then
      compute T2 := T2 - A12^T * T1
   else if hermitian then
      compute T2 := T2 - A12^H * T1
   else if nonsymmetric then
      compute T2 := T2 - A21 * T1
   end
   --------------------------------
*/
MARKTIME(t1) ;
switch ( semimtx->symmetryflag ) {
case SPOOLES_SYMMETRIC :
   InpMtx_nonsym_mmm_T(semimtx->A12, T2, alpha, T1) ;
   break ;
case SPOOLES_HERMITIAN :
   InpMtx_nonsym_mmm_H(semimtx->A12, T2, alpha, T1) ;
   break ;
case SPOOLES_NONSYMMETRIC :
   InpMtx_nonsym_mmm(semimtx->A21, T2, alpha, T1) ;
   break ;
}
MARKTIME(t2) ;
cpus[3] = t2 - t1 ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n schur rhs") ;
   DenseMtx_writeForHumanEye(T2, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(T2, "T2", msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   if symmetric then
      solve (U22^T+I)D22(I+U22) T2 = B2
   else if hermitian then
      solve (U22^H+I)D22(I+U22) T2 = B2
   else if nonsymmetric then
      solve (L22+I)D22(I+U22) T2 = B2
   end
   note: X2 = T2
   ------------------------------------
*/
MARKTIME(t1) ;
DVzero(6, localcpus) ;
FrontMtx_solve(semimtx->schurMtx, T2, T2, mtxmanager,
               localcpus, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[4] = t2 - t1 ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n schur solution") ;
   DenseMtx_writeForHumanEye(T2, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(T2, "X2", msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   compute T1 := B1 - A21 * T2
   ---------------------------
*/
MARKTIME(t1) ;
for ( irow = 0 ; irow < ndomeqns ; irow++ ) {
   DenseMtx_copyRow(T1, irow, mtxB, domRows[irow]) ;
}
InpMtx_nonsym_mmm(semimtx->A12, T1, alpha, T2) ;
MARKTIME(t2) ;
cpus[5] = t2 - t1 ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n domain rhs") ;
   DenseMtx_writeForHumanEye(T1, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(T1, "W1", msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   if symmetric then
      solve (U11^T+I)D11(I+U11) X1 = T1
   else if hermitian then
      solve (U11^H+I)D11(I+U11) X1 = T1
   else if nonsymmetric then
      solve (L11+I)D11(I+U11) X1 = T1
   end
   ------------------------------------
*/
MARKTIME(t1) ;
DVzero(6, localcpus) ;
FrontMtx_solve(semimtx->domainMtx, T1, T1, mtxmanager,
               localcpus, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[6] = t2 - t1 ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n domain solution") ;
   DenseMtx_writeForHumanEye(T1, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(T1, "X1", msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   load mtxX with the domains' solution
   and the schur complement's solution
   ------------------------------------
*/
MARKTIME(t1) ;
IV_sizeAndEntries(semimtx->domColsIV, &ndomeqns, &domCols) ;
IV_sizeAndEntries(semimtx->schurColsIV, &nschureqns, &schurCols) ;
for ( irow = 0 ; irow < ndomeqns ; irow++ ) {
   DenseMtx_copyRow(mtxX, domCols[irow], T1, irow) ;
}
for ( irow = 0 ; irow < nschureqns ; irow++ ) {
   DenseMtx_copyRow(mtxX, schurCols[irow], T2, irow) ;
}
MARKTIME(t2) ;
cpus[7] = t2 - t1 ;
if ( msglvl > 4 ) {
   DenseMtx_writeForMatlab(mtxX, "Xcomp", msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DenseMtx_free(T1) ;
DenseMtx_free(T2) ;

MARKTIME(t2) ;
cpus[9] = t2 - t0 ;
cpus[8] = cpus[9] - cpus[0] - cpus[1] - cpus[2] - cpus[3] 
        - cpus[4] - cpus[5] - cpus[6] - cpus[7] ;
/*
fprintf(msgFile, "\n solve cpus") ;
DVfprintf(msgFile, 10, cpus) ;
*/

return(1) ; }

/*--------------------------------------------------------------------*/
