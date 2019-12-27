/*  util.c  */

#include "../SemiImplMtx.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   fill a statistics array for a semi-implicit factorization
     stats[0]  -- # of equations
     stats[1]  -- # of equations in the (1,1) block
     stats[2]  -- # of equations in the (2,2) block
     stats[3]  -- # of entries in L11
     stats[4]  -- # of entries in D11
     stats[5]  -- # of entries in U11
     stats[6]  -- # of entries in L22
     stats[7]  -- # of entries in D22
     stats[8]  -- # of entries in U22
     stats[9]  -- # of entries in A12
     stats[10] -- # of entries in A21
     stats[11] -- total # of entries
     stats[12] -- # of operations for a solve

   return value ---
      1 -- normal return
     -1 -- semimtx is NULL
     -2 -- stats is NULL

   created -- 98oct22, cca
   ---------------------------------------------------------
*/
int
SemiImplMtx_stats (
   SemiImplMtx   *semimtx,
   int           stats[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( semimtx == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_stats()"
           "\n semimtx is NULL\n") ;
   return(-1) ;
}
if ( stats == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_stats()"
           "\n stats is NULL\n") ;
   return(-2) ;
}
stats[0] = semimtx->neqns ;
stats[1] = semimtx->ndomeqns ;
stats[2] = semimtx->nschureqns ;
if ( semimtx->domainMtx != NULL ) {
   stats[3] = semimtx->domainMtx->nentL ;
   stats[4] = semimtx->domainMtx->nentD ;
   stats[5] = semimtx->domainMtx->nentU ;
} else {
   stats[3] = 0 ;
   stats[4] = 0 ;
   stats[5] = 0 ;
}
if ( semimtx->schurMtx != NULL ) {
   stats[6] = semimtx->schurMtx->nentL ;
   stats[7] = semimtx->schurMtx->nentD ;
   stats[8] = semimtx->schurMtx->nentU ;
} else {
   stats[6] = 0 ;
   stats[7] = 0 ;
   stats[8] = 0 ;
}
if ( semimtx->A12 != NULL ) {
   stats[9] = semimtx->A12->nent ;
} else {
   stats[9] = 0 ;
}
if ( semimtx->A21 != NULL ) {
   stats[10] = semimtx->A21->nent ;
} else {
   stats[10] = 0 ;
}
stats[11] = stats[3] + stats[4] + stats[5] + stats[6] + stats[7]
          + stats[8] + stats[9] + stats[10] ;
stats[12] = 0.0 ;
if ( semimtx->domainMtx != NULL ) {
   stats[12] += 2*FrontMtx_nSolveOps(semimtx->domainMtx) ;
}
if ( semimtx->schurMtx != NULL ) {
   stats[12] += FrontMtx_nSolveOps(semimtx->schurMtx) ;
}
if ( semimtx->A12 != NULL ) {
   if ( semimtx->type == SPOOLES_REAL ) {
      stats[12] += 2*semimtx->A12->nent ;
   } else if ( semimtx->type == SPOOLES_COMPLEX ) {
      stats[12] += 8*semimtx->A12->nent ;
   }
}
if ( semimtx->A21 != NULL ) {
   if ( semimtx->type == SPOOLES_REAL ) {
      stats[12] += 2*semimtx->A21->nent ;
   } else if ( semimtx->type == SPOOLES_COMPLEX ) {
      stats[12] += 8*semimtx->A21->nent ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
