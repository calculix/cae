/*  IO.c  */

#include "../SemiImplMtx.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to write a SemiImplMtx to a file 
              in a human readable format

   return values ---
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- type is invalid
     -3 -- symmetry flag is invalid
     -4 -- fp is NULL

   created -- 98oct16, cca
   -------------------------------------------
*/
int
SemiImplMtx_writeForHumanEye (
   SemiImplMtx   *mtx,
   FILE          *fp
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_writeForHumanEye()"
           "\n mtx is NULL\n") ;
   return(-1) ;
}
switch ( mtx->type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n error in SemiImplMtx_writeForHumanEye()"
           "\n invalid type %d\n", mtx->type) ;
   return(-2) ;
   break ;
}
switch ( mtx->symmetryflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   fprintf(stderr, "\n error in SemiImplMtx_writeForHumanEye()"
           "\n invalid symmetry flag %d\n", mtx->symmetryflag) ;
   return(-3) ;
   break ;
}
if ( fp == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_writeForHumanEye()"
           "\n fp is NULL\n") ;
   return(-4) ;
}

fprintf(fp, "\n\n Semi-Implicit Matrix") ;
fprintf(fp, 
        "\n %d equations, %d in the domain, %d in the schur complement",
        mtx->neqns, mtx->ndomeqns, mtx->nschureqns) ;
switch ( mtx->type ) {
case SPOOLES_REAL :
   fprintf(fp, "\n real entries") ;
   break ;
case SPOOLES_COMPLEX :
   fprintf(fp, "\n complex entries") ;
   break ;
}
switch ( mtx->symmetryflag ) {
case SPOOLES_SYMMETRIC :
   fprintf(fp, ", symmetric matrix") ;
   break ;
case SPOOLES_HERMITIAN :
   fprintf(fp, ", Hermitian matrix") ;
   break ;
case SPOOLES_NONSYMMETRIC :
   fprintf(fp, ", nonsymmetric matrix") ;
   break ;
}
if ( mtx->domColsIV != NULL ) {
   fprintf(fp, "\n\n domain columns") ;
   IV_writeForHumanEye(mtx->domColsIV, fp) ;
}
if ( mtx->symmetryflag == SPOOLES_NONSYMMETRIC ) {
   if ( mtx->domRowsIV != NULL ) {
      fprintf(fp, "\n\n domain rows") ;
      IV_writeForHumanEye(mtx->domRowsIV, fp) ;
   }
}
if ( mtx->schurColsIV != NULL ) {
   fprintf(fp, "\n\n schur complement columns") ;
   IV_writeForHumanEye(mtx->schurColsIV, fp) ;
}
if ( mtx->symmetryflag == SPOOLES_NONSYMMETRIC ) {
   if ( mtx->schurRowsIV != NULL ) {
      fprintf(fp, "\n\n schur complement rows") ;
      IV_writeForHumanEye(mtx->schurRowsIV, fp) ;
   }
}
if ( mtx->domainMtx != NULL ) {
   fprintf(fp, "\n\n domain FrontMtx object") ;
   FrontMtx_writeForHumanEye(mtx->domainMtx, fp) ;
}
if ( mtx->schurMtx != NULL ) {
   fprintf(fp, "\n\n schur complement FrontMtx object") ;
   FrontMtx_writeForHumanEye(mtx->schurMtx, fp) ;
}
if ( mtx->A12 != NULL ) {
   fprintf(fp, "\n\n original (1,2) matrix") ;
   InpMtx_writeForHumanEye(mtx->A12, fp) ;
}
if ( mtx->symmetryflag == SPOOLES_NONSYMMETRIC ) {
   if ( mtx->A21 != NULL ) {
      fprintf(fp, "\n\n original (2,1) matrix") ;
      InpMtx_writeForHumanEye(mtx->A21, fp) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
