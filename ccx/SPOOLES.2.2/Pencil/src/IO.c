/*  IO.c  */

#include "../Pencil.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- to read in a Pencil object from a file

   input --

      inpmtxAfileName -- filename for A, must be *.inpmtxb or *.inpmtxf
      inpmtxBfileName -- filename for B, must be *.inpmtxb or *.inpmtxf

   return value -- 1 if success, 0 if failure

   created -- 97jul18, cca
   -------------------------------------------------------------------
*/
int
Pencil_readFromFiles ( 
   Pencil   *pencil, 
   char     *inpmtxAfileName, 
   char     *inpmtxBfileName 
) {
int    rc = 1 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  pencil == NULL || inpmtxAfileName == NULL 
   || inpmtxBfileName == NULL ) {
   fprintf(stderr, "\n error in Pencil_readFromFile(%p,%s,%s)"
           "\n bad input\n", pencil, inpmtxAfileName, inpmtxBfileName) ;
   return(0) ;
}
if ( strcmp(inpmtxAfileName, "none") != 0 ) {
   rc = InpMtx_readFromFile(pencil->inpmtxA, inpmtxAfileName) ;
   if ( rc != 1 ) {
      return(rc) ;
   } 
}
if ( strcmp(inpmtxBfileName, "none") != 0 ) {
   rc = InpMtx_readFromFile(pencil->inpmtxB, inpmtxBfileName) ;
   if ( rc != 1 ) {
      return(rc) ;
   } 
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to write a Pencil object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 97jul18, cca
   ----------------------------------------------------
*/
int
Pencil_writeForHumanEye ( 
   Pencil   *pencil, 
   FILE     *fp 
) {
if ( pencil == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Pencil_writeForHumanEye(%p,%p)"
           "\n bad input\n", pencil, fp) ;
   exit(-1) ;
}
/*
   ------------------------
   write out the statistics
   ------------------------
*/
Pencil_writeStats(pencil, fp) ;
if ( pencil->inpmtxA != NULL ) {
   fprintf(fp, "\n\n inpmtxA") ;
   InpMtx_writeForHumanEye(pencil->inpmtxA, fp) ;
}
if ( pencil->inpmtxB != NULL ) {
   fprintf(fp, "\n\n inpmtxB") ;
   InpMtx_writeForHumanEye(pencil->inpmtxB, fp) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to write out the statistics for the Pencil object

   return value -- 1 if success, 0 otherwise

   created -- 97jul18, cca
   -------------------------------------------------------------
*/
int
Pencil_writeStats ( 
   Pencil   *pencil, 
   FILE     *fp 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( pencil == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Pencil_writeStats(%p,%p)"
           "\n bad input\n", pencil, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n\n Pencil : matrix pencil object :") ;
if ( PENCIL_IS_REAL(pencil) ) {
   fprintf(fp, " real entries") ;
   fprintf(fp, "\n sigma = %20.12e ", pencil->sigma[0]) ;
} else if ( PENCIL_IS_COMPLEX(pencil) ) {
   fprintf(fp, " complex entries") ;
   fprintf(fp, "\n sigma = %20.12e + %20.12e*i",
           pencil->sigma[0], pencil->sigma[1]) ;
}
if ( pencil->inpmtxA != NULL ) {
   fprintf(fp, "\n\n inpmtxA") ;
   InpMtx_writeStats(pencil->inpmtxA, fp) ;
}
if ( pencil->inpmtxB != NULL ) {
   fprintf(fp, "\n\n inpmtxB") ;
   InpMtx_writeStats(pencil->inpmtxB, fp) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
