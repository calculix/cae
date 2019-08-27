/*  init.c  */

#include "../Pencil.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   initialize the object

   created -- 98may02, cca
   -----------------------
*/
void
Pencil_init (
  Pencil   *pencil,
  int      type,
  int      symflag,
  InpMtx   *inpmtxA,
  double   sigma[],
  InpMtx   *inpmtxB
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  pencil == NULL || sigma == NULL ) {
   fprintf(stderr, "\n fatal error in Pencil_init(%p,%d,%d,%p,%p,%p)"
           "\n bad input\n", 
           pencil, type, symflag, inpmtxA, sigma, inpmtxB) ;
   exit(-1) ;
}
if ( !(type == SPOOLES_REAL || type == SPOOLES_COMPLEX) ) {
   fprintf(stderr, "\n fatal error in Pencil_init(%p,%d,%d,%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           pencil, type, symflag, inpmtxA, sigma, inpmtxB, type) ;
   exit(-1) ;
}
if ( !(symflag == SPOOLES_SYMMETRIC 
   ||  symflag == SPOOLES_HERMITIAN
   ||  symflag == SPOOLES_NONSYMMETRIC) ) { 
   fprintf(stderr, "\n fatal error in Pencil_init(%p,%d,%d,%p,%p,%p)"
           "\n bad symflag %d, must be SPOOLES_SYMMETRIC,"
           "\n SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC\n", 
           pencil, type, symflag, inpmtxA, sigma, inpmtxB, symflag) ;
   exit(-1) ;
}
/*
   ------------------------
   clear the data structure
   ------------------------
*/
Pencil_clearData(pencil) ;
/*
   --------------
   set the fields
   --------------
*/
pencil->type     = type     ;
pencil->symflag  = symflag  ;
pencil->inpmtxA  = inpmtxA  ;
pencil->sigma[0] = sigma[0] ;
pencil->sigma[1] = sigma[1] ;
pencil->inpmtxB  = inpmtxB  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   change the coordinate type

   created -- 98may02, cca
   --------------------------
*/
void
Pencil_changeCoordType (
   Pencil   *pencil,
   int       newType
) {
if ( pencil->inpmtxA != NULL ) {
   InpMtx_changeCoordType(pencil->inpmtxA, newType) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_changeCoordType(pencil->inpmtxB, newType) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   change the storage mode

   created -- 98may02, cca
   -----------------------
*/
void
Pencil_changeStorageMode (
   Pencil   *pencil,
   int       newMode
) {
if ( pencil->inpmtxA != NULL ) {
   InpMtx_changeStorageMode(pencil->inpmtxA, newMode) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_changeStorageMode(pencil->inpmtxB, newMode) ;
}
return ; }

/*--------------------------------------------------------------------*/
