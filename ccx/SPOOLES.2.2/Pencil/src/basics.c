/*  basics.c  */

#include "../Pencil.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may02, cca
   -----------------------
*/
Pencil *
Pencil_new ( 
   void 
) {
Pencil   *pencil ;

ALLOCATE(pencil, struct _Pencil, 1) ;
Pencil_setDefaultFields(pencil) ;

return(pencil) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may02, cca
   -----------------------
*/
void
Pencil_setDefaultFields (
   Pencil   *pencil
) {
if ( pencil == NULL ) {
   fprintf(stderr, "\n fatal error in Pencil_setDefaultFields(%p)"
           "\n bad input", pencil) ;
   exit(-1) ;
}
pencil->type     = SPOOLES_REAL ;
pencil->symflag  = SPOOLES_SYMMETRIC ;
pencil->sigma[0] = 0.0  ;
pencil->sigma[1] = 0.0  ;
pencil->inpmtxA  = NULL ;
pencil->inpmtxB  = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may02, cca
   --------------------------------------------------
*/
void
Pencil_clearData ( 
   Pencil   *pencil 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( pencil == NULL ) {
   fprintf(stderr, "\n fatal error in Pencil_clearData(%p)"
           "\n bad input\n", pencil) ;
   exit(-1) ;
}
/*
   -----------------
   free the matrices
   -----------------
*/
if ( pencil->inpmtxA != NULL ) {
   InpMtx_free(pencil->inpmtxA) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_free(pencil->inpmtxB) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
Pencil_setDefaultFields(pencil) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may02, cca
   ------------------------------------------
*/
Pencil *
Pencil_free ( 
   Pencil   *pencil 
) {
if ( pencil == NULL ) {
   fprintf(stderr, "\n fatal error in Pencil_free(%p)"
           "\n bad input\n", pencil) ;
   exit(-1) ;
}
Pencil_clearData(pencil) ;
FREE(pencil) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
