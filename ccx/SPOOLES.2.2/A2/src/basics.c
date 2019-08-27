/*  basics.c  */

#include "../A2.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may01, cca
   -----------------------
*/
A2 *
A2_new ( 
   void 
) {
A2   *mtx ;

ALLOCATE(mtx, struct _A2, 1) ;
A2_setDefaultFields(mtx) ;

return(mtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may01, cca
   -----------------------
*/
void
A2_setDefaultFields (
   A2   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_setDefaultFields(%p)"
           "\n bad input", mtx) ;
   exit(-1) ;
}
mtx->type    =   SPOOLES_REAL  ;
mtx->n1      =   0  ;
mtx->inc1    =   0  ;
mtx->n2      =   0  ;
mtx->inc2    =   0  ;
mtx->nowned  =   0  ;
mtx->entries = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may01, cca
   --------------------------------------------------
*/
void
A2_clearData ( 
   A2   *mtx 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_clearData(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
/*
   -------------------------------------
   free the entries if present and owned
   -------------------------------------
*/
if ( mtx->nowned > 0 && mtx->entries != NULL ) {
   DVfree(mtx->entries) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
A2_setDefaultFields(mtx) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may01, cca
   ------------------------------------------
*/
void
A2_free ( 
   A2   *mtx 
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in A2_free(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
A2_clearData(mtx) ;
FREE(mtx) ;

return ; }

/*--------------------------------------------------------------------*/
