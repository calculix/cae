/*  basics.c  */

#include "../SubMtx.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may01, cca
   -----------------------
*/
SubMtx *
SubMtx_new ( 
   void 
) {
SubMtx   *mtx ;

ALLOCATE(mtx, struct _SubMtx, 1) ;
SubMtx_setDefaultFields(mtx) ;

return(mtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may01, cca
   -----------------------
*/
void
SubMtx_setDefaultFields (
   SubMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_setDefaultFields(%p)"
           "\n bad input", mtx) ;
   exit(-1) ;
}
mtx->type    = SPOOLES_REAL ;
mtx->mode    = SUBMTX_DENSE_COLUMNS ;
mtx->rowid   =  -1  ;
mtx->colid   =  -1  ;
mtx->nrow    =   0  ;
mtx->ncol    =   0  ;
mtx->nent    =   0  ;
mtx->next    = NULL ;
DV_setDefaultFields(&mtx->wrkDV) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may01, cca
   --------------------------------------------------
*/
void
SubMtx_clearData ( 
   SubMtx   *mtx 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_clearData(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_clearData(&mtx->wrkDV) ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
SubMtx_setDefaultFields(mtx) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may01, cca
   ------------------------------------------
*/
void
SubMtx_free ( 
   SubMtx   *mtx 
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_free(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
SubMtx_clearData(mtx) ;
FREE(mtx) ;

return ; }

/*--------------------------------------------------------------------*/
