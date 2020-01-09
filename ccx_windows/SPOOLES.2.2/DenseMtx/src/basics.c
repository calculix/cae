/*  basics.c  */

#include "../DenseMtx.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may02, cca
   -----------------------
*/
DenseMtx *
DenseMtx_new ( 
   void 
) {
DenseMtx   *mtx ;

ALLOCATE(mtx, struct _DenseMtx, 1) ;
DenseMtx_setDefaultFields(mtx) ;

return(mtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may02, cca
   -----------------------
*/
void
DenseMtx_setDefaultFields (
   DenseMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setDefaultFields(%p)"
           "\n bad input", mtx) ;
   exit(-1) ;
}
mtx->type    =  SPOOLES_REAL ;
mtx->rowid   =  -1  ;
mtx->colid   =  -1  ;
mtx->nrow    =   0  ;
mtx->ncol    =   0  ;
mtx->inc1    =   0  ;
mtx->inc2    =   0  ;
mtx->rowind  = NULL ;
mtx->colind  = NULL ;
mtx->entries = NULL ;
DV_setDefaultFields(&mtx->wrkDV) ;
mtx->next    = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may02, cca
   --------------------------------------------------
*/
void
DenseMtx_clearData ( 
   DenseMtx   *mtx 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_clearData(%p)"
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
DenseMtx_setDefaultFields(mtx) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may02, cca
   ------------------------------------------
*/
void
DenseMtx_free ( 
   DenseMtx   *mtx 
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_free(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
DenseMtx_clearData(mtx) ;
FREE(mtx) ;

return ; }

/*--------------------------------------------------------------------*/
