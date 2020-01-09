/*  basics.c  */

#include "../InpMtx.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98jan28, cca
   -----------------------
*/
InpMtx *
InpMtx_new ( 
   void 
) {
InpMtx   *inpmtx ;

ALLOCATE(inpmtx, struct _InpMtx, 1) ;
InpMtx_setDefaultFields(inpmtx) ;

return(inpmtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_setDefaultFields (
   InpMtx   *inpmtx
) {
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_setDefaultFields(%p)"
           "\n bad input", inpmtx) ;
   exit(-1) ;
}
inpmtx->coordType      =   INPMTX_BY_ROWS  ;
inpmtx->storageMode    =   INPMTX_RAW_DATA ;
inpmtx->inputMode      =   SPOOLES_REAL    ;
inpmtx->maxnent        =   0  ;
inpmtx->nent           =   0  ;
inpmtx->resizeMultiple = 1.25 ;
inpmtx->maxnvector     =   0  ;
inpmtx->nvector        =   0  ;
IV_setDefaultFields(&inpmtx->ivec1IV)   ;
IV_setDefaultFields(&inpmtx->ivec2IV)   ;
DV_setDefaultFields(&inpmtx->dvecDV)    ;
IV_setDefaultFields(&inpmtx->vecidsIV)  ;
IV_setDefaultFields(&inpmtx->sizesIV)   ;
IV_setDefaultFields(&inpmtx->offsetsIV) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98jan28, cca
   --------------------------------------------------
*/
void
InpMtx_clearData ( 
   InpMtx   *inpmtx 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_clearData(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
/*
   -----------------------------------------------------
   free any storage held in the IV and DV vector objects
   -----------------------------------------------------
*/
IV_clearData(&inpmtx->ivec1IV)   ;
IV_clearData(&inpmtx->ivec2IV)   ;
DV_clearData(&inpmtx->dvecDV)    ;
IV_clearData(&inpmtx->vecidsIV)  ;
IV_clearData(&inpmtx->sizesIV)   ;
IV_clearData(&inpmtx->offsetsIV) ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
InpMtx_setDefaultFields(inpmtx) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98jan28, cca
   ------------------------------------------
*/
InpMtx *
InpMtx_free ( 
   InpMtx   *inpmtx 
) {
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_free(%p)"
           "\n bad input\n", inpmtx) ;
   exit(-1) ;
}
InpMtx_clearData(inpmtx) ;
FREE(inpmtx) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
