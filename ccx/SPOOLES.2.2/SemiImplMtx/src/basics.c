/*  basics.c  */

#include "../SemiImplMtx.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98oct16, cca
   -----------------------
*/
SemiImplMtx *
SemiImplMtx_new ( 
   void 
) {
SemiImplMtx   *mtx ;

ALLOCATE(mtx, struct _SemiImplMtx, 1) ;
SemiImplMtx_setDefaultFields(mtx) ;

return(mtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   return code --
      1 -- normal return
     -1 -- mtx is NULL

   created -- 98oct16, cca
   -----------------------
*/
int
SemiImplMtx_setDefaultFields (
   SemiImplMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SemiImplMtx_setDefaultFields(%p)"
           "\n bad input", mtx) ;
   return(-1) ;
}
mtx->neqns        =   0  ;
mtx->type         = SPOOLES_REAL       ;
mtx->symmetryflag = SPOOLES_SYMMETRIC  ;
mtx->ndomeqns     =   0  ;
mtx->nschureqns   =   0  ;
mtx->domainMtx    = NULL ;
mtx->schurMtx     = NULL ;
mtx->A21          = NULL ;
mtx->A12          = NULL ;
mtx->domRowsIV    = NULL ;
mtx->schurRowsIV  = NULL ;
mtx->domColsIV    = NULL ;
mtx->schurColsIV  = NULL ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   return code --
      1 -- normal return
     -1 -- mtx is NULL

   created -- 98oct16, cca
   --------------------------------------------------
*/
int
SemiImplMtx_clearData ( 
   SemiImplMtx   *mtx 
) {
int   ieqn, neqns ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SemiImplMtx_clearData(%p)"
           "\n bad input\n", mtx) ;
   return(-1) ;
}
if ( (neqns = mtx->neqns) <= 0 ) {
   return(1) ;
}
if ( mtx->domainMtx != NULL ) {
   ETree   *etree      = mtx->domainMtx->frontETree ;
   IVL     *symbfacIVL = mtx->domainMtx->symbfacIVL ;
   FrontMtx_free(mtx->domainMtx) ;
   ETree_free(etree) ;
   IVL_free(symbfacIVL) ;
}
if ( mtx->schurMtx != NULL ) {
   ETree   *etree      = mtx->schurMtx->frontETree ;
   IVL     *symbfacIVL = mtx->schurMtx->symbfacIVL ;
   FrontMtx_free(mtx->schurMtx) ;
   ETree_free(etree) ;
   IVL_free(symbfacIVL) ;
}
if ( mtx->A12 != NULL ) {
   InpMtx_free(mtx->A12) ;
}
if ( mtx->domRowsIV != NULL ) {
   IV_free(mtx->domRowsIV) ;
}
if ( mtx->domColsIV != NULL ) {
   IV_free(mtx->domColsIV) ;
}
if ( mtx->schurRowsIV != NULL ) {
   IV_free(mtx->schurRowsIV) ;
}
if ( mtx->schurColsIV != NULL ) {
   IV_free(mtx->schurColsIV) ;
}
if ( mtx->symmetryflag == SPOOLES_NONSYMMETRIC ) {
   if ( mtx->A21 != NULL ) {
      InpMtx_free(mtx->A21) ;
   }
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
SemiImplMtx_setDefaultFields(mtx) ;

return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   return code --
      1 -- normal return
     -1 -- mtx is NULL

   created -- 98oct16, cca
   ------------------------------------------
*/
int
SemiImplMtx_free ( 
   SemiImplMtx   *mtx 
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SemiImplMtx_free(%p)"
           "\n bad input\n", mtx) ;
   return(-1) ;
}
SemiImplMtx_clearData(mtx) ;
FREE(mtx) ;

return(1) ; }

/*--------------------------------------------------------------------*/
