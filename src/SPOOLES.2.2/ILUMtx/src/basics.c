/*  basics.c  */

#include "../ILUMtx.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98oct03, cca
   -----------------------
*/
ILUMtx *
ILUMtx_new ( 
   void 
) {
ILUMtx   *mtx ;

ALLOCATE(mtx, struct _ILUMtx, 1) ;
ILUMtx_setDefaultFields(mtx) ;

return(mtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   return code --
      1 -- normal return
     -1 -- mtx is NULL

   created -- 98oct03, cca
   -----------------------
*/
int
ILUMtx_setDefaultFields (
   ILUMtx   *mtx
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in ILUMtx_setDefaultFields(%p)"
           "\n bad input", mtx) ;
   return(-1) ;
}
mtx->neqns        =   0  ;
mtx->type         = SPOOLES_REAL       ;
mtx->symmetryflag = SPOOLES_SYMMETRIC  ;
mtx->UstorageMode = SPOOLES_BY_ROWS    ;
mtx->LstorageMode = SPOOLES_BY_COLUMNS ;
mtx->sizesL       = NULL ;
mtx->sizesU       = NULL ;
mtx->p_indL       = NULL ;
mtx->p_indU       = NULL ;
mtx->entD         = NULL ;
mtx->p_entU       = NULL ;
mtx->p_entU       = NULL ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   return code --
      1 -- normal return
     -1 -- mtx is NULL

   created -- 98oct03, cca
   --------------------------------------------------
*/
int
ILUMtx_clearData ( 
   ILUMtx   *mtx 
) {
int   ieqn, neqns ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in ILUMtx_clearData(%p)"
           "\n bad input\n", mtx) ;
   return(-1) ;
}
if ( (neqns = mtx->neqns) < 0 ) {
   return(1) ;
}
if ( mtx->entD != NULL ) {
   DVfree(mtx->entD) ;
}
IVfree(mtx->sizesU) ;
for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
   if ( mtx->p_indU[ieqn] != NULL ) {
      IVfree(mtx->p_indU[ieqn]) ;
   }
   if ( mtx->p_entU[ieqn] != NULL ) {
      DVfree(mtx->p_entU[ieqn]) ;
   }
}
PIVfree(mtx->p_indU) ;
PDVfree(mtx->p_entU) ;
if ( mtx->symmetryflag == SPOOLES_NONSYMMETRIC ) {
   IVfree(mtx->sizesL) ;
   for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
      if ( mtx->p_indL[ieqn] != NULL ) {
         IVfree(mtx->p_indL[ieqn]) ;
      }
      if ( mtx->p_entL[ieqn] != NULL ) {
         DVfree(mtx->p_entL[ieqn]) ;
      }
   }
   PIVfree(mtx->p_indL) ;
   PDVfree(mtx->p_entL) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
ILUMtx_setDefaultFields(mtx) ;

return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   return code --
      1 -- normal return
     -1 -- mtx is NULL

   created -- 98oct03, cca
   ------------------------------------------
*/
int
ILUMtx_free ( 
   ILUMtx   *mtx 
) {
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in ILUMtx_free(%p)"
           "\n bad input\n", mtx) ;
   return(-1) ;
}
ILUMtx_clearData(mtx) ;
FREE(mtx) ;

return(1) ; }

/*--------------------------------------------------------------------*/
