/*  basics.c  */

#include "../FrontMtx.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may04, cca
   -----------------------
*/
FrontMtx *
FrontMtx_new ( 
   void 
) {
FrontMtx   *frontmtx ;

ALLOCATE(frontmtx, struct _FrontMtx, 1) ;
FrontMtx_setDefaultFields(frontmtx) ;

return(frontmtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may04, cca
   -----------------------
*/
void
FrontMtx_setDefaultFields (
   FrontMtx   *frontmtx
) {
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_setDefaultFields(%p)"
           "\n bad input", frontmtx) ;
   exit(-1) ;
}
frontmtx->nfront        =   0  ;
frontmtx->neqns         =   0  ;
frontmtx->type          =   SPOOLES_REAL          ;
frontmtx->symmetryflag  =   SPOOLES_SYMMETRIC     ;
frontmtx->sparsityflag  =   FRONTMTX_DENSE_FRONTS ;
frontmtx->pivotingflag  =   SPOOLES_NO_PIVOTING   ;
frontmtx->dataMode      =   FRONTMTX_1D_MODE      ;
frontmtx->nentD         =   0  ;
frontmtx->nentL         =   0  ;
frontmtx->nentU         =   0  ;
frontmtx->tree          = NULL ;
frontmtx->frontETree    = NULL ;
frontmtx->frontsizesIV  = NULL ;
frontmtx->symbfacIVL    = NULL ;
frontmtx->rowadjIVL     = NULL ;
frontmtx->coladjIVL     = NULL ;
frontmtx->lowerblockIVL = NULL ;
frontmtx->upperblockIVL = NULL ;
frontmtx->p_mtxDJJ      = NULL ;
frontmtx->p_mtxUJJ      = NULL ;
frontmtx->p_mtxUJN      = NULL ;
frontmtx->p_mtxLJJ      = NULL ;
frontmtx->p_mtxLNJ      = NULL ;
frontmtx->lowerhash     = NULL ;
frontmtx->upperhash     = NULL ;
frontmtx->lock          = NULL ;
frontmtx->nlocks        =   0  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may04, cca
   --------------------------------------------------
*/
void
FrontMtx_clearData ( 
   FrontMtx   *frontmtx 
) {
SubMtx   *mtx ;
int      ii, J, K, nadj, nfront ;
int      *adj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_clearData(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
/*
   ----------------------
   free the owned storage
   ----------------------
*/
if ( frontmtx->frontsizesIV != NULL ) {
   IV_free(frontmtx->frontsizesIV) ;
}
if ( frontmtx->rowadjIVL != NULL ) {
   IVL_free(frontmtx->rowadjIVL) ;
}
if ( frontmtx->coladjIVL != NULL ) {
   IVL_free(frontmtx->coladjIVL) ;
}
if ( frontmtx->p_mtxDJJ != NULL ) {
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (mtx = frontmtx->p_mtxDJJ[J]) != NULL ) {
         SubMtx_free(mtx) ;
      }
   }
   FREE(frontmtx->p_mtxDJJ) ;
}
if ( frontmtx->tree != NULL ) {
   if (  frontmtx->frontETree == NULL 
      || frontmtx->frontETree->tree != frontmtx->tree ) {
      Tree_free(frontmtx->tree) ;
   }
   frontmtx->tree = NULL ;
}
if ( frontmtx->dataMode == FRONTMTX_1D_MODE ) {
   if ( frontmtx->p_mtxUJJ != NULL ) {
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = frontmtx->p_mtxUJJ[J]) != NULL ) {
            SubMtx_free(mtx) ;
         }
      }
      FREE(frontmtx->p_mtxUJJ) ;
   }
   if ( frontmtx->p_mtxUJN != NULL ) {
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = frontmtx->p_mtxUJN[J]) != NULL ) {
            SubMtx_free(mtx) ;
         }
      }
      FREE(frontmtx->p_mtxUJN) ;
   }
   if ( frontmtx->p_mtxLJJ != NULL ) {
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = frontmtx->p_mtxLJJ[J]) != NULL ) {
            SubMtx_free(mtx) ;
         }
      }
      FREE(frontmtx->p_mtxLJJ) ;
   }
   if ( frontmtx->p_mtxLNJ != NULL ) {
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = frontmtx->p_mtxLNJ[J]) != NULL ) {
            SubMtx_free(mtx) ;
         }
      }
      FREE(frontmtx->p_mtxLNJ) ;
   }
} else if ( frontmtx->dataMode == FRONTMTX_2D_MODE ) {
   for ( J = 0 ; J < nfront ; J++ ) {
      FrontMtx_upperAdjFronts(frontmtx, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         K = adj[ii] ;
         if ( (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            SubMtx_free(mtx) ;
         }
      }
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      for ( J = 0 ; J < nfront ; J++ ) {
         FrontMtx_lowerAdjFronts(frontmtx, J, &nadj, &adj) ;
         for ( ii = 0 ; ii < nadj ; ii++ ) {
            K = adj[ii] ;
            if ( (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
               SubMtx_free(mtx) ;
            }
         }
      }
   }
   if ( frontmtx->lowerblockIVL != NULL ) {
      IVL_free(frontmtx->lowerblockIVL) ;
   }
   if ( frontmtx->upperblockIVL != NULL ) {
      IVL_free(frontmtx->upperblockIVL) ;
   }
   if ( frontmtx->lowerhash != NULL ) {
      I2Ohash_free(frontmtx->lowerhash) ;
   }
   if ( frontmtx->upperhash != NULL ) {
      I2Ohash_free(frontmtx->upperhash) ;
   }
}
if ( frontmtx->lock != NULL ) {
/*
   -------------------------
   destroy and free the lock
   -------------------------
*/
   Lock_free(frontmtx->lock) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
FrontMtx_setDefaultFields(frontmtx) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may04, cca
   ------------------------------------------
*/
void
FrontMtx_free ( 
   FrontMtx   *frontmtx 
) {
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_free(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
FrontMtx_clearData(frontmtx) ;
FREE(frontmtx) ;

return ; }

/*--------------------------------------------------------------------*/
