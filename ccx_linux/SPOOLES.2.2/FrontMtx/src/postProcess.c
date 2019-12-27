/*  postProcess.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- post-process the factorization
      (1) permute row and column adjacency objects if necessary
      (2) permute lower and upper matrices if necessary
      (3) update the block adjacency objects if necessary
      (4) split the chevron submatrices into submatrices
          and make the submatrix indices local w.r.t their fronts

   created -- 98mar05, cca
   --------------------------------------------------------------
*/
void
FrontMtx_postProcess (
   FrontMtx   *frontmtx,
   int        msglvl,
   FILE       *msgFile
) {
int   nfront ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_postProcess(%p,%d,%p)"
           "\n bad input\n", frontmtx, msglvl, msgFile) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   IV   *colmapIV, *rowmapIV ;
/*
   -------------------------------
   permute the adjacency object(s)
   -------------------------------
*/
   FrontMtx_permuteUpperAdj(frontmtx, msglvl, msgFile) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n new column adjacency object") ;
      IVL_writeForHumanEye(frontmtx->coladjIVL, msgFile) ;
      fflush(msgFile) ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      FrontMtx_permuteLowerAdj(frontmtx, msglvl, msgFile) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n new row adjacency object") ;
         IVL_writeForHumanEye(frontmtx->rowadjIVL, msgFile) ;
         fflush(msgFile) ;
      }
   }
/*
   -------------------------------------------------------------
   permute the U_{J,bnd{J}} and L_{bnd{J},J} triangular matrices
   -------------------------------------------------------------
*/
   FrontMtx_permuteUpperMatrices(frontmtx, msglvl, msgFile) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      FrontMtx_permuteLowerMatrices(frontmtx, msglvl, msgFile) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n front factor matrix after pivoting") ;
     FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   }
/*
   -----------------------------------------------
   get the map from columns to owning fronts
   and create the new upper block adjacency object
   -----------------------------------------------
*/
   colmapIV = FrontMtx_colmapIV(frontmtx) ;
   frontmtx->upperblockIVL = FrontMtx_makeUpperBlockIVL(frontmtx,
                                                        colmapIV) ;
   IV_free(colmapIV) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      -------------------------------------------
      get the map from rows to owning fronts and
      create the new lower block adjacency object
      -------------------------------------------
*/
      rowmapIV = FrontMtx_rowmapIV(frontmtx) ;
      frontmtx->lowerblockIVL 
                   = FrontMtx_makeLowerBlockIVL(frontmtx, rowmapIV) ;
      IV_free(rowmapIV) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n new upper block adjacency object") ;
      IVL_writeForHumanEye(frontmtx->upperblockIVL, msgFile) ;
      if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
         fprintf(msgFile, "\n\n new lower block adjacency object") ;
         IVL_writeForHumanEye(frontmtx->lowerblockIVL, msgFile) ;
      }
      fflush(msgFile) ;
   }
} else {
/*
   ---------------------------------------
   get the upper block adjacency structure
   ---------------------------------------
*/
   IV *vtxToFrontIV = ETree_vtxToFrontIV(frontmtx->frontETree) ;
   frontmtx->upperblockIVL 
                 = FrontMtx_makeUpperBlockIVL(frontmtx, vtxToFrontIV) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      ---------------------------------------
      get the lower block adjacency structure
      ---------------------------------------
*/
      frontmtx->lowerblockIVL 
                 = FrontMtx_makeLowerBlockIVL(frontmtx, vtxToFrontIV) ;
   }
}
/*
   ------------------------
   allocate the hash tables
   ------------------------
*/
frontmtx->upperhash = I2Ohash_new() ;
I2Ohash_init(frontmtx->upperhash, nfront, nfront, nfront) ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   frontmtx->lowerhash = I2Ohash_new() ;
   I2Ohash_init(frontmtx->lowerhash, nfront, nfront, nfront) ;
} else {
   frontmtx->lowerhash = NULL ;
}
/*
   --------------------------------------------------------
   split the U_{J,bnd{J}} and L_{bnd{J},J} into submatrices
   put the U_{J,K} and L_{K,J} matrices into hash tables,
   free the p_mtx*[] vectors.
   --------------------------------------------------------
*/
FrontMtx_splitUpperMatrices(frontmtx, msglvl, msgFile) ;
FREE(frontmtx->p_mtxUJJ) ; frontmtx->p_mtxUJJ = NULL ;
FREE(frontmtx->p_mtxUJN) ; frontmtx->p_mtxUJN = NULL ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   FrontMtx_splitLowerMatrices(frontmtx, msglvl, msgFile) ;
   FREE(frontmtx->p_mtxLJJ) ; frontmtx->p_mtxLJJ = NULL ;
   FREE(frontmtx->p_mtxLNJ) ; frontmtx->p_mtxLNJ = NULL ;
}
frontmtx->dataMode    = FRONTMTX_2D_MODE ;

return ; }

/*--------------------------------------------------------------------*/
