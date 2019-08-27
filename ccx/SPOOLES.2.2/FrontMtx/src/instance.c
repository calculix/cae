/*  instance.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- return the number of fronts
 
   created -- 98may04, cca
   --------------------------------------
*/
int
FrontMtx_nfront (
   FrontMtx   *frontmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_nfront(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
return(frontmtx->nfront) ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- return the number of equations
 
   created -- 98may04, cca
   -----------------------------------------
*/
int
FrontMtx_neqns (
   FrontMtx   *frontmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_neqns(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
return(frontmtx->neqns) ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- return a pointer to the front Tree object
 
   created -- 98may04, cca
   ----------------------------------------------------
*/
Tree *
FrontMtx_frontTree (
   FrontMtx   *frontmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_frontTree(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
return(frontmtx->tree) ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   simple method to return the dimensions of front J and the number 
   of bytes necessary for the Chv object to hold the front.

   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_initialFrontDimensions (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnD,
   int         *pnL,
   int         *pnU,
   int         *pnbytes
) {
int   nbytes, nD, nL, nU ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || J < 0 || J >= frontmtx->nfront
   || pnD == NULL || pnL == NULL || pnU == NULL || pnbytes == NULL ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_initialFrontDimensions()"
           "\n frontmtx = %p, J = %d, pnD = %p, "
           "pnL = %p, pnU = %p, pnbytes = %p",
           frontmtx, J, pnD, pnL, pnU, pnbytes) ;
   exit(-1) ;
}
switch ( frontmtx->type ) {
case SPOOLES_REAL :
   switch ( frontmtx->symmetryflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, 
              "\n fatal error in FrontMtx_initialFrontDimensions()"
              "\n real type, must be symmetric or nonsymmetric\n") ;
      exit(-1) ;
      break ;
   }
  break ;
case SPOOLES_COMPLEX :
   switch ( frontmtx->symmetryflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
   case SPOOLES_NONSYMMETRIC :
      break ;
      fprintf(stderr, 
              "\n fatal error in FrontMtx_initialFrontDimensions()"
              "\n complex type, must be symmetric,"
              "\n hermitian or nonsymmetric\n") ;
      exit(-1) ;
      break ;
   }
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in FrontMtx_initialFrontDimensions()"
           "\n bad type, must be real or complex") ;
   exit(-1) ;
   break ;
}
nD = frontmtx->frontETree->nodwghtsIV->vec[J] ;
nL = nU = frontmtx->frontETree->bndwghtsIV->vec[J] ;
nbytes = Chv_nbytesNeeded(nD, nL, nU, 
                          frontmtx->type, frontmtx->symmetryflag) ;
*pnD = nD ;
*pnL = nL ;
*pnU = nU ;
*pnbytes = nbytes ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   return the number of internal rows and columns in front J

   created -- 98may04, cca
   ---------------------------------------------------------
*/
int
FrontMtx_frontSize (
   FrontMtx   *frontmtx,
   int         J
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || frontmtx->frontsizesIV == NULL 
   || J < 0 || J >= frontmtx->nfront ) {
   fprintf(stderr, "\n fatal error in FrontMtx_frontSize(%p,%d)"
           "\n bad input\n", frontmtx, J) ;
   exit(-1) ;
}
return(IV_entry(frontmtx->frontsizesIV, J)) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   set the number of internal rows and columns in front J

   created -- 98may04, cca
   ------------------------------------------------------
*/
void
FrontMtx_setFrontSize (
   FrontMtx   *frontmtx,
   int         J,
   int         size
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || frontmtx->frontsizesIV == NULL 
   || J < 0 || J >= frontmtx->nfront || size < 0 ) {
   fprintf(stderr, "\n fatal error in FrontMtx_setFrontSize(%p,%d,%d)"
           "\n bad input\n", frontmtx, J, size) ;
   exit(-1) ;
}
IV_setEntry(frontmtx->frontsizesIV, J, size) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   fill *pncol with the number of columns and 
   *pcolind with a pointer to the column indices

   created -- 98may04, cca
   ---------------------------------------------
*/
void
FrontMtx_columnIndices (
   FrontMtx   *frontmtx,
   int         J,
   int         *pncol,
   int         **pcolind
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || J < 0 || J >= frontmtx->nfront 
   || pncol == NULL || pcolind == NULL ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_columnIndices(%p,%d,%p,%p)"
           "\n bad input\n", frontmtx, J, pncol, pcolind) ;
   exit(-1) ;
}
if ( ! FRONTMTX_IS_PIVOTING(frontmtx) ) {
   IVL_listAndSize(frontmtx->symbfacIVL, J, pncol, pcolind) ;
} else {
   IVL_listAndSize(frontmtx->coladjIVL, J, pncol, pcolind) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   fill *pnrow with the number of rows and 
   *prowind with a pointer to the rows indices

   created -- 98may04, cca
   -------------------------------------------
*/
void
FrontMtx_rowIndices (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnrow,
   int         **prowind
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || J < 0 || J >= frontmtx->nfront 
   || pnrow == NULL || prowind == NULL ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_rowIndices(%p,%d,%p,%p)"
           "\n bad input\n", frontmtx, J, pnrow, prowind) ;
   exit(-1) ;
}
if ( ! FRONTMTX_IS_PIVOTING(frontmtx) ) {
   IVL_listAndSize(frontmtx->symbfacIVL, J, pnrow, prowind) ;
} else if ( FRONTMTX_IS_SYMMETRIC(frontmtx) 
         || FRONTMTX_IS_HERMITIAN(frontmtx) ) {
   IVL_listAndSize(frontmtx->coladjIVL, J, pnrow, prowind) ;
} else if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   IVL_listAndSize(frontmtx->rowadjIVL, J, pnrow, prowind) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- return a pointer to the (J,J) diagonal submatrix
 
   created -- 98may04, cca
   -----------------------------------------------------------
*/
SubMtx *
FrontMtx_diagMtx (
   FrontMtx   *frontmtx,
   int         J
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || J < 0 || J >= frontmtx->nfront ) {
   fprintf(stderr, "\n fatal error in FrontMtx_diagMtx(%p,%d)"
           "\n bad input\n", frontmtx, J) ;
   exit(-1) ;
}
return(frontmtx->p_mtxDJJ[J]) ; }
 
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- return a pointer to the (J,K) upper submatrix
 
   created -- 98may04, cca
   --------------------------------------------------------
*/
SubMtx *
FrontMtx_upperMtx (
   FrontMtx   *frontmtx,
   int         J,
   int         K
) {
int    rc ;
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL
   || J < 0 || J >= frontmtx->nfront
   || K < J || K > frontmtx->nfront ) {
   fprintf(stderr, "\n fatal error in FrontMtx_upperMtx(%p,%d,%d)"
           "\n bad input\n", frontmtx, J, K) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   if ( K == frontmtx->nfront ) {
      mtx = frontmtx->p_mtxUJN[J] ;
   } else if ( K == J ) {
      mtx = frontmtx->p_mtxUJJ[J] ;
   }
} else if ( frontmtx->upperhash == NULL ) {
   mtx = NULL ;
} else {
   rc = I2Ohash_locate(frontmtx->upperhash, J, K, (void *) &mtx) ;
   if ( rc == 0 ) {
      mtx = NULL ;
   }
}
return(mtx) ; }
 
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- return a pointer to the (K,J) lower submatrix
 
   created -- 98may04, cca
   --------------------------------------------------------
*/
SubMtx *
FrontMtx_lowerMtx (
   FrontMtx   *frontmtx,
   int         K,
   int         J
) {
int    rc ;
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL
   || J < 0 || J >= frontmtx->nfront
   || K < J || K > frontmtx->nfront ) {
   fprintf(stderr, "\n fatal error in FrontMtx_lowerMtx(%p,%d,%d)"
           "\n bad input\n", frontmtx, K, J) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   if ( K == frontmtx->nfront ) {
      mtx = frontmtx->p_mtxLNJ[J] ;
   } else if ( K == J ) {
      mtx = frontmtx->p_mtxLJJ[J] ;
   }
} else if ( frontmtx->lowerhash == NULL ) {
   mtx = NULL ;
} else {
   rc = I2Ohash_locate(frontmtx->lowerhash, K, J, (void *) &mtx) ;
   if ( rc == 0 ) {
      mtx = NULL ;
   }
}
return(mtx) ; }
 
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- fill *pnadj with the number of fronts K
              such that L_{K,J} != 0 and *padj with a
              pointer to a list of those fronts
 
   created -- 98may04, cca
   --------------------------------------------------
*/
void
FrontMtx_lowerAdjFronts (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnadj,
   int         **padj
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL
   || J < 0 || J >= frontmtx->nfront
   || pnadj == NULL || padj == NULL ) {
   fprintf(stderr,
           "\n fatal error in FrontMtx_lowerAdjFronts(%p,%d,%p,%p)"
          "\n bad input\n", frontmtx, J, pnadj, padj) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_lowerAdjFronts()"
           "\n data mode is 1-D, not 2-D\n") ;
   exit(-1) ;
} else if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   IVL_listAndSize(frontmtx->lowerblockIVL, J, pnadj, padj) ;
} else {
   IVL_listAndSize(frontmtx->upperblockIVL, J, pnadj, padj) ;
}
return ; }
 
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- fill *pnadj with the number of fronts K
              such that U_{J,K} != 0 and *padj with a
              pointer to a list of those fronts
 
   created -- 98may04, cca
   --------------------------------------------------
*/
void
FrontMtx_upperAdjFronts (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnadj,
   int         **padj
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL
   || J < 0 || J >= frontmtx->nfront
   || pnadj == NULL || padj == NULL ) {
   fprintf(stderr,
           "\n fatal error in FrontMtx_upperAdjFronts(%p,%d,%p,%p)"
          "\n bad input\n", frontmtx, J, pnadj, padj) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_upperAdjFronts()"
           "\n data mode is 1, not 2\n") ;
   exit(-1) ;
}
IVL_listAndSize(frontmtx->upperblockIVL, J, pnadj, padj) ;
 
return ; }
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- return the number of nonzero L_{K,J} blocks
 
   created -- 98may04, cca
   ------------------------------------------------------
*/
int
FrontMtx_nLowerBlocks (
   FrontMtx   *frontmtx
) {
int   nblocks ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_nLowerBlocks(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_nLowerBlocks()"
           "\n data mode is 1, not 2\n") ;
   exit(-1) ;
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   nblocks = frontmtx->lowerblockIVL->tsize ;
} else {
   nblocks = frontmtx->upperblockIVL->tsize ;
}
return(nblocks) ; }
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- return the number of nonzero U_{K,J} blocks
 
   created -- 98may04, cca
   ------------------------------------------------------
*/
int
FrontMtx_nUpperBlocks (
   FrontMtx   *frontmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_nUpperBlocks(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_nUpperBlocks()"
           "\n data mode is 1, not 2\n") ;
   exit(-1) ;
}
return(frontmtx->upperblockIVL->tsize) ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- return a pointer to the upper block IVL object
 
   created -- 98jun13, cca
   ---------------------------------------------------------
*/
IVL *
FrontMtx_upperBlockIVL (
   FrontMtx   *frontmtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_upperBlockIVL(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_upperBlockIVL()"
           "\n data mode is 1, not 2\n") ;
   exit(-1) ;
}
return(frontmtx->upperblockIVL) ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- return a pointer to the lower block IVL object
 
   created -- 98jun13, cca
   ---------------------------------------------------------
*/
IVL *
FrontMtx_lowerBlockIVL (
   FrontMtx   *frontmtx
) {
IVL    *ivl ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_lowerBlockIVL(%p)"
           "\n bad input\n", frontmtx) ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_lowerBlockIVL()"
           "\n data mode is 1, not 2\n") ;
   exit(-1) ;
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   ivl = frontmtx->lowerblockIVL ;
} else {
   ivl = frontmtx->upperblockIVL ;
}
return(ivl) ; }
 
/*--------------------------------------------------------------------*/
