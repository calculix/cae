/*  permute.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   permute the entries

   created -- 96jul05, cca
   -----------------------
*/
void
InpMtx_permute (
   InpMtx   *inpmtx,
   int       rowOldToNew[],
   int       colOldToNew[]
) {
int      col, ii, nent, row ;
int      *ivec1, *ivec2 ;
/*
   --------------
   check the data
   --------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_permute(%p,%p,%p)"
           "\n bad input\n", inpmtx, rowOldToNew, colOldToNew);
   exit(-1) ;
}
if ( inpmtx->coordType <= 0 || inpmtx->coordType >= 4 ) {
   fprintf(stderr, "\n fatal error in InpMtx_permute(%p,%p,%p)"
           "\n coordType = %d, must be 1, 2 or 3\n", 
           inpmtx, rowOldToNew, colOldToNew, inpmtx->coordType);
   exit(-1) ;
}
/*
   ----------------------
   check for quick return
   ----------------------
*/
if ( rowOldToNew == NULL && colOldToNew == NULL ) {
   return ; 
}
if ( (nent = inpmtx->nent) == 0 ) {
   return ;
}
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
if ( ivec1 == NULL || ivec2 == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_permute(%p,%p,%p)"
           "\n nent = %d, ivec1 = %p, ivec2 = %p",
           inpmtx, rowOldToNew, colOldToNew, nent, ivec1, ivec2) ;
   exit(-1) ;
}
/*
   --------------------------------------
   convert coordinates to new permutation
   --------------------------------------
*/
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec1[ii] ; col = ivec2[ii] ;
      if ( 0 <= row && rowOldToNew != NULL ) {
         row = rowOldToNew[row] ;
      }
      if ( 0 <= col && colOldToNew != NULL ) {
         col = colOldToNew[col] ;
      }
      ivec1[ii] = row ; ivec2[ii] = col ;
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      col = ivec1[ii] ; row = ivec2[ii] ;
      if ( 0 <= row && rowOldToNew != NULL ) {
         row = rowOldToNew[row] ;
      }
      if ( 0 <= col && colOldToNew != NULL ) {
         col = colOldToNew[col] ;
      }
      ivec1[ii] = col ; ivec2[ii] = row ;
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   int   chv, off ;

   for ( ii = 0 ; ii < nent ; ii++ ) {
      chv = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chv ; col = chv + off ;
      } else {
         col = chv ; row = chv - off ;
      }
      if ( 0 <= row && rowOldToNew != NULL ) {
         row = rowOldToNew[row] ;
      }
      if ( 0 <= col && colOldToNew != NULL ) {
         col = colOldToNew[col] ;
      }
      ivec1[ii] = (col <= row) ? col : row ;
      ivec2[ii] = col - row ;
   }
} 
/*
   -----------------------------------
   set the storage mode to raw triples
   -----------------------------------
*/
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }
   
/*--------------------------------------------------------------------*/
