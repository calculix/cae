/*  map.c  */

#include "../../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- to map the coordinates of the entries of the matrix 
      into new coordinates. this method is used during the distributed
      matrix-vector multiply where a matrix local to a processor is
      mapped into a local coordinate system.

   (row,col) --> (rowmap[row],colmap[col])

   we check that row is in [0,nrow) and col is in [0,ncol),
   where nrow is the size of rowmapIV and ncol is the size of colmapIV.

   note, the storage mode is not changed. i.e., if the data is
   stored by vectors, it may be invalid after the indices have 
   been mapped. on the other hand, it may not, so it is the user's
   responsibility to reset the storage mode if necessary.

   created -- 98aug02, cca
   --------------------------------------------------------------------
*/
void
InpMtx_mapEntries (
   InpMtx   *inpmtx,
   IV       *rowmapIV,
   IV       *colmapIV
) {
int   chv, col, ii, ncol, nent, nrow, off, row ;
int   *colmap, *ivec1, *ivec2, *rowmap ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || rowmapIV == NULL || colmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( !(INPMTX_IS_BY_ROWS(inpmtx)
   ||  INPMTX_IS_BY_COLUMNS(inpmtx)
   ||  INPMTX_IS_BY_CHEVRONS(inpmtx)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
           "\n bad coordinate type\n") ;
   exit(-1) ;
}
IV_sizeAndEntries(rowmapIV, &nrow, &rowmap) ;
if ( rowmap == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
           "\n rowmap is NULL\n") ;
   exit(-1) ;
}
IV_sizeAndEntries(colmapIV, &ncol, &colmap) ;
if ( colmap == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
           "\n colmap is NULL\n") ;
   exit(-1) ;
}
nent  = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec1[ii] ; col = ivec2[ii] ;
      if ( row < 0 || row >= nrow ) {
         fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
                 "\n entry (%d,%d), nrow = %d\n", row, col, nrow) ;
         exit(-1) ;
      }
      ivec1[ii] = rowmap[row] ;
      if ( col < 0 || col >= ncol ) {
         fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
                 "\n entry (%d,%d), ncol = %d\n", row, col, ncol) ;
         exit(-1) ;
      }
      ivec2[ii] = colmap[col] ;
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec2[ii] ; col = ivec1[ii] ;
      if ( row < 0 || row >= nrow ) {
         fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
                 "\n entry (%d,%d), nrow = %d\n", row, col, nrow) ;
         exit(-1) ;
      }
      ivec2[ii] = rowmap[row] ;
      if ( col < 0 || col >= ncol ) {
         fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
                 "\n entry (%d,%d), ncol = %d\n", row, col, ncol) ;
         exit(-1) ;
      }
      ivec1[ii] = colmap[col] ;
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      chv = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chv ; col = chv + off ;
      } else {
         row = chv - off ; col = chv ;
      }
      if ( row < 0 || row >= nrow ) {
         fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
                 "\n entry (%d,%d), nrow = %d\n", row, col, nrow) ;
         exit(-1) ;
      }
      row = rowmap[row] ;
      if ( col < 0 || col >= ncol ) {
         fprintf(stderr, "\n fatal error in InpMtx_mapEntries()"
                 "\n entry (%d,%d), ncol = %d\n", row, col, ncol) ;
         exit(-1) ;
      }
      col = colmap[col] ;
      ivec1[ii] = (row >= col) ? col : row ;
      ivec2[ii] = col - row ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
