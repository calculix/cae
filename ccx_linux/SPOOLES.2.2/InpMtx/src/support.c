/*  support.c  */

#include "../../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- 
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      nonsymmetric matrix.

      rowsupIV -- filled with row indices of y[] which will be updated.
      colsupIV -- filled with row indices of x[] which will be used.

   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportNonsym (
   InpMtx   *A,
   IV       *rowsupIV,
   IV       *colsupIV
) {
char   *colmark, *rowmark ;
int    chev, col, colcount, ii, loc, maxcol, maxrow, nent, off, row,
       rowcount ;
int    *colsup, *ivec1, *ivec2, *rowsup ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || rowsupIV == NULL || colsupIV == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportNonsym(%p,%p,%p)"
           "\n bad input\n", A, rowsupIV, colsupIV) ;
   exit(-1) ;
}
if (  !INPMTX_IS_BY_ROWS(A) 
   && !INPMTX_IS_BY_COLUMNS(A) 
   && !INPMTX_IS_BY_CHEVRONS(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportNonsym(%p,%p,%p)"
           "\n coordinate type\n", A, rowsupIV, colsupIV) ;
   exit(-1) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
nent  = A->nent ;
/*
   -----------------------------------------------------------------
   (1) determine the maximum row and column numbers in these entries
   (2) allocate marking vectors for rows and columns
   (3) fill marking vectors for rows and columns
   (4) fill support vectors 
   -----------------------------------------------------------------
*/
if ( INPMTX_IS_BY_ROWS(A) ) {
   maxrow   = IVmax(nent, ivec1, &loc) ;
   maxcol   = IVmax(nent, ivec2, &loc) ;
   rowmark  = CVinit(1+maxrow, 'O') ;
   colmark  = CVinit(1+maxcol, 'O') ;
   rowcount = colcount = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec1[ii] ; col = ivec2[ii] ;
      if ( rowmark[row] == 'O' ) {
         rowcount++ ;
      }
      rowmark[row] = 'X' ;
      if ( colmark[col] == 'O' ) {
         colcount++ ;
      }
      colmark[col] = 'X' ;
   }
} else if ( INPMTX_IS_BY_COLUMNS(A) ) {
   maxrow   = IVmax(nent, ivec2, &loc) ;
   maxcol   = IVmax(nent, ivec1, &loc) ;
   rowmark  = CVinit(1+maxrow, 'O') ;
   colmark  = CVinit(1+maxcol, 'O') ;
   rowcount = colcount = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec2[ii] ; col = ivec1[ii] ;
      if ( rowmark[row] == 'O' ) {
         rowcount++ ;
      }
      rowmark[row] = 'X' ;
      if ( colmark[col] == 'O' ) {
         colcount++ ;
      }
      colmark[col] = 'X' ;
   }
} else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
   maxrow = maxcol = -1 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      chev = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chev ; col = chev + off ;
      } else {
         col = chev ; row = chev - off ;
      }
      if ( maxrow < row ) {
         maxrow = row ;
      }
      if ( maxcol < col ) {
         maxcol = col ;
      }
   }
   rowmark  = CVinit(1+maxrow, 'O') ;
   colmark  = CVinit(1+maxcol, 'O') ;
   rowcount = colcount = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      chev = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chev ; col = chev + off ;
      } else {
         col = chev ; row = chev - off ;
      }
      if ( rowmark[row] == 'O' ) {
         rowcount++ ;
      }
      rowmark[row] = 'X' ;
      if ( colmark[col] == 'O' ) {
         colcount++ ;
      }
      colmark[col] = 'X' ;
   }
}
IV_setSize(rowsupIV, rowcount) ;
rowsup = IV_entries(rowsupIV) ;
for ( row = rowcount = 0 ; row <= maxrow ; row++ ) {
   if ( rowmark[row] == 'X' ) {
      rowsup[rowcount++] = row ;
   }
}
IV_setSize(colsupIV, colcount) ;
colsup = IV_entries(colsupIV) ;
for ( col = colcount = 0 ; col <= maxcol ; col++ ) {
   if ( colmark[col] == 'X' ) {
      colsup[colcount++] = col ;
   }
}
CVfree(colmark) ;
CVfree(rowmark) ;

return ; }
   
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- 
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      nonsymmetric matrix.

      rowsupIV -- filled with row indices of y[] which will be updated.
      colsupIV -- filled with row indices of x[] which will be used.

   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportNonsymT (
   InpMtx   *A,
   IV       *rowsupIV,
   IV       *colsupIV
) {
char   *colmark, *rowmark ;
int    chev, col, colcount, ii, loc, maxcol, maxrow, nent, off, row,
       rowcount ;
int    *colsup, *ivec1, *ivec2, *rowsup ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || rowsupIV == NULL || colsupIV == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportNonsymT(%p,%p,%p)"
           "\n bad input\n", A, rowsupIV, colsupIV) ;
   exit(-1) ;
}
if (  !INPMTX_IS_BY_ROWS(A) 
   && !INPMTX_IS_BY_COLUMNS(A) 
   && !INPMTX_IS_BY_CHEVRONS(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportNonsymT(%p,%p,%p)"
           "\n coordinate type\n", A, rowsupIV, colsupIV) ;
   exit(-1) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
nent  = A->nent ;
/*
   -----------------------------------------------------------------
   (1) determine the maximum row and column numbers in these entries
   (2) allocate marking vectors for rows and columns
   (3) fill marking vectors for rows and columns
   (4) fill support vectors 
   -----------------------------------------------------------------
*/
if ( INPMTX_IS_BY_ROWS(A) ) {
   maxrow   = IVmax(nent, ivec1, &loc) ;
   maxcol   = IVmax(nent, ivec2, &loc) ;
   rowmark  = CVinit(1+maxcol, 'O') ;
   colmark  = CVinit(1+maxrow, 'O') ;
   rowcount = colcount = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec1[ii] ; col = ivec2[ii] ;
      if ( colmark[row] == 'O' ) {
         colcount++ ;
      }
      colmark[row] = 'X' ;
      if ( rowmark[col] == 'O' ) {
         rowcount++ ;
      }
      rowmark[col] = 'X' ;
   }
} else if ( INPMTX_IS_BY_COLUMNS(A) ) {
   maxrow   = IVmax(nent, ivec2, &loc) ;
   maxcol   = IVmax(nent, ivec1, &loc) ;
   rowmark  = CVinit(1+maxcol, 'O') ;
   colmark  = CVinit(1+maxrow, 'O') ;
   rowcount = colcount = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec2[ii] ; col = ivec1[ii] ;
      if ( colmark[row] == 'O' ) {
         colcount++ ;
      }
      colmark[row] = 'X' ;
      if ( rowmark[col] == 'O' ) {
         rowcount++ ;
      }
      rowmark[col] = 'X' ;
   }
} else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
   maxrow = maxcol = -1 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      chev = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chev ; col = chev + off ;
      } else {
         col = chev ; row = chev - off ;
      }
      if ( maxrow < row ) {
         maxrow = row ;
      }
      if ( maxcol < col ) {
         maxcol = col ;
      }
   }
   rowmark  = CVinit(1+maxcol, 'O') ;
   colmark  = CVinit(1+maxrow, 'O') ;
   rowcount = colcount = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      chev = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chev ; col = chev + off ;
      } else {
         col = chev ; row = chev - off ;
      }
      if ( colmark[row] == 'O' ) {
         colcount++ ;
      }
      colmark[row] = 'X' ;
      if ( rowmark[col] == 'O' ) {
         rowcount++ ;
      }
      rowmark[col] = 'X' ;
   }
}
IV_setSize(rowsupIV, rowcount) ;
rowsup = IV_entries(rowsupIV) ;
for ( col = rowcount = 0 ; col <= maxcol ; col++ ) {
   if ( rowmark[col] == 'X' ) {
      rowsup[rowcount++] = col ;
   }
}
IV_setSize(colsupIV, colcount) ;
colsup = IV_entries(colsupIV) ;
for ( row = colcount = 0 ; row <= maxrow ; row++ ) {
   if ( colmark[row] == 'X' ) {
      colsup[colcount++] = row ;
   }
}
CVfree(colmark) ;
CVfree(rowmark) ;

return ; }
   
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- 
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      nonsymmetric matrix.

      rowsupIV -- filled with row indices of y[] which will be updated.
      colsupIV -- filled with row indices of x[] which will be used.

   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportNonsymH (
   InpMtx   *A,
   IV       *rowsupIV,
   IV       *colsupIV
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || rowsupIV == NULL || colsupIV == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportNonsymH(%p,%p,%p)"
           "\n bad input\n", A, rowsupIV, colsupIV) ;
   exit(-1) ;
}
if (  !INPMTX_IS_BY_ROWS(A) 
   && !INPMTX_IS_BY_COLUMNS(A) 
   && !INPMTX_IS_BY_CHEVRONS(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportNonsymH(%p,%p,%p)"
           "\n coordinate type\n", A, rowsupIV, colsupIV) ;
   exit(-1) ;
}
InpMtx_supportNonsymT(A, rowsupIV, colsupIV) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- 
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a 
      symmetric matrix.

      supIV -- filled with row indices of y[] which will be updated
               and row indices of x[] which will be used.

   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportSym (
   InpMtx   *A,
   IV       *supIV
) {
char   *mark ;
int    chev, col, count, ii, loc, maxcol, maxrow, maxv, nent, off, row ;
int    *ivec1, *ivec2, *sup ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || supIV == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportSym(%p,%p)"
           "\n bad input\n", A, supIV) ;
   exit(-1) ;
}
if (  !INPMTX_IS_BY_ROWS(A) 
   && !INPMTX_IS_BY_COLUMNS(A) 
   && !INPMTX_IS_BY_CHEVRONS(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportSym(%p,%p)"
           "\n coordinate type\n", A, supIV) ;
   exit(-1) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
nent  = A->nent ;
/*
   -----------------------------------------------------------------
   (1) determine the maximum row and column numbers in these entries
   (2) allocate marking vectors for rows and columns
   (3) fill marking vectors for rows and columns
   (4) fill support vectors 
   -----------------------------------------------------------------
*/
if ( INPMTX_IS_BY_ROWS(A) ) {
   maxrow = IVmax(nent, ivec1, &loc) ;
   maxcol = IVmax(nent, ivec2, &loc) ;
   maxv   = (maxrow >= maxcol) ? maxrow : maxcol ;
   mark   = CVinit(1+maxv, 'O') ;
   count  = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec1[ii] ; col = ivec2[ii] ;
      if ( mark[row] == 'O' ) {
         count++ ;
      }
      mark[row] = 'X' ;
      if ( mark[col] == 'O' ) {
         count++ ;
      }
      mark[col] = 'X' ;
   }
} else if ( INPMTX_IS_BY_COLUMNS(A) ) {
   maxrow = IVmax(nent, ivec2, &loc) ;
   maxcol = IVmax(nent, ivec1, &loc) ;
   maxv   = (maxrow >= maxcol) ? maxrow : maxcol ;
   mark   = CVinit(1+maxv, 'O') ;
   count  = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      row = ivec2[ii] ; col = ivec1[ii] ;
      if ( mark[row] == 'O' ) {
         count++ ;
      }
      mark[row] = 'X' ;
      if ( mark[col] == 'O' ) {
         count++ ;
      }
      mark[col] = 'X' ;
   }
} else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
   maxv = -1 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      chev = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chev ; col = chev + off ;
         if ( maxv < col ) {
            maxv = col ;
         }
      } else {
         col = chev ; row = chev - off ;
         if ( maxv < row ) {
            maxv = row ;
         }
      }
   }
   mark = CVinit(1+maxv, 'O') ;
   count = 0 ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      chev = ivec1[ii] ; off = ivec2[ii] ;
      if ( off >= 0 ) {
         row = chev ; col = chev + off ;
      } else {
         col = chev ; row = chev - off ;
      }
      if ( mark[row] == 'O' ) {
         count++ ;
      }
      mark[row] = 'X' ;
      if ( mark[col] == 'O' ) {
         count++ ;
      }
      mark[col] = 'X' ;
   }
}
IV_setSize(supIV, count) ;
sup = IV_entries(supIV) ;
for ( row = count = 0 ; row <= maxv ; row++ ) {
   if ( mark[row] == 'X' ) {
      sup[count++] = row ;
   }
}
CVfree(mark) ;

return ; }
   
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- 
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a 
      Hermitian matrix.

      supIV -- filled with row indices of y[] which will be updated
               and row indices of x[] which will be used.

   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportHerm (
   InpMtx   *A,
   IV       *supIV
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || supIV == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportHerm(%p,%p)"
           "\n bad input\n", A, supIV) ;
   exit(-1) ;
}
if (  !INPMTX_IS_BY_ROWS(A) 
   && !INPMTX_IS_BY_COLUMNS(A) 
   && !INPMTX_IS_BY_CHEVRONS(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_supportHerm(%p,%p)"
           "\n coordinate type\n", A, supIV) ;
   exit(-1) ;
}
InpMtx_supportSym(A, supIV) ;

return ; }

/*--------------------------------------------------------------------*/
