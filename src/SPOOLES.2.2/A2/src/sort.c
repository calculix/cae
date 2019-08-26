/*  sort.c  */

#include "../A2.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   permute the rows of the matrix
   A(*,*) = A(index(*),*)
   this method calls A2_sortRowsUp
   but does not overwrite the index[] vector

   created -- 98apr15, cca
   -----------------------------------------
*/
void
A2_permuteRows (
   A2    *mtx,
   int   nrow,
   int   index[]
) {
int   *rowids ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || nrow < 0 || nrow > mtx->n1 || index == NULL ) {
   fprintf(stderr, "\n fatal error in A2_permuteRows(%p,%d,%p)"
           "\n bad input\n", mtx, nrow, index) ;
   exit(-1) ;
}
rowids = IVinit(nrow, -1) ;
IVcopy(nrow, rowids, index) ;
A2_sortRowsUp(mtx, nrow, rowids) ;
IVfree(rowids) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   permute the columns of the matrix
   A(*,*) = A(*,index(*))
   this method calls A2_sortColumnsUp
   but does not overwrite the index[] vector

   created -- 98apr15, cca
   -----------------------------------------
*/
void
A2_permuteColumns (
   A2   *mtx,
   int   ncol,
   int   index[]
) {
int   *colids ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || ncol < 0 || ncol > mtx->n2 || index == NULL ) {
   fprintf(stderr, "\n fatal error in A2_permuteColumns(%p,%d,%p)"
           "\n bad input\n", mtx, ncol, index) ;
   exit(-1) ;
}
colids = IVinit(ncol, -1) ;
IVcopy(ncol, colids, index) ;
A2_sortColumnsUp(mtx, ncol, colids) ;
IVfree(colids) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   sort the rows of the matrix in ascending order
   of the rowids[] vector. on return, rowids is
   in asending order. return value is the number
   of row swaps made.

   created -- 98apr15, cca
   ----------------------------------------------
*/
int
A2_sortRowsUp (
   A2    *mtx,
   int   nrow,
   int   rowids[]
) {
int   ii, minrow, minrowid, nswap, target ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || mtx->n1 < nrow || nrow < 0 || rowids == NULL ) {
   fprintf(stderr, "\n fatal error in A2_sortRowsUp(%p,%d,%p)"
           "\n bad input\n", mtx, nrow, rowids) ;
   if ( mtx != NULL ) {
      A2_writeStats(mtx, stderr) ;
   }
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_sortRowsUp(%p,%d,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, nrow, rowids, mtx->type) ;
   exit(-1) ;
}
nswap = 0 ;
if ( mtx->inc1 == 1 ) {
   double   *dvtmp ;
   int      jcol, ncol ;
   int      *ivtmp ;
/*
   ---------------------------------------------------
   matrix is stored by columns, so permute each column
   ---------------------------------------------------
*/
   ivtmp = IVinit(nrow, -1) ;
   if ( A2_IS_REAL(mtx) ) {
      dvtmp = DVinit(nrow, 0.0) ;
   } else if ( A2_IS_COMPLEX(mtx) ) {
      dvtmp = DVinit(2*nrow, 0.0) ;
   }
   IVramp(nrow, ivtmp, 0, 1) ;
   IV2qsortUp(nrow, rowids, ivtmp) ;
   ncol = mtx->n2 ;
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      if ( A2_IS_REAL(mtx) ) {
         DVcopy(nrow, dvtmp, A2_column(mtx, jcol)) ;
         DVgather(nrow, A2_column(mtx, jcol), dvtmp, ivtmp) ;
      } else if ( A2_IS_COMPLEX(mtx) ) {
         ZVcopy(nrow, dvtmp, A2_column(mtx, jcol)) ;
         ZVgather(nrow, A2_column(mtx, jcol), dvtmp, ivtmp) ;
      }
   }
   IVfree(ivtmp) ;
   DVfree(dvtmp) ;
} else {
/*
   ----------------------------------------
   use a simple insertion sort to swap rows
   ----------------------------------------
*/
   for ( target = 0 ; target < nrow ; target++ ) {
      minrow   = target ;
      minrowid = rowids[target] ;
      for ( ii = target + 1 ; ii < nrow ; ii++ ) {
         if ( minrowid > rowids[ii] ) {
            minrow   = ii ;
            minrowid = rowids[ii] ;
         }
      }
      if ( minrow != target ) {
         rowids[minrow] = rowids[target] ;
         rowids[target] = minrowid ;
         A2_swapRows(mtx, target, minrow) ;
         nswap++ ;
      }
   }
}

return(nswap) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   sort the columns of the matrix in ascending order
   of the colids[] vector. on return, colids is
   in asending order. return value is the number
   of column swaps made.

   created -- 98apr15, cca
   -------------------------------------------------
*/
int
A2_sortColumnsUp (
   A2   *mtx,
   int   ncol,
   int   colids[]
) {
int   ii, mincol, mincolid, nswap, target ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || mtx->n2 < ncol || ncol < 0 || colids == NULL ) {
   fprintf(stderr, "\n fatal error in A2_sortColumnsUp(%p,%d,%p)"
           "\n bad input\n", mtx, ncol, colids) ;
   if ( mtx != NULL ) {
      A2_writeStats(mtx, stderr) ;
   }
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtx) || A2_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in A2_sortColumnsUp(%p,%d,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, ncol, colids, mtx->type) ;
   exit(-1) ;
}
nswap = 0 ;
if ( mtx->inc2 == 1 ) {
   double   *dvtmp ;
   int      irow, nrow ;
   int      *ivtmp ;
/*
   ---------------------------------------------------
   matrix is stored by rows, so permute each row
   ---------------------------------------------------
*/
   ivtmp = IVinit(ncol, -1) ;
   if ( A2_IS_REAL(mtx) ) {
      dvtmp = DVinit(ncol, 0.0) ;
   } else if ( A2_IS_COMPLEX(mtx) ) {
      dvtmp = DVinit(2*ncol, 0.0) ;
   }
   IVramp(ncol, ivtmp, 0, 1) ;
   IV2qsortUp(ncol, colids, ivtmp) ;
   nrow = mtx->n1 ;
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      if ( A2_IS_REAL(mtx) ) {
         DVcopy(ncol, dvtmp, A2_row(mtx, irow)) ;
         DVgather(ncol, A2_row(mtx, irow), dvtmp, ivtmp) ;
      } else if ( A2_IS_COMPLEX(mtx) ) {
         ZVcopy(ncol, dvtmp, A2_row(mtx, irow)) ;
         ZVgather(ncol, A2_row(mtx, irow), dvtmp, ivtmp) ;
      }
   }
   IVfree(ivtmp) ;
   DVfree(dvtmp) ;
} else {
/*
   ----------------------------------------
   use a simple insertion sort to swap cols
   ----------------------------------------
*/
   for ( target = 0 ; target < ncol ; target++ ) {
      mincol   = target ;
      mincolid = colids[target] ;
      for ( ii = target + 1 ; ii < ncol ; ii++ ) {
         if ( mincolid > colids[ii] ) {
            mincol   = ii ;
            mincolid = colids[ii] ;
         }
      }
      if ( mincol != target ) {
         colids[mincol] = colids[target] ;
         colids[target] = mincolid ;
         A2_swapColumns(mtx, target, mincol) ;
         nswap++ ;
      }
   }
}
return(nswap) ; }

/*--------------------------------------------------------------------*/
