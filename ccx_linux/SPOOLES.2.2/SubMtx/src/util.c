/*  util.c  */

#include "../SubMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to fill rowDV with the entries 
              in row irow of the matrix

   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillRowDV (
   SubMtx   *mtx,
   int      irow,
   DV       *rowDV
) {
double   *rowvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || irow < 0 || rowDV == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillRowDV(%p,%d,%p)"
           "\n bad input\n", mtx, irow, rowDV) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_REAL(mtx) ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillRowDV(%p,%d,%p)"
           "\n type = %d, must be SPOOLES_REAL\n", 
           mtx, irow, rowDV, mtx->type) ;
   exit(-1) ;
}
DV_setSize(rowDV, mtx->ncol) ;
rowvec = DV_entries(rowDV) ;
DVzero(mtx->ncol, rowvec) ;
/*
   --------------------------------------
   switch over the different matrix types
   --------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, jcol, loc, ncol, nrow ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      loc = irow*inc1 + jcol*inc2 ;
      rowvec[jcol] = entries[loc] ;
   }
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jrow, kk, nent, nrow, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, 
                         &nrow, &nent, &sizes, &indices, &entries) ;
   for ( jrow = offset = 0 ; jrow < irow ; jrow++ ) {
       offset += sizes[jrow] ;
   }
   for ( ii = 0, kk = offset ; ii < sizes[irow] ; ii++, kk++ ) {
      rowvec[indices[kk]] = entries[kk] ;
   }
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jcol, kk, nent, ncol, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   for ( jcol = offset = 0 ; jcol < ncol ; jcol++ ) {
      for ( ii = 0, kk = offset ; ii < sizes[jcol] ; ii++, kk++ ) {
         if ( indices[kk] == irow ) {
            rowvec[jcol] = entries[kk] ;
            break ;
         }
      }
      offset += sizes[jcol] ;
   }
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent ;
   int      *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( rowids[ii] == irow ) {
         rowvec[colids[ii]] = entries[ii] ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      first, ii, jrow, kk, last, nent, nrow, offset ;
   int      *firstlocs, *sizes ;


   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n %% entries(%d) :", nent) ;
   DVfprintf(stdout, nent, entries) ;
#endif
   for ( jrow = offset = 0 ; jrow < irow ; jrow++ ) {
      offset += sizes[jrow] ;
   }
#if MYDEBUG > 0
fprintf(stdout, "\n %% irow = %d, offset = %d", irow, offset) ;
fprintf(stdout, "\n %% first = %d, size = %d", 
        firstlocs[irow], sizes[irow]) ;
#endif
   if ( sizes[irow] > 0 ) {
      first = firstlocs[irow] ;
      last  = first + sizes[irow] - 1  ;
      for ( kk = first, ii = offset ; kk <= last ; kk++, ii++ ) {
         rowvec[kk] = entries[ii] ;
#if MYDEBUG > 0
fprintf(stdout, "\n %% rowvec[%d] = entries[%d]", kk, ii) ;
#endif
      }
   }
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      first, jcol, last, loc, nent, ncol, offset ;
   int      *firstlocs, *sizes ;

   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   for ( jcol = offset = 0 ; jcol < ncol ; jcol++ ) {
      if ( sizes[jcol] > 0 ) {
         first = firstlocs[jcol] ;
         last  = first + sizes[jcol] - 1 ;
         if ( first <= irow && irow <= last ) {
            loc = offset + irow - first ;
            rowvec[jcol] = entries[loc] ;
         }
         offset += sizes[jcol] ;
      }
   }
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      nent ;

   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   rowvec[irow] = entries[irow] ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jrow <= irow && irow < jrow + m ) {
         stride = m - 1 ;
         kk += irow - jrow ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jrow ; ii <= irow ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[ii] = entries[kk] ;
            kk += stride-- ;
         }
         for (    ; ii < jrow + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[ii] = entries[kk] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jrow += m ;
   }
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jrow <= irow && irow < jrow + m ) {
         stride = m - 1 ;
         kk += irow - jrow ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jrow ; ii < irow ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[ii] = entries[kk] ;
            kk += stride-- ;
         }
         for (    ; ii < jrow + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[ii] = entries[kk] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jrow += m ;
   }
   } break ;
default :
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to fill rowZV with the entries 
              in row irow of the matrix

   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillRowZV (
   SubMtx   *mtx,
   int      irow,
   ZV       *rowZV
) {
double   *rowvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || irow < 0 || rowZV == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillRowZV(%p,%d,%p)"
           "\n bad input\n", mtx, irow, rowZV) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_COMPLEX(mtx) ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillRowZV(%p,%d,%p)"
           "\n type = %d, must be SPOOLES_COMPLEX\n", 
           mtx, irow, rowZV, mtx->type) ;
   exit(-1) ;
}
ZV_setSize(rowZV, mtx->ncol) ;
rowvec = ZV_entries(rowZV) ;
DVzero(2*mtx->ncol, rowvec) ;
/*
   --------------------------------------
   switch over the different matrix types
   --------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, jcol, loc, ncol, nrow ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      loc = irow*inc1 + jcol*inc2 ;
      rowvec[2*jcol]   = entries[2*loc]   ;
      rowvec[2*jcol+1] = entries[2*loc+1] ;
   }
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jrow, kk, nent, nrow, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries) ;
   for ( jrow = offset = 0 ; jrow < irow ; jrow++ ) {
       offset += sizes[jrow] ;
   }
   for ( ii = 0, kk = offset ; ii < sizes[irow] ; ii++, kk++ ) {
      rowvec[2*indices[kk]]   = entries[2*kk]   ;
      rowvec[2*indices[kk]+1] = entries[2*kk+1] ;
   }
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jcol, kk, nent, ncol, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   for ( jcol = offset = 0 ; jcol < ncol ; jcol++ ) {
      for ( ii = 0, kk = offset ; ii < sizes[jcol] ; ii++, kk++ ) {
         if ( indices[kk] == irow ) {
            rowvec[2*jcol]   = entries[2*kk]   ;
            rowvec[2*jcol+1] = entries[2*kk+1] ;
            break ;
         }
      }
      offset += sizes[jcol] ;
   }
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent ;
   int      *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( rowids[ii] == irow ) {
         rowvec[2*colids[ii]]   = entries[2*ii]   ;
         rowvec[2*colids[ii]+1] = entries[2*ii+1] ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      first, ii, jrow, kk, last, nent, nrow, offset ;
   int      *firstlocs, *sizes ;


   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                           &firstlocs, &sizes, &entries) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n %% entries(%d) :", nent) ;
   ZVfprintf(stdout, nent, entries) ;
#endif
   for ( jrow = offset = 0 ; jrow < irow ; jrow++ ) {
      offset += sizes[jrow] ;
   }
#if MYDEBUG > 0
fprintf(stdout, "\n %% irow = %d, offset = %d", irow, offset) ;
fprintf(stdout, "\n %% first = %d, size = %d", 
        firstlocs[irow], sizes[irow]) ;
#endif
   if ( sizes[irow] > 0 ) {
      first = firstlocs[irow] ;
      last  = first + sizes[irow] - 1  ;
      for ( kk = first, ii = offset ; kk <= last ; kk++, ii++ ) {
         rowvec[2*kk]   = entries[2*ii]   ;
         rowvec[2*kk+1] = entries[2*ii+1] ;
#if MYDEBUG > 0
fprintf(stdout, "\n %% rowvec[%d] = entries[%d]", kk, ii) ;
#endif
      }
   }
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      first, jcol, last, loc, nent, ncol, offset ;
   int      *firstlocs, *sizes ;

   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   for ( jcol = offset = 0 ; jcol < ncol ; jcol++ ) {
      if ( sizes[jcol] > 0 ) {
         first = firstlocs[jcol] ;
         last  = first + sizes[jcol] - 1 ;
         if ( first <= irow && irow <= last ) {
            loc = offset + irow - first ;
            rowvec[2*jcol]   = entries[2*loc]   ;
            rowvec[2*jcol+1] = entries[2*loc+1] ;
         }
         offset += sizes[jcol] ;
      }
   }
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      nent ;

   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   rowvec[2*irow]   = entries[2*irow]   ;
   rowvec[2*irow+1] = entries[2*irow+1] ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jrow <= irow && irow < jrow + m ) {
         stride = m - 1 ;
         kk += irow - jrow ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jrow ; ii <= irow ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[2*ii]   = entries[2*kk]   ;
            rowvec[2*ii+1] = entries[2*kk+1] ;
            kk += stride-- ;
         }
         for (    ; ii < jrow + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[2*ii]   = entries[2*kk]   ;
            rowvec[2*ii+1] = entries[2*kk+1] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jrow += m ;
   }
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jrow <= irow && irow < jrow + m ) {
         stride = m - 1 ;
         kk += irow - jrow ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jrow ; ii < irow ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[2*ii]   =  entries[2*kk]   ;
            rowvec[2*ii+1] = -entries[2*kk+1] ;
            kk += stride-- ;
         }
         for (    ; ii < jrow + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            rowvec[2*ii]   = entries[2*kk]   ;
            rowvec[2*ii+1] = entries[2*kk+1] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jrow += m ;
   }
   } break ;
default :
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to fill colDV with the entries 
              in column icol of the matrix

   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillColumnDV (
   SubMtx   *mtx,
   int      icol,
   DV       *colDV
) {
double   *colvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || icol < 0 || colDV == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillColumnDV(%p,%d,%p)"
           "\n bad input\n", mtx, icol, colDV) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_REAL(mtx) ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillColumnDV(%p,%d,%p)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, icol, colDV, mtx->type) ;
   exit(-1) ;
}
DV_setSize(colDV, mtx->nrow) ;
colvec = DV_entries(colDV) ;
DVzero(mtx->nrow, colvec) ;
/*
   --------------------------------------
   switch over the different matrix types
   --------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, jrow, loc, ncol, nrow ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   for ( jrow = 0 ; jrow < nrow ; jrow++ ) {
      loc = jrow*inc1 + icol*inc2 ;
      colvec[jrow] = entries[loc] ;
   }
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jcol, kk, ncol, nent, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   for ( jcol = offset = 0 ; jcol < icol ; jcol++ ) {
       offset += sizes[jcol] ;
   }
   for ( ii = 0, kk = offset ; ii < sizes[icol] ; ii++, kk++ ) {
      colvec[indices[kk]] = entries[kk] ;
   }
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jrow, kk, nent, nrow, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries) ;
   for ( jrow = offset = 0 ; jrow < nrow ; jrow++ ) {
      for ( ii = 0, kk = offset ; ii < sizes[jrow] ; ii++, kk++ ) {
         if ( indices[kk] == icol ) {
            colvec[jrow] = entries[kk] ;
            break ;
         }
      }
      offset += sizes[jrow] ;
   }
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent ;
   int      *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( colids[ii] == icol ) {
         colvec[rowids[ii]] = entries[ii] ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      first, ii, jcol, kk, last, ncol, nent, offset ;
   int      *firstlocs, *sizes ;


   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n %% entries(%d) :", nent) ;
   DVfprintf(stdout, nent, entries) ;
#endif
   for ( jcol = offset = 0 ; jcol < icol ; jcol++ ) {
      offset += sizes[jcol] ;
   }
#if MYDEBUG > 0
fprintf(stdout, "\n %% icol = %d, offset = %d", icol, offset) ;
fprintf(stdout, "\n %% first = %d, size = %d", 
        firstlocs[icol], sizes[icol]) ;
#endif
   if ( sizes[icol] > 0 ) {
      first = firstlocs[icol] ;
      last  = first + sizes[icol] - 1 ;
      for ( kk = first, ii = offset ; kk <= last ; kk++, ii++ ) {
         colvec[kk] = entries[ii] ;
#if MYDEBUG > 0
fprintf(stdout, "\n %% colvec[%d] = entries[%d]", kk, ii) ;
#endif
      }
   }
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      first, jrow, last, loc, nent, nrow, offset ;
   int      *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   for ( jrow = offset = 0 ; jrow < nrow ; jrow++ ) {
      if ( sizes[jrow] > 0 ) {
         first = firstlocs[jrow] ;
         last  = first + sizes[jrow] - 1 ;
         if ( first <= icol && icol <= last ) {
            loc = offset + icol - first ;
            colvec[jrow] = entries[loc] ;
         }
         offset += sizes[jrow] ;
      }
   }
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      nent ;

   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   colvec[icol] = entries[icol] ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   double   *entries ;
   int      ii, ipivot, jcol, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jcol = ipivot = kk = 0 ; jcol <= icol ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jcol <= icol && icol < jcol + m ) {
         stride = m - 1 ;
         kk += icol - jcol ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jcol ; ii <= icol ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[ii] = entries[kk] ;
            kk += stride-- ;
         }
         for (    ; ii < jcol + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[ii] = entries[kk] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jcol += m ;
   }
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   *entries ;
   int      ii, ipivot, jcol, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jcol = ipivot = kk = 0 ; jcol <= icol ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jcol <= icol && icol < jcol + m ) {
         stride = m - 1 ;
         kk += icol - jcol ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jcol ; ii <= icol ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[ii] = entries[kk] ;
            kk += stride-- ;
         }
         for (    ; ii < jcol + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[ii] = entries[kk] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jcol += m ;
   }
   } break ;
default :
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to fill colZV with the entries 
              in column icol of the matrix

   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillColumnZV (
   SubMtx   *mtx,
   int      icol,
   ZV       *colZV
) {
double   *colvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || icol < 0 || colZV == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillColumnZV(%p,%d,%p)"
           "\n bad input\n", mtx, icol, colZV) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_COMPLEX(mtx) ) {
   fprintf(stderr, "\n fatal error in SubMtx_fillColumnZV(%p,%d,%p)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, icol, colZV, mtx->type) ;
   exit(-1) ;
}
ZV_setSize(colZV, mtx->nrow) ;
colvec = ZV_entries(colZV) ;
DVzero(2*mtx->nrow, colvec) ;
/*
   --------------------------------------
   switch over the different matrix types
   --------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, jrow, loc, ncol, nrow ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   for ( jrow = 0 ; jrow < nrow ; jrow++ ) {
      loc = jrow*inc1 + icol*inc2 ;
      colvec[2*jrow]   = entries[2*loc]   ;
      colvec[2*jrow+1] = entries[2*loc+1] ;
   }
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jcol, kk, ncol, nent, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   for ( jcol = offset = 0 ; jcol < icol ; jcol++ ) {
       offset += sizes[jcol] ;
   }
   for ( ii = 0, kk = offset ; ii < sizes[icol] ; ii++, kk++ ) {
      colvec[2*indices[kk]]   = entries[2*kk] ;
      colvec[2*indices[kk]+1] = entries[2*kk+1] ;
   }
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jrow, kk, nent, nrow, offset ;
   int      *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries) ;
   for ( jrow = offset = 0 ; jrow < nrow ; jrow++ ) {
      for ( ii = 0, kk = offset ; ii < sizes[jrow] ; ii++, kk++ ) {
         if ( indices[kk] == icol ) {
            colvec[2*jrow]   = entries[2*kk]   ;
            colvec[2*jrow+1] = entries[2*kk+1] ;
            break ;
         }
      }
      offset += sizes[jrow] ;
   }
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent ;
   int      *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( colids[ii] == icol ) {
         colvec[2*rowids[ii]]   = entries[2*ii]   ;
         colvec[2*rowids[ii]+1] = entries[2*ii+1] ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      first, ii, jcol, kk, last, ncol, nent, offset ;
   int      *firstlocs, *sizes ;


   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n %% entries(%d) :", nent) ;
   ZVfprintf(stdout, nent, entries) ;
#endif
   for ( jcol = offset = 0 ; jcol < icol ; jcol++ ) {
      offset += sizes[jcol] ;
   }
#if MYDEBUG > 0
fprintf(stdout, "\n %% icol = %d, offset = %d", icol, offset) ;
fprintf(stdout, "\n %% first = %d, size = %d", 
        firstlocs[icol], sizes[icol]) ;
#endif
   if ( sizes[icol] > 0 ) {
      first = firstlocs[icol] ;
      last  = first + sizes[icol] - 1 ;
      for ( kk = first, ii = offset ; kk <= last ; kk++, ii++ ) {
         colvec[2*kk]   = entries[2*ii]   ;
         colvec[2*kk+1] = entries[2*ii+1] ;
#if MYDEBUG > 0
fprintf(stdout, "\n %% colvec[%d] = entries[%d]", kk, ii) ;
#endif
      }
   }
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      first, jrow, last, loc, nent, nrow, offset ;
   int      *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   for ( jrow = offset = 0 ; jrow < nrow ; jrow++ ) {
      if ( sizes[jrow] > 0 ) {
         first = firstlocs[jrow] ;
         last  = first + sizes[jrow] - 1 ;
         if ( first <= icol && icol <= last ) {
            loc = offset + icol - first ;
            colvec[2*jrow]   = entries[2*loc]   ;
            colvec[2*jrow+1] = entries[2*loc+1] ;
         }
         offset += sizes[jrow] ;
      }
   }
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      nent ;

   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   colvec[2*icol]   = entries[2*icol]   ;
   colvec[2*icol+1] = entries[2*icol+1] ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   double   *entries ;
   int      ii, ipivot, jcol, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jcol = ipivot = kk = 0 ; jcol <= icol ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jcol <= icol && icol < jcol + m ) {
         stride = m - 1 ;
         kk += icol - jcol ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jcol ; ii <= icol ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[2*ii]   = entries[2*kk]   ;
            colvec[2*ii+1] = entries[2*kk+1] ;
            kk += stride-- ;
         }
         for (    ; ii < jcol + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[2*ii]   = entries[2*kk]   ;
            colvec[2*ii+1] = entries[2*kk+1] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jcol += m ;
   }
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   *entries ;
   int      ii, ipivot, jcol, kk, m, nrow, nent, stride ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   for ( jcol = ipivot = kk = 0 ; jcol <= icol ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
/*
fprintf(stdout, "\n jrow %d, m %d, kk %d", jrow, m, kk) ;
fflush(stdout) ;
*/
      if ( jcol <= icol && icol < jcol + m ) {
         stride = m - 1 ;
         kk += icol - jcol ;
/*
fprintf(stdout, "\n    0. kk %d, stride %d", kk, stride) ;
fflush(stdout) ;
*/
         for ( ii = jcol ; ii <= icol ; ii++ ) {
/*
fprintf(stdout, "\n    1. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[2*ii]   = entries[2*kk]   ;
            colvec[2*ii+1] = entries[2*kk+1] ;
            kk += stride-- ;
         }
         for (    ; ii < jcol + m ; ii++ ) {
/*
fprintf(stdout, "\n    2. ii %d, kk %d", ii, kk) ;
fflush(stdout) ;
*/
            colvec[2*ii]   =  entries[2*kk]   ;
            colvec[2*ii+1] = -entries[2*kk+1] ;
            kk++ ;
         }
      } else {
         kk += (m*(m+1))/2 ;
      }
      jcol += m ;
   }
   } break ;
default :
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- return the magnitude of the largest element

   created -- 98may01, cca
   ------------------------------------------------------
*/
double
SubMtx_maxabs (
   SubMtx   *mtx
) {
double   maxabs ;
double   *entries ;
int      loc, nent ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_maxabs(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in SubMtx_maxabs(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
/*
   --------------------------------------
   switch over the different matrix types
   --------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   int      inc1, inc2, ncol, nrow ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   nent = nrow*ncol ;
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   int      ncol ;
   int      *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   } break ;
case SUBMTX_SPARSE_ROWS : {
   int      nrow ;
   int      *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries) ;
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   int      *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   int      ncol ;
   int      *firstlocs, *sizes ;


   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   int      nrow ;
   int      *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   } break ;
case SUBMTX_DIAGONAL : {
   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   int      nrow ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   int      nrow ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   } break ;
default :
   break ;
}
if ( SUBMTX_IS_REAL(mtx) ) {
   maxabs = DVmaxabs(nent, entries, &loc) ;
} else if ( SUBMTX_IS_COMPLEX(mtx) ) {
   maxabs = ZVmaxabs(nent, entries) ;
}
return(maxabs) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- zero the entries in the matrix

   created -- 98may04, cca
   -----------------------------------------
*/
void
SubMtx_zero (
   SubMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_zero(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( SUBMTX_IS_REAL(mtx) ) {
   DVzero(mtx->nent, mtx->entries) ;
} else if ( SUBMTX_IS_COMPLEX(mtx) ) {
   DVzero(2*mtx->nent, mtx->entries) ;
}
return ; }

/*--------------------------------------------------------------------*/
