/*  sort.c  */

#include "../SubMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- sort the rows of the matrix into ascending order

   created -- 98mar02, cca
   -----------------------------------------------------------
*/
void
SubMtx_sortRowsUp (
   SubMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_sortRowsUp(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n fatal error in SubMtx_sortRowsUp(%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtx, mtx->type) ;
   exit(-1) ;
}
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   A2       a2 ;
   double   *entries ;
   int      inc1, inc2, ncol, nrow ;
   int      *rowind ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   A2_setDefaultFields(&a2) ;
   A2_init(&a2, mtx->type, nrow, ncol, inc1, inc2, entries) ;
   SubMtx_rowIndices(mtx, &nrow, &rowind) ;
   A2_sortRowsUp(&a2, nrow, rowind) ;
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, irow, jrow, kk, nent, nrow, offset, rowid, size ;
   int      *indices, *ivtemp, *rowind, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries);
   SubMtx_rowIndices(mtx, &nrow, &rowind) ;
/*
   ----------------------------------------------------------------
   get a companion vector and fill with the row id's of the entries
   ----------------------------------------------------------------
*/
   ivtemp = IVinit(nent, -1) ;
   for ( irow = kk = 0 ; irow < nrow ; irow++ ) {
      rowid = rowind[irow] ;
      size  = sizes[irow] ;
#if MYDEBUG > 0
      fprintf(stdout, "\n rowid %d, size %d, kk %d", rowid, size, kk) ;
      fflush(stdout) ;
#endif
      for ( ii = 0 ; ii < size ; ii++, kk++ ) {
         ivtemp[kk] = rowid ;
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n ivtemp[%d]", nent) ;
   IVfprintf(stdout, nent, ivtemp) ;
   fflush(stdout) ;
#endif
/*
   -----------------------------------------------------------
   zero the sizes vector, sort the (rowid,colid,entry) triples 
   into ascending order of rowids, and sort the row indices
   -----------------------------------------------------------
*/
   IVzero(nrow, sizes) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      IV2DVqsortUp(nent, ivtemp, indices, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      IV2ZVqsortUp(nent, ivtemp, indices, entries) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n after sort, ivtemp[%d]", nent) ;
   IVfprintf(stdout, nent, ivtemp) ;
   fflush(stdout) ;
#endif
   IVqsortUp(nrow, rowind) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n after sort, rowind") ;
   IVfprintf(stdout, nrow, rowind) ;
   fflush(stdout) ;
#endif
/*
   ----------------------------------------------------------------
   sort each row in ascending order and fill the sizes[nrow] vector
   ----------------------------------------------------------------
*/
   rowid = ivtemp[0] ;
   jrow  = offset = 0 ;
   kk    = size = 1 ;
   while ( kk < nent ) {
      if ( ivtemp[kk] == rowid ) {
         size++ ;
      } else {
         while ( rowid != rowind[jrow] ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n rowid %d, rowind[%d] %d",
                    rowid, jrow, rowind[jrow]) ;
            fflush(stdout) ;
#endif
            jrow++ ;
         }
         sizes[jrow] = size ;
         if ( SUBMTX_IS_REAL(mtx) ) {
            IVDVqsortUp(size, indices + offset, entries + offset) ;
         } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
            IVZVqsortUp(size, indices + offset, entries + 2*offset) ;
         }
         rowid = ivtemp[kk] ;
         jrow++ ;
         offset += size ;
         size = 1 ;
      }
      kk++ ;
   }
   while ( rowid != rowind[jrow] ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n rowid %d, rowind[%d] %d",
              rowid, jrow, rowind[jrow]) ;
      fflush(stdout) ;
#endif
      jrow++ ;
   }
   sizes[jrow] = size ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      IVDVqsortUp(size, indices + offset, entries + offset) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      IVZVqsortUp(size, indices + offset, entries + 2*offset) ;
   }
   IVfree(ivtemp) ;
   } break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_sortRowsUp(%p)"
           "\n bad type = %d", mtx, mtx->type) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- sort the columns of the matrix into ascending order

   created -- 98mar02, cca
   --------------------------------------------------------------
*/
void
SubMtx_sortColumnsUp (
   SubMtx   *mtx
) {
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   A2       a2 ;
   double   *entries ;
   int      inc1, inc2, ncol, nrow ;
   int      *colind ;

   A2_setDefaultFields(&a2) ;
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   A2_init(&a2, mtx->type, nrow, ncol, inc1, inc2, entries) ;
   SubMtx_columnIndices(mtx, &ncol, &colind) ;
   A2_sortColumnsUp(&a2, ncol, colind) ;
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      colid, ii, jcol, kk, ncol, nent, offset, size ;
   int      *colind, *indices, *ivtemp, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, 
                          &ncol, &nent, &sizes, &indices, &entries) ;
   SubMtx_columnIndices(mtx, &ncol, &colind) ;
/*
   -------------------------------------------------------------------
   get a companion vector and fill with the column id's of the entries
   -------------------------------------------------------------------
*/
   ivtemp = IVinit(nent, -1) ;
   for ( jcol = kk = 0 ; jcol < ncol ; jcol++ ) {
      colid = colind[jcol] ;
      size  = sizes[jcol] ;
      for ( ii = 0 ; ii < size ; ii++, kk++ ) {
         ivtemp[kk] = colid ;
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n ivtemp[%d]", nent) ;
   IVfprintf(stdout, nent, ivtemp) ;
   fflush(stdout) ;
#endif
/*
   -----------------------------------------------------------
   zero the sizes vector, sort the (colid,rowid,entry) triples 
   into ascending order of colids, and sort the column indices
   -----------------------------------------------------------
*/
   IVzero(ncol, sizes) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      IV2DVqsortUp(nent, ivtemp, indices, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      IV2ZVqsortUp(nent, ivtemp, indices, entries) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n after sort, ivtemp[%d]", nent) ;
   IVfprintf(stdout, nent, ivtemp) ;
   fflush(stdout) ;
#endif
   IVqsortUp(ncol, colind) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n after sort, colind") ;
   IVfprintf(stdout, ncol, colind) ;
   fflush(stdout) ;
#endif
/*
   -------------------------------------------------------------------
   sort each column in ascending order and fill the sizes[ncol] vector
   -------------------------------------------------------------------
*/
   colid = ivtemp[0] ;
   jcol  = offset = 0 ;
   kk    = size = 1 ;
   while ( kk < nent ) {
      if ( ivtemp[kk] == colid ) {
         size++ ;
      } else {
         while ( colid != colind[jcol] ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n colid %d, colind[%d] %d",
                    colid, jcol, colind[jcol]) ;
            fflush(stdout) ;
#endif
            jcol++ ;
         }
         sizes[jcol] = size ;
         if ( SUBMTX_IS_REAL(mtx) ) {
            IVDVqsortUp(size, indices + offset, entries + offset) ;
         } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
            IVZVqsortUp(size, indices + offset, entries + 2*offset) ;
         }
         colid = ivtemp[kk] ;
         jcol++ ;
         offset += size ;
         size = 1 ;
      }
      kk++ ;
   }
   while ( colid != colind[jcol] ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n colid %d, colind[%d] %d",
              colid, jcol, colind[jcol]) ;
      fflush(stdout) ;
#endif
      jcol++ ;
   }
   sizes[jcol] = size ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      IVDVqsortUp(size, indices + offset, entries + offset) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      IVZVqsortUp(size, indices + offset, entries + 2*offset) ;
   }
   IVfree(ivtemp) ;
   } break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_sortColumnsUp(%p)"
           "\n bad type = %d", mtx, mtx->type) ;
   SubMtx_writeForHumanEye(mtx, stderr) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
