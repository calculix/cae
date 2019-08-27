/*  initRandom.c  */

#include "../../Drand.h"
#include "../SubMtx.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to initialize a matrix object with 
      random entries and possibly random structure.

   created -- 98feb16, cca
   ------------------------------------------------
*/
void
SubMtx_initRandom (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent,
   int      seed
) {
Drand   *drand ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || type < 1 || type > 2 || mode < 0 || mode > 9 
|| nrow <= 0 || ncol <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_initRandom()"
           "\n bad input\n") ;
   exit(-1) ;
}
SubMtx_clearData(mtx) ;
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
if ( seed > 0 ) {
   Drand_setSeed(drand, seed) ;
}
Drand_setUniform(drand, 0, 1) ;
/*
   ---------------------------
   switch over the matrix mode
   ---------------------------
*/
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2 ;
/*
   ----------------------------------------
   this case is simple, fill entries vector
   ----------------------------------------
*/
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nrow*ncol) ;
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nrow*ncol, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nrow*ncol, entries) ;
   }
   } break ;
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS :
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii ;
   int      *colids, *indices, *ivec, *ivec1, *ivec2, *rowids, *sizes ;
/*
   ----------------------------------------------------------
   (1) fill ivec[nrow*ncol] with 0, 1, ..., nrow*ncol-1
   (2) shuffle ivec[nrow*ncol] 
   (3) set rowids[nent] to be the row id 
       of the entries of ivec[0:nent-1]
   (4) set colids[nent] to be the column id 
       of the entries of ivec[0:nent-1]
   (5) for sparse rows and columns, set sizes[] and indices[]
   ----------------------------------------------------------
*/ 
   ivec = IVinit(nrow*ncol, -1) ;
   IVramp(nrow*ncol, ivec, 0, 1) ;
   IVshuffle(nrow*ncol, ivec, seed + 3) ;
   if ( nent > nrow*ncol ) {
      nent = nrow*ncol ;
   }
   rowids = IVinit(nent, -1) ;
   colids = IVinit(nent, -1) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      rowids[ii] = ivec[ii] % nrow ;
      colids[ii] = ivec[ii] / nrow ;
   }
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   switch ( mode ) {
   case SUBMTX_SPARSE_ROWS :
      IV2qsortUp(nent, rowids, colids) ;
      SubMtx_sparseRowsInfo(mtx, &nrow, &nent, 
                          &sizes, &indices, &entries) ;
      IVzero(nrow, sizes) ;
      for ( ii = 0 ; ii < nent ; ii++ ) {
         sizes[rowids[ii]]++ ;
      }
      IVcopy(nent, indices, colids) ;
      break ;
   case SUBMTX_SPARSE_COLUMNS :
      IV2qsortUp(nent, colids, rowids) ;
      SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                             &sizes, &indices, &entries) ;
      IVzero(ncol, sizes) ;
      for ( ii = 0 ; ii < nent ; ii++ ) {
         sizes[colids[ii]]++ ;
      }
      IVcopy(nent, indices, rowids) ;
      break ;
   case SUBMTX_SPARSE_TRIPLES :
      IV2qsortUp(nent, colids, rowids) ;
      SubMtx_sparseTriplesInfo(mtx, &nent, &ivec1, &ivec2, &entries) ;
      IVcopy(nent, ivec1, rowids) ;
      IVcopy(nent, ivec2, colids) ;
      break ;
   }
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec) ;
   IVfree(rowids) ;
   IVfree(colids) ;
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      irow, size, *firstlocs, *ivec1, *ivec2, *sizes ;
 
   ivec1 = IVinit(nrow, -1) ;
   ivec2 = IVinit(nrow,  0) ;
   Drand_setUniform(drand, 0, ncol) ;
   Drand_fillIvector(drand, nrow, ivec1) ;
   Drand_fillIvector(drand, nrow, ivec2) ;
   for ( irow = nent = 0 ; irow < nrow ; irow++ ) {
      if ( (size = ivec2[irow] - ivec1[irow] + 1) <= 0 ) {
         ivec1[irow] = -1 ;
         ivec2[irow] = 0 ;
      } else {
         nent += size ;
         ivec2[irow] = size ;
      }
   }
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent,
                         &firstlocs, &sizes, &entries) ;
   IVcopy(nrow, firstlocs, ivec1) ;
   IVcopy(nrow, sizes,     ivec2) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec1) ;
   IVfree(ivec2) ;
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      jcol, size, *firstlocs, *ivec1, *ivec2, *sizes ;
 
   ivec1 = IVinit(ncol, -1) ;
   ivec2 = IVinit(ncol, -1) ;
   Drand_setUniform(drand, 0, nrow) ;
   Drand_fillIvector(drand, ncol, ivec1) ;
   Drand_fillIvector(drand, ncol, ivec2) ;
   for ( jcol = nent = 0 ; jcol < ncol ; jcol++ ) {
      if ( (size = ivec2[jcol] - ivec1[jcol] + 1) <= 0 ) {
         ivec1[jcol] = -1 ;
         ivec2[jcol] =  0 ;
      } else {
         nent += size ;
         ivec2[jcol] = size ;
      }
   }
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent,
                            &firstlocs, &sizes, &entries) ;
   IVcopy(ncol, firstlocs, ivec1) ;
   IVcopy(ncol, sizes,     ivec2) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec1) ;
   IVfree(ivec2) ;
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
 
   ncol = nrow ;
   nent = nrow ;
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_diagonalInfo(mtx, &ncol, &entries) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   double   *entries ;
   int      ipivot, irow, m ;
   int      *pivots, *pivotsizes ;
 
   ncol = nrow ;
   pivotsizes = IVinit(nrow, -1) ;
   Drand_setUniform(drand, 1,3) ;
   Drand_fillIvector(drand, ncol, pivotsizes) ;
   for ( ipivot = irow = nent = 0 ; irow < nrow ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( irow + m > nrow ) {
         m = pivotsizes[ipivot] = nrow - irow ;
      }
      nent += (m*(m+1))/2 ;
      irow += m ;
   }
   IVzero(nrow - ipivot, pivotsizes + ipivot) ;
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivots, &entries) ;
   IVcopy(nrow, pivots, pivotsizes) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(pivotsizes) ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   *entries ;
   int      ipivot, irow, kk, m ;
   int      *pivots, *pivotsizes ;
 
   ncol = nrow ;
   pivotsizes = IVinit(nrow, -1) ;
   Drand_setUniform(drand, 1,3) ;
   Drand_fillIvector(drand, ncol, pivotsizes) ;
   for ( ipivot = irow = nent = 0 ; irow < nrow ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( irow + m > nrow ) {
         m = pivotsizes[ipivot] = nrow - irow ;
      }
      nent += (m*(m+1))/2 ;
      irow += m ;
   }
   IVzero(nrow - ipivot, pivotsizes + ipivot) ;
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivots, &entries) ;
   IVcopy(nrow, pivots, pivotsizes) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   if ( SUBMTX_IS_COMPLEX(mtx) ) {
      for ( ipivot = irow = kk = 0 ; irow < nrow ; ipivot++ ) {
         m = pivotsizes[ipivot] ;
         if ( m == 1 ) {
            entries[2*kk+1] = 0.0 ;
            kk++ ;
         } else {
            entries[2*kk+1] = 0.0 ;
            entries[2*kk+5] = 0.0 ;
            kk += 3 ;
         }
         irow += m ;
      }
   }
   IVfree(pivotsizes) ;
   } break ;
default :
   fprintf(stderr, "\n\n %% type %d not yet supported", type) ;
   exit(0) ;
   break ;
}
Drand_free(drand) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to initialize a matrix object with 
      random entries in the upper triangle 

   strict = 1 --> strict upper triangle

   created -- 98feb16, cca
   ------------------------------------------------
*/
void
SubMtx_initRandomUpperTriangle (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent,
   int      seed,
   int      strict
) {
Drand   *drand ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || type < 1 || type > 2 || mode < 0 || mode > 9 
|| nrow <= 0 || ncol <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_initRandomUpperTriangle()"
           "\n bad input\n") ;
   exit(-1) ;
}
SubMtx_clearData(mtx) ;
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
if ( seed > 0 ) {
   Drand_setSeed(drand, seed) ;
}
Drand_setUniform(drand, 0, 1) ;
/*
   ---------------------------
   switch over the matrix type
   ---------------------------
*/
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      ij, inc1, inc2, irow, jcol ;
/*
   -------------------------------------------------------------------
   (1) initialize the matrix and fill all entries with random numbers.
   (2) put zeros in (strict) lower triangle
   -------------------------------------------------------------------
*/
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nrow*ncol) ;
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
      if ( strict == 1 ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = 0 ; jcol <= irow ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[ij] = 0.0 ;
            }
         }
      } else {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = 0 ; jcol < irow ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[ij] = 0.0 ;
            }
         }
      }
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
      if ( strict == 1 ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = 0 ; jcol <= irow ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[2*ij] = entries[2*ij+1] = 0.0 ;
            }
         }
      } else {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = 0 ; jcol < irow ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[2*ij] = entries[2*ij+1] = 0.0 ;
            }
         }
      }
   }
   } break ;
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      count, ii, irow, jcol ;
   int      *colids, *indices, *ivec, *rowids, *sizes ;
/*
   ------------------------------------------
   fill ivec[*] with upper triangular entries
   ------------------------------------------
*/ 
   ivec = IVinit(nrow*ncol, -1) ;
   if ( strict == 1 ) {
      for ( irow = count = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = irow + 1 ; jcol < ncol ; jcol++ ) {
            ivec[count++] = irow + jcol*nrow ;
         }
      }
   } else {
      for ( irow = count = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = irow ; jcol < ncol ; jcol++ ) {
            ivec[count++] = irow + jcol*nrow ;
         }
      }
   }
   IVshuffle(count, ivec, seed) ;
   if ( nent > count ) {
      nent = count ;
   }
   rowids = IVinit(nent, -1) ;
   colids = IVinit(nent, -1) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      rowids[ii] = ivec[ii] % nrow ;
      colids[ii] = ivec[ii] / nrow ;
   }
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   if ( mode == SUBMTX_SPARSE_ROWS ) {
      SubMtx_sparseRowsInfo(mtx, &nrow, &nent, 
                          &sizes, &indices, &entries) ;
      IVzero(nrow, sizes) ;
      IV2qsortUp(nent, rowids, colids) ;
      for ( ii = 0 ; ii < nent ; ii++ ) {
         sizes[rowids[ii]]++ ;
      }
      IVcopy(nent, indices, colids) ;
   } else {
      SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                             &sizes, &indices, &entries) ;
      IV2qsortUp(nent, colids, rowids) ;
      IVzero(ncol, sizes) ;
      for ( ii = 0 ; ii < nent ; ii++ ) {
         sizes[colids[ii]]++ ;
      }
      IVcopy(nent, indices, rowids) ;
   }
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec) ;
   IVfree(rowids) ;
   IVfree(colids) ;
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      irow, size ;
   int      *firstlocs, *ivec1, *ivec2, *sizes ;
/* 
   ------------------------------------------------------------
   fill ivec1[] and ivec2[] with start and end column locations
   ------------------------------------------------------------
*/
   ivec1 = IVinit(nrow, -1) ;
   ivec2 = IVinit(nrow, -1) ;
   Drand_setUniform(drand, 0, ncol) ;
   Drand_fillIvector(drand, nrow, ivec1) ;
   Drand_fillIvector(drand, nrow, ivec2) ;
/*
   -------------------------------------------------
   modify ivec1[] to be on or just past the diagonal
   -------------------------------------------------
*/
   if ( strict == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         if ( ivec1[irow] <= irow ) {
            ivec1[irow] = irow + 1 ;
         }
      }
   } else {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         if ( ivec1[irow] < irow ) {
            ivec1[irow] = irow ;
         }
      }
   }
/*
   -----------------------------------------------------------
   set ivec2[] to be the sizes and count the number of entries
   -----------------------------------------------------------
*/
   for ( irow = nent = 0 ; irow < nrow ; irow++ ) {
      if ( (size = ivec2[irow] - ivec1[irow] + 1) <= 0 ) {
         ivec1[irow] = -1 ;
         ivec2[irow] = 0 ;
      } else {
         nent += size ;
         ivec2[irow] = size ;
      }
   }
/*
   -------------------------------------------------
   initialize, set first location, sizes and entries
   -------------------------------------------------
*/
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent,
                         &firstlocs, &sizes, &entries) ;
   IVcopy(nrow, firstlocs, ivec1) ;
   IVcopy(nrow, sizes,     ivec2) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec1) ;
   IVfree(ivec2) ;
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      jcol, size ;
   int      *firstlocs, *ivec1, *ivec2, *sizes ;
/* 
   ---------------------------------------------------------
   fill ivec1[] and ivec2[] with start and end row locations
   ---------------------------------------------------------
*/
   ivec1 = IVinit(ncol, -1) ;
   ivec2 = IVinit(ncol, -1) ;
   Drand_setUniform(drand, 0, nrow) ;
   Drand_fillIvector(drand, ncol, ivec1) ;
   Drand_fillIvector(drand, ncol, ivec2) ;
/*
   -------------------------------------------------
   modify ivec2[] to be on or just past the diagonal
   -------------------------------------------------
*/
   if ( strict == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         if ( ivec2[jcol] >= jcol ) {
            ivec2[jcol] = jcol - 1 ;
         }
      }
   } else {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         if ( ivec2[jcol] > jcol ) {
            ivec2[jcol] = jcol ;
         }
      }
   }
/*
   -----------------------------------------------------------
   set ivec2[] to be the sizes and count the number of entries
   -----------------------------------------------------------
*/
   for ( jcol = nent = 0 ; jcol < ncol ; jcol++ ) {
      if ( (size = ivec2[jcol] - ivec1[jcol] + 1) <= 0 ) {
         ivec1[jcol] = -1 ;
         ivec2[jcol] = 0 ;
      } else {
         nent += size ;
         ivec2[jcol] = size ;
      }
   }
/*
   -------------------------------------------------
   initialize, set first location, sizes and entries
   -------------------------------------------------
*/
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent,
                            &firstlocs, &sizes, &entries) ;
   IVcopy(nrow, firstlocs, ivec1) ;
   IVcopy(nrow, sizes,     ivec2) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec1) ;
   IVfree(ivec2) ;
   } break ;
default :
   fprintf(stderr, "\n\n %% type %d not yet supported", type) ;
   exit(0) ;
   break ;
}
/*
   --------------------------------
   free the random number generator
   --------------------------------
*/
Drand_free(drand) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to initialize a matrix object with 
      random entries in the lower triangle 

   strict = 1 --> strict lower triangle

   created -- 98feb16, cca
   ------------------------------------------------
*/
void
SubMtx_initRandomLowerTriangle (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent,
   int      seed,
   int      strict
) {
Drand   *drand ;
int     *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || type < 0 || type > 9 || nrow <= 0 || ncol <= 0 ) {
   fprintf(stderr, "\n fatal error in SubMtx_initRandomLowerTriangle()"
           "\n bad input\n") ;
   exit(-1) ;
}
SubMtx_clearData(mtx) ;
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
if ( seed > 0 ) {
   Drand_setSeed(drand, seed) ;
}
Drand_setUniform(drand, 0, 1) ;
/*
   ---------------------------
   switch over the matrix type
   ---------------------------
*/
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      ij, inc1, inc2, irow, jcol ;
/*
   -------------------------------------------------------------------
   (1) initialize the matrix and fill all entries with random numbers.
   (2) put zeros in (strict) upper triangle
   -------------------------------------------------------------------
*/
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nrow*ncol) ;
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
      if ( strict == 1 ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = irow ; jcol < ncol ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[ij] = 0.0 ;
            }
         }
      } else {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = irow + 1 ; jcol < ncol ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[ij] = 0.0 ;
            }
         }
      }
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
      if ( strict == 1 ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = irow ; jcol < ncol ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[2*ij] = entries[2*ij+1] = 0.0 ;
            }
         }
      } else {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            for ( jcol = irow + 1 ; jcol < ncol ; jcol++ ) {
               ij = irow*inc1 + jcol*inc2 ;
               entries[2*ij] = entries[2*ij+1] = 0.0 ;
            }
         }
      }
   }
   } break ;
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      count, ii, irow, jcol ;
   int      *colids, *indices, *ivec, *rowids, *sizes ;
 
   ivec = IVinit(nrow*ncol, -1) ;
   if ( strict == 1 ) {
      for ( irow = count = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < irow ; jcol++ ) {
            ivec[count++] = irow + jcol*nrow ;
         }
      }
   } else {
      for ( irow = count = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol <= irow ; jcol++ ) {
            ivec[count++] = irow + jcol*nrow ;
         }
      }
   }
   IVshuffle(count, ivec, seed) ;
   if ( nent > count ) {
      nent = count ;
   }
   rowids = IVinit(nent, -1) ;
   colids = IVinit(nent, -1) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      rowids[ii] = ivec[ii] % nrow ;
      colids[ii] = ivec[ii] / nrow ;
   }
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   if ( mode == SUBMTX_SPARSE_ROWS ) {
      SubMtx_sparseRowsInfo(mtx, &nrow, &nent, 
                          &sizes, &indices, &entries) ;
      IVzero(nrow, sizes) ;
      IV2qsortUp(nent, rowids, colids) ;
      for ( ii = 0 ; ii < nent ; ii++ ) {
         sizes[rowids[ii]]++ ;
      }
      IVcopy(nent, indices, colids) ;
   } else {
      SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                             &sizes, &indices, &entries) ;
      IV2qsortUp(nent, colids, rowids) ;
      IVzero(ncol, sizes) ;
      for ( ii = 0 ; ii < nent ; ii++ ) {
         sizes[colids[ii]]++ ;
      }
      IVcopy(nent, indices, rowids) ;
   }
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec) ;
   IVfree(rowids) ;
   IVfree(colids) ;
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      irow, size ;
   int      *firstlocs, *ivec1, *ivec2, *sizes ;
/* 
   ------------------------------------------------------------
   fill ivec1[] and ivec2[] with start and end column locations
   ------------------------------------------------------------
*/
   ivec1 = IVinit(nrow, -1) ;
   ivec2 = IVinit(nrow, -1) ;
   Drand_setUniform(drand, 0, ncol) ;
   Drand_fillIvector(drand, nrow, ivec1) ;
   Drand_fillIvector(drand, nrow, ivec2) ;
/*
   -------------------------------------------------
   modify ivec2[] to be on or just past the diagonal
   -------------------------------------------------
*/
   if ( strict == 1 ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         if ( ivec2[irow] >= irow ) {
            ivec2[irow] = irow - 1 ;
         }
      }
   } else {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         if ( ivec1[irow] > irow ) {
            ivec1[irow] = irow ;
         }
      }
   }
/*
   -----------------------------------------------------------
   set ivec2[] to be the sizes and count the number of entries
   -----------------------------------------------------------
*/
   for ( irow = nent = 0 ; irow < nrow ; irow++ ) {
      if ( (size = ivec2[irow] - ivec1[irow] + 1) <= 0 ) {
         ivec1[irow] = -1 ;
         ivec2[irow] = 0 ;
      } else {
         nent += size ;
         ivec2[irow] = size ;
      }
   }
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent,
                         &firstlocs, &sizes, &entries) ;
   IVcopy(nrow, firstlocs, ivec1) ;
   IVcopy(nrow, sizes,     ivec2) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec1) ;
   IVfree(ivec2) ;
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      jcol, size ;
   int      *firstlocs, *ivec1, *ivec2, *sizes ;
/* 
   ---------------------------------------------------------
   fill ivec1[] and ivec2[] with start and end row locations
   ---------------------------------------------------------
*/
   ivec1 = IVinit(ncol, -1) ;
   ivec2 = IVinit(ncol, -1) ;
   Drand_setUniform(drand, 0, nrow) ;
   Drand_fillIvector(drand, ncol, ivec1) ;
   Drand_fillIvector(drand, ncol, ivec2) ;
/*
   -------------------------------------------------
   modify ivec1[] to be on or just past the diagonal
   -------------------------------------------------
*/
   if ( strict == 1 ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         if ( ivec1[jcol] <= jcol ) {
            ivec1[jcol] = jcol + 1 ;
         }
      }
   } else {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         if ( ivec1[jcol] < jcol ) {
            ivec1[jcol] = jcol ;
         }
      }
   }
/*
   -----------------------------------------------------------
   set ivec2[] to be the sizes and count the number of entries
   -----------------------------------------------------------
*/
   for ( jcol = nent = 0 ; jcol < ncol ; jcol++ ) {
      if ( (size = ivec2[jcol] - ivec1[jcol] + 1) <= 0 ) {
         ivec1[jcol] = -1 ;
         ivec2[jcol] = 0 ;
      } else {
         nent += size ;
      }
   }
   SubMtx_init(mtx, type, mode, rowid, colid, nrow, ncol, nent) ;
   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent,
                            &firstlocs, &sizes, &entries) ;
   IVcopy(ncol, firstlocs, ivec1) ;
   IVcopy(ncol, sizes,     ivec2) ;
   if ( SUBMTX_IS_REAL(mtx) ) {
      Drand_fillDvector(drand, nent, entries) ;
   } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
      Drand_fillDvector(drand, 2*nent, entries) ;
   }
   IVfree(ivec1) ;
   IVfree(ivec2) ;
   } break ;
default :
   fprintf(stderr, "\n\n %% type %d not yet supported", type) ;
   exit(0) ;
   break ;
}
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
IVramp(nrow, rowind, 0, 1) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
IVramp(ncol, colind, 0, 1) ;
/*
   --------------------------------
   free the random number generator
   --------------------------------
*/
Drand_free(drand) ;

return ; }

/*--------------------------------------------------------------------*/
