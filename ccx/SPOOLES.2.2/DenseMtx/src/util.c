/*  util.c  */

#include "../DenseMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   sort the rows so the row ids are in ascending order
   sort the columns so the column ids are in ascending order

   created -- 98may02, cca
   ---------------------------------------------------------
*/
void
DenseMtx_sort (
   DenseMtx   *mtx
) {
A2    a2 ;
int   ii, ncol, nrow, sortColumns, sortRows ;
int   *colind, *rowind ;
/*
   ----------------
   check the output
   ----------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_sort(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
if ( nrow <= 0 || ncol <= 0 ) {
   return ;
}
sortRows = sortColumns = 0 ;
for ( ii = 1 ; ii < nrow ; ii++ ) {
   if ( rowind[ii-1] > rowind[ii] ) {
      sortRows = 1 ;
      break ;
   }
}
for ( ii = 1 ; ii < ncol ; ii++ ) {
   if ( colind[ii-1] > colind[ii] ) {
      sortColumns = 1 ;
      break ;
   }
}
if ( sortRows == 0 && sortColumns == 0 ) {
   return ;
}
A2_setDefaultFields(&a2) ;
DenseMtx_setA2(mtx, &a2) ;
if ( sortRows == 1 ) {
   A2_sortRowsUp(&a2, nrow, rowind) ;
}
if ( sortColumns == 1 ) {
   A2_sortColumnsUp(&a2, ncol, colind) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   copy row irowA from mtxA into row irowB in mtxB

   created -- 98may02, cca
   -----------------------------------------------
*/
void
DenseMtx_copyRow (
   DenseMtx   *mtxB,
   int        irowB,
   DenseMtx   *mtxA,
   int        irowA
) {
double   *rowA, *rowB ;
int      ii, inc2A, inc2B, iA, iB, ncol ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtxB == NULL || irowB < 0 || irowB >= mtxB->nrow 
   || mtxA == NULL || irowA < 0 || irowA >= mtxA->nrow 
   || (ncol = mtxA->ncol) != mtxB->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_copyRow(%p,%d,%p,%d)"
           "\n bad input\n", mtxB, irowB, mtxA, irowA) ;
   exit(-1) ;
}
inc2A = mtxA->inc2 ;
inc2B = mtxB->inc2 ;
/*
mtxB->rowind[irowB] = mtxA->rowind[irowA] ; 
*/
if ( DENSEMTX_IS_REAL(mtxB) && DENSEMTX_IS_REAL(mtxA) ) {
   rowA  = mtxA->entries + irowA*mtxA->inc1 ;
   rowB  = mtxB->entries + irowB*mtxB->inc1 ;
   for ( ii = iA = iB = 0 ; ii < ncol ; ii++, iA += inc2A, iB += inc2B){
      rowB[iB] = rowA[iA] ;
   }
} else if ( DENSEMTX_IS_COMPLEX(mtxB) && DENSEMTX_IS_COMPLEX(mtxA) ) {
   rowA  = mtxA->entries + 2*irowA*mtxA->inc1 ;
   rowB  = mtxB->entries + 2*irowB*mtxB->inc1 ;
   for ( ii = iA = iB = 0 ; ii < ncol ; ii++, iA += inc2A, iB += inc2B){
      rowB[2*iB]   = rowA[2*iA] ;
      rowB[2*iB+1] = rowA[2*iA+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   copy row irowA from mtxA into row irowB in mtxB
   and copy row index irowB from mtxB into row index irowA of mtxA

   created -- 98aug12, cca
   ---------------------------------------------------------------
*/
void
DenseMtx_copyRowAndIndex (
   DenseMtx   *mtxB,
   int        irowB,
   DenseMtx   *mtxA,
   int        irowA
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtxB == NULL || irowB < 0 || irowB >= mtxB->nrow 
   || mtxA == NULL || irowA < 0 || irowA >= mtxA->nrow 
   || mtxA->ncol != mtxB->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_copyRow(%p,%d,%p,%d)"
           "\n bad input\n", mtxB, irowB, mtxA, irowA) ;
   exit(-1) ;
}
DenseMtx_copyRow(mtxB, irowB, mtxA, irowA) ;
mtxB->rowind[irowB] = mtxA->rowind[irowA] ; 

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   add row irowA from mtxA into row irowB in mtxB

   created -- 98aug12, cca
   ----------------------------------------------
*/
void
DenseMtx_addRow (
   DenseMtx   *mtxB,
   int        irowB,
   DenseMtx   *mtxA,
   int        irowA
) {
double   *rowA, *rowB ;
int      ii, inc2A, inc2B, iA, iB, ncol ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtxB == NULL || irowB < 0 || irowB >= mtxB->nrow 
   || mtxA == NULL || irowA < 0 || irowA >= mtxA->nrow 
   || (ncol = mtxA->ncol) != mtxB->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_addRow(%p,%d,%p,%d)"
           "\n bad input\n", mtxB, irowB, mtxA, irowA) ;
   exit(-1) ;
}
inc2A = mtxA->inc2 ;
inc2B = mtxB->inc2 ;
if ( DENSEMTX_IS_REAL(mtxB) && DENSEMTX_IS_REAL(mtxA) ) {
   rowA  = mtxA->entries + irowA*mtxA->inc1 ;
   rowB  = mtxB->entries + irowB*mtxB->inc1 ;
   for ( ii = iA = iB = 0 ; ii < ncol ; ii++, iA += inc2A, iB += inc2B){
      rowB[iB] += rowA[iA] ;
   }
} else if ( DENSEMTX_IS_COMPLEX(mtxB) && DENSEMTX_IS_COMPLEX(mtxA) ) {
   rowA  = mtxA->entries + 2*irowA*mtxA->inc1 ;
   rowB  = mtxB->entries + 2*irowB*mtxB->inc1 ;
   for ( ii = iA = iB = 0 ; ii < ncol ; ii++, iA += inc2A, iB += inc2B){
      rowB[2*iB]   += rowA[2*iA] ;
      rowB[2*iB+1] += rowA[2*iA+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   zero the entries

   created -- 98may16, cca
   -----------------------
*/
void
DenseMtx_zero (
   DenseMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_zero(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( DENSEMTX_IS_REAL(mtx) ) {
   DVzero(mtx->nrow*mtx->ncol, mtx->entries) ;
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   DVzero(2*mtx->nrow*mtx->ncol, mtx->entries) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------
   fill with random entries

   created -- 98may16, cca
   ------------------------
*/
void
DenseMtx_fillRandomEntries (
   DenseMtx   *mtx,
   Drand      *drand
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || drand == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_fillRandomEntries(%p,%p)"
           "\n bad input\n", mtx, drand) ;
   exit(-1) ;
}
if ( DENSEMTX_IS_REAL(mtx) ) {
   Drand_fillDvector(drand, mtx->nrow*mtx->ncol, mtx->entries) ;
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   Drand_fillDvector(drand, 2*mtx->nrow*mtx->ncol, mtx->entries) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   compute three checksums
     sums[0] = sum of row indices
     sums[1] = sum of columns indices
     sums[2] = sum of entry magnitudes

   created -- 98may16, cca
   -----------------------------------
*/
void
DenseMtx_checksums (
   DenseMtx   *mtx,
   double     sums[]
) {
double   *entries ;
int      ii, ncol, nent, nrow ;
int      *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || sums == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_checksums(%p,%p)"
           "\n bad input\n", mtx, sums) ;
   exit(-1) ;
}
sums[0] = sums[1] = sums[2] = 0.0 ;
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
for ( ii = 0 ; ii < nrow ; ii++ ) {
   sums[0] += rowind[ii] ;
}
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
for ( ii = 0 ; ii < ncol ; ii++ ) {
   sums[1] += colind[ii] ;
}
entries = DenseMtx_entries(mtx) ;
nent    = nrow*ncol ;
if ( DENSEMTX_IS_REAL(mtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      sums[2] += fabs(entries[ii]) ;
   }
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   for ( ii = 0 ; ii < nent ; ii++ ) {
      sums[2] += Zabs(entries[2*ii], entries[2*ii+1]) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   return the maximum magnitude of the entries

   created -- 98may15, cca
   -------------------------------------------
*/
double
DenseMtx_maxabs (
   DenseMtx   *mtx
) {
double   maxabs ;
int      loc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_maxabs(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
if ( DENSEMTX_IS_REAL(mtx) ) {
   maxabs = DVmaxabs(mtx->nrow*mtx->ncol, mtx->entries, &loc) ;
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   maxabs = ZVmaxabs(mtx->nrow*mtx->ncol, mtx->entries) ;
} else {
   fprintf(stderr, "\n fatal error in DenseMtx_maxabs(%p)"
           "\n bad type %d\n", mtx, mtx->type) ;
   exit(-1) ;
}
return(maxabs) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   subtract one matrix from another, B := B - A

   created -- 98may25, cca
   --------------------------------------------
*/
void
DenseMtx_sub (
   DenseMtx   *mtxB,
   DenseMtx   *mtxA
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtxB == NULL || mtxA == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_sub(%p,%p)"
           "\n bad input\n", mtxB, mtxA) ;
   exit(-1) ;
}
if ( mtxB->type != mtxA->type ) {
   fprintf(stderr, "\n fatal error in DenseMtx_sub(%p,%p)"
           "\n mtxB->type = %d, mtxA->type = %d\n", 
           mtxB, mtxA, mtxB->type, mtxA->type) ;
   exit(-1) ;
}
if ( mtxB->nrow != mtxA->nrow ) {
   fprintf(stderr, "\n fatal error in DenseMtx_sub(%p,%p)"
           "\n mtxB->nrow = %d, mtxA->nrow = %d\n", 
           mtxB, mtxA, mtxB->nrow, mtxA->nrow) ;
   exit(-1) ;
}
if ( mtxB->ncol != mtxA->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_sub(%p,%p)"
           "\n mtxB->ncol = %d, mtxA->ncol = %d\n", 
           mtxB, mtxA, mtxB->ncol, mtxA->ncol) ;
   exit(-1) ;
}
if ( mtxB->entries == NULL || mtxA->entries == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_sub(%p,%p)"
           "\n mtxB->entries = %p, mtxA->entries = %p\n", 
           mtxB, mtxA, mtxB->entries, mtxA->entries) ;
   exit(-1) ;
}
if ( DENSEMTX_IS_REAL(mtxB) ) {
   DVsub(mtxB->nrow*mtxB->ncol, mtxB->entries, mtxA->entries) ;
} else if ( DENSEMTX_IS_COMPLEX(mtxB) ) {
   ZVsub(mtxB->nrow*mtxB->ncol, mtxB->entries, mtxA->entries) ;
} else {
   fprintf(stderr, "\n fatal error in DenseMtx_sub(%p,%p)"
           "\n mtxB->type = %d, mtxA->type = %d\n", 
           mtxB, mtxA, mtxB->type, mtxA->type) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to copy a row of the matrix into a vector

   irow -- local row id
   vec  -- double vector to receive the row entries

   created -- 98jul31, cca
   ----------------------------------------------------
*/
void
DenseMtx_copyRowIntoVector (
   DenseMtx   *mtx,
   int        irow,
   double     *vec
) {
double   *entries ;
int      inc1, inc2, jcol, jj, kk, nrow, ncol ;
int      *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || irow < 0 || vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_copyRowIntoVector()"
           "\n bad input\n") ;
   exit(-1) ;
}
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
if ( irow >= nrow ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_copyRowIntoVector()"
           "\n irow = %d, nrow = %d\n", irow, nrow) ;
   exit(-1) ;
}
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
inc1    = DenseMtx_rowIncrement(mtx) ;
inc2    = DenseMtx_columnIncrement(mtx) ;
entries = DenseMtx_entries(mtx) ;
if ( DENSEMTX_IS_REAL(mtx) ) {
   for ( jcol = jj = 0, kk = irow*inc1 ; 
         jcol < ncol ; 
         jcol++, jj++, kk += inc2 ) {
      vec[jj] = entries[kk] ;
   }
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   for ( jcol = jj = 0, kk = irow*inc1 ; 
         jcol < ncol ; 
         jcol++, jj++, kk += inc2 ) {
      vec[2*jj]   = entries[2*kk]   ;
      vec[2*jj+1] = entries[2*kk+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to copy a row of the matrix into a vector

   irow -- local row id
   vec  -- double vector to supply the row entries

   created -- 98jul31, cca
   ----------------------------------------------------
*/
void
DenseMtx_copyVectorIntoRow (
   DenseMtx   *mtx,
   int        irow,
   double     *vec
) {
double   *entries ;
int      inc1, inc2, jcol, jj, kk, nrow, ncol ;
int      *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || irow < 0 || vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_copyVectorIntoRow()"
           "\n bad input, mtx %p, irow %d, vec %p\n",
           mtx, irow, vec) ;
   exit(-1) ;
}
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
if ( irow >= nrow ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_copyVectorIntoRow()"
           "\n irow = %d, nrow = %d\n", irow, nrow) ;
   exit(-1) ;
}
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
inc1    = DenseMtx_rowIncrement(mtx) ;
inc2    = DenseMtx_columnIncrement(mtx) ;
entries = DenseMtx_entries(mtx) ;
if ( DENSEMTX_IS_REAL(mtx) ) {
   for ( jcol = jj = 0, kk = irow*inc1 ; 
         jcol < ncol ; 
         jcol++, jj++, kk += inc2 ) {
      entries[kk] = vec[jj] ;
   }
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   for ( jcol = jj = 0, kk = irow*inc1 ; 
         jcol < ncol ; 
         jcol++, jj++, kk += inc2 ) {
      entries[2*kk]   = vec[2*jj]   ;
      entries[2*kk+1] = vec[2*jj+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to add a row of the matrix into a vector

   irow -- local row id
   vec  -- double vector to supply the row entries

   created -- 98aug12, cca
   ----------------------------------------------------
*/
void
DenseMtx_addVectorIntoRow (
   DenseMtx   *mtx,
   int        irow,
   double     *vec
) {
double   *entries ;
int      inc1, inc2, jcol, jj, kk, nrow, ncol ;
int      *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || irow < 0 || vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_addVectorIntoRow()"
           "\n bad input, mtx %p, irow %d, vec %p\n",
           mtx, irow, vec) ;
   exit(-1) ;
}
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
if ( irow >= nrow ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_addVectorIntoRow()"
           "\n irow = %d, nrow = %d\n", irow, nrow) ;
   exit(-1) ;
}
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
inc1    = DenseMtx_rowIncrement(mtx) ;
inc2    = DenseMtx_columnIncrement(mtx) ;
entries = DenseMtx_entries(mtx) ;
if ( DENSEMTX_IS_REAL(mtx) ) {
   for ( jcol = jj = 0, kk = irow*inc1 ; 
         jcol < ncol ; 
         jcol++, jj++, kk += inc2 ) {
      entries[kk] += vec[jj] ;
   }
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   for ( jcol = jj = 0, kk = irow*inc1 ; 
         jcol < ncol ; 
         jcol++, jj++, kk += inc2 ) {
      entries[2*kk]   += vec[2*jj]   ;
      entries[2*kk+1] += vec[2*jj+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
