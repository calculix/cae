/*  solveH.c  */

#include "../SubMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
static void solveDenseSubrows ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void solveDenseSubcolumns ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void solveSparseRows ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void solveSparseColumns ( SubMtx *mtxA, SubMtx *mtxB ) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- solve (A^H + I) X = B, 
      where 
        (1) X overwrites B
        (2) A must be strict lower or upper triangular
        (2) A, B and X are complex
        (4) columns(A) = rows(X)
        (5) rows(A)    = rows(B)
        (6) B has mode SUBMTX_DENSE_COLUMNS
        (7) if A is SUBMTX_DENSE_SUBROWS or SUBMTX_SPARSE_ROWS
            then A must be strict lower triangular
        (8) if A is SUBMTX_DENSE_SUBCOLUMNS or SUBMTX_SPARSE_COLUMNS
            then A must be strict upper triangular

   created -- 98may01, cca
   -------------------------------------------------------------
*/
void
SubMtx_solveH (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtxA == NULL || mtxB == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveH(%p,%p)"
           "\n bad input\n", mtxA, mtxB) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_COMPLEX(mtxB) ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveH(%p,%p)"
           "\n mtxB has bad type %d\n", mtxA, mtxB, mtxB->type) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_DENSE_COLUMNS(mtxB) ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveH(%p,%p)"
           "\n mtxB has bad mode %d\n", mtxA, mtxB, mtxB->mode) ;
   exit(-1) ;
}
if ( mtxA->nrow != mtxB->nrow ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveH(%p,%p)"
           "\n mtxA->nrow = %d, mtxB->nrwo = %d\n", 
           mtxA, mtxB, mtxA->nrow, mtxB->nrow) ;
   exit(-1) ;
}
/*
   -------------------------
   switch over the mode of A
   -------------------------
*/
switch ( mtxA->mode ) {
case SUBMTX_DENSE_SUBROWS :
   solveDenseSubrows(mtxA, mtxB) ;
   break ;
case SUBMTX_SPARSE_ROWS :
   solveSparseRows(mtxA, mtxB) ;
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   solveDenseSubcolumns(mtxA, mtxB) ;
   break ;
case SUBMTX_SPARSE_COLUMNS :
   solveSparseColumns(mtxA, mtxB) ;
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_solveH(%p,%p)"
           "\n bad mode %d\n", mtxA, mtxB, mtxA->mode) ;
   exit(-1) ;
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- solve (A^T + I) X = B, where 
     (1) A is strictly upper triangular
     (2) X overwrites B
     (B) B has mode SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   ---------------------------------------
*/
static void
solveDenseSubcolumns (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ai, ar, bi0, bi1, bi2, br0, br1, br2, 
         isum0, isum1, isum2, rsum0, rsum1, rsum2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      first, ii, iloc, inc1, inc2, irowA, jcolB, kk, last, 
         ncolB, nentA, nrowA, nrowB, rloc ;
int      *firstlocsA, *sizesA ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_denseSubcolumnsInfo(mtxA, &nrowA, &nentA, 
                         &firstlocsA, &sizesA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n nentA = %d", nentA) ;
   fflush(stdout) ;
#endif
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + 2*nrowB ;
   colB2 = colB1 + 2*nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n %% jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n %% irowA %d, size %d", irowA, sizesA[irowA]) ;
      fflush(stdout) ;
#endif
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
#if MYDEBUG > 0
         fprintf(stdout, ", first %d, last %d", first, last) ;
         fflush(stdout) ;
#endif
         rsum0 = isum0 = 0.0 ;
         rsum1 = isum1 = 0.0 ;
         rsum2 = isum2 = 0.0 ;
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            ar = entriesA[2*kk] ; ai = entriesA[2*kk+1] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %%   A(%d,%d) = (%12.4e,%12.4e)", 
                    irowA+1, ii+1, ar, ai) ;
            fflush(stdout) ;
#endif
            rloc = 2*ii ; iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
            br2 = colB2[rloc] ; bi2 = colB2[iloc] ;
            rsum0 += ar*br0 + ai*bi0 ; isum0 += ar*bi0 - ai*br0 ;
            rsum1 += ar*br1 + ai*bi1 ; isum1 += ar*bi1 - ai*br1 ;
            rsum2 += ar*br2 + ai*bi2 ; isum2 += ar*bi2 - ai*br2 ;
         }
         rloc = 2*irowA ; iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
         colB1[rloc] -= rsum1 ; colB1[iloc] -= isum1 ;
         colB2[rloc] -= rsum2 ; colB2[iloc] -= isum2 ;
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n %% kk = %d", kk) ;
   fflush(stdout) ;
#endif
   colB0 = colB2 + 2*nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + 2*nrowB ;
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n %% irowA %d, size %d", irowA, sizesA[irowA]) ;
      fflush(stdout) ;
#endif
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
#if MYDEBUG > 0
         fprintf(stdout, ", first %d, last %d", first, last) ;
         fflush(stdout) ;
#endif
         rsum0 = isum0 = 0.0 ;
         rsum1 = isum1 = 0.0 ;
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            ar = entriesA[2*kk] ; ai = entriesA[2*kk+1] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %%   A(%d,%d) = (%12.4e,%12.4e)", 
                    irowA+1, ii+1, ar, ai) ;
            fflush(stdout) ;
#endif
            rloc = 2*ii ; iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
            rsum0 += ar*br0 + ai*bi0 ; isum0 += ar*bi0 - ai*br0 ;
            rsum1 += ar*br1 + ai*bi1 ; isum1 += ar*bi1 - ai*br1 ;
         }
         rloc = 2*irowA ; iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
         colB1[rloc] -= rsum1 ; colB1[iloc] -= isum1 ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n %% kk = %d", kk) ;
      fflush(stdout) ;
#endif
   }
} else if ( jcolB == ncolB - 1 ) {
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n %% irowA %d, size %d", irowA, sizesA[irowA]) ;
      fflush(stdout) ;
#endif
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
#if MYDEBUG > 0
         fprintf(stdout, ", first %d, last %d", first, last) ;
         fflush(stdout) ;
#endif
         rsum0 = isum0 = 0.0 ;
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            ar = entriesA[2*kk] ; ai = entriesA[2*kk+1] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %%   A(%d,%d) = (%12.4e,%12.4e)", 
                    irowA+1, ii+1, ar, ai) ;
            fflush(stdout) ;
#endif
            rloc = 2*ii ; iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            rsum0 += ar*br0 + ai*bi0 ; isum0 += ar*bi0 - ai*br0 ;
         }
         rloc = 2*irowA ; iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n %% kk = %d", kk) ;
      fflush(stdout) ;
#endif
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- solve (A^T + I) X = B, where 
     (1) A is strictly upper triangular
     (2) X overwrites B
     (B) B has mode SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   ---------------------------------------
*/
static void
solveSparseColumns (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ai, ar, bi0, bi1, bi2, br0, br1, br2,
         isum0, isum1, isum2, rsum0, rsum1, rsum2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      ii, iloc, inc1, inc2, irowA, jcolB, jj, kk, 
         ncolB, nentA, nrowA, nrowB, rloc, size ;
int      *indicesA, *sizesA ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_sparseColumnsInfo(mtxA, &nrowA, &nentA, 
                       &sizesA, &indicesA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + 2*nrowB ;
   colB2 = colB1 + 2*nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizesA[irowA]) > 0 ) {
         rsum0 = isum0 = 0.0 ;
         rsum1 = isum1 = 0.0 ;
         rsum2 = isum2 = 0.0 ;
         for ( ii = 0 ; ii < size ; ii++, kk++ ) {
            ar = entriesA[2*kk] ;
            ai = entriesA[2*kk+1] ;
            jj  = indicesA[kk] ;
            if ( jj < 0 || jj >= irowA ) {
               fprintf(stderr, 
            "\n fatal error, irowA = %d, kk =%d, ii = %d, jj = %d",
            irowA, kk, ii, jj) ;
               exit(-1) ;
            }
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
            br2 = colB2[rloc] ; bi2 = colB2[iloc] ;
            rsum0 += ar*br0 + ai*bi0 ; isum0 += ar*bi0 - ai*br0 ;
            rsum1 += ar*br1 + ai*bi1 ; isum1 += ar*bi1 - ai*br1 ;
            rsum2 += ar*br2 + ai*bi2 ; isum2 += ar*bi2 - ai*br2 ;
         }
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
         colB1[rloc] -= rsum1 ; colB1[iloc] -= isum1 ;
         colB2[rloc] -= rsum2 ; colB2[iloc] -= isum2 ;
      }
   }
   colB0 = colB2 + 2*nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + 2*nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizesA[irowA]) > 0 ) {
         rsum0 = isum0 = 0.0 ;
         rsum1 = isum1 = 0.0 ;
         for ( ii = 0 ; ii < size ; ii++, kk++ ) {
            ar = entriesA[2*kk] ;
            ai = entriesA[2*kk+1] ;
            jj  = indicesA[kk] ;
            if ( jj < 0 || jj >= irowA ) {
               fprintf(stderr, 
            "\n fatal error, irowA = %d, kk =%d, ii = %d, jj = %d",
            irowA, kk, ii, jj) ;
               exit(-1) ;
            }
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
            rsum0 += ar*br0 + ai*bi0 ; isum0 += ar*bi0 - ai*br0 ;
            rsum1 += ar*br1 + ai*bi1 ; isum1 += ar*bi1 - ai*br1 ;
         }
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
         colB1[rloc] -= rsum1 ; colB1[iloc] -= isum1 ;
      }
   }
} else if ( jcolB == ncolB - 1 ) {
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizesA[irowA]) > 0 ) {
         rsum0 = isum0 = 0.0 ;
         for ( ii = 0 ; ii < size ; ii++, kk++ ) {
            ar = entriesA[2*kk] ;
            ai = entriesA[2*kk+1] ;
            jj  = indicesA[kk] ;
            if ( jj < 0 || jj >= irowA ) {
               fprintf(stderr, 
            "\n fatal error, irowA = %d, kk =%d, ii = %d, jj = %d",
            irowA, kk, ii, jj) ;
               exit(-1) ;
            }
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            rsum0 += ar*br0 + ai*bi0 ; isum0 += ar*bi0 - ai*br0 ;
         }
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- solve (I + A^T) X = B, where 
     (1) A is strictly lower triangular
     (2) X overwrites B
     (B) B has mode SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   ---------------------------------------
*/
static void
solveDenseSubrows (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ai, ar, bi0, bi1, bi2, br0, br1, br2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      colstart, first, iloc, inc1, inc2, irowA, jcolB, 
         jj, kk, last, ncolB, nentA, nrowA, nrowB, rloc ;
int      *firstlocsA, *sizesA ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_denseSubrowsInfo(mtxA, &nrowA, &nentA, 
                      &firstlocsA, &sizesA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n nrowA = %d, ncolA = %d", nrowA, nentA) ;
   fflush(stdout) ;
#endif
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + 2*nrowB ;
   colB2 = colB1 + 2*nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = nrowA - 1, colstart = nentA ; 
         irowA >= 0 ; 
         irowA-- ) {
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
         colstart -= last - first + 1 ;
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
         br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
         br2 = colB2[rloc] ; bi2 = colB2[iloc] ;
         for ( jj = first, kk = colstart ; jj <= last ; jj++, kk++ ) {
            ar = entriesA[2*kk] ;
            ai = entriesA[2*kk+1] ;
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            colB0[rloc] -= ar*br0 + ai*bi0 ;
            colB0[iloc] -= ar*bi0 - ai*br0 ;
            colB1[rloc] -= ar*br1 + ai*bi1 ;
            colB1[iloc] -= ar*bi1 - ai*br1 ;
            colB2[rloc] -= ar*br2 + ai*bi2 ;
            colB2[iloc] -= ar*bi2 - ai*br2 ;
         }
      }
   }
   colB0 = colB2 + 2*nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + 2*nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = nrowA - 1, colstart = nentA ; 
         irowA >= 0 ; 
         irowA-- ) {
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
         colstart -= last - first + 1 ;
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
         br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
         for ( jj = first, kk = colstart ; jj <= last ; jj++, kk++ ) {
            ar = entriesA[2*kk] ;
            ai = entriesA[2*kk+1] ;
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            colB0[rloc] -= ar*br0 + ai*bi0 ;
            colB0[iloc] -= ar*bi0 - ai*br0 ;
            colB1[rloc] -= ar*br1 + ai*bi1 ;
            colB1[iloc] -= ar*bi1 - ai*br1 ;
         }
      }
   }
} else if ( jcolB == ncolB - 1 ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = nrowA - 1, colstart = nentA ; 
         irowA >= 0 ; 
         irowA-- ) {
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
         colstart -= last - first + 1 ;
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
         for ( jj = first, kk = colstart ; jj <= last ; jj++, kk++ ) {
            ar = entriesA[2*kk] ;
            ai = entriesA[2*kk+1] ;
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            colB0[rloc] -= ar*br0 + ai*bi0 ;
            colB0[iloc] -= ar*bi0 - ai*br0 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- solve (I + A^T) X = B, where 
     (1) A is strictly lower triangular
     (2) X overwrites B
     (B) B has mode SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   ---------------------------------------
*/
static void
solveSparseRows (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ai, ar, bi0, bi1, bi2, br0, br1, br2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      colstart, ii, iloc, inc1, inc2, jcolA, jcolB, 
         jj, kk, ncolB, nentA, nrowA, nrowB, rloc, size ;
int      *indicesA, *sizesA ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_sparseRowsInfo(mtxA, &nrowA, &nentA, 
                    &sizesA, &indicesA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n nrowA = %d, ncolA = %d", nrowA, nentA) ;
   fflush(stdout) ;
#endif
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + 2*nrowB ;
   colB2 = colB1 + 2*nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( jcolA = nrowA - 1, colstart = nentA ; 
         jcolA >= 0 ; 
         jcolA-- ) {
      if ( (size = sizesA[jcolA]) > 0 ) {
         colstart -= size ;
         rloc = 2*jcolA ;
         iloc = rloc + 1 ;
         br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
         br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
         br2 = colB2[rloc] ; bi2 = colB2[iloc] ;
         for ( ii = 0, kk = colstart ; ii < size ; ii++, kk++ ) {
            ar = entriesA[2*kk] ; ai = entriesA[2*kk+1] ;
            jj  = indicesA[kk] ;
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            colB0[rloc] -= ar*br0 + ai*bi0 ;
            colB0[iloc] -= ar*bi0 - ai*br0 ;
            colB1[rloc] -= ar*br1 + ai*bi1 ;
            colB1[iloc] -= ar*bi1 - ai*br1 ;
            colB2[rloc] -= ar*br2 + ai*bi2 ;
            colB2[iloc] -= ar*bi2 - ai*br2 ;
         }
      }
   }
   colB0 = colB2 + 2*nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + 2*nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( jcolA = nrowA - 1, colstart = nentA ; 
         jcolA >= 0 ; 
         jcolA-- ) {
      if ( (size = sizesA[jcolA]) > 0 ) {
         colstart -= size ;
         rloc = 2*jcolA ;
         iloc = rloc + 1 ;
         br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
         br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
         for ( ii = 0, kk = colstart ; ii < size ; ii++, kk++ ) {
            ar = entriesA[2*kk] ; ai = entriesA[2*kk+1] ;
            jj  = indicesA[kk] ;
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            colB0[rloc] -= ar*br0 + ai*bi0 ;
            colB0[iloc] -= ar*bi0 - ai*br0 ;
            colB1[rloc] -= ar*br1 + ai*bi1 ;
            colB1[iloc] -= ar*bi1 - ai*br1 ;
         }
      }
   }
} else if ( jcolB == ncolB - 1 ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( jcolA = nrowA - 1, colstart = nentA ; 
         jcolA >= 0 ; 
         jcolA-- ) {
      if ( (size = sizesA[jcolA]) > 0 ) {
         colstart -= size ;
         rloc = 2*jcolA ;
         iloc = rloc + 1 ;
         br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
         for ( ii = 0, kk = colstart ; ii < size ; ii++, kk++ ) {
            ar = entriesA[2*kk] ; ai = entriesA[2*kk+1] ;
            jj  = indicesA[kk] ;
            rloc = 2*jj ;
            iloc = rloc + 1 ;
            colB0[rloc] -= ar*br0 + ai*bi0 ;
            colB0[iloc] -= ar*bi0 - ai*br0 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
