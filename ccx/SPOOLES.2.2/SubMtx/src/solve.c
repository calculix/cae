/*  solve.c  */

#include "../SubMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
static void real_solveDenseSubrows ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void real_solveDenseSubcolumns ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void real_solveSparseRows ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void real_solveSparseColumns ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void real_solveDiagonal ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void real_solveBlockDiagonalSym ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void complex_solveDenseSubrows ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void complex_solveDenseSubcolumns ( SubMtx *mtxA, SubMtx *mtxB );
static void complex_solveSparseRows ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void complex_solveSparseColumns ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void complex_solveDiagonal ( SubMtx *mtxA, SubMtx *mtxB ) ;
static void complex_solveBlockDiagonalSym ( SubMtx *mtxA, SubMtx *mtxB);
static void complex_solveBlockDiagonalHerm (SubMtx *mtxA, SubMtx *mtxB);
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- solve A X = B, 
      where 
        (1) X overwrites B
        (2) A must be lower or upper triangular, 
            diagonal or block diagonal
        (3) if A is strict lower or upper triangular,
            then we solve (I + A) X = B
        (4) columns(A) = rows(X)
        (5) rows(A)    = rows(B)
        (6) B has type SUBMTX_DENSE_COLUMNS
        (7) if A is SUBMTX_DENSE_SUBROWS or SUBMTX_SPARSE_ROWS
            then A must be strict lower triangular
        (8) if A is SUBMTX_DENSE_SUBCOLUMNS or SUBMTX_SPARSE_COLUMNS
            then A must be strict upper triangular
        (9) A can be SUBMTX_DIAGONAL, SUBMTX_BLOCK_DIAGONAL_SYM
            or SUBMTX_BLOCK_DIAGONAL_HERM

   created -- 98may01, cca
   -----------------------------------------------------------------
*/
void
SubMtx_solve (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtxA == NULL || mtxB == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
           "\n bad input\n", mtxA, mtxB) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_DENSE_COLUMNS(mtxB) ) {
   fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
           "\n mtxB has bad type %d\n", mtxA, mtxB, mtxB->type) ;
   exit(-1) ;
}
if ( mtxA->nrow != mtxB->nrow ) {
   fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
           "\n mtxA->nrow = %d, mtxB->nrwo = %d\n", 
           mtxA, mtxB, mtxA->nrow, mtxB->nrow) ;
   exit(-1) ;
}
/*
   -------------------------
   switch over the type of A
   -------------------------
*/
switch ( mtxA->type ) {
case SPOOLES_REAL :
/*
   -------------------------
   switch over the mode of A
   -------------------------
*/
   switch ( mtxA->mode ) {
   case SUBMTX_DENSE_SUBROWS :
      real_solveDenseSubrows(mtxA, mtxB) ;
      break ;
   case SUBMTX_SPARSE_ROWS :
      real_solveSparseRows(mtxA, mtxB) ;
      break ;
   case SUBMTX_DENSE_SUBCOLUMNS :
      real_solveDenseSubcolumns(mtxA, mtxB) ;
      break ;
   case SUBMTX_SPARSE_COLUMNS :
      real_solveSparseColumns(mtxA, mtxB) ;
      break ;
   case SUBMTX_DIAGONAL :
      real_solveDiagonal(mtxA, mtxB) ;
      break ;
   case SUBMTX_BLOCK_DIAGONAL_SYM :
      real_solveBlockDiagonalSym(mtxA, mtxB) ;
      break ;
   case SUBMTX_DENSE_ROWS :
   case SUBMTX_DENSE_COLUMNS :
   case SUBMTX_SPARSE_TRIPLES :
   default :
      fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
              "\n bad mode %d\n", mtxA, mtxB, mtxA->type) ;
      exit(-1) ;
      break ;
   }
   break ;
case SPOOLES_COMPLEX :
/*
   -------------------------
   switch over the mode of A
   -------------------------
*/
   switch ( mtxA->mode ) {
   case SUBMTX_DENSE_SUBROWS :
      complex_solveDenseSubrows(mtxA, mtxB) ;
      break ;
   case SUBMTX_SPARSE_ROWS :
      complex_solveSparseRows(mtxA, mtxB) ;
      break ;
   case SUBMTX_DENSE_SUBCOLUMNS :
      complex_solveDenseSubcolumns(mtxA, mtxB) ;
      break ;
   case SUBMTX_SPARSE_COLUMNS :
      complex_solveSparseColumns(mtxA, mtxB) ;
      break ;
   case SUBMTX_DIAGONAL :
      complex_solveDiagonal(mtxA, mtxB) ;
      break ;
   case SUBMTX_BLOCK_DIAGONAL_SYM :
      complex_solveBlockDiagonalSym(mtxA, mtxB) ;
      break ;
   case SUBMTX_BLOCK_DIAGONAL_HERM :
      complex_solveBlockDiagonalHerm(mtxA, mtxB) ;
      break ;
   case SUBMTX_DENSE_ROWS :
   case SUBMTX_DENSE_COLUMNS :
   case SUBMTX_SPARSE_TRIPLES :
   default :
      fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
              "\n bad mode %d\n", mtxA, mtxB, mtxA->type) ;
      exit(-1) ;
      break ;
   }
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
           "\n bad type %d\n", mtxA, mtxB, mtxA->type) ;
   exit(-1) ;
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- solve (A + I) X = B, where 
     (1) A is strictly lower triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
real_solveDenseSubrows (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   Aki, sum0, sum1, sum2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      first, ii, inc1, inc2, irowA, jcolB, kk, last, 
         ncolB, nentA, nrowA, nrowB ;
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
   fprintf(stdout, "\n nentA = %d", nentA) ;
   fflush(stdout) ;
#endif
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + nrowB ;
   colB2 = colB1 + nrowB ;
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
         sum0 = sum1 = sum2 = 0.0 ;
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            Aki = entriesA[kk] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %%   Aki A(%d,%d) = %12.4e", 
                    irowA+1, ii+1, Aki) ;
            fflush(stdout) ;
#endif
            sum0 += Aki * colB0[ii] ;
            sum1 += Aki * colB1[ii] ;
            sum2 += Aki * colB2[ii] ;
         }
         colB0[irowA] -= sum0 ;
         colB1[irowA] -= sum1 ;
         colB2[irowA] -= sum2 ;
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n %% kk = %d", kk) ;
   fflush(stdout) ;
#endif
   colB0 = colB2 + nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n %% jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
         sum0 = sum1 = 0.0 ;
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            Aki = entriesA[kk] ;
            sum0 += Aki * colB0[ii] ;
            sum1 += Aki * colB1[ii] ;
         }
         colB0[irowA] -= sum0 ;
         colB1[irowA] -= sum1 ;
      }
   }
} else if ( jcolB == ncolB - 1 ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n %% jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n %% irowA = %d, kk = %d, sizesA[%d] = %d",
           irowA, kk, irowA, sizesA[irowA]) ;
   fflush(stdout) ;
#endif
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
#if MYDEBUG > 0
         fprintf(stdout, "\n %% first = %d, last = %d", first, last) ;
         fflush(stdout) ;
#endif
         sum0 = 0.0 ;
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            Aki = entriesA[kk] ;
#if MYDEBUG > 0
         fprintf(stdout, "\n %% Aki = %12.4e, colB0[%d] = %12.4e",
                 Aki, ii, colB0[ii]) ;
         fflush(stdout) ;
#endif
            sum0 += Aki * colB0[ii] ;
         }
         colB0[irowA] -= sum0 ;
#if MYDEBUG > 0
         fprintf(stdout, "\n %% colB0[%d] -= %12.4e",
                 irowA, sum0) ;
         fflush(stdout) ;
#endif
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- solve (A + I) X = B, where 
     (1) A is strictly lower triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
real_solveSparseRows (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   Aki, sum0, sum1, sum2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      ii, inc1, inc2, irowA, jcolB, jj, kk, 
         ncolB, nentA, nrowA, nrowB, size ;
int      *indicesA, *sizesA ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_sparseRowsInfo(mtxA, &nrowA, &nentA, 
                    &sizesA, &indicesA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + nrowB ;
   colB2 = colB1 + nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizesA[irowA]) > 0 ) {
         sum0 = sum1 = sum2 = 0.0 ;
         for ( ii = 0 ; ii < size ; ii++, kk++ ) {
            Aki = entriesA[kk] ;
            jj  = indicesA[kk] ;
if ( jj < 0 || jj >= irowA ) {
   fprintf(stderr, 
"\n fatal error, irowA = %d, kk =%d, ii = %d, jj = %d",
irowA, kk, ii, jj) ;
   exit(-1) ;
}
            sum0 += Aki * colB0[jj] ;
            sum1 += Aki * colB1[jj] ;
            sum2 += Aki * colB2[jj] ;
         }
         colB0[irowA] -= sum0 ;
         colB1[irowA] -= sum1 ;
         colB2[irowA] -= sum2 ;
      }
   }
   colB0 = colB2 + nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizesA[irowA]) > 0 ) {
         sum0 = sum1 = 0.0 ;
         for ( ii = 0 ; ii < size ; ii++, kk++ ) {
            Aki = entriesA[kk] ;
            jj  = indicesA[kk] ;
            sum0 += Aki * colB0[jj] ;
            sum1 += Aki * colB1[jj] ;
         }
         colB0[irowA] -= sum0 ;
         colB1[irowA] -= sum1 ;
      }
   }
} else if ( jcolB == ncolB - 1 ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizesA[irowA]) >= 0 ) {
         sum0 = 0.0 ;
         for ( ii = 0 ; ii < size ; ii++, kk++ ) {
            Aki = entriesA[kk] ;
            jj  = indicesA[kk] ;
if ( jj < 0 || jj >= irowA ) {
   fprintf(stderr, 
"\n fatal error, irowA = %d, kk =%d, ii = %d, jj = %d",
irowA, kk, ii, jj) ;
   exit(-1) ;
}
            sum0 += Aki * colB0[jj] ;
         }
         colB0[irowA] -= sum0 ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- solve (I + A) X = B, where 
     (1) A is strictly upper triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
real_solveDenseSubcolumns (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   Aji, Bi0, Bi1, Bi2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      colstart, first, inc1, inc2, irowA, jcolB, 
         jj, kk, last, ncolB, nentA, nrowA, nrowB ;
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
   fprintf(stdout, "\n nrowA = %d, ncolA = %d", nrowA, nentA) ;
   fflush(stdout) ;
#endif
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + nrowB ;
   colB2 = colB1 + nrowB ;
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
         Bi0 = colB0[irowA] ;
         Bi1 = colB1[irowA] ;
         Bi2 = colB2[irowA] ;
         for ( jj = first, kk = colstart ; jj <= last ; jj++, kk++ ) {
            Aji = entriesA[kk] ;
            colB0[jj] -= Aji * Bi0 ;
            colB1[jj] -= Aji * Bi1 ;
            colB2[jj] -= Aji * Bi2 ;
         }
      }
   }
   colB0 = colB2 + nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + nrowB ;
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
         Bi0 = colB0[irowA] ;
         Bi1 = colB1[irowA] ;
         for ( jj = first, kk = colstart ; jj <= last ; jj++, kk++ ) {
            Aji = entriesA[kk] ;
            colB0[jj] -= Aji * Bi0 ;
            colB1[jj] -= Aji * Bi1 ;
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
#if MYDEBUG > 0
      fprintf(stdout, "\n %% irowA = %d, sizesA[%d] = %d", 
              irowA, irowA, sizesA[irowA]) ;
      fflush(stdout) ;
#endif
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
         colstart -= last - first + 1 ;
         Bi0 = colB0[irowA] ;
#if MYDEBUG > 0
         fprintf(stdout, 
                 "\n %% first %d, last %d, colstart %d, Bi0 = %12.4e",
                 first, last, colstart, Bi0) ;
         fflush(stdout) ;
#endif
         for ( jj = first, kk = colstart ; jj <= last ; jj++, kk++ ) {
            Aji = entriesA[kk] ;
#if MYDEBUG > 0
            fprintf(stdout, 
                 "\n %% jj %d, kk %d, Aji %12.4e", jj, kk, Aji) ;
            fflush(stdout) ;
#endif
            colB0[jj] -= Aji * Bi0 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- solve (I + A) X = B, where 
     (1) A is strictly upper triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
real_solveSparseColumns (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   Aji, Bi0, Bi1, Bi2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      colstart, ii, inc1, inc2, jcolA, jcolB, 
         jj, kk, ncolB, nentA, nrowA, nrowB, size ;
int      *indicesA, *sizesA ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_sparseColumnsInfo(mtxA, &nrowA, &nentA, 
                       &sizesA, &indicesA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n nrowA = %d, ncolA = %d", nrowA, nentA) ;
   fflush(stdout) ;
#endif
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + nrowB ;
   colB2 = colB1 + nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( jcolA = nrowA - 1, colstart = nentA ; 
         jcolA >= 0 ; 
         jcolA-- ) {
      if ( (size = sizesA[jcolA]) > 0 ) {
         colstart -= size ;
         Bi0 = colB0[jcolA] ;
         Bi1 = colB1[jcolA] ;
         Bi2 = colB2[jcolA] ;
         for ( ii = 0, kk = colstart ; ii < size ; ii++, kk++ ) {
            Aji = entriesA[kk] ;
            jj  = indicesA[kk] ;
            colB0[jj] -= Aji * Bi0 ;
            colB1[jj] -= Aji * Bi1 ;
            colB2[jj] -= Aji * Bi2 ;
         }
      }
   }
   colB0 = colB2 + nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + nrowB ;
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( jcolA = nrowA - 1, colstart = nentA ; 
         jcolA >= 0 ; 
         jcolA-- ) {
      if ( (size = sizesA[jcolA]) > 0 ) {
         colstart -= size ;
         Bi0 = colB0[jcolA] ;
         Bi1 = colB1[jcolA] ;
         for ( ii = 0, kk = colstart ; ii < size ; ii++, kk++ ) {
            Aji = entriesA[kk] ;
            jj  = indicesA[kk] ;
            colB0[jj] -= Aji * Bi0 ;
            colB1[jj] -= Aji * Bi1 ;
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
         Bi0 = colB0[jcolA] ;
         for ( ii = 0, kk = colstart ; ii < size ; ii++, kk++ ) {
            Aji = entriesA[kk] ;
            jj  = indicesA[kk] ;
            colB0[jj] -= Aji * Bi0 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- solve A X = B, where 
     (1) A is diagonal
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -----------------------------------
*/
static void
real_solveDiagonal (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   Aii ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      inc1, inc2, irowA, jcolB, ncolB, nrowA, nrowB ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_diagonalInfo(mtxA, &nrowA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + nrowB ;
   colB2 = colB1 + nrowB ;
   for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
      Aii = entriesA[irowA] ;
      colB0[irowA] /= Aii ;
      colB1[irowA] /= Aii ;
      colB2[irowA] /= Aii ;
   }
   colB0 = colB2 + nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + nrowB ;
   for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
      Aii = entriesA[irowA] ;
      colB0[irowA] /= Aii ;
      colB1[irowA] /= Aii ;
   }
} else if ( jcolB == ncolB - 1 ) {
   for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
      Aii = entriesA[irowA] ;
      colB0[irowA] /= Aii ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- solve A X = B, where 
     (1) A is block diagonal symmetric
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -----------------------------------
*/
static void
real_solveBlockDiagonalSym (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   Aii, Arr, Ars, Ass, recip, t1, t2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      inc1, inc2, ipivot, irowA, jcolB, kk, m, 
         ncolB, nentA, nrowA, nrowB ;
int      *pivotsizes ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entriesA);
for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
   m = pivotsizes[ipivot] ;
   if ( m != 1 && m != 2 ) {
      fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
              "\n mtxA is block diagonal symmetric"
              "\n pivot %d is %d, not supported",
              mtxA, mtxB, ipivot, m) ;
      exit(-1) ;
   }
   irowA += m ;
}
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + nrowB ;
   colB2 = colB1 + nrowB ;
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
         Aii = entriesA[kk++] ;
         colB0[irowA] /= Aii ;
         colB1[irowA] /= Aii ;
         colB2[irowA] /= Aii ;
      } else if ( m == 2 ) {
         Arr = entriesA[kk++] ;
         Ars = entriesA[kk++] ;
         Ass = entriesA[kk++] ;
         recip = 1./(Arr*Ass - Ars*Ars) ;
         t1 = colB0[irowA] ;
         t2 = colB0[irowA+1] ;
         colB0[irowA]   = recip*( Ass*t1 - Ars*t2) ;
         colB0[irowA+1] = recip*(-Ars*t1 + Arr*t2) ;
         t1 = colB1[irowA] ;
         t2 = colB1[irowA+1] ;
         colB1[irowA]   = recip*( Ass*t1 - Ars*t2) ;
         colB1[irowA+1] = recip*(-Ars*t1 + Arr*t2) ;
         t1 = colB2[irowA] ;
         t2 = colB2[irowA+1] ;
         colB2[irowA]   = recip*( Ass*t1 - Ars*t2) ;
         colB2[irowA+1] = recip*(-Ars*t1 + Arr*t2) ;
      }
      irowA += m ;
   }
   colB0 = colB2 + nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + nrowB ;
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
         Aii = entriesA[kk++] ;
         colB0[irowA] /= Aii ;
         colB1[irowA] /= Aii ;
      } else if ( m == 2 ) {
         Arr = entriesA[kk++] ;
         Ars = entriesA[kk++] ;
         Ass = entriesA[kk++] ;
         recip = 1./(Arr*Ass - Ars*Ars) ;
         t1 = colB0[irowA] ;
         t2 = colB0[irowA+1] ;
         colB0[irowA]   = recip*( Ass*t1 - Ars*t2) ;
         colB0[irowA+1] = recip*(-Ars*t1 + Arr*t2) ;
         t1 = colB1[irowA] ;
         t2 = colB1[irowA+1] ;
         colB1[irowA]   = recip*( Ass*t1 - Ars*t2) ;
         colB1[irowA+1] = recip*(-Ars*t1 + Arr*t2) ;
      }
      irowA += m ;
   }
} else if ( jcolB == ncolB - 1 ) {
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
         Aii = entriesA[kk++] ;
         colB0[irowA] /= Aii ;
      } else if ( m == 2 ) {
         Arr = entriesA[kk++] ;
         Ars = entriesA[kk++] ;
         Ass = entriesA[kk++] ;
         recip = 1./(Arr*Ass - Ars*Ars) ;
         t1 = colB0[irowA] ;
         t2 = colB0[irowA+1] ;
         colB0[irowA]   = recip*( Ass*t1 - Ars*t2) ;
         colB0[irowA+1] = recip*(-Ars*t1 + Arr*t2) ;
      }
      irowA += m ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- solve (A + I) X = B, where 
     (1) A is strictly lower triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
complex_solveDenseSubrows (
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
SubMtx_denseSubrowsInfo(mtxA, &nrowA, &nentA, 
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
            rloc = 2*kk ;
            iloc = rloc + 1 ;
            ar = entriesA[rloc] ;
            ai = entriesA[iloc] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %%   A(%d,%d) = (%12.4e,%12.4e)", 
                    irowA+1, ii+1, ar, ai) ;
            fflush(stdout) ;
#endif
            rloc = 2*ii ;
            iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
            br2 = colB2[rloc] ; bi2 = colB2[iloc] ;
            rsum0 += ar*br0 - ai*bi0 ; isum0 += ar*bi0 + ai*br0 ;
            rsum1 += ar*br1 - ai*bi1 ; isum1 += ar*bi1 + ai*br1 ;
            rsum2 += ar*br2 - ai*bi2 ; isum2 += ar*bi2 + ai*br2 ;
         }
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
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
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            rloc = 2*kk ;
            iloc = rloc + 1 ;
            ar = entriesA[rloc] ;
            ai = entriesA[iloc] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %%   A(%d,%d) = (%12.4e,%12.4e)", 
                    irowA+1, ii+1, ar, ai) ;
            fflush(stdout) ;
#endif
            rloc = 2*ii ;
            iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            br1 = colB1[rloc] ; bi1 = colB1[iloc] ;
            rsum0 += ar*br0 - ai*bi0 ; isum0 += ar*bi0 + ai*br0 ;
            rsum1 += ar*br1 - ai*bi1 ; isum1 += ar*bi1 + ai*br1 ;
         }
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
         colB1[rloc] -= rsum1 ; colB1[iloc] -= isum1 ;
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n %% kk = %d", kk) ;
   fflush(stdout) ;
#endif
} else if ( jcolB == ncolB - 1 ) {
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
         for ( ii = first ; ii <= last ; ii++, kk++ ) {
            rloc = 2*kk ;
            iloc = rloc + 1 ;
            ar = entriesA[rloc] ;
            ai = entriesA[iloc] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n %%   A(%d,%d) = (%12.4e,%12.4e)", 
                    irowA+1, ii+1, ar, ai) ;
            fflush(stdout) ;
#endif
            rloc = 2*ii ;
            iloc = rloc + 1 ;
            br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
            rsum0 += ar*br0 - ai*bi0 ; isum0 += ar*bi0 + ai*br0 ;
         }
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n %% kk = %d", kk) ;
   fflush(stdout) ;
#endif
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- solve (A + I) X = B, where 
     (1) A is strictly lower triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
complex_solveSparseRows (
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
SubMtx_sparseRowsInfo(mtxA, &nrowA, &nentA, 
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
            rloc = 2*kk ;
            iloc = rloc + 1 ;
            ar = entriesA[rloc] ;
            ai = entriesA[iloc] ;
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
            rsum0 += ar*br0 - ai*bi0 ; isum0 += ar*bi0 + ai*br0 ;
            rsum1 += ar*br1 - ai*bi1 ; isum1 += ar*bi1 + ai*br1 ;
            rsum2 += ar*br2 - ai*bi2 ; isum2 += ar*bi2 + ai*br2 ;
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
            rloc = 2*kk ;
            iloc = rloc + 1 ;
            ar = entriesA[rloc] ;
            ai = entriesA[iloc] ;
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
            rsum0 += ar*br0 - ai*bi0 ; isum0 += ar*bi0 + ai*br0 ;
            rsum1 += ar*br1 - ai*bi1 ; isum1 += ar*bi1 + ai*br1 ;
         }
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
         colB0[rloc] -= rsum0 ; colB0[iloc] -= isum0 ;
         colB1[rloc] -= rsum1 ; colB1[iloc] -= isum1 ;
      }
   }
} else if ( jcolB == ncolB - 1 ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n jcolB = %d", jcolB) ;
   fflush(stdout) ;
#endif
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizesA[irowA]) > 0 ) {
         rsum0 = isum0 = 0.0 ;
         for ( ii = 0 ; ii < size ; ii++, kk++ ) {
            rloc = 2*kk ;
            iloc = rloc + 1 ;
            ar = entriesA[rloc] ;
            ai = entriesA[iloc] ;
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
            rsum0 += ar*br0 - ai*bi0 ; isum0 += ar*bi0 + ai*br0 ;
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
   -------------------------------------
   purpose -- solve (I + A) X = B, where 
     (1) A is strictly upper triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
complex_solveDenseSubcolumns (
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
SubMtx_denseSubcolumnsInfo(mtxA, &nrowA, &nentA, 
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
            colB0[rloc] -= ar*br0 - ai*bi0 ;
            colB0[iloc] -= ar*bi0 + ai*br0 ;
            colB1[rloc] -= ar*br1 - ai*bi1 ;
            colB1[iloc] -= ar*bi1 + ai*br1 ;
            colB2[rloc] -= ar*br2 - ai*bi2 ;
            colB2[iloc] -= ar*bi2 + ai*br2 ;
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
            colB0[rloc] -= ar*br0 - ai*bi0 ;
            colB0[iloc] -= ar*bi0 + ai*br0 ;
            colB1[rloc] -= ar*br1 - ai*bi1 ;
            colB1[iloc] -= ar*bi1 + ai*br1 ;
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
/*
fprintf(stdout, "\n %% irowA %d, size %d", irowA, sizesA[irowA]) ;
*/
      if ( sizesA[irowA] > 0 ) {
         first = firstlocsA[irowA] ;
         last  = first + sizesA[irowA] - 1 ;
         colstart -= last - first + 1 ;
         rloc = 2*irowA ;
         iloc = rloc + 1 ;
/*
fprintf(stdout, 
        "\n %% first %d, last %d, colstart %d, rloc %d, iloc %d",
        first, last, colstart, rloc, iloc) ;
*/
         br0 = colB0[rloc] ; bi0 = colB0[iloc] ;
/*
fprintf(stdout, "\n %% br0 %12.4e, bi0 %12.4e", br0, bi0) ;
*/
         for ( jj = first, kk = colstart ; jj <= last ; jj++, kk++ ) {
            ar = entriesA[2*kk] ;
            ai = entriesA[2*kk+1] ;
/*
fprintf(stdout, "\n %% ar %12.4e, ai %12.4e", ar, ai) ;
*/
            rloc = 2*jj ;
            iloc = rloc + 1 ;
/*
fprintf(stdout, "\n %% rloc %d, iloc %d", rloc, iloc) ;
*/
            colB0[rloc] -= ar*br0 - ai*bi0 ;
            colB0[iloc] -= ar*bi0 + ai*br0 ;
/*
fprintf(stdout, "\n %% colB[%d] = %12.5e, colB[%d] = %12.5e",
        rloc, colB0[rloc], iloc, colB0[iloc]) ;
*/
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- solve (I + A) X = B, where 
     (1) A is strictly upper triangular
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -------------------------------------
*/
static void
complex_solveSparseColumns (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ar, ai, bi0, bi1, bi2, br0, br1, br2 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      colstart, ii, iloc, inc1, inc2, jcolA, jcolB, 
         jj, kk, ncolB, nentA, nrowA, nrowB, rloc, size ;
int      *indicesA, *sizesA ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_sparseColumnsInfo(mtxA, &nrowA, &nentA, 
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
            colB0[rloc] -= ar*br0 - ai*bi0 ;
            colB0[iloc] -= ar*bi0 + ai*br0 ;
            colB1[rloc] -= ar*br1 - ai*bi1 ;
            colB1[iloc] -= ar*bi1 + ai*br1 ;
            colB2[rloc] -= ar*br2 - ai*bi2 ;
            colB2[iloc] -= ar*bi2 + ai*br2 ;
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
            colB0[rloc] -= ar*br0 - ai*bi0 ;
            colB0[iloc] -= ar*bi0 + ai*br0 ;
            colB1[rloc] -= ar*br1 - ai*bi1 ;
            colB1[iloc] -= ar*bi1 + ai*br1 ;
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
            colB0[rloc] -= ar*br0 - ai*bi0 ;
            colB0[iloc] -= ar*bi0 + ai*br0 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- solve A X = B, where 
     (1) A is diagonal
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -----------------------------------
*/
static void
complex_solveDiagonal (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ai, ar, bi, br, ci, cr ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      inc1, inc2, irowA, jcolB, kk, ncolB, nrowA, nrowB ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_diagonalInfo(mtxA, &nrowA, &entriesA) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + 2*nrowB ;
   colB2 = colB1 + 2*nrowB ;
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++, kk += 2 ) {
      ar = entriesA[kk] ;
      ai = entriesA[kk+1] ;
      Zrecip(ar, ai, &cr, &ci) ;
      br = colB0[kk] ; bi = colB0[kk+1] ;
      colB0[kk]   = br*cr - bi*ci ; 
      colB0[kk+1] = br*ci + bi*cr ;
      br = colB1[kk] ; bi = colB1[kk+1] ;
      colB1[kk]   = br*cr - bi*ci ; 
      colB1[kk+1] = br*ci + bi*cr ;
      br = colB2[kk] ; bi = colB2[kk+1] ;
      colB2[kk]   = br*cr - bi*ci ; 
      colB2[kk+1] = br*ci + bi*cr ;
   }
   colB0 = colB2 + 2*nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + 2*nrowB ;
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++, kk += 2 ) {
      ar = entriesA[kk] ;
      ai = entriesA[kk+1] ;
      Zrecip(ar, ai, &cr, &ci) ;
      br = colB0[kk] ; bi = colB0[kk+1] ;
      colB0[kk]   = br*cr - bi*ci ; 
      colB0[kk+1] = br*ci + bi*cr ;
      br = colB1[kk] ; bi = colB1[kk+1] ;
      colB1[kk]   = br*cr - bi*ci ; 
      colB1[kk+1] = br*ci + bi*cr ;
   }
} else if ( jcolB == ncolB - 1 ) {
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++, kk += 2 ) {
      ar = entriesA[kk] ;
      ai = entriesA[kk+1] ;
      Zrecip(ar, ai, &cr, &ci) ;
      br = colB0[kk] ; bi = colB0[kk+1] ;
      colB0[kk]   = br*cr - bi*ci ; 
      colB0[kk+1] = br*ci + bi*cr ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- solve A X = B, where 
     (1) A is block diagonal symmetric
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -----------------------------------
*/
static void
complex_solveBlockDiagonalSym (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ai00, ai01, ai11, ar00, ar01, ar11, bi0, bi1, bi2, 
         br0, br1, br2, ci00, ci01, ci11, cr00, cr01, cr11 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      inc1, inc2, ipivot, irowA, jcolB, kk, m, 
         ncolB, nentA, nrowA, nrowB ;
int      *pivotsizes ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entriesA);
for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
   m = pivotsizes[ipivot] ;
   if ( m != 1 && m != 2 ) {
      fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
              "\n mtxA is block diagonal symmetric"
              "\n pivot %d is %d, not supported",
              mtxA, mtxB, ipivot, m) ;
      exit(-1) ;
   }
   irowA += m ;
}
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + 2*nrowB ;
   colB2 = colB1 + 2*nrowB ;
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         Zrecip(ar00, ai00, &cr00, &ci00) ;
         br0 = colB0[2*irowA] ; bi0 = colB0[2*irowA+1] ;
         colB0[2*irowA]   = br0*cr00 - bi0*ci00 ;
         colB0[2*irowA+1] = br0*ci00 + bi0*cr00 ;
         br1 = colB1[2*irowA] ; bi1 = colB1[2*irowA+1] ;
         colB1[2*irowA]   = br1*cr00 - bi1*ci00 ;
         colB1[2*irowA+1] = br1*ci00 + bi1*cr00 ;
         br2 = colB2[2*irowA] ; bi2 = colB2[2*irowA+1] ;
         colB2[2*irowA]   = br2*cr00 - bi2*ci00 ;
         colB2[2*irowA+1] = br2*ci00 + bi2*cr00 ;
      } else if ( m == 2 ) {
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = entriesA[2*kk+1] ; kk++ ;
         Zrecip2(ar00, ai00, ar01, ai01, ar01, ai01, ar11, ai11,
                 &cr00, &ci00, &cr01, &ci01, NULL, NULL, &cr11, &ci11) ;
         br0 = colB0[2*irowA]   ; bi0 = colB0[2*irowA+1] ;
         br1 = colB0[2*irowA+2] ; bi1 = colB0[2*irowA+3] ;
         colB0[2*irowA]   = cr00*br0 - ci00*bi0 + cr01*br1 - ci01*bi1 ;
         colB0[2*irowA+1] = cr00*bi0 + ci00*br0 + cr01*bi1 + ci01*br1 ;
         colB0[2*irowA+2] = cr01*br0 - ci01*bi0 + cr11*br1 - ci11*bi1 ;
         colB0[2*irowA+3] = cr01*bi0 + ci01*br0 + cr11*bi1 + ci11*br1 ;
         br0 = colB1[2*irowA]   ; bi0 = colB1[2*irowA+1] ;
         br1 = colB1[2*irowA+2] ; bi1 = colB1[2*irowA+3] ;
         colB1[2*irowA]   = cr00*br0 - ci00*bi0 + cr01*br1 - ci01*bi1 ;
         colB1[2*irowA+1] = cr00*bi0 + ci00*br0 + cr01*bi1 + ci01*br1 ;
         colB1[2*irowA+2] = cr01*br0 - ci01*bi0 + cr11*br1 - ci11*bi1 ;
         colB1[2*irowA+3] = cr01*bi0 + ci01*br0 + cr11*bi1 + ci11*br1 ;
         br0 = colB2[2*irowA]   ; bi0 = colB2[2*irowA+1] ;
         br1 = colB2[2*irowA+2] ; bi1 = colB2[2*irowA+3] ;
         colB2[2*irowA]   = cr00*br0 - ci00*bi0 + cr01*br1 - ci01*bi1 ;
         colB2[2*irowA+1] = cr00*bi0 + ci00*br0 + cr01*bi1 + ci01*br1 ;
         colB2[2*irowA+2] = cr01*br0 - ci01*bi0 + cr11*br1 - ci11*bi1 ;
         colB2[2*irowA+3] = cr01*bi0 + ci01*br0 + cr11*bi1 + ci11*br1 ;
      }
      irowA += m ;
   }
   colB0 = colB2 + 2*nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + 2*nrowB ;
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         Zrecip(ar00, ai00, &cr00, &ci00) ;
         br0 = colB0[2*irowA] ; bi0 = colB0[2*irowA+1] ;
         colB0[2*irowA]   = br0*cr00 - bi0*ci00 ;
         colB0[2*irowA+1] = br0*ci00 + bi0*cr00 ;
         br1 = colB1[2*irowA] ; bi1 = colB1[2*irowA+1] ;
         colB1[2*irowA]   = br1*cr00 - bi1*ci00 ;
         colB1[2*irowA+1] = br1*ci00 + bi1*cr00 ;
      } else if ( m == 2 ) {
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = entriesA[2*kk+1] ; kk++ ;
         Zrecip2(ar00, ai00, ar01, ai01, ar01, ai01, ar11, ai11,
                 &cr00, &ci00, &cr01, &ci01, NULL, NULL, &cr11, &ci11) ;
         br0 = colB0[2*irowA]   ; bi0 = colB0[2*irowA+1] ;
         br1 = colB0[2*irowA+2] ; bi1 = colB0[2*irowA+3] ;
         colB0[2*irowA]   = cr00*br0 - ci00*bi0 + cr01*br1 - ci01*bi1 ;
         colB0[2*irowA+1] = cr00*bi0 + ci00*br0 + cr01*bi1 + ci01*br1 ;
         colB0[2*irowA+2] = cr01*br0 - ci01*bi0 + cr11*br1 - ci11*bi1 ;
         colB0[2*irowA+3] = cr01*bi0 + ci01*br0 + cr11*bi1 + ci11*br1 ;
         br0 = colB1[2*irowA]   ; bi0 = colB1[2*irowA+1] ;
         br1 = colB1[2*irowA+2] ; bi1 = colB1[2*irowA+3] ;
         colB1[2*irowA]   = cr00*br0 - ci00*bi0 + cr01*br1 - ci01*bi1 ;
         colB1[2*irowA+1] = cr00*bi0 + ci00*br0 + cr01*bi1 + ci01*br1 ;
         colB1[2*irowA+2] = cr01*br0 - ci01*bi0 + cr11*br1 - ci11*bi1 ;
         colB1[2*irowA+3] = cr01*bi0 + ci01*br0 + cr11*bi1 + ci11*br1 ;
      }
      irowA += m ;
   }
} else if ( jcolB == ncolB - 1 ) {
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         Zrecip(ar00, ai00, &cr00, &ci00) ;
         br0 = colB0[2*irowA] ; bi0 = colB0[2*irowA+1] ;
         colB0[2*irowA]   = br0*cr00 - bi0*ci00 ;
         colB0[2*irowA+1] = br0*ci00 + bi0*cr00 ;
      } else if ( m == 2 ) {
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = entriesA[2*kk+1] ; kk++ ;
         Zrecip2(ar00, ai00, ar01, ai01, ar01, ai01, ar11, ai11,
                 &cr00, &ci00, &cr01, &ci01, NULL, NULL, &cr11, &ci11) ;
         br0 = colB0[2*irowA]   ; bi0 = colB0[2*irowA+1] ;
         br1 = colB0[2*irowA+2] ; bi1 = colB0[2*irowA+3] ;
         colB0[2*irowA]   = cr00*br0 - ci00*bi0 + cr01*br1 - ci01*bi1 ;
         colB0[2*irowA+1] = cr00*bi0 + ci00*br0 + cr01*bi1 + ci01*br1 ;
         colB0[2*irowA+2] = cr01*br0 - ci01*bi0 + cr11*br1 - ci11*bi1 ;
         colB0[2*irowA+3] = cr01*bi0 + ci01*br0 + cr11*bi1 + ci11*br1 ;
      }
      irowA += m ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- solve A X = B, where 
     (1) A is block diagonal hermitian
     (2) X overwrites B
     (B) B has type SUBMTX_DENSE_COLUMNS

   created -- 98may01, cca
   -----------------------------------
*/
static void
complex_solveBlockDiagonalHerm (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) {
double   ai00, ai01, ai11, ar00, ar01, ar11, bi0, bi1, bi2,
         br0, br1, br2, ci00, ci01, ci11, cr00, cr01, cr11 ;
double   *colB0, *colB1, *colB2, *entriesA, *entriesB ;
int      inc1, inc2, ipivot, irowA, jcolB, kk, m, 
         ncolB, nentA, nrowA, nrowB ;
int      *pivotsizes ;
/*
   ----------------------------------------------------
   extract the pointer and dimensions from two matrices
   ----------------------------------------------------
*/
SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entriesA);
for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
   m = pivotsizes[ipivot] ;
   if ( m != 1 && m != 2 ) {
      fprintf(stderr, "\n fatal error in SubMtx_solve(%p,%p)"
              "\n mtxA is block diagonal hermitian"
              "\n pivot %d is %d, not supported",
              mtxA, mtxB, ipivot, m) ;
      exit(-1) ;
   }
   irowA += m ;
}
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entriesB) ;
colB0 = entriesB ;
for ( jcolB = 0 ; jcolB < ncolB - 2 ; jcolB += 3 ) {
   colB1 = colB0 + 2*nrowB ;
   colB2 = colB1 + 2*nrowB ;
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
/*
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
*/
         ar00 = entriesA[2*kk] ; ai00 = 0.0 ; kk++ ;
         Zrecip(ar00, ai00, &cr00, &ci00) ;
         br0 = colB0[2*irowA] ; bi0 = colB0[2*irowA+1] ;
         colB0[2*irowA]   = br0*cr00 ;
         colB0[2*irowA+1] = bi0*cr00 ;
         br1 = colB1[2*irowA] ; bi1 = colB1[2*irowA+1] ;
         colB1[2*irowA]   = br1*cr00 ;
         colB1[2*irowA+1] = bi1*cr00 ;
         br2 = colB2[2*irowA] ; bi2 = colB2[2*irowA+1] ;
         colB2[2*irowA]   = br2*cr00 ;
         colB2[2*irowA+1] = bi2*cr00 ;
      } else if ( m == 2 ) {
/*
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = entriesA[2*kk+1] ; kk++ ;
*/
         ar00 = entriesA[2*kk] ; ai00 = 0.0 ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = 0.0 ; kk++ ;
         Zrecip2(ar00, ai00, ar01, ai01, ar01, -ai01, ar11, ai11,
               &cr00, &ci00, &cr01, &ci01, NULL, NULL, &cr11, &ci11) ;
         br0 = colB0[2*irowA]   ; bi0 = colB0[2*irowA+1] ;
         br1 = colB0[2*irowA+2] ; bi1 = colB0[2*irowA+3] ;
         colB0[2*irowA]   = cr00*br0 + cr01*br1 - ci01*bi1 ;
         colB0[2*irowA+1] = cr00*bi0 + cr01*bi1 + ci01*br1 ;
         colB0[2*irowA+2] = cr01*br0 + ci01*bi0 + cr11*br1 ;
         colB0[2*irowA+3] = cr01*bi0 - ci01*br0 + cr11*bi1 ;
         br0 = colB1[2*irowA]   ; bi0 = colB1[2*irowA+1] ;
         br1 = colB1[2*irowA+2] ; bi1 = colB1[2*irowA+3] ;
         colB1[2*irowA]   = cr00*br0 + cr01*br1 - ci01*bi1 ;
         colB1[2*irowA+1] = cr00*bi0 + cr01*bi1 + ci01*br1 ;
         colB1[2*irowA+2] = cr01*br0 + ci01*bi0 + cr11*br1 ;
         colB1[2*irowA+3] = cr01*bi0 - ci01*br0 + cr11*bi1 ;
         br0 = colB2[2*irowA]   ; bi0 = colB2[2*irowA+1] ;
         br1 = colB2[2*irowA+2] ; bi1 = colB2[2*irowA+3] ;
         colB2[2*irowA]   = cr00*br0 + cr01*br1 - ci01*bi1 ;
         colB2[2*irowA+1] = cr00*bi0 + cr01*bi1 + ci01*br1 ;
         colB2[2*irowA+2] = cr01*br0 + ci01*bi0 + cr11*br1 ;
         colB2[2*irowA+3] = cr01*bi0 - ci01*br0 + cr11*bi1 ;
      }
      irowA += m ;
   }
   colB0 = colB2 + 2*nrowB ;
}
if ( jcolB == ncolB - 2 ) {
   colB1 = colB0 + 2*nrowB ;
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
/*
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
*/
         ar00 = entriesA[2*kk] ; ai00 = 0.0 ; kk++ ;
         Zrecip(ar00, ai00, &cr00, &ci00) ;
         br0 = colB0[2*irowA] ; bi0 = colB0[2*irowA+1] ;
         colB0[2*irowA]   = br0*cr00 ;
         colB0[2*irowA+1] = bi0*cr00 ;
         br1 = colB1[2*irowA] ; bi1 = colB1[2*irowA+1] ;
         colB1[2*irowA]   = br1*cr00 ;
         colB1[2*irowA+1] = bi1*cr00 ;
      } else if ( m == 2 ) {
/*
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = entriesA[2*kk+1] ; kk++ ;
*/
         ar00 = entriesA[2*kk] ; ai00 = 0.0 ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = 0.0 ; kk++ ;
         Zrecip2(ar00, ai00, ar01, ai01, ar01, -ai01, ar11, ai11,
               &cr00, &ci00, &cr01, &ci01, NULL, NULL, &cr11, &ci11) ;
         br0 = colB0[2*irowA]   ; bi0 = colB0[2*irowA+1] ;
         br1 = colB0[2*irowA+2] ; bi1 = colB0[2*irowA+3] ;
         colB0[2*irowA]   = cr00*br0 + cr01*br1 - ci01*bi1 ;
         colB0[2*irowA+1] = cr00*bi0 + cr01*bi1 + ci01*br1 ;
         colB0[2*irowA+2] = cr01*br0 + ci01*bi0 + cr11*br1 ;
         colB0[2*irowA+3] = cr01*bi0 - ci01*br0 + cr11*bi1 ;
         br0 = colB1[2*irowA]   ; bi0 = colB1[2*irowA+1] ;
         br1 = colB1[2*irowA+2] ; bi1 = colB1[2*irowA+3] ;
         colB1[2*irowA]   = cr00*br0 + cr01*br1 - ci01*bi1 ;
         colB1[2*irowA+1] = cr00*bi0 + cr01*bi1 + ci01*br1 ;
         colB1[2*irowA+2] = cr01*br0 + ci01*bi0 + cr11*br1 ;
         colB1[2*irowA+3] = cr01*bi0 - ci01*br0 + cr11*bi1 ;
      }
      irowA += m ;
   }
} else if ( jcolB == ncolB - 1 ) {
   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      if ( m == 1 ) {
/*
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
*/
         ar00 = entriesA[2*kk] ; ai00 = 0.0 ; kk++ ;
         Zrecip(ar00, ai00, &cr00, &ci00) ;
         br0 = colB0[2*irowA] ; bi0 = colB0[2*irowA+1] ;
         colB0[2*irowA]   = br0*cr00 ;
         colB0[2*irowA+1] = bi0*cr00 ;
      } else if ( m == 2 ) {
/*
         ar00 = entriesA[2*kk] ; ai00 = entriesA[2*kk+1] ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = entriesA[2*kk+1] ; kk++ ;
*/
         ar00 = entriesA[2*kk] ; ai00 = 0.0 ; kk++ ;
         ar01 = entriesA[2*kk] ; ai01 = entriesA[2*kk+1] ; kk++ ;
         ar11 = entriesA[2*kk] ; ai11 = 0.0 ; kk++ ;
         Zrecip2(ar00, ai00, ar01, ai01, ar01, -ai01, ar11, ai11,
               &cr00, &ci00, &cr01, &ci01, NULL, NULL, &cr11, &ci11) ;
         br0 = colB0[2*irowA]   ; bi0 = colB0[2*irowA+1] ;
         br1 = colB0[2*irowA+2] ; bi1 = colB0[2*irowA+3] ;
         colB0[2*irowA]   = cr00*br0 + cr01*br1 - ci01*bi1 ;
         colB0[2*irowA+1] = cr00*bi0 + cr01*bi1 + ci01*br1 ;
         colB0[2*irowA+2] = cr01*br0 + ci01*bi0 + cr11*br1 ;
         colB0[2*irowA+3] = cr01*bi0 - ci01*br0 + cr11*bi1 ;
      }
      irowA += m ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
