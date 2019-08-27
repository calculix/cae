/*  extract.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- to extract a submatrix, B = A(BrowsIV, BcolsIV)

   return values ---
      1 -- normal return
     -1 -- B is NULL
     -2 -- BcolsIV is NULL
     -3 -- BrowsIV is NULL
     -4 -- B is NULL
     -5 -- invalid input mode for A
     -6 -- invalid coordinate type for A
     -7 -- invalid symmetryflag
     -8 -- hermitian flag but not complex
     -9 -- msglvl > 0 and msgFile = NULL

   created -- 98oct15, cca
   ----------------------------------------------------------
*/
int
InpMtx_initFromSubmatrix (
   InpMtx   *B,
   InpMtx   *A,
   IV       *BrowsIV,
   IV       *BcolsIV,
   int      symmetryflag,
   int      msglvl,
   FILE     *msgFile
) {
double   *dbuf, *entA ;
int      colA, colB, ii, jj, jjfirst, jjlast, kk, maxcolA, maxn,
         maxrowA, nbuf, ncolB, nent, nrowB, nvector, rowA, rowB, 
         oldCoordType, oldStorageMode, rowsAndColumnsAreIdentical ;
int      *colmap, *colsA, *colsB, *ibuf1, *ibuf2, *offsets, *rowmap, 
         *rowsB, *sizes, *vecids ;
/*
   ---------------
   check the input
   ---------------
*/
if ( B == NULL ) {
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n B is NULL\n") ;
   return(-1) ;
}
if ( BrowsIV == NULL ) {
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n BrowsIV is NULL\n") ;
   return(-2) ;
}
if ( BcolsIV == NULL ) {
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n BcolsIV is NULL\n") ;
   return(-3) ;
}
if ( A == NULL ) {
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n A is NULL\n") ;
   return(-4) ;
}
switch ( A->inputMode ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
case INPMTX_INDICES_ONLY :
   break ;
default :
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n invalid inputMode %d for A\n", A->inputMode) ;
   return(-5) ;
   break ;
}
switch ( A->coordType ) {
case INPMTX_BY_ROWS     :
case INPMTX_BY_COLUMNS  :
case INPMTX_BY_CHEVRONS :
   break ;
default :
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n invalid coordType %d for A\n", A->coordType) ;
   return(-6) ;
   break ;
}
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
case SPOOLES_NONSYMMETRIC :
   break ;
default :
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n invalid symmetryflag %d \n", symmetryflag) ;
   return(-7) ;
   break ;
}
if (  A->inputMode != SPOOLES_COMPLEX 
   && symmetryflag == SPOOLES_HERMITIAN ) {
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n Hermitian symmetry flag but A is not complex\n") ;
   return(-8) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n error in InpMtx_initFromSubmatrix()"
           "\n msglvl = %d, msgFile = NULL\n", msglvl) ;
   return(-9) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------
   initialize B
   ------------
*/
InpMtx_init(B, INPMTX_BY_ROWS, A->inputMode, 0, 0) ;
/*
   -----------------------------
   get the range of entries in A
   -----------------------------
*/
InpMtx_range(A, NULL, &maxcolA, NULL, &maxrowA) ;
maxn = (maxcolA >= maxrowA) ? maxcolA : maxrowA ;
/*
   --------------------------------------------------
   get the # and indices of the rows and columns of B
   --------------------------------------------------
*/
IV_sizeAndEntries(BrowsIV, &nrowB, &rowsB) ;
IV_sizeAndEntries(BcolsIV, &ncolB, &colsB) ;
if ( nrowB != ncolB ) {
   rowsAndColumnsAreIdentical = 0 ;
} else {
   rowsAndColumnsAreIdentical = 1 ;
   for ( ii = 0 ; ii < nrowB ; ii++ ) {
      if ( rowsB[ii] != colsB[ii] ) {
         rowsAndColumnsAreIdentical = 0 ;
         break ;
      }
   }
}
/*
   ---------------------------
   set up the local column map
   ---------------------------
*/
colmap = IVinit(1+maxn, -1) ;
for ( ii = 0 ; ii < ncolB ; ii++ ) {
   colmap[colsB[ii]] = ii ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n colmap") ;
   IVfprintf(msgFile, 1+maxn, colmap) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   change the coordinate type of A to rows
   and storage mode to vectors
   ---------------------------------------
*/
if ( (oldCoordType = A->coordType) != INPMTX_BY_ROWS ) {
   InpMtx_changeCoordType(A, INPMTX_BY_ROWS) ;
}
if ( (oldStorageMode = A->storageMode) != INPMTX_BY_VECTORS ) {
   InpMtx_changeStorageMode(A, INPMTX_BY_VECTORS) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n A") ;
   InpMtx_writeForHumanEye(A, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   switch over the input mode. search the rows of A that will
   belong to B, load into the buffer all entries that belong 
   in columns of B.  when the buffer is full, add to B matrix.
   -----------------------------------------------------------
*/
nbuf  = 100 ;
ibuf1 = IVinit(nbuf, -1)  ;
ibuf2 = IVinit(nbuf, -1)  ;
kk    = 0 ;
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   dbuf = DVinit(nbuf, 0.0) ;
   for ( rowB = 0 ; rowB < nrowB ; rowB++ ) {
      rowA = rowsB[rowB] ;
      InpMtx_realVector(A, rowA, &nent, &colsA, &entA) ;
      for ( jj = 0 ; jj < nent ; jj++ ) {
         if ( (colB = colmap[colsA[jj]]) != -1 ) {
            ibuf1[kk] = rowB ; ibuf2[kk] = colB ;
            dbuf[kk] = entA[jj] ;
            if ( ++kk == nbuf ) {
               InpMtx_inputRealTriples(B, kk, ibuf1, ibuf2, dbuf) ;
               kk = 0 ;
            }
         }
      }
   } 
} else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   dbuf = DVinit(2*nbuf, 0.0) ;
   for ( rowB = 0 ; rowB < nrowB ; rowB++ ) {
      rowA = rowsB[rowB] ;
      InpMtx_complexVector(A, rowA, &nent, &colsA, &entA) ;
      for ( jj = 0 ; jj < nent ; jj++ ) {
         if ( (colB = colmap[colsA[jj]]) != -1 ) {
            ibuf1[kk] = rowB ; ibuf2[kk] = colB ;
            dbuf[2*kk] = entA[2*jj] ; dbuf[2*kk+1] = entA[2*jj+1] ;
            if ( ++kk == nbuf ) {
               InpMtx_inputComplexTriples(B, kk, ibuf1, ibuf2, dbuf) ;
               kk = 0 ;
            }
         }
      }
   }
} else if ( INPMTX_IS_INDICES_ONLY(A) ) {
   dbuf = NULL ;
   for ( rowB = 0 ; rowB < nrowB ; rowB++ ) {
      rowA = rowsB[rowB] ;
      InpMtx_vector(A, rowA, &nent, &colsA) ;
      for ( jj = 0 ; jj < nent ; jj++ ) {
         if ( (colB = colmap[colsA[jj]]) != -1 ) {
            ibuf1[kk] = rowB ; ibuf2[kk] = colB ;
            if ( ++kk == nbuf ) {
               InpMtx_inputTriples(B, kk, ibuf1, ibuf2) ;
               kk = 0 ;
            }
         }
      }
   }
}
if (  (   symmetryflag == SPOOLES_SYMMETRIC
       || symmetryflag == SPOOLES_HERMITIAN) 
   && !rowsAndColumnsAreIdentical ) {
/*
   ----------------------------------------------
   matrix is symmetric or hermitian, search lower 
   triangle of A for entries that belong to B.
   ----------------------------------------------
*/
   rowmap = IVinit(1+maxn, -1) ;
   for ( ii = 0 ; ii < nrowB ; ii++ ) {
      rowmap[rowsB[ii]] = ii ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n rowmap") ;
      IVfprintf(msgFile, 1+maxn, rowmap) ;
      fflush(msgFile) ;
   }
   if ( INPMTX_IS_REAL_ENTRIES(A) ) {
      for ( colB = 0 ; colB < ncolB ; colB++ ) {
         colA = colsB[colB] ;
         InpMtx_realVector(A, colA, &nent, &colsA, &entA) ;
         for ( jj = 0 ; jj < nent ; jj++ ) {
            if ( (rowB = rowmap[colsA[jj]]) != -1 ) {
               ibuf1[kk] = rowB ; ibuf2[kk] = colB ;
               dbuf[kk] = entA[jj] ;
               if ( ++kk == nbuf ) {
                  InpMtx_inputRealTriples(B, kk, ibuf1, ibuf2, dbuf) ;
                  kk = 0 ;
               }
            }
         }
      } 
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
      for ( colB = 0 ; colB < ncolB ; colB++ ) {
         colA = colsB[colB] ;
         InpMtx_complexVector(A, colA, &nent, &colsA, &entA) ;
         for ( jj = 0 ; jj < nent ; jj++ ) {
            if ( (rowB = rowmap[colsA[jj]]) != -1 ) {
               ibuf1[kk] = rowB ; ibuf2[kk] = colB ;
               dbuf[2*kk] = entA[2*jj] ; dbuf[2*kk+1] = -entA[2*jj+1] ;
               if ( ++kk == nbuf ) {
                  InpMtx_inputComplexTriples(B, kk, ibuf1, ibuf2, dbuf);
                  kk = 0 ;
               }
            }
         }
      }
   } else if ( INPMTX_IS_INDICES_ONLY(A) ) {
      for ( colB = 0 ; colB < ncolB ; colB++ ) {
         colA = colsB[colB] ;
         InpMtx_vector(A, colA, &nent, &colsA) ;
         for ( jj = 0 ; jj < nent ; jj++ ) {
            if ( (rowB = rowmap[colsA[jj]]) != -1 ) {
               ibuf1[kk] = rowB ; ibuf2[kk] = colB ;
               if ( ++kk == nbuf ) {
                  InpMtx_inputTriples(B, kk, ibuf1, ibuf2) ;
                  kk = 0 ;
               }
            }
         }
      }
   }
   IVfree(rowmap) ;
}
if ( kk > 0 ) {
/*
   -------------------------------------------
   load remaining entries in the buffer into B
   -------------------------------------------
*/
   if ( INPMTX_IS_REAL_ENTRIES(A) ) {
      InpMtx_inputRealTriples(B, kk, ibuf1, ibuf2, dbuf) ;
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
      InpMtx_inputComplexTriples(B, kk, ibuf1, ibuf2, dbuf) ;
   } else {
      InpMtx_inputTriples(B, kk, ibuf1, ibuf2) ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n B") ;
   InpMtx_writeForHumanEye(B, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------
   set B's coordinate type and storage 
   mode to be the same as A's on input
   -----------------------------------
*/
InpMtx_changeCoordType(B, oldCoordType) ;
InpMtx_changeStorageMode(B, oldStorageMode) ;
/*
   -------------------------------------------
   change back to the old coordinate type of A
   -------------------------------------------
*/
if ( (oldCoordType = A->coordType) != INPMTX_BY_ROWS ) {
   InpMtx_changeCoordType(A, oldCoordType) ;
}
if ( oldStorageMode != A->storageMode ) {
   InpMtx_changeStorageMode(A, oldStorageMode) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(colmap) ;
IVfree(ibuf1) ;
IVfree(ibuf2) ;
if ( dbuf != NULL ) {
   DVfree(dbuf) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
