/*  solveupd.c  */

#include "../SubMtx.h"

/*--------------------------------------------------------------------*/
static void 
real_updDenseColumns  ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
static void 
real_updDenseRows     ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
static void 
real_updSparseRows    ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
static void 
real_updSparseColumns ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
static void 
complex_updDenseColumns  ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
static void 
complex_updDenseRows     ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
static void 
complex_updSparseRows    ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
static void 
complex_updSparseColumns ( SubMtx *mtxY, SubMtx *mtxA, SubMtx *mtxX ) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- perform the matrix-matrix multiply 
      Y := Y - A * X used in the forward and backsolves
      where
        (1) rows(A) \subseteq rows(Y)
        (2) rows(A) are local w.r.t. rows(Y)
        (3) cols(A) \subseteq rows(X)
        (4) cols(A) are local w.r.t. rows(X)
        (5) cols(Y) = cols(X)
        (6) Y and X have mode SUBMTX_DENSE_COLUMNS
 
   created -- 98may02, cca
   ----------------------------------------------------
*/
void
SubMtx_solveupd (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtxY == NULL || mtxA == NULL || mtxX == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveupd(%p,%p,%p)"
           "\n bad input\n", mtxY, mtxA, mtxX) ;
   exit(-1) ;
}
if ( mtxY->mode != SUBMTX_DENSE_COLUMNS ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveupd(%p,%p,%p)"
           "\n Y must have mode SUBMTX_DENSE_COLUMNS\n", 
           mtxY, mtxA, mtxX) ;
   exit(-1) ;
}
if ( mtxX->mode != SUBMTX_DENSE_COLUMNS ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveupd(%p,%p,%p)"
           "\n X must have mode SUBMTX_DENSE_COLUMNS\n", 
           mtxY, mtxA, mtxX) ;
   exit(-1) ;
}
switch ( mtxA->type ) {
case SPOOLES_REAL :
   switch ( mtxA->mode ) {
   case SUBMTX_DENSE_COLUMNS :
      real_updDenseColumns(mtxY, mtxA, mtxX) ;
      break ;
   case SUBMTX_DENSE_ROWS :
      real_updDenseRows(mtxY, mtxA, mtxX) ;
      break ;
   case SUBMTX_SPARSE_ROWS :
      real_updSparseRows(mtxY, mtxA, mtxX) ;
      break ;
   case SUBMTX_SPARSE_COLUMNS :
      real_updSparseColumns(mtxY, mtxA, mtxX) ;
      break ;
   default :
      fprintf(stderr, "\n fatal error in SubMtx_solveupd(%p,%p,%p)"
              "\n unsupported mode %d for A\n",
              mtxY, mtxA, mtxX, mtxA->mode) ;
      exit(-1) ;
      break ;
   }
   break ;
case SPOOLES_COMPLEX :
   switch ( mtxA->mode ) {
   case SUBMTX_DENSE_COLUMNS :
      complex_updDenseColumns(mtxY, mtxA, mtxX) ;
      break ;
   case SUBMTX_DENSE_ROWS :
      complex_updDenseRows(mtxY, mtxA, mtxX) ;
      break ;
   case SUBMTX_SPARSE_ROWS :
      complex_updSparseRows(mtxY, mtxA, mtxX) ;
      break ;
   case SUBMTX_SPARSE_COLUMNS :
      complex_updSparseColumns(mtxY, mtxA, mtxX) ;
      break ;
   default :
      fprintf(stderr, "\n fatal error in SubMtx_solveupd(%p,%p,%p)"
              "\n unsupported mode %d for A\n",
              mtxY, mtxA, mtxX, mtxA->mode) ;
      SubMtx_writeForHumanEye(mtxA, stderr) ;
      exit(-1) ;
      break ;
   }
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_solveupd(%p,%p,%p)"
           "\n unsupported type %d for A\n",
           mtxY, mtxA, mtxX, mtxA->type) ;
   exit(-1) ;
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------
   A has dense columns
   -------------------
*/
static void
real_updDenseColumns (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   Ak0, Ak1, Ak2, x00, x01, x02, x10, x11, x12,
         x20, x21, x22 ;
double   *colA0, *colA1, *colA2, *colX0, *colX1, *colX2, 
         *colY0, *colY1, *colY2, *entA, *entX, *entY ;
int      icolA, inc1, inc2, irowX, jcolX, krowA, krowY, 
         ncolA, ncolX, ncolY, nrowA, nrowX, nrowY ;
int      *colindA, *rowindA ;

SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &nrowA, &ncolA, &inc1, &inc2, &entA) ;
colX0 = entX ;
colY0 = entY ;
if ( ncolA != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( nrowA != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + nrowX ;
   colX2 = colX1 + nrowX ;
   colY1 = colY0 + nrowY ;
   colY2 = colY1 + nrowY ;
   colA0 = entA ;
   for ( icolA = 0 ; icolA < ncolA - 2 ; icolA += 3 ) {
      colA1 = colA0 + nrowA ;
      colA2 = colA1 + nrowA ;
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x01 = colX1[icolA] ;
         x02 = colX2[icolA] ;
         x10 = colX0[icolA+1] ;
         x11 = colX1[icolA+1] ;
         x12 = colX2[icolA+1] ;
         x20 = colX0[icolA+2] ;
         x21 = colX1[icolA+2] ;
         x22 = colX2[icolA+2] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         x02 = colX2[irowX] ;
         irowX = colindA[icolA+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
         x12 = colX2[irowX] ;
         irowX = colindA[icolA+2] ;
         x20 = colX0[irowX] ;
         x21 = colX1[irowX] ;
         x22 = colX2[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            Ak2 = colA2[krowA] ;
            colY0[krowA] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowA] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
            colY2[krowA] -= Ak0 * x02 + Ak1 * x12 + Ak2 * x22 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            Ak2 = colA2[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
            colY2[krowY] -= Ak0 * x02 + Ak1 * x12 + Ak2 * x22 ;
         }
      }
      colA0 = colA2 + nrowA ;
   }
   if ( icolA == ncolA - 2 ) {
      colA1 = colA0 + nrowA ;
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x01 = colX1[icolA] ;
         x02 = colX2[icolA] ;
         x10 = colX0[icolA+1] ;
         x11 = colX1[icolA+1] ;
         x12 = colX2[icolA+1] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         x02 = colX2[irowX] ;
         irowX = colindA[icolA+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
         x12 = colX2[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            colY0[krowA] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowA] -= Ak0 * x01 + Ak1 * x11 ;
            colY2[krowA] -= Ak0 * x02 + Ak1 * x12 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 ;
            colY2[krowY] -= Ak0 * x02 + Ak1 * x12 ;
         }
      }
   } else if ( icolA == ncolA - 1 ) {
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x01 = colX1[icolA] ;
         x02 = colX2[icolA] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         x02 = colX2[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            colY0[krowA] -= Ak0 * x00 ;
            colY1[krowA] -= Ak0 * x01 ;
            colY2[krowA] -= Ak0 * x02 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 ;
            colY1[krowY] -= Ak0 * x01 ;
            colY2[krowY] -= Ak0 * x02 ;
         }
      }
   }
   colX0 = colX2 + nrowX ;
   colY0 = colY2 + nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + nrowX ;
   colY1 = colY0 + nrowY ;
   colA0 = entA ;
   for ( icolA = 0 ; icolA < ncolA - 2 ; icolA += 3 ) {
      colA1 = colA0 + nrowA ;
      colA2 = colA1 + nrowA ;
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x01 = colX1[icolA] ;
         x10 = colX0[icolA+1] ;
         x11 = colX1[icolA+1] ;
         x20 = colX0[icolA+2] ;
         x21 = colX1[icolA+2] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         irowX = colindA[icolA+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
         irowX = colindA[icolA+2] ;
         x20 = colX0[irowX] ;
         x21 = colX1[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            Ak2 = colA2[krowA] ;
            colY0[krowA] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowA] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            Ak2 = colA2[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
         }
      }
      colA0 = colA2 + nrowA ;
   }
   if ( icolA == ncolA - 2 ) {
      colA1 = colA0 + nrowA ;
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x01 = colX1[icolA] ;
         x10 = colX0[icolA+1] ;
         x11 = colX1[icolA+1] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         irowX = colindA[icolA+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            colY0[krowA] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowA] -= Ak0 * x01 + Ak1 * x11 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 ;
         }
      }
   } else if ( icolA == ncolA - 1 ) {
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x01 = colX1[icolA] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            colY0[krowA] -= Ak0 * x00 ;
            colY1[krowA] -= Ak0 * x01 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 ;
            colY1[krowY] -= Ak0 * x01 ;
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   colA0 = entA ;
   for ( icolA = 0 ; icolA < ncolA - 2 ; icolA += 3 ) {
      colA1 = colA0 + nrowA ;
      colA2 = colA1 + nrowA ;
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x10 = colX0[icolA+1] ;
         x20 = colX0[icolA+2] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         irowX = colindA[icolA+1] ;
         x10 = colX0[irowX] ;
         irowX = colindA[icolA+2] ;
         x20 = colX0[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            Ak2 = colA2[krowA] ;
            colY0[krowA] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            Ak2 = colA2[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
         }
      }
      colA0 = colA2 + nrowA ;
   }
   if ( icolA == ncolA - 2 ) {
      colA1 = colA0 + nrowA ;
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
         x10 = colX0[icolA+1] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
         irowX = colindA[icolA+1] ;
         x10 = colX0[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            colY0[krowA] -= Ak0 * x00 + Ak1 * x10 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            Ak1 = colA1[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 ;
         }
      }
   } else if ( icolA == ncolA - 1 ) {
      if ( ncolA == nrowX ) {
         x00 = colX0[icolA] ;
      } else {
         irowX = colindA[icolA] ;
         x00 = colX0[irowX] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            colY0[krowA] -= Ak0 * x00 ;
         }
      } else {
         for ( krowA = 0 ; krowA < nrowA ; krowA++ ) {
            Ak0 = colA0[krowA] ;
            krowY = rowindA[krowA] ;
            colY0[krowY] -= Ak0 * x00 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------
   A has dense rows
   ----------------
*/
static void
real_updDenseRows (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2, 
         *rowA0, *rowA1, *rowA2, *entA, *entX, *entY ;
int      inc1, inc2, irowA, irowY, jcolX, kcolA, krowX, 
         ncolA, ncolX, ncolY, nrowA, nrowX, nrowY ;
int      *colindA, *rowindA ;

SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &nrowA, &ncolA, &inc1, &inc2, &entA) ;
if ( ncolA != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( nrowA != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + nrowX ;
   colX2 = colX1 + nrowX ;
   colY1 = colY0 + nrowY ;
   colY2 = colY1 + nrowY ;
   rowA0 = entA ;
   for ( irowA = 0 ; irowA < nrowA - 2 ; irowA += 3 ) {
      double   A0k, A1k, A2k, Xk0, Xk1, Xk2 ;
      double   sum00, sum01, sum02, 
               sum10, sum11, sum12, 
               sum20, sum21, sum22 ;

      sum00 = sum01 = sum02 
            = sum10 = sum11 = sum12 = sum20 = sum21 = sum22 = 0.0 ;
      rowA1 = rowA0 + ncolA ;
      rowA2 = rowA1 + ncolA ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            A2k = rowA2[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            Xk1 = colX1[kcolA] ;
            Xk2 = colX2[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum12 += A1k * Xk2 ;
            sum20 += A2k * Xk0 ;
            sum21 += A2k * Xk1 ;
            sum22 += A2k * Xk2 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            A2k = rowA2[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            Xk2 = colX2[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum12 += A1k * Xk2 ;
            sum20 += A2k * Xk0 ;
            sum21 += A2k * Xk1 ;
            sum22 += A2k * Xk2 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY1[irowA]   -= sum01 ;
         colY2[irowA]   -= sum02 ;
         colY0[irowA+1] -= sum10 ;
         colY1[irowA+1] -= sum11 ;
         colY2[irowA+1] -= sum12 ;
         colY0[irowA+2] -= sum20 ;
         colY1[irowA+2] -= sum21 ;
         colY2[irowA+2] -= sum22 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         colY2[irowY] -= sum02 ;
         irowY = rowindA[irowA+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
         colY2[irowY] -= sum12 ;
         irowY = rowindA[irowA+2] ;
         colY0[irowY] -= sum20 ;
         colY1[irowY] -= sum21 ;
         colY2[irowY] -= sum22 ;
      }
      rowA0 = rowA2 + ncolA ;
   }
   if ( irowA == nrowA - 2 ) {
      double   A0k, A1k, Xk0, Xk1, Xk2 ;
      double   sum00, sum01, sum02, sum10, sum11, sum12 ;

      sum00 = sum01 = sum02 = sum10 = sum11 = sum12 = 0.0 ;
      rowA1 = rowA0 + ncolA ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            Xk1 = colX1[kcolA] ;
            Xk2 = colX2[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum12 += A1k * Xk2 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            Xk2 = colX2[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum12 += A1k * Xk2 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY1[irowA]   -= sum01 ;
         colY2[irowA]   -= sum02 ;
         colY0[irowA+1] -= sum10 ;
         colY1[irowA+1] -= sum11 ;
         colY2[irowA+1] -= sum12 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         colY2[irowY] -= sum02 ;
         irowY = rowindA[irowA+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
         colY2[irowY] -= sum12 ;
      }
   } else if ( irowA == nrowA - 1 ) {
      double   A0k, Xk0, Xk1, Xk2 ;
      double   sum00, sum01, sum02 ;

      sum00 = sum01 = sum02 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            Xk1 = colX1[kcolA] ;
            Xk2 = colX2[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            Xk2 = colX2[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY1[irowA]   -= sum01 ;
         colY2[irowA]   -= sum02 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         colY2[irowY] -= sum02 ;
      }
   }
   colX0 = colX2 + nrowX ;
   colY0 = colY2 + nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + nrowX ;
   colY1 = colY0 + nrowY ;
   rowA0 = entA ;
   for ( irowA = 0 ; irowA < nrowA - 2 ; irowA += 3 ) {
      double   A0k, A1k, A2k, Xk0, Xk1 ;
      double   sum00, sum01, sum10, sum11, sum20, sum21 ;

      sum00 = sum01 = sum10 = sum11 = sum20 = sum21 = 0.0 ;
      rowA1 = rowA0 + ncolA ;
      rowA2 = rowA1 + ncolA ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            A2k = rowA2[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            Xk1 = colX1[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum20 += A2k * Xk0 ;
            sum21 += A2k * Xk1 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            A2k = rowA2[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum20 += A2k * Xk0 ;
            sum21 += A2k * Xk1 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY1[irowA]   -= sum01 ;
         colY0[irowA+1] -= sum10 ;
         colY1[irowA+1] -= sum11 ;
         colY0[irowA+2] -= sum20 ;
         colY1[irowA+2] -= sum21 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         irowY = rowindA[irowA+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
         irowY = rowindA[irowA+2] ;
         colY0[irowY] -= sum20 ;
         colY1[irowY] -= sum21 ;
      }
      rowA0 = rowA2 + ncolA ;
   }
   if ( irowA == nrowA - 2 ) {
      double   A0k, A1k, Xk0, Xk1 ;
      double   sum00, sum01, sum10, sum11 ;

      sum00 = sum01 = sum10 = sum11 = 0.0 ;
      rowA1 = rowA0 + ncolA ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            Xk1 = colX1[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY1[irowA]   -= sum01 ;
         colY0[irowA+1] -= sum10 ;
         colY1[irowA+1] -= sum11 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         irowY = rowindA[irowA+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
      }
   } else if ( irowA == nrowA - 1 ) {
      double   A0k, Xk0, Xk1 ;
      double   sum00, sum01 ;

      sum00 = sum01 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            Xk1 = colX1[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY1[irowA]   -= sum01 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   rowA0 = entA ;
   for ( irowA = 0 ; irowA < nrowA - 2 ; irowA += 3 ) {
      double   A0k, A1k, A2k, Xk0 ;
      double   sum00, sum10, sum20 ;

      sum00 = sum10 = sum20 = 0.0 ;
      rowA1 = rowA0 + ncolA ;
      rowA2 = rowA1 + ncolA ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            A2k = rowA2[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
            sum20 += A2k * Xk0 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            A2k = rowA2[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
            sum20 += A2k * Xk0 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY0[irowA+1] -= sum10 ;
         colY0[irowA+2] -= sum20 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         irowY = rowindA[irowA+1] ;
         colY0[irowY] -= sum10 ;
         irowY = rowindA[irowA+2] ;
         colY0[irowY] -= sum20 ;
      }
      rowA0 = rowA2 + ncolA ;
   }
   if ( irowA == nrowA - 2 ) {
      double   A0k, A1k, Xk0 ;
      double   sum00, sum10 ;

      sum00 = sum10 = 0.0 ;
      rowA1 = rowA0 + ncolA ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            A1k = rowA1[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
         colY0[irowA+1] -= sum10 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
         irowY = rowindA[irowA+1] ;
         colY0[irowY] -= sum10 ;
      }
   } else if ( irowA == nrowA - 1 ) {
      double   A0k, Xk0 ;
      double   sum00 ;

      sum00 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            Xk0 = colX0[kcolA] ;
            sum00 += A0k * Xk0 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            A0k = rowA0[kcolA] ; 
            krowX = colindA[kcolA] ;
            Xk0 = colX0[krowX] ; 
            sum00 += A0k * Xk0 ;
         }
      }
      if ( nrowY == nrowA ) {
         colY0[irowA]   -= sum00 ;
      } else {
         irowY = rowindA[irowA] ;
         colY0[irowY] -= sum00 ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------
   A has sparse rows
   -----------------
*/
static void
real_updSparseRows (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   Aik, sum0, sum1, sum2 ;
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, inc1, inc2, irowA, irowY, jcolX, kk, krowX, 
         ncolA, ncolX, ncolY, nentA, nrowA, nrowX, nrowY, size ;
int      *colindA, *indices, *rowindA, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_ROWS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseRowsInfo(mtxA, &nrowA, &nentA, &sizes, &indices, &entA) ;
if ( (ncolA = mtxA->ncol) != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( nrowA != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + nrowX ;
   colX2 = colX1 + nrowX ;
   colY1 = colY0 + nrowY ;
   colY2 = colY1 + nrowY ;
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizes[irowA]) > 0 ) {
         sum0 = sum1 = sum2 = 0.0 ;
         if ( ncolA == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = indices[kk] ;
               sum0 += Aik * colX0[krowX] ;
               sum1 += Aik * colX1[krowX] ;
               sum2 += Aik * colX2[krowX] ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = colindA[indices[kk]] ;
               sum0 += Aik * colX0[krowX] ;
               sum1 += Aik * colX1[krowX] ;
               sum2 += Aik * colX2[krowX] ;
            }
         }
         if ( nrowA == nrowY ) {
            colY0[irowA] -= sum0 ;
            colY1[irowA] -= sum1 ;
            colY2[irowA] -= sum2 ;
         } else {
            irowY = rowindA[irowA] ;
            colY0[irowY] -= sum0 ;
            colY1[irowY] -= sum1 ;
            colY2[irowY] -= sum2 ;
         }
      }
   }
   colX0 = colX2 + nrowX ;
   colY0 = colY2 + nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + nrowX ;
   colY1 = colY0 + nrowY ;
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizes[irowA]) > 0 ) {
         sum0 = sum1 = 0.0 ;
         if ( ncolA == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = indices[kk] ;
               sum0 += Aik * colX0[krowX] ;
               sum1 += Aik * colX1[krowX] ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = colindA[indices[kk]] ;
               sum0 += Aik * colX0[krowX] ;
               sum1 += Aik * colX1[krowX] ;
            }
         }
         if ( nrowA == nrowY ) {
            colY0[irowA] -= sum0 ;
            colY1[irowA] -= sum1 ;
         } else {
            irowY = rowindA[irowA] ;
            colY0[irowY] -= sum0 ;
            colY1[irowY] -= sum1 ;
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizes[irowA]) > 0 ) {
         sum0 = 0.0 ;
         if ( ncolA == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = indices[kk] ;
               sum0 += Aik * colX0[krowX] ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = colindA[indices[kk]] ;
               sum0 += Aik * colX0[krowX] ;
            }
         }
         if ( nrowA == nrowY ) {
            colY0[irowA] -= sum0 ;
         } else {
            irowY = rowindA[irowA] ;
            colY0[irowY] -= sum0 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------
   A has sparse columns
   --------------------
*/
static void
real_updSparseColumns (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   Aij, Xj0, Xj1, Xj2 ;
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, inc1, inc2, irowY, jcolA, jcolX, jrowX, kk, 
         ncolA, ncolX, ncolY, nentA, nrowA, nrowX, nrowY, size ;
int      *colindA, *indices, *rowindA, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_COLUMNS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseColumnsInfo(mtxA, &ncolA, &nentA, &sizes, &indices, &entA) ;
if ( ncolA != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( (nrowA = mtxA->nrow) != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + nrowX ;
   colX2 = colX1 + nrowX ;
   colY1 = colY0 + nrowY ;
   colY2 = colY1 + nrowY ;
   for ( jcolA = kk = 0 ; jcolA < ncolA ; jcolA++ ) {
      if ( (size = sizes[jcolA]) > 0 ) {
         if ( ncolA == nrowX ) {
            jrowX = jcolA ;
         } else {
            jrowX = colindA[jcolA] ;
         }
         Xj0 = colX0[jrowX] ;
         Xj1 = colX1[jrowX] ;
         Xj2 = colX2[jrowX] ;
         if ( nrowA == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = indices[kk] ;
               colY0[irowY] -= Aij * Xj0 ;
               colY1[irowY] -= Aij * Xj1 ;
               colY2[irowY] -= Aij * Xj2 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = rowindA[indices[kk]] ;
               colY0[irowY] -= Aij * Xj0 ;
               colY1[irowY] -= Aij * Xj1 ;
               colY2[irowY] -= Aij * Xj2 ;
            }
         }
      }
   }
   colX0 = colX2 + nrowX ;
   colY0 = colY2 + nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + nrowX ;
   colY1 = colY0 + nrowY ;
   for ( jcolA = kk = 0 ; jcolA < ncolA ; jcolA++ ) {
      if ( (size = sizes[jcolA]) > 0 ) {
         if ( ncolA == nrowX ) {
            jrowX = jcolA ;
         } else {
            jrowX = colindA[jcolA] ;
         }
         Xj0 = colX0[jrowX] ;
         Xj1 = colX1[jrowX] ;
         if ( nrowA == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = indices[kk] ;
               colY0[irowY] -= Aij * Xj0 ;
               colY1[irowY] -= Aij * Xj1 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = rowindA[indices[kk]] ;
               colY0[irowY] -= Aij * Xj0 ;
               colY1[irowY] -= Aij * Xj1 ;
            }
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   for ( jcolA = kk = 0 ; jcolA < ncolA ; jcolA++ ) {
      if ( (size = sizes[jcolA]) > 0 ) {
         if ( ncolA == nrowX ) {
            jrowX = jcolA ;
         } else {
            jrowX = colindA[jcolA] ;
         }
         Xj0 = colX0[jrowX] ;
         if ( nrowA == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = indices[kk] ;
               colY0[irowY] -= Aij * Xj0 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = rowindA[indices[kk]] ;
               colY0[irowY] -= Aij * Xj0 ;
            }
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------
   A has dense columns
   -------------------
*/
static void
complex_updDenseColumns (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   ai0, ai1, ai2, ar0, ar1, ar2,
         xi00, xi01, xi02, xi10, xi11, xi12, xi20, xi21, xi22, 
         xr00, xr01, xr02, xr10, xr11, xr12, xr20, xr21, xr22 ;
double   *colA0, *colA1, *colA2, *colX0, *colX1, *colX2, 
         *colY0, *colY1, *colY2, *entA, *entX, *entY ;
int      icolA, iloc, inc1, inc2, iyloc, jcolX, krowA, krowY, 
         ncolA, ncolX, ncolY, nrowA, nrowX, nrowY, rloc, ryloc ;
int      *colindA, *rowindA ;

SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &nrowA, &ncolA, &inc1, &inc2, &entA) ;
colX0 = entX ;
colY0 = entY ;
if ( ncolA != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( nrowA != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + 2*nrowX ;
   colX2 = colX1 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colY2 = colY1 + 2*nrowY ;
   colA0 = entA ;
   for ( icolA = 0 ; icolA < ncolA - 2 ; icolA += 3 ) {
      colA1 = colA0 + 2*nrowA ;
      colA2 = colA1 + 2*nrowA ;
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc += 2, iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc += 2, iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc = 2*colindA[icolA+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc = 2*colindA[icolA+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            ar2 = colA2[rloc] ; ai2 = colA2[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 
                         + ar1*xr10 - ai1*xi10 
                         + ar2*xr20 - ai2*xi20 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 
                         + ar1*xi10 + ai1*xr10 
                         + ar2*xi20 + ai2*xr20 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 
                         + ar1*xr11 - ai1*xi11 
                         + ar2*xr21 - ai2*xi21 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 
                         + ar1*xi11 + ai1*xr11 
                         + ar2*xi21 + ai2*xr21 ;
            colY2[rloc] -= ar0*xr02 - ai0*xi02 
                         + ar1*xr12 - ai1*xi12 
                         + ar2*xr22 - ai2*xi22 ;
            colY2[iloc] -= ar0*xi02 + ai0*xr02 
                         + ar1*xi12 + ai1*xr12 
                         + ar2*xi22 + ai2*xr22 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            ar2 = colA2[rloc] ; ai2 = colA2[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 
                          + ar1*xr10 - ai1*xi10 
                          + ar2*xr20 - ai2*xi20 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 
                          + ar1*xi10 + ai1*xr10 
                          + ar2*xi20 + ai2*xr20 ;
            colY1[ryloc] -= ar0*xr01 - ai0*xi01 
                          + ar1*xr11 - ai1*xi11 
                          + ar2*xr21 - ai2*xi21 ;
            colY1[iyloc] -= ar0*xi01 + ai0*xr01 
                          + ar1*xi11 + ai1*xr11 
                          + ar2*xi21 + ai2*xr21 ;
            colY2[ryloc] -= ar0*xr02 - ai0*xi02 
                          + ar1*xr12 - ai1*xi12 
                          + ar2*xr22 - ai2*xi22 ;
            colY2[iyloc] -= ar0*xi02 + ai0*xr02 
                          + ar1*xi12 + ai1*xr12 
                          + ar2*xi22 + ai2*xr22 ;
         }
      }
      colA0 = colA2 + 2*nrowA ;
   }
   if ( icolA == ncolA - 2 ) {
      colA1 = colA0 + 2*nrowA ;
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc += 2, iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc = 2*colindA[icolA+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
            colY2[rloc] -= ar0*xr02 - ai0*xi02 + ar1*xr12 - ai1*xi12 ;
            colY2[iloc] -= ar0*xi02 + ai0*xr02 + ar1*xi12 + ai1*xr12 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[ryloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iyloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
            colY2[ryloc] -= ar0*xr02 - ai0*xi02 + ar1*xr12 - ai1*xi12 ;
            colY2[iyloc] -= ar0*xi02 + ai0*xr02 + ar1*xi12 + ai1*xr12 ;
         }
      }
   } else if ( icolA == ncolA - 1 ) {
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 ;
            colY2[rloc] -= ar0*xr02 - ai0*xi02 ;
            colY2[iloc] -= ar0*xi02 + ai0*xr02 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[ryloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iyloc] -= ar0*xi01 + ai0*xr01 ;
            colY2[ryloc] -= ar0*xr02 - ai0*xi02 ;
            colY2[iyloc] -= ar0*xi02 + ai0*xr02 ;
         }
      }
   }
   colX0 = colX2 + 2*nrowX ;
   colY0 = colY2 + 2*nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colA0 = entA ;
   for ( icolA = 0 ; icolA < ncolA - 2 ; icolA += 3 ) {
      colA1 = colA0 + 2*nrowA ;
      colA2 = colA1 + 2*nrowA ;
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc += 2, iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc += 2, iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc = 2*colindA[icolA+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc = 2*colindA[icolA+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            ar2 = colA2[rloc] ; ai2 = colA2[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 
                         + ar1*xr10 - ai1*xi10 
                         + ar2*xr20 - ai2*xi20 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 
                         + ar1*xi10 + ai1*xr10 
                         + ar2*xi20 + ai2*xr20 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 
                         + ar1*xr11 - ai1*xi11 
                         + ar2*xr21 - ai2*xi21 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 
                         + ar1*xi11 + ai1*xr11 
                         + ar2*xi21 + ai2*xr21 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            ar2 = colA2[rloc] ; ai2 = colA2[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 
                          + ar1*xr10 - ai1*xi10 
                          + ar2*xr20 - ai2*xi20 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 
                          + ar1*xi10 + ai1*xr10 
                          + ar2*xi20 + ai2*xr20 ;
            colY1[ryloc] -= ar0*xr01 - ai0*xi01 
                          + ar1*xr11 - ai1*xi11 
                          + ar2*xr21 - ai2*xi21 ;
            colY1[iyloc] -= ar0*xi01 + ai0*xr01 
                          + ar1*xi11 + ai1*xr11 
                          + ar2*xi21 + ai2*xr21 ;
         }
      }
      colA0 = colA2 + 2*nrowA ;
   }
   if ( icolA == ncolA - 2 ) {
      colA1 = colA0 + 2*nrowA ;
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc += 2, iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc = 2*colindA[icolA+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[ryloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iyloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
         }
      }
   } else if ( icolA == ncolA - 1 ) {
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[ryloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iyloc] -= ar0*xi01 + ai0*xr01 ;
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   colA0 = entA ;
   for ( icolA = 0 ; icolA < ncolA - 2 ; icolA += 3 ) {
      colA1 = colA0 + 2*nrowA ;
      colA2 = colA1 + 2*nrowA ;
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc += 2, iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc += 2, iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc = 2*colindA[icolA+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc = 2*colindA[icolA+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            ar2 = colA2[rloc] ; ai2 = colA2[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 
                         + ar1*xr10 - ai1*xi10 
                         + ar2*xr20 - ai2*xi20 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 
                         + ar1*xi10 + ai1*xr10 
                         + ar2*xi20 + ai2*xr20 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            ar2 = colA2[rloc] ; ai2 = colA2[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 
                          + ar1*xr10 - ai1*xi10 
                          + ar2*xr20 - ai2*xi20 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 
                          + ar1*xi10 + ai1*xr10 
                          + ar2*xi20 + ai2*xr20 ;
         }
      }
      colA0 = colA2 + 2*nrowA ;
   }
   if ( icolA == ncolA - 2 ) {
      colA1 = colA0 + 2*nrowA ;
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc += 2, iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc = 2*colindA[icolA+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            ar1 = colA1[rloc] ; ai1 = colA1[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
         }
      }
   } else if ( icolA == ncolA - 1 ) {
      if ( ncolA == nrowX ) {
         rloc = 2*icolA ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
      } else {
         rloc = 2*colindA[icolA] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
      }
      if ( nrowY == nrowA ) {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
         }
      } else {
         for ( krowA = 0, rloc = 0, iloc = 1 ; 
               krowA < nrowA ; 
               krowA++, rloc += 2, iloc += 2 ) {
            ar0 = colA0[rloc] ; ai0 = colA0[iloc] ;
            krowY = rowindA[krowA] ;
            ryloc = 2*rowindA[krowA] ; iyloc = ryloc + 1 ;
            colY0[ryloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iyloc] -= ar0*xi00 + ai0*xr00 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------
   A has dense rows
   ----------------
*/
static void
complex_updDenseRows (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2, 
         *rowA0, *rowA1, *rowA2, *entA, *entX, *entY ;
int      inc1, inc2, irowA, jcolX, kcolA, 
         ncolA, ncolX, ncolY, nrowA, nrowX, nrowY ;
int      *colindA, *rowindA ;

SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &nrowA, &ncolA, &inc1, &inc2, &entA) ;
if ( ncolA != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( nrowA != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + 2*nrowX ;
   colX2 = colX1 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colY2 = colY1 + 2*nrowY ;
   rowA0 = entA ;
   for ( irowA = 0 ; irowA < nrowA - 2 ; irowA += 3 ) {
      double   ai0, ai1, ai2, ar0, ar1, ar2, 
               xi0, xi1, xi2, xr0, xr1, xr2,
               isum00, isum01, isum02, isum10, isum11, isum12,
               isum20, isum21, isum22, rsum00, rsum01, rsum02,
               rsum10, rsum11, rsum12, rsum20, rsum21, rsum22 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      rowA1 = rowA0 + 2*ncolA ;
      rowA2 = rowA1 + 2*ncolA ;
      isum00 = isum01 = isum02 =
      isum10 = isum11 = isum12 =
      isum20 = isum21 = isum22 = 0.0 ;
      rsum00 = rsum01 = rsum02 =
      rsum10 = rsum11 = rsum12 =
      rsum20 = rsum21 = rsum22 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            ar1 = rowA1[rloc] ; ai1 = rowA1[iloc] ;
            ar2 = rowA2[rloc] ; ai2 = rowA2[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
            xr2 = colX2[rloc] ; xi2 = colX2[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum12 += ar1*xr2 - ai1*xi2 ; isum12 += ar1*xi2 + ai1*xr2 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
            rsum21 += ar2*xr1 - ai2*xi1 ; isum21 += ar2*xi1 + ai2*xr1 ;
            rsum22 += ar2*xr2 - ai2*xi2 ; isum22 += ar2*xi2 + ai2*xr2 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            ar1 = rowA1[raloc] ; ai1 = rowA1[ialoc] ;
            ar2 = rowA2[raloc] ; ai2 = rowA2[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc] ; xi1 = colX1[ixloc] ;
            xr2 = colX2[rxloc] ; xi2 = colX2[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum12 += ar1*xr2 - ai1*xi2 ; isum12 += ar1*xi2 + ai1*xr2 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
            rsum21 += ar2*xr1 - ai2*xi1 ; isum21 += ar2*xi1 + ai2*xr1 ;
            rsum22 += ar2*xr2 - ai2*xi2 ; isum22 += ar2*xi2 + ai2*xr2 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         colY2[ryloc] -= rsum02 ; colY2[iyloc] -= isum02 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
         colY2[ryloc] -= rsum12 ; colY2[iyloc] -= isum12 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum20 ; colY0[iyloc] -= isum20 ;
         colY1[ryloc] -= rsum21 ; colY1[iyloc] -= isum21 ;
         colY2[ryloc] -= rsum22 ; colY2[iyloc] -= isum22 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         colY2[ryloc] -= rsum02 ; colY2[iyloc] -= isum02 ;
         ryloc = 2*rowindA[irowA+1] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
         colY2[ryloc] -= rsum12 ; colY2[iyloc] -= isum12 ;
         ryloc = 2*rowindA[irowA+2] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum20 ; colY0[iyloc] -= isum20 ;
         colY1[ryloc] -= rsum21 ; colY1[iyloc] -= isum21 ;
         colY2[ryloc] -= rsum22 ; colY2[iyloc] -= isum22 ;
      }
      rowA0 = rowA2 + 2*ncolA ;
   }
   if ( irowA == nrowA - 2 ) {
      double   ai0, ai1, ar0, ar1, xi0, xi1, xi2, xr0, xr1, xr2,
               isum00, isum01, isum02, isum10, isum11, isum12,
               rsum00, rsum01, rsum02, rsum10, rsum11, rsum12 ; 
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      rowA1 = rowA0 + 2*ncolA ;
      isum00 = isum01 = isum02 =
      isum10 = isum11 = isum12 = 0.0 ;
      rsum00 = rsum01 = rsum02 =
      rsum10 = rsum11 = rsum12 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            ar1 = rowA1[rloc] ; ai1 = rowA1[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
            xr2 = colX2[rloc] ; xi2 = colX2[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum12 += ar1*xr2 - ai1*xi2 ; isum12 += ar1*xi2 + ai1*xr2 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            ar1 = rowA1[raloc] ; ai1 = rowA1[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc] ; xi1 = colX1[ixloc] ;
            xr2 = colX2[rxloc] ; xi2 = colX2[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum12 += ar1*xr2 - ai1*xi2 ; isum12 += ar1*xi2 + ai1*xr2 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         colY2[ryloc] -= rsum02 ; colY2[iyloc] -= isum02 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
         colY2[ryloc] -= rsum12 ; colY2[iyloc] -= isum12 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         colY2[ryloc] -= rsum02 ; colY2[iyloc] -= isum02 ;
         ryloc = 2*rowindA[irowA+1] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
         colY2[ryloc] -= rsum12 ; colY2[iyloc] -= isum12 ;
      }
   } else if ( irowA == nrowA - 1 ) {
      double   ai0, ar0, xi0, xi1, xi2, xr0, xr1, xr2,
               isum00, isum01, isum02, rsum00, rsum01, rsum02 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      isum00 = isum01 = isum02 = rsum00 = rsum01 = rsum02 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
            xr2 = colX2[rloc] ; xi2 = colX2[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc] ; xi1 = colX1[ixloc] ;
            xr2 = colX2[rxloc] ; xi2 = colX2[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         colY2[ryloc] -= rsum02 ; colY2[iyloc] -= isum02 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         colY2[ryloc] -= rsum02 ; colY2[iyloc] -= isum02 ;
      }
   }
   colX0 = colX2 + 2*nrowX ;
   colY0 = colY2 + 2*nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   rowA0 = entA ;
   for ( irowA = 0 ; irowA < nrowA - 2 ; irowA += 3 ) {
      double   ai0, ai1, ai2, ar0, ar1, ar2, xi0, xi1, xr0, xr1, 
               isum00, isum01, isum10, isum11, 
               isum20, isum21, rsum00, rsum01, 
               rsum10, rsum11, rsum20, rsum21 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      rowA1 = rowA0 + 2*ncolA ;
      rowA2 = rowA1 + 2*ncolA ;
      isum00 = isum01 = isum10 = isum11 = isum20 = isum21 = 0.0 ;
      rsum00 = rsum01 = rsum10 = rsum11 = rsum20 = rsum21 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            ar1 = rowA1[rloc] ; ai1 = rowA1[iloc] ;
            ar2 = rowA2[rloc] ; ai2 = rowA2[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
            rsum21 += ar2*xr1 - ai2*xi1 ; isum21 += ar2*xi1 + ai2*xr1 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            ar1 = rowA1[raloc] ; ai1 = rowA1[ialoc] ;
            ar2 = rowA2[raloc] ; ai2 = rowA2[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc] ; xi1 = colX1[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
            rsum21 += ar2*xr1 - ai2*xi1 ; isum21 += ar2*xi1 + ai2*xr1 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum20 ; colY0[iyloc] -= isum20 ;
         colY1[ryloc] -= rsum21 ; colY1[iyloc] -= isum21 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         ryloc = 2*rowindA[irowA+1] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
         ryloc = 2*rowindA[irowA+2] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum20 ; colY0[iyloc] -= isum20 ;
         colY1[ryloc] -= rsum21 ; colY1[iyloc] -= isum21 ;
      }
      rowA0 = rowA2 + 2*ncolA ;
   }
   if ( irowA == nrowA - 2 ) {
      double   ai0, ai1, ar0, ar1, xi0, xi1, xr0, xr1, 
               isum00, isum01, isum10, isum11, rsum00, rsum01, 
               rsum10, rsum11 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      rowA1 = rowA0 + 2*ncolA ;
      isum00 = isum01 = isum10 = isum11 = 0.0 ;
      rsum00 = rsum01 = rsum10 = rsum11 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            ar1 = rowA1[rloc] ; ai1 = rowA1[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            ar1 = rowA1[raloc] ; ai1 = rowA1[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc] ; xi1 = colX1[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
         ryloc = 2*rowindA[irowA+1] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         colY1[ryloc] -= rsum11 ; colY1[iyloc] -= isum11 ;
      }
      rowA0 = rowA2 + 2*ncolA ;
   } else if ( irowA == nrowA - 1 ) {
      double   ai0, ar0, xi0, xi1, xr0, xr1, 
               isum00, isum01, rsum00, rsum01 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      isum00 = isum01 = 0.0 ;
      rsum00 = rsum01 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc] ; xi1 = colX1[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         colY1[ryloc] -= rsum01 ; colY1[iyloc] -= isum01 ;
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   rowA0 = entA ;
   for ( irowA = 0 ; irowA < nrowA - 2 ; irowA += 3 ) {
      double   ai0, ai1, ai2, ar0, ar1, ar2, xi0, xr0,
               isum00, isum10, isum20, rsum00, rsum10, rsum20 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      rowA1 = rowA0 + 2*ncolA ;
      rowA2 = rowA1 + 2*ncolA ;
      isum00 = isum10 = isum20 = 0.0 ;
      rsum00 = rsum10 = rsum20 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            ar1 = rowA1[rloc] ; ai1 = rowA1[iloc] ;
            ar2 = rowA2[rloc] ; ai2 = rowA2[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            ar1 = rowA1[raloc] ; ai1 = rowA1[ialoc] ;
            ar2 = rowA2[raloc] ; ai2 = rowA2[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum20 ; colY0[iyloc] -= isum20 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         ryloc = 2*rowindA[irowA+1] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
         ryloc = 2*rowindA[irowA+2] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum20 ; colY0[iyloc] -= isum20 ;
      }
      rowA0 = rowA2 + 2*ncolA ;
   }
   if ( irowA == nrowA - 2 ) {
      double   ai0, ai1, ar0, ar1, xi0, xr0,
               isum00, isum10, rsum00, rsum10 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      rowA1 = rowA0 + 2*ncolA ;
      isum00 = isum10 = 0.0 ;
      rsum00 = rsum10 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            ar1 = rowA1[rloc] ; ai1 = rowA1[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            ar1 = rowA1[raloc] ; ai1 = rowA1[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         ryloc += 2, iyloc += 2 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
         ryloc = 2*rowindA[irowA+1] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum10 ; colY0[iyloc] -= isum10 ;
      }
   } else if ( irowA == nrowA - 1 ) {
      double   ai0, ar0, xi0, xr0, isum00, rsum00 ;
      int      ialoc, iloc, ixloc, iyloc, raloc, rloc, rxloc, ryloc ;

      isum00 = 0.0 ;
      rsum00 = 0.0 ;
      if ( ncolA == nrowX ) {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            rloc = 2*kcolA ; iloc = rloc + 1 ;
            ar0 = rowA0[rloc] ; ai0 = rowA0[iloc] ;
            xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
         }
      } else {
         for ( kcolA = 0 ; kcolA < ncolA ; kcolA++ ) {
            raloc = 2*kcolA ; ialoc = raloc + 1 ;
            ar0 = rowA0[raloc] ; ai0 = rowA0[ialoc] ;
            rxloc = 2*colindA[kcolA] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc] ; xi0 = colX0[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
         }
      }
      if ( nrowY == nrowA ) {
         ryloc = 2*irowA ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
      } else {
         ryloc = 2*rowindA[irowA] ; iyloc = ryloc + 1 ;
         colY0[ryloc] -= rsum00 ; colY0[iyloc] -= isum00 ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------
   A has sparse rows
   -----------------
*/
static void
complex_updSparseRows (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   ai, ar, xi0, isum0, isum1, isum2, 
         xi1, xi2, xr0, xr1, xr2, rsum0, rsum1, rsum2 ;
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, iloc, inc1, inc2, irowA, jcolX, kk, krowX, 
         ncolA, ncolX, ncolY, nentA, nrowA, nrowX, nrowY, rloc, size ;
int      *colindA, *indices, *rowindA, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_ROWS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseRowsInfo(mtxA, &nrowA, &nentA, &sizes, &indices, &entA) ;
if ( (ncolA = mtxA->ncol) != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( nrowA != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + 2*nrowX ;
   colX2 = colX1 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colY2 = colY1 + 2*nrowY ;
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizes[irowA]) > 0 ) {
         isum0 = isum1 = isum2 = rsum0 = rsum1 = rsum2 = 0.0 ;
         if ( ncolA == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               krowX = indices[kk] ;
               xr0 = colX0[2*krowX] ; xi0 = colX0[2*krowX+1] ;
               xr1 = colX1[2*krowX] ; xi1 = colX1[2*krowX+1] ;
               xr2 = colX2[2*krowX] ; xi2 = colX2[2*krowX+1] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
               rsum2 += ar*xr2 - ai*xi2 ; isum2 += ar*xi2 + ai*xr2 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               krowX = colindA[indices[kk]] ;
               xr0 = colX0[2*krowX] ; xi0 = colX0[2*krowX+1] ;
               xr1 = colX1[2*krowX] ; xi1 = colX1[2*krowX+1] ;
               xr2 = colX2[2*krowX] ; xi2 = colX2[2*krowX+1] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
               rsum2 += ar*xr2 - ai*xi2 ; isum2 += ar*xi2 + ai*xr2 ;
            }
         }
         if ( nrowA == nrowY ) {
            rloc = 2*irowA ; iloc = rloc + 1 ;
         } else {
            rloc = 2*rowindA[irowA] ; iloc = rloc + 1 ;
         }
         colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
         colY1[rloc] -= rsum1 ; colY1[iloc] -= isum1 ;
         colY2[rloc] -= rsum2 ; colY2[iloc] -= isum2 ;
      }
   }
   colX0 = colX2 + 2*nrowX ;
   colY0 = colY2 + 2*nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizes[irowA]) > 0 ) {
         isum0 = isum1 = rsum0 = rsum1 = 0.0 ;
         if ( ncolA == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               krowX = indices[kk] ;
               xr0 = colX0[2*krowX] ; xi0 = colX0[2*krowX+1] ;
               xr1 = colX1[2*krowX] ; xi1 = colX1[2*krowX+1] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               krowX = colindA[indices[kk]] ;
               xr0 = colX0[2*krowX] ; xi0 = colX0[2*krowX+1] ;
               xr1 = colX1[2*krowX] ; xi1 = colX1[2*krowX+1] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
            }
         }
         if ( nrowA == nrowY ) {
            rloc = 2*irowA ; iloc = rloc + 1 ;
         } else {
            rloc = 2*rowindA[irowA] ; iloc = rloc + 1 ;
         }
         colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
         colY1[rloc] -= rsum1 ; colY1[iloc] -= isum1 ;
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   for ( irowA = kk = 0 ; irowA < nrowA ; irowA++ ) {
      if ( (size = sizes[irowA]) > 0 ) {
         isum0 = rsum0 = 0.0 ;
         if ( ncolA == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               krowX = indices[kk] ;
               xr0 = colX0[2*krowX] ; xi0 = colX0[2*krowX+1] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               krowX = colindA[indices[kk]] ;
               xr0 = colX0[2*krowX] ; xi0 = colX0[2*krowX+1] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
            }
         }
         if ( nrowA == nrowY ) {
            rloc = 2*irowA ; iloc = rloc + 1 ;
         } else {
            rloc = 2*rowindA[irowA] ; iloc = rloc + 1 ;
         }
         colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------
   A has sparse columns
   --------------------
*/
static void
complex_updSparseColumns (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) {
double   ai, ar, xi0, xi1, xi2, xr0, xr1, xr2 ;
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, iloc, inc1, inc2, jcolA, jcolX, jrowX, kk, 
         ncolA, ncolX, ncolY, nentA, nrowA, nrowX, nrowY, rloc, size ;
int      *colindA, *indices, *rowindA, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_COLUMNS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseColumnsInfo(mtxA, &ncolA, &nentA, &sizes, &indices, &entA) ;
if ( ncolA != nrowX ) {
   SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
} else {
   colindA = NULL ;
}
if ( (nrowA = mtxA->nrow) != nrowY ) {
   SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
} else {
   rowindA = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + 2*nrowX ;
   colX2 = colX1 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colY2 = colY1 + 2*nrowY ;
   for ( jcolA = kk = 0 ; jcolA < ncolA ; jcolA++ ) {
      if ( (size = sizes[jcolA]) > 0 ) {
         if ( ncolA == nrowX ) {
            jrowX = jcolA ;
         } else {
            jrowX = colindA[jcolA] ;
         }
         rloc = 2*jrowX ; iloc = rloc + 1 ;
         xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
         xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
         xr2 = colX2[rloc] ; xi2 = colX2[iloc] ;
         if ( nrowA == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*indices[kk] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
               colY1[rloc] -= ar*xr1 - ai*xi1 ;
               colY1[iloc] -= ar*xi1 + ai*xr1 ;
               colY2[rloc] -= ar*xr2 - ai*xi2 ;
               colY2[iloc] -= ar*xi2 + ai*xr2 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*rowindA[indices[kk]] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
               colY1[rloc] -= ar*xr1 - ai*xi1 ;
               colY1[iloc] -= ar*xi1 + ai*xr1 ;
               colY2[rloc] -= ar*xr2 - ai*xi2 ;
               colY2[iloc] -= ar*xi2 + ai*xr2 ;
            }
         }
      }
   }
   colX0 = colX2 + 2*nrowX ;
   colY0 = colY2 + 2*nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   for ( jcolA = kk = 0 ; jcolA < ncolA ; jcolA++ ) {
      if ( (size = sizes[jcolA]) > 0 ) {
         if ( ncolA == nrowX ) {
            jrowX = jcolA ;
         } else {
            jrowX = colindA[jcolA] ;
         }
         rloc = 2*jrowX ; iloc = rloc + 1 ;
         xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
         xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
         if ( nrowA == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*indices[kk] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
               colY1[rloc] -= ar*xr1 - ai*xi1 ;
               colY1[iloc] -= ar*xi1 + ai*xr1 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*rowindA[indices[kk]] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
               colY1[rloc] -= ar*xr1 - ai*xi1 ;
               colY1[iloc] -= ar*xi1 + ai*xr1 ;
            }
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   for ( jcolA = kk = 0 ; jcolA < ncolA ; jcolA++ ) {
      if ( (size = sizes[jcolA]) > 0 ) {
         if ( ncolA == nrowX ) {
            jrowX = jcolA ;
         } else {
            jrowX = colindA[jcolA] ;
         }
         rloc = 2*jrowX ; iloc = rloc + 1 ;
         xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
         if ( nrowA == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*indices[kk] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*rowindA[indices[kk]] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
            }
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
