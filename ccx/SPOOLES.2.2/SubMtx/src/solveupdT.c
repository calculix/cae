/*  solveupdT.c  */

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
   ------------------------------------------------------
   purpose -- perform the matrix-matrix multiply 
      Y := Y - A^T * X used in the forward and backsolves
      where
        (1) rows(A) \subseteq rows(Y)
        (2) rows(A) are local w.r.t. rows(Y)
        (3) cols(A) \subseteq rows(X)
        (4) cols(A) are local w.r.t. rows(X)
        (5) cols(Y) = cols(X)
        (6) Y and X have type SUBMTX_DENSE_COLUMNS
 
   created -- 98may02, cca
   -----------------------------------------------------
*/
void
SubMtx_solveupdT (
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtxY == NULL || mtxA == NULL || mtxX == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveupdT(%p,%p,%p)"
           "\n bad input\n", mtxY, mtxA, mtxX) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_DENSE_COLUMNS(mtxY) ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveupdT(%p,%p,%p)"
           "\n Y must have mode SUBMTX_DENSE_COLUMNS\n", 
           mtxY, mtxA, mtxX) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_DENSE_COLUMNS(mtxX) ) {
   fprintf(stderr, "\n fatal error in SubMtx_solveupdT(%p,%p,%p)"
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
      fprintf(stderr, "\n fatal error in SubMtx_solveupdT(%p,%p,%p)"
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
      fprintf(stderr, "\n fatal error in SubMtx_solveupdT(%p,%p,%p)"
              "\n unsupported mode %d for A\n",
              mtxY, mtxA, mtxX, mtxA->mode) ;
      exit(-1) ;
      break ;
   }
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_solveupdT(%p,%p,%p)"
           "\n unsupported type %d for A\n",
           mtxY, mtxA, mtxX, mtxA->type) ;
   exit(-1) ;
   break ;
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
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   Ak0, Ak1, Ak2, x00, x01, x02, x10, x11, x12,
         x20, x21, x22 ;
double   *colAT0, *colAT1, *colAT2, *colX0, *colX1, *colX2, 
         *colY0, *colY1, *colY2, *entA, *entX, *entY ;
int      icolAT, inc1, inc2, irowX, jcolX, krowAT, krowY, 
         ncolAT, ncolX, ncolY, nrowAT, nrowX, nrowY ;
int      *colindAT, *rowindAT ;

SubMtx_denseInfo(mtxY, &nrowY,  &ncolY,  &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX,  &ncolX,  &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &ncolAT, &nrowAT, &inc1, &inc2, &entA) ;
colX0 = entX ;
colY0 = entY ;
if ( ncolAT != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( nrowAT != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1  = colX0 + nrowX ;
   colX2  = colX1 + nrowX ;
   colY1  = colY0 + nrowY ;
   colY2  = colY1 + nrowY ;
   colAT0 = entA ;
   for ( icolAT = 0 ; icolAT < ncolAT - 2 ; icolAT += 3 ) {
      colAT1 = colAT0 + nrowAT ;
      colAT2 = colAT1 + nrowAT ;
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x01 = colX1[icolAT] ;
         x02 = colX2[icolAT] ;
         x10 = colX0[icolAT+1] ;
         x11 = colX1[icolAT+1] ;
         x12 = colX2[icolAT+1] ;
         x20 = colX0[icolAT+2] ;
         x21 = colX1[icolAT+2] ;
         x22 = colX2[icolAT+2] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         x02 = colX2[irowX] ;
         irowX = colindAT[icolAT+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
         x12 = colX2[irowX] ;
         irowX = colindAT[icolAT+2] ;
         x20 = colX0[irowX] ;
         x21 = colX1[irowX] ;
         x22 = colX2[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            Ak2 = colAT2[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowAT] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
            colY2[krowAT] -= Ak0 * x02 + Ak1 * x12 + Ak2 * x22 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            Ak2 = colAT2[krowAT] ;
            krowY = rowindAT[krowAT] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
            colY2[krowY] -= Ak0 * x02 + Ak1 * x12 + Ak2 * x22 ;
         }
      }
      colAT0 = colAT2 + nrowAT ;
   }
   if ( icolAT == ncolAT - 2 ) {
      colAT1 = colAT0 + nrowAT ;
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x01 = colX1[icolAT] ;
         x02 = colX2[icolAT] ;
         x10 = colX0[icolAT+1] ;
         x11 = colX1[icolAT+1] ;
         x12 = colX2[icolAT+1] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         x02 = colX2[irowX] ;
         irowX = colindAT[icolAT+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
         x12 = colX2[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowAT] -= Ak0 * x01 + Ak1 * x11 ;
            colY2[krowAT] -= Ak0 * x02 + Ak1 * x12 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            krowY = rowindAT[krowAT] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 ;
            colY2[krowY] -= Ak0 * x02 + Ak1 * x12 ;
         }
      }
   } else if ( icolAT == ncolAT - 1 ) {
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x01 = colX1[icolAT] ;
         x02 = colX2[icolAT] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         x02 = colX2[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 ;
            colY1[krowAT] -= Ak0 * x01 ;
            colY2[krowAT] -= Ak0 * x02 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            krowY = rowindAT[krowAT] ;
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
   colAT0 = entA ;
   for ( icolAT = 0 ; icolAT < ncolAT - 2 ; icolAT += 3 ) {
      colAT1 = colAT0 + nrowAT ;
      colAT2 = colAT1 + nrowAT ;
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x01 = colX1[icolAT] ;
         x10 = colX0[icolAT+1] ;
         x11 = colX1[icolAT+1] ;
         x20 = colX0[icolAT+2] ;
         x21 = colX1[icolAT+2] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         irowX = colindAT[icolAT+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
         irowX = colindAT[icolAT+2] ;
         x20 = colX0[irowX] ;
         x21 = colX1[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            Ak2 = colAT2[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowAT] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            Ak2 = colAT2[krowAT] ;
            krowY = rowindAT[krowAT] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 + Ak2 * x21 ;
         }
      }
      colAT0 = colAT2 + nrowAT ;
   }
   if ( icolAT == ncolAT - 2 ) {
      colAT1 = colAT0 + nrowAT ;
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x01 = colX1[icolAT] ;
         x10 = colX0[icolAT+1] ;
         x11 = colX1[icolAT+1] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
         irowX = colindAT[icolAT+1] ;
         x10 = colX0[irowX] ;
         x11 = colX1[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowAT] -= Ak0 * x01 + Ak1 * x11 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            krowY = rowindAT[krowAT] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 ;
            colY1[krowY] -= Ak0 * x01 + Ak1 * x11 ;
         }
      }
   } else if ( icolAT == ncolAT - 1 ) {
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x01 = colX1[icolAT] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         x01 = colX1[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 ;
            colY1[krowAT] -= Ak0 * x01 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            krowY = rowindAT[krowAT] ;
            colY0[krowY] -= Ak0 * x00 ;
            colY1[krowY] -= Ak0 * x01 ;
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   colAT0 = entA ;
   for ( icolAT = 0 ; icolAT < ncolAT - 2 ; icolAT += 3 ) {
      colAT1 = colAT0 + nrowAT ;
      colAT2 = colAT1 + nrowAT ;
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x10 = colX0[icolAT+1] ;
         x20 = colX0[icolAT+2] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         irowX = colindAT[icolAT+1] ;
         x10 = colX0[irowX] ;
         irowX = colindAT[icolAT+2] ;
         x20 = colX0[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            Ak2 = colAT2[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            Ak2 = colAT2[krowAT] ;
            krowY = rowindAT[krowAT] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 + Ak2 * x20 ;
         }
      }
      colAT0 = colAT2 + nrowAT ;
   }
   if ( icolAT == ncolAT - 2 ) {
      colAT1 = colAT0 + nrowAT ;
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
         x10 = colX0[icolAT+1] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
         irowX = colindAT[icolAT+1] ;
         x10 = colX0[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 + Ak1 * x10 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            Ak1 = colAT1[krowAT] ;
            krowY = rowindAT[krowAT] ;
            colY0[krowY] -= Ak0 * x00 + Ak1 * x10 ;
         }
      }
   } else if ( icolAT == ncolAT - 1 ) {
      if ( ncolAT == nrowX ) {
         x00 = colX0[icolAT] ;
      } else {
         irowX = colindAT[icolAT] ;
         x00 = colX0[irowX] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            colY0[krowAT] -= Ak0 * x00 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            Ak0 = colAT0[krowAT] ;
            krowY = rowindAT[krowAT] ;
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
real_updDenseColumns (
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2, 
         *rowAT0, *rowAT1, *rowAT2, *entA, *entX, *entY ;
int      inc1, inc2, irowAT, irowY, jcolX, kcolAT, krowX, 
         ncolAT, ncolX, ncolY, nrowAT, nrowX, nrowY ;
int      *colindAT, *rowindAT ;

SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &ncolAT, &nrowAT, &inc1, &inc2, &entA) ;
if ( ncolAT != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( nrowAT != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + nrowX ;
   colX2 = colX1 + nrowX ;
   colY1 = colY0 + nrowY ;
   colY2 = colY1 + nrowY ;
   rowAT0 = entA ;
   for ( irowAT = 0 ; irowAT < nrowAT - 2 ; irowAT += 3 ) {
      double   A0k, A1k, A2k, Xk0, Xk1, Xk2 ;
      double   sum00, sum01, sum02, 
               sum10, sum11, sum12, 
               sum20, sum21, sum22 ;

      sum00 = sum01 = sum02 
            = sum10 = sum11 = sum12 = sum20 = sum21 = sum22 = 0.0 ;
      rowAT1 = rowAT0 + ncolAT ;
      rowAT2 = rowAT1 + ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            A2k = rowAT2[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            Xk1 = colX1[kcolAT] ;
            Xk2 = colX2[kcolAT] ;
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
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            A2k = rowAT2[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
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
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY1[irowAT]   -= sum01 ;
         colY2[irowAT]   -= sum02 ;
         colY0[irowAT+1] -= sum10 ;
         colY1[irowAT+1] -= sum11 ;
         colY2[irowAT+1] -= sum12 ;
         colY0[irowAT+2] -= sum20 ;
         colY1[irowAT+2] -= sum21 ;
         colY2[irowAT+2] -= sum22 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         colY2[irowY] -= sum02 ;
         irowY = rowindAT[irowAT+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
         colY2[irowY] -= sum12 ;
         irowY = rowindAT[irowAT+2] ;
         colY0[irowY] -= sum20 ;
         colY1[irowY] -= sum21 ;
         colY2[irowY] -= sum22 ;
      }
      rowAT0 = rowAT2 + ncolAT ;
   }
   if ( irowAT == nrowAT - 2 ) {
      double   A0k, A1k, Xk0, Xk1, Xk2 ;
      double   sum00, sum01, sum02, sum10, sum11, sum12 ;

      sum00 = sum01 = sum02 = sum10 = sum11 = sum12 = 0.0 ;
      rowAT1 = rowAT0 + ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            Xk1 = colX1[kcolAT] ;
            Xk2 = colX2[kcolAT] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum12 += A1k * Xk2 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
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
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY1[irowAT]   -= sum01 ;
         colY2[irowAT]   -= sum02 ;
         colY0[irowAT+1] -= sum10 ;
         colY1[irowAT+1] -= sum11 ;
         colY2[irowAT+1] -= sum12 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         colY2[irowY] -= sum02 ;
         irowY = rowindAT[irowAT+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
         colY2[irowY] -= sum12 ;
      }
   } else if ( irowAT == nrowAT - 1 ) {
      double   A0k, Xk0, Xk1, Xk2 ;
      double   sum00, sum01, sum02 ;

      sum00 = sum01 = sum02 = 0.0 ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            Xk1 = colX1[kcolAT] ;
            Xk2 = colX2[kcolAT] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            Xk2 = colX2[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum02 += A0k * Xk2 ;
         }
      }
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY1[irowAT]   -= sum01 ;
         colY2[irowAT]   -= sum02 ;
      } else {
         irowY = rowindAT[irowAT] ;
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
   rowAT0 = entA ;
   for ( irowAT = 0 ; irowAT < nrowAT - 2 ; irowAT += 3 ) {
      double   A0k, A1k, A2k, Xk0, Xk1 ;
      double   sum00, sum01, sum10, sum11, sum20, sum21 ;

      sum00 = sum01 = sum10 = sum11 = sum20 = sum21 = 0.0 ;
      rowAT1 = rowAT0 + ncolAT ;
      rowAT2 = rowAT1 + ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            A2k = rowAT2[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            Xk1 = colX1[kcolAT] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
            sum20 += A2k * Xk0 ;
            sum21 += A2k * Xk1 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            A2k = rowAT2[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
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
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY1[irowAT]   -= sum01 ;
         colY0[irowAT+1] -= sum10 ;
         colY1[irowAT+1] -= sum11 ;
         colY0[irowAT+2] -= sum20 ;
         colY1[irowAT+2] -= sum21 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         irowY = rowindAT[irowAT+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
         irowY = rowindAT[irowAT+2] ;
         colY0[irowY] -= sum20 ;
         colY1[irowY] -= sum21 ;
      }
      rowAT0 = rowAT2 + ncolAT ;
   }
   if ( irowAT == nrowAT - 2 ) {
      double   A0k, A1k, Xk0, Xk1 ;
      double   sum00, sum01, sum10, sum11 ;

      sum00 = sum01 = sum10 = sum11 = 0.0 ;
      rowAT1 = rowAT0 + ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            Xk1 = colX1[kcolAT] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
            sum10 += A1k * Xk0 ;
            sum11 += A1k * Xk1 ;
         }
      }
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY1[irowAT]   -= sum01 ;
         colY0[irowAT+1] -= sum10 ;
         colY1[irowAT+1] -= sum11 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
         irowY = rowindAT[irowAT+1] ;
         colY0[irowY] -= sum10 ;
         colY1[irowY] -= sum11 ;
      }
   } else if ( irowAT == nrowAT - 1 ) {
      double   A0k, Xk0, Xk1 ;
      double   sum00, sum01 ;

      sum00 = sum01 = 0.0 ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            Xk1 = colX1[kcolAT] ;
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
            Xk0 = colX0[krowX] ; 
            Xk1 = colX1[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum01 += A0k * Xk1 ;
         }
      }
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY1[irowAT]   -= sum01 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
         colY1[irowY] -= sum01 ;
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   rowAT0 = entA ;
   for ( irowAT = 0 ; irowAT < nrowAT - 2 ; irowAT += 3 ) {
      double   A0k, A1k, A2k, Xk0 ;
      double   sum00, sum10, sum20 ;

      sum00 = sum10 = sum20 = 0.0 ;
      rowAT1 = rowAT0 + ncolAT ;
      rowAT2 = rowAT1 + ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            A2k = rowAT2[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
            sum20 += A2k * Xk0 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            A2k = rowAT2[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
            Xk0 = colX0[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
            sum20 += A2k * Xk0 ;
         }
      }
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY0[irowAT+1] -= sum10 ;
         colY0[irowAT+2] -= sum20 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
         irowY = rowindAT[irowAT+1] ;
         colY0[irowY] -= sum10 ;
         irowY = rowindAT[irowAT+2] ;
         colY0[irowY] -= sum20 ;
      }
      rowAT0 = rowAT2 + ncolAT ;
   }
   if ( irowAT == nrowAT - 2 ) {
      double   A0k, A1k, Xk0 ;
      double   sum00, sum10 ;

      sum00 = sum10 = 0.0 ;
      rowAT1 = rowAT0 + ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            A1k = rowAT1[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
            Xk0 = colX0[krowX] ; 
            sum00 += A0k * Xk0 ;
            sum10 += A1k * Xk0 ;
         }
      }
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
         colY0[irowAT+1] -= sum10 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
         irowY = rowindAT[irowAT+1] ;
         colY0[irowY] -= sum10 ;
      }
   } else if ( irowAT == nrowAT - 1 ) {
      double   A0k, Xk0 ;
      double   sum00 ;

      sum00 = 0.0 ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            Xk0 = colX0[kcolAT] ;
            sum00 += A0k * Xk0 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            A0k = rowAT0[kcolAT] ; 
            krowX = colindAT[kcolAT] ;
            Xk0 = colX0[krowX] ; 
            sum00 += A0k * Xk0 ;
         }
      }
      if ( nrowY == nrowAT ) {
         colY0[irowAT]   -= sum00 ;
      } else {
         irowY = rowindAT[irowAT] ;
         colY0[irowY] -= sum00 ;
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
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   Aik, sum0, sum1, sum2 ;
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, inc1, inc2, irowAT, irowY, jcolX, kk, krowX, 
         ncolAT, ncolX, ncolY, nentA, nrowAT, nrowX, nrowY, size ;
int      *colindAT, *indices, *rowindAT, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_ROWS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseColumnsInfo(mtxA, &nrowAT, &nentA, &sizes, &indices, &entA) ;
if ( (ncolAT = mtxA->nrow) != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( (nrowAT = mtxA->ncol) != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + nrowX ;
   colX2 = colX1 + nrowX ;
   colY1 = colY0 + nrowY ;
   colY2 = colY1 + nrowY ;
   for ( irowAT = kk = 0 ; irowAT < nrowAT ; irowAT++ ) {
      if ( (size = sizes[irowAT]) > 0 ) {
         sum0 = sum1 = sum2 = 0.0 ;
         if ( ncolAT == nrowX ) {
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
               krowX = colindAT[indices[kk]] ;
               sum0 += Aik * colX0[krowX] ;
               sum1 += Aik * colX1[krowX] ;
               sum2 += Aik * colX2[krowX] ;
            }
         }
         if ( nrowAT == nrowY ) {
            colY0[irowAT] -= sum0 ;
            colY1[irowAT] -= sum1 ;
            colY2[irowAT] -= sum2 ;
         } else {
            irowY = rowindAT[irowAT] ;
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
   for ( irowAT = kk = 0 ; irowAT < nrowAT ; irowAT++ ) {
      if ( (size = sizes[irowAT]) > 0 ) {
         sum0 = sum1 = 0.0 ;
         if ( ncolAT == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = indices[kk] ;
               sum0 += Aik * colX0[krowX] ;
               sum1 += Aik * colX1[krowX] ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = colindAT[indices[kk]] ;
               sum0 += Aik * colX0[krowX] ;
               sum1 += Aik * colX1[krowX] ;
            }
         }
         if ( nrowAT == nrowY ) {
            colY0[irowAT] -= sum0 ;
            colY1[irowAT] -= sum1 ;
         } else {
            irowY = rowindAT[irowAT] ;
            colY0[irowY] -= sum0 ;
            colY1[irowY] -= sum1 ;
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   for ( irowAT = kk = 0 ; irowAT < nrowAT ; irowAT++ ) {
      if ( (size = sizes[irowAT]) > 0 ) {
         sum0 = 0.0 ;
         if ( ncolAT == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = indices[kk] ;
               sum0 += Aik * colX0[krowX] ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aik = entA[kk] ;
               krowX = colindAT[indices[kk]] ;
               sum0 += Aik * colX0[krowX] ;
            }
         }
         if ( nrowAT == nrowY ) {
            colY0[irowAT] -= sum0 ;
         } else {
            irowY = rowindAT[irowAT] ;
            colY0[irowY] -= sum0 ;
         }
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
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   Aij, Xj0, Xj1, Xj2 ;
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, inc1, inc2, irowY, jcolAT, jcolX, jrowX, kk, 
         ncolAT, ncolX, ncolY, nentA, nrowAT, nrowX, nrowY, size ;
int      *colindAT, *indices, *rowindAT, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_COLUMNS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseRowsInfo(mtxA, &ncolAT, &nentA, &sizes, &indices, &entA) ;
if ( (ncolAT = mtxA->nrow) != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( (nrowAT = mtxA->ncol) != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + nrowX ;
   colX2 = colX1 + nrowX ;
   colY1 = colY0 + nrowY ;
   colY2 = colY1 + nrowY ;
   for ( jcolAT = kk = 0 ; jcolAT < ncolAT ; jcolAT++ ) {
      if ( (size = sizes[jcolAT]) > 0 ) {
         if ( ncolAT == nrowX ) {
            jrowX = jcolAT ;
         } else {
            jrowX = colindAT[jcolAT] ;
         }
         Xj0 = colX0[jrowX] ;
         Xj1 = colX1[jrowX] ;
         Xj2 = colX2[jrowX] ;
         if ( nrowAT == nrowY ) {
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
               irowY = rowindAT[indices[kk]] ;
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
   for ( jcolAT = kk = 0 ; jcolAT < ncolAT ; jcolAT++ ) {
      if ( (size = sizes[jcolAT]) > 0 ) {
         if ( ncolAT == nrowX ) {
            jrowX = jcolAT ;
         } else {
            jrowX = colindAT[jcolAT] ;
         }
         Xj0 = colX0[jrowX] ;
         Xj1 = colX1[jrowX] ;
         if ( nrowAT == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = indices[kk] ;
               colY0[irowY] -= Aij * Xj0 ;
               colY1[irowY] -= Aij * Xj1 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = rowindAT[indices[kk]] ;
               colY0[irowY] -= Aij * Xj0 ;
               colY1[irowY] -= Aij * Xj1 ;
            }
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   for ( jcolAT = kk = 0 ; jcolAT < ncolAT ; jcolAT++ ) {
      if ( (size = sizes[jcolAT]) > 0 ) {
         if ( ncolAT == nrowX ) {
            jrowX = jcolAT ;
         } else {
            jrowX = colindAT[jcolAT] ;
         }
         Xj0 = colX0[jrowX] ;
         if ( nrowAT == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = indices[kk] ;
               colY0[irowY] -= Aij * Xj0 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               Aij = entA[kk] ;
               irowY = rowindAT[indices[kk]] ;
               colY0[irowY] -= Aij * Xj0 ;
            }
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
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   ai0, ai1, ai2, ar0, ar1, ar2,
         xi00, xi01, xi02, xi10, xi11, xi12, xi20, xi21, xi22,
         xr00, xr01, xr02, xr10, xr11, xr12, xr20, xr21, xr22 ;
double   *colAT0, *colAT1, *colAT2, *colX0, *colX1, *colX2, 
         *colY0, *colY1, *colY2, *entA, *entX, *entY ;
int      icolAT, ialoc, iloc, inc1, inc2, jcolX, krowAT, 
         ncolAT, ncolX, ncolY, nrowAT, nrowX, nrowY, raloc, rloc ;
int      *colindAT, *rowindAT ;

SubMtx_denseInfo(mtxY, &nrowY,  &ncolY,  &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX,  &ncolX,  &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &ncolAT, &nrowAT, &inc1, &inc2, &entA) ;
colX0 = entX ;
colY0 = entY ;
if ( ncolAT != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( nrowAT != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1  = colX0 + 2*nrowX ;
   colX2  = colX1 + 2*nrowX ;
   colY1  = colY0 + 2*nrowY ;
   colY2  = colY1 + 2*nrowY ;
   colAT0 = entA ;
   for ( icolAT = 0 ; icolAT < ncolAT - 2 ; icolAT += 3 ) {
      colAT1 = colAT0 + 2*nrowAT ;
      colAT2 = colAT1 + 2*nrowAT ;
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            ar1 = colAT1[rloc] ; ai1 = colAT1[iloc] ;
            ar2 = colAT2[rloc] ; ai2 = colAT2[iloc] ;
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
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            ar1 = colAT1[raloc] ; ai1 = colAT1[ialoc] ;
            ar2 = colAT2[raloc] ; ai2 = colAT2[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
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
      }
      colAT0 = colAT2 + 2*nrowAT ;
   }
   if ( icolAT == ncolAT - 2 ) {
      colAT1 = colAT0 + 2*nrowAT ;
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            ar1 = colAT1[rloc] ; ai1 = colAT1[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
            colY2[rloc] -= ar0*xr02 - ai0*xi02 + ar1*xr12 - ai1*xi12 ;
            colY2[iloc] -= ar0*xi02 + ai0*xr02 + ar1*xi12 + ai1*xr12 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            ar1 = colAT1[raloc] ; ai1 = colAT1[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
            colY2[rloc] -= ar0*xr02 - ai0*xi02 + ar1*xr12 - ai1*xi12 ;
            colY2[iloc] -= ar0*xi02 + ai0*xr02 + ar1*xi12 + ai1*xr12 ;
         }
      }
   } else if ( icolAT == ncolAT - 1 ) {
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         xr02 = colX2[rloc] ; xi02 = colX2[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         xr12 = colX2[rloc] ; xi12 = colX2[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
         xr22 = colX2[rloc] ; xi22 = colX2[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 ;
            colY2[rloc] -= ar0*xr02 - ai0*xi02 ;
            colY2[iloc] -= ar0*xi02 + ai0*xr02 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 ;
            colY2[rloc] -= ar0*xr02 - ai0*xi02 ;
            colY2[iloc] -= ar0*xi02 + ai0*xr02 ;
         }
      }
   }
   colX0 = colX2 + 2*nrowX ;
   colY0 = colY2 + 2*nrowY ;
}
/*
fprintf(stdout, "\n %% jcolX = %d", jcolX) ;
*/
if ( jcolX == ncolX - 2 ) {
   colX1  = colX0 + 2*nrowX ;
   colY1  = colY0 + 2*nrowY ;
   colAT0 = entA ;
   for ( icolAT = 0 ; icolAT < ncolAT - 2 ; icolAT += 3 ) {
      colAT1 = colAT0 + 2*nrowAT ;
      colAT2 = colAT1 + 2*nrowAT ;
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            ar1 = colAT1[rloc] ; ai1 = colAT1[iloc] ;
            ar2 = colAT2[rloc] ; ai2 = colAT2[iloc] ;
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
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            ar1 = colAT1[raloc] ; ai1 = colAT1[ialoc] ;
            ar2 = colAT2[raloc] ; ai2 = colAT2[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
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
      }
      colAT0 = colAT2 + 2*nrowAT ;
   }
   if ( icolAT == ncolAT - 2 ) {
      colAT1 = colAT0 + 2*nrowAT ;
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            ar1 = colAT1[rloc] ; ai1 = colAT1[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            ar1 = colAT1[raloc] ; ai1 = colAT1[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 + ar1*xr11 - ai1*xi11 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 + ar1*xi11 + ai1*xr11 ;
         }
      }
   } else if ( icolAT == ncolAT - 1 ) {
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         xr01 = colX1[rloc] ; xi01 = colX1[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         xr11 = colX1[rloc] ; xi11 = colX1[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
         xr21 = colX1[rloc] ; xi21 = colX1[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
            colY1[rloc] -= ar0*xr01 - ai0*xi01 ;
            colY1[iloc] -= ar0*xi01 + ai0*xr01 ;
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   colAT0 = entA ;
   for ( icolAT = 0 ; icolAT < ncolAT - 2 ; icolAT += 3 ) {
/*
fprintf(stdout, "\n %% icolAT = %d", icolAT) ;
*/
      colAT1 = colAT0 + 2*nrowAT ;
      colAT2 = colAT1 + 2*nrowAT ;
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      }
/*
fprintf(stdout, "\n %% x00 = (%12.4e,%12.4e)", xr00, xi00) ;
fprintf(stdout, "\n %% x10 = (%12.4e,%12.4e)", xr10, xi10) ;
fprintf(stdout, "\n %% x20 = (%12.4e,%12.4e)", xr20, xi20) ;
*/
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            ar1 = colAT1[rloc] ; ai1 = colAT1[iloc] ;
            ar2 = colAT2[rloc] ; ai2 = colAT2[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00
                         + ar1*xr10 - ai1*xi10
                         + ar2*xr20 - ai2*xi20 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00
                         + ar1*xi10 + ai1*xr10
                         + ar2*xi20 + ai2*xr20 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            ar1 = colAT1[raloc] ; ai1 = colAT1[ialoc] ;
            ar2 = colAT2[raloc] ; ai2 = colAT2[ialoc] ;
/*
fprintf(stdout, "\n %% raloc = %d, ialoc = %d", raloc, ialoc) ;
fprintf(stdout, "\n %% a0 = (%12.4e,%12.4e)", ar0, ai0) ;
fprintf(stdout, "\n %% a1 = (%12.4e,%12.4e)", ar1, ai1) ;
fprintf(stdout, "\n %% a2 = (%12.4e,%12.4e)", ar2, ai2) ;
*/
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
/*
fprintf(stdout, "\n %% rloc = %d, iloc = %d", rloc, iloc) ;
*/
            colY0[rloc] -= ar0*xr00 - ai0*xi00
                         + ar1*xr10 - ai1*xi10
                         + ar2*xr20 - ai2*xi20 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00
                         + ar1*xi10 + ai1*xr10
                         + ar2*xi20 + ai2*xr20 ;
         }
      }
      colAT0 = colAT2 + 2*nrowAT ;
   }
   if ( icolAT == ncolAT - 2 ) {
      colAT1 = colAT0 + 2*nrowAT ;
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            ar1 = colAT1[rloc] ; ai1 = colAT1[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            ar1 = colAT1[raloc] ; ai1 = colAT1[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 + ar1*xr10 - ai1*xi10 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 + ar1*xi10 + ai1*xr10 ;
         }
      }
   } else if ( icolAT == ncolAT - 1 ) {
      if ( ncolAT == nrowX ) {
         rloc = 2*icolAT ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc += 2 ; iloc += 2 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      } else {
         rloc = 2*colindAT[icolAT] ; iloc = rloc + 1 ;
         xr00 = colX0[rloc] ; xi00 = colX0[iloc] ;
         rloc = 2*colindAT[icolAT+1] ; iloc = rloc + 1 ;
         xr10 = colX0[rloc] ; xi10 = colX0[iloc] ;
         rloc = 2*colindAT[icolAT+2] ; iloc = rloc + 1 ;
         xr20 = colX0[rloc] ; xi20 = colX0[iloc] ;
      }
      if ( nrowY == nrowAT ) {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            rloc = 2*krowAT ; iloc = rloc + 1 ;
            ar0 = colAT0[rloc] ; ai0 = colAT0[iloc] ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
         }
      } else {
         for ( krowAT = 0 ; krowAT < nrowAT ; krowAT++ ) {
            raloc = 2*krowAT ; ialoc = raloc + 1 ;
            ar0 = colAT0[raloc] ; ai0 = colAT0[ialoc] ;
            rloc = 2*rowindAT[krowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= ar0*xr00 - ai0*xi00 ;
            colY0[iloc] -= ar0*xi00 + ai0*xr00 ;
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
complex_updDenseColumns (
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2, 
         *rowAT0, *rowAT1, *rowAT2, *entA, *entX, *entY ;
int      inc1, inc2, irowAT, jcolX, kcolAT, 
         ncolAT, ncolX, ncolY, nrowAT, nrowX, nrowY ;
int      *colindAT, *rowindAT ;

SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_denseInfo(mtxA, &ncolAT, &nrowAT, &inc1, &inc2, &entA) ;
if ( ncolAT != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( nrowAT != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   colX1 = colX0 + 2*nrowX ;
   colX2 = colX1 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colY2 = colY1 + 2*nrowY ;
   rowAT0 = entA ;
   for ( irowAT = 0 ; irowAT < nrowAT - 2 ; irowAT += 3 ) {
      double   ai0, ai1, ai2, ar0, ar1, ar2, isum00, isum01, isum02, 
               isum10, isum11, isum12, isum20, isum21, isum22, 
               rsum00, rsum01, rsum02, rsum10, rsum11, rsum12, 
               rsum20, rsum21, rsum22, xi0, xi1, xi2, xr0, xr1, xr2 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum01 = isum02 = isum10 = isum11 = isum12 
             = isum20 = isum21 = isum22 = 0.0 ;
      rsum00 = rsum01 = rsum02 = rsum10 = rsum11 = rsum12 
             = rsum20 = rsum21 = rsum22 = 0.0 ;
      rowAT1 = rowAT0 + 2*ncolAT ;
      rowAT2 = rowAT1 + 2*ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            ar1 = rowAT1[rloc] ; ai1 = rowAT1[iloc] ;
            ar2 = rowAT2[rloc] ; ai2 = rowAT2[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc]  ; xi1 = colX1[iloc] ;
            xr2 = colX2[rloc]  ; xi2 = colX2[iloc] ;
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
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            ar1 = rowAT1[raloc] ; ai1 = rowAT1[ialoc] ;
            ar2 = rowAT2[raloc] ; ai2 = rowAT2[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc]  ; xi1 = colX1[ixloc] ;
            xr2 = colX2[rxloc]  ; xi2 = colX2[ixloc] ;
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
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         colY2[rloc] -= rsum02 ; colY2[iloc] -= isum02 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
         colY2[rloc] -= rsum12 ; colY2[iloc] -= isum12 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum20 ; colY0[iloc] -= isum20 ;
         colY1[rloc] -= rsum21 ; colY1[iloc] -= isum21 ;
         colY2[rloc] -= rsum22 ; colY2[iloc] -= isum22 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         colY2[rloc] -= rsum02 ; colY2[iloc] -= isum02 ;
         rloc = 2*rowindAT[irowAT+1] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
         colY2[rloc] -= rsum12 ; colY2[iloc] -= isum12 ;
         rloc = 2*rowindAT[irowAT+2] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum20 ; colY0[iloc] -= isum20 ;
         colY1[rloc] -= rsum21 ; colY1[iloc] -= isum21 ;
         colY2[rloc] -= rsum22 ; colY2[iloc] -= isum22 ;
      }
      rowAT0 = rowAT2 + 2*ncolAT ;
   }
   if ( irowAT == nrowAT - 2 ) {
      double   ai0, ai1, ar0, ar1, isum00, isum01, isum02, 
               isum10, isum11, isum12, rsum00, rsum01, rsum02, 
               rsum10, rsum11, rsum12, xi0, xi1, xi2, xr0, xr1, xr2 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum01 = isum02 = isum10 = isum11 = isum12 = 0.0 ;
      rsum00 = rsum01 = rsum02 = rsum10 = rsum11 = rsum12 = 0.0 ;
      rowAT1 = rowAT0 + 2*ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            ar1 = rowAT1[rloc] ; ai1 = rowAT1[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc]  ; xi1 = colX1[iloc] ;
            xr2 = colX2[rloc]  ; xi2 = colX2[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum12 += ar1*xr2 - ai1*xi2 ; isum12 += ar1*xi2 + ai1*xr2 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            ar1 = rowAT1[raloc] ; ai1 = rowAT1[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc]  ; xi1 = colX1[ixloc] ;
            xr2 = colX2[rxloc]  ; xi2 = colX2[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum12 += ar1*xr2 - ai1*xi2 ; isum12 += ar1*xi2 + ai1*xr2 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         colY2[rloc] -= rsum02 ; colY2[iloc] -= isum02 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
         colY2[rloc] -= rsum12 ; colY2[iloc] -= isum12 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         colY2[rloc] -= rsum02 ; colY2[iloc] -= isum02 ;
         rloc = 2*rowindAT[irowAT+1] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
         colY2[rloc] -= rsum12 ; colY2[iloc] -= isum12 ;
      }
   } else if ( irowAT == nrowAT - 1 ) {
      double   ai0, ar0, isum00, isum01, isum02, 
               rsum00, rsum01, rsum02, xi0, xi1, xi2, xr0, xr1, xr2 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum01 = isum02 = 0.0 ;
      rsum00 = rsum01 = rsum02 = 0.0 ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc]  ; xi1 = colX1[iloc] ;
            xr2 = colX2[rloc]  ; xi2 = colX2[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc]  ; xi1 = colX1[ixloc] ;
            xr2 = colX2[rxloc]  ; xi2 = colX2[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum02 += ar0*xr2 - ai0*xi2 ; isum02 += ar0*xi2 + ai0*xr2 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         colY2[rloc] -= rsum02 ; colY2[iloc] -= isum02 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         colY2[rloc] -= rsum02 ; colY2[iloc] -= isum02 ;
      }
   }
   colX0 = colX2 + 2*nrowX ;
   colY0 = colY2 + 2*nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   colX1 = colX0 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   rowAT0 = entA ;
   for ( irowAT = 0 ; irowAT < nrowAT - 2 ; irowAT += 3 ) {
      double   ai0, ai1, ai2, ar0, ar1, ar2, isum00, isum01, 
               isum10, isum11, isum20, isum21, rsum00, rsum01, rsum10, 
               rsum11, rsum20, rsum21, xi0, xi1, xr0, xr1 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum01 = isum10 = isum11 = isum20 = isum21 = 0.0 ;
      rsum00 = rsum01 = rsum10 = rsum11 = rsum20 = rsum21 = 0.0 ;
      rowAT1 = rowAT0 + 2*ncolAT ;
      rowAT2 = rowAT1 + 2*ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            ar1 = rowAT1[rloc] ; ai1 = rowAT1[iloc] ;
            ar2 = rowAT2[rloc] ; ai2 = rowAT2[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc]  ; xi1 = colX1[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
            rsum21 += ar2*xr1 - ai2*xi1 ; isum21 += ar2*xi1 + ai2*xr1 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            ar1 = rowAT1[raloc] ; ai1 = rowAT1[ialoc] ;
            ar2 = rowAT2[raloc] ; ai2 = rowAT2[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc]  ; xi1 = colX1[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
            rsum21 += ar2*xr1 - ai2*xi1 ; isum21 += ar2*xi1 + ai2*xr1 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum20 ; colY0[iloc] -= isum20 ;
         colY1[rloc] -= rsum21 ; colY1[iloc] -= isum21 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         rloc = 2*rowindAT[irowAT+1] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
         rloc = 2*rowindAT[irowAT+2] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum20 ; colY0[iloc] -= isum20 ;
         colY1[rloc] -= rsum21 ; colY1[iloc] -= isum21 ;
      }
      rowAT0 = rowAT2 + 2*ncolAT ;
   }
   if ( irowAT == nrowAT - 2 ) {
      double   ai0, ai1, ar0, ar1, isum00, isum01, isum10, isum11, 
               rsum00, rsum01, rsum10, rsum11, xi0, xi1, xr0, xr1 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum01 = isum10 = isum11 = 0.0 ;
      rsum00 = rsum01 = rsum10 = rsum11 = 0.0 ;
      rowAT1 = rowAT0 + 2*ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            ar1 = rowAT1[rloc] ; ai1 = rowAT1[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc]  ; xi1 = colX1[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            ar1 = rowAT1[raloc] ; ai1 = rowAT1[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc]  ; xi1 = colX1[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum11 += ar1*xr1 - ai1*xi1 ; isum11 += ar1*xi1 + ai1*xr1 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
         rloc = 2*rowindAT[irowAT+1] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         colY1[rloc] -= rsum11 ; colY1[iloc] -= isum11 ;
      }
   } else if ( irowAT == nrowAT - 1 ) {
      double   ai0, ar0, isum00, isum01, 
               rsum00, rsum01, xi0, xi1, xr0, xr1 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum01 = 0.0 ;
      rsum00 = rsum01 = 0.0 ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            xr1 = colX1[rloc]  ; xi1 = colX1[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            xr1 = colX1[rxloc]  ; xi1 = colX1[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum01 += ar0*xr1 - ai0*xi1 ; isum01 += ar0*xi1 + ai0*xr1 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         colY1[rloc] -= rsum01 ; colY1[iloc] -= isum01 ;
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   rowAT0 = entA ;
   for ( irowAT = 0 ; irowAT < nrowAT - 2 ; irowAT += 3 ) {
      double   ai0, ai1, ai2, ar0, ar1, ar2, isum00, isum10, isum20, 
               rsum00, rsum10, rsum20, xi0, xr0 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum10 = isum20 = 0.0 ;
      rsum00 = rsum10 = rsum20 = 0.0 ;
      rowAT1 = rowAT0 + 2*ncolAT ;
      rowAT2 = rowAT1 + 2*ncolAT ;
/*
fprintf(stdout, "\n %% irowAT %d", irowAT) ;
*/
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
/*
fprintf(stdout, "\n %%    rloc %d, iloc %d", rloc, iloc) ;
*/
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            ar1 = rowAT1[rloc] ; ai1 = rowAT1[iloc] ;
            ar2 = rowAT2[rloc] ; ai2 = rowAT2[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
/*
fprintf(stdout, "\n %%    a0 = (%12.4e,%12.4e)", ar0, ai0) ;
fprintf(stdout, "\n %%    a1 = (%12.4e,%12.4e)", ar1, ai1) ;
fprintf(stdout, "\n %%    a2 = (%12.4e,%12.4e)", ar2, ai2) ;
fprintf(stdout, "\n %%    x0 = (%12.4e,%12.4e)", xr0, xi0) ;
*/
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            ar1 = rowAT1[raloc] ; ai1 = rowAT1[ialoc] ;
            ar2 = rowAT2[raloc] ; ai2 = rowAT2[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
            rsum20 += ar2*xr0 - ai2*xi0 ; isum20 += ar2*xi0 + ai2*xr0 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum20 ; colY0[iloc] -= isum20 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         rloc = 2*rowindAT[irowAT+1] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
         rloc = 2*rowindAT[irowAT+2] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum20 ; colY0[iloc] -= isum20 ;
      }
      rowAT0 = rowAT2 + 2*ncolAT ;
   }
   if ( irowAT == nrowAT - 2 ) {
      double   ai0, ai1, ar0, ar1, isum00, isum10, rsum00, rsum10, 
               xi0, xr0 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = isum10 = 0.0 ;
      rsum00 = rsum10 = 0.0 ;
      rowAT1 = rowAT0 + 2*ncolAT ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            ar1 = rowAT1[rloc] ; ai1 = rowAT1[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            ar1 = rowAT1[raloc] ; ai1 = rowAT1[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
            rsum10 += ar1*xr0 - ai1*xi0 ; isum10 += ar1*xi0 + ai1*xr0 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         rloc+= 2 ; iloc += 2 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
         rloc = 2*rowindAT[irowAT+1] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum10 ; colY0[iloc] -= isum10 ;
      }
   } else if ( irowAT == nrowAT - 1 ) {
      double   ai0, ar0, isum00, rsum00, xi0, xr0 ;
      int      ialoc, iloc, ixloc, raloc, rloc, rxloc ;

      isum00 = 0.0 ;
      rsum00 = 0.0 ;
      if ( ncolAT == nrowX ) {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            rloc = 2*kcolAT ; iloc = rloc + 1 ;
            ar0 = rowAT0[rloc] ; ai0 = rowAT0[iloc] ;
            xr0 = colX0[rloc]  ; xi0 = colX0[iloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
         }
      } else {
         for ( kcolAT = 0 ; kcolAT < ncolAT ; kcolAT++ ) {
            raloc = 2*kcolAT ; ialoc = raloc + 1 ;
            ar0 = rowAT0[raloc] ; ai0 = rowAT0[ialoc] ;
            rxloc = 2*colindAT[kcolAT] ; ixloc = rxloc + 1 ;
            xr0 = colX0[rxloc]  ; xi0 = colX0[ixloc] ;
            rsum00 += ar0*xr0 - ai0*xi0 ; isum00 += ar0*xi0 + ai0*xr0 ;
         }
      }
      if ( nrowY == nrowAT ) {
         rloc = 2*irowAT ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
      } else {
         rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
         colY0[rloc] -= rsum00 ; colY0[iloc] -= isum00 ;
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
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, iloc, inc1, inc2, irowAT, jcolX, kk, ncolAT, ncolX, 
         ncolY, nentA, nrowAT, nrowX, nrowY, rloc, size ;
int      *colindAT, *indices, *rowindAT, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_ROWS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseColumnsInfo(mtxA, &nrowAT, &nentA, &sizes, &indices, &entA) ;
if ( (ncolAT = mtxA->nrow) != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( (nrowAT = mtxA->ncol) != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   double   ai, ar, isum0, isum1, isum2, rsum0, rsum1, rsum2, 
            xi0, xi1, xi2, xr0, xr1, xr2 ;

   colX1 = colX0 + 2*nrowX ;
   colX2 = colX1 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colY2 = colY1 + 2*nrowY ;
   for ( irowAT = kk = 0 ; irowAT < nrowAT ; irowAT++ ) {
      if ( (size = sizes[irowAT]) > 0 ) {
         isum0 = isum1 = isum2 = rsum0 = rsum1 = rsum2 = 0.0 ;
         if ( ncolAT == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*indices[kk] ; iloc = rloc + 1 ;
               xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
               xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
               xr2 = colX2[rloc] ; xi2 = colX2[iloc] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
               rsum2 += ar*xr2 - ai*xi2 ; isum2 += ar*xi2 + ai*xr2 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*colindAT[indices[kk]] ; iloc = rloc + 1 ;
               xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
               xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
               xr2 = colX2[rloc] ; xi2 = colX2[iloc] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
               rsum2 += ar*xr2 - ai*xi2 ; isum2 += ar*xi2 + ai*xr2 ;
            }
         }
         if ( nrowAT == nrowY ) {
            rloc = 2*irowAT ; iloc = rloc + 1 ;
            colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
            colY1[rloc] -= rsum1 ; colY1[iloc] -= isum1 ;
            colY2[rloc] -= rsum2 ; colY2[iloc] -= isum2 ;
         } else {
            rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
            colY1[rloc] -= rsum1 ; colY1[iloc] -= isum1 ;
            colY2[rloc] -= rsum2 ; colY2[iloc] -= isum2 ;
         }
      }
   }
   colX0 = colX2 + 2*nrowX ;
   colY0 = colY2 + 2*nrowY ;
}
if ( jcolX == ncolX - 2 ) {
   double   ai, ar, isum0, isum1, rsum0, rsum1, xi0, xi1, xr0, xr1 ;

   colX1 = colX0 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   for ( irowAT = kk = 0 ; irowAT < nrowAT ; irowAT++ ) {
      if ( (size = sizes[irowAT]) > 0 ) {
         isum0 = isum1 = rsum0 = rsum1 = 0.0 ;
         if ( ncolAT == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*indices[kk] ; iloc = rloc + 1 ;
               xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
               xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*colindAT[indices[kk]] ; iloc = rloc + 1 ;
               xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
               xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
               rsum1 += ar*xr1 - ai*xi1 ; isum1 += ar*xi1 + ai*xr1 ;
            }
         }
         if ( nrowAT == nrowY ) {
            rloc = 2*irowAT ; iloc = rloc + 1 ;
            colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
            colY1[rloc] -= rsum1 ; colY1[iloc] -= isum1 ;
         } else {
            rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
            colY1[rloc] -= rsum1 ; colY1[iloc] -= isum1 ;
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   double   ai, ar, isum0, rsum0, xi0, xr0 ;

   for ( irowAT = kk = 0 ; irowAT < nrowAT ; irowAT++ ) {
      if ( (size = sizes[irowAT]) > 0 ) {
         isum0 = rsum0 = 0.0 ;
         if ( ncolAT == nrowX ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*indices[kk] ; iloc = rloc + 1 ;
               xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*colindAT[indices[kk]] ; iloc = rloc + 1 ;
               xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
               rsum0 += ar*xr0 - ai*xi0 ; isum0 += ar*xi0 + ai*xr0 ;
            }
         }
         if ( nrowAT == nrowY ) {
            rloc = 2*irowAT ; iloc = rloc + 1 ;
            colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
         } else {
            rloc = 2*rowindAT[irowAT] ; iloc = rloc + 1 ;
            colY0[rloc] -= rsum0 ; colY0[iloc] -= isum0 ;
         }
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
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) {
double   *colX0, *colX1, *colX2, *colY0, *colY1, *colY2,
         *entA, *entX, *entY ;
int      ii, inc1, inc2, jcolAT, jcolX, jrowX, kk, 
         ncolAT, ncolX, ncolY, nentA, nrowAT, nrowX, nrowY, size ;
int      *colindAT, *indices, *rowindAT, *sizes ;
/*
fprintf(stdout, "\n UPDATE_SPARSE_COLUMNS(%d,%d)", 
        mtxA->rowid, mtxA->colid) ;
*/
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
SubMtx_sparseRowsInfo(mtxA, &ncolAT, &nentA, &sizes, &indices, &entA) ;
if ( (ncolAT = mtxA->nrow) != nrowX ) {
   SubMtx_rowIndices(mtxA, &ncolAT, &colindAT) ;
} else {
   colindAT = NULL ;
}
if ( (nrowAT = mtxA->ncol) != nrowY ) {
   SubMtx_columnIndices(mtxA, &nrowAT, &rowindAT) ;
} else {
   rowindAT = NULL ;
}
colX0 = entX ;
colY0 = entY ;
for ( jcolX = 0 ; jcolX < ncolX - 2 ; jcolX += 3 ) {
   double   ai, ar, xi0, xi1, xi2, xr0, xr1, xr2 ;
   int      iloc, rloc ;

   colX1 = colX0 + 2*nrowX ;
   colX2 = colX1 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   colY2 = colY1 + 2*nrowY ;
   for ( jcolAT = kk = 0 ; jcolAT < ncolAT ; jcolAT++ ) {
      if ( (size = sizes[jcolAT]) > 0 ) {
         if ( ncolAT == nrowX ) {
            jrowX = jcolAT ;
         } else {
            jrowX = colindAT[jcolAT] ;
         }
         rloc = 2*jrowX ; iloc = rloc + 1 ;
         xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
         xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
         xr2 = colX2[rloc] ; xi2 = colX2[iloc] ;
         if ( nrowAT == nrowY ) {
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
               rloc = 2*rowindAT[indices[kk]] ; iloc = rloc + 1 ;
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
   double   ai, ar, xi0, xi1, xr0, xr1 ;
   int      iloc, rloc ;

   colX1 = colX0 + 2*nrowX ;
   colY1 = colY0 + 2*nrowY ;
   for ( jcolAT = kk = 0 ; jcolAT < ncolAT ; jcolAT++ ) {
      if ( (size = sizes[jcolAT]) > 0 ) {
         if ( ncolAT == nrowX ) {
            jrowX = jcolAT ;
         } else {
            jrowX = colindAT[jcolAT] ;
         }
         rloc = 2*jrowX ; iloc = rloc + 1 ;
         xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
         xr1 = colX1[rloc] ; xi1 = colX1[iloc] ;
         if ( nrowAT == nrowY ) {
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
               rloc = 2*rowindAT[indices[kk]] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
               colY1[rloc] -= ar*xr1 - ai*xi1 ;
               colY1[iloc] -= ar*xi1 + ai*xr1 ;
            }
         }
      }
   }
} else if ( jcolX == ncolX - 1 ) {
   double   ai, ar, xi0, xr0 ;
   int      iloc, rloc ;

   for ( jcolAT = kk = 0 ; jcolAT < ncolAT ; jcolAT++ ) {
      if ( (size = sizes[jcolAT]) > 0 ) {
         if ( ncolAT == nrowX ) {
            jrowX = jcolAT ;
         } else {
            jrowX = colindAT[jcolAT] ;
         }
         rloc = 2*jrowX ; iloc = rloc + 1 ;
         xr0 = colX0[rloc] ; xi0 = colX0[iloc] ;
         if ( nrowAT == nrowY ) {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*indices[kk] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
            }
         } else {
            for ( ii = 0 ; ii < size ; ii++, kk++ ) {
               ar = entA[2*kk] ; ai = entA[2*kk+1] ;
               rloc = 2*rowindAT[indices[kk]] ; iloc = rloc + 1 ;
               colY0[rloc] -= ar*xr0 - ai*xi0 ;
               colY0[iloc] -= ar*xi0 + ai*xr0 ;
            }
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
