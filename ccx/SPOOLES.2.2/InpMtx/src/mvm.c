/*  mvm.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X

   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_nonsym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || Y == NULL || alpha == NULL || X == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
           "\n bad input\n", A, Y, alpha, X) ;
   exit(-1) ;
}
if ( ! (INPMTX_IS_REAL_ENTRIES(A) || INPMTX_IS_COMPLEX_ENTRIES(A)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, Y, alpha, X, A->inputMode) ;
   exit(-1) ;
}
if ( ! (DENSEMTX_IS_REAL(Y) || DENSEMTX_IS_COMPLEX(Y)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n bad type %d for Y\n", A, Y, alpha, X, Y->type) ;
   exit(-1) ;
}
if ( ! (DENSEMTX_IS_REAL(X) || DENSEMTX_IS_COMPLEX(X)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n bad type %d for X\n", A, Y, alpha, X, X->type) ;
   exit(-1) ;
}
if ( DENSEMTX_IS_REAL(Y) != DENSEMTX_IS_REAL(X) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n X's type %d, Y's type %d \n", A, Y, alpha, X, 
          X->type, Y->type) ;
   exit(-1) ;
}
if (  (INPMTX_IS_REAL_ENTRIES(A) && !DENSEMTX_IS_REAL(Y))
   || (INPMTX_IS_COMPLEX_ENTRIES(A) && !DENSEMTX_IS_COMPLEX(Y)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n A's inputMode %d, Y's type %d \n", A, Y, alpha, X, 
          A->inputMode, Y->type) ;
   exit(-1) ;
}
if ( DenseMtx_rowIncrement(Y) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n Y's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(Y)) ;
   exit(-1) ;
}
incY = DenseMtx_columnIncrement(Y) ;
if ( DenseMtx_rowIncrement(X) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n X's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(X)) ;
   exit(-1) ;
}
incX = DenseMtx_columnIncrement(X) ;
x    = DenseMtx_entries(X) ;
y    = DenseMtx_entries(Y) ;
DenseMtx_dimensions(Y, &nrowY, &ncolY) ;
DenseMtx_dimensions(X, &nrowX, &ncolX) ;
if ( (nrhs = ncolY) != ncolX ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
          "\n Y's nrhs = %d, X's nrhs = %d",
          A, Y, alpha, X, nrhs, ncolX) ;
   exit(-1) ;
}
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A) ;
if ( ivec1 == NULL || ivec2 == NULL || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, Y, alpha, X, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
nent  = A->nent ;
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   double   rfac ;
   int      chev, col, ii, jrhs, off, row ;

   rfac = alpha[0] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               row = ivec1[ii] ; col = ivec2[ii] ;
               y[row] += dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               row = ivec1[ii] ; col = ivec2[ii] ;
               y[row] += rfac * dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               col = ivec1[ii] ; row = ivec2[ii] ;
               y[row] += dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               col = ivec1[ii] ; row = ivec2[ii] ;
               y[row] += rfac * dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               y[row] += dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               y[row] += rfac * dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      }
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, jrhs, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row   = ivec1[ii] ; col   = ivec2[ii] ;
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec1[ii]  ; col = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec1[ii]  ; col = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^T*X

   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmm_T (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || Y == NULL || alpha == NULL || X == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
           "\n bad input\n", A, Y, alpha, X) ;
   exit(-1) ;
}
if ( ! (INPMTX_IS_REAL_ENTRIES(A) || INPMTX_IS_COMPLEX_ENTRIES(A)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, Y, alpha, X, A->inputMode) ;
   exit(-1) ;
}
if ( ! (DENSEMTX_IS_REAL(Y) || DENSEMTX_IS_COMPLEX(Y)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n bad type %d for Y\n", A, Y, alpha, X, Y->type) ;
   exit(-1) ;
}
if ( ! (DENSEMTX_IS_REAL(X) || DENSEMTX_IS_COMPLEX(X)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n bad type %d for X\n", A, Y, alpha, X, X->type) ;
   exit(-1) ;
}
if ( DENSEMTX_IS_REAL(Y) != DENSEMTX_IS_REAL(X) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n X's type %d, Y's type %d \n", A, Y, alpha, X, 
          X->type, Y->type) ;
   exit(-1) ;
}
if (  (INPMTX_IS_REAL_ENTRIES(A) && !DENSEMTX_IS_REAL(Y))
   || (INPMTX_IS_COMPLEX_ENTRIES(A) && !DENSEMTX_IS_COMPLEX(Y)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n A's inputMode %d, Y's type %d \n", A, Y, alpha, X, 
          A->inputMode, Y->type) ;
   exit(-1) ;
}
if ( DenseMtx_rowIncrement(Y) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n Y's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(Y)) ;
   exit(-1) ;
}
incY = DenseMtx_columnIncrement(Y) ;
if ( DenseMtx_rowIncrement(X) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n X's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(X)) ;
   exit(-1) ;
}
incX = DenseMtx_columnIncrement(X) ;
x    = DenseMtx_entries(X) ;
y    = DenseMtx_entries(Y) ;
DenseMtx_dimensions(Y, &nrowY, &ncolY) ;
DenseMtx_dimensions(X, &nrowX, &ncolX) ;
if ( (nrhs = ncolY) != ncolX ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
          "\n Y's nrhs = %d, X's nrhs = %d",
          A, Y, alpha, X, nrhs, ncolX) ;
   exit(-1) ;
}
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A) ;
if ( ivec1 == NULL || ivec2 == NULL || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_T(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, Y, alpha, X, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
nent  = A->nent ;
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   double   rfac ;
   int      chev, col, ii, jrhs, off, row ;

   rfac = alpha[0] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               row = ivec2[ii] ; col = ivec1[ii] ;
               y[row] += dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               row = ivec2[ii] ; col = ivec1[ii] ;
               y[row] += rfac * dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               col = ivec2[ii] ; row = ivec1[ii] ;
               y[row] += dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               col = ivec2[ii] ; row = ivec1[ii] ;
               y[row] += rfac * dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               y[row] += dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               y[row] += rfac * dvec[ii]*x[col] ;
            }
            x += incX ; y += incY ;
         }
      }
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, jrhs, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row   = ivec2[ii] ; col   = ivec1[ii] ;
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec2[ii]  ; col = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec2[ii]  ; col = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec2[ii]  ; row = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec2[ii]  ; row = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec2[ii]  ; row = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^H*X

   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmm_H (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || Y == NULL || alpha == NULL || X == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
           "\n bad input\n", A, Y, alpha, X) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, Y, alpha, X, A->inputMode) ;
   exit(-1) ;
}
if ( ! DENSEMTX_IS_COMPLEX(Y) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
          "\n bad type %d for Y\n", A, Y, alpha, X, Y->type) ;
   exit(-1) ;
}
if ( ! DENSEMTX_IS_COMPLEX(X) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
          "\n bad type %d for X\n", A, Y, alpha, X, X->type) ;
   exit(-1) ;
}
if ( DenseMtx_rowIncrement(Y) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
          "\n Y's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(Y)) ;
   exit(-1) ;
}
incY = DenseMtx_columnIncrement(Y) ;
if ( DenseMtx_rowIncrement(X) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
          "\n X's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(X)) ;
   exit(-1) ;
}
incX = DenseMtx_columnIncrement(X) ;
x    = DenseMtx_entries(X) ;
y    = DenseMtx_entries(Y) ;
DenseMtx_dimensions(Y, &nrowY, &ncolY) ;
DenseMtx_dimensions(X, &nrowX, &ncolX) ;
if ( (nrhs = ncolY) != ncolX ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
          "\n Y's nrhs = %d, X's nrhs = %d",
          A, Y, alpha, X, nrhs, ncolX) ;
   exit(-1) ;
}
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A) ;
if ( ivec1 == NULL || ivec2 == NULL || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_mmm_H(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, Y, alpha, X, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
nent  = A->nent ;
if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, jrhs, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row   = ivec2[ii] ; col   = ivec1[ii] ;
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal + aimag*ximag ;
               y[2*row+1] += areal*ximag - aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec2[ii]  ; col = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal + aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag - aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec2[ii]  ; col = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal + aimag*ximag ;
               t2 = areal*ximag - aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec2[ii]  ; row = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal + aimag*ximag ;
               y[2*row+1] += areal*ximag - aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec2[ii]  ; row = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal + aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag - aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec2[ii]  ; row = ivec1[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal + aimag*ximag ;
               t2 = areal*ximag - aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal + aimag*ximag ;
               y[2*row+1] += areal*ximag - aimag*xreal ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal + aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag - aimag*xreal) ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  col = chev ; row = chev + off ;
               } else {
                  row = chev ; col = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               t1 = areal*xreal + aimag*ximag ;
               t2 = areal*ximag - aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
              where A is symmetric

   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_sym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || Y == NULL || alpha == NULL || X == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
           "\n bad input\n", A, Y, alpha, X) ;
   exit(-1) ;
}
if ( ! (INPMTX_IS_REAL_ENTRIES(A) || INPMTX_IS_COMPLEX_ENTRIES(A)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, Y, alpha, X, A->inputMode) ;
   exit(-1) ;
}
if ( ! (DENSEMTX_IS_REAL(Y) || DENSEMTX_IS_COMPLEX(Y)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n bad type %d for Y\n", A, Y, alpha, X, Y->type) ;
   exit(-1) ;
}
if ( ! (DENSEMTX_IS_REAL(X) || DENSEMTX_IS_COMPLEX(X)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n bad type %d for X\n", A, Y, alpha, X, X->type) ;
   exit(-1) ;
}
if ( DENSEMTX_IS_REAL(Y) != DENSEMTX_IS_REAL(X) ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n X's type %d, Y's type %d \n", A, Y, alpha, X, 
          X->type, Y->type) ;
   exit(-1) ;
}
if (  (INPMTX_IS_REAL_ENTRIES(A) && !DENSEMTX_IS_REAL(Y))
   || (INPMTX_IS_COMPLEX_ENTRIES(A) && !DENSEMTX_IS_COMPLEX(Y)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n A's inputMode %d, Y's type %d \n", A, Y, alpha, X, 
          A->inputMode, Y->type) ;
   exit(-1) ;
}
if ( DenseMtx_rowIncrement(Y) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n Y's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(Y)) ;
   exit(-1) ;
}
incY = DenseMtx_columnIncrement(Y) ;
if ( DenseMtx_rowIncrement(X) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n X's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(X)) ;
   exit(-1) ;
}
incX = DenseMtx_columnIncrement(X) ;
x    = DenseMtx_entries(X) ;
y    = DenseMtx_entries(Y) ;
DenseMtx_dimensions(Y, &nrowY, &ncolY) ;
DenseMtx_dimensions(X, &nrowX, &ncolX) ;
if ( (nrhs = ncolY) != ncolX ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
          "\n Y's nrhs = %d, X's nrhs = %d",
          A, Y, alpha, X, nrhs, ncolX) ;
   exit(-1) ;
}
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A) ;
if ( ivec1 == NULL || ivec2 == NULL || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmm(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, Y, alpha, X, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
nent  = A->nent ;
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   double   rfac ;
   int      chev, col, ii, jrhs, off, row ;

   rfac = alpha[0] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               row = ivec1[ii] ; col = ivec2[ii] ;
               y[row] += dvec[ii]*x[col] ;
               if ( row != col ) {
                  y[col] += dvec[ii]*x[row] ;
               }
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               row = ivec1[ii] ; col = ivec2[ii] ;
               y[row] += rfac * dvec[ii]*x[col] ;
               if ( row != col ) {
                  y[col] += rfac * dvec[ii]*x[row] ;
               }
            }
            x += incX ; y += incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               col = ivec1[ii] ; row = ivec2[ii] ;
               y[row] += dvec[ii]*x[col] ;
               if ( row != col ) {
                  y[col] += dvec[ii]*x[row] ;
               }
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               col = ivec1[ii] ; row = ivec2[ii] ;
               y[row] += rfac * dvec[ii]*x[col] ;
               if ( row != col ) {
                  y[col] += rfac * dvec[ii]*x[row] ;
               }
            }
            x += incX ; y += incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               y[row] += dvec[ii]*x[col] ;
               if ( row != col ) {
                  y[col] += dvec[ii]*x[row] ;
               }
            }
            x += incX ; y += incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = 0 ; ii < nent ; ii++ ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               y[row] += rfac * dvec[ii]*x[col] ;
               if ( row != col ) {
                  y[col] += rfac * dvec[ii]*x[row] ;
               }
            }
            x += incX ; y += incY ;
         }
      }
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, jrhs, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row   = ivec1[ii] ; col   = ivec2[ii] ;
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  y[2*col]   += areal*xreal - aimag*ximag ;
                  y[2*col+1] += areal*ximag + aimag*xreal ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec1[ii]  ; col = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  y[2*col]   += rfac*(areal*xreal - aimag*ximag) ;
                  y[2*col+1] += rfac*(areal*ximag + aimag*xreal) ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec1[ii]  ; col = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  t1 = areal*xreal - aimag*ximag ;
                  t2 = areal*ximag + aimag*xreal ;
                  y[2*col]   += rfac*t1 - ifac*t2 ;
                  y[2*col+1] += rfac*t2 + ifac*t1 ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  y[2*col]   += areal*xreal - aimag*ximag ;
                  y[2*col+1] += areal*ximag + aimag*xreal ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  y[2*col]   += rfac*(areal*xreal - aimag*ximag) ;
                  y[2*col+1] += rfac*(areal*ximag + aimag*xreal) ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  t1 = areal*xreal - aimag*ximag ;
                  t2 = areal*ximag + aimag*xreal ;
                  y[2*col]   += rfac*t1 - ifac*t2 ;
                  y[2*col+1] += rfac*t2 + ifac*t1 ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  y[2*col]   += areal*xreal - aimag*ximag ;
                  y[2*col+1] += areal*ximag + aimag*xreal ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  y[2*col]   += rfac*(areal*xreal - aimag*ximag) ;
                  y[2*col+1] += rfac*(areal*ximag + aimag*xreal) ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  t1 = areal*xreal - aimag*ximag ;
                  t2 = areal*ximag + aimag*xreal ;
                  y[2*col]   += rfac*t1 - ifac*t2 ;
                  y[2*col+1] += rfac*t2 + ifac*t1 ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
              where A is hermitian

   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_herm_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || Y == NULL || alpha == NULL || X == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
           "\n bad input\n", A, Y, alpha, X) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, Y, alpha, X, A->inputMode) ;
   exit(-1) ;
}
if ( ! DENSEMTX_IS_COMPLEX(Y) ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
          "\n bad type %d for Y\n", A, Y, alpha, X, Y->type) ;
   exit(-1) ;
}
if ( ! DENSEMTX_IS_COMPLEX(X) ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
          "\n bad type %d for X\n", A, Y, alpha, X, X->type) ;
   exit(-1) ;
}
if ( DenseMtx_rowIncrement(Y) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
          "\n Y's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(Y)) ;
   exit(-1) ;
}
incY = DenseMtx_columnIncrement(Y) ;
if ( DenseMtx_rowIncrement(X) != 1 ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
          "\n X's row increment = %d\n",
          A, Y, alpha, X, DenseMtx_rowIncrement(X)) ;
   exit(-1) ;
}
incX = DenseMtx_columnIncrement(X) ;
x    = DenseMtx_entries(X) ;
y    = DenseMtx_entries(Y) ;
DenseMtx_dimensions(Y, &nrowY, &ncolY) ;
DenseMtx_dimensions(X, &nrowX, &ncolX) ;
if ( (nrhs = ncolY) != ncolX ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
          "\n Y's nrhs = %d, X's nrhs = %d",
          A, Y, alpha, X, nrhs, ncolX) ;
   exit(-1) ;
}
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A) ;
if ( ivec1 == NULL || ivec2 == NULL || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_mmm(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, Y, alpha, X, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
nent  = A->nent ;
if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, jrhs, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row   = ivec1[ii] ; col   = ivec2[ii] ;
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  y[2*col]   += areal*xreal + aimag*ximag ;
                  y[2*col+1] += areal*ximag - aimag*xreal ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec1[ii]  ; col = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  y[2*col]   += rfac*(areal*xreal + aimag*ximag) ;
                  y[2*col+1] += rfac*(areal*ximag - aimag*xreal) ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               row = ivec1[ii]  ; col = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  t1 = areal*xreal + aimag*ximag ;
                  t2 = areal*ximag - aimag*xreal ;
                  y[2*col]   += rfac*t1 - ifac*t2 ;
                  y[2*col+1] += rfac*t2 + ifac*t1 ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  y[2*col]   += areal*xreal + aimag*ximag ;
                  y[2*col+1] += areal*ximag - aimag*xreal ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  y[2*col]   += rfac*(areal*xreal + aimag*ximag) ;
                  y[2*col+1] += rfac*(areal*ximag - aimag*xreal) ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               col = ivec1[ii]  ; row = ivec2[ii] ;
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  t1 = areal*xreal + aimag*ximag ;
                  t2 = areal*ximag - aimag*xreal ;
                  y[2*col]   += rfac*t1 - ifac*t2 ;
                  y[2*col+1] += rfac*t2 + ifac*t1 ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj] ; aimag = dvec[jj+1] ;
               xreal = x[2*col] ; ximag = x[2*col+1] ;
               y[2*row]   += areal*xreal - aimag*ximag ;
               y[2*row+1] += areal*ximag + aimag*xreal ;
               if ( row != col ) {
                  xreal = x[2*row] ; ximag = x[2*row+1] ;
                  y[2*col]   += areal*xreal + aimag*ximag ;
                  y[2*col+1] += areal*ximag - aimag*xreal ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else if ( ifac == 0.0 ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
               y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  y[2*col]   += rfac*(areal*xreal + aimag*ximag) ;
                  y[2*col+1] += rfac*(areal*ximag - aimag*xreal) ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      } else {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
               chev = ivec1[ii] ; off  = ivec2[ii] ;
               if ( off >= 0 ) {
                  row = chev ; col = chev + off ;
               } else {
                  col = chev ; row = chev - off ;
               }
               areal = dvec[jj]  ; aimag = dvec[jj+1] ;
               xreal = x[2*col]  ; ximag = x[2*col+1] ;
               t1 = areal*xreal - aimag*ximag ;
               t2 = areal*ximag + aimag*xreal ;
               y[2*row]   += rfac*t1 - ifac*t2 ;
               y[2*row+1] += rfac*t2 + ifac*t1 ;
               if ( row != col ) {
                  xreal = x[2*row]  ; ximag = x[2*row+1] ;
                  t1 = areal*xreal + aimag*ximag ;
                  t2 = areal*ximag - aimag*xreal ;
                  y[2*col]   += rfac*t1 - ifac*t2 ;
                  y[2*col+1] += rfac*t2 + ifac*t1 ;
               }
            }
            x += 2*incX ; y += 2*incY ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
