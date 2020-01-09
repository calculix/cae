/*  gmmm.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
static int checkInput ( InpMtx *A, double beta[], DenseMtx *Y,
   double alpha[], DenseMtx *X, char *methodname ) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
   where X and Y are DenseMtx objects, and X and Y
   must be column major.

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- alpha is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal

   created -- 98may02, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmmm (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY, rc ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, Y, alpha, X, "InpMtx_nonsym_gmmm") ;
if ( rc != 1 ) {
   return(rc) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
incY  = Y->inc2    ;
y     = Y->entries ;
nrhs  = Y->ncol    ;
incX  = X->inc2    ;
x     = X->entries ;
/*
   ----------------
   scale Y by beta
   ----------------
*/
DenseMtx_scale(Y, beta) ;
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
nent = A->nent ;
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
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A^T*X

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- alpha is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal

   created -- 98may02, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmmm_T (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY, rc ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, Y, alpha, X, "InpMtx_nonsym_gmmm_T") ;
if ( rc != 1 ) {
   return(rc) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
incY  = Y->inc2    ;
y     = Y->entries ;
nrhs  = Y->ncol    ;
incX  = X->inc2    ;
x     = X->entries ;
/*
   ----------------
   scale Y by beta
   ----------------
*/
DenseMtx_scale(Y, beta) ;
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
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
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A^H*X

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- alpha is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal
    -16 -- A, X and Y are real

   created -- 98may02, cca
   ---------------------------------------------------
*/
int
InpMtx_nonsym_gmmm_H (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY, rc ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, Y, alpha, X, "InpMtx_nonsym_gmmm_H") ;
if ( rc != 1 ) {
   return(rc) ;
}
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_gmmm_H()"
           "\n A, X and Y are real\n") ;
   return(-16) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
incY  = Y->inc2    ;
y     = Y->entries ;
nrhs  = Y->ncol    ;
incX  = X->inc2    ;
x     = X->entries ;
/*
   ----------------
   scale Y by beta
   ----------------
*/
DenseMtx_scale(Y, beta) ;
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
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
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
              where A is symmetric

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- alpha is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal

   created -- 98nov06, cca
   ----------------------------------------------------
*/
int
InpMtx_sym_gmmm (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY, rc ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, Y, alpha, X, "InpMtx_sym_gmmm") ;
if ( rc != 1 ) {
   return(rc) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
incY  = Y->inc2    ;
y     = Y->entries ;
nrhs  = Y->ncol    ;
incX  = X->inc2    ;
x     = X->entries ;
/*
   ----------------
   scale Y by beta
   ----------------
*/
DenseMtx_scale(Y, beta) ;
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
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
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
              where A is hermitian

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- alpha is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal
    -16 -- A, X and Y are real

   created -- 98nov06, cca
   --------------------------------------------------
*/
int
InpMtx_herm_gmmm (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) {
int      incX, incY, ncolX, ncolY, nent, nrhs, nrowX, nrowY, rc ;
int      *ivec1, *ivec2 ;
double   *dvec, *x, *y ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, Y, alpha, X, "InpMtx_herm_gmmm") ;
if ( rc != 1 ) {
   return(rc) ;
}
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_gmmm()"
           "\n A, X and Y are real\n") ;
   return(-16) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
incY  = Y->inc2    ;
y     = Y->entries ;
nrhs  = Y->ncol    ;
incX  = X->inc2    ;
x     = X->entries ;
/*
   ----------------
   scale Y by beta
   ----------------
*/
DenseMtx_scale(Y, beta) ;
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
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
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to check the input

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- alpha is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal

   created -- 98may02, cca
   --------------------------------------------------
*/
static int
checkInput (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   char       *methodname
) {
double   *dvec, *x, *y ;
int      colincX, colincY, rowincX, rowincY, ncolX, ncolY, 
         nrowX, nrowY, typeA, typeX, typeY ;
int      *ivec1, *ivec2 ;

if ( A == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n A is NULL\n", methodname) ;
   return(-1) ;
}
typeA = A->inputMode ;
if ( typeA != SPOOLES_REAL && typeA != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n type of A is %d, invalid\n", methodname, typeA) ;
   return(-2) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A) ;
if ( ivec1 == NULL || ivec2 == NULL || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           methodname, ivec1, ivec2, dvec) ;
   return(-3) ;
}
if ( beta == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n beta is NULL\n", methodname) ;
   return(-4) ;
}
if ( Y == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n Y is NULL\n", methodname) ;
   return(-5) ;
}
typeY = Y->type ;
if ( typeY != SPOOLES_REAL && typeY != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n type of Y is %d, invalid\n", methodname, typeY) ;
   return(-6) ;
}
DenseMtx_dimensions(Y, &nrowY, &ncolY) ;
colincY = DenseMtx_columnIncrement(Y) ;
rowincY = DenseMtx_rowIncrement(Y) ;
if ( nrowY <= 0 || ncolY <= 0 || rowincY != 1 || colincY != nrowY ) {
   fprintf(stderr, "\n fatal error in %s()"
          "\n nrowY %d, ncolY %d, rowincY %d, colincY %d\n",
          methodname, nrowY, ncolY, rowincY, colincY) ;
   return(-7) ;
}
y = DenseMtx_entries(Y) ;
if ( y == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n Y's entries are NULL\n", methodname) ;
   return(-8) ;
}
if ( alpha == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n alpha is NULL\n", methodname) ;
   return(-9) ;
}
if ( X == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n X is NULL\n", methodname) ;
   return(-10) ;
}
typeX = X->type ;
if ( typeX != SPOOLES_REAL && typeX != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n type of X is %d, invalid\n", methodname, typeX) ;
   return(-11) ;
}
DenseMtx_dimensions(X, &nrowX, &ncolX) ;
colincX = DenseMtx_columnIncrement(X) ;
rowincX = DenseMtx_rowIncrement(X) ;
if ( nrowX <= 0 || ncolX <= 0 || rowincX != 1 || colincX != nrowX ) {
   fprintf(stderr, "\n fatal error in %s()"
          "\n nrowX %d, ncolX %d, rowincX %d, colincX %d\n",
          methodname, nrowX, ncolX, rowincX, colincX) ;
   return(-12) ;
}
x = DenseMtx_entries(X) ;
if ( x == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n entries of X are NULL\n", methodname) ;
   return(-13) ;
}
if ( typeA != typeX || typeA != typeY || typeX != typeY ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n types do not match, typeA %d, typeX %d, typeY %d\n",
           methodname, typeA, typeX, typeY) ;
   return(-14) ;
}
if ( ncolY != ncolX ) {
   fprintf(stderr, "\n fatal error in %s()"
          "\n ncolY %d, ncolX %d\n", methodname, ncolY, ncolX) ;
   return(-15) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
