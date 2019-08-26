/*  gmvm.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
static int checkInput ( InpMtx *A, double beta[], int ny, double y[],
   double alpha[], int nx, double x[], char *methodname ) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to compute y := beta*y + alpha*A*x
   where x and y are double vectors.

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL

   created -- 98nov14, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmvm (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) {
int      nent, rc ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, ny, y, alpha, nx, x, "InpMtx_nonsym_gmvm") ;
if ( rc != 1 ) {
   return(rc) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
/*
   ----------------
   scale y by beta
   ----------------
*/
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   if ( beta[0] == 0.0 ) {
      DVzero(ny, y) ;
   } else if ( beta[0] != 0.0 ) {
      DVscale(ny, y, beta[0]) ;
   }
} else {
   if ( beta[0] == 0.0 && beta[1] == 0.0 ) {
      ZVzero(ny, y) ;
   } else if ( beta[0] != 1.0 || beta[1] != 0.0 ) {
      ZVscale(ny, y, beta[0], beta[1]) ;
   }
}
/*
   --------------------------------
   data is stored as triples
   (deal with vector storage later)
   --------------------------------
*/
nent = A->nent ;
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   double   rfac ;
   int      chev, col, ii, off, row ;

   rfac = alpha[0] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 ) {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            row = ivec1[ii] ; col = ivec2[ii] ;
            y[row] += dvec[ii]*x[col] ;
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            row = ivec1[ii] ; col = ivec2[ii] ;
            y[row] += rfac * dvec[ii]*x[col] ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 ) {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            col = ivec1[ii] ; row = ivec2[ii] ;
            y[row] += dvec[ii]*x[col] ;
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            col = ivec1[ii] ; row = ivec2[ii] ;
            y[row] += rfac * dvec[ii]*x[col] ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 ) {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            chev = ivec1[ii] ; off = ivec2[ii] ;
            if ( off >= 0 ) {
               row = chev ; col = chev + off ;
            } else {
               col = chev ; row = chev - off ;
            }
            y[row] += dvec[ii]*x[col] ;
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            chev = ivec1[ii] ; off  = ivec2[ii] ;
            if ( off >= 0 ) {
               row = chev ; col = chev + off ;
            } else {
               col = chev ; row = chev - off ;
            }
            y[row] += rfac * dvec[ii]*x[col] ;
         }
      }
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row   = ivec1[ii] ; col   = ivec2[ii] ;
            areal = dvec[jj]  ; aimag = dvec[jj+1] ;
            xreal = x[2*col]  ; ximag = x[2*col+1] ;
            y[2*row]   += areal*xreal - aimag*ximag ;
            y[2*row+1] += areal*ximag + aimag*xreal ;
         }
      } else if ( ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row = ivec1[ii]  ; col = ivec2[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
            y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
         }
      } else {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row = ivec1[ii]  ; col = ivec2[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            t1 = areal*xreal - aimag*ximag ;
            t2 = areal*ximag + aimag*xreal ;
            y[2*row]   += rfac*t1 - ifac*t2 ;
            y[2*row+1] += rfac*t2 + ifac*t1 ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec1[ii]  ; row = ivec2[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += areal*xreal - aimag*ximag ;
            y[2*row+1] += areal*ximag + aimag*xreal ;
         }
      } else if ( ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec1[ii]  ; row = ivec2[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
            y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
         }
      } else {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec1[ii]  ; row = ivec2[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            t1 = areal*xreal - aimag*ximag ;
            t2 = areal*ximag + aimag*xreal ;
            y[2*row]   += rfac*t1 - ifac*t2 ;
            y[2*row+1] += rfac*t2 + ifac*t1 ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL

   created -- 98may02, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmvm_T (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) {
int      nent, rc ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, ny, y, alpha, nx, x, "InpMtx_nonsym_gmvm_T") ;
if ( rc != 1 ) {
   return(rc) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
/*
   ----------------
   scale y by beta
   ----------------
*/
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   if ( beta[0] == 0.0 ) {
      DVzero(ny, y) ;
   } else if ( beta[0] != 0.0 ) {
      DVscale(ny, y, beta[0]) ;
   }
} else {
   if ( beta[0] == 0.0 && beta[1] == 0.0 ) {
      ZVzero(ny, y) ;
   } else if ( beta[0] != 1.0 || beta[1] != 0.0 ) {
      ZVscale(ny, y, beta[0], beta[1]) ;
   }
}
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
         for ( ii = 0 ; ii < nent ; ii++ ) {
            row = ivec2[ii] ; col = ivec1[ii] ;
            y[row] += dvec[ii]*x[col] ;
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            row = ivec2[ii] ; col = ivec1[ii] ;
            y[row] += rfac * dvec[ii]*x[col] ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 ) {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            col = ivec2[ii] ; row = ivec1[ii] ;
            y[row] += dvec[ii]*x[col] ;
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            col = ivec2[ii] ; row = ivec1[ii] ;
            y[row] += rfac * dvec[ii]*x[col] ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 ) {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            chev = ivec1[ii] ; off = ivec2[ii] ;
            if ( off >= 0 ) {
               col = chev ; row = chev + off ;
            } else {
               row = chev ; col = chev - off ;
            }
            y[row] += dvec[ii]*x[col] ;
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            chev = ivec1[ii] ; off  = ivec2[ii] ;
            if ( off >= 0 ) {
               col = chev ; row = chev + off ;
            } else {
               row = chev ; col = chev - off ;
            }
            y[row] += rfac * dvec[ii]*x[col] ;
         }
      }
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, jrhs, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row   = ivec2[ii] ; col   = ivec1[ii] ;
            areal = dvec[jj]  ; aimag = dvec[jj+1] ;
            xreal = x[2*col]  ; ximag = x[2*col+1] ;
            y[2*row]   += areal*xreal - aimag*ximag ;
            y[2*row+1] += areal*ximag + aimag*xreal ;
         }
      } else if ( ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row = ivec2[ii]  ; col = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
            y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
         }
      } else {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row = ivec2[ii]  ; col = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            t1 = areal*xreal - aimag*ximag ;
            t2 = areal*ximag + aimag*xreal ;
            y[2*row]   += rfac*t1 - ifac*t2 ;
            y[2*row+1] += rfac*t2 + ifac*t1 ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec2[ii]  ; row = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += areal*xreal - aimag*ximag ;
            y[2*row+1] += areal*ximag + aimag*xreal ;
         }
      } else if ( ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec2[ii]  ; row = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += rfac*(areal*xreal - aimag*ximag) ;
            y[2*row+1] += rfac*(areal*ximag + aimag*xreal) ;
         }
      } else {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec2[ii]  ; row = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            t1 = areal*xreal - aimag*ximag ;
            t2 = areal*ximag + aimag*xreal ;
            y[2*row]   += rfac*t1 - ifac*t2 ;
            y[2*row+1] += rfac*t2 + ifac*t1 ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL
    -10 -- A is real

   created -- 98may02, cca
   ---------------------------------------------------
*/
int
InpMtx_nonsym_gmvm_H (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) {
int      nent, rc ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, ny, y, alpha, nx, x, "InpMtx_nonsym_gmvm_H") ;
if ( rc != 1 ) {
   return(rc) ;
}
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_nonsym_gmvm_H()"
           "\n A, X and Y are real\n") ;
   return(-10) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
/*
   ----------------
   scale y by beta
   ----------------
*/
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   if ( beta[0] == 0.0 ) {
      DVzero(ny, y) ;
   } else if ( beta[0] != 0.0 ) {
      DVscale(ny, y, beta[0]) ;
   }
} else {
   if ( beta[0] == 0.0 && beta[1] == 0.0 ) {
      ZVzero(ny, y) ;
   } else if ( beta[0] != 1.0 || beta[1] != 0.0 ) {
      ZVscale(ny, y, beta[0], beta[1]) ;
   }
}
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
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row   = ivec2[ii] ; col   = ivec1[ii] ;
            areal = dvec[jj]  ; aimag = dvec[jj+1] ;
            xreal = x[2*col]  ; ximag = x[2*col+1] ;
            y[2*row]   += areal*xreal + aimag*ximag ;
            y[2*row+1] += areal*ximag - aimag*xreal ;
         }
      } else if ( ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row = ivec2[ii]  ; col = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += rfac*(areal*xreal + aimag*ximag) ;
            y[2*row+1] += rfac*(areal*ximag - aimag*xreal) ;
         }
      } else {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            row = ivec2[ii]  ; col = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            t1 = areal*xreal + aimag*ximag ;
            t2 = areal*ximag - aimag*xreal ;
            y[2*row]   += rfac*t1 - ifac*t2 ;
            y[2*row+1] += rfac*t2 + ifac*t1 ;
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec2[ii]  ; row = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += areal*xreal + aimag*ximag ;
            y[2*row+1] += areal*ximag - aimag*xreal ;
         }
      } else if ( ifac == 0.0 ) {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec2[ii]  ; row = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            y[2*row]   += rfac*(areal*xreal + aimag*ximag) ;
            y[2*row+1] += rfac*(areal*ximag - aimag*xreal) ;
         }
      } else {
         for ( ii = jj = 0 ; ii < nent ; ii++, jj += 2 ) {
            col = ivec2[ii]  ; row = ivec1[ii] ;
            areal = dvec[jj] ; aimag = dvec[jj+1] ;
            xreal = x[2*col] ; ximag = x[2*col+1] ;
            t1 = areal*xreal + aimag*ximag ;
            t2 = areal*ximag - aimag*xreal ;
            y[2*row]   += rfac*t1 - ifac*t2 ;
            y[2*row+1] += rfac*t2 + ifac*t1 ;
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL

   created -- 98nov06, cca
   ----------------------------------------------------
*/
int
InpMtx_sym_gmvm (
   InpMtx     *A,
   double     beta[],
   int        ny, 
   double     y[],
   double     alpha[],
   int        nx, 
   double     x[]
) {
int      nent, rc ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, ny, y, alpha, nx, x, "InpMtx_sym_gmvm") ;
if ( rc != 1 ) {
   return(rc) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
/*
   ----------------
   scale y by beta
   ----------------
*/
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   if ( beta[0] == 0.0 ) {
      DVzero(ny, y) ;
   } else if ( beta[0] != 0.0 ) {
      DVscale(ny, y, beta[0]) ;
   }
} else {
   if ( beta[0] == 0.0 && beta[1] == 0.0 ) {
      ZVzero(ny, y) ;
   } else if ( beta[0] != 1.0 || beta[1] != 0.0 ) {
      ZVscale(ny, y, beta[0], beta[1]) ;
   }
}
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
         for ( ii = 0 ; ii < nent ; ii++ ) {
            row = ivec1[ii] ; col = ivec2[ii] ;
            y[row] += dvec[ii]*x[col] ;
            if ( row != col ) {
               y[col] += dvec[ii]*x[row] ;
            }
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            row = ivec1[ii] ; col = ivec2[ii] ;
            y[row] += rfac * dvec[ii]*x[col] ;
            if ( row != col ) {
               y[col] += rfac * dvec[ii]*x[row] ;
            }
         }
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 ) {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            col = ivec1[ii] ; row = ivec2[ii] ;
            y[row] += dvec[ii]*x[col] ;
            if ( row != col ) {
               y[col] += dvec[ii]*x[row] ;
            }
         }
      } else {
         for ( ii = 0 ; ii < nent ; ii++ ) {
            col = ivec1[ii] ; row = ivec2[ii] ;
            y[row] += rfac * dvec[ii]*x[col] ;
            if ( row != col ) {
               y[col] += rfac * dvec[ii]*x[row] ;
            }
         }
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 ) {
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
      } else {
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
      }
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   double   aimag, areal, ifac, rfac, t1, t2, ximag, xreal ;
   int      chev, col, ii, jj, jrhs, off, row ;

   rfac = alpha[0] ; ifac = alpha[1] ;
   if ( INPMTX_IS_BY_ROWS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL
    -10 -- A, X and Y are real

   created -- 98nov06, cca
   --------------------------------------------------
*/
int
InpMtx_herm_gmvm (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) {
int      nent, rc ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
rc = checkInput(A, beta, ny, y, alpha, nx, x, "InpMtx_herm_gmvm") ;
if ( rc != 1 ) {
   return(rc) ;
}
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   fprintf(stderr, "\n fatal error in InpMtx_herm_gmvm()"
           "\n A, X and Y are real\n") ;
   return(-10) ;
}
ivec1 = InpMtx_ivec1(A) ;
ivec2 = InpMtx_ivec2(A) ;
dvec  = InpMtx_dvec(A)  ;
/*
   ----------------
   scale y by beta
   ----------------
*/
if ( INPMTX_IS_REAL_ENTRIES(A) ) {
   if ( beta[0] == 0.0 ) {
      DVzero(ny, y) ;
   } else if ( beta[0] != 0.0 ) {
      DVscale(ny, y, beta[0]) ;
   }
} else {
   if ( beta[0] == 0.0 && beta[1] == 0.0 ) {
      ZVzero(ny, y) ;
   } else if ( beta[0] != 1.0 || beta[1] != 0.0 ) {
      ZVscale(ny, y, beta[0], beta[1]) ;
   }
}
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
      }
   } else if ( INPMTX_IS_BY_COLUMNS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
      }
   } else if ( INPMTX_IS_BY_CHEVRONS(A) ) {
      if ( rfac == 1.0 && ifac == 0.0 ) {
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
      } else if ( ifac == 0.0 ) {
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
      } else {
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
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL

   created -- 98may02, cca
   --------------------------------------------------
*/
static int
checkInput (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[],
   char       *methodname
) {
double   *dvec ;
int      typeX, typeY ;
int      *ivec1, *ivec2 ;

if ( A == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n A is NULL\n", methodname) ;
   return(-1) ;
}
if ( ! INPMTX_IS_REAL_ENTRIES(A) && ! INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n type of A is %d, invalid\n", methodname, A->inputMode) ;
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
if ( ny <= 0 ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n ny = %d\n", methodname, ny) ;
   return(-5) ;
}
if ( y == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n y is NULL\n", methodname) ;
   return(-6) ;
}
if ( alpha == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n alpha is NULL\n", methodname) ;
   return(-7) ;
}
if ( nx <= 0 ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n nx = %d\n", methodname, nx) ;
   return(-8) ;
}
if ( x == NULL ) {
   fprintf(stderr, "\n fatal error in %s()"
           "\n x is NULL\n", methodname) ;
   return(-9) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
