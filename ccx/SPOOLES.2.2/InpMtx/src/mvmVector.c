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
InpMtx_nonsym_mmmVector (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) {
int      nent ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || y == NULL || alpha == NULL || x == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector(%p,%p,%p,%p)"
           "\n bad input\n", A, y, alpha, x) ;
   exit(-1) ;
}
if ( ! (INPMTX_IS_REAL_ENTRIES(A) || INPMTX_IS_COMPLEX_ENTRIES(A)) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, y, alpha, x, A->inputMode) ;
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
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, y, alpha, x, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
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
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^T*X

   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmmVector_T (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) {
int      nent ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || y == NULL || alpha == NULL || x == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector_T(%p,%p,%p,%p)"
           "\n bad input\n", A, y, alpha, x) ;
   exit(-1) ;
}
if ( ! (INPMTX_IS_REAL_ENTRIES(A) || INPMTX_IS_COMPLEX_ENTRIES(A)) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector_T(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, y, alpha, x, A->inputMode) ;
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
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector_T(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, y, alpha, x, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
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
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^H*X

   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmmVector_H (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) {
int      nent ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || y == NULL || alpha == NULL || x == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector_H(%p,%p,%p,%p)"
           "\n bad input\n", A, y, alpha, x) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector_H(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, y, alpha, x, A->inputMode) ;
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
   fprintf(stderr, 
           "\n fatal error in InpMtx_nonsym_mmmVector_H(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, y, alpha, x, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
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
InpMtx_sym_mmmVector (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) {
int      nent ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || y == NULL || alpha == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmmVector(%p,%p,%p,%p)"
           "\n bad input\n", A, y, alpha, x) ;
   exit(-1) ;
}
if ( ! (INPMTX_IS_REAL_ENTRIES(A) || INPMTX_IS_COMPLEX_ENTRIES(A)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmmVector(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, y, alpha, x, A->inputMode) ;
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
   fprintf(stderr, "\n fatal error in InpMtx_sym_mmmVector(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, y, alpha, x, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
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
InpMtx_herm_mmmVector (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) {
int      nent ;
int      *ivec1, *ivec2 ;
double   *dvec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL || y == NULL || alpha == NULL || x == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_herm_mmmVector(%p,%p,%p,%p)"
           "\n bad input\n", A, y, alpha, x) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_herm_mmmVector(%p,%p,%p,%p)"
          "\n bad inputMode %d for A\n", A, y, alpha, x, A->inputMode) ;
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
   fprintf(stderr, 
           "\n fatal error in InpMtx_herm_mmmVector(%p,%p,%p,%p)"
           "\n ivec1 %p, ivec2 %p, dvec %p\n", 
           A, y, alpha, x, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
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
return ; }

/*--------------------------------------------------------------------*/
