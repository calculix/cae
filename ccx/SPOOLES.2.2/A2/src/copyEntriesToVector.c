/*  copyEntriesToVector.c  */

#include "../A2.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- copy entries to a vector. the portion copied
              can be a union of the strict lower portion,
              the diagonal portion, and the strict upper
              portion. there is one restriction, if the strict
              lower and strict upper are to be copied, the
              diagonal will also be copied.

   length -- length of dvec[]
   dvec[] -- vector to receive matrix entries
   copyflag  -- flag to denote what part of the entries to move
      A2_STRICT_LOWER --> move strict lower entries
      A2_LOWER        --> move lower entries (includes the diagonal)
      A2_DIAGONAL     --> move diagonal entries
      A2_UPPER        --> move upper entries (includes the diagonal)
      A2_STRICT_UPPER --> move strict upper entries
      A2_ALL_ENTRIES  --> move all entries
   storeflag -- flag to denote how to store entries in dvec[]
      A2_BY_ROWS    --> store by rows
      A2_BY_COLUMNS --> store by columns

   return value -- number of entries copied

   created  -- 97jun03, cca, dkw
   modified -- 98may25, cca
   ----------------------------------------------------------------
*/
int
A2_copyEntriesToVector (
   A2       *mtx,
   int      length,
   double   *dvec,
   int      copyflag, 
   int      storeflag
) {
int      inc1, inc2, kk, ncol, ndiag, nent, nrow ;
/*
   --------------------------------------------
   check the input, get dimensions and pointers
   and check that length is large enough
   --------------------------------------------
*/
if (  mtx == NULL || length < 0 || dvec == NULL ) {
   fprintf(stderr,
           "\n fatal error in A2_copyEntriesToVector(%p,%d,%p,,%d,%d)"
           "\n bad input\n", mtx, length, dvec, copyflag, storeflag) ;
   exit(-1) ;
}
if ( copyflag  < 1 || copyflag > 6 ) {
   fprintf(stderr,
           "\n fatal error in A2_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n bad input\n" 
           "\n copyflag = %d, must be\n" 
           "\n    1 --> strictly lower entries"
           "\n    2 --> lower entries"
           "\n    3 --> diagonal entries"
           "\n    4 --> strictly upper entries"
           "\n    5 --> upper entries"
           "\n    6 --> all entries",
           mtx, length, dvec, copyflag, storeflag, copyflag) ;
   exit(-1) ;
}
if ( storeflag  < 0 || storeflag > 1 ) {
   fprintf(stderr,
           "\n fatal error in A2_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n bad input\n" 
           "\n storeflag = %d, must be\n" 
           "\n    0 --> store by rows"
           "\n    1 --> store by columns",
           mtx, length, dvec, copyflag, storeflag, storeflag) ;
   exit(-1) ;
}
nrow = mtx->n1 ;
ncol = mtx->n2 ;
inc1 = mtx->inc1 ;
inc2 = mtx->inc2 ;
if ( nrow >= ncol ) {
   ndiag = ncol ;
} else {
   ndiag = nrow ;
}
/*
   ------------------------------------------
   compute the number of entries to be copied
   ------------------------------------------
*/
switch ( copyflag ) {
case A2_STRICT_LOWER : /* strictly lower entries  */
   nent = (ndiag*(ndiag - 1))/2 + (nrow - ndiag)*ncol ;
   break ;
case A2_LOWER : /* lower entries  */
   nent = (ndiag*(ndiag + 1))/2 + (nrow - ndiag)*ncol ;
   break ;
case A2_DIAGONAL : /* diagonal entries  */
   nent = ndiag ;
   break ;
case A2_UPPER : /* upper entries  */
   nent = (ndiag*(ndiag + 1))/2 + (ncol - ndiag)*nrow ;
   break ;
case A2_STRICT_UPPER : /* strictly upper entries  */
   nent = (ndiag*(ndiag - 1))/2 + (ncol - ndiag)*nrow ;
   break ;
case A2_ALL_ENTRIES : /* all entries  */
   nent = nrow*ncol ;
   break ;
default :
   break ;
}
if ( nent > length ) {
   fprintf(stderr,
           "\n fatal error in A2_copyEntriesToVector(%p,%d,%p,%d,%d)"
           "\n nent = %d, buffer length = %d", 
           mtx, length, dvec, copyflag, storeflag, nent, length) ;
   exit(-1) ;
}
/*
   --------------------------------------------
   make life simple, unit stride through dvec[]
   --------------------------------------------
*/
kk = 0 ;
if ( storeflag == A2_BY_ROWS ) {
   int      irow, jstart, jcol, jend, jj ;
   double   *row ;
/*
   --------------
   loop over rows
   --------------
*/
   switch ( copyflag ) {
   case A2_STRICT_LOWER :
/*
      ----------------------
      strictly lower entries
      ----------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += inc1 ) {
            jstart = 0 ;
            jend   = (irow < ndiag) ? irow - 1 : ndiag - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[kk] = row[jj] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += 2*inc1 ) {
            jstart = 0 ;
            jend   = (irow < ndiag) ? irow - 1 : ndiag - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[2*kk]   = row[2*jj]   ;
               dvec[2*kk+1] = row[2*jj+1] ;
            }
         }
      }
      break ;
   case A2_LOWER :
/*
      ------------------------------------
      lower entries including the diagonal
      ------------------------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += inc1 ) {
            jstart = 0    ;
            jend   = (irow < ndiag) ? irow : ndiag - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[kk] = row[jj] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += 2*inc1 ) {
            jstart = 0    ;
            jend   = (irow < ndiag) ? irow : ndiag - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[2*kk]   = row[2*jj]   ;
               dvec[2*kk+1] = row[2*jj+1] ;
            }
         }
      }
      break ;
   case A2_DIAGONAL :
/*
      -----------------
      just the diagonal
      -----------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < ndiag ; 
               irow++, row += inc1, kk++ ) {
            dvec[kk] = row[irow*inc2] ;
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < ndiag ; 
               irow++, row += 2*inc1, kk++ ) {
            dvec[2*kk] = row[2*irow*inc2] ;
            dvec[2*kk+1] = row[2*irow*inc2+1] ;
         }
      }
      break ;
   case A2_UPPER :
/*
      ------------------------------------
      upper entries including the diagonal
      ------------------------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += inc1 ) {
            jstart = irow ;
            jend   = ncol - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[kk] = row[jj] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += 2*inc1 ) {
            jstart = irow ;
            jend   = ncol - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[2*kk]   = row[2*jj]   ;
               dvec[2*kk+1] = row[2*jj+1] ;
            }
         }
      }
      break ;
   case A2_STRICT_UPPER :
/*
      --------------------------
      strictly the upper entries
      --------------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += inc1 ) {
            jstart = irow + 1 ;
            jend   = ncol - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[kk] = row[jj] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += 2*inc1 ) {
            jstart = irow + 1 ;
            jend   = ncol - 1 ;
            for ( jcol = jstart, jj = jcol*inc2 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[2*kk]   = row[2*jj]   ;
               dvec[2*kk+1] = row[2*jj+1] ;
            }
         }
      }
      break ;
   case A2_ALL_ENTRIES :
/*
      -----------
      all entries
      -----------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += inc1 ) {
            jstart = 0 ;
            jend   = ncol - 1 ;
            for ( jcol = jstart, jj = 0 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[kk] = row[jj] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( irow = 0, row = mtx->entries, kk = 0 ; 
               irow < nrow ; 
               irow++, row += 2*inc1 ) {
            jstart = 0 ;
            jend   = ncol - 1 ;
            for ( jcol = jstart, jj = 0 ;
                  jcol <= jend ; 
                  jcol++, jj += inc2, kk++ ) {
               dvec[2*kk]   = row[2*jj]   ;
               dvec[2*kk+1] = row[2*jj+1] ;
            }
         }
      }
      break ;
   default :
      break ;
   }
} else {
   int      iend, ii, irow, istart, jcol ;
   double   *col ;
/*
   -----------------
   loop over columns
   -----------------
*/
   kk = 0 ;
   switch ( copyflag ) {
   case A2_STRICT_LOWER :
/*
      ----------------------
      strictly lower entries
      ----------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += inc2 ) {
            istart = jcol + 1 ;
            iend   = nrow - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[kk] = col[ii] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += 2*inc2 ) {
            istart = jcol + 1 ;
            iend   = nrow - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[2*kk]   = col[2*ii]   ;
               dvec[2*kk+1] = col[2*ii+1] ;
            }
         }
      }
      break ;
   case A2_LOWER :
/*
      ------------------------------------
      lower entries including the diagonal
      ------------------------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += inc2 ) {
            istart = jcol ;
            iend   = nrow - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[kk] = col[ii] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += 2*inc2 ) {
            istart = jcol ;
            iend   = nrow - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[2*kk]   = col[2*ii]   ;
               dvec[2*kk+1] = col[2*ii+1] ;
            }
         }
      }
      break ;
   case A2_DIAGONAL :
/*
      -----------------
      just the diagonal
      -----------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ndiag ; 
               jcol++, col += inc2, kk++ ) {
            dvec[kk] = col[jcol*inc1] ;
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ndiag ; 
               jcol++, col += 2*inc2, kk++ ) {
            dvec[2*kk]   = col[2*jcol*inc1]   ;
            dvec[2*kk+1] = col[2*jcol*inc1+1] ;
         }
      }
      break ;
   case A2_UPPER :
/*
      ------------------------------------
      upper entries including the diagonal
      ------------------------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += inc2 ) {
            istart = 0 ;
            iend   = (jcol < ndiag) ? jcol : ndiag - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[kk] = col[ii] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += 2*inc2 ) {
            istart = 0 ;
            iend   = (jcol < ndiag) ? jcol : ndiag - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[2*kk]   = col[2*ii]   ;
               dvec[2*kk+1] = col[2*ii+1] ;
            }
         }
      }
      break ;
   case A2_STRICT_UPPER :
/*
      --------------------------
      strictly the upper entries
      --------------------------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += inc2 ) {
            istart = 0 ;
            iend   = (jcol < ndiag) ? jcol - 1 : ndiag - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[kk] = col[ii] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += 2*inc2 ) {
            istart = 0 ;
            iend   = (jcol < ndiag) ? jcol - 1 : ndiag - 1 ;
            for ( irow = istart, ii = irow*inc1 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[2*kk]   = col[2*ii]   ;
               dvec[2*kk+1] = col[2*ii+1] ;
            }
         }
      }
      break ;
   case A2_ALL_ENTRIES :
/*
      -----------
      all entries
      -----------
*/
      if ( A2_IS_REAL(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += inc2 ) {
            istart = 0 ;
            iend   = nrow - 1 ;
            for ( irow = istart, ii = 0 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[kk] = col[ii] ;
            }
         }
      } else if ( A2_IS_COMPLEX(mtx) ) {
         for ( jcol = 0, col = mtx->entries, kk = 0 ; 
               jcol < ncol ; 
               jcol++, col += 2*inc2 ) {
            istart = 0 ;
            iend   = nrow - 1 ;
            for ( irow = istart, ii = 0 ;
                  irow <= iend ; 
                  irow++, ii += inc1, kk++ ) {
               dvec[2*kk]   = col[2*ii]   ;
               dvec[2*kk+1] = col[2*ii+1] ;
            }
         }
      }
      break ;
   default :
      break ;
   }
}
return(kk) ; }

/*--------------------------------------------------------------------*/
