/*  search.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   find the first unmarked entry in 
   the diagonal with largest magnitude
   if ( mark[jj] == tag ) then
      we can compare this entry
   endif

   created -- 98apr30, cca
   -----------------------------------
*/
int
Chv_maxabsInDiagonal11 (
   Chv      *chv,
   int      mark[],
   int      tag,
   double   *pmaxval
) {
double   maxval, val ;
double   *entries ;
int      jcol, jj, nD, nL, nU, off, stride ;
/*
   --------------
   check the data
   --------------
*/
if ( chv == NULL || mark == NULL || pmaxval == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInDiagonal11(%p,%p,%d,%p)"
           "\n bad input\n", chv, mark, tag, pmaxval) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 2 ;
      for ( jj = 0 ; jj < nD ; jj++ ) {
         if ( mark[jj] == tag ) {
            val = fabs(entries[off]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride -= 2 ;
      }
   } else if ( CHV_IS_SYMMETRIC(chv) ) {
/*
      -----------------
      symmetric chevron
      -----------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = 0 ;
      stride = nD + nU ;
      for ( jj = 0 ; jj < nD ; jj++ ) {
         if ( mark[jj] == tag ) {
            val = fabs(entries[off]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride-- ;
      }
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_maxabsInDiagonal11(%p,%p,%d,%p)"
              "\n type = SPOOLES_REAL, bad symflag %d \n", 
              chv, mark, tag, pmaxval, chv->symflag) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = nD + nL - 1 ;
      stride = 2*nD + nL + nU - 2 ;
      for ( jj = 0 ; jj < nD ; jj++ ) {
         if ( mark[jj] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride -= 2 ;
      }
   } else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
      ------------------------------
      hermitian or symmetric chevron
      ------------------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = 0 ;
      stride = nD + nU ;
      for ( jj = 0 ; jj < nD ; jj++ ) {
         if ( mark[jj] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride-- ;
      }
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_maxabsInDiagonal11(%p,%p,%d,%p)"
              "\n type = SPOOLES_COMPLEX, bad symflag %d \n", 
              chv, mark, tag, pmaxval, chv->symflag) ;
      exit(-1) ;
   }
} else {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInDiagonal11(%p,%p,%d,%p)"
           "\n bad type, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, mark, tag, pmaxval) ;
   exit(-1) ;
}
*pmaxval = maxval ;
   
return(jcol) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   find the first unmarked entry in 
   row irow with largest magnitude
   if ( colmark[jj] == tag ) then
      we can examined this entry
   endif
   only entries in the (1,1) block are examined

   created -- 98apr30, cca
   --------------------------------------------
*/
int
Chv_maxabsInRow11 (
   Chv      *chv,
   int      irow,
   int      colmark[],
   int      tag,
   double   *pmaxval
) {
double   maxval, val ;
double   *entries ;
int      jcol, jj, nD, nL, nU, off, stride ;
/*
   --------------
   check the data
   --------------
*/
if ( chv == NULL || irow < 0 || colmark == NULL || pmaxval == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInRow11(%p,%d,%p,%d,%p)"
           "\n bad input\n", chv, irow, colmark, tag, pmaxval) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = nD + nL - 1 - irow ;
      stride = 2*nD + nL + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         if ( colmark[jj] == tag ) {
            val = fabs(entries[off]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( jj = irow ; jj < nD ; jj++, off++ ) {
         if ( colmark[jj] == tag ) {
            val = fabs(entries[off]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) ) {
/*
      -----------------
      symmetric chevron
      -----------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = irow ;
      stride = nD + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         if ( colmark[jj] == tag ) {
            val = fabs(entries[off]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride-- ;
      }
      for ( jj = irow ; jj < nD ; jj++, off++ ) {
         if ( colmark[jj] == tag ) {
            val = fabs(entries[off]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
      }
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_maxabsInRow11(%p,%d,%p,%d,%p)"
              "\n type is SPOOLES_REAL, bad symflag %d \n", 
              chv, irow, colmark, tag, pmaxval, chv->symflag) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = nD + nL - 1 - irow ;
      stride = 2*nD + nL + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         if ( colmark[jj] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( jj = irow ; jj < nD ; jj++, off++ ) {
         if ( colmark[jj] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
      ------------------------------
      hermitian or symmetric chevron
      ------------------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = irow ;
      stride = nD + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         if ( colmark[jj] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
         off += stride ;
         stride-- ;
      }
      for ( jj = irow ; jj < nD ; jj++, off++ ) {
         if ( colmark[jj] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( jcol == -1 || maxval < val ) {
               jcol = jj ; maxval = val ;
            }
         }
      }
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_maxabsInRow11(%p,%d,%p,%d,%p)"
              "\n type is SPOOLES_COMPLEX, bad symflag %d \n", 
              chv, irow, colmark, tag, pmaxval, chv->symflag) ;
      exit(-1) ;
   }
} else {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInRow11(%p,%d,%p,%d,%p)"
           "\n bad type, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, irow, colmark, tag, pmaxval) ;
   exit(-1) ;
}
*pmaxval = maxval ;
   
return(jcol) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   find the first unmarked entry in 
   column jcol with largest magnitude
   if ( rowmark[ii] == tag ) then
      we can examined this entry
   endif
   only entries in the (1,1) block are examined

   created -- 98apr30, cca
   --------------------------------------------
*/
int
Chv_maxabsInColumn11 (
   Chv      *chv,
   int      jcol,
   int      rowmark[],
   int      tag,
   double   *pmaxval
) {
double   maxval, val ;
double   *entries ;
int      irow, ii, nD, nL, nU, off, stride ;
/*
   --------------
   check the data
   --------------
*/
if ( chv == NULL || jcol < 0 || rowmark == NULL || pmaxval == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInColumn11(%p,%d,%p,%d,%p)"
           "\n bad input\n", chv, jcol, rowmark, tag, pmaxval) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
irow    = -1  ;
maxval  = 0.0 ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      maxval = 0.0 ;
      off    = nD + nL + jcol - 1 ;
      stride = 2*nD + nL + nU - 3 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         if ( rowmark[ii] == tag ) {
            val = fabs(entries[off]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( ii = jcol ; ii < nD ; ii++, off-- ) {
         if ( rowmark[ii] == tag ) {
            val = fabs(entries[off]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) ) {
/*
      -----------------
      symmetric chevron
      -----------------
*/
      maxval = 0.0 ;
      off    = jcol ;
      stride = nD + nU - 1 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         if ( rowmark[ii] == tag ) {
            val = fabs(entries[off]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
         off += stride ;
         stride-- ;
      }
      for ( ii = jcol ; ii < nD ; ii++, off++ ) {
         if ( rowmark[ii] == tag ) {
            val = fabs(entries[off]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
      }
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      maxval = 0.0 ;
      off    = nD + nL + jcol - 1 ;
      stride = 2*nD + nL + nU - 3 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         if ( rowmark[ii] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( ii = jcol ; ii < nD ; ii++, off-- ) {
         if ( rowmark[ii] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
      ------------------------------
      hermitian or symmetric chevron
      ------------------------------
*/
      maxval = 0.0 ;
      off    = jcol ;
      stride = nD + nU - 1 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         if ( rowmark[ii] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
         off += stride ;
         stride-- ;
      }
      for ( ii = jcol ; ii < nD ; ii++, off++ ) {
         if ( rowmark[ii] == tag ) {
            val = Zabs(entries[2*off], entries[2*off+1]) ;
            if ( irow == -1 || maxval < val ) {
               irow = ii ; maxval = val ;
            }
         }
      }
   }
} else {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInColumn11(%p,%d,%p,%d,%p)"
           "\n bad type, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, jcol, rowmark, tag, pmaxval) ;
   exit(-1) ;
}
*pmaxval = maxval ;

return(irow) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   return the location of the first entry 
   with largest magnitude in row irow. 
   *pmaxval is filled with its magnitude.

   created -- 98apr30, cca
   --------------------------------------
*/
int
Chv_maxabsInRow (
   Chv      *chv,
   int      irow,
   double   *pmaxval
) {
double   maxval, val ;
double   *entries ;
int      jcol, jj, ncol, nD, nL, nU, off, stride ;
/*
   --------------
   check the data
   --------------
*/
if ( chv == NULL || irow < 0 || pmaxval == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_maxabsInRow(%p,%d,%p)"
           "\n bad input\n", chv, irow, pmaxval) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
ncol    = nD + nU ;
jcol    = -1  ;
maxval  = 0.0 ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = nD + nL - 1 - irow ;
      stride = 2*nD + nL + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         val = fabs(entries[off]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( jj = irow ; jj < ncol ; jj++, off++ ) {
         val = fabs(entries[off]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) ) {
/*
      -----------------
      symmetric chevron
      -----------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = irow ;
      stride = nD + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         val = fabs(entries[off]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
         off += stride ;
         stride-- ;
      }
      for ( jj = irow ; jj < ncol ; jj++, off++ ) {
         val = fabs(entries[off]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
      }
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = nD + nL - 1 - irow ;
      stride = 2*nD + nL + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( jj = irow ; jj < ncol ; jj++, off++ ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
      ------------------------------
      hermitian or symmetric chevron
      ------------------------------
*/
      jcol   = -1  ;
      maxval = 0.0 ;
      off    = irow ;
      stride = nD + nU - 1 ;
      for ( jj = 0 ; jj < irow ; jj++ ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
         off += stride ;
         stride-- ;
      }
      for ( jj = irow ; jj < ncol ; jj++, off++ ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( jcol == -1 || maxval < val ) {
            jcol = jj ; maxval = val ;
         }
      }
   }
} else {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInRow(%p,%d,%p)"
          "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX \n", 
           chv, irow, pmaxval, chv->symflag) ;
   exit(-1) ;
}
*pmaxval = maxval ;
   
return(jcol) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   return the location of the first entry 
   with largest magnitude in column jcol. 
   *pmaxval is filled with its magnitude.

   created -- 98apr30, cca
   --------------------------------------
*/
int
Chv_maxabsInColumn (
   Chv      *chv,
   int      jcol,
   double   *pmaxval
) {
double   maxval, val ;
double   *entries ;
int      irow, ii, nD, nL, nrow, nU, off, stride ;
/*
   --------------
   check the data
   --------------
*/
if ( chv == NULL || jcol < 0 || pmaxval == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_maxabsInColumn(%p,%d,%p)"
           "\n bad input\n", chv, jcol, pmaxval) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
nrow    = nD + nL ;
irow    = -1  ;
maxval  = 0.0 ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      maxval = 0.0 ;
      off    = nD + nL + jcol - 1 ;
      stride = 2*nD + nL + nU - 3 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         val = fabs(entries[off]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( ii = jcol ; ii < nrow ; ii++, off-- ) {
         val = fabs(entries[off]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) ) {
/*
      -----------------
      symmetric chevron
      -----------------
*/
      maxval = 0.0 ;
      off    = jcol ;
      stride = nD + nU - 1 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         val = fabs(entries[off]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
         off += stride ;
         stride-- ;
      }
      for ( ii = jcol ; ii < nrow ; ii++, off++ ) {
         val = fabs(entries[off]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
      }
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      --------------------
      nonsymmetric chevron
      --------------------
*/
      maxval = 0.0 ;
      off    = nD + nL + jcol - 1 ;
      stride = 2*nD + nL + nU - 3 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
         off += stride ;
         stride -= 2 ;
      }
      for ( ii = jcol ; ii < nrow ; ii++, off-- ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
      }
   } else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
      ------------------------------
      hermitian or symmetric chevron
      ------------------------------
*/
      maxval = 0.0 ;
      off    = jcol ;
      stride = nD + nU - 1 ;
      for ( ii = 0 ; ii < jcol ; ii++ ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
         off += stride ;
         stride-- ;
      }
      for ( ii = jcol ; ii < nrow ; ii++, off++ ) {
         val = Zabs(entries[2*off], entries[2*off+1]) ;
         if ( irow == -1 || maxval < val ) {
            irow = ii ; maxval = val ;
         }
      }
   }
} else {
   fprintf(stderr, 
           "\n fatal error in Chv_maxabsInColumn(%p,%d,%p)"
           "\n bad symflag %d \n", chv, jcol, pmaxval, chv->symflag) ;
   exit(-1) ;
}
*pmaxval = maxval ;

return(irow) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   return the magnitude of a quasimax entry from the unmarked 
   rows and columns and fill *pirow and *pjcol with its location

   created -- 98apr30, cca
   -------------------------------------------------------------
*/
double
Chv_quasimax (
   Chv   *chv,
   int   rowmark[],
   int   colmark[],
   int   tag,
   int   *pirow,
   int   *pjcol
) {
double   maxval ;
int      irow, jcol, nD, qcol, qrow ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || rowmark == NULL || colmark == NULL
  || pirow == NULL || pjcol == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_quasimax(%p,%p,%p,%d,%p,%p)"
           "\n bad input\n", chv, rowmark, colmark, tag, pirow, pjcol) ;
   exit(-1) ;
}
if ( ! CHV_IS_NONSYMMETRIC(chv) ) {
   fprintf(stderr, 
           "\n fatal error in Chv_quasimax(%p,%p,%p,%d,%p,%p)"
           "\n chv->symflag =  %d"
           "\n chevron is not symmetric or hermitian"
           "\n method cannot be used \n", 
           chv, rowmark, colmark, tag, pirow, pjcol, chv->symflag) ;
   exit(-1) ;
}
nD = chv->nD ;
/*
   ----------------------
   set the default values
   ----------------------
*/
*pirow = *pjcol = -1 ;
maxval = 0.0 ;
/*
   -----------------------
   find an unmarked column
   -----------------------
*/
for ( jcol = 0 ; jcol < nD ; jcol++ ) {
   if ( colmark[jcol] == tag ) {
      break ;
   }
}
if ( jcol == nD ) {
/*
   ---------------------------
   no unmarked columns, return
   ---------------------------
*/
   return(maxval) ;
}
/*
   ----------------------------------------
   find a maxabs entry in the unmarked rows
   ----------------------------------------
*/
irow = Chv_maxabsInColumn11(chv, jcol, rowmark, tag, &maxval) ;
#if MYDEBUG > 0
fprintf(stdout, "\n first maxabs = %12.4e in (%d,%d)",
        Chv_entry(chv, irow, jcol), irow, jcol) ;
fflush(stdout) ;
#endif
if ( irow == -1 ) {
/*
   ------------------------
   no unmarked rows, return
   ------------------------
*/
   return(maxval) ; 
}
while ( 1 ) {
/*
   ----------------------------------
   find a new maxabs entry in the row
   ----------------------------------
*/
   qcol = Chv_maxabsInRow11(chv, irow, colmark, tag, &maxval) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n new maxabs   = %12.4e in (%d,%d)",
           Chv_entry(chv, irow, qcol), irow, qcol) ;
   fflush(stdout) ;
#endif
   if ( qcol == jcol ) {
/*
      --------------------------------------------
      same as before, break out of the search loop
      --------------------------------------------
*/
      break ;
   }
   jcol = qcol ;
/*
   -------------------------------------
   find a new maxabs entry in the column
   -------------------------------------
*/
   qrow = Chv_maxabsInColumn11(chv, jcol, rowmark, tag, &maxval) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n new maxabs   = %12.4e in (%d,%d)",
           Chv_entry(chv, qrow, jcol), qrow, jcol) ;
   fflush(stdout) ;
#endif
   if ( qrow == irow ) {
/*
      --------------------------------------------
      same as before, break out of the search loop
      --------------------------------------------
*/
      break ;
   }
   irow = qrow ;
}
/*
   --------------------------------------------------------
   set the row and column where the quasimax entry is found
   --------------------------------------------------------
*/
*pjcol = jcol ;
*pirow = irow ;

return(maxval) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   find a 1x1 or 2x2 pivot using the fast Bunch-Parlett algorithm.
   used only with symmetric chevrons.

   created -- 98apr30, cca
   ---------------------------------------------------------------
*/
void
Chv_fastBunchParlettPivot (
   Chv   *chv,
   int   mark[],
   int   tag,
   int   *pirow,
   int   *pjcol
) {
double   maxdiag, gamma_r, gamma_s ;
double   *entries ;
int      nD, nL, nU, r, s, t ;
/*
   --------------
   check the data
   --------------
*/
if ( chv == NULL || mark == NULL || pirow == NULL || pjcol == NULL ) {
   fprintf(stderr, 
          "\n fatal error in Chv_fastBunchParlettPivot(%p,%p,%d,%p,%p)"
          "\n bad input\n",
          chv, mark, tag, pirow, pjcol) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
/*
   ----------------------
   set the default values
   ----------------------
*/
*pirow = *pjcol = -1 ;
/*
   ------------------------------------------
   find an unmarked entry of maximum magitude
   ------------------------------------------
*/
r = Chv_maxabsInDiagonal11(chv, mark, tag, &maxdiag) ;
if ( r == -1 ) {
/*
   -------------------------------------------------------
   all rows and columns are marked, return without success
   -------------------------------------------------------
*/
   *pirow = *pjcol = -1 ;
   return ;
}
/*
   -------------------------------------------------------------------
   find the offdiagonal entry of maximum magnitude in row and column r
   -------------------------------------------------------------------
*/
s = -1 ;
gamma_r = 0.0 ;
s = Chv_maxabsInRow11(chv, r, mark, tag, &gamma_r) ;
if ( s == -1 ) {
/*
   ---------------------------------------------
   r is the only unmarked row and column, return
   ---------------------------------------------
*/
   *pirow = *pjcol = r ;
   return ; 
}
if ( maxdiag >= 0.6404 * gamma_r ) {
/*
   -------------------------
   1 x 1 pivot is acceptable
   -------------------------
*/
   *pirow = *pjcol = r ;
   return ;
} else {
/*
   ---------------
   loop until done
   ---------------
*/
   while ( 1 ) {
/*
   ----------------------------------------------
   find
      t = index of max off diag entry in column s
      gamma_s is its magnitude
   ----------------------------------------------
*/
      t = Chv_maxabsInRow11(chv, s, mark, tag, &gamma_s) ;
      if ( t == r || gamma_s == gamma_r ) {
/*
         -------------------------
         2 x 2 pivot is acceptable
         -------------------------
*/
         *pirow = r ;
         *pjcol = s ;
         return ;
      } else {
/*
         --------------------------------
         keep looking for a local maximum
         --------------------------------
*/
         r = s ;
         gamma_r = gamma_s ;
         s = t ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
