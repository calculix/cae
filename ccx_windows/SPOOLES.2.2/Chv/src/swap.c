/*  swap.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   swap rows irow and jrow

   created -- 98apr30, cca
   -----------------------
*/
void
Chv_swapRows (
   Chv   *chv,
   int   irow,
   int   jrow
) {
double   dtmp ;
double   *entries ;
int      ii, ioff, itmp, joff, nD, nL, nrow, nU, stride ;
int      *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || irow < 0 || jrow < 0 ) {
   fprintf(stderr, "\n fatal error in Chv_swapRows(%p,%d,%d)"
           "\n bad input\n", chv, irow, jrow) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout,"\n %% Chv_swapRows(%p,%d,%d)", chv, irow, jrow) ;
fprintf(stdout,"\n %% chv->symflag = %d", chv->symflag) ;
fflush(stdout) ;
#endif
if ( irow == jrow ) {
   return ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
if ( irow >= nD || jrow >= nD ) {
   fprintf(stderr, "\n fatal error in Chv_swapRows(%p,%d,%d)"
           "\n rows must be less than nD = %d", chv, irow, irow, nD) ;
   exit(-1) ;
}
entries = Chv_entries(chv) ;
#if MYDEBUG > 0
fprintf(stdout,"\n %% nD = %d, nL = %d, nU = %d, entries = %p",
        nD, nL, nU, entries) ;
fflush(stdout) ;
#endif
if ( entries == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_swapRows(%p,%d,%d)"
           "\n bad input, entries = %p, nD = %d\n", 
           chv, irow, jrow, entries, nD) ;
   exit(-1) ;
}
if ( ! (CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_swapRows(%p,%d,%d)"
           "\n type = %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           chv, irow, jrow, chv->type) ;
   exit(-1) ;
}
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
   ------------------------------------------------
   call method for hermitian and symmetric chevrons
   ------------------------------------------------
*/
   Chv_swapRowsAndColumns(chv, irow, jrow) ;
} else if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
   ------------------------
   swap the two row indices
   ------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n %% getting ready to swap row indices") ;
   fflush(stdout) ;
#endif
   Chv_rowIndices(chv, &nrow, &rowind) ;
#if MYDEBUG > 0
   fprintf(stdout,
         "\n %% before: rowind = %p, rowind[%d] = %d, rowind[%d] = %d",
         rowind, irow, rowind[irow], jrow, rowind[jrow]) ;
   IVfprintf(stdout, nrow, rowind) ;
   fflush(stdout) ;
#endif
   itmp         = rowind[irow] ;
   rowind[irow] = rowind[jrow] ;
   rowind[jrow] = itmp         ;
#if MYDEBUG > 0
   fprintf(stdout,
           "\n %% after: rowind = %p, rowind[%d] = %d, rowind[%d] = %d",
           rowind, irow, rowind[irow], jrow, rowind[jrow]) ;
   IVfprintf(stdout, nrow, rowind) ;
   fflush(stdout) ;
#endif
/*
   --------------------------------------------------
   set irow = min(irow, jrow), jrow = max(irow, jrow)
   --------------------------------------------------
*/
   if ( irow > jrow ) {
      itmp = irow ;
      irow = jrow ;
      jrow = itmp ;
   }
/*
   --------------------------------
   swap the entries in the two rows
   --------------------------------
*/
   ioff   = nD + nL - 1 - irow ;
   joff   = nD + nL - 1 - jrow ;
   stride = 2*nD + nL + nU - 1 ;
   if ( CHV_IS_REAL(chv) ) {
      for ( ii = 0 ; ii < irow ; ii++ ) {
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff += stride, joff += stride, stride -= 2 ;
      }
      for ( ii = irow ; ii < jrow ; ii++ ) {
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff += 1, joff += stride, stride -= 2 ;
      }
      for ( ii = jrow ; ii < nD + nU ; ii++ ) {
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff += 1, joff += 1 ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( ii = 0 ; ii < irow ; ii++ ) {
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff += stride, joff += stride, stride -= 2 ;
      }
      for ( ii = irow ; ii < jrow ; ii++ ) {
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff += 1, joff += stride, stride -= 2 ;
      }
      for ( ii = jrow ; ii < nD + nU ; ii++ ) {
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff += 1, joff += 1 ;
      }
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_swapRows(%p,%d,%d)"
           "\n bad symmetryflag %d\n", chv, irow, jrow, chv->symflag) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   swap columns icol and jcol

   created -- 98apr30, cca
   --------------------------
*/
void
Chv_swapColumns (
   Chv   *chv,
   int   icol,
   int   jcol
) {
double   dtmp ;
double   *entries ;
int      ii, ioff, itmp, joff, ncol, nD, nL, nU, stride ;
int      *colind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || icol < 0 || jcol < 0 ) {
   fprintf(stderr, "\n fatal error in Chv_swapColumns(%p,%d,%d)"
           "\n bad input\n", chv, icol, jcol) ;
   exit(-1) ;
}
if ( icol == jcol ) {
   return ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
if ( entries == NULL || icol >= nD || jcol >= nD ) {
   fprintf(stderr, "\n fatal error in Chv_swapColumns(%p,%d,%d)"
           "\n bad input, entries = %p, nD = %d\n", 
           chv, icol, jcol, entries, nD) ;
   exit(-1) ;
}
if ( ! (CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_swapColumns(%p,%d,%d)"
           "\n type = %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           chv, icol, jcol, chv->type) ;
   exit(-1) ;
}
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
   ------------------------------------------------
   call method for symmetric and hermitian chevrons
   ------------------------------------------------
*/
   Chv_swapRowsAndColumns(chv, icol, jcol) ;
} else if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
   ---------------------------
   swap the two column indices
   ---------------------------
*/
   Chv_columnIndices(chv, &ncol, &colind) ;
   itmp         = colind[icol] ;
   colind[icol] = colind[jcol] ;
   colind[jcol] = itmp         ;
/*
   --------------------------------------------------
   set icol = min(icol, jcol), jcol = max(icol, jcol)
   --------------------------------------------------
*/
   if ( icol > jcol ) {
      itmp = icol ;
      icol = jcol ;
      jcol = itmp ;
   }
/*
   -----------------------------------
   swap the entries in the two columns
   -----------------------------------
*/
   ioff   = nD + nL - 1 + icol ;
   joff   = nD + nL - 1 + jcol ;
   stride = 2*nD + nL + nU - 3 ;
   if ( CHV_IS_REAL(chv) ) {
      for ( ii = 0 ; ii < icol ; ii++ ) {
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff += stride, joff += stride, stride -= 2 ;
      }
      for ( ii = icol ; ii < jcol ; ii++ ) {
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff -= 1, joff += stride, stride -= 2 ;
      }
      for ( ii = jcol ; ii < nD + nL ; ii++ ) {
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff -= 1, joff -= 1 ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( ii = 0 ; ii < icol ; ii++ ) {
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff += stride, joff += stride, stride -= 2 ;
      }
      for ( ii = icol ; ii < jcol ; ii++ ) {
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff -= 1, joff += stride, stride -= 2 ;
      }
      for ( ii = jcol ; ii < nD + nL ; ii++ ) {
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff -= 1, joff -= 1 ;
      }
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_swapColumns(%p,%d,%d)"
           "\n bad symmetryflag %d\n", chv, icol, jcol, chv->symflag) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   swap rows and columns ii and jj

   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_swapRowsAndColumns (
   Chv   *chv,
   int   ii,
   int   jj
) {
double   dtmp ;
double   *entries ;
int      iiloc, ioff, itmp, jjloc, joff, kk, ncol, nD, nL, nU, stride ;
int      *colind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || ii < 0 || jj < 0 ) {
   fprintf(stderr, 
           "\n fatal error in Chv_swapRowsAndColumns(%p,%d,%d)"
           "\n bad input\n", chv, ii, jj) ;
   exit(-1) ;
}
if ( ii == jj ) {
   return ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
if ( entries == NULL || ii >= nD || jj >= nD ) {
   fprintf(stderr, 
           "\n fatal error in Chv_swapRowsAndColumns(%p,%d,%d)"
           "\n bad input, entries = %p, nD = %d\n", 
           chv, ii, jj, entries, nD) ;
   exit(-1) ;
}
if ( ! (CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_swapRowsAndColumns(%p,%d,%d)"
           "\n type = %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           chv, ii, jj, chv->type) ;
   exit(-1) ;
}
if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
   --------------------------------------
   call methods for nonsymmetric chevrons
   --------------------------------------
*/
   Chv_swapRows(chv, ii, jj) ;
   Chv_swapColumns(chv, ii, jj) ;
} else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
   ---------------------------
   swap the two column indices
   ---------------------------
*/
   Chv_columnIndices(chv, &ncol, &colind) ;
   itmp       = colind[ii] ;
   colind[ii] = colind[jj] ;
   colind[jj] = itmp       ;
/*
   --------------------------------------
   set ii = min(ii, jj), jj = max(ii, jj)
   --------------------------------------
*/
   if ( ii > jj ) {
      itmp = ii   ;
      ii   = jj   ;
      jj   = itmp ;
   }
/*
   ------------------------------------------
   step 1, swap A(0:ii-1,ii) and A(0:ii-1,jj)
   ------------------------------------------
*/
   ioff   = ii ;
   joff   = jj ;
   stride = nD + nU - 1 ;
   if ( CHV_IS_REAL(chv) ) {
      for ( kk = 0 ; kk < ii ; kk++ ) {
#if MYDEBUG > 0
      fprintf(stdout, 
      "\n\n %% first step, kk = %d, ioff = %d, joff = %d, stride = %d", 
         kk, ioff, joff, stride) ;
#endif
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff += stride, joff += stride, stride -= 1 ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( kk = 0 ; kk < ii ; kk++ ) {
#if MYDEBUG > 0
      fprintf(stdout, 
      "\n\n %% first step, kk = %d, ioff = %d, joff = %d, stride = %d", 
         kk, ioff, joff, stride) ;
#endif
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff += stride, joff += stride, stride -= 1 ;
      }
   }
   iiloc = ioff ;
/*
   -------------------------------------------
   2. swap A(ii,ii+1:jj-1) and A(ii+1:jj-1,jj)
      tricky part if hermitian
   -------------------------------------------
*/
   ioff++, joff += stride, stride -= 1 ;
   if ( CHV_IS_REAL(chv) ) {
      for ( kk = ii + 1 ; kk < jj ; kk++ ) {
         double   aiikk, akkjj ;
#if MYDEBUG > 0
      fprintf(stdout, 
     "\n\n second step, %% kk = %d, ioff = %d, joff = %d, stride = %d", 
        kk, ioff, joff, stride) ;
#endif
         aiikk = entries[ioff] ;
         akkjj = entries[joff] ;
         entries[ioff] = akkjj ;
         entries[joff] = aiikk ;
         ioff += 1, joff += stride, stride -= 1 ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( kk = ii + 1 ; kk < jj ; kk++ ) {
         double   imagiikk, imagkkjj, realiikk, realkkjj ;
#if MYDEBUG > 0
      fprintf(stdout, 
     "\n\n second step, %% kk = %d, ioff = %d, joff = %d, stride = %d", 
        kk, ioff, joff, stride) ;
#endif
         realiikk = entries[2*ioff] ;
         imagiikk = entries[2*ioff+1] ;
         realkkjj = entries[2*joff] ;
         imagkkjj = entries[2*joff+1] ;
         entries[2*ioff] = realkkjj ;
         entries[2*joff] = realiikk ;
         if ( CHV_IS_SYMMETRIC(chv) ) {
            entries[2*ioff+1] = imagkkjj ;
            entries[2*joff+1] = imagiikk ;
         } else {
            entries[2*ioff+1] = -imagkkjj ;
            entries[2*joff+1] = -imagiikk ;
         }
         ioff += 1, joff += stride, stride -= 1 ;
      }
      if ( CHV_IS_HERMITIAN(chv) ) {
/*
         ----------------------------------
         set (ii,jj) entry to its conjugate
         ----------------------------------
*/
         entries[2*ioff+1] = -entries[2*ioff+1] ;
      }
   }
   jjloc = joff ;
/*
   -----------------------------
   4. swap A(ii,ii) and A(jj,jj)
   -----------------------------
*/
#if MYDEBUG > 0
      fprintf(stdout, 
     "\n\n %% third step, iiloc = %d, jjloc = %d", iiloc, jjloc) ;
#endif
   if ( CHV_IS_REAL(chv) ) {
      dtmp             = entries[iiloc]   ;
      entries[iiloc]   = entries[jjloc]   ;
      entries[jjloc]   = dtmp             ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      dtmp               = entries[2*iiloc]   ;
      entries[2*iiloc]   = entries[2*jjloc]   ;
      entries[2*jjloc]   = dtmp               ;
      dtmp               = entries[2*iiloc+1] ;
      entries[2*iiloc+1] = entries[2*jjloc+1] ;
      entries[2*jjloc+1] = dtmp               ;
   }
/*
   -------------------------------------------------
   swap
   4. swap A(ii,jj+1:nD+nU-1) and A(jj,jj+1:nD+nU-1)
   -------------------------------------------------
*/
   ioff++, joff++ ;
   if ( CHV_IS_REAL(chv) ) {
      for ( kk = jj + 1 ; kk < nD + nU ; kk++ ) {
#if MYDEBUG > 0
         fprintf(stdout, 
      "\n\n %% fourth step, kk = %d, ioff = %d, joff = %d, stride = %d",
           kk, ioff, joff, stride) ;
#endif
         dtmp            = entries[ioff]   ;
         entries[ioff]   = entries[joff]   ;
         entries[joff]   = dtmp            ;
         ioff += 1, joff += 1 ;
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( kk = jj + 1 ; kk < nD + nU ; kk++ ) {
#if MYDEBUG > 0
         fprintf(stdout, 
      "\n\n %% fourth step, kk = %d, ioff = %d, joff = %d, stride = %d",
           kk, ioff, joff, stride) ;
#endif
         dtmp              = entries[2*ioff]   ;
         entries[2*ioff]   = entries[2*joff]   ;
         entries[2*joff]   = dtmp              ;
         dtmp              = entries[2*ioff+1] ;
         entries[2*ioff+1] = entries[2*joff+1] ;
         entries[2*joff+1] = dtmp              ;
         ioff += 1, joff += 1 ;
      }
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_swapRowsAndColumns(%p,%d,%d)"
           "\n bad symmetryflag %d\n", chv, ii, jj, chv->symflag) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
