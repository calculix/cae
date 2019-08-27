/*  util.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   shift the indices, entries and adjust the nD dimension.
   note: shift can be positive or negative

   created -- 98apr30, cca
   -------------------------------------------------------
*/
void
Chv_shift (
   Chv   *chv,
   int   shift
) {
int   ii, stride ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_shift(%p,%d)"
           "\n bad input\n", chv, shift) ;
   exit(-1) ;
}
if ( shift == 0 ) {
   return ;
}
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) ) {
   /*
      --------------------
      chevron is symmetric
      --------------------
   */
      chv->colind += shift ;
      if ( shift < 0 ) {
         stride = chv->nD + chv->nU + 1 ;
         for ( ii = shift ; ii < 0 ; ii++ ) {
            chv->entries -= stride ;
            stride++ ;
         }
      } else {
         stride = chv->nD + chv->nU ;
         for ( ii = 0 ; ii < shift ; ii++ ) {
            chv->entries += stride ;
            stride-- ;
         }
      }
      chv->nD -= shift ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      -----------------------
      chevron is nonsymmetric
      -----------------------
*/
      chv->rowind += shift ;
      chv->colind += shift ;
      if ( shift < 0 ) {
         stride = 2*chv->nD + chv->nL + chv->nU + 1 ;
         for ( ii = shift ; ii < 0 ; ii++ ) {
            chv->entries -= stride ;
            stride += 2 ;
         }
      } else {
         stride = 2*chv->nD + chv->nL + chv->nU - 1 ;
         for ( ii = 0 ; ii < shift ; ii++ ) {
            chv->entries += stride ;
            stride -= 2 ;
         }
      }
      chv->nD -= shift ;
   } else {
      fprintf(stderr, "\n fatal error in Chv_shift(%p,%d)"
              "\n type is SPOOLES_REAL, symflag = %d" 
              "\n must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n", 
              chv, shift, chv->symflag) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   /*
      --------------------
      chevron is symmetric
      --------------------
   */
      chv->colind += shift ;
      if ( shift < 0 ) {
         stride = chv->nD + chv->nU + 1 ;
         for ( ii = shift ; ii < 0 ; ii++ ) {
            chv->entries -= 2*stride ;
            stride++ ;
         }
      } else {
         stride = chv->nD + chv->nU ;
         for ( ii = 0 ; ii < shift ; ii++ ) {
            chv->entries += 2*stride ;
            stride-- ;
         }
      }
      chv->nD -= shift ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
      -----------------------
      chevron is nonsymmetric
      -----------------------
*/
      chv->rowind += shift ;
      chv->colind += shift ;
      if ( shift < 0 ) {
         stride = 2*chv->nD + chv->nL + chv->nU + 1 ;
         for ( ii = shift ; ii < 0 ; ii++ ) {
            chv->entries -= 2*stride ;
            stride += 2 ;
         }
      } else {
         stride = 2*chv->nD + chv->nL + chv->nU - 1 ;
         for ( ii = 0 ; ii < shift ; ii++ ) {
            chv->entries += 2*stride ;
            stride -= 2 ;
         }
      }
      chv->nD -= shift ;
   } else {
      fprintf(stderr, "\n fatal error in Chv_shift(%p,%d)"
        "\n type is SPOOLES_COMPLEX, symflag = %d" 
        "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
        "\n or SPOOLES_NONSYMMETRIC\n",
        chv, shift, chv->symflag) ;
      exit(-1) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   return the maximum magnitude of the entries in the chevron

   created -- 98apr30, cca
   ----------------------------------------------------------
*/
double
Chv_maxabs (
   Chv   *chv 
) {
double   maxabs ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_maxabs(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
maxabs = 0.0 ;
if ( CHV_IS_REAL(chv) ) {
   int   loc ;
   maxabs = DVmaxabs(Chv_nent(chv), Chv_entries(chv), &loc) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   maxabs = ZVmaxabs(Chv_nent(chv), Chv_entries(chv)) ;
} else {
   fprintf(stderr, "\n fatal error in Chv_maxabs(%p)"
           "\n type is %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, chv->type) ;
   exit(-1) ;
}
return(maxabs) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   return the frobenius norm of the entries in the chevron

   created -- 98apr30, cca
   -------------------------------------------------------
*/
double
Chv_frobNorm (
   Chv   *chv 
) {
double   sum ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_frobNorm(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
if ( CHV_IS_REAL(chv) ) {
   double   value, *entries ;
   int      ii, nent ;
  
   nent = Chv_nent(chv) ;
   entries = Chv_entries(chv) ;
   for ( ii = 0, sum = 0.0 ; ii < nent ; ii++ ) {
      value = entries[ii] ;
      sum += value*value ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   double   imag, real, *entries ;
   int      ii, nent ;
  
   nent = Chv_nent(chv) ;
   entries = Chv_entries(chv) ;
   for ( ii = 0, sum = 0.0 ; ii < nent ; ii++ ) {
      real = entries[2*ii]   ;
      imag = entries[2*ii+1] ;
      sum += real*real + imag*imag ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_frobNorm(%p)"
           "\n type is %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, chv->type) ;
   exit(-1) ;
}
return(sqrt(sum)) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   subtract chvI from chvJ

   created -- 98apr30, cca
   -----------------------
*/
void
Chv_sub (
   Chv   *chvJ,
   Chv   *chvI 
) {
double   *entriesI, *entriesJ ;
int      ii, nDI, nDJ, nent, nLI, nLJ, nUI, nUJ ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvI == NULL || chvJ == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_sub(%p,%p)"
           "\n bad input\n", chvI, chvJ) ;
   exit(-1) ;
}
Chv_dimensions(chvJ, &nDJ, &nLJ, &nUJ) ;
Chv_dimensions(chvI, &nDI, &nLI, &nUI) ;
if ( nDJ != nDI || nLJ != nLI || nUJ != nUI ) {
   fprintf(stderr, "\n fatal error in Chv_sub(%p,%p)"
           "\n dimensions do not match\n", chvJ, chvI) ;
   exit(-1) ;
}
entriesJ = Chv_entries(chvJ) ;
entriesI = Chv_entries(chvI) ;
if ( entriesJ == NULL || entriesI == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_sub(%p,%p)"
           "\n entriesJ = %p, entriesI = %p\n", 
           chvJ, chvI, entriesJ, entriesI) ;
   exit(-1) ;
}
if ( CHV_IS_REAL(chvJ) && CHV_IS_REAL(chvI) ) {
   nent = Chv_nent(chvJ) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      entriesJ[ii] -= entriesI[ii] ;
   }
} else if ( CHV_IS_COMPLEX(chvJ) && CHV_IS_COMPLEX(chvI) ) {
   nent = Chv_nent(chvJ) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      entriesJ[2*ii]   -= entriesI[2*ii]   ;
      entriesJ[2*ii+1] -= entriesI[2*ii+1] ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_sub(%p,%p)"
           "\n typeJ = %d, typeI = %d"
           "\n both must be SPOOLES_REAL or SPOOLES_COMPLEX",
           chvJ, chvI, chvJ->type, chvI->type) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   zero the entries in the chevron

   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_zero (
   Chv   *chv
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_zero(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
if ( CHV_IS_REAL(chv) ) {
   DVzero(Chv_nent(chv), Chv_entries(chv)) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   ZVzero(Chv_nent(chv), Chv_entries(chv)) ;
} else {
   fprintf(stderr, "\n fatal error in Chv_zero(%p)"
           "\n type = %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, chv->type) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   fill A2 object with (1,1) block

   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_fill11block (
   Chv   *chv,
   A2    *mtx
) {
double   *entries ;
int      nD, nL, nU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || mtx == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_fill11block(%p,%p)"
           "\n bad input\n", chv, mtx) ;
   exit(-1) ;
}
if ( ! (CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_fill11block(%p,%p)"
           "\n type = %d, must be real or complex\n", 
           chv, mtx, chv->type) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
if ( CHV_IS_REAL(chv) ) {
   A2_init(mtx, SPOOLES_REAL, nD, nD, 1, nD, NULL) ;
   A2_zero(mtx) ;
   if ( CHV_IS_SYMMETRIC(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
/*
      ----------------
      chv is symmetric
      ----------------
*/
      ioff = 0 ;
      istride = nD + nU ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
         A2_setRealEntry(mtx, ii, ii, entries[ioff]) ;
         for ( jj = ii + 1, iioff = ioff + 1 ; 
               jj < nD ; 
               jj++, iioff++ ) {
            A2_setRealEntry(mtx, ii, jj, entries[iioff]) ;
            A2_setRealEntry(mtx, jj, ii, 0.0) ;
         }
         ioff += istride ;
         istride-- ;
      }
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
/*
      ---------------------
      chv is nonsymmetric
      ---------------------
*/
      ioff = nD + nL - 1 ;
      istride = 2*nD + nU + nL - 2 ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
         A2_setRealEntry(mtx, ii, ii, entries[ioff]) ;
         for ( jj = ii + 1, iioff = ioff + 1 ; 
               jj < nD ; 
               jj++, iioff++ ) {
            A2_setRealEntry(mtx, ii, jj, entries[iioff]) ;
         }
         for ( jj = ii + 1, iioff = ioff - 1 ; 
               jj < nD ; 
               jj++, iioff-- ) {
            A2_setRealEntry(mtx, jj, ii, entries[iioff]) ;
         }
         ioff += istride ;
         istride -= 2 ;
      }
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   A2_init(mtx, SPOOLES_COMPLEX, nD, nD, 1, nD, NULL) ;
   A2_zero(mtx) ;
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
/*
      ----------------
      chv is symmetric
      ----------------
*/
      ioff = 0 ;
      istride = nD + nU ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
         A2_setComplexEntry(mtx, ii, ii, 
                            entries[2*ioff], entries[2*ioff+1]) ;
         for ( jj = ii + 1, iioff = ioff + 1 ; 
               jj < nD ; 
               jj++, iioff++ ) {
            A2_setComplexEntry(mtx, ii, jj, 
                               entries[2*iioff], entries[2*iioff+1]) ;
            A2_setComplexEntry(mtx, jj, ii, 0.0, 0.0) ;
         }
         ioff += istride ;
         istride-- ;
      }
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
   /*
      ---------------------
      chv is nonsymmetric
      ---------------------
   */
      ioff = nD + nL - 1 ;
      istride = 2*nD + nU + nL - 2 ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
   #if MYDEBUG > 0
         fprintf(stdout, 
              "\n ii %d, ioff %d, setting (%d,%d) = (%20.12e,%20.12e)",
              ii, ioff, ii, ii, entries[2*ioff], entries[2*ioff+1]) ;
   #endif
         A2_setComplexEntry(mtx, ii, ii, 
                            entries[2*ioff], entries[2*ioff+1]) ;
         for ( jj = ii + 1, iioff = ioff + 1 ; 
               jj < nD ; 
               jj++, iioff++ ) {
   #if MYDEBUG > 0
            fprintf(stdout, 
              "\n jj %d, iioff %d, setting (%d,%d) = (%20.12e,%20.12e)",
              jj, iioff, ii, ii, entries[2*iioff], entries[2*iioff+1]) ;
   #endif
            A2_setComplexEntry(mtx, ii, jj, 
                               entries[2*iioff], entries[2*iioff+1]) ;
         }
         for ( jj = ii + 1, iioff = ioff - 1 ; 
               jj < nD ; 
               jj++, iioff-- ) {
   #if MYDEBUG > 0
            fprintf(stdout, 
              "\n jj %d, iioff %d, setting (%d,%d) = (%20.12e,%20.12e)",
              jj, iioff, ii, ii, entries[2*iioff], entries[2*iioff+1]) ;
   #endif
            A2_setComplexEntry(mtx, jj, ii, 
                               entries[2*iioff], entries[2*iioff+1]) ;
         }
         ioff += istride ;
         istride -= 2 ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   fill A2 object with (1,2) block

   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_fill12block (
   Chv   *chv,
   A2    *mtx
) {
double   *entries ;
int      nD, nL, nU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || mtx == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_fill12block(%p,%p)"
           "\n bad input\n", chv, mtx) ;
   exit(-1) ;
}
if ( ! (CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_fill12block(%p,%p)"
           "\n type = %d, must be real or complex\n", 
           chv, mtx, chv->type) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
/*
   ---------------------------------
   resize the A2 object as necessary
   ---------------------------------
*/
if ( CHV_IS_REAL(chv) ) {
   A2_init(mtx, SPOOLES_REAL, nD, nU, 1, nD, NULL) ;
   A2_zero(mtx) ;
   if ( CHV_IS_SYMMETRIC(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
/*
      ------------------
      chv is symmetric
      ------------------
*/
      ioff = 0 ;
      istride = nD + nU ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
         for ( jj = 0, iioff = ioff + nD - ii ; 
               jj < nU ; 
               jj++, iioff++ ) {
            A2_setRealEntry(mtx, ii, jj, entries[iioff]) ;
         }
         ioff += istride ;
         istride-- ;
      }
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
/*
      ---------------------
      chv is nonsymmetric
      ---------------------
*/
      ioff = nD + nL - 1 ;
      istride = 2*nD + nU + nL - 2 ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
         for ( jj = 0, iioff = ioff + nD - ii ; 
               jj < nU ; 
               jj++, iioff++ ) {
            A2_setRealEntry(mtx, ii, jj, entries[iioff]) ;
         }
         ioff += istride ;
         istride -= 2 ;
      }
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   A2_init(mtx, SPOOLES_COMPLEX, nD, nU, 1, nD, NULL) ;
   A2_zero(mtx) ;
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
/*
      ------------------
      chv is symmetric
      ------------------
*/
      ioff = 0 ;
      istride = nD + nU ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
         for ( jj = 0, iioff = ioff + nD - ii ; 
               jj < nU ; 
               jj++, iioff++ ) {
            A2_setComplexEntry(mtx, ii, jj, 
                               entries[2*iioff], entries[2*iioff+1]) ;
#if MYDEBUG > 0
            fprintf(stdout, 
   "\n 21: ii %d, jj %d, iioff %d, setting (%d,%d) = (%20.12e,%20.12e)",
   ii, jj, iioff, ii, ii, entries[2*iioff], entries[2*iioff+1]) ;
#endif
         }
         ioff += istride ;
         istride-- ;
      }
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      int   ii, iioff, ioff, istride, jj ;
/*
      ---------------------
      chv is nonsymmetric
      ---------------------
*/
      ioff = nD + nL - 1 ;
      istride = 2*nD + nU + nL - 2 ;
      for ( ii = 0 ; ii < nD ; ii++ ) {
         for ( jj = 0, iioff = ioff + nD - ii ; 
               jj < nU ; 
               jj++, iioff++ ) {
#if MYDEBUG > 0
            fprintf(stdout, 
   "\n 21: ii %d, jj %d, iioff %d, setting (%d,%d) = (%20.12e,%20.12e)",
   ii, jj, iioff, ii, ii, entries[2*iioff], entries[2*iioff+1]) ;
#endif
            A2_setComplexEntry(mtx, ii, jj, 
                               entries[2*iioff], entries[2*iioff+1]) ;
         }
         ioff += istride ;
         istride -= 2 ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   fill A2 object with (2,1) block

   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_fill21block (
   Chv   *chv,
   A2    *mtx
) {
double   *entries ;
int      nD, nL, nU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || mtx == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_fillReal21block(%p,%p)"
           "\n bad input\n", chv, mtx) ;
   exit(-1) ;
}
if ( ! (CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_fill21block(%p,%p)"
           "\n type = %d, must be real or complex\n", 
           chv, mtx, chv->type) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
if ( CHV_IS_REAL(chv) ) {
   int   ii, iioff, ioff, istride, jj ;

   A2_init(mtx, SPOOLES_REAL, nL, nD, nD, 1, NULL) ;
   A2_zero(mtx) ;
   ioff = nD + nL - 1 ;
   istride = 2*nD + nU + nL - 2 ;
   for ( ii = 0 ; ii < nD ; ii++ ) {
      for ( jj = 0, iioff = ioff - nD + ii ; 
            jj < nL ; 
            jj++, iioff-- ) {
         A2_setRealEntry(mtx, jj, ii, entries[iioff]) ;
      }
      ioff += istride ;
      istride -= 2 ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   int   ii, iioff, ioff, istride, jj ;

   A2_init(mtx, SPOOLES_COMPLEX, nL, nD, nD, 1, NULL) ;
   A2_zero(mtx) ;
   ioff = nD + nL - 1 ;
   istride = 2*nD + nU + nL - 2 ;
   for ( ii = 0 ; ii < nD ; ii++ ) {
      for ( jj = 0, iioff = ioff - nD + ii ; 
            jj < nL ; 
            jj++, iioff-- ) {
#if MYDEBUG > 0
         fprintf(stdout, 
   "\n 12: ii %d, jj %d, iioff %d, setting (%d,%d) = (%20.12e,%20.12e)",
          ii, jj, iioff, ii, ii, entries[2*iioff], entries[2*iioff+1]) ;
#endif
         A2_setComplexEntry(mtx, jj, ii, 
                            entries[2*iioff], entries[2*iioff+1]) ;
      }
      ioff += istride ;
      istride -= 2 ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
