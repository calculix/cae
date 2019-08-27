/*  assemble.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   add a scaled multiple of a simple chevron to a Chv object.
   the indices are offsets.
   note: for this purpose, (assembling original entries into the
   matrix), the row and column indices of the chevron are identical.
   also, the indices of both the Chv object and the chvind[]
   vector are assumed to be in ascending order.

   created -- 98apr30, cca
   -----------------------------------------------------------------
*/
void
Chv_addChevron (
   Chv      *chv,
   double   alpha[],
   int      ichv,
   int      chvsize,
   int      chvind[],
   double   chvent[]
) {
int      ii, iloc, jcol, jj, jjfirst, jjlast, 
         ncol, nD, nL, nU, offset ;
int      *colind ;
double   *diag ;
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL || ichv < 0 || chvsize < 0 
   || chvind == NULL || chvent == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_addChevron(%p,%p,%d,%d,%p,%p)"
           "\n bad input\n", 
           chv, alpha, ichv, chvsize, chvind, chvent) ;
   exit(-1) ;
}
switch ( chv->type ) {
case SPOOLES_REAL :
   switch ( chv->symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, "\n fatal error in Chv_addChevron()"
              "\n type is SPOOLES_REAL, symflag = %d"
              "\n must be SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n",
              chv->symflag) ;
      exit(-1) ;
      break ;
   }
case SPOOLES_COMPLEX :
   switch ( chv->symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
   case SPOOLES_NONSYMMETRIC :
      break ;
   default :
      fprintf(stderr, "\n fatal error in Chv_addChevron()"
        "\n type is SPOOLES_REAL, symflag = %d"
        "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
        "\n or SPOOLES_NONSYMMETRIC\n",
        chv->symflag) ;
      exit(-1) ;
      break ;
   }
   break ;
default :
   fprintf(stderr, "\n fatal error in Chv_addChevron()"
     "\n type is %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
     chv->type) ;
   exit(-1) ;
   break ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n alpha = %f, ichv = %d, chvsize = %d", 
        alpha, ichv, chvsize) ;
#endif
if ( chvsize == 0 
    || (CHV_IS_REAL(chv) && alpha[0] == 0.0)
    || (CHV_IS_COMPLEX(chv) && (alpha[0] == 0.0 && alpha[1] == 0.0)) ) {
/*
   ----------------------------
   quick return, nothing to add
   ----------------------------
*/
   return ;
}
#if MYDEBUG > 0 
fprintf(stdout, "\n\n Chv_addChevron(%d): ", ichv) ;
IVfprintf(stdout, chvsize, chvind) ;
DVfprintf(stdout, chvsize, chvent) ;
fflush(stdout) ;
#endif
Chv_dimensions(chv, &nD, &nL, &nU) ;
Chv_columnIndices(chv, &ncol, &colind) ;
/*
   -------------------------------------
   locate the chevron in the Chv object 
   that will accumulate these entries
   -------------------------------------
*/
for ( iloc = 0 ; iloc < nD ; iloc++ ) {
   if ( colind[iloc] == ichv ) {
      break ;
   }
}
if ( iloc == nD ) {
/*
   --------------------------------------------
   unable to assemble these entries, error exit
   --------------------------------------------
*/
   fprintf(stderr, "\n fatal error in Chv_addChevron(%p,%d,%d,%p,%p)"
           "\n chevron id %d not found in colind[]",
           chv, ichv, chvsize, chvind, chvent, ichv) ;
   exit(-1) ;
}
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
/*
   ---------------------
   symmetric chevron
   get the local indices
   ---------------------
*/
   jjlast = nD + nU - 1 ;
   for ( ii = 0, jj = iloc ; ii < chvsize ; ii++ ) {
      if ( (offset = chvind[ii]) < 0 ) {
         fprintf(stderr, 
                 "\n fatal error in Chv_addChevron(%p,%d,%d,%p,%p)"
                 "\n ii %d, negative offset %d\n", 
                 chv, ichv, chvsize, chvind, chvent, ii, chvind[ii]) ;
         IVfprintf(stderr, chvsize, chvind) ;
         exit(-1) ;
      }
      jcol = ichv + offset ;
#if MYDEBUG > 0 
      fprintf(stdout, "\n ii = %d, offset = %d, jcol = %d",
              ii, offset, jcol) ;
      fflush(stdout) ;
#endif
      while ( jj <= jjlast && jcol != colind[jj] ) {
         jj++ ;
      }
#if MYDEBUG > 0 
      fprintf(stdout, ", jj = %d", jj) ;
      fflush(stdout) ;
#endif
      if ( jj > jjlast ) {
         fprintf(stderr, 
                 "\n fatal error in Chv_addChevron(%p,%d,%d,%p,%p)"
                 "\n jcol %d not found in colind[]\n", 
                 chv, ichv, chvsize, chvind, chvent, jcol) ;
         fprintf(stderr, "\n colind") ;
         IVfprintf(stderr, ncol, colind) ;
         fprintf(stderr, "\n chvind") ;
         IVfprintf(stderr, chvsize, chvind) ;
         exit(-1) ;
      }
      chvind[ii] = jj ;
   }
#if MYDEBUG > 0 
   fprintf(stdout, "\n local indices") ;
   IVfprintf(stdout, chvsize, chvind) ;
   fflush(stdout) ;
#endif
/*
   --------------------
   assemble the chevron
   --------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      diag = Chv_diagLocation(chv, iloc) - iloc ;
#if MYDEBUG > 0 
      fprintf(stdout, "\n ichv = %d, iloc = %d, diag = %p" 
              "\n chv->entries = %p, diag - chv->entries = %d",
              ichv, iloc, diag, chv->entries,
              diag - chv->entries) ;
      fflush(stdout) ;
#endif
      if ( alpha[0] == 1.0 ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
#if MYDEBUG > 0 
         fprintf(stdout, "\n location %d", 
                 &diag[chvind[ii]] - chv->entries) ;
         fflush(stdout) ;
#endif
            diag[chvind[ii]] += chvent[ii] ;
         }
      } else {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            diag[chvind[ii]] += alpha[0]*chvent[ii] ;
         }
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      diag = Chv_diagLocation(chv, iloc) - 2*iloc ;
#if MYDEBUG > 0 
      fprintf(stdout, "\n ichv = %d, iloc = %d, diag = %p" 
              "\n chv->entries = %p, diag - chv->entries = %d",
              ichv, iloc, diag, chv->entries,
              diag - chv->entries) ;
      fflush(stdout) ;
#endif
      if ( alpha[0] == 1.0 && alpha[1] == 0.0 ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
#if MYDEBUG > 0 
         fprintf(stdout, "\n location %d", 
                 &diag[chvind[ii]] - chv->entries) ;
         fflush(stdout) ;
#endif
            diag[2*chvind[ii]]   += chvent[2*ii] ;
            diag[2*chvind[ii]+1] += chvent[2*ii+1] ;
         }
      } else if ( alpha[0] != 0.0 && alpha[1] == 0.0 ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            diag[2*chvind[ii]]   += alpha[0]*chvent[2*ii] ;
            diag[2*chvind[ii]+1] += alpha[0]*chvent[2*ii+1] ;
         }
      } else if ( CHV_IS_SYMMETRIC(chv) ) {
         double   alphareal, alphaimag, xreal, ximag ; 
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            alphareal = alpha[0] ;     alphaimag = alpha[1] ;
            xreal     = chvent[2*ii] ; ximag     = chvent[2*ii+1] ; 
            diag[2*chvind[ii]]   += alphareal*xreal - alphaimag*ximag ;
            diag[2*chvind[ii]+1] += alphareal*ximag + alphaimag*xreal ;
         }
      } else {
         fprintf(stderr, "\n fatal error in Chv_addChevron()"
    "\n chevron is hermitian, but the scalar has nonzero imaginary part"
    "\n sum is no longer hermitian\n") ;
         exit(-1) ;
      }
   }
/*
   ------------------------------------
   restore the indices to their offsets
   ------------------------------------
*/
   for ( ii = 0 ; ii < chvsize ; ii++ ) {
      chvind[ii] = colind[chvind[ii]] - ichv ;
   }
#if MYDEBUG > 0 
   fprintf(stdout, "\n restored indices") ;
   IVfprintf(stdout, chvsize, chvind) ;
   fflush(stdout) ;
#endif
} else if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
   -----------------------------------------
   nonsymmetric chevron, symmetric structure
   overwrite chvind[] with local indices
   -----------------------------------------
*/
   jjfirst = iloc ;
   jjlast  = nD + nU - 1 ;
   for ( ii = 0, jj = jjlast ; ii < chvsize ; ii++ ) {
      if ( (offset = chvind[ii]) >= 0 ) {
         break ;
      } 
      jcol = ichv - offset ;
      while ( jj >= jjfirst && jcol != colind[jj] ) {
         jj-- ;
      }
      if ( jj < jjfirst ) {
         fprintf(stderr, 
                 "\n fatal error in Chv_addChevron(%p,%d,%d,%p,%p)"
                 "\n jcol %d not found in colind[]\n", 
                 chv, ichv, chvsize, chvind, chvent, jcol) ;
         exit(-1) ;
      }
      chvind[ii] = -jj + iloc ;
   }
   for ( jj = jjfirst ; ii < chvsize ; ii++ ) {
      jcol = ichv + chvind[ii] ;
      while ( jj <= jjlast && jcol != colind[jj] ) {
         jj++ ;
      }
      if ( jj > jjlast ) {
         fprintf(stderr, 
                 "\n fatal error in Chv_addChevron(%p,%d,%d,%p,%p)"
                 "\n jcol %d not found in colind[]\n", 
                 chv, ichv, chvsize, chvind, chvent, jcol) ;
         exit(-1) ;
      }
      chvind[ii] = jj - iloc ;
   }
#if MYDEBUG > 0 
   fprintf(stdout, "\n local indices") ;
   IVfprintf(stdout, chvsize, chvind) ;
   fflush(stdout) ;
#endif
/*
   --------------------
   assemble the chevron
   --------------------
*/
   diag = Chv_diagLocation(chv, iloc) ;
#if MYDEBUG > 0 
   fprintf(stdout, "\n ichv = %d, iloc = %d, diag = %p" 
           "\n chv->entries = %p, diag - chv->entries = %d",
           ichv, iloc, diag, chv->entries,
           diag - chv->entries) ;
   fflush(stdout) ;
#endif
   if ( CHV_IS_REAL(chv) ) {
      if ( alpha[0] == 1.0 ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
#if MYDEBUG > 0 
            fprintf(stdout, "\n  diag[%d] += %12.4e, ii = %d", 
                    chvind[ii], chvent[ii], ii) ;
            fflush(stdout) ;
#endif
            diag[chvind[ii]] += chvent[ii] ;
         }
      } else {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
#if MYDEBUG > 0 
            fprintf(stdout, "\n  diag[%d] += %12.4e, ii = %d", 
                    chvind[ii], chvent[ii], ii) ;
            fflush(stdout) ;
#endif
            diag[chvind[ii]] += alpha[0] * chvent[ii] ;
         }
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      if ( alpha[0] == 1.0 && alpha[1] == 0.0 ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
#if MYDEBUG > 0 
            fprintf(stdout, "\n  diag[%d] += %12.4e, ii = %d", 
                    chvind[ii], chvent[ii], ii) ;
            fflush(stdout) ;
#endif
            diag[2*chvind[ii]]   += chvent[2*ii] ;
            diag[2*chvind[ii]+1] += chvent[2*ii+1] ;
         }
      } else if ( alpha[0] != 1.0 && alpha[1] == 0.0 ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            diag[2*chvind[ii]]   += alpha[0] * chvent[2*ii] ;
            diag[2*chvind[ii]+1] += alpha[0] * chvent[2*ii+1] ;
         }
      } else {
         double   alphareal, alphaimag, xreal, ximag ; 
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            alphareal = alpha[0]     ; alphaimag = alpha[1] ;
            xreal     = chvent[2*ii] ; ximag     = chvent[2*ii+1] ; 
            diag[2*chvind[ii]]   += alphareal*xreal - alphaimag*ximag ;
            diag[2*chvind[ii]+1] += alphareal*ximag + alphaimag*xreal ;
         }
      }
   }
/*
   ------------------------------------
   restore the indices to their offsets
   ------------------------------------
*/
   for ( ii = 0 ; ii < chvsize ; ii++ ) {
      if ( chvind[ii] < 0 ) {
         chvind[ii] = ichv - colind[iloc - chvind[ii]] ;
      } else {
         chvind[ii] = colind[chvind[ii] + iloc] - ichv ;
      }
   }
#if MYDEBUG > 0 
   fprintf(stdout, "\n restored indices") ;
   IVfprintf(stdout, chvsize, chvind) ;
   fflush(stdout) ;
#endif
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   assemble Chv object chvI into Chv object chvJ.
   note: the two objects must be of the same symmetry type,
         the row indices of chvI must nest into those of chvJ,
         the column indices of chvI must nest into those of chvJ.

   created -- 98apr30, cca
   --------------------------------------------------------------
*/
void
Chv_assembleChv (
   Chv   *chvJ,
   Chv   *chvI
) {
double   *diagI, *diagJ ;
int      ii, ichvI, ichvJ, jj, ncolI, ncolJ, nDI, nDJ, 
         nLI, nLJ, nrowI, nrowJ, nUI, nUJ, offset ;
int      *colindJ, *colindI, *rowindI, *rowindJ ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvJ == NULL || chvI == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_assembleChv(%p,%p)"
           "\n bad input\n", chvJ, chvI) ;
   exit(-1) ;
}
if ( !(CHV_IS_SYMMETRIC(chvI) 
       || CHV_IS_HERMITIAN(chvI) 
       || CHV_IS_NONSYMMETRIC(chvI) ) ) {
   fprintf(stderr, "\n fatal error in Chv_assembleChv(%p,%p)"
           "\n bad symflag %d\n", chvJ, chvI, chvI->symflag) ;
   exit(-1) ;
}
if ( chvI->symflag != chvJ->symflag ) {
   fprintf(stderr, "\n fatal error in Chv_assembleChv(%p,%p)"
           "\n chvI->symflag = %d, chvJ->symflag = %d\n", 
           chvJ, chvI, chvI->symflag, chvJ->symflag) ;
   exit(-1) ;
}
/*
   -------------------------------------
   get the dimensions of the two objects
   -------------------------------------
*/
Chv_dimensions(chvJ, &nDJ, &nLJ, &nUJ) ;
Chv_dimensions(chvI, &nDI, &nLI, &nUI) ;
if (  nDI + nLI > nDJ + nLJ ||  nDI + nUI > nDJ + nUJ ) {
   fprintf(stderr, "\n fatal error in Chv_assembleChv(%p,%p)"
           "\n bad dimensions"
           "\n nDI = %d, nLI = %d, nUI = %d"
           "\n nDI = %d, nLI = %d, nUI = %d",
           chvJ, chvI, nDI, nLI, nUI, nDJ, nLJ, nUJ) ;
   exit(-1) ;
}
/*
   -----------------
   get the local ids
   -----------------
*/
Chv_columnIndices(chvJ, &ncolJ, &colindJ) ;
Chv_columnIndices(chvI, &ncolI, &colindI) ;
#if MYDEBUG > 0
fprintf(stdout, "\n colindI") ;
IVfprintf(stdout, ncolI, colindI) ;
fprintf(stdout, "\n colindJ") ;
IVfprintf(stdout, ncolJ, colindJ) ;
#endif
for ( ii = 0, jj = 0 ; ii < ncolI ; ii++ ) {
   while ( jj < ncolJ && colindI[ii] != colindJ[jj] ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n colindI[%d] = %d, colindJ[%d] = %d",
              ii, colindI[ii], jj, colindJ[jj]) ;
#endif
      jj++ ;
   }
   if ( jj == ncolJ ) {
      break ;
   }
   colindI[ii] = jj ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n local column indices") ;
IVfprintf(stdout, ncolI, colindI) ;
#endif
if ( jj == ncolJ ) {
   fprintf(stderr, "\n fatal error in Chv_assembleChv(%p,%p)"
           "\n column indicesI do not nest in indicesJ\n", chvJ, chvI) ;
   fprintf(stderr, "\n colindI") ;
   IVfprintf(stderr, ncolI, colindI) ;
   fprintf(stderr, "\n colindJ") ;
   IVfprintf(stderr, ncolJ, colindJ) ;
   exit(-1) ;
}
if ( CHV_IS_SYMMETRIC(chvJ) || CHV_IS_HERMITIAN(chvJ) ) {
/*
   -------------------
   symmetric structure
   -------------------
*/
   nrowI   = ncolI   ;
   rowindI = colindI ;
} else if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
/*
   ----------------------
   nonsymmetric structure
   ----------------------
*/
   Chv_rowIndices(chvJ, &nrowJ, &rowindJ) ;
   Chv_rowIndices(chvI, &nrowI, &rowindI) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n rowindI") ;
   IVfprintf(stdout, nrowI, rowindI) ;
   fprintf(stdout, "\n rowindJ") ;
   IVfprintf(stdout, nrowJ, rowindJ) ;
#endif
   for ( ii = 0, jj = 0 ; ii < nrowI ; ii++ ) {
      while ( jj < nrowJ && rowindI[ii] != rowindJ[jj] ) {
         jj++ ;
      }
      if ( jj == nrowJ ) {
         break ;
      }
      rowindI[ii] = jj ;
   }
   if ( jj == nrowJ ) {
      fprintf(stderr, "\n fatal error in Chv_assembleChv(%p,%p)"
              "\n row indicesI do not nest in indicesJ\n", chvJ, chvI) ;
      fprintf(stderr, "\n rowindI") ;
      IVfprintf(stderr, nrowI, rowindI) ;
      fprintf(stderr, "\n rowindJ") ;
      IVfprintf(stderr, nrowJ, rowindJ) ;
      exit(-1) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n local row indices") ;
   IVfprintf(stdout, nrowI, rowindI) ;
#endif
}
#if MYDEBUG > 0
fprintf(stdout, "\n local column indices") ;
IVfprintf(stdout, ncolI, colindI) ;
fprintf(stdout, "\n local row indices") ;
IVfprintf(stdout, nrowI, rowindI) ;
#endif
/*
   ---------------------------
   loop over the chevrons in I
   ---------------------------
*/
for ( ichvI = 0 ; ichvI < nDI ; ichvI++ ) {
   ichvJ = colindI[ichvI] ;
   if ( ichvJ != rowindI[ichvI] ) {
      fprintf(stderr, "\n fatal error in Chv_assembleChv(%p,%p)"
              "\n ichvI = %d, ichvJ = %d, rowindI[ichvI] = %d",
              chvJ, chvI, ichvI, ichvJ, rowindI[ichvI]) ;
      exit(-1) ;
   }
   diagI = Chv_diagLocation(chvI, ichvI) ;
   diagJ = Chv_diagLocation(chvJ, ichvJ) ;
#if MYDEBUG > 0
   fprintf(stdout, 
          "\n ichvI = %d, diagI - entriesI = %d"
          "\n ichvJ = %d, diagJ - entriesJ = %d",
          ichvI, diagI - chvI->entries, ichvJ, diagJ - chvJ->entries) ;
#endif
   if ( CHV_IS_REAL(chvJ) ) {
      diagJ[0] += diagI[0] ;
   } else if ( CHV_IS_COMPLEX(chvJ) ) {
      diagJ[0] += diagI[0] ;
      diagJ[1] += diagI[1] ;
   }
   if ( CHV_IS_SYMMETRIC(chvJ) || CHV_IS_HERMITIAN(chvJ) ) {
/*
      ----------------------
      symmetric or hermitian
      ----------------------
*/
      if ( CHV_IS_REAL(chvJ) ) {
         for ( ii = ichvI + 1, jj = 1 ; ii < ncolI ; ii++, jj++ ) {
            offset = colindI[ii] - ichvJ ;
#if MYDEBUG > 0
            fprintf(stdout, 
                   "\n ii = %d, jj = %d, offset = %d", ii, jj, offset) ;
#endif
            diagJ[offset] += diagI[jj] ;
         }
      } else if ( CHV_IS_COMPLEX(chvJ) ) {
         for ( ii = ichvI + 1, jj = 1 ; ii < ncolI ; ii++, jj++ ) {
            offset = colindI[ii] - ichvJ ;
#if MYDEBUG > 0
            fprintf(stdout, 
                   "\n ii = %d, jj = %d, offset = %d", ii, jj, offset) ;
#endif
            diagJ[2*offset]   += diagI[2*jj] ;
            diagJ[2*offset+1] += diagI[2*jj+1] ;
         }
      }
   } else if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
/*
      ------------
      nonsymmetric
      ------------
*/
      if ( CHV_IS_REAL(chvJ) ) {
         for ( ii = ichvI + 1, jj = 1 ; ii < ncolI ; ii++, jj++ ) {
            offset = colindI[ii] - ichvJ ;
#if MYDEBUG > 0
            fprintf(stdout, "\n diagJ[%d] += diagI[%d]", offset, jj) ;
#endif
            diagJ[offset] += diagI[jj] ;
         }
         for ( ii = ichvI + 1, jj = -1 ; ii < nrowI ; ii++, jj-- ) {
            offset = rowindI[ii] - ichvJ ;
#if MYDEBUG > 0
         fprintf(stdout, "\n diagJ[%d] += diagI[%d]", -offset, jj) ;
#endif
            diagJ[-offset] += diagI[jj] ;
         }
      } else if ( CHV_IS_COMPLEX(chvJ) ) {
         for ( ii = ichvI + 1, jj = 1 ; ii < ncolI ; ii++, jj++ ) {
            offset = colindI[ii] - ichvJ ;
#if MYDEBUG > 0
            fprintf(stdout, "\n diagJ[%d] += diagI[%d]", offset, jj) ;
#endif
            diagJ[2*offset]   += diagI[2*jj] ;
            diagJ[2*offset+1] += diagI[2*jj+1] ;
         }
         for ( ii = ichvI + 1, jj = -1 ; ii < nrowI ; ii++, jj-- ) {
            offset = rowindI[ii] - ichvJ ;
#if MYDEBUG > 0
         fprintf(stdout, "\n diagJ[%d] += diagI[%d]", -offset, jj) ;
#endif
            diagJ[-2*offset]   += diagI[2*jj] ;
            diagJ[-2*offset+1] += diagI[2*jj+1] ;
         }
      }
   }
}
/*
   -------------------
   restore the indices
   -------------------
*/
for ( ii = 0 ; ii < ncolI ; ii++ ) {
   colindI[ii] = colindJ[colindI[ii]] ;
}
if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
   for ( ii = 0 ; ii < nrowI ; ii++ ) {
      rowindI[ii] = rowindJ[rowindI[ii]] ;
   }
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- assemble the postponed data from the children


   newchv     -- Chv object to contain fully assembled front
   oldchv     -- Chv object that contains former front
   firstchild -- pointer to first child in the list of children
                 Chv objects to be merged into the new front

   return value -- # of delayed rows and columns added to the front

   created -- 98apr30, cca
   ----------------------------------------------------------------
*/
int
Chv_assemblePostponedData (
   Chv   *newchv,
   Chv   *oldchv,
   Chv   *firstchild
) {
Chv     *child ;
int      ierr, ndelay ;
/*
   ---------------
   check the input
   ---------------
*/
if ( newchv == NULL || oldchv == NULL || firstchild == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_assemblePostponedData(%p,%p,%p)"
           "\n bad input\n", newchv, oldchv, firstchild) ;
   exit(-1) ;
}
Chv_zero(newchv) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n inside Chv_assemblePostponedData, J = %d", 
        oldchv->id) ;
fprintf(stdout, "\n\n old chevron %d", oldchv->id) ;
Chv_writeForHumanEye(oldchv, stdout) ;
for ( child = firstchild ; child != NULL ; child = child->next ) {
   fprintf(stdout, "\n\n child %d", child->id) ;
   Chv_writeForHumanEye(child, stdout) ;
}
fprintf(stdout, "\n\n new chevron %d", newchv->id) ;
Chv_writeForHumanEye(newchv, stdout) ;
fflush(stdout) ;
#endif
/*
   ------------------------------------
   copy the indices into the new object
   ------------------------------------
*/
ndelay = 0 ;
for ( child = firstchild ; child != NULL ; child = child->next ) {
   IVcopy(child->nD, newchv->colind + ndelay, child->colind) ;
   ndelay += child->nD ;
}
IVcopy(oldchv->nD + oldchv->nU, 
       newchv->colind + ndelay, oldchv->colind) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ndelay = %d, new column indices", ndelay) ;
IVfp80(stdout, newchv->nD + newchv->nU, newchv->colind, 80, &ierr) ;
fflush(stdout) ;
#endif
if ( CHV_IS_NONSYMMETRIC(newchv) ) {
   ndelay = 0 ;
   for ( child = firstchild ; child != NULL ; child = child->next ) {
      IVcopy(child->nD, newchv->rowind + ndelay, child->rowind) ;
      ndelay += child->nD ;
   }
   IVcopy(oldchv->nD + oldchv->nU, 
          newchv->rowind + ndelay, oldchv->rowind) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n ndelay = %d, new row indices", ndelay) ;
   IVfp80(stdout, newchv->nD + newchv->nL, newchv->rowind, 80, &ierr) ;
   fflush(stdout) ;
#endif
}
/*
   ------------------------
   assemble the old chevron
   ------------------------
*/
Chv_assembleChv(newchv, oldchv) ;
#if MYDEBUG > 0
fprintf(stdout, "\n after merging old chevron") ;
Chv_writeForHumanEye(newchv, stdout) ;
#endif
/*
   ------------------------------
   assemble the children chevrons
   ------------------------------
*/
for ( child = firstchild ; child != NULL ; child = child->next ) {
   Chv_assembleChv(newchv, child) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n after merging child %d", child->id) ;
   Chv_writeForHumanEye(newchv, stdout) ;
#endif
}

return(ndelay) ; }

/*--------------------------------------------------------------------*/
