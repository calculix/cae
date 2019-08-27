/*  factor.c  */

#include "../Chv.h"

#define MYDEBUG 0
#define MYCHECK 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to factor the front without pivoting

   return value -- # of eliminated rows and columns

   created -- 98aug27, cca
   ------------------------------------------------
*/
int
Chv_factorWithNoPivoting (
   Chv              *chv,
   PatchAndGoInfo   *info
) {
Chv   wrkChv ;
int   ncol, nD, nelim, nrow ;
int   *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_factorWithNoPivoting()"
           "\n bad input\n") ;
   exit(-1) ;
}
nD = chv->nD ;
/*
   --------------------------
   set up the working chevron
   --------------------------
*/
Chv_setDefaultFields(&wrkChv) ;
Chv_rowIndices(chv, &nrow, &rowind) ;
Chv_columnIndices(chv, &ncol, &colind) ;
Chv_initWithPointers(&wrkChv, chv->id, nD, chv->nL, chv->nU, chv->type, 
                     chv->symflag, rowind, colind, Chv_entries(chv)) ;
/*
   -------------------------------------
   switch over the patch-and-go strategy
   -------------------------------------
*/
if ( info == NULL ) {
/*
   ------------------
   simple no pivoting
   ------------------
*/
   nelim = 0 ;
   while ( nelim < nD ) {
      if ( Chv_r1upd(&wrkChv) == 0 ) {
         break ;
      }
      Chv_shift(&wrkChv, 1) ;
      nelim++ ;
   }
} else if ( info->strategy == 1 ) {
   double   colmaxabs, diagabs, offmaxabs, rowmaxabs ;
/*
   ----------------------------------------
   Patch-and-go for optimization matrices.
   if |diag| < toosmall * max|offdiag| then
      diag = 1.0
      offdiag = 0.0
   endif
   ----------------------------------------
*/ 
   for ( nelim = 0 ; nelim < nD ; nelim++ ) {
      Chv_maxabsInChevron(&wrkChv, 0, &diagabs, &rowmaxabs, &colmaxabs);
      offmaxabs = (rowmaxabs >= colmaxabs) ? rowmaxabs : colmaxabs ;
      if ( diagabs <= info->toosmall * offmaxabs ) {
         if ( CHV_IS_REAL(chv) ) {
            wrkChv.entries[0] = 1.0 ;
         } else {
            wrkChv.entries[0] = 1.0 ;
            wrkChv.entries[1] = 0.0 ;
         }
         Chv_zeroOffdiagonalOfChevron(chv, 0) ;
         if ( info->fudgeIV != NULL ) {
            IV_push(info->fudgeIV, chv->colind[0]) ;
         }
      }
      Chv_r1upd(&wrkChv) ;
      Chv_shift(&wrkChv, 1) ;
   }
} else if ( info->strategy == 2 ) {
   double   colmaxabs, diagabs, olddiag, newdiag, offmaxabs, rowmaxabs ;
/*
   ----------------------------------------------
   Patch-and-go for structural analysis matrices.
   if |diag| < fudge then
      diag = fudge * max(max|offdiag|, 1.0)
   endif
   ----------------------------------------------
*/ 
   for ( nelim = 0 ; nelim < nD ; nelim++ ) {
      Chv_maxabsInChevron(&wrkChv, 0, &diagabs, &rowmaxabs, &colmaxabs);
      offmaxabs = (rowmaxabs >= colmaxabs) ? rowmaxabs : colmaxabs ;
      if ( diagabs <= info->fudge ) {
         olddiag = diagabs ;
         if ( offmaxabs < 1.0 ) {
            offmaxabs = 1.0 ;
         }
         if ( CHV_IS_REAL(chv) ) {
            wrkChv.entries[0] = newdiag = info->fudge * offmaxabs ;
         } else {
            wrkChv.entries[0] = newdiag = info->fudge * offmaxabs ;
            wrkChv.entries[1] = 0.0 ;
         }
         if ( info->fudgeIV != NULL ) {
            IV_push(info->fudgeIV, chv->colind[0]) ;
         }
         if ( info->fudgeDV != NULL ) {
            DV_push(info->fudgeDV, fabs(olddiag - newdiag)) ;
         }
      }
      Chv_r1upd(&wrkChv) ;
      Chv_shift(&wrkChv, 1) ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_factorWithNoPivoting()"
           "\n patch-and-go info != NULL, strategy = %d",
           info->strategy) ;
   exit(-1) ;
}
return(nelim) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- factor the pivot chevron with pivoting

   ndelay -- number of delayed rows and columns
   pivotflag -- enable pivoting or not
      0 --> no pivoting
      1 --> enable pivoting
   pivotsizesIV -- IV object that holds the sizes of the pivots,
      used only when the front is symmetric or hermitian
      and pivoting is enabled
   workDV -- DV object used for working storage, resized as necessary
   tau    -- upper bound on the magnitude of the entries 
      in the factors, used only when pivoting is enabled
   pntest -- pointer to be incremented with the number of pivot tests

   return value -- # of eliminated rows and columns

   created -- 98aug27, cca
   ------------------------------------------------------------------
*/
int
Chv_factorWithPivoting (
   Chv     *chv,
   int      ndelay,
   int      pivotflag,
   IV       *pivotsizesIV,
   DV       *workDV,
   double   tau,
   int      *pntest
) {
Chv   wrkChv ;
int   irow, jcol, ncol, nD, nelim, nrow, pivotsize ;
int   *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || pivotflag != 1 || ndelay < 0 ) {
   fprintf(stderr, "\n fatal error in Chv_factorWithPivoting()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( workDV == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_factorWithPivoting()"
           "\n workDV is NULL \n") ;
   exit(-1) ;
}
if ( tau < 1.0 ) {
   fprintf(stderr, "\n fatal error in Chv_factorWithPivoting()"
           "\n tau = %f is invalid \n", tau) ; 
   exit(-1) ;
}
if ( CHV_IS_REAL(chv) && CHV_IS_SYMMETRIC(chv)
   && pivotsizesIV == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_factorWithPivoting()"
           "\n real symmetric front"
           "\n pivoting enabled, pivotsizesIV is NULL\n") ;
   exit(-1) ;
}
if ( CHV_IS_COMPLEX(chv) 
   && (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv))
   && pivotsizesIV == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_factorWithPivoting()"
           "\n complex symmetric or hermitian front"
           "\n pivoting enabled, pivotsizesIV is NULL\n") ;
   exit(-1) ;
}
nD = chv->nD ;
/*
   --------------------------
   set up the working chevron
   --------------------------
*/
Chv_setDefaultFields(&wrkChv) ;
Chv_rowIndices(chv, &nrow, &rowind) ;
Chv_columnIndices(chv, &ncol, &colind) ;
Chv_initWithPointers(&wrkChv, chv->id, nD, chv->nL, chv->nU, chv->type,
                      chv->symflag, rowind, colind, Chv_entries(chv)) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n after initializing wrkChv") ;
Chv_writeForHumanEye(&wrkChv, stdout) ;
fflush(stdout) ;
#endif
/*
   -----------------------------
   switch over the symmetry flag
   -----------------------------
*/
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n\n pivoting, hermitian or symmetric front %d", chv->id) ;
   fflush(stdout) ;
#endif
/*
   -------------------------
   symmetric structure front
   -------------------------
*/
   IV_setSize(pivotsizesIV, 0) ;
   nelim = 0 ;
   while ( nelim < nD ) {
/*
      -------------------------
      find the 1x1 or 2x2 pivot
      -------------------------
*/
#if MYDEBUG > 0
      fprintf(stdout, 
            "\n trying to find pivot, nelim = %d, nD = %d, ndelay = %d",
            nelim, nD, ndelay) ;
      fflush(stdout) ;
#endif
      pivotsize = Chv_findPivot(&wrkChv, workDV, tau, ndelay, 
                                &irow, &jcol, pntest) ;
      if ( irow > jcol ) {
         int itemp = irow ;
         irow = jcol ;
         jcol = itemp ;
      }
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n pivotsize = %d, local irow = %d, local jcol = %d",
              pivotsize, irow, jcol) ;
      fflush(stdout) ;
#endif
      irow += nelim ;
      jcol += nelim ;
      if ( pivotsize == 0 ) {
/*
         ---------------------------------
         no pivot found, break out of loop
         ---------------------------------
*/
         break ;
      } else {
         ndelay = 0 ;
         if ( irow == jcol ) {
/*
            ------------------------------------------------------
            1x1 pivot found, swap row and column, update and shift
            ------------------------------------------------------
*/
#if MYDEBUG > 0
            fprintf(stdout, "\n\n before swaps") ;
            Chv_writeForHumanEye(chv, stdout) ;
            fflush(stdout) ;
#endif
            Chv_swapRowsAndColumns(chv, nelim, irow) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n\n after swaps") ;
            Chv_writeForHumanEye(chv, stdout) ;
            fflush(stdout) ;
#endif
#if MYDEBUG > 0
            fprintf(stdout, "\n\n nelim = %d, before update", nelim) ;
            Chv_writeForHumanEye(&wrkChv, stdout) ;
            fflush(stdout) ;
#endif
            if ( Chv_r1upd(&wrkChv) == 0 ) {
               break ;
            }
#if MYDEBUG > 0
            fprintf(stdout, "\n\n nelim = %d, after update", nelim) ;
            Chv_writeForHumanEye(&wrkChv, stdout) ;
            fflush(stdout) ;
#endif
            Chv_shift(&wrkChv, 1) ;
            nelim++ ;
            IV_push(pivotsizesIV, 1) ;
         } else {
/*
            --------------------------------------------------------
            2x2 pivot found, swap rows and columns, update and shift
            --------------------------------------------------------
*/
#if MYDEBUG > 0
            fprintf(stdout, "\n\n before swaps") ;
            Chv_writeForHumanEye(chv, stdout) ;
            fflush(stdout) ;
#endif
            Chv_swapRowsAndColumns(chv, nelim, irow) ;
            Chv_swapRowsAndColumns(chv, nelim+1, jcol) ;
#if MYDEBUG > 0
            fprintf(stdout, "\n\n after swaps") ;
            Chv_writeForHumanEye(chv, stdout) ;
            fflush(stdout) ;
#endif
#if MYDEBUG > 0
            fprintf(stdout, "\n\n irow = %d, jcol = %d", irow, jcol) ;
            fprintf(stdout, "\n\n nelim = %d, before update", nelim) ;
#endif
#if MYDEBUG > 0
            Chv_writeForHumanEye(&wrkChv, stdout) ;
            fflush(stdout) ;
#endif
            if ( Chv_r2upd(&wrkChv) == 0 ) {
               break ;
            }
#if MYDEBUG > 0
            fprintf(stdout, "\n\n nelim = %d, after update", nelim) ;
            Chv_writeForHumanEye(&wrkChv, stdout) ;
            fflush(stdout) ;
#endif
            Chv_shift(&wrkChv, 2) ;
            nelim += 2 ;
            IV_push(pivotsizesIV, 2) ;
         }
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n\n ok, done with this pivot") ;
      fflush(stdout) ;
#endif
   }
} else {
/*
   ------------------
   nonsymmetric front
   ------------------
*/
   nelim = 0 ;
   while ( nelim < nD ) {
/*
      ------------------
      find the 1x1 pivot
      ------------------
*/
      pivotsize = Chv_findPivot(&wrkChv, workDV, tau, ndelay, 
                                 &irow, &jcol, pntest) ;
      irow += nelim ;
      jcol += nelim ;
#if MYDEBUG > 0
      fprintf(stdout, "\n\n irow = %d, jcol = %d", irow, jcol) ;
      fflush(stdout) ;
#endif
      if ( pivotsize == 0 ) {
/*
         ---------------------------------
         no pivot found, break out of loop
         ---------------------------------
*/
         break ;
      } else {
         ndelay = 0 ;
/*
         ------------------------------------------------------
         1x1 pivot found, swap row and column, update and shift
         ------------------------------------------------------
*/
#if MYDEBUG > 1
         fprintf(stdout, "\n\n before swaps") ;
         Chv_writeForHumanEye(chv, stdout) ;
         fflush(stdout) ;
#endif
         Chv_swapRows(chv, nelim, irow) ;
         Chv_swapColumns(chv, nelim, jcol) ;
#if MYDEBUG > 1
         fprintf(stdout, "\n\n after swaps") ;
         Chv_writeForHumanEye(chv, stdout) ;
         fflush(stdout) ;
#endif
#if MYDEBUG > 0
         fprintf(stdout, "\n\n nelim = %d, before update", nelim) ;
         fflush(stdout) ;
#endif
#if MYDEBUG > 1
         Chv_writeForHumanEye(&wrkChv, stdout) ;
         fflush(stdout) ;
#endif
         if ( Chv_r1upd(&wrkChv) == 0 ) {
            break ;
         }
#if MYDEBUG > 0
         fprintf(stdout, "\n\n nelim = %d, after update", nelim) ;
#endif
#if MYDEBUG > 1
         Chv_writeForHumanEye(&wrkChv, stdout) ;
         fflush(stdout) ;
#endif
         Chv_shift(&wrkChv, 1) ;
         nelim++ ;
      }
   }
}
return(nelim) ; }

/*--------------------------------------------------------------------*/
static int symmetric1x1 ( Chv *chv ) ;
static int nonsym1x1 ( Chv *chv) ;
static int symmetric2x2 ( Chv *chv ) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   perform a rank one update using the first row and column.
   this is used in the (L + I)D(I + U) factorization

   return code ---
      0 if the pivot was zero
      1 if the pivot was nonzero

   created -- 98jan23, cca
   ---------------------------------------------------------
*/
int
Chv_r1upd ( 
   Chv   *chv
) {
int   rc = 0 ;

if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_r1upd(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
} 
if ( CHV_IS_NONSYMMETRIC(chv) ) {
   rc = nonsym1x1(chv) ;
} else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   rc = symmetric1x1(chv) ;
} else {
   fprintf(stderr, "\n fatal error in Chv_r1upd(%p)"
           "\n symflag = %d\n", chv, chv->symflag) ;
   exit(-1) ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   perform a rank two update using the first two rows.
   used in the (U^T + I)D(I + U) and (U^H + I)D(I + U) factorizations

   return code ---
      0 if the pivot was zero
      1 if the pivot was nonzero

   created -- 98jan23, cca
   ------------------------------------------------------------------
*/
int
Chv_r2upd ( 
   Chv   *chv
) {
int   rc = 0 ;

if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_r2upd(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
} 
if ( !(CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_r2upd(%p)"
           "\n symflag = %d\n", chv, chv->symflag) ;
   exit(-1) ;
} 
rc = symmetric2x2(chv) ;

return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   perform an internal rank-1 update for a symmetric chevron

   return code ---
      0 if the pivot was zero
      1 if the pivot was nonzero

   created -- 98jan23, cca
   ------------------------------------------------------------
*/
static int
symmetric1x1 (
  Chv   *chv
) {
double   dimag, dreal, fac1, fac2, limag, lreal, uimag, ureal ;
double   *entries ;
int      dloc, dstride, kchv, nD, nL, nU, uijloc, ukbeg, usize ;
/*
   -------------------------------------
   get dimensions and pointer to entries
   -------------------------------------
*/
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
#if MYDEBUG > 0
fprintf(stdout, "\n nD = %d, nL = %d, nU = %d, entries = %p",
        nD, nL, nU, entries) ;
#endif
/*
   ----------------------------------------------
   dloc    : offset to the first diagonal element
   dstride : stride to next diagonal element
   usize   : size of first row in upper part
   ----------------------------------------------
*/
dloc    = 0 ;
dstride = nD + nU ;
usize   = nD + nU - 1 ;
#if MYDEBUG > 0
fprintf(stdout, "\n dloc = %d, dstride = %d, usize = %d",
        dloc, dstride, usize) ;
#endif
/*
   ----------------------
   check for a zero pivot
   ----------------------
*/
if ( CHV_IS_REAL(chv) ) {
   dreal = entries[dloc] ;
   if ( dreal == 0.0 ) {
      return(0) ;
   }
   fac1 = 1./dreal ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   dreal = entries[2*dloc] ;
   dimag = entries[2*dloc+1] ;
   if ( dreal == 0.0 && dimag == 0.0 ) {
      return(0) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n chv->id = %d : (dreal,dimag) = <%12.4e,%12.4e>", 
           chv->id, dreal, dimag) ;
   fprintf(stdout, "\n (dreal,dimag) = <%12.4e,%12.4e>", 
           dreal, dimag) ;
#endif
/*
   ------------------------------
   compute (fac1,fac2) = 1/d(0,0)
   ------------------------------
*/
   if ( CHV_IS_HERMITIAN(chv) ) {
      fac1 = 1./dreal ; fac2 = 0.0 ;
      entries[2*dloc+1] = 0.0 ;
   } else {
      Zrecip(dreal, dimag, &fac1, &fac2) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n 1/(dreal,dimag) = <%12.4e,%12.4e>", 
           fac1, fac2) ;
#endif
}
/*
   ------------------------
   scale the first row of U
   ------------------------
*/
if ( CHV_IS_REAL(chv) ) {
   DVscale(usize, &entries[1], fac1) ;
} else if ( CHV_IS_HERMITIAN(chv) ) {
   DVscale(2*usize, &entries[2], fac1) ;
} else {
   ZVscale(usize, &entries[2], fac1, fac2) ;
}
/*
   ------------------------------------
   loop over the following chevrons
   uijloc -- offset into uij multiplier
   ------------------------------------
*/
uijloc  = dloc + 1 ;
for ( kchv = 1 ; kchv < nD ; kchv++ ) {
/*
   --------------------------------------------------
   dloc now points to next diagonal location
   ukbeg -- offset into start of row in upper part
   --------------------------------------------------
*/
   dloc += dstride ;
   ukbeg = dloc + 1 ;
#if MYDEBUG > 0
   fprintf(stdout, "\n kchv = %5d, dloc = %5d"
           "\n uijloc = %5d, usize = %5d, ukbeg = %5d",
           kchv, dloc, uijloc, usize, ukbeg) ;
#endif
/*
   ------------------------------------
   pull out the multiplier coefficients
   ------------------------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      ureal = entries[uijloc] ;
      lreal =  dreal*ureal ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      ureal = entries[2*uijloc] ;
      uimag = entries[2*uijloc+1] ;
      if (CHV_IS_HERMITIAN(chv) ) {
         lreal =  dreal*ureal ;
         limag = -dreal*uimag ;
      } else {
         lreal = dreal*ureal - dimag*uimag ;
         limag = dreal*uimag + dimag*ureal ;
      }
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n (lreal,limag) = <%12.4e,%12.4e>"
              "\n (ureal,uimag) = <%12.4e,%12.4e>", 
              lreal, limag, ureal, uimag) ;
#endif
   }
/*
   -------------------------------------------------
   update the upper row including the diagonal entry
   -------------------------------------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      DVaxpy(usize, &entries[dloc], -lreal, &entries[uijloc]) ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      ZVaxpy(usize, &entries[2*dloc], -lreal, -limag, 
             &entries[2*uijloc]) ;
   }
/*
   ----------------------------------
   adjust offsets and diagonal stride
   ----------------------------------
*/
   uijloc++ ; dstride-- ; usize-- ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   perform an internal rank-1 update for a nonsymmetric chevron

   return code ---
      0 if the pivot was zero
      1 if the pivot was nonzero

   created -- 98jan23, cca
   ------------------------------------------------------------
*/
static int
nonsym1x1 (
  Chv   *chv
) {
double   dimag, dreal, fac1, fac2, limag, lreal, uimag, ureal ;
double   *entries ;
int      dloc, dstride, kchv, ljiloc, lkbeg, lsize, nD, nL, nU,
         uijloc, ukbeg, usize ;
/*
   -------------------------------------
   get dimensions and pointer to entries
   -------------------------------------
*/
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
#if MYDEBUG > 0
fprintf(stdout, "\n nD = %d, nL = %d, nU = %d, entries = %p",
        nD, nL, nU, entries) ;
#endif
/*
   ----------------------------------------------
   dloc    : offset to the first diagonal element
   dstride : stride to next diagonal element
   lsize   : size of first column in lower part
   usize   : size of first row in upper part
   ----------------------------------------------
*/
dloc    = nD + nL - 1 ;
dstride = 2*nD + nL + nU - 2 ;
lsize   = nD + nL - 1 ;
usize   = nD + nU - 1 ;
#if MYDEBUG > 0
fprintf(stdout, "\n dloc = %d, dstride = %d, lsize = %d, usize = %d",
        dloc, dstride, lsize, usize) ;
#endif
/*
   ----------------------
   check for a zero pivot
   ----------------------
*/
if ( CHV_IS_REAL(chv) ) {
   dreal = entries[dloc] ;
   if (  dreal == 0.0 ) {
      return(0) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   dreal = entries[2*dloc] ;
   dimag = entries[2*dloc+1] ;
   if (  dreal == 0.0 && dimag == 0.0 ) {
      return(0) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n (dreal,dimag) = <%12.4e,%12.4e>", 
           dreal, dimag) ;
#endif
}
/*
   -----------------------------------------
   compute the inverse of the diagonal pivot
   real:    fac1 = 1/d(0,0)
   complex: (fac1,fac2) = 1/d(0,0)
   -----------------------------------------
*/
if ( CHV_IS_REAL(chv) ) {
   fac1 = 1./dreal ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   Zrecip(dreal, dimag, &fac1, &fac2) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n 1/(dreal,dimag) = <%12.4e,%12.4e>", 
           fac1, fac2) ;
#endif
}
/*
   ---------------------------
   scale the first column of L
   (fac1,fac2) = 1/d(0,0)
   ---------------------------
*/
if ( CHV_IS_REAL(chv) ) {
   DVscale(lsize, entries, fac1) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   ZVscale(lsize, entries, fac1, fac2) ;
#if MYDEBUG > 2
   { double   real, imag ;
     int      irow ;
      fprintf(stdout, "\n entries in L after scaling") ;
      for ( irow = 1 ; irow < nD + nL ; irow++ ) {
         Chv_entry(chv, irow, 0, &real, &imag) ;
         fprintf(stdout, "\n %% A(%d,%d) = %20.12e + %20.12ei", 
                 irow, 0, real, imag) ;
      }
   }
#endif
}
/*
   ------------------------------------
   loop over the following chevrons
   ljiloc -- offset into lij multiplier
   uijloc -- offset into uij multiplier
   ------------------------------------
*/
ljiloc  = dloc - 1 ;
uijloc  = dloc + 1 ;
for ( kchv = 1 ; kchv < nD ; kchv++ ) {
/*
   --------------------------------------------------
   dloc now points to next diagonal location
   lsize and usize decremented
   lkbeg -- offset into start of column in lower part
   ukbeg -- offset into start of row in upper part
   --------------------------------------------------
*/
   dloc += dstride ;
   lsize-- ;
   usize-- ;
   lkbeg = dloc - lsize ;
   ukbeg = dloc + 1 ;
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n kchv   = %5d, dloc   = %5d"
           "\n ljiloc = %5d, uijloc = %5d"
           "\n lsize  = %5d, usize  = %5d"
           "\n lkbeg  = %5d, ukbeg  = %5d",
           kchv, dloc, ljiloc, uijloc, lsize, usize, lkbeg, ukbeg) ;
#endif
/*
   ------------------------------------
   pull out the multiplier coefficients
   ------------------------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      lreal = entries[ljiloc] ;
      ureal = entries[uijloc] ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      lreal = entries[2*ljiloc] ;
      limag = entries[2*ljiloc+1] ;
      ureal = entries[2*uijloc] ;
      uimag = entries[2*uijloc+1] ;
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n (lreal,limag) = <%12.4e,%12.4e>"
              "\n (ureal,uimag) = <%12.4e,%12.4e>", 
              lreal, limag, ureal, uimag) ;
#endif
   }
/*
   -------------------------
   update the diagonal entry
   -------------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      entries[dloc] -= lreal*ureal ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      entries[2*dloc]   -= lreal*ureal - limag*uimag ;
      entries[2*dloc+1] -= lreal*uimag + limag*ureal ;
   }
/*
   -----------------------
   update the lower column
   -----------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      DVaxpy(lsize, &entries[lkbeg], -ureal, entries) ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      ZVaxpy(lsize, &entries[2*lkbeg], -ureal, -uimag, entries) ;
   }
/*
   --------------------
   update the upper row
   --------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      DVaxpy(usize, &entries[ukbeg], -lreal, &entries[uijloc+1]) ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      ZVaxpy(usize, &entries[2*ukbeg], 
             -lreal, -limag, &entries[2*uijloc+2]) ;
   }
/*
   ----------------------------------
   adjust offsets and diagonal stride
   ----------------------------------
*/
   ljiloc-- ; uijloc++ ; dstride -= 2 ;
}
/*
   ------------------
   scale the row of U
   ------------------
*/
usize = nD + nU - 1 ;
if ( CHV_IS_REAL(chv) ) {
   DVscale(usize, &entries[nD+nL], fac1) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   ZVscale(usize, &entries[2*(nD+nL)], fac1, fac2) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   perform an internal rank-2 update for 
   a hermitian or symmetric chevron


   return code ---
      0 if the pivot was zero
      1 if the pivot was nonzero

   created -- 98jan23, cca
   -------------------------------------
*/
static int
symmetric2x2 (
  Chv   *chv
) {
double   areal, aimag, breal, bimag, creal, cimag, 
         ereal, eimag, freal, fimag, greal, gimag,
         l0imag, l1imag, l0real, l1real,
         u0imag, u1imag, u0real, u1real ;
double   *entries ;
int      dloc, dstride, kchv, nD, nL, nU, rc, 
         u0jloc, u1jloc, ukbeg, usize ;
/*
   -------------------------------------
   get dimensions and pointer to entries
   -------------------------------------
*/
Chv_dimensions(chv, &nD, &nL, &nU) ;
entries = Chv_entries(chv) ;
#if MYDEBUG > 0
fprintf(stdout, "\n nD = %d, nL = %d, nU = %d, entries = %p",
        nD, nL, nU, entries) ;
#endif
/*
   ----------------------------------------
   check for a zero pivot
   D = [    a    b ] for hermitian
       [ conj(b) c ]
   D = [ a b ] for symmetric
       [ b c ]
   compute the inverse of D
   E = inv(D) = [    e    f ] for hermitian
                [ conj(f) g ]
   E = inv(D) = [ e f ] for symmetric
                [ f g ]
   ----------------------------------------
*/
if ( CHV_IS_REAL(chv) ) {
   double   denom ;

   areal = entries[0] ;
   breal = entries[1] ;
   creal = entries[nD+nU] ;
   if ( (denom = areal*creal - breal*breal) == 0.0 ) {
      rc = 0 ;
   } else {
      rc = 1 ;
      denom = 1./denom ;
      ereal =  creal*denom ;
      freal = -breal*denom ;
      greal =  areal*denom ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   areal = entries[0] ;
   aimag = entries[1] ;
   breal = entries[2] ;
   bimag = entries[3] ;
   creal = entries[2*(nD+nU)] ;
   cimag = entries[2*(nD+nU)+1] ;
#if MYDEBUG > 0
   if ( CHV_IS_HERMITIAN(chv) ) {
      fprintf(stdout, 
              "\n hermitian D = [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]"
              "\n               [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]", 
              areal,  aimag, breal, bimag,
              breal, -bimag, creal, cimag) ;
   } else {
      fprintf(stdout, 
              "\n symmetric D = [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]" 
              "\n               [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]", 
              areal, aimag, breal, bimag,
              breal, bimag, creal, cimag) ;
   }
#endif
   if ( CHV_IS_HERMITIAN(chv) ) {
      rc = Zrecip2(areal,  0.0,    breal,  bimag, 
                   breal,  -bimag, creal,  0.0,
                   &ereal, NULL,   &freal, &fimag,
                   NULL,   NULL,   &greal, NULL) ;
      eimag = gimag = 0.0 ;
   } else {
      rc = Zrecip2(areal,  aimag,  breal,  bimag, 
                   breal,  bimag,  creal,  cimag,
                   &ereal, &eimag, &freal, &fimag,
                   NULL,   NULL,   &greal, &gimag) ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_symmetric2x2"
           "\n chevron must be real or complex") ;
   exit(-1) ;
}
if ( rc == 0 ) {
/*
   -------------------------
   pivot is singular, return
   -------------------------
*/
   return(0) ;
}
#if MYDEBUG > 0
if ( CHV_IS_HERMITIAN(chv) ) {
   fprintf(stdout, 
           "\n hermitian DINV = [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]"
           "\n                  [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]", 
           ereal,  eimag, freal, fimag,
           freal, -fimag, greal, gimag) ;
} else {
   fprintf(stdout, 
           "\n symmetric DINV = [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]" 
           "\n                  [ <%12.4e,%12.4e> <%12.4e,%12.4e> ]", 
           ereal, eimag, freal, fimag,
           freal, fimag, greal, gimag) ;
}
#endif
/*
   -----------------------------
   scale the first two rows of U
   -----------------------------
*/
u0jloc = 2 ;
u1jloc = nD + nU + 1 ;
usize  = nD + nU - 2 ;
if ( CHV_IS_REAL(chv) ) {
   DVscale2(usize, &entries[u0jloc], &entries[u1jloc],
           ereal, freal, freal, greal) ;
} else if ( CHV_IS_HERMITIAN(chv) ) {
   ZVscale2(usize, &entries[2*u0jloc], &entries[2*u1jloc],
           ereal, 0.0, freal, fimag, freal, -fimag, greal, 0.0) ;
} else {
   ZVscale2(usize, &entries[2*u0jloc], &entries[2*u1jloc],
           ereal, eimag, freal, fimag, freal, fimag, greal, gimag) ;
}
#if MYDEBUG > 2
{ double   real, imag ;
  int      irow, jcol ;
  fprintf(stdout, "\n entries in U after scaling") ;
  for ( irow = 0 ; irow <= 1 ; irow++ ) {
     for ( jcol = 2 ; jcol < nD + nU ; jcol++ ) {
        Chv_entry(chv, 0, jcol, &real, &imag) ;
        fprintf(stdout, "\n %% A(%d,%d) = %20.12e + %20.12ei", 
                0, jcol, real, imag) ;
     }
  }
}
#endif
/*
   ------------------------------------
   loop over the following chevrons
   u0jloc -- offset into u0j multiplier
   u1jloc -- offset into u1j multiplier
   ------------------------------------
*/
usize   = nD + nU - 2 ;
dloc    = nD + nU ;
dstride = nD + nU - 1 ;
for ( kchv = 2 ; kchv < nD ; kchv++ ) {
/*
   --------------------------------------------------
   dloc now points to next diagonal location
   ukbeg -- offset into start of row in upper part
   --------------------------------------------------
*/
   dloc += dstride ;
   ukbeg = dloc + 1 ;
#if MYDEBUG > 0
   fprintf(stdout, "\n kchv = %5d, dloc = %5d"
           "\n u0jloc = %5d, u1jloc = %d, usize = %5d, ukbeg = %5d",
           kchv, dloc, u0jloc, u1jloc, usize, ukbeg) ;
#endif
/*
   ------------------------------------
   pull out the multiplier coefficients
   ------------------------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      u0real = entries[u0jloc] ;
      u1real = entries[u1jloc] ;
      l0real = u0real*areal + u1real*breal ;
      l1real = u0real*breal + u1real*creal ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      u0real = entries[2*u0jloc] ;
      u0imag = entries[2*u0jloc+1] ;
      u1real = entries[2*u1jloc] ;
      u1imag = entries[2*u1jloc+1] ;
      if ( CHV_IS_HERMITIAN(chv) ) {
         l0real = u0real*areal + u1real*breal - u1imag*bimag ;
         l0imag = -u0imag*areal - u1real*bimag - u1imag*breal ;
         l1real = u0real*breal + u0imag*bimag + u1real*creal ;
         l1imag = u0real*bimag - u0imag*breal - u1imag*creal ;
      } else {
         l0real = u0real*areal - u0imag*aimag
                + u1real*breal - u1imag*bimag ;
         l0imag = u0real*aimag + u0imag*areal
                + u1real*bimag + u1imag*breal ;
         l1real = u0real*breal - u0imag*bimag
                + u1real*creal - u1imag*cimag ;
         l1imag = u0real*bimag + u0imag*breal
                + u1real*cimag + u1imag*creal ;
      }
#if MYDEBUG > 0
      fprintf(stdout, 
              "\n (l0real,l0imag) = <%12.4e,%12.4e>"
              "\n (l1real,l1imag) = <%12.4e,%12.4e>", 
              l0real, l0imag, l1real, l1imag) ;
#endif
   }
/*
   -------------------------------------------------
   update the upper row including the diagonal entry
   -------------------------------------------------
*/
   if ( CHV_IS_REAL(chv) ) {
      DVaxpy2(usize, &entries[dloc], -l0real, &entries[u0jloc],
                                     -l1real, &entries[u1jloc]) ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      ZVaxpy2(usize, &entries[2*dloc], 
              -l0real, -l0imag, &entries[2*u0jloc],
              -l1real, -l1imag, &entries[2*u1jloc]) ;
   }
/*
   ----------------------------------
   adjust offsets and diagonal stride
   ----------------------------------
*/
   u0jloc++ ; u1jloc++ ; dstride-- ; usize-- ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- looking at just a single chevron inside the Chv object,
              find the absolute value of the diagonal element, and
              the maximum absolute values of the offdiagonal elements 
              in the chevron's row and column.

   created -- 98aug26, cca
   ------------------------------------------------------------------
*/
void
Chv_maxabsInChevron (
   Chv      *chv,
   int      ichv,
   double   *pdiagmaxabs,
   double   *prowmaxabs,
   double   *pcolmaxabs
) {
double   *pdiag ;
int      length, loc, nD, nL, nU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || ichv < 0 || ichv >= chv->nD
  || pdiagmaxabs == NULL || prowmaxabs == NULL || pcolmaxabs == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_maxabsInChevron()"
           "\n bad input\n") ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
pdiag  = Chv_diagLocation(chv, ichv) ;
length = nD - ichv - 1 + nU ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) ) {
      *pdiagmaxabs = fabs(*pdiag) ;
      *prowmaxabs  = DVmaxabs(length, pdiag + 1, &loc) ;
      *pcolmaxabs  = *prowmaxabs ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      *pdiagmaxabs = fabs(*pdiag) ;
      *prowmaxabs  = DVmaxabs(length, pdiag + 1, &loc) ;
      *pcolmaxabs  = DVmaxabs(length, pdiag - length, &loc) ;
   } else {
      fprintf(stderr, "\n fatal error in Chv_maxabsInChevron()"
              "\n chv is real, chv->symflag = %d"
              "\n must be symmetric or nonsymmetric\n", chv->symflag) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      *pdiagmaxabs = Zabs(*pdiag, *(pdiag+1)) ;
      *prowmaxabs  = ZVmaxabs(length, pdiag + 2) ;
      *pcolmaxabs  = *prowmaxabs ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      *pdiagmaxabs = Zabs(*pdiag, *(pdiag+1)) ;
      *prowmaxabs  = ZVmaxabs(length, pdiag + 2) ;
      *pcolmaxabs  = ZVmaxabs(length, pdiag - 2*length) ;
   } else {
      fprintf(stderr, "\n fatal error in Chv_maxabsInChevron()"
              "\n chv is complex, chv->symflag = %d"
              "\n must be symmetric or nonsymmetric\n", chv->symflag) ;
      exit(-1) ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_maxabsInChevron()"
           "\n chv->type = %d, must be real or complex\n", chv->type) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- zero the offdiagonal entries of chevron ichv

   created -- 98aug26, cca
   -------------------------------------------------------
*/
void
Chv_zeroOffdiagonalOfChevron (
   Chv   *chv,
   int   ichv
) {
double   *pdiag ;
int      length, nD, nL, nU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || ichv < 0 || ichv >= chv->nD ) {
   fprintf(stderr, "\n fatal error in Chv_zeroOffdiagonalOfChevron()"
           "\n bad input\n") ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
pdiag  = Chv_diagLocation(chv, ichv) ;
length = nD - ichv - 1 + nU ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) ) {
      DVzero(length, pdiag+1) ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      DVzero(length, pdiag + 1) ;
      DVzero(length, pdiag - length) ;
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_zeroOffdiagonalOfChevron()"
              "\n chv is real, chv->symflag = %d"
              "\n must be symmetric or nonsymmetric\n", chv->symflag) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      ZVzero(length, pdiag+2) ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      ZVzero(length, pdiag+2) ;
      ZVzero(length, pdiag-2*length) ;
   } else {
      fprintf(stderr, 
              "\n fatal error in Chv_zeroOffdiagonalOfChevron()"
              "\n chv is complex, chv->symflag = %d"
              "\n must be symmetric or nonsymmetric\n", chv->symflag) ;
      exit(-1) ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_zeroOffdiagonalOfChevron()"
           "\n chv->type = %d, must be real or complex\n", chv->type) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
