/*  findPivot.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
static int findPivotSH ( Chv *chv, DV *workDV, double tau,
   int ndelay, int *pirow, int *pjcol, int *pntest ) ;
static int findPivotN ( Chv *chv, DV *workDV, double tau,
   int ndelay, int *pirow, int *pjcol, int *pntest ) ;
static int nonsym1x1 ( Chv *chv, int irow, int jcol, double tau,
                       double rowmaxes[], double colmaxes[] ) ;
static int sym1x1 ( Chv *chv, int irow, 
                    double tau, double rowmaxes[] );
static int sym2x2 ( Chv *chv, int irow, int jcol, 
                    double tau, double rowmaxes[] ) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- find and test a pivot

   workDV -- object that contains work vectors
   tau    -- upper bound on magnitude of factor entries
   ndelay -- number of delayed rows and columns on input
   pirow  -- pointer to be filled with pivot row
   pjcol  -- pointer to be filled with pivot column
   pntest -- pointer to be incremented with the number of pivot tests

   return value -- size of pivot
     0 --> pivot not found
     1 --> 1x1 pivot in row *pirow and column *pjcol
     2 --> 2x2 pivot in rows and columns *pirow and *pjcol,
           symmetric front only
   
   created -- 98jan24, cca
   ------------------------------------------------------------------
*/
int
Chv_findPivot (
   Chv     *chv,
   DV       *workDV,
   double   tau,
   int      ndelay,
   int      *pirow,
   int      *pjcol,
   int      *pntest
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL || workDV == NULL || tau < 1.0 || ndelay < 0
   || pirow == NULL || pjcol == NULL || pntest == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_findPivot(%p,%p,%f,%d,%p,%p,%p)"
           "\n bad input\n", 
           chv, workDV, tau, ndelay, pirow, pjcol, pntest) ;
   exit(-1) ;
}
if ( !(CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, 
           "\n fatal error in Chv_findPivot(%p,%p,%f,%d,%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           chv, workDV, tau, ndelay, pirow, pjcol, pntest, chv->type) ;
   exit(-1) ;
}
if ( !(CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) 
        || CHV_IS_NONSYMMETRIC(chv)) ) {
   fprintf(stderr, 
        "\n fatal error in Chv_findPivot(%p,%p,%f,%d,%p,%p,%p)"
        "\n bad symflag %d"
        "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or CHV_NONSYMMETRIC\n",
        chv, workDV, tau, ndelay, pirow, pjcol, pntest, chv->symflag) ;
   exit(-1) ;
}
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   rc = findPivotSH(chv, workDV, tau, ndelay, pirow, pjcol, pntest) ;
} else if ( CHV_IS_NONSYMMETRIC(chv) ) {
   rc = findPivotN(chv, workDV, tau, ndelay, pirow, pjcol, pntest) ;
} else {
   fprintf(stderr, 
           "\n fatal error in Chv_findPivot(%p,%p,%f,%d,%p,%p,%p)"
           "\n bad symflag %d\n", chv, workDV, tau, ndelay, pirow, 
           pjcol, pntest, chv->symflag) ;
   exit(-1) ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- find and test a pivot for a symmetric or hermitian matrix

   workDV -- object that contains work vectors
   tau    -- upper bound on magnitude of factor entries
   ndelay -- number of delayed rows and columns on input
   pirow  -- pointer to be filled with pivot row
   pjcol  -- pointer to be filled with pivot column
   pntest -- pointer to be incremented with the number of pivot tests

   return value -- size of pivot
     0 --> pivot not found
     1 --> 1x1 pivot in row *pirow and column *pjcol
     2 --> 2x2 pivot in rows and columns *pirow and *pjcol,
           symmetric front only
   
   created -- 98apr18, cca
   --------------------------------------------------------------------
*/
static int
findPivotSH (
   Chv     *chv,
   DV       *workDV,
   double   tau,
   int      ndelay,
   int      *pirow,
   int      *pjcol,
   int      *pntest
) {
double   maxval ;
double   *rowmaxes ;
int      ii, irow, jrow, krow, ncand, nD, 
         ndouble, ntest, pivotsize, tag, untag ;
int      *rowids, *rowmark ;

untag  = 0 ;
tag    = 1 ;
nD     = chv->nD ;
#if MYDEBUG > 0
fprintf(stdout, 
"\n %% findPivotSH, id = %d, nD = %d, nL = %d, nU = %d, ndelay = %d",
chv->id, chv->nD, chv->nL, chv->nU, ndelay) ;
fflush(stdout) ;
#endif
*pirow = *pjcol = -1 ;
ntest  = *pntest ;
/*
   ------------------------------------
   symmetric front, set up work vectors
   ------------------------------------
*/
if ( sizeof(int) == sizeof(double) ) {
   ndouble = 3*nD ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   ndouble = 2*nD ;
}
DV_setSize(workDV, ndouble) ;
rowmaxes = DV_entries(workDV) ;
DVfill(nD, rowmaxes, 0.0) ;
rowmark  = (int *) (rowmaxes + nD) ;
rowids   = rowmark + nD ;
if ( ndelay > 0 ) {
   IVfill(ndelay, rowmark, untag) ;
   IVfill(nD - ndelay, rowmark + ndelay, tag) ;
} else {
   IVfill(nD, rowmark, tag) ;
}
ncand = 0 ;
do {
   pivotsize = 0 ;
/*
   --------------------
   find candidate pivot
   --------------------
*/
   Chv_fastBunchParlettPivot(chv, rowmark, tag, &irow, &jrow) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n %% FBP: irow = %d, jrow = %d",
           irow, jrow) ;
   if ( irow != -1 ) {
      double   imag, real ;
      Chv_entry(chv, irow, irow, &real, &imag) ;
      fprintf(stdout, "\n%%  entry(%d,%d) = %20.12e + %20.12e*i",
              irow, irow, real, imag) ;
      if ( jrow != irow ) {
         Chv_entry(chv, irow, jrow, &real, &imag) ;
         fprintf(stdout, "\n%%  entry(%d,%d) = %20.12e + %20.12e*i",
                 irow, jrow, real, imag) ;
         Chv_entry(chv, jrow, jrow, &real, &imag) ;
         fprintf(stdout, "\n%%  entry(%d,%d) = %20.12e + %20.12e*i",
                 jrow, jrow, real, imag) ;
      }
   }
   fflush(stdout) ;
#endif
   if ( irow == -1 ) {
/*
      ----------------------------------------------
      unable to find pivot, break out of search loop
      ----------------------------------------------
*/
      pivotsize = 0 ; break ;
   } else {
/*
      -------------------------------
      (irow,jrow) is a possible pivot
      mark as visited and get row max
      -------------------------------
*/
      Chv_maxabsInRow(chv, irow, &maxval) ;
      rowmaxes[irow] = maxval ;
      rowmark[irow]  = untag ;
      if ( irow != jrow ) {
         Chv_maxabsInRow(chv, jrow, &maxval) ;
         rowmaxes[jrow] = maxval ;
         rowmark[jrow] = untag ;
      }
      if ( irow == jrow ) {
/*
         ------------------
         test the 1x1 pivot
         ------------------
*/
         pivotsize = sym1x1(chv, irow, tau, rowmaxes) ;
         ntest++ ;
#if MYDEBUG > 0
         fprintf(stdout, 
                 "\n\n %% pivotsize from sym1x1 = %d", pivotsize) ;
#endif
         if ( pivotsize == 1 ) {
            *pirow = irow ; *pjcol = jrow ;
         } else {
            for ( ii = 0 ; ii < ncand ; ii++ ) {
/*
               ----------------------------------
               test the 2x2 pivot (irow, krow)
               where krow is a previous candidate
               ----------------------------------
*/
               krow = rowids[ii] ;
               pivotsize = sym2x2(chv, irow, krow, tau, rowmaxes) ; 
               ntest++ ;
               if ( pivotsize == 2 ) {
                  *pirow = irow ; *pjcol = krow ; break ;
               }
            }
         }
      } else {
/*
         -------------------------------
         test the 2x2 pivot (irow, jrow)
         -------------------------------
*/
         pivotsize = sym2x2(chv, irow, jrow, tau, rowmaxes) ; 
         ntest++ ;
         if ( pivotsize == 2 ) {
            *pirow = irow ; *pjcol = jrow ; 
         } else {
            for ( ii = 0 ; ii < ncand ; ii++ ) {
               krow = rowids[ii] ;
/*
               ----------------------------------
               test the 2x2 pivot (irow, krow)
               where krow is a previous candidate
               ----------------------------------
*/
               pivotsize = sym2x2(chv, irow, krow, tau, rowmaxes) ;
               ntest++ ;
               if ( pivotsize == 2 ) {
                  *pirow = irow ; *pjcol = krow ; break ;
               }
/*
               ----------------------------------
               test the 2x2 pivot (jrow, krow)
               where krow is a previous candidate
               ----------------------------------
*/
               pivotsize = sym2x2(chv, jrow, krow, tau, rowmaxes) ;
               ntest++ ;
               if ( pivotsize == 2 ) {
                  *pirow = jrow ; *pjcol = krow ; break ;
               }
            }
         }
      }
      if ( pivotsize == 0 ) {
/*
         ------------------------
         add new candidate row(s)
         ------------------------
*/
         rowids[ncand++] = irow ;
         if ( irow != jrow ) {
            rowids[ncand++] = jrow ;
         }
      }
   }
} while ( pivotsize == 0 ) ;
*pntest = ntest ;

return(pivotsize) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- find and test a pivot for a nonsymmetric matrix

   workDV -- object that contains work vectors
   tau    -- upper bound on magnitude of factor entries
   ndelay -- number of delayed rows and columns on input
   pirow  -- pointer to be filled with pivot row
   pjcol  -- pointer to be filled with pivot column
   pntest -- pointer to be incremented with the number of pivot tests

   return value -- size of pivot
     0 --> pivot not found
     1 --> 1x1 pivot in row *pirow and column *pjcol
     2 --> 2x2 pivot in rows and columns *pirow and *pjcol,
           symmetric front only
   
   created -- 98jan24, cca
   ------------------------------------------------------------------
*/
static int
findPivotN (
   Chv     *chv,
   DV       *workDV,
   double   tau,
   int      ndelay,
   int      *pirow,
   int      *pjcol,
   int      *pntest
) {
double   maxval ;
double   *colmaxes, *rowmaxes ;
int      icol, ii, irow, jcol, jrow, ncand, nD, 
         ndouble, ntest, pivotsize, tag, untag ;
int      *colids, *colmark, *rowids, *rowmark ;

untag  = 0 ;
tag    = 1 ;
nD     = chv->nD ;
#if MYDEBUG > 0
fprintf(stdout, 
"\n %% Chv_findPivot, id = %d, nD = %d, nL = %d, nU = %d, ndelay = %d",
chv->id, chv->nD, chv->nL, chv->nU, ndelay) ;
fflush(stdout) ;
#endif
*pirow = *pjcol = -1 ;
ntest  = *pntest ;
/*
   -------------------
   set up work vectors
   -------------------
*/
if ( sizeof(int) == sizeof(double) ) {
   ndouble = 6*nD ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   ndouble = 4*nD ;
}
DV_setSize(workDV, ndouble) ;
rowmaxes = DV_entries(workDV) ;
colmaxes = rowmaxes + nD ;
DVfill(nD, rowmaxes, 0.0) ;
DVfill(nD, colmaxes, 0.0) ;
rowmark  = (int *) (colmaxes + nD) ;
colmark  = rowmark + nD ;
rowids   = colmark + nD ;
colids   = rowids + nD ;
if ( ndelay > 0 ) {
   IVfill(ndelay, rowmark, untag) ;
   IVfill(nD - ndelay, rowmark + ndelay, tag) ;
   IVfill(ndelay, colmark, untag) ;
   IVfill(nD - ndelay, colmark + ndelay, tag) ;
} else {
   IVfill(nD, rowmark, tag) ;
   IVfill(nD, colmark, tag) ;
}
ncand = 0 ;
do {
   pivotsize = 0 ;
/*
   --------------------
   find candidate pivot
   --------------------
*/
   Chv_quasimax(chv, rowmark, colmark, tag, &irow, &jcol) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n %% quasimax: irow = %d, jcol = %d",
           irow, jcol) ;
   if ( irow != -1 ) {
      double   imag, real ;
      Chv_entry(chv, irow, jcol, &real, &imag) ;
      fprintf(stdout, "\n%%  entry(%d,%d) = %20.12e + %20.12e*i",
              irow, jcol, real, imag) ;
   }
   fflush(stdout) ;
#endif
   if ( irow == -1 ) {
/*
      ----------------------------------------------
      unable to find pivot, break out of search loop
      ----------------------------------------------
*/
      break ;
   } else {
/*
      ------------------------------------------------------------
      find the row max for row irow and column max for column jcol
      ------------------------------------------------------------
*/
      Chv_maxabsInRow(chv, irow, &maxval) ;
      rowmaxes[irow] = maxval ;
      Chv_maxabsInColumn(chv, jcol, &maxval) ;
      colmaxes[jcol] = maxval ;
      rowmark[irow]  = untag ;
      colmark[jcol]  = untag ;
/*
      -------------------------------------
      test the (irow,jcol) entry as a pivot
      -------------------------------------
*/
      pivotsize = nonsym1x1(chv, irow, jcol, tau, rowmaxes, colmaxes) ;
      ntest++ ;
      if ( pivotsize == 1 ) {
         *pirow = irow ; *pjcol = jcol ; 
      } else {
/*
         ---------------------------------------
         test the other matrix entries as pivots
         ---------------------------------------
*/
         for ( ii = 0 ; ii < ncand ; ii++ ) {
            jrow = rowids[ii] ;
            icol = colids[ii] ;
/*
            --------------------------
            test the (irow,icol) entry
            --------------------------
*/
            pivotsize = nonsym1x1(chv, irow, icol, tau, 
                                  rowmaxes, colmaxes) ;
            ntest++ ;
            if ( pivotsize == 1 ) {
               *pirow = irow ; *pjcol = icol ; break ;
            }
/*
            --------------------------
            test the (jrow,jcol) entry
            --------------------------
*/
            pivotsize = nonsym1x1(chv, jrow, jcol, tau, 
                                  rowmaxes, colmaxes) ;
            ntest++ ;
            if ( pivotsize == 1 ) {
               *pirow = jrow ; *pjcol = jcol ; break ;
            }
         }
/*
         ----------------------------------------------
         no pivots found, add irow to candidate row ids
         and add jcol to candidate column ids
         ----------------------------------------------
*/
         rowids[ncand] = irow ;
         colids[ncand] = jcol ;
         ncand++ ;
      }
   }
} while ( pivotsize == 0 ) ;
*pntest = ntest ;

return(pivotsize) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   return 1 if the nonsymmetric 1x1 pivot passes
   return 0 otherwise

   created -- 98jan24, cca
   ---------------------------------------------
*/
static int
nonsym1x1 (
   Chv   *chv,
   int    irow,
   int    jcol,
   double   tau,
   double   rowmaxes[],
   double   colmaxes[]
) {
double   cutoff, magn ;
int      rc ;

if ( CHV_IS_REAL(chv) ) {
   double   value ;
   Chv_realEntry(chv, irow, jcol, &value) ;
   magn = fabs(value) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   double   imag, real ;
   Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
   magn = Zabs(real, imag) ;
}
cutoff = tau * magn ;
#if MYDEBUG > 0
fprintf(stdout, "\n %% magn = %12.4e, cutoff = %12.4e", magn, cutoff) ;
fprintf(stdout, "\n %% rowmaxes[%d] = %12.4e, colmaxes[%d] = %12.4e",
        irow, rowmaxes[irow], jcol, colmaxes[jcol]) ;
#endif
if ( rowmaxes[irow] <= cutoff && colmaxes[jcol] <= cutoff ) {
   rc = 1 ;
} else {
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return 1 if the symmetric 1x1 pivot passes
   return 0 otherwise

   created -- 98jan24, cca
   ------------------------------------------
*/
static int
sym1x1 (
   Chv     *chv,
   int      irow,
   double   tau,
   double   rowmaxes[]
) {
double   cutoff ;
int      rc ;

if ( CHV_IS_REAL(chv) ) {
   double   value ;
   Chv_realEntry(chv, irow, irow, &value) ;
   cutoff = tau * fabs(value) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   double   imag, real ;
   Chv_complexEntry(chv, irow, irow, &real, &imag) ;
   cutoff = tau * Zabs(real, imag) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n %% cutoff = %12.4e, rowmaxes[%d] = %12.4e", 
        cutoff, irow, rowmaxes[irow]) ;
#endif
if ( rowmaxes[irow] <= cutoff ) {
   rc = 1 ;
} else {
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return 2 if the symmetric 2x2 pivot passes
   return 0 otherwise

   created -- 98jan24, cca
   ------------------------------------------
*/
static int
sym2x2 (
   Chv     *chv,
   int      irow,
   int      jcol,
   double   tau,
   double   rowmaxes[]
) {
double   amag, bmag, cmag, denom, val1, val2 ;
int      rc ;

if ( CHV_IS_REAL(chv) ) {
   double  a, b, c ;

   Chv_realEntry(chv, irow, irow, &a) ;
   Chv_realEntry(chv, irow, jcol, &b) ;
   Chv_realEntry(chv, jcol, jcol, &c) ;
   amag  = fabs(a) ;
   bmag  = fabs(b) ;
   cmag  = fabs(c) ;
   denom = fabs(a*c - b*b) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   double   aimag, areal, bimag, breal, cimag, creal, imag, real ;

   Chv_complexEntry(chv, irow, irow, &areal, &aimag) ;
   Chv_complexEntry(chv, irow, jcol, &breal, &bimag) ;
   Chv_complexEntry(chv, jcol, jcol, &creal, &cimag) ;
   if ( CHV_IS_HERMITIAN(chv) ) {
      amag  = fabs(areal) ;
      bmag  = Zabs(breal, bimag) ;
      cmag  = fabs(creal) ;
      denom = fabs(areal*creal - breal*breal - bimag*bimag) ;
   } else if ( CHV_IS_SYMMETRIC(chv) ) {
      amag  = Zabs(areal, aimag) ;
      bmag  = Zabs(breal, bimag) ;
      cmag  = Zabs(creal, cimag) ;
      real  = areal*creal - aimag*cimag - breal*breal + bimag*bimag ;
      imag  = areal*cimag + aimag*creal - 2*breal*bimag ;
      denom = Zabs(real, imag) ;
   }
}
#if MYDEBUG > 0
      fprintf(stdout,
              "\n amag = %20.12e ; "
              "\n bmag = %20.12e ; "
              "\n cmag = %20.12e ; ", amag, bmag, cmag) ;
#endif
if ( denom == 0.0 ) {
   return(0) ;
}
val1 = (cmag*rowmaxes[irow] + bmag*rowmaxes[jcol])/denom ;
val2 = (bmag*rowmaxes[irow] + amag*rowmaxes[jcol])/denom ;
#if MYDEBUG > 0
fprintf(stdout, "\n %% sym2x2"
        "\n rowmax1 = %20.12e"
        "\n rowmax2 = %20.12e"
        "\n val1 = %20.12e"
        "\n val2 = %20.12e"
        "\n denom = %20.12e", 
        rowmaxes[irow], rowmaxes[jcol], val1, val2, denom) ;
#endif
if (  val1 <= tau && val2 <= tau ) {
   rc = 2 ;
} else {
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
