/*  QRreduce.c  */

#include "../A2.h"

/*--------------------------------------------------------------------*/
static int copyIntoVec1 ( A2 *mtxA, double H0[], 
   int msglvl, FILE *msgFile ) ;
static double getHouseholderVector1 ( int type, int n, double H0[],
   double *pbeta0, int msglvl, FILE *msgFile ) ;
static double computeW1 ( A2 *mtxA, double H0[], double W0[],
   int msglvl, FILE *msgFile ) ;
static double updateA1 ( A2 *mtxA, double H0[], double beta0,
   double W0[], int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- compute A = QR, where Q is a product of householder
              vectors, (I - beta_j v_j v_j^T). on return, v_j is 
              found in the lower triangle of A, v_j(j) = 1.0.

   return value -- # of floating point operations

   created -- 98may25, cca
   --------------------------------------------------------------
*/
double
A2_QRreduce (
   A2       *mtxA,
   DV       *workDV,
   int      msglvl,
   FILE     *msgFile
) {
A2       tempA ;
double   nops ;
double   beta0 ;
double   *colA, *H0, *W0 ;
int      inc1, inc2, jcol, lastrow, length, ncolA, nrowA, nstep ;
/*
   ---------------
   check the input
   ---------------
*/
if (   mtxA == NULL || workDV == NULL
    || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in A2_QRreduce()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtxA) || A2_IS_COMPLEX(mtxA)) ) {
   fprintf(stderr, "\n fatal error in A2_QRreduce()"
           "\n matrix must be real or complex\n") ;
   exit(-1) ;
}
nrowA = A2_nrow(mtxA) ; 
ncolA = A2_ncol(mtxA) ;
inc1  = A2_inc1(mtxA) ;
inc2  = A2_inc2(mtxA) ;
if ( A2_IS_REAL(mtxA) ) {
   DV_setSize(workDV, nrowA + ncolA) ;
   H0 = DV_entries(workDV) ;
   W0 = H0 + nrowA ;
} else if ( A2_IS_COMPLEX(mtxA) ) {
   DV_setSize(workDV, 2*(nrowA + ncolA)) ;
   H0 = DV_entries(workDV) ;
   W0 = H0 + 2*nrowA ;
}
/*
   -------------------------------------------------
   determine the number of steps = min(ncolA, nrowA)
   -------------------------------------------------
*/
nstep = (ncolA <= nrowA) ? ncolA : nrowA ;
/*
   -------------------
   loop over the steps
   -------------------
*/
nops = 0.0 ; 
for ( jcol = 0 ; jcol < nstep ; jcol++ ) {
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n %% jcol = %d", jcol) ;
   }
/*
   ----------------------------------
   copy the column of A into a vector
   and find the last nonzero element
   ----------------------------------
*/
   A2_subA2(&tempA, mtxA, jcol, nrowA-1, jcol, ncolA-1) ;
   length = 1 + copyIntoVec1(&tempA, H0, msglvl, msgFile) ;
   lastrow = jcol + length - 1 ;
   if ( msglvl > 5 ) {
      fprintf(msgFile, 
            "\n %% return from copyIntoVec1, length = %d, lastrow = %d",
            length, lastrow) ;
   }
/*
   ------------------------------
   compute the Householder vector
   and place into the column of A
   ------------------------------
*/
   colA = A2_column(mtxA, jcol) ;
   if ( A2_IS_REAL(mtxA) ) {
      nops += getHouseholderVector1(SPOOLES_REAL, length, H0, 
                                    &beta0, msglvl, msgFile) ;
      A2_subA2(&tempA, mtxA, jcol, lastrow, jcol, jcol) ;
      A2_setColumn(&tempA, H0, 0) ;
      H0[0] = 1.0 ;
   } else if ( A2_IS_COMPLEX(mtxA) ) {
      nops += getHouseholderVector1(SPOOLES_COMPLEX, length, H0, 
                                    &beta0, msglvl, msgFile) ;
      A2_subA2(&tempA, mtxA, jcol, lastrow, jcol, jcol) ;
      A2_setColumn(&tempA, H0, 0) ;
      H0[0] = 1.0 ; H0[1] = 0.0 ;
   }
   if ( msglvl > 5 && jcol == 0 ) {
      fprintf(msgFile, "\n %% beta0 = %12.4e;", beta0) ;
   }
   if ( beta0 != 0.0 && jcol + 1 < ncolA ) {
      A2_subA2(&tempA, mtxA, jcol, lastrow, jcol+1, ncolA-1) ;
/*
      ------------------------------------------------
      compute w = v^T * A(jcol:lastrow,jcol+1:nrowA-1)
      ------------------------------------------------
*/
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n %% compute w") ;
      }
      nops += computeW1(&tempA, H0, W0, msglvl, msgFile) ;
/*
      -------------------------------------------------
      update A(jcol:lastrow,jcol+1:nrowA-1) -= beta*v*w
      -------------------------------------------------
*/
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n %% update A") ;
      }
      nops += updateA1(&tempA, H0, beta0, W0, msglvl, msgFile) ;
   }
}
return(nops) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   copy the first column of mtxA into the vector H0[]

   created -- 98may30, cca
   --------------------------------------------------
*/
static int
copyIntoVec1 (
   A2       *mtxA,
   double   H0[],
   int      msglvl,
   FILE     *msgFile
) {
double   ival, rval ;
double   *colA ;
int      ii, inc1, irow, jj, lastrow, ncolA, nrowA ;
/*
   ----------------------------------
   copy the column of A into a vector
   and find the last nonzero element
   ----------------------------------
*/
nrowA   = mtxA->n1 ;
ncolA   = mtxA->n2 ;
inc1    = mtxA->inc1 ;
lastrow = -1 ;
colA    = A2_column(mtxA, 0) ;
if ( A2_IS_REAL(mtxA) ) {
   for ( irow = ii = jj = 0 ;
         irow < nrowA ;
         irow++, ii += inc1, jj++ ) {
      rval = colA[ii] ; 
      if ( rval != 0.0 ) {
         H0[jj] = rval ; 
         lastrow = irow ;
      }
   }
} else if ( A2_IS_COMPLEX(mtxA) ) {
   for ( irow = ii = jj = 0 ;
         irow < nrowA ;
         irow++, ii += 2*inc1, jj += 2 ) {
      rval = colA[ii] ; ival = colA[ii+1] ;
      if ( rval != 0.0 || ival != 0.0 ) {
         H0[jj] = rval ; H0[jj+1] = ival ;
         lastrow = irow ;
      }
   }
}
return(lastrow) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   compute a householder transformation
    (I - beta*v*v^H)x = alpha*e_1
   where v[0] = 1.0
   on input, 
      H0 contains x,
   on output,
      beta0 contains beta
      H0[0] = alpha 
      H0[1:n-1] = v[1:n] ;

   created -- 98may30, cca
   ----------------------------------------------------------------
*/
static double
getHouseholderVector1 (
   int      type,
   int      n,
   double   H0[],
   double   *pbeta0,
   int      msglvl,
   FILE     *msgFile
) {
double   beta0, ifac, ival, nops, normx, rfac, rval, sigma, 
         sum, v0imag, v0real, y0imag, y0real ;
int      ii, jj ;
/*
   --------------------------------------------
   compute ||H0(1:n-1)||_2^2 and the
   row that contains the last nonzero entry
   --------------------------------------------
*/
sigma   = 0.0 ; 
beta0   = 0.0 ;
nops    = 0.0 ;
if ( type == SPOOLES_REAL ) {
   for ( ii = 1 ; ii < n ; ii++ ) {
      rval = H0[ii] ;
      sigma   += rval*rval ;
   }
   nops += 2*(n-1) ;
} else if ( type == SPOOLES_COMPLEX ) {
   for ( ii = 1, jj = 2 ; ii < n ; ii++, jj += 2 ) {
      rval = H0[jj] ; ival = H0[jj+1] ;
      sigma   += rval*rval + ival*ival ;
   }
   nops += 4*(n-1) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% sigma = %12.4e", sigma) ;
}
if ( sigma != 0.0 ) {
/*
   --------------------------------------------
   there are nonzero entries below the diagonal
   --------------------------------------------
*/
   if ( type == SPOOLES_REAL ) {
      rval = H0[0] ;
      if ( rval == 0.0 ) {
         normx  = sqrt(sigma) ;
         v0real =   normx ;
         y0real = - normx ;
         nops++ ;
      } else {
         normx  = sqrt(sigma + rval*rval) ;
         rfac   = normx/fabs(rval) ;
         v0real = rval*(1 + rfac) ;
         y0real = -rfac*rval ;
         nops += 7 ;
      }
   } else if ( type == SPOOLES_COMPLEX ) {
      rval = H0[0] ; ival = H0[1] ;
      if ( rval == 0.0 && ival == 0.0 ) {
         normx  = sqrt(sigma) ;
         v0real =   normx ; v0imag = 0.0 ;
         y0real = - normx ; y0imag = 0.0 ;
         nops += 2 ;
      } else {
         normx  = sqrt(sigma + rval*rval + ival*ival) ;
         rfac   = normx/Zabs(rval, ival) ;
         v0real = rval + rfac*rval ; v0imag = ival + rfac*ival ;
         y0real = -rfac*rval ;       y0imag = -rfac*ival ;
         nops += 16 ;
      }
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n %% normx = %12.4e", normx) ;
   }
/*
   -------------------------------------
   scale u so u1 = 1.0 and compute beta0
   -------------------------------------
*/
   if ( type == SPOOLES_REAL ) {
      rfac = 1./v0real ;
      for ( ii = 1 ; ii < n ; ii++ ) {
         H0[ii] *= rfac ;
      }
      sum = 1.0 ;
      for ( ii = 1 ; ii < n ; ii++ ) {
         rval = H0[ii] ;
         sum += rval*rval ;
      }
      nops += 3*(n-1) ;
      beta0 = 2./sum ; 
/*
      rfac = 1./v0real ;
      sum = 1.0 ;
      for ( ii = 1 ; ii < n ; ii++ ) {
         rval = H0[ii] = rfac*H0[ii] ;
         sum += rval*rval ;
      }
      nops += 3*(n-1) ;
      beta0 = 2./sum ; 
*/
   } else if ( type == SPOOLES_COMPLEX ) {
      Zrecip(v0real, v0imag, &rfac, &ifac) ;
      sum = 1.0 ;
      for ( ii = 1, jj = 2 ; ii < n ; ii++, jj += 2 ) {
         rval = H0[jj] ; ival = H0[jj+1] ;
         H0[jj]   = rfac*rval - ifac*ival ;
         H0[jj+1] = rfac*ival + ifac*rval ;
         rval = H0[jj] ; ival = H0[jj+1] ;
         sum += rval*rval + ival*ival ;
      }
      nops += 10*(n-1) + 5 ;
      beta0 = 2./sum ; 
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n %% sum = %12.4e, beta0 = %12.4e", 
              sum, beta0) ;
   }
/*
   ---------------------------------------------
   set the first entry of the transformed column
   ---------------------------------------------
*/
   if ( type == SPOOLES_REAL ) {
      H0[0] = y0real ;
   } else if ( type == SPOOLES_COMPLEX ) {
      H0[0] = y0real ; H0[1] = y0imag ;
   }
}
*pbeta0 = beta0 ;
/*
fprintf(msgFile, "\n H0") ;
DVfprintf(msgFile, n, H0) ;
*/

return(nops) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   compute W0 = v^H * A

   created -- 98may30, cca
   -----------------------
*/
static double
computeW1 (
   A2       *mtxA,
   double   H0[],
   double   W0[],
   int      msglvl,
   FILE     *msgFile
) {
double   nops ;
int      inc1, inc2, ncolA, nrowA ;

if ( msglvl > 5 ) {
   fprintf(msgFile, "\n %% inside computeW1, nrow %d, ncol %d",
           mtxA->n1, mtxA->n2) ;
}

nrowA = mtxA->n1 ;
ncolA = mtxA->n2 ;
inc1  = mtxA->inc1 ;
inc2  = mtxA->inc2 ;
if ( inc1 != 1 && inc2 != 1 ) {
   fprintf(stderr, "\n error in computeW1"
           "\n inc1 = %d, inc2 = %d\n", inc1, inc2) ;
   exit(-1) ;
}
nops  = 0.0 ;
if ( A2_IS_REAL(mtxA) ) {
   int      irow, jcol ;

   if ( inc1 == 1 ) {
      double   sums[3] ;
      double   *colA0, *colA1, *colA2 ;
/*
      ----------------------------
      A is column major, 
      compute W(j) = H0^T * A(*,j)
      ----------------------------
*/
      for ( jcol = 0 ; jcol < ncolA - 2 ; jcol += 3 ) {
         colA0 = A2_column(mtxA, jcol)   ;
         colA1 = A2_column(mtxA, jcol+1) ;
         colA2 = A2_column(mtxA, jcol+2) ;
         DVdot13(nrowA, H0, colA0, colA1, colA2, sums) ;
         W0[jcol]   = sums[0] ;
         W0[jcol+1] = sums[1] ;
         W0[jcol+2] = sums[2] ;
         nops += 6*nrowA ;
      }
      if ( jcol == ncolA - 2 ) {
         colA0 = A2_column(mtxA, jcol)   ;
         colA1 = A2_column(mtxA, jcol+1) ;
         DVdot12(nrowA, H0, colA0, colA1, sums) ;
         W0[jcol]   = sums[0] ;
         W0[jcol+1] = sums[1] ;
         nops += 4*nrowA ;
      } else if ( jcol == ncolA - 1 ) {
         colA0 = A2_column(mtxA, jcol)   ;
         DVdot11(nrowA, H0, colA0, sums) ;
         W0[jcol] = sums[0] ;
         nops += 2*nrowA ;
      }
   } else {
      double   alpha[3] ;
      double   *rowA0, *rowA1, *rowA2 ;
/*
      -------------------------------
      A is row major
      compute W := W + H0(j) * A(j,*)
      -------------------------------
*/
      DVzero(ncolA, W0) ;
      for ( irow = 0 ; irow < nrowA - 2 ; irow += 3 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         rowA2 = A2_row(mtxA, irow+2) ;
         alpha[0] = H0[irow]   ; 
         alpha[1] = H0[irow+1] ; 
         alpha[2] = H0[irow+2] ; 
         DVaxpy13(ncolA, W0, alpha, rowA0, rowA1, rowA2) ;
         nops += 6*ncolA ;
      }
      if ( irow == nrowA - 2 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         alpha[0] = H0[irow]   ; 
         alpha[1] = H0[irow+1] ; 
         DVaxpy12(ncolA, W0, alpha, rowA0, rowA1) ;
         nops += 4*ncolA ;
      } else if ( irow == nrowA - 1 ) {
         rowA0 = A2_row(mtxA, irow) ;
         alpha[0] = H0[irow]   ; 
         DVaxpy11(ncolA, W0, alpha, rowA0) ;
         nops += 2*ncolA ;
      }
   }
} else if ( A2_IS_COMPLEX(mtxA) ) {
   int      irow, jcol ;

   if ( inc1 == 1 ) {
      double   sums[6] ;
      double   *colA0, *colA1, *colA2 ;
/*
      ----------------------------
      A is column major
      compute W(j) = H0^H * A(*,j)
      ----------------------------
*/
      for ( jcol = 0 ; jcol < ncolA - 2 ; jcol += 3 ) {
         colA0 = A2_column(mtxA, jcol)   ;
         colA1 = A2_column(mtxA, jcol+1) ;
         colA2 = A2_column(mtxA, jcol+2) ;
         ZVdotC13(nrowA, H0, colA0, colA1, colA2, sums) ;
         W0[2*jcol]     = sums[0] ; W0[2*jcol+1]     = sums[1] ;
         W0[2*(jcol+1)] = sums[2] ; W0[2*(jcol+1)+1] = sums[3] ;
         W0[2*(jcol+2)] = sums[4] ; W0[2*(jcol+2)+1] = sums[5] ;
         nops += 24*nrowA ;
      }
      if ( jcol == ncolA - 2 ) {
         colA0 = A2_column(mtxA, jcol)   ;
         colA1 = A2_column(mtxA, jcol+1) ;
         ZVdotC12(nrowA, H0, colA0, colA1, sums) ;
         W0[2*jcol]     = sums[0] ; W0[2*jcol+1]     = sums[1] ;
         W0[2*(jcol+1)] = sums[2] ; W0[2*(jcol+1)+1] = sums[3] ;
         nops += 16*nrowA ;
      } else if ( jcol == ncolA - 1 ) {
         colA0 = A2_column(mtxA, jcol)   ;
         ZVdotC11(nrowA, H0, colA0, sums) ;
         W0[2*jcol]     = sums[0] ; W0[2*jcol+1]     = sums[1] ;
         nops += 8*nrowA ;
      }
   } else {
      double   alpha[6] ;
      double   *rowA0, *rowA1, *rowA2 ;
/*
      ---------------------------------
      A is row major
      compute W := W + H0(j)^H * A(j,*)
      ---------------------------------
*/
      DVzero(2*ncolA, W0) ;
      for ( irow = 0 ; irow < nrowA - 2 ; irow += 3 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         rowA2 = A2_row(mtxA, irow+2) ;
         alpha[0] = H0[2*irow]     ; alpha[1] = -H0[2*irow+1]   ; 
         alpha[2] = H0[2*(irow+1)] ; alpha[3] = -H0[2*(irow+1)+1] ;
         alpha[4] = H0[2*(irow+2)] ; alpha[5] = -H0[2*(irow+2)+1] ;
         ZVaxpy13(ncolA, W0, alpha, rowA0, rowA1, rowA2) ;
         nops += 24*ncolA ;
      }
      if ( irow == nrowA - 2 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         alpha[0] = H0[2*irow]     ; alpha[1] = -H0[2*irow+1]   ; 
         alpha[2] = H0[2*(irow+1)] ; alpha[3] = -H0[2*(irow+1)+1] ;
         ZVaxpy12(ncolA, W0, alpha, rowA0, rowA1) ;
         nops += 16*ncolA ;
      } else if ( irow == nrowA - 1 ) {
         rowA0 = A2_row(mtxA, irow) ;
         alpha[0] = H0[2*irow] ; alpha[1] = -H0[2*irow+1] ; 
         ZVaxpy11(ncolA, W0, alpha, rowA0) ;
         nops += 8*ncolA ;
      }
   }
}
return(nops) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   compute A := A - H0 * beta0 * W0

   created -- 98may30, cca
   ---------------------------------
*/
static double
updateA1 (
   A2       *mtxA,
   double   H0[],
   double   beta0,
   double   W0[],
   int      msglvl,
   FILE     *msgFile
) {
double   nops ;
int      inc1, inc2, ncolA, nrowA ;

if ( msglvl > 5 ) {
   fprintf(msgFile, "\n %% inside updateA1, nrow %d, ncol %d",
           mtxA->n1, mtxA->n2) ;
}

nrowA = mtxA->n1 ;
ncolA = mtxA->n2 ;
inc1  = mtxA->inc1 ;
inc2  = mtxA->inc2 ;
nops  = 0.0 ;
if ( A2_IS_REAL(mtxA) ) {
   int      irow, jcol ;

   if ( inc1 == 1 ) {
      double   alpha[3] ;
      double   *colA0, *colA1, *colA2 ;
/*
      -----------------------------------------
      A is column major
      compute A(:,jcol) -= beta * W0(jcol) * H0
      -----------------------------------------
*/
      for ( jcol = 0 ; jcol < ncolA - 2 ; jcol += 3 ) {
         colA0 = A2_column(mtxA, jcol) ;
         colA1 = A2_column(mtxA, jcol+1) ;
         colA2 = A2_column(mtxA, jcol+2) ;
         alpha[0] = -beta0 * W0[jcol] ;
         alpha[1] = -beta0 * W0[jcol+1] ;
         alpha[2] = -beta0 * W0[jcol+2] ;
         DVaxpy31(nrowA, colA0, colA1, colA2, alpha, H0) ;
         nops += 6*nrowA ;
      }
      if ( jcol == ncolA - 2 ) {
         colA0 = A2_column(mtxA, jcol) ;
         colA1 = A2_column(mtxA, jcol+1) ;
         alpha[0] = -beta0 * W0[jcol] ;
         alpha[1] = -beta0 * W0[jcol+1] ;
         DVaxpy21(nrowA, colA0, colA1, alpha, H0) ;
         nops += 4*nrowA ;
      } else if ( jcol == ncolA - 1 ) {
         colA0 = A2_column(mtxA, jcol) ;
         alpha[0] = -beta0 * W0[jcol] ;
         DVaxpy11(nrowA, colA0, alpha, H0) ;
         nops += 2*nrowA ;
      }
   } else {
      double   alpha[3] ;
      double   *rowA0, *rowA1, *rowA2 ;
/*
      -----------------------------------------
      A is row major
      compute A(irow,:) -= H0[irow]*beta0*W0(:)
      -----------------------------------------
*/
      for ( irow = 0 ; irow < nrowA - 2 ; irow += 3 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         rowA2 = A2_row(mtxA, irow+2) ;
         alpha[0] = -beta0 * H0[irow] ;
         alpha[1] = -beta0 * H0[irow+1] ;
         alpha[2] = -beta0 * H0[irow+2] ;
         DVaxpy31(ncolA, rowA0, rowA1, rowA2, alpha, W0) ;
         nops += 6*ncolA + 3 ;
      }
      if ( irow == nrowA - 2 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         alpha[0] = -beta0 * H0[irow] ;
         alpha[1] = -beta0 * H0[irow+1] ;
         DVaxpy21(ncolA, rowA0, rowA1, alpha, W0) ;
         nops += 4*ncolA + 2 ;
      } else if ( irow == nrowA - 1 ) {
         rowA0 = A2_row(mtxA, irow) ;
         alpha[0] = -beta0 * H0[irow] ;
         DVaxpy11(ncolA, rowA0, alpha, W0) ;
         nops += 2*ncolA + 1 ;
      }
   }
} else if ( A2_IS_COMPLEX(mtxA) ) {
   int      irow, jcol ;

   if ( inc1 == 1 ) {
      double   alpha[6] ;
      double   *colA0, *colA1, *colA2 ;
/*
      -----------------------------------------
      A is column major
      compute A(:,jcol) -= beta * W0(jcol) * H0
      -----------------------------------------
*/
      for ( jcol = 0 ; jcol < ncolA - 2 ; jcol += 3 ) {
         colA0 = A2_column(mtxA, jcol) ;
         colA1 = A2_column(mtxA, jcol+1) ;
         colA2 = A2_column(mtxA, jcol+2) ;
         alpha[0] = -beta0 * W0[2*jcol] ;
         alpha[1] = -beta0 * W0[2*jcol+1] ;
         alpha[2] = -beta0 * W0[2*(jcol+1)] ;
         alpha[3] = -beta0 * W0[2*(jcol+1)+1] ;
         alpha[4] = -beta0 * W0[2*(jcol+2)] ;
         alpha[5] = -beta0 * W0[2*(jcol+2)+1] ;
         ZVaxpy31(nrowA, colA0, colA1, colA2, alpha, H0) ;
         nops += 24*nrowA ;
      }
      if ( jcol == ncolA - 2 ) {
         colA0 = A2_column(mtxA, jcol) ;
         colA1 = A2_column(mtxA, jcol+1) ;
         alpha[0] = -beta0 * W0[2*jcol] ;
         alpha[1] = -beta0 * W0[2*jcol+1] ;
         alpha[2] = -beta0 * W0[2*(jcol+1)] ;
         alpha[3] = -beta0 * W0[2*(jcol+1)+1] ;
         ZVaxpy21(nrowA, colA0, colA1, alpha, H0) ;
         nops += 16*nrowA ;
      } else if ( jcol == ncolA - 1 ) {
         colA0 = A2_column(mtxA, jcol) ;
         alpha[0] = -beta0 * W0[2*jcol] ;
         alpha[1] = -beta0 * W0[2*jcol+1] ;
         ZVaxpy11(nrowA, colA0, alpha, H0) ;
         nops += 8*nrowA ;
      }
   } else {
      double   alpha[6] ;
      double   *rowA0, *rowA1, *rowA2 ;
/*
      -----------------------------------------
      A is row major
      compute A(irow,:) -= H0[irow]*beta0*W0(:)
      -----------------------------------------
*/
      for ( irow = 0 ; irow < nrowA - 2 ; irow += 3 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         rowA2 = A2_row(mtxA, irow+2) ;
         alpha[0] = -beta0 * H0[2*irow] ;
         alpha[1] = -beta0 * H0[2*irow+1] ;
         alpha[2] = -beta0 * H0[2*(irow+1)] ;
         alpha[3] = -beta0 * H0[2*(irow+1)+1] ;
         alpha[4] = -beta0 * H0[2*(irow+2)] ;
         alpha[5] = -beta0 * H0[2*(irow+2)+1] ;
         ZVaxpy31(ncolA, rowA0, rowA1, rowA2, alpha, W0) ;
         nops += 24*ncolA + 12 ;
      }
      if( irow == nrowA - 2 ) {
         rowA0 = A2_row(mtxA, irow) ;
         rowA1 = A2_row(mtxA, irow+1) ;
         alpha[0] = -beta0 * H0[2*irow] ;
         alpha[1] = -beta0 * H0[2*irow+1] ;
         alpha[2] = -beta0 * H0[2*(irow+1)] ;
         alpha[3] = -beta0 * H0[2*(irow+1)+1] ;
         ZVaxpy21(ncolA, rowA0, rowA1, alpha, W0) ;
         nops += 16*ncolA + 8 ;
      } else if( irow == nrowA - 1 ) {
         rowA0 = A2_row(mtxA, irow) ;
         alpha[0] = -beta0 * H0[2*irow] ;
         alpha[1] = -beta0 * H0[2*irow+1] ;
         ZVaxpy11(ncolA, rowA0, alpha, W0) ;
         nops += 8*ncolA + 4 ;
      }
   }
}
return(nops) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   A contains the following data from the A = QR factorization

   A(1:ncolA,1:ncolA) = R
   A(j+1:nrowA,j) is v_j, the j-th householder vector, 
       where v_j[j] = 1.0

   NOTE: A and Q must be column major

   created -- 98dec10, cca
   -----------------------------------------------------------
*/
void
A2_computeQ (
   A2     *Q,
   A2     *A,
   DV     *workDV,
   int    msglvl,
   FILE   *msgFile
) {
double   *betas ;
int      irowA, jcolA, ncolA, nrowA ;
/*
   ---------------
   check the input
   ---------------
*/
if ( Q == NULL ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n Q is NULL\n") ;
   exit(-1) ;
}
if ( A == NULL ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n A is NULL\n") ;
   exit(-1) ;
}
if ( workDV == NULL ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n workDV is NULL\n") ;
   exit(-1) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n msglvl > 0 and msgFile is NULL\n") ;
   exit(-1) ;
}
nrowA = A2_nrow(A) ;
ncolA = A2_ncol(A) ;
if ( nrowA <= 0 ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n nrowA = %d\n", nrowA) ;
   exit(-1) ;
}
if ( ncolA <= 0 ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n ncolA = %d\n", nrowA) ;
   exit(-1) ;
}
if ( nrowA != A2_nrow(Q) ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n nrowA = %d, nrowQ = %d\n", nrowA, A2_nrow(Q)) ;
   exit(-1) ;
}
if ( ncolA != A2_ncol(Q) ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n ncolA = %d, ncolQ = %d\n", ncolA, A2_ncol(Q)) ;
   exit(-1) ;
}
switch ( A->type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n invalid type for A\n") ;
   exit(-1) ;
}
if ( A->type != Q->type ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n A->type = %d, Q->type = %d\n", A->type, Q->type) ;
   exit(-1) ;
}
if ( A2_inc1(A) != 1 ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n A->inc1 = %d \n", A2_inc1(A)) ; 
   exit(-1) ;
}
if ( A2_inc1(Q) != 1 ) {
   fprintf(stderr, "\n fatal error in A2_computeQ()"
           "\n Q->inc1 = %d, \n", A2_inc1(Q)) ;
   exit(-1) ;
}
/*
   --------------------------------------------------
   compute the beta values, beta_j = 2./(V_j^H * V_j)
   --------------------------------------------------
*/
DV_setSize(workDV, ncolA) ;
betas = DV_entries(workDV) ;
if ( A2_IS_REAL(A) ) {
   int   irowA, jcolA ;
   double   sum ;
   double   *colA ;

   for ( jcolA = 0 ; jcolA < ncolA ; jcolA++ ) {
      sum = 1.0 ;
      colA = A2_column(A, jcolA) ;
      for ( irowA = jcolA + 1 ; irowA < nrowA ; irowA++ ) {
         sum += colA[irowA] * colA[irowA] ;
      }
      betas[jcolA] = 2./sum ;
   }
} else {
   double   ival, rval, sum ;
   double   *colA ;

   for ( jcolA = 0 ; jcolA < ncolA ; jcolA++ ) {
      sum = 1.0 ;
      colA = A2_column(A, jcolA) ;
      for ( irowA = jcolA + 1 ; irowA < nrowA ; irowA++ ) {
         rval = colA[2*irowA] ; ival = colA[2*irowA+1] ;
         sum += rval*rval + ival*ival ;
      }
      betas[jcolA] = 2./sum ;
   }
}
/*
   -------------------------------------------
   loop over the number of householder vectors
   -------------------------------------------
*/
for ( jcolA = 0 ; jcolA < ncolA ; jcolA++ ) {
   double   *V, *X ;
   int      jcolV ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n %% jcolA = %d", jcolA) ;
      fflush(msgFile) ;
   }
/*
   ------------------
   set X[] to e_jcolA
   ------------------
*/
   X = A2_column(Q, jcolA) ;
   if ( A2_IS_REAL(Q) ) {
      DVzero(nrowA, X) ;
      X[jcolA] = 1.0 ;
   } else {
      DVzero(2*nrowA, X) ;
      X[2*jcolA] = 1.0 ;
   }
   for ( jcolV = jcolA ; jcolV >= 0 ; jcolV-- ) {
      double   beta ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n   %% jcolV = %d", jcolV) ;
         fflush(msgFile) ;
      }
/*
      -----------------------------------------------------
      update X = (I - beta_jcolV * V_jcolV * V_jcolV^T)X
               = X - beta_jcolV * V_jcolV * V_jcolV^T * X
               = X - (beta_jcolV * V_jcolV^T * X) * V_jcolV 
      -----------------------------------------------------
*/
      V = A2_column(A, jcolV) ;
      beta = betas[jcolV] ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n   %% beta = %12.4e", beta) ;
         fflush(msgFile) ;
      }
      if ( A2_IS_REAL(Q) ) {
         double   fac, sum = X[jcolV] ;
         int      irow ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, 
                       "\n      %% V[%d] = %12.4e, X[%d] = %12.4e",
                       irow, V[irow], irow, X[irow]) ;
               fflush(msgFile) ;
            }
            sum += V[irow] * X[irow] ;
         }
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n   %% sum = %12.4e", sum) ;
            fflush(msgFile) ;
         }
         fac = beta * sum ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n   %% fac = %12.4e", fac) ;
            fflush(msgFile) ;
         }
         X[jcolV] -= fac ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            X[irow] -= fac * V[irow] ;
         }
      } else {
         double   rfac, ifac, rsum = X[2*jcolV], isum = X[2*jcolV+1] ;
         int      irow ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            double   Vi, Vr, Xi, Xr ;
            Vr = V[2*irow] ; Vi = V[2*irow+1] ;
            Xr = X[2*irow] ; Xi = X[2*irow+1] ;
            rsum += Vr*Xr + Vi*Xi ;
            isum += Vr*Xi - Vi*Xr ;
         }
         rfac = beta * rsum ;
         ifac = beta * isum ;
         X[2*jcolV]   -= rfac ;
         X[2*jcolV+1] -= ifac ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            double   Vi, Vr ;
            Vr = V[2*irow] ; Vi = V[2*irow+1] ;
            X[2*irow]   -= rfac*Vr - ifac*Vi ;
            X[2*irow+1] -= rfac*Vi + ifac*Vr ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   A contains the following data from the A = QR factorization

   A(1:ncolA,1:ncolA) = R
   A(j+1:nrowA,j) is v_j, the j-th householder vector, 
       where v_j[j] = 1.0

   we compute Y = Q^T X when A is real
          and Y = Q^H X when A is complex

   NOTE: A, Y and X must be column major.
   NOTE: Y and X can be the same object,
         in which case X is overwritten with Y

   created -- 98dec10, cca
   -----------------------------------------------------------
*/
void
A2_applyQT (
   A2     *Y,
   A2     *A,
   A2     *X,
   DV     *workDV,
   int    msglvl,
   FILE   *msgFile
) {
double   *betas ;
int      irowA, jcolA, jcolX, ncolA, ncolX, nrowA ;
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n A is NULL\n") ;
   exit(-1) ;
}
if ( X == NULL ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n X is NULL\n") ;
   exit(-1) ;
}
if ( workDV == NULL ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n workDV is NULL\n") ;
   exit(-1) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n msglvl > 0 and msgFile is NULL\n") ;
   exit(-1) ;
}
nrowA = A2_nrow(A) ;
ncolA = A2_ncol(A) ;
ncolX = A2_ncol(X) ;
if ( nrowA <= 0 ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n nrowA = %d\n", nrowA) ;
   exit(-1) ;
}
if ( ncolA <= 0 ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n ncolA = %d\n", nrowA) ;
   exit(-1) ;
}
if ( nrowA != A2_nrow(X) ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n nrowA = %d, nrowX = %d\n", nrowA, A2_nrow(X)) ;
   exit(-1) ;
}
switch ( A->type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n invalid type for A\n") ;
   exit(-1) ;
}
if ( A->type != X->type ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n A->type = %d, X->type = %d\n", A->type, X->type) ;
   exit(-1) ;
}
if ( A2_inc1(A) != 1 ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n A->inc1 = %d \n", A2_inc1(A)) ; 
   exit(-1) ;
}
if ( A2_inc1(X) != 1 ) {
   fprintf(stderr, "\n fatal error in A2_applyQT()"
           "\n X->inc1 = %d, \n", A2_inc1(X)) ;
   exit(-1) ;
}
/*
   --------------------------------------------------
   compute the beta values, beta_j = 2./(V_j^H * V_j)
   --------------------------------------------------
*/
DV_setSize(workDV, ncolA) ;
betas = DV_entries(workDV) ;
if ( A2_IS_REAL(A) ) {
   int   irowA, jcolA ;
   double   sum ;
   double   *colA ;

   for ( jcolA = 0 ; jcolA < ncolA ; jcolA++ ) {
      sum = 1.0 ;
      colA = A2_column(A, jcolA) ;
      for ( irowA = jcolA + 1 ; irowA < nrowA ; irowA++ ) {
         sum += colA[irowA] * colA[irowA] ;
      }
      betas[jcolA] = 2./sum ;
   }
} else {
   double   ival, rval, sum ;
   double   *colA ;

   for ( jcolA = 0 ; jcolA < ncolA ; jcolA++ ) {
      sum = 1.0 ;
      colA = A2_column(A, jcolA) ;
      for ( irowA = jcolA + 1 ; irowA < nrowA ; irowA++ ) {
         rval = colA[2*irowA] ; ival = colA[2*irowA+1] ;
         sum += rval*rval + ival*ival ;
      }
      betas[jcolA] = 2./sum ;
   }
}
/*
   ------------------------------------------
   loop over the number of columns in X and Y
   ------------------------------------------
*/
for ( jcolX = 0 ; jcolX < ncolX ; jcolX++ ) {
   double   *V, *colX, *colY ;
   int      jcolV ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n %% jcolX = %d", jcolX) ;
      fflush(msgFile) ;
   }
/*
   -------------------------------
   copy X(:,jcolX) into Y(:,jcolX)
   -------------------------------
*/
   colY = A2_column(Y, jcolX) ;
   colX = A2_column(X, jcolX) ;
   if ( A2_IS_REAL(A) ) {
      DVcopy(nrowA, colY, colX) ;
   } else {
      DVcopy(2*nrowA, colY, colX) ;
   }
   for ( jcolV = 0 ; jcolV < ncolA ; jcolV++ ) {
      double   beta ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n   %% jcolV = %d", jcolV) ;
         fflush(msgFile) ;
      }
/*
      ------------------------------------------------------------
      update colY = (I - beta_jcolV * V_jcolV * V_jcolV^T)colY
                  = colY - beta_jcolV * V_jcolV * V_jcolV^T * colY
                  = colY - (beta_jcolV * V_jcolV^T * Y) * V_jcolV 
      ------------------------------------------------------------
*/
      V = A2_column(A, jcolV) ;
      beta = betas[jcolV] ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n   %% beta = %12.4e", beta) ;
         fflush(msgFile) ;
      }
      if ( A2_IS_REAL(A) ) {
         double   fac, sum = colY[jcolV] ;
         int      irow ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, 
                       "\n      %% V[%d] = %12.4e, X[%d] = %12.4e",
                       irow, V[irow], irow, colY[irow]) ;
               fflush(msgFile) ;
            }
            sum += V[irow] * colY[irow] ;
         }
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n   %% sum = %12.4e", sum) ;
            fflush(msgFile) ;
         }
         fac = beta * sum ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n   %% fac = %12.4e", fac) ;
            fflush(msgFile) ;
         }
         colY[jcolV] -= fac ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            colY[irow] -= fac * V[irow] ;
         }
      } else {
         double   rfac, ifac, 
                  rsum = colY[2*jcolV], isum = colY[2*jcolV+1] ;
         int      irow ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            double   Vi, Vr, Yi, Yr ;
            Vr = V[2*irow] ; Vi = V[2*irow+1] ;
            Yr = colY[2*irow] ; Yi = colY[2*irow+1] ;
            rsum += Vr*Yr + Vi*Yi ;
            isum += Vr*Yi - Vi*Yr ;
         }
         rfac = beta * rsum ;
         ifac = beta * isum ;
         colY[2*jcolV]   -= rfac ;
         colY[2*jcolV+1] -= ifac ;
         for ( irow = jcolV + 1 ; irow < nrowA ; irow++ ) {
            double   Vi, Vr ;
            Vr = V[2*irow] ; Vi = V[2*irow+1] ;
            colY[2*irow]   -= rfac*Vr - ifac*Vi ;
            colY[2*irow+1] -= rfac*Vi + ifac*Vr ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
