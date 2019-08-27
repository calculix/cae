/*  ZV.c  */

#include "../Utilities.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to return the absolute value of (areal,aimag)

   created -- 98jan24, cca
   --------------------------------------------------------
*/
double
Zabs (
   double  real,
   double  imag
) {
double   abs, val ;

if ( real == 0.0 ) {
   abs =  fabs(imag) ;
} else if ( imag == 0.0 ) {
   abs =  fabs(real) ;
} else if ( real >= imag ) {
   val = imag/real ;
   abs = fabs(real)*sqrt(1 + val*val) ;
} else {
   val = real/imag ;
   abs = fabs(imag)*sqrt(1 + val*val) ;
}
return(abs) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- given (areal,aimag), 
              compute (breal,bimag) = 1/(areal,aimag), 

   put breal into *pbreal
   put bimag into *pbimag

   created -- 98jan23, cca
   ---------------------------------------------------
*/
int
Zrecip ( 
   double   areal,
   double   aimag,
   double   *pbreal,
   double   *pbimag
) {
double   bimag, breal, fac ;
if ( areal == 0.0 && aimag == 0.0 ) {
   return(0) ;
}
if ( fabs(areal) >= fabs(aimag) ) {
   fac  = aimag/areal ;
   breal = 1./(areal + fac*aimag) ;
   bimag = -fac * breal ;
} else {
   fac = areal/aimag ;
   bimag = -1./(aimag + fac*areal) ;
   breal = -fac * bimag ;
}
*pbreal = breal ;
*pbimag = bimag ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   given [ (areal,aimag) (breal,bimag) ]
         [ (creal,cimag) (dreal,dimag) ]
   compute [ (ereal,eimag) (freal,fimag) ]
           [ (greal,gimag) (hreal,himag) ]
   where
   I = [ (areal,aimag) (breal,bimag) ] * [ (ereal,eimag) (freal,fimag) ]
       [ (creal,cimag) (dreal,dimag) ]   [ (greal,gimag) (hreal,himag) ]

   note, any of {pereal, peimag, pfreal, pfimag, pgreal, pgimag,
   phreal, phimag} can be NULL. if not NULL, their value is filled.
   this is useful when the input matrix is symmetric or hermitian.

   created -- 98jan23, cca
   ---------------------------------------------------------------------
*/
int
Zrecip2 (
   double areal,
   double aimag,
   double breal,
   double bimag,
   double creal,
   double cimag,
   double dreal,
   double dimag,
   double *pereal,
   double *peimag,
   double *pfreal,
   double *pfimag,
   double *pgreal,
   double *pgimag,
   double *phreal,
   double *phimag
) {
double   xreal, ximag, yreal, yimag ;
/*
   -------------------------------
   step one, compute x = a*d - b*c
   -------------------------------
*/
xreal = areal*dreal - aimag*dimag - breal*creal + bimag*cimag ;
ximag = areal*dimag + aimag*dreal - breal*cimag - bimag*creal ;
/*
   -------------------------
   step two, compute y = 1/x
   -------------------------
*/
Zrecip(xreal, ximag, &yreal, &yimag) ;
/*
   -----------------------
   step three, 
   [ e f ] = y * [  d -b ]
   [ g h ]       [ -c  a ]
   -----------------------
*/
if ( pereal != NULL ) {
   *pereal =  dreal*yreal - dimag*yimag ;
}
if ( peimag != NULL ) {
   *peimag =  dreal*yimag + dimag*yreal ;
}
if ( pfreal != NULL ) {
   *pfreal = -breal*yreal + bimag*yimag ;
}
if ( pfimag != NULL ) {
   *pfimag = -breal*yimag - bimag*yreal ;
}
if ( pgreal != NULL ) {
   *pgreal = -creal*yreal + cimag*yimag ;
}
if ( pgimag != NULL ) {
   *pgimag = -creal*yimag - cimag*yreal ;
}
if ( phreal != NULL ) {
   *phreal =  areal*yreal - aimag*yimag ;
}
if ( phimag != NULL ) {
   *phimag =  areal*yimag + aimag*yreal ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to allocate, initialize and
              return a complex vector x[]

   n -- length of the complex vector (2*n double's)
   
   x[ii] = (real, imag) for 0 <= ii < n

   created -- 98jan23, cca
   ------------------------------------------------
*/
double *
ZVinit (
   int      n,
   double   real,
   double   imag
) {
double   *x ;
int      ii, jj ;

if ( n <= 0 ) {
   fprintf(stderr, "\n fatal error in ZVinit(%d,%f,%f)"
           "\n bad input\n", n, real, imag) ;
   exit(-1) ;
}
ALLOCATE(x, double, 2*n) ;
for ( ii = jj = 0 ; ii < n ; ii++, jj += 2 ) {
   x[jj]   = real ;
   x[jj+1] = imag ;
}
return(x) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to perform a complex dot product

   (*prdot,*pidot) = y^T x

   where y and x are complex

   created -- 98apr15, cca
   -------------------------------------------
*/
void
ZVdotU ( 
   int      size, 
   double   y[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) {
double   isum, rsum, ximag, xreal, yimag, yreal ;
int      ii, jj ;

if (  size < 0 || y == NULL || x == NULL 
   || prdot == NULL || pidot == NULL ) {
   fprintf(stderr, "\n fatal error in ZVdotU(%d,%p,%p,%p,%p)"
           "\n bad input\n", size, y, x, prdot, pidot) ;
   exit(-1) ;
} 
isum = rsum = 0.0 ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   xreal = x[jj] ;
   ximag = x[jj+1] ;
   yreal = y[jj] ;
   yimag = y[jj+1] ;
   rsum += xreal*yreal - ximag*yimag ;
   isum += xreal*yimag + ximag*yreal ;
}
*prdot = rsum ;
*pidot = isum ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to perform a conjugated complex dot product

   (*prdot,*pidot) = (conjugate(y))^T x = y^H x

   where y and x are complex

   created -- 98apr15, cca
   ------------------------------------------------------
*/
void
ZVdotC ( 
   int      size, 
   double   y[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) {
double   isum, rsum, ximag, xreal, yimag, yreal ;
int      ii, jj ;

if (  size < 0 || y == NULL || x == NULL 
   || prdot == NULL || pidot == NULL ) {
   fprintf(stderr, "\n fatal error in ZVdotC(%d,%p,%p,%p,%p)"
           "\n bad input\n", size, y, x, prdot, pidot) ;
   exit(-1) ;
} 
isum = rsum = 0.0 ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   xreal = x[jj] ;
   ximag = x[jj+1] ;
   yreal = y[jj] ;
   yimag = y[jj+1] ;
   rsum +=   xreal*yreal + ximag*yimag ;
   isum += - xreal*yimag + ximag*yreal ;
}
*prdot = rsum ;
*pidot = isum ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to perform a indexed complex dot product

   (*prdot,*pidot) = sum_k y[index[k]]*x[k]

   where y and x are complex

   created -- 98apr15, cca
   ---------------------------------------------------
*/
void
ZVdotiU ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) {
double   isum, rsum, ximag, xreal, yimag, yreal ;
int      ii, jj ;

if (  size < 0 || y == NULL || index == NULL || x == NULL 
   || prdot == NULL || pidot == NULL ) {
   fprintf(stderr, "\n fatal error in ZVdotiU(%d,%p,%p,%p,%p,%p)"
           "\n bad input\n", size, y, index, x, prdot, pidot) ;
   exit(-1) ;
} 
isum = rsum = 0.0 ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
/*
   fprintf(stdout, 
           "\n %% ii = %d, jj = %d, kk = %d", ii, jj, index[ii]) ;
   fflush(stdout) ;
*/
   xreal = x[jj] ;
   ximag = x[jj+1] ;
   yreal = y[2*index[ii]] ;
   yimag = y[2*index[ii]+1] ;
   rsum += xreal*yreal - ximag*yimag ;
   isum += xreal*yimag + ximag*yreal ;
}
*prdot = rsum ;
*pidot = isum ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to perform a indexed conjugate complex dot product

   (*prdot,*pidot) = sum_k conjugate(y[index[k]])*x[k]

   where y and x are complex

   created -- 98apr15, cca
   -------------------------------------------------------------
*/
void
ZVdotiC ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[],
   double   *prdot, 
   double   *pidot
) {
double   isum, rsum, ximag, xreal, yimag, yreal ;
int      ii, jj ;

if (  size < 0 || y == NULL || index == NULL || x == NULL 
   || prdot == NULL || pidot == NULL ) {
   fprintf(stderr, "\n fatal error in ZVdotiU(%d,%p,%p,%p,%p,%p)"
           "\n bad input\n", size, y, index, x, prdot, pidot) ;
   exit(-1) ;
} 
isum = rsum = 0.0 ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   xreal = x[jj] ;
   ximag = x[jj+1] ;
   yreal = y[2*index[ii]] ;
   yimag = y[2*index[ii]+1] ;
   rsum +=   xreal*yreal + ximag*yimag ;
   isum += - xreal*yimag + ximag*yreal ;
}
*prdot = rsum ;
*pidot = isum ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- to perform a complex axpy

   y := y + (areal, aimag) * x

   where y and x are complex

   created -- 98jan22, cca
   ------------------------------------
*/
void
ZVaxpy ( 
   int      size, 
   double   y[], 
   double   areal, 
   double   aimag, 
   double   x[]
) {
double   ximag, xreal ;
int      ii, jj ;

if ( size < 0 || y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in ZVaxpy(%d,%p,%f,%f,%p)"
           "\n bad input\n", size, y, areal, aimag, x) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   xreal = x[jj] ;
   ximag = x[jj+1] ;
/*
fprintf(stdout, "\n ii = %d, xreal = %20.12e, ximag = %20.12e",
        ii, xreal, ximag) ;
fprintf(stdout, "\n before   yreal = %20.12e, yimag = %20.12e",
        y[jj], y[jj+1]) ;
*/
   y[jj]   += areal*xreal - aimag*ximag ;
   y[jj+1] += areal*ximag + aimag*xreal ;
/*
fprintf(stdout, "\n after    yreal = %20.12e, yimag = %20.12e",
        y[jj], y[jj+1]) ;
*/
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to scale a double complex vector by (xreal,ximag)

   y := y * (areal, aimag)

   created -- 98jan22, cca
   ------------------------------------------------------------
*/
void
ZVscale ( 
   int      size, 
   double   y[], 
   double   areal, 
   double   aimag 
) {
double   yimag, yreal ;
int      ii, jj ;

if ( size < 0 || y == NULL ) {
   fprintf(stderr, "\n fatal error in ZVscale(%d,%p,%f,%f)"
           "\n bad input\n", size, y, areal, aimag) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   yreal = y[jj] ;
   yimag = y[jj+1] ;
   y[jj]   = areal*yreal - aimag*yimag ;
   y[jj+1] = areal*yimag + aimag*yreal ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose --- to print a complex vector to a file

   created -- 98jan23, cca
   -----------------------------------------------
*/
void
ZVfprintf (
   FILE     *fp,
   int      size,
   double   y[]
) {
int      ii, jj ;

if ( size < 0 || y == NULL ) {
   fprintf(stderr, "\n fatal error in ZVfprintf(%p,%d,%p)"
           "\n bad input\n", fp, size, y) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
/*
   if ( ii % 2 == 0 ) {
      fprintf(fp, "\n") ;
   }
*/
   fprintf(fp, "\n < %12.4e, %12.4e >", y[jj], y[jj+1]) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return the minimum absolute value 
   of the entries in a complex vector

   created -- 98jan23, cca
   ----------------------------------
*/
double
ZVminabs (
   int      size,
   double   x[]
) {
double   abs, imag, minabs, real, val ;
int      ii, jj ;

if ( size < 0 || x == NULL ) {
   fprintf(stderr, "\n fatal error in ZVminabs(%d,%p)"
           "\n bad input\n", size, x) ;
   exit(-1) ;
} 
minabs = 0.0 ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   real = fabs(x[jj]) ;
   imag = fabs(x[jj+1]) ;
/*
   fprintf(stdout, "\n i = %d, <%12.4e, %12.4e>", ii, real, imag) ;
*/
   if ( real == 0.0 ) {
      abs =  imag ;
   } else if ( imag == 0.0 ) {
      abs =  real ;
   } else if ( real >= imag ) {
      val = imag/real ;
      abs = real*sqrt(1 + val*val) ;
   } else {
      val = real/imag ;
      abs = imag*sqrt(1 + val*val) ;
   }
/*
   fprintf(stdout, " abs = %12.4e", abs) ;
*/
   if ( ii == 0 || minabs > abs ) {
      minabs = abs ;
   }
}
/*
fprintf(stdout, "\n minabs = %12.4e", minabs) ;
*/

return(minabs) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return the maximum absolute value 
   of the entries in a complex vector

   created -- 98jan23, cca
   ----------------------------------
*/
double
ZVmaxabs (
   int      size,
   double   x[]
) {
double   abs, imag, maxabs, real, val ;
int      ii, jj ;

if ( size < 0 || x == NULL ) {
   fprintf(stderr, "\n fatal error in ZVmaxabs(%d,%p)"
           "\n bad input\n", size, x) ;
   exit(-1) ;
} 
maxabs = 0.0 ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   real = fabs(x[jj]) ;
   imag = fabs(x[jj+1]) ;
/*
   fprintf(stdout, "\n i = %d, <%12.4e, %12.4e>", ii, real, imag) ;
*/
   if ( real == 0.0 ) {
      abs =  imag ;
   } else if ( imag == 0.0 ) {
      abs =  real ;
   } else if ( real >= imag ) {
      val = imag/real ;
      abs = real*sqrt(1 + val*val) ;
   } else {
      val = real/imag ;
      abs = imag*sqrt(1 + val*val) ;
   }
/*
   fprintf(stdout, " abs = %12.4e", abs) ;
*/
   if ( ii == 0 || maxabs < abs ) {
      maxabs = abs ;
   }
}
/*
fprintf(stdout, "\n maxabs = %12.4e", maxabs) ;
*/

return(maxabs) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   copy a complex vector into another
   y[] := x[]

   created -- 98jan23, cca
   ----------------------------------
*/
void
ZVcopy (
   int      size,
   double   y[],
   double   x[]
) {
int      ii, jj ;

if ( size < 0 || y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in ZVcopy(%d,%p,%p)"
           "\n bad input\n", size, y, x) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   y[jj]   = x[jj]   ;
   y[jj+1] = x[jj+1] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   subtract a complex vector from another
   y[] := y[] - x[]

   created -- 98may25, cca
   --------------------------------------
*/
void
ZVsub (
   int      size,
   double   y[],
   double   x[]
) {
int      ii, jj ;

if ( size < 0 || y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in ZVsub(%d,%p,%p)"
           "\n bad input\n", size, y, x) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   y[jj]   -= x[jj]   ;
   y[jj+1] -= x[jj+1] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to perform a complex axpy with two vectors

   z := z + (areal, aimag)*x + (breal, bimag)*y

   where y and x are complex

   created -- 98jan23, cca
   ----------------------------------------------------
*/
void
ZVaxpy2 ( 
   int      size, 
   double   z[], 
   double   areal, 
   double   aimag, 
   double   x[],
   double   breal, 
   double   bimag, 
   double   y[]
) {
double   ximag, xreal, yimag, yreal ;
int      ii, jj ;

if ( size < 0 || y == NULL || x == NULL ) {
   fprintf(stderr, "\n fatal error in ZVaxpy(%d,%p,%f,%f,%p)"
           "\n bad input\n", size, y, areal, aimag, x) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   xreal = x[jj] ;
   ximag = x[jj+1] ;
   yreal = y[jj] ;
   yimag = y[jj+1] ;
/*
fprintf(stdout, "\n ii = %d, xreal = %20.12e, ximag = %20.12e",
        ii, xreal, ximag) ;
fprintf(stdout, "\n          yreal = %20.12e, yimag = %20.12e",
        y[jj], y[jj+1]) ;
fprintf(stdout, "\n before   zreal = %20.12e, zimag = %20.12e",
        z[jj], z[jj+1]) ;
*/
   z[jj]   += areal*xreal - aimag*ximag + breal*yreal - bimag*yimag ;
   z[jj+1] += areal*ximag + aimag*xreal + breal*yimag + bimag*yreal ;
/*
fprintf(stdout, "\n after    zreal = %20.12e, zimag = %20.12e",
        z[jj], z[jj+1]) ;
*/
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to scale a double complex vector by a 2x2 matrix

   [ x ] := [ a b ] [ x ]
   [ y ]    [ c d ] [ y ]

   created -- 98jan23, cca
   ------------------------------------------------------------
*/
void
ZVscale2 ( 
   int      size, 
   double   x[], 
   double   y[], 
   double   areal,
   double   aimag,
   double   breal,
   double   bimag,
   double   creal,
   double   cimag,
   double   dreal,
   double   dimag
) {
double   ximag, xreal, yimag, yreal ;
int      ii, jj ;

if ( size < 0 || x == NULL || y == NULL ) {
   fprintf(stderr, "\n fatal error in ZVscale2(%d,%p,%p,...)"
           "\n bad input\n", size, x, y) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   xreal = x[jj] ;
   ximag = x[jj+1] ;
   yreal = y[jj] ;
   yimag = y[jj+1] ;
   x[jj]   = areal*xreal - aimag*ximag + breal*yreal - bimag*yimag ;
   x[jj+1] = areal*ximag + aimag*xreal + breal*yimag + bimag*yreal ;
   y[jj]   = creal*xreal - cimag*ximag + dreal*yreal - dimag*yimag ;
   y[jj+1] = creal*ximag + cimag*xreal + dreal*yimag + dimag*yreal ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- to gather y[*] = x[index[*]]
 
   created -- 98apr15, cca
   ---------------------------------------
*/
void
ZVgather ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in ZVgather, invalid input"
              "\n size = %d, y = %p, x = %p, index = %p\n",
              size, y, x, index) ;
      exit(-1) ;
   } else {
      int   i, j, k ;
      for ( i = j = 0 ; i < size ; i++, j += 2 ) {
         k      = 2*index[i] ; 
         y[j]   = x[k] ;
         y[j+1] = x[k+1] ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to scatter y[index[*]] = x[*]
 
   created -- 98apr15, cca
   ----------------------------------------
*/
void
ZVscatter ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) {
if ( size > 0 ) {
   if ( y == NULL || x == NULL || index == NULL ) {
      fprintf(stderr, "\n fatal error in ZVscatter, invalid data"
              "\n size = %d, y = %p, index = %p, x = %p\n",
              size, y, index, x) ;
      exit(-1) ; 
   } else {
      int   i, j, k ;
      for ( i = j = 0 ; i < size ; i++, j += 2 ) {
         k      = 2*index[i] ;
         y[k]   = x[j]       ; 
         y[k+1] = x[j+1]     ;
      }
   }
}
return ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 y2 ]^T [ x0 x1 x2]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) {
double   i00, i01, i02, i10, i11, i12, i20, i21, i22,
         r00, r01, r02, r10, r11, r12, r20, r21, r22,
         xi0, xi1, xi2, xr0, xr1, xr2, yi0, yi1, yi2, yr0, yr1, yr2 ;
int      ii, iloc, rloc ;

i00 = i01 = i02 = i10 = i11 = i12 = i20 = i21 = i22 
    = r00 = r01 = r02 = r10 = r11 = r12 = r20 = r21 = r22 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   yr2 = y2[rloc] ; yi2 = y2[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   xr2 = x2[rloc] ; xi2 = x2[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r01 += yr0*xr1 - yi0*xi1 ; i01 += yr0*xi1 + yi0*xr1 ;
   r02 += yr0*xr2 - yi0*xi2 ; i02 += yr0*xi2 + yi0*xr2 ;
   r10 += yr1*xr0 - yi1*xi0 ; i10 += yr1*xi0 + yi1*xr0 ;
   r11 += yr1*xr1 - yi1*xi1 ; i11 += yr1*xi1 + yi1*xr1 ;
   r12 += yr1*xr2 - yi1*xi2 ; i12 += yr1*xi2 + yi1*xr2 ;
   r20 += yr2*xr0 - yi2*xi0 ; i20 += yr2*xi0 + yi2*xr0 ;
   r21 += yr2*xr1 - yi2*xi1 ; i21 += yr2*xi1 + yi2*xr1 ;
   r22 += yr2*xr2 - yi2*xi2 ; i22 += yr2*xi2 + yi2*xr2 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r02 ; sums[ 5] = i02 ;
sums[ 6] = r10 ; sums[ 7] = i10 ;
sums[ 8] = r11 ; sums[ 9] = i11 ;
sums[10] = r12 ; sums[11] = i12 ;
sums[12] = r20 ; sums[13] = i20 ;
sums[14] = r21 ; sums[15] = i21 ;
sums[16] = r22 ; sums[17] = i22 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 y2 ]^T [ x0 x1 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   sums[]
) {
double   i00, i01, i10, i11, i20, i21, 
         r00, r01, r10, r11, r20, r21, 
         xi0, xi1, xr0, xr1, yi0, yi1, yi2, yr0, yr1, yr2 ;
int      ii, iloc, rloc ;

i00 = i01 = i10 = i11 = i20 = i21 
    = r00 = r01 = r10 = r11 = r20 = r21 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   yr2 = y2[rloc] ; yi2 = y2[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r01 += yr0*xr1 - yi0*xi1 ; i01 += yr0*xi1 + yi0*xr1 ;
   r10 += yr1*xr0 - yi1*xi0 ; i10 += yr1*xi0 + yi1*xr0 ;
   r11 += yr1*xr1 - yi1*xi1 ; i11 += yr1*xi1 + yi1*xr1 ;
   r20 += yr2*xr0 - yi2*xi0 ; i20 += yr2*xi0 + yi2*xr0 ;
   r21 += yr2*xr1 - yi2*xi1 ; i21 += yr2*xi1 + yi2*xr1 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r10 ; sums[ 5] = i10 ;
sums[ 6] = r11 ; sums[ 7] = i11 ;
sums[ 8] = r20 ; sums[ 9] = i20 ;
sums[10] = r21 ; sums[11] = i21 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 y2 ]^T [ x0 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   sums[]
) {
double   i00, i10, i20, r00, r10, r20, 
         xi0, xr0, yi0, yi1, yi2, yr0, yr1, yr2 ;
int      ii, iloc, rloc ;

i00 = i10 = i20 = r00 = r10 = r20 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   yr2 = y2[rloc] ; yi2 = y2[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r10 += yr1*xr0 - yi1*xi0 ; i10 += yr1*xi0 + yi1*xr0 ;
   r20 += yr2*xr0 - yi2*xi0 ; i20 += yr2*xi0 + yi2*xr0 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r10 ; sums[ 3] = i10 ;
sums[ 4] = r20 ; sums[ 5] = i20 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 ]^T [ x0 x1 x2]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU23 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) {
double   i00, i01, i02, i10, i11, i12, 
         r00, r01, r02, r10, r11, r12, 
         xi0, xi1, xi2, xr0, xr1, xr2, yi0, yi1, yr0, yr1 ;
int      ii, iloc, rloc ;

i00 = i01 = i02 = i10 = i11 = i12 
    = r00 = r01 = r02 = r10 = r11 = r12 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   xr2 = x2[rloc] ; xi2 = x2[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r01 += yr0*xr1 - yi0*xi1 ; i01 += yr0*xi1 + yi0*xr1 ;
   r02 += yr0*xr2 - yi0*xi2 ; i02 += yr0*xi2 + yi0*xr2 ;
   r10 += yr1*xr0 - yi1*xi0 ; i10 += yr1*xi0 + yi1*xr0 ;
   r11 += yr1*xr1 - yi1*xi1 ; i11 += yr1*xi1 + yi1*xr1 ;
   r12 += yr1*xr2 - yi1*xi2 ; i12 += yr1*xi2 + yi1*xr2 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r02 ; sums[ 5] = i02 ;
sums[ 6] = r10 ; sums[ 7] = i10 ;
sums[ 8] = r11 ; sums[ 9] = i11 ;
sums[10] = r12 ; sums[11] = i12 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 ]^T [ x0 x1 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU22 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   sums[]
) {
double   i00, i01, i10, i11, r00, r01, r10, r11, 
         xi0, xi1, xr0, xr1, yi0, yi1, yr0, yr1 ;
int      ii, iloc, rloc ;

i00 = i01 = i10 = i11 = r00 = r01 = r10 = r11 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r01 += yr0*xr1 - yi0*xi1 ; i01 += yr0*xi1 + yi0*xr1 ;
   r10 += yr1*xr0 - yi1*xi0 ; i10 += yr1*xi0 + yi1*xr0 ;
   r11 += yr1*xr1 - yi1*xi1 ; i11 += yr1*xi1 + yi1*xr1 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r10 ; sums[ 5] = i10 ;
sums[ 6] = r11 ; sums[ 7] = i11 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 ]^T [ x0 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU21 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   sums[]
) {
double   i00, i10, r00, r10, xi0, xr0, yi0, yi1, yr0, yr1 ;
int      ii, iloc, rloc ;

i00 = i10 = r00 = r10 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r10 += yr1*xr0 - yi1*xi0 ; i10 += yr1*xi0 + yi1*xr0 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r10 ; sums[ 3] = i10 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 ]^T [ x0 x1 x2]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU13 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) {
double   i00, i01, i02, r00, r01, r02, 
         xi0, xi1, xi2, xr0, xr1, xr2, yi0, yr0 ;
int      ii, iloc, rloc ;

i00 = i01 = i02 = r00 = r01 = r02 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   xr2 = x2[rloc] ; xi2 = x2[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r01 += yr0*xr1 - yi0*xi1 ; i01 += yr0*xi1 + yi0*xr1 ;
   r02 += yr0*xr2 - yi0*xi2 ; i02 += yr0*xi2 + yi0*xr2 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r02 ; sums[ 5] = i02 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 ]^T [ x0 x1 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU12 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   sums[]
) {
double   i00, i01, r00, r01, xi0, xi1, xr0, xr1, yi0, yr0 ;
int      ii, iloc, rloc ;

i00 = i01 = r00 = r01 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
   r01 += yr0*xr1 - yi0*xi1 ; i01 += yr0*xi1 + yi0*xr1 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 ]^T [ x0 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU11 (
   int      n,
   double   y0[],
   double   x0[],
   double   sums[]
) {
double   i00, r00, xi0, xr0, yi0, yr0 ;
int      ii, iloc, rloc ;

i00 = r00 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   r00 += yr0*xr0 - yi0*xi0 ; i00 += yr0*xi0 + yi0*xr0 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 y2 ]^H [ x0 x1 x2]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) {
double   i00, i01, i02, i10, i11, i12, i20, i21, i22,
         r00, r01, r02, r10, r11, r12, r20, r21, r22,
         xi0, xi1, xi2, xr0, xr1, xr2, yi0, yi1, yi2, yr0, yr1, yr2 ;
int      ii, iloc, rloc ;

i00 = i01 = i02 = i10 = i11 = i12 = i20 = i21 = i22 
    = r00 = r01 = r02 = r10 = r11 = r12 = r20 = r21 = r22 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   yr2 = y2[rloc] ; yi2 = y2[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   xr2 = x2[rloc] ; xi2 = x2[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r01 += yr0*xr1 + yi0*xi1 ; i01 += yr0*xi1 - yi0*xr1 ;
   r02 += yr0*xr2 + yi0*xi2 ; i02 += yr0*xi2 - yi0*xr2 ;
   r10 += yr1*xr0 + yi1*xi0 ; i10 += yr1*xi0 - yi1*xr0 ;
   r11 += yr1*xr1 + yi1*xi1 ; i11 += yr1*xi1 - yi1*xr1 ;
   r12 += yr1*xr2 + yi1*xi2 ; i12 += yr1*xi2 - yi1*xr2 ;
   r20 += yr2*xr0 + yi2*xi0 ; i20 += yr2*xi0 - yi2*xr0 ;
   r21 += yr2*xr1 + yi2*xi1 ; i21 += yr2*xi1 - yi2*xr1 ;
   r22 += yr2*xr2 + yi2*xi2 ; i22 += yr2*xi2 - yi2*xr2 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r02 ; sums[ 5] = i02 ;
sums[ 6] = r10 ; sums[ 7] = i10 ;
sums[ 8] = r11 ; sums[ 9] = i11 ;
sums[10] = r12 ; sums[11] = i12 ;
sums[12] = r20 ; sums[13] = i20 ;
sums[14] = r21 ; sums[15] = i21 ;
sums[16] = r22 ; sums[17] = i22 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 y2 ]^H [ x0 x1 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   sums[]
) {
double   i00, i01, i10, i11, i20, i21, 
         r00, r01, r10, r11, r20, r21, 
         xi0, xi1, xr0, xr1, yi0, yi1, yi2, yr0, yr1, yr2 ;
int      ii, iloc, rloc ;

i00 = i01 = i10 = i11 = i20 = i21 
    = r00 = r01 = r10 = r11 = r20 = r21 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   yr2 = y2[rloc] ; yi2 = y2[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r01 += yr0*xr1 + yi0*xi1 ; i01 += yr0*xi1 - yi0*xr1 ;
   r10 += yr1*xr0 + yi1*xi0 ; i10 += yr1*xi0 - yi1*xr0 ;
   r11 += yr1*xr1 + yi1*xi1 ; i11 += yr1*xi1 - yi1*xr1 ;
   r20 += yr2*xr0 + yi2*xi0 ; i20 += yr2*xi0 - yi2*xr0 ;
   r21 += yr2*xr1 + yi2*xi1 ; i21 += yr2*xi1 - yi2*xr1 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r10 ; sums[ 5] = i10 ;
sums[ 6] = r11 ; sums[ 7] = i11 ;
sums[ 8] = r20 ; sums[ 9] = i20 ;
sums[10] = r21 ; sums[11] = i21 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 y2 ]^H [ x0 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   sums[]
) {
double   i00, i10, i20, r00, r10, r20, 
         xi0, xr0, yi0, yi1, yi2, yr0, yr1, yr2 ;
int      ii, iloc, rloc ;

i00 = i10 = i20 = r00 = r10 = r20 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   yr2 = y2[rloc] ; yi2 = y2[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r10 += yr1*xr0 + yi1*xi0 ; i10 += yr1*xi0 - yi1*xr0 ;
   r20 += yr2*xr0 + yi2*xi0 ; i20 += yr2*xi0 - yi2*xr0 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r10 ; sums[ 3] = i10 ;
sums[ 4] = r20 ; sums[ 5] = i20 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 ]^H [ x0 x1 x2]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC23 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) {
double   i00, i01, i02, i10, i11, i12, 
         r00, r01, r02, r10, r11, r12, 
         xi0, xi1, xi2, xr0, xr1, xr2, yi0, yi1, yr0, yr1 ;
int      ii, iloc, rloc ;

i00 = i01 = i02 = i10 = i11 = i12 
    = r00 = r01 = r02 = r10 = r11 = r12 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   xr2 = x2[rloc] ; xi2 = x2[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r01 += yr0*xr1 + yi0*xi1 ; i01 += yr0*xi1 - yi0*xr1 ;
   r02 += yr0*xr2 + yi0*xi2 ; i02 += yr0*xi2 - yi0*xr2 ;
   r10 += yr1*xr0 + yi1*xi0 ; i10 += yr1*xi0 - yi1*xr0 ;
   r11 += yr1*xr1 + yi1*xi1 ; i11 += yr1*xi1 - yi1*xr1 ;
   r12 += yr1*xr2 + yi1*xi2 ; i12 += yr1*xi2 - yi1*xr2 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r02 ; sums[ 5] = i02 ;
sums[ 6] = r10 ; sums[ 7] = i10 ;
sums[ 8] = r11 ; sums[ 9] = i11 ;
sums[10] = r12 ; sums[11] = i12 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 ]^H [ x0 x1 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC22 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   sums[]
) {
double   i00, i01, i10, i11, r00, r01, r10, r11, 
         xi0, xi1, xr0, xr1, yi0, yi1, yr0, yr1 ;
int      ii, iloc, rloc ;

i00 = i01 = i10 = i11 = r00 = r01 = r10 = r11 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r01 += yr0*xr1 + yi0*xi1 ; i01 += yr0*xi1 - yi0*xr1 ;
   r10 += yr1*xr0 + yi1*xi0 ; i10 += yr1*xi0 - yi1*xr0 ;
   r11 += yr1*xr1 + yi1*xi1 ; i11 += yr1*xi1 - yi1*xr1 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r10 ; sums[ 5] = i10 ;
sums[ 6] = r11 ; sums[ 7] = i11 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 y1 ]^H [ x0 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC21 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   sums[]
) {
double   i00, i10, r00, r10, xi0, xr0, yi0, yi1, yr0, yr1 ;
int      ii, iloc, rloc ;

i00 = i10 = r00 = r10 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   yr1 = y1[rloc] ; yi1 = y1[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r10 += yr1*xr0 + yi1*xi0 ; i10 += yr1*xi0 - yi1*xr0 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r10 ; sums[ 3] = i10 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 ]^H [ x0 x1 x2]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC13 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) {
double   i00, i01, i02, r00, r01, r02, 
         xi0, xi1, xi2, xr0, xr1, xr2, yi0, yr0 ;
int      ii, iloc, rloc ;

i00 = i01 = i02 = r00 = r01 = r02 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   xr2 = x2[rloc] ; xi2 = x2[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r01 += yr0*xr1 + yi0*xi1 ; i01 += yr0*xi1 - yi0*xr1 ;
   r02 += yr0*xr2 + yi0*xi2 ; i02 += yr0*xi2 - yi0*xr2 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;
sums[ 4] = r02 ; sums[ 5] = i02 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 ]^H [ x0 x1 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC12 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   sums[]
) {
double   i00, i01, r00, r01, xi0, xi1, xr0, xr1, yi0, yr0 ;
int      ii, iloc, rloc ;

i00 = i01 = r00 = r01 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   xr1 = x1[rloc] ; xi1 = x1[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
   r01 += yr0*xr1 + yi0*xi1 ; i01 += yr0*xi1 - yi0*xr1 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;
sums[ 2] = r01 ; sums[ 3] = i01 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product

     [ y0 ]^H [ x0 ]

   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC11 (
   int      n,
   double   y0[],
   double   x0[],
   double   sums[]
) {
double   i00, r00, xi0, xr0, yi0, yr0 ;
int      ii, iloc, rloc ;

i00 = r00 = 0.0 ;
for ( ii = rloc = 0, iloc = 1 ; ii < n ; ii++, rloc += 2, iloc += 2 ) {
   yr0 = y0[rloc] ; yi0 = y0[iloc] ;
   xr0 = x0[rloc] ; xi0 = x0[iloc] ;
   r00 += yr0*xr0 + yi0*xi0 ; i00 += yr0*xi0 - yi0*xr0 ;
}
sums[ 0] = r00 ; sums[ 1] = i00 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   purpose -- to zero the vector

   y := 0

   where y is complex

   created -- 98apr25, cca
   -----------------------------
*/
void
ZVzero ( 
   int      size, 
   double   y[] 
) {
int      ii, jj ;

if ( size < 0 || y == NULL ) {
   fprintf(stderr, "\n fatal error in ZVzero(%d,%p)"
           "\n bad input\n", size, y) ;
   exit(-1) ;
} 
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   y[jj] = y[jj+1] = 0.0 ;
}
return ; }

/*--------------------------------------------------------------------*/
