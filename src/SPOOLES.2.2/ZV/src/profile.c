/*  profile  */

#include "../ZV.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   to fill xDV and yDV with a log10 profile of the magnitudes of
   the entries in the ZV object. tausmall and tau big provide
   cutoffs within which to examine the entries. pnzero, pnsmall 
   and pnbig are addresses to hold the number of entries zero,
   smaller than tausmall and larger than taubig, respectively.

   created -- 97feb14, cca
   ------------------------------------------------------------------
*/
void
ZV_log10profile (
   ZV      *zv,
   int      npts,
   DV       *xDV,
   DV       *yDV,
   double   tausmall,
   double   taubig,
   int      *pnzero,
   int      *pnsmall,
   int      *pnbig
) {
double   deltaVal, imag, maxval, minval, real, val ;
double   *zvec, *sums, *x, *y ;
int      ii, ipt, jj, nbig, nsmall, nzero, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || npts <= 0 || xDV == NULL || yDV == NULL
   || tausmall < 0.0 || taubig < 0.0 || tausmall > taubig
   || pnzero == NULL || pnsmall == NULL || pnbig == NULL ) {
   fprintf(stderr, 
       "\n fatal error in ZV_log10profile(%p,%d,%p,%p,%f,%f,%p,%p,%p)"
       "\n bad input\n",
       zv, npts, xDV, yDV, tausmall, taubig, pnzero, pnsmall, pnbig) ;
   exit(-1) ;
}
/*
   -------------------------------------
   find the largest and smallest entries 
   in the range [tausmall, taubig]
   -------------------------------------
*/
nbig = nsmall = nzero = 0 ;
minval = maxval = 0.0 ;
ZV_sizeAndEntries(zv, &size, &zvec) ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   real = fabs(zvec[jj]) ;
   imag = fabs(zvec[jj+1]) ;
   if ( real == 0.0 ) {
      val =  imag ;
   } else if ( imag == 0.0 ) {
      val =  real ;
   } else if ( real >= imag ) {
      val = imag/real ;
      val = real*sqrt(1 + val*val) ;
   } else {
      val = real/imag ;
      val = imag*sqrt(1 + val*val) ;
   }
   if ( val == 0.0 ) {
      nzero++ ;
   } else if ( val <= tausmall ) {
      nsmall++ ;
   } else if ( val >= taubig ) {
      nbig++ ;
   } else {
      if ( minval == 0.0 || minval > val ) {
         minval = val ;
      }
      if ( maxval < val ) {
         maxval = val ;
      }
   }
}
*pnzero  = nzero  ;
*pnsmall = nsmall ;
*pnbig   = nbig   ;
#if MYDEBUG > 0
fprintf(stdout, 
        "\n nzero = %d, minval = %e, nsmall = %d, maxval = %e, nbig = %d",
        nzero, minval, nsmall, maxval, nbig) ;
#endif
/*
   ------------------
   set up the buckets
   ------------------
*/
DV_setSize(xDV, npts) ;
DV_setSize(yDV, npts) ;
x = DV_entries(xDV) ;
y = DV_entries(yDV) ;
sums = DVinit(npts, 0.0) ;
minval = log10(minval) ;
maxval = log10(maxval) ;
/*
minval = log10(tausmall) ;
maxval = log10(taubig) ;
*/
deltaVal = (maxval - minval)/(npts - 1) ;
DVfill(npts, x, 0.0) ;
DVfill(npts, y, 0.0) ;
/*
   --------------------------------
   fill the sums and counts vectors
   --------------------------------
*/
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   real = fabs(zvec[jj]) ;
   imag = fabs(zvec[jj+1]) ;
   if ( real == 0.0 ) {
      val =  imag ;
   } else if ( imag == 0.0 ) {
      val =  real ;
   } else if ( real >= imag ) {
      val = imag/real ;
      val = real*sqrt(1 + val*val) ;
   } else {
      val = real/imag ;
      val = imag*sqrt(1 + val*val) ;
   }
   if ( tausmall < val && val < taubig ) {
      ipt = (log10(val) - minval) / deltaVal ;
      sums[ipt] += val ;
      y[ipt]++ ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n sum(y) = %.0f", DV_sum(yDV)) ;
#endif
/*
   ---------------------------
   set the x-coordinate vector
   ---------------------------
*/
for ( ipt = 0 ; ipt < npts ; ipt++ ) {
   if ( sums[ipt] == 0.0 ) {
      x[ipt] = minval + ipt*deltaVal ;
   } else {
      x[ipt] = log10(sums[ipt]/y[ipt]) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DVfree(sums) ;

return ; }

/*--------------------------------------------------------------------*/
