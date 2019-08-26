/*  profile  */
 
#include "../InpMtx.h"
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   to fill xDV and yDV with a log10 profile of the magnitudes of
   the entries in the InpMtx object. tausmall and tau big provide
   cutoffs within which to examine the entries. pnsmall and pnbig 
   are address to hold the number of entries smaller than tausmall,

   and larger than taubig, respectively.
 
   created -- 97feb14, cca
   ------------------------------------------------------------------
*/
void
InpMtx_log10profile (
   InpMtx    *inpmtx,
   int        npts,
   DV         *xDV,
   DV         *yDV,
   double     tausmall,
   double     taubig,
   int        *pnzero,
   int        *pnsmall,
   int        *pnbig
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL 
   || npts <= 0 || xDV == NULL || yDV == NULL 
   || tausmall < 0.0 || taubig < 0.0 || tausmall > taubig
   || pnzero == NULL || pnsmall == NULL || pnbig == NULL ) {
   fprintf(stderr, 
   "\n fatal error in InpMtx_log10profile(%p,%d,%p,%p,%f,%f,%p,%p,%p)"
      "\n bad input\n",
       inpmtx, npts, xDV, yDV, tausmall, taubig, 
       pnzero, pnsmall, pnbig) ;
   exit(-1) ;
}
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   DV   *dv = DV_new();
   DV_init(dv, inpmtx->nent, InpMtx_dvec(inpmtx)) ;
   DV_log10profile(dv, npts, xDV, yDV, 
                   tausmall, taubig, pnzero, pnsmall, pnbig) ;
   DV_free(dv) ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   ZV   *zv = ZV_new();
   ZV_init(zv, inpmtx->nent, InpMtx_dvec(inpmtx)) ;
   ZV_log10profile(zv, npts, xDV, yDV, 
                   tausmall, taubig, pnzero, pnsmall, pnbig) ;
   ZV_free(zv) ;
}
return ; }

/*--------------------------------------------------------------------*/
