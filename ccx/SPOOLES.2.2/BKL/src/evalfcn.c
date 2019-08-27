/*  evalfcn.c  */

#include "../BKL.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   evaluate the partition

   created -- 95oct07, cca
   -----------------------
*/
float
BKL_evalfcn (
   BKL   *bkl
) {
float   cost ;
int     wmax, wmin ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_evalfcn(%p)"
           "\n bad input\n",  bkl) ;
   exit(-1) ;
}
if ( bkl->cweights[1] <= bkl->cweights[2] ) {
   wmin = bkl->cweights[1] ;
   wmax = bkl->cweights[2] ;
} else {
   wmin = bkl->cweights[2] ;
   wmax = bkl->cweights[1] ;
}
if ( wmin == 0 ) {
   cost = ((float) bkl->totweight) * bkl->totweight ;
} else {
   cost = bkl->cweights[0] * (1. + (bkl->alpha * wmax)/wmin) ;
}
return(cost) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   evaluate the partition

   created -- 95oct07, cca
   -----------------------
*/
float
BKL_eval (
   BKL   *bkl,
   int   Sweight,
   int   Bweight,
   int   Wweight
) {
float   cost ;
int     wmax, wmin ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_evalfcn(%p)"
           "\n bad input\n",  bkl) ;
   exit(-1) ;
}
if ( Bweight <= Wweight ) {
   wmin = Bweight ;
   wmax = Wweight ;
} else {
   wmin = Wweight ;
   wmax = Bweight ;
}
if ( wmin == 0 ) {
   cost = ((float) bkl->totweight) * bkl->totweight ;
} else {
   cost = Sweight * (1. + (bkl->alpha * wmax)/wmin) ;
}
return(cost) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   evaluate the (deltaS, deltaB and deltaW) of a domain flip

   created -- 950ct11, cca
   ---------------------------------------------------------
*/
void
BKL_evalgain (
   BKL   *bkl,
   int   dom,
   int   *pdeltaS,
   int   *pdeltaB,
   int   *pdeltaW
) {
int   ii, newc, oldc, seg, size ;
int   *adj, *colors, *regwghts ;
int   stats[3] ;
/*
   ---------------
   check the input
   ---------------
*/
if (  bkl == NULL || dom < 0 || dom >= bkl->ndom
   || pdeltaS == NULL || pdeltaB == NULL || pdeltaW == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_evalGain(%p,%d,%p,%p,%p)"
           "\n bad input\n", bkl, dom, pdeltaS, pdeltaB, pdeltaW) ;
   exit(-1) ;
}
colors   = bkl->colors   ;
regwghts = bkl->regwghts ;
stats[0] = stats[1] = stats[2] = 0 ;
/*
   ---------------
   flip the domain
   ---------------
*/
if ( colors[dom] == 1 ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n domain %d, old color = 1, new color = 2", dom) ;
   fflush(stdout) ;
#endif
   stats[1] -= regwghts[dom] ;
   stats[2] += regwghts[dom] ;
   colors[dom] = 2 ;
} else {
#if MYDEBUG > 0
   fprintf(stdout, "\n domain %d, old color = 2, new color = 1", dom) ;
   fflush(stdout) ;
#endif
   stats[2] -= regwghts[dom] ;
   stats[1] += regwghts[dom] ;
   colors[dom] = 1 ;
}
/*
   -------------------------------
   loop over the adjacent segments
   -------------------------------
*/
Graph_adjAndSize(bkl->bpg->graph, dom, &size, &adj) ;
for ( ii = 0 ; ii < size ; ii++ ) {
   seg = adj[ii] ;
   oldc = colors[seg] ;
   newc = BKL_segColor(bkl, seg) ;
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n segment %d, weight = %d, old color = %d, new color = %d",
           seg, regwghts[seg], oldc, newc) ;
   fflush(stdout) ;
#endif
   if ( oldc != newc ) {
      stats[oldc] -= regwghts[seg] ;
      stats[newc] += regwghts[seg] ;
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n stats = < %d %d %d >",
           stats[0], stats[1], stats[2]) ;
   fflush(stdout) ;
#endif
   }
}
#if MYDEBUG > 0
   fprintf(stdout, "\n stats = < %d %d %d > ",
           stats[0], stats[1], stats[2]) ;
   fflush(stdout) ;
#endif
/*
   ------------------------
   set the output variables
   ------------------------
*/
*pdeltaS = stats[0] ;
*pdeltaB = stats[1] ;
*pdeltaW = stats[2] ;
/*
   --------------------
   flip the domain back
   --------------------
*/
if ( colors[dom] == 1 ) {
   colors[dom] = 2 ;
} else {
   colors[dom] = 1 ;
}
/*
   ---------------------------------
   increment the number of gainevals
   ---------------------------------
*/
bkl->ngaineval++ ;

return ; }

/*--------------------------------------------------------------------*/
