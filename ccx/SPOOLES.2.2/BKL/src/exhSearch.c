/*  exhSearch.c  */

#include "../BKL.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   perform an exhaustive search over a subspace of domains
   
   mdom    -- number of domains in the subspace
   domids  -- vector of domain ids, size mdom
   tcolors -- temporary vector to hold active domain colors, size mdom

   note : region colors and component weights of the best
          partition are left in bkl->colors[] and bkl->cweights[].

   return value -- cost of best partition

   created -- 95oct07, cca
   --------------------------------------------------------------------
*/
float
BKL_exhSearch (
   BKL   *bkl,
   int   mdom,
   int   domids[],
   int   tcolors[]
) {
float   bestcost, newcost ;
int     idom, ierr, iflip, iloc, jloc, nflip ;
int     *colors ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL || mdom < 1 || domids == NULL || tcolors == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_exhaustiveSearch(%p,%d,%p,%p)"
           "\n bad input\n",  bkl, mdom, domids, tcolors) ;
   exit(-1) ;
}
colors = bkl->colors ;
bkl->nflips = 0 ;
#if MYDEBUG > 0
fprintf(stdout, "\n inside BKL_exhSearch(%p,%d,%p,%p)",
        bkl, mdom, domids, tcolors) ;
fprintf(stdout, "\n bkl->nflips = %d", bkl->nflips) ;
fflush(stdout) ;
#endif
/*
   ---------------------------------------------
   copy the present colors and component weights
   ---------------------------------------------
*/
for ( iloc = 0 ; iloc < mdom ; iloc++ ) {
   idom = domids[iloc] ;
   tcolors[iloc] = colors[idom] ;
}
/*
   ---------------------
   compute the best cost
   ---------------------
*/
bestcost = BKL_evalfcn(bkl) ;
#if MYDEBUG > 0
fprintf(stdout, "\n inside BKL_exhSearch(%p,%d,%p,%p)",
        bkl, mdom, domids, tcolors) ;
fprintf(stdout, "\n %d domain ids : ", mdom) ;
IVfp80(stdout, mdom, domids, 20, &ierr) ;
fprintf(stdout, "\n color weights < %6d %6d %6d >, cost %9.2f",
        bkl->cweights[0], bkl->cweights[1], bkl->cweights[2],
        bestcost) ;
fflush(stdout) ;
#endif
#if MYDEBUG > 2
fprintf(stdout, "\n colors ") ;
IVfp80(stdout, bkl->nreg, colors, 80, &ierr) ;
fflush(stdout) ;
#endif
/*
   ---------------------------------
   count the number of flips to make
   ---------------------------------
*/
for ( idom = 0, nflip = 1 ; idom < mdom ; idom++ ) {
   nflip *= 2 ;
}
/*
   --------------------------
   loop over the 2^mdom flips
   --------------------------
*/
for ( iflip = 1 ; iflip < nflip ; iflip++ ) {
   iloc = BKL_greyCodeDomain(bkl, iflip) ;
   idom = domids[iloc] ;
#if MYDEBUG > 1
   fprintf(stdout, "\n FLIP %4d domain %4d", bkl->nflips, idom) ;
   fflush(stdout) ;
#endif
#if MYDEBUG > 2
   fprintf(stdout, "\n colors before flip") ;
   IVfp80(stdout, bkl->nreg, colors, 80, &ierr) ;
   fflush(stdout) ;
#endif
   BKL_flipDomain(bkl, idom) ;
#if MYDEBUG > 2
fprintf(stdout, "\n colors after flip") ;
IVfp80(stdout, bkl->nreg, colors, 80, &ierr) ;
fprintf(stdout, "\n cweights : < %9d %9d %9d > ",
        bkl->cweights[0], bkl->cweights[1], bkl->cweights[2]) ;
fflush(stdout) ;
#endif
   newcost = BKL_evalfcn(bkl) ;
#if MYDEBUG > 1
fprintf(stdout, ", < %6d %6d %6d >, cost %9.2f",
        bkl->cweights[0], bkl->cweights[1], bkl->cweights[2],
           newcost) ;
   fflush(stdout) ;
#endif
   if ( newcost < bestcost ) {
#if MYDEBUG > 1
      fprintf(stdout, ", better") ;
      fflush(stdout) ;
#endif
      bkl->nimprove++ ;
      for ( jloc = 0 ; jloc < mdom ; jloc++ ) {
         tcolors[jloc] = colors[domids[jloc]] ;
      }
      bestcost = newcost ;
   }
}
/*
   -----------------------------------------------------
   restore the best colors and update the segment colors
   -----------------------------------------------------
*/
for ( iloc = 0 ; iloc < mdom ; iloc++ ) {
   idom = domids[iloc] ;
   if ( colors[idom] != tcolors[iloc] ) {
      BKL_flipDomain(bkl, idom) ;
   }
}

return(bestcost) ; }
  
/*--------------------------------------------------------------------*/
