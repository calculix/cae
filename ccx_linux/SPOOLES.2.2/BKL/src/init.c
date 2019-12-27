/*  init.c  */

#include "../BKL.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   initialize the object

   created -- 95oct07, cca
   -----------------------
*/
void
BKL_init (
   BKL     *bkl,
   BPG     *bpg,
   float   alpha
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( bkl == NULL || bpg == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_init(%p,%p,%f)"
           "\n bad input\n", bkl, bpg, alpha) ;
   exit(-1) ;
}
/*
   --------------
   clear the data
   --------------
*/
BKL_clearData(bkl) ;
/*
   ---------------------
   initialize the fields
   ---------------------
*/
bkl->bpg  = bpg ;
bkl->ndom = bpg->nX ;
bkl->nseg = bpg->nY ;
bkl->nreg = bpg->nX + bpg->nY ;
if ( bpg->graph->vwghts == NULL ) {
   bkl->totweight = bkl->nreg ;
   bkl->regwghts  = IVinit(bkl->nreg, 1) ;
} else {
   bkl->regwghts  = bpg->graph->vwghts ;
   bkl->totweight = IVsum(bkl->nreg, bkl->regwghts) ; 
}
bkl->colors = IVinit(bkl->nreg, 0) ;
bkl->alpha  = alpha ;

return ; }

/*--------------------------------------------------------------------*/
