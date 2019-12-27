/*  basics.c  */

#include "../BKL.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor

   created -- 95oct07, cca
   -----------------------
*/
BKL *
BKL_new (
   void
) {
BKL   *bkl ;

ALLOCATE(bkl, struct _BKL, 1) ;

BKL_setDefaultFields(bkl) ;

return(bkl) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 95oct07, cca
   -----------------------
*/
void
BKL_setDefaultFields (
   BKL   *bkl
) {
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_setDefaultFields(%p)"
           "\n bad input\n", bkl) ;
   exit(-1) ;
}
bkl->bpg       = NULL ;
bkl->ndom      =   0  ;
bkl->nseg      =   0  ;
bkl->nreg      =   0  ;
bkl->totweight =   0  ;
bkl->npass     =   0  ;
bkl->npatch    =   0  ;
bkl->nflips    =   0  ;
bkl->nimprove  =   0  ;
bkl->ngaineval =   0  ;
bkl->colors    = NULL ;
bkl->alpha     =  0.0 ;
bkl->cweights[0] = bkl->cweights[1] = bkl->cweights[2] = 0 ;
bkl->regwghts  = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   clear the data fields

   created  -- 95oct07, cca
   modified -- 95dec07, cca
      memory leak (bkl->regwghts) fixed
   ------------------------------------
*/
void
BKL_clearData (
   BKL   *bkl
) {
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_clearData(%p)"
           "\n bad input\n", bkl) ;
   exit(-1) ;
}
if ( bkl->colors != NULL ) {
   IVfree(bkl->colors) ;
}
if ( bkl->bpg != NULL
   && bkl->bpg->graph != NULL
   && bkl->bpg->graph->vwghts == NULL 
   && bkl->regwghts != NULL ) {
   IVfree(bkl->regwghts) ;
}
BKL_setDefaultFields(bkl) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 95oct07, cca
   -----------------------
*/
void
BKL_free (
   BKL   *bkl
) {
if ( bkl == NULL ) {
   fprintf(stderr, "\n fatal error in BKL_free(%p)"
           "\n bad input\n", bkl) ;
   exit(-1) ;
}
BKL_clearData(bkl) ;
FREE(bkl) ;

return ; }

/*--------------------------------------------------------------------*/
