/*  basics.c  */

#include "../GPart.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   construct a new instance of the GPart object

   created -- 95oct05, cca
   --------------------------------------------
*/
GPart *
GPart_new (
   void
) {
GPart   *gpart ;

ALLOCATE(gpart, struct _GPart, 1) ;

GPart_setDefaultFields(gpart) ;

return(gpart) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   set the default fields of the GPart object

   created  -- 95oct05, cca
   modified -- 95nov29, cca
      par, fch, sib and vtxMap fields included
   ---------------------------------------------
*/
void
GPart_setDefaultFields (
   GPart   *gpart
) {
if ( gpart == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_setDefaultFields(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
gpart->id       =  -1  ;
gpart->g        = NULL ;
gpart->nvtx     =   0  ;
gpart->nvbnd    =   0  ;
gpart->ncomp    =   0  ;
gpart->par      = NULL ;
gpart->fch      = NULL ;
gpart->sib      = NULL ;
IV_setDefaultFields(&gpart->compidsIV)  ;
IV_setDefaultFields(&gpart->cweightsIV) ;
IV_setDefaultFields(&gpart->vtxMapIV)   ;
gpart->msglvl  =   0  ;
gpart->msgFile = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   clear the data fields for a GPart object

   created  -- 95oct05, cca
   modified -- 95nov29, cca
      par, fch, sib and vtxMap fields included
   ---------------------------------------------
*/
void
GPart_clearData (
   GPart   *gpart
) {
if ( gpart == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_clearData(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
IV_clearData(&gpart->compidsIV)  ;
IV_clearData(&gpart->cweightsIV) ;
IV_clearData(&gpart->vtxMapIV)   ;
GPart_setDefaultFields(gpart) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------
   free the GPart object

   created  -- 95oct05, cca
   modified -- 95nov29, cca
      gpart now free'd
   ------------------------
*/
void
GPart_free (
   GPart   *gpart
) {
if ( gpart == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_free(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
GPart_clearData(gpart) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n trying to free gpart") ;
   fflush(stdout) ;
#endif
FREE(gpart) ;

return ; }

/*--------------------------------------------------------------------*/
