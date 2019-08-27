/*  init.c  */

#include "../GPart.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   initialize the GPart object

   created -- 95oct05, cca
   ---------------------------
*/
void
GPart_init (
   GPart   *gpart,
   Graph   *g
) {
if ( gpart == NULL || g == NULL || g->nvtx <= 0 ) {
   fprintf(stderr, "\n fatal error in GPart_init(%p,%p)"
           "\n bad input\n", gpart, g) ;
   exit(-1) ;
}
GPart_clearData(gpart) ;
gpart->nvtx     = g->nvtx  ;
gpart->nvbnd    = g->nvbnd ;
gpart->g        = g ;
gpart->ncomp    = 1 ;
IV_setSize(&gpart->compidsIV, g->nvtx) ;
IV_fill(&gpart->compidsIV, 1) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the message fields

   created -- 96oct21, cca
   -----------------------
*/
void 
GPart_setMessageInfo (
   GPart   *gpart,
   int     msglvl,
   FILE    *msgFile
) {
if ( gpart == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_setMessageInfo(%p,%d,%p)"
           "\n bad input\n", gpart, msglvl, msgFile) ;
   exit(-1) ;
}
gpart->msglvl = msglvl ;
if ( msgFile != NULL ) {
   gpart->msgFile = msgFile ;
} else {
   gpart->msgFile = stdout ;
}
return ; }

/*--------------------------------------------------------------------*/
