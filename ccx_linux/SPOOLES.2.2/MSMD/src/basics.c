/*  basics.c  */

#include "../MSMD.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor

   created -- 96feb25, cca
   -----------------------
*/
MSMD *
MSMD_new ( 
   void 
) {
MSMD   *msmd ;

ALLOCATE(msmd, struct _MSMD, 1) ;
MSMD_setDefaultFields(msmd) ;

return(msmd) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   set the default data fields
   
   created -- 96feb25, cca
   ---------------------------
*/
void
MSMD_setDefaultFields(
   MSMD   *msmd
) {
msmd->nvtx      =   0    ;
msmd->heap      = NULL   ;
msmd->incrIP    =   0    ;
msmd->baseIP    = NULL   ;
msmd->freeIP    = NULL   ;
msmd->vertices  = NULL   ;
IV_setDefaultFields(&msmd->ivtmpIV) ;
IV_setDefaultFields(&msmd->reachIV) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 96feb25, cca
   -----------------------
*/
void
MSMD_clearData ( 
   MSMD   *msmd
) {
IP        *ip ;
MSMDvtx   *first, *last, *v ;
/*
   --------------
   check the data
   --------------
*/
if ( msmd == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_clearData(%p)"
           "\n bad input\n", msmd) ;
   exit(-1) ;
}
if ( msmd->heap != NULL ) { 
#if MYDEBUG > 0
   fprintf(stdout, "\n trying to free the heap") ;
   fflush(stdout) ;
#endif
   IIheap_free(msmd->heap) ; 
}
if ( msmd->vertices != NULL ) { 
#if MYDEBUG > 0
   fprintf(stdout, "\n trying to free the vertices") ;
   fflush(stdout) ;
#endif
   first = msmd->vertices ;
   last  = first + msmd->nvtx - 1 ;
   for ( v = first ; v <= last ; v++ ) {
      if ( v->status == 'E' && v->adj != NULL ) {
         IVfree(v->adj) ;
      }
   }
   FREE(msmd->vertices) ;
}
IV_clearData(&msmd->ivtmpIV) ;
IV_clearData(&msmd->reachIV) ;
while ( (ip = msmd->baseIP) != NULL ) {
   msmd->baseIP = ip->next ;
   IP_free(ip) ;
}
MSMD_setDefaultFields(msmd) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 96feb25, cca
   -----------------------
*/
void
MSMD_free (
   MSMD   *msmd
) {
MSMD_clearData(msmd) ;
FREE(msmd) ;

return ; }

/*--------------------------------------------------------------------*/
