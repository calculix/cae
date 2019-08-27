/*  basics.C  */

#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor method

   created -- 95oct06, cca
   -----------------------
*/
DV *
DV_new ( 
   void 
) {
DV   *dv ;

ALLOCATE(dv, struct _DV, 1) ;

DV_setDefaultFields(dv) ;

return(dv) ; }
/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 95oct06, cca
   -----------------------
*/
void
DV_setDefaultFields ( 
   DV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_setDefaultFields(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
dv->size    =   0  ;
dv->maxsize =   0  ;
dv->owned   =   0  ;
dv->vec     = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 95oct06, cca
   -----------------------
*/
void
DV_clearData ( 
   DV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_clearData(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
if ( dv->vec != NULL && dv->owned == 1 ) {
   DVfree(dv->vec) ;
}
DV_setDefaultFields(dv) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 95oct06, cca
   -----------------------
*/
void
DV_free ( 
   DV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_free(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
DV_clearData(dv) ;
FREE(dv) ;

return ; }

/*--------------------------------------------------------------------*/
