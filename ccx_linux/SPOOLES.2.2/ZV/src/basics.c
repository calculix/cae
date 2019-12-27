/*  basics.C  */

#include "../ZV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor method

   created -- 98jan22, cca
   -----------------------
*/
ZV *
ZV_new ( 
   void 
) {
ZV   *zv ;

ALLOCATE(zv, struct _ZV, 1) ;

ZV_setDefaultFields(zv) ;

return(zv) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98jan22, cca
   -----------------------
*/
void
ZV_setDefaultFields ( 
   ZV   *zv
) {
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_setDefaultFields(%p)"
           "\n bad input\n", zv) ;
   exit(-1) ;
}
zv->size    =   0  ;
zv->maxsize =   0  ;
zv->owned   =   0  ;
zv->vec     = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 98jan22, cca
   -----------------------
*/
void
ZV_clearData ( 
   ZV   *zv
) {
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_clearData(%p)"
           "\n bad input\n", zv) ;
   exit(-1) ;
}
if ( zv->vec != NULL && zv->owned == 1 ) {
   DVfree(zv->vec) ;
}
ZV_setDefaultFields(zv) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 98jan22, cca
   -----------------------
*/
void
ZV_free ( 
   ZV   *zv
) {
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_free(%p)"
           "\n bad input\n", zv) ;
   exit(-1) ;
}
ZV_clearData(zv) ;
FREE(zv) ;

return ; }

/*--------------------------------------------------------------------*/
