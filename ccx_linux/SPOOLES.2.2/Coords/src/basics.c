/*  basics.c  */

#include "../Coords.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 95dec17, cca
   -----------------------
*/
Coords *
Coords_new ( 
   void 
) {
Coords   *coords ;

ALLOCATE(coords, struct _Coords, 1) ;
Coords_setDefaultFields(coords) ;

return(coords) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 95dec17, cca
   -----------------------
*/
void
Coords_setDefaultFields (
   Coords   *coords
) {
if ( coords == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_setDefaultFields(%p)"
           "\n bad input", coords) ;
   exit(-1) ;
}
coords->type   = COORDS_BY_TUPLE ;
coords->ndim   =   0  ;
coords->ncoor  =   0  ;
coords->coors  = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 95dec17, cca
   --------------------------------------------------
*/
void
Coords_clearData ( 
   Coords   *coords 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_clearData(%p)"
           "\n bad input\n", coords) ;
   exit(-1) ;
}
if ( coords->coors != NULL ) {
   FVfree(coords->coors) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
Coords_setDefaultFields(coords) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 95dec17, cca
   ------------------------------------------
*/
Coords *
Coords_free ( 
   Coords   *coords 
) {
if ( coords == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_free(%p)"
           "\n bad input\n", coords) ;
   exit(-1) ;
}
Coords_clearData(coords) ;
FREE(coords) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
