/*  basics.c  */

#include "../EGraph.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor

   created -- 95nov03, cca
   -----------------------
*/
EGraph *
EGraph_new ( 
   void 
) {
EGraph   *eg ;

ALLOCATE(eg, struct _EGraph, 1) ;
EGraph_setDefaultFields(eg) ;

return(eg) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 95nov03, cca
   -----------------------
*/
void
EGraph_setDefaultFields ( 
   EGraph   *eg 
) {
if ( eg == NULL ) {
   fprintf(stderr, "\n fatal error in Egraph_setDefaultFields(%p)"
           "\n bad input\n", eg) ;
   exit(-1) ;
}
eg->type   =   0  ;
eg->nelem  =   0  ;
eg->nvtx   =   0  ;
eg->adjIVL = NULL ;
eg->vwghts = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 95nov03, cca
   -----------------------
*/
void
EGraph_clearData ( 
   EGraph   *eg 
) {
if ( eg == NULL ) {
   fprintf(stderr, "\n fatal error in Egraph_clearData(%p)"
           "\n bad input\n", eg) ;
   exit(-1) ;
}
if ( eg->adjIVL != NULL ) {
   IVL_free(eg->adjIVL) ;
}
if ( eg->vwghts != NULL ) {
   IVfree(eg->vwghts) ;
}
EGraph_setDefaultFields(eg) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 95nov03, cca
   -----------------------
*/
void
EGraph_free ( 
   EGraph   *eg 
) {
if ( eg == NULL ) {
   fprintf(stderr, "\n fatal error in Egraph_free(%p)"
           "\n bad input\n", eg) ;
   exit(-1) ;
}
EGraph_clearData(eg) ;
FREE(eg) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object
 
   created -- 95nov03, cca
   ----------------------------------------------
*/
int
EGraph_sizeOf ( 
   EGraph   *eg 
) {
int   bytes ;

bytes = sizeof(struct _EGraph) ;
if ( eg->adjIVL != NULL ) {
   bytes += IVL_sizeOf(eg->adjIVL) ;
}
if ( eg->vwghts != NULL ) {
   bytes += eg->nvtx * sizeof(int) ;
}
return(bytes) ; }

/*--------------------------------------------------------------------*/
