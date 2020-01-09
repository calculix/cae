/*  basics.c  */

#include "../Graph.h"

#define MYTRACE 0
#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- create and return a new Graph object

   created -- 95sep27, cca
   -----------------------------------------------
*/
Graph *
Graph_new ( 
   void
) {
Graph   *g ;

#if MYTRACE > 0
fprintf(stdout, "\n just inside Graph_new()") ;
fflush(stdout) ;
#endif

ALLOCATE(g, struct _Graph, 1) ;

Graph_setDefaultFields(g) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Graph_new()") ;
fflush(stdout) ;
#endif

return(g) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- set the default fields for the Graph object

   created -- 95sep27, cca
   ------------------------------------------------------
*/
void
Graph_setDefaultFields (
   Graph   *g
) {

#if MYTRACE > 0
fprintf(stdout, "\n just inside Graph_setDefaultFields(%p)", g) ;
fflush(stdout) ;
#endif

if ( g == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_setDefaultFields(%p)"
           "\n graph is NULL\n", g) ;
   exit(-1) ;
}
g->type     =  0   ;
g->nvtx     =  0   ;
g->nvbnd    =  0   ;
g->nedges   =  0   ;
g->totvwght =  0   ;
g->totewght =  0   ;
g->adjIVL   = NULL ;
g->vwghts   = NULL ;
g->ewghtIVL = NULL ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Graph_setDefaultFields(%p)", g) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 95sep27, cca
   --------------------------------
*/
void
Graph_clearData (
   Graph   *g
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside Graph_clearData(%p)", g) ;
fflush(stdout) ;
#endif

if ( g == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_clearData(%p)"
           "\n graph is NULL\n", g) ;
   exit(-1) ;
}

if ( g->adjIVL != NULL ) {
   IVL_free(g->adjIVL) ;
}
if ( g->vwghts != NULL ) {
   IVfree(g->vwghts) ;
}
if ( g->ewghtIVL != NULL ) {
   IVL_free(g->ewghtIVL) ;
}
Graph_setDefaultFields(g) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Graph_clearData(%p)", g) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- free the Graph object

   created -- 95sep27, cca
   --------------------------------
*/
void
Graph_free (
   Graph   *g
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside Graph_free(%p)", g) ;
fflush(stdout) ;
#endif

if ( g == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_free(%p)"
           "\n graph is NULL\n", g) ;
   exit(-1) ;
}

Graph_clearData(g) ;

FREE(g) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Graph_free(%p)", g) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
