/*  init.c  */

#include "../Graph.h"

#define MYTRACE 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   basic initializer for the Graph object

   type      -- graph type
   nvtx      -- # of vertices
   nvbnd     -- # of boundary vertices
   nedges    -- # of edges
   adjType   -- IVL type for adjacency object
   ewghtType -- IVL type for edge weights object

   created -- 95sep27, cca
   ---------------------------------------------
*/
void
Graph_init1 (
   Graph   *g,
   int     type,
   int     nvtx,
   int     nvbnd,
   int     nedges,
   int     adjType,
   int     ewghtType
) {
int   nvtot ;
#if MYTRACE > 0
fprintf(stdout, "\n just inside Graph_init1(%p,%d,%d,%d,%d,%d,%d)",
        g, type, nvtx, nvbnd, nedges, adjType, ewghtType) ;
fflush(stdout) ;
#endif
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_init1(%p,%d,%d,%d,%d,%d,%d)"
           "\n graph is NULL\n",
           g, type, nvtx, nvbnd, nedges, adjType, ewghtType) ;
   exit(-1) ;
}
if ( type < 0 || type >= 4 ) {
   fprintf(stderr, "\n fatal error in Graph_init1(%p,%d,%d,%d,%d,%d,%d)"
           "\n invalid type = %d, must be in [0,3]\n",
           g, type, nvtx, nvbnd, nedges, adjType, ewghtType, type) ;
   exit(-1) ;
}
if ( nvtx <= 0 ) {
   fprintf(stderr, "\n fatal error in Graph_init1(%p,%d,%d,%d,%d,%d,%d)"
           "\n nvtx = %d, must be positive\n",
           g, type, nvtx, nvbnd, nedges, adjType, ewghtType, nvtx) ;
   exit(-1) ;
}
if ( nvbnd < 0 ) {
   fprintf(stderr, "\n fatal error in Graph_init1(%p,%d,%d,%d,%d,%d,%d)"
           "\n nvbnd = %d, must be nonnegative\n",
           g, type, nvtx, nvbnd, nedges, adjType, ewghtType, nvbnd) ;
   exit(-1) ;
}
if ( nedges < 0 ) {
   fprintf(stderr, "\n fatal error in Graph_init1(%p,%d,%d,%d,%d,%d,%d)"
           "\n nedges = %d, must be nonnegative\n",
           g, type, nvtx, nvbnd, nedges, adjType, ewghtType, nedges) ;
   exit(-1) ;
}
#if MYTRACE > 0
fprintf(stdout, "\n input checks out") ;
fflush(stdout) ;
#endif
switch ( adjType ) {
case IVL_CHUNKED :
case IVL_SOLO    :
case IVL_UNKNOWN :
   break ;
default :
   fprintf(stderr, "\n fatal error in Graph_init1(%p,%d,%d,%d,%d,%d,%d)"
           "\n invalid adjType = %d\n",
           g, type, nvtx, nvbnd, nedges, adjType, ewghtType, adjType) ;
   exit(-1) ;
}
switch ( ewghtType ) {
case IVL_CHUNKED :
case IVL_SOLO    :
case IVL_UNKNOWN :
   break ;
default :
   fprintf(stderr, "\n fatal error in Graph_init1(%p,%d,%d,%d,%d,%d,%d)"
           "\n invalid ewghtType = %d\n",
        g, type, nvtx, nvbnd, nedges, adjType, ewghtType, ewghtType) ;
   exit(-1) ;
}
/*
   -------------------------
   clear the data structures
   -------------------------
*/
Graph_clearData(g) ;
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
g->type   = type   ;
g->nvtx   = nvtx   ;
g->nvbnd  = nvbnd  ;
g->nedges = nedges ;
nvtot     = nvtx + nvbnd ;
/*
   ----------------------------
   set the adjacency IVL object
   ----------------------------
*/
g->adjIVL = IVL_new() ;
if ( nedges == 0 || adjType == IVL_UNKNOWN ) {
   IVL_init1(g->adjIVL, adjType, nvtot) ;
} else {
   IVL_init2(g->adjIVL, adjType, nvtot, nedges) ;
}
if ( type % 2 == 1 ) {
/*
   -----------------------------
   set the vertex weights vector
   -----------------------------
*/
   g->vwghts = IVinit(nvtot, 0) ;
}
if ( type >= 2 ) {
/*
   -------------------------------
   set the edge weights IVL object
   -------------------------------
*/
   g->ewghtIVL = IVL_new() ;
   if ( nedges == 0 || ewghtType == IVL_UNKNOWN ) {
      IVL_init1(g->ewghtIVL, ewghtType, nvtot) ;
   } else {
      IVL_init2(g->ewghtIVL, ewghtType, nvtot, nedges) ;
   }
}

#if MYTRACE > 0
fprintf(stdout, "\n leaving Graph_init1(%p,%d,%d,%d,%d,%d,%d)",
        g, type, nvtx, nvbnd, nedges, adjType, ewghtType) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   second initializer for the Graph object.
   this function is used in the I/O routines
   Graph_readFromFormattedFile(Graph *g, FILE *fp) and
   Graph_readFromBinaryFile(Graph *g, FILE *fp) where
   the IVL object(s) and vwghts[] vector are created
   independently.

   type     -- graph type
   nvtx     -- # of vertices
   nvbnd    -- # of boundary vertices
   nedges   -- # of edges
   totvwght -- total vertex weight
   totewght -- total edge weight
   adjIVL   -- IVL object for adjacency structure
   vwghts   -- pointer to integer vector for vertex weights
   ewghtIVL -- IVL object for edge weights 

   created -- 95sep27, cca
   --------------------------------------------------------
*/
void
Graph_init2 (
   Graph   *g,
   int     type,
   int     nvtx,
   int     nvbnd,
   int     nedges,
   int     totvwght,
   int     totewght,
   IVL     *adjIVL,
   int     *vwghts,
   IVL     *ewghtIVL
) {
#if MYTRACE > 0
fprintf(stdout, 
        "\n just inside Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght, 
        adjIVL, vwghts, ewghtIVL) ;
   fflush(stdout) ;
#endif
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n graph is NULL\n",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL) ;
   exit(-1) ;
}
if ( type < 0 || type >= 4 ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n invalid type = %d, must be in [0,3]\n",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL, type) ;
   exit(-1) ;
}
if ( nvtx <= 0 ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n nvtx = %d, must be positive\n",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL, nvtx) ;
   exit(-1) ;
}
if ( nvbnd < 0 ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n nvbnd = %d, must be nonnegative\n",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL, nvbnd) ;
   exit(-1) ;
}
if ( nedges < 0 ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n nedges = %d, must be nonnegative\n",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL, nedges) ;
   exit(-1) ;
}
if ( adjIVL == NULL ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n adjIVL is NULL\n",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL) ;
   exit(-1) ;
}
if ( (type % 2 == 1) && (vwghts == NULL) ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n type = %d, vwghts is NULL",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL, type) ;
   exit(-1) ;
}
if ( (type >= 2) && (ewghtIVL == NULL) ) {
   fprintf(stdout, 
        "\n fatal error in Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n type = %d, ewghtIVL is NULL",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL, type) ;
   exit(-1) ;
}
/*
   -------------------------
   clear the data structures
   -------------------------
*/
Graph_clearData(g) ;
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
g->type     = type     ;
g->nvtx     = nvtx     ;
g->nvbnd    = nvbnd    ;
g->nedges   = nedges   ;
g->totvwght = totvwght ;
g->totewght = totewght ;
/*
   ---------------------------------------------
   set the IVL objects and vertex weights vector
   ---------------------------------------------
*/
g->adjIVL = adjIVL ;
if ( type % 2 == 1 ) {
   g->vwghts = vwghts ;
}
if ( type >= 2 ) {
   g->ewghtIVL = ewghtIVL ;
}

#if MYTRACE > 0
fprintf(stdout, 
        "\n leaving Graph_init2(%p,%d,%d,%d,%d,%d,%d,%p,%p,%p)"
        "\n type = %d, ewghtIVL is NULL",
        g, type, nvtx, nvbnd, nedges, totvwght, totewght,
        adjIVL, vwghts, ewghtIVL, type) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
