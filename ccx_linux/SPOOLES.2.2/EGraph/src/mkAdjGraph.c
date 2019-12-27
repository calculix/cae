/*  mkAdjGraph.c.c  */

#include "../EGraph.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   create a Graph object that holds the adjacency graph 
   of the assembled elements. 

   created -- 95nov03, cca
   ----------------------------------------------------
*/
Graph *
EGraph_mkAdjGraph ( 
   EGraph   *egraph 
) {
int     elem, esize, i, nelem, nvtx, v, vsize, w ;
int     *eind, *head, *link, *marker, *offsets, *vind ;
IVL     *eadjIVL, *gadjIVL ;
Graph   *graph ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || (eadjIVL = egraph->adjIVL) == NULL ) {
   fprintf(stderr, "\n fatal error in EGraph_mkAdjGraph(%p)"
           "\n bad input\n", egraph) ;
   exit(-1) ;
}
nelem  = egraph->nelem  ;
nvtx   = egraph->nvtx   ;
/*
   --------------------------------
   set up the linked list structure
   --------------------------------
*/
head    = IVinit(nvtx,  -1) ;
link    = IVinit(nelem, -1) ;
offsets = IVinit(nelem,  0) ;
/*
   -----------------------------------------------------------
   sort the vertices in each element list into ascending order
   and link them into their first vertex
   -----------------------------------------------------------
*/
for ( elem = 0 ; elem < nelem ; elem++ ) {
   IVL_listAndSize(eadjIVL, elem, &esize, &eind) ;
   if ( esize > 0 ) {
      IVqsortUp(esize, eind) ;
      v          = eind[0] ;
      link[elem] = head[v] ;
      head[v]    = elem    ;
   }
}
/*
   ---------------------------
   create the new Graph object
   ---------------------------
*/
graph = Graph_new() ;
Graph_init1(graph, egraph->type, nvtx, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
gadjIVL = graph->adjIVL ;
/*
   ----------------------
   loop over the vertices
   ----------------------
*/
vind   = IVinit(nvtx, -1) ;
marker = IVinit(nvtx, -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
/*
   ---------------------------------
   loop over the supporting elements
   ---------------------------------
*/
   vsize = 0 ;
   vind[vsize++] = v ;
   marker[v]     = v ;
   while ( (elem = head[v]) != -1 ) {
/*
fprintf(stdout, "\n    checking out element %d :", jelem) ;
*/
      head[v] = link[elem] ;
      IVL_listAndSize(eadjIVL, elem, &esize, &eind) ;
      for ( i = 0 ; i < esize ; i++ ) {
         w = eind[i] ;
         if ( marker[w] != v ) {
            marker[w]     = v ;
            vind[vsize++] = w ;
         }
      }
      if ( (i = ++offsets[elem]) < esize ) {
         w          = eind[i] ;
         link[elem] = head[w] ;
         head[w]    = elem    ;
      }
   }
   IVqsortUp(vsize, vind) ;
   IVL_setList(gadjIVL, v, vsize, vind) ;
}
graph->nedges = gadjIVL->tsize ;
if ( egraph->type == 0 ) {
   graph->totvwght = nvtx ;
} else if ( egraph->type == 1 ) {
/*
   ------------------------------
   fill the vertex weights vector
   ------------------------------
*/
   IVcopy(nvtx, graph->vwghts, egraph->vwghts) ;
   graph->totvwght = IVsum(nvtx, graph->vwghts) ;
}
graph->totewght = graph->nedges ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(head)    ;
IVfree(link)    ;
IVfree(marker)  ;
IVfree(vind)    ;
IVfree(offsets) ;

return(graph) ; }

/*--------------------------------------------------------------------*/
