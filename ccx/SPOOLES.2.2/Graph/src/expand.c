/*  expand.c  */

#include "../Graph.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- take a Graph object and a map to expand it, create and
              return a bigger unit weight Graph object. this is useful 
              for expanding a compressed graph into a unit weight graph.

   created -- 96mar02, cca
   ---------------------------------------------------------------------
*/
Graph *
Graph_expand2 (
   Graph   *g,
   IV      *mapIV
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  g == NULL 
   || g->nvtx <= 0
   || mapIV == NULL 
   || IV_entries(mapIV) == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_expand2(%p,%p)"
           "\n bad input\n", g, mapIV) ;
   exit(-1) ;
}
return(Graph_expand(g, IV_size(mapIV), IV_entries(mapIV))) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- take a Graph object and a map to expand it, create and
              return a bigger unit weight Graph object. this is useful 
              for expanding a compressed graph into a unit weight graph.

   created -- 96mar02, cca
   ---------------------------------------------------------------------
*/
Graph *
Graph_expand (
   Graph   *g,
   int     nvtxbig,
   int     map[]
) {
Graph   *gbig ;
int     count, ii, nedge, nvtx, v, vbig, vsize, w ;
int     *head, *indices, *link, *mark, *vadj ;
IVL     *adjIVL, *adjbigIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL || nvtxbig <= 0 || map == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_expand(%p,%d,%p)"
           "\n bad input\n", g, nvtxbig, map) ;
   exit(-1) ;
}
nvtx   = g->nvtx   ;
adjIVL = g->adjIVL ;
/*
   ----------------------------------------
   set up the linked lists for the vertices
   ----------------------------------------
*/
head = IVinit(nvtx,    -1) ;
link = IVinit(nvtxbig, -1) ;
for ( vbig = 0 ; vbig < nvtxbig ; vbig++ ) {
   v          = map[vbig] ;
   link[vbig] = head[v]   ;
   head[v]    = vbig      ;
}
/*
   --------------------------------
   create the expanded Graph object
   --------------------------------
*/
gbig = Graph_new() ;
Graph_init1(gbig, 0, nvtxbig, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
adjbigIVL = gbig->adjIVL ;
/*
   -------------------------------------------
   fill the lists in the expanded Graph object
   -------------------------------------------
*/
indices = IVinit(nvtxbig, -1) ;
mark    = IVinit(nvtx,    -1) ;
nedge   = 0 ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( head[v] != -1 ) {
/*
      ------------------------------
      load the indices that map to v
      ------------------------------
*/
      mark[v] = v ;
      count   = 0 ;
      for ( vbig = head[v] ; vbig != -1 ; vbig = link[vbig] ) {
         indices[count++] = vbig ;
      }
/*
      ---------------------------------------------------
      load the indices that map to vertices adjacent to v
      ---------------------------------------------------
*/
      IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w < nvtx && mark[w] != v ) {
            mark[w] = v ;
            for ( vbig = head[w] ; vbig != -1 ; vbig = link[vbig] ) {
               indices[count++] = vbig ;
            }
         }
      }
/*
      --------------------------------------
      sort the index list in ascending order
      --------------------------------------
*/
      IVqsortUp(count, indices) ;
/*
      -------------------------------------------------------
      each vertex in the big IVL object has its own list.
      -------------------------------------------------------
*/
      for ( vbig = head[v] ; vbig != -1 ; vbig = link[vbig] ) {
         IVL_setList(adjbigIVL, vbig, count, indices) ;
         nedge += count ;
      }
   }
}
gbig->nedges = nedge ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(head)    ;
IVfree(link)    ;
IVfree(indices) ;
IVfree(mark)    ;

return(gbig) ; }

/*--------------------------------------------------------------------*/
