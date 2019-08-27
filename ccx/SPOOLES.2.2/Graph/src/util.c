/*  util.c  */

#include "../Graph.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   return the external degree (in terms of vertex weight) of vertex v

   created -- 95oct05, cca
   ---------------------------------------------------------------------
*/
int
Graph_externalDegree (
      Graph   *g,
      int     v
) {
int   ii, vsize, w, weight ;
int   *vadj, *vwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if (  g == NULL || v < 0 || g->nvtx + g->nvbnd <= v ) {
   fprintf(stderr, "\n fatal error in Graph_externalDegree(%p,%d)"
           "\n bad input\n", g, v) ;
   exit(-1) ;
}
vwghts = g->vwghts ;
Graph_adjAndSize(g, v, &vsize, &vadj) ;
for ( ii = 0, weight = 0 ; ii < vsize ; ii++ ) {
   if ( (w = vadj[ii]) != v ) {
      weight += (vwghts != NULL) ? vwghts[w] : 1 ;
   }
}
return(weight) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   method to access the adjacency list

   created -- 95oct05, cca
   -----------------------------------
*/
void
Graph_adjAndSize (
   Graph   *g,
   int     jvtx,
   int     *psize,
   int     **padj
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  g == NULL || jvtx < 0 || g->nvtx + g->nvbnd <= jvtx
   || psize == NULL || padj == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Graph_adjAndSize(%p,%d,%p,%p)"
           "\n bad input\n", g, jvtx, psize, padj) ;
   exit(-1) ;
}
if ( g->adjIVL == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_adjAndSize(%p,%d,%p,%p)"
           "\n g->adjIVL == NULL\n", g, jvtx, psize, padj) ;
   exit(-1) ;
}
IVL_listAndSize(g->adjIVL, jvtx, psize, padj) ;
/*
   here is fast code but not safe code
*/
/*
*psize = g->adjIVL->sizes[jvtx] ;
*padj  = g->adjIVL->p_vec[jvtx] ;
*/

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   method to access the adjacency list 
   and possibly the edge weight list for a vertex.

   created -- 95sep29, cca
   -----------------------------------------------
*/
void
Graph_adjAndEweights (
   Graph   *g,
   int     jvtx,
   int     *psize,
   int     **padj,
   int     **pewghts
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  g == NULL || jvtx < 0 || g->nvtx + g->nvbnd <= jvtx
   || psize == NULL || padj == NULL || pewghts == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Graph_adjAndEwghts(%p,%d,%p,%p,%p)"
           "\n bad input\n",
           g, jvtx, psize, padj, pewghts) ;
   exit(-1) ;
}
if ( g->adjIVL == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Graph_adjAndEwghts(%p,%d,%p,%p,%p)"
           "\n g->adjIVL == NULL\n",
           g, jvtx, psize, padj, pewghts) ;
   exit(-1) ;
}
if ( g->type >= 2 && g->ewghtIVL == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Graph_adjAndEwghts(%p,%d,%p,%p,%p)"
           "\n g->type = %d and g->ewghtIVL == NULL\n",
           g, jvtx, psize, padj, pewghts, g->type) ;
   exit(-1) ;
}
IVL_listAndSize(g->adjIVL, jvtx, psize, padj) ;
if ( g->type >= 2 ) {
   IVL_listAndSize(g->ewghtIVL, jvtx, psize, pewghts) ;
} else {
   *pewghts = NULL ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   return the number of bytes taken by the Graph object

   created -- 95oct05, cca
   ----------------------------------------------------
*/
int
Graph_sizeOf (
   Graph   *g
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_sizeOf(%p)"
           "\n bad input\n", g) ;
   exit(-1) ;
}
nbytes = sizeof(struct _Graph) ;
if ( g->vwghts != NULL ) {
   nbytes += (g->nvtx + g->nvbnd)*sizeof(int) ;
}
if ( g->adjIVL != NULL ) {
   nbytes += IVL_sizeOf(g->adjIVL) ;
}
if ( g->ewghtIVL != NULL ) {
   nbytes += IVL_sizeOf(g->ewghtIVL) ;
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   create and return an IV object filled 
   with a map from vertices to components

   created -- 96feb25, cca
   --------------------------------------
*/
IV *
Graph_componentMap (
   Graph   *g
) {
int   ii, last, ncomp, now, nvtx, u, usize, v, w ;
int   *list, *map, *uadj ;
IV    *mapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_componentMap(%p)"
           "\n bad input\n", g) ;
   exit(-1) ;
}
if ( (nvtx = g->nvtx) <= 0 ) {
   return(NULL) ;
}
mapIV = IV_new() ;
IV_init(mapIV, nvtx, NULL) ;
map = IV_entries(mapIV) ;

list = IVinit(nvtx, -1) ;
for ( v = 0, ncomp = 0 ; v < nvtx ; v++ ) {
   if ( map[v] == -1 ) {
/*
      -------------------------------
      seed vertex for a new component
      -------------------------------
*/
      map[v] = ncomp ;
      now = last = 0 ;
      list[0] = v ;
      while ( now <= last ) {
         u = list[now++] ;
         Graph_adjAndSize(g, u, &usize, &uadj) ;
         for ( ii = 0 ; ii < usize ; ii++ ) {
            w = uadj[ii] ;
            if ( map[w] == -1 ) {
/*
               ------------------
               add w to component
               ------------------
*/
               list[++last] = w ;
               map[w] = ncomp ;
            }
         }
      }
      ncomp++ ;
   }
}
IVfree(list) ;

return(mapIV) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   given a Graph g and map from vertices to components,
   fill counts[icomp] with the number of vertices in component icomp
   and fill weight[icomp] with their weight

   created -- 96feb25, cca
   -----------------------------------------------------------------
*/
void
Graph_componentStats (
   Graph   *g,
   int     map[],
   int     counts[],
   int     weights[]
) {
int   ncomp, nvtx, v, vcomp ;
int   *vwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL || map == NULL || counts == NULL || weights == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_componentStats(%p,%p,%p,%p)"
           "\n bad input\n", g, map, counts, weights) ;
   exit(-1) ;
}
/*
   --------------------
   fill the two vectors
   --------------------
*/
nvtx = g->nvtx ;
ncomp = 1 + IVmax(nvtx, map, &v) ;
if ( (vwghts = g->vwghts) != NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      vcomp = map[v] ;
      counts[vcomp]++ ;
      weights[vcomp] += vwghts[v] ;
   }
} else {
   for ( v = 0 ; v < nvtx ; v++ ) {
      vcomp = map[v] ;
      counts[vcomp]++ ;
   }
   IVcopy(ncomp, weights, counts) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   create and return a subgraph and a map from 
   its vertices to the vertices of the graph.

   g       -- graph from which to extract the subgraph
   icomp   -- component from which comes the vertices of the subgraph,
              icomp > 0
   compids -- component ids vector of graph
   pmap    -- pointer to hold address of map vector, the map from
              the subgraph's vertices to the graph's vertices

   return value -- pointer to subgraph Graph object

   created -- 95nov10, cca
   -------------------------------------------------------------------
*/
Graph *
Graph_subGraph (
   Graph   *g,
   int     icomp,
   int     compids[],
   int     **pmap
) {
Graph   *gsub ;
int     count, ii, iv, nvbnd, nvbndsub, nvtx, nvtxsub, nvtot, nvtotsub,
        v, vsub, vsize, w, wsub ;
int     *invmap, *map, *vadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL || icomp <= 0 || compids == NULL || pmap == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_subGraph(%p,%d,%p,%p)"
           "\n bad input\n", g, icomp, compids, pmap) ;
   exit(-1) ;
}
if ( g->type < 0 || g->type >= 2 ) {
   fprintf(stderr, "\n fatal error in Graph_subGraph(%p,%d,%p,%p)"
           "\n g->type = 0 or 1 currently supported\n", 
           g, icomp, compids, pmap) ;
   exit(-1) ;
}
nvtx  = g->nvtx  ;
nvbnd = g->nvbnd ;
nvtot = nvtx + nvbnd ;
/*
   ------------------------------------------------
   generate the map vectors
   map    : subgraph's vertices to graph's vertices
   invmap : graph's vertices to subgraph's vertices
   ------------------------------------------------
*/
map    = IVinit(nvtot, -1) ;
invmap = IVinit(nvtot, -1) ;
for ( v = 0, count = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] == icomp ) {
      map[count] = v ;
      invmap[v]  = count++ ;
   }
}
nvtxsub = count ;
/*
   ----------------------------------------------
   now get the boundary vertices for the subgraph
   ----------------------------------------------
*/
for ( iv = 0 ; iv < nvtxsub ; iv++ ) {
   v = map[iv] ;
   Graph_adjAndSize(g, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      w = vadj[ii] ;
      if ( w < nvtx ) {
         if ( compids[w] == 0 && invmap[w] == -1 ) {
            map[count] = w       ;
            invmap[w]  = count++ ;
         }
      } else if ( invmap[w] == -1 ) {
         map[count] = w       ;
         invmap[w]  = count++ ;
      }
   }
}
nvbndsub = count - nvtxsub;
nvtotsub = count ;
IVqsortUp(nvbndsub, &map[nvtxsub]) ;
for ( ii = nvtxsub ; ii < nvtotsub ; ii++ ) {
   v         = map[ii] ;
   invmap[v] = ii      ;
}
/*
   -----------------------
   initialize the subgraph
   -----------------------
*/
gsub = Graph_new() ;
Graph_init1(gsub, g->type, nvtxsub, nvbndsub, 
            0, IVL_CHUNKED, IVL_UNKNOWN) ;
/*
   ---------------------------------------------
   fill the adjacency structure of the subgraph
   note: the pointers of the subgraph point into
         the adjacency lists of the parent graph
         and the indices are overwritten.
   ---------------------------------------------
*/
for ( vsub = 0 ; vsub < nvtxsub ; vsub++ ) {
   v = map[vsub] ;
   Graph_adjAndSize(g, v, &vsize, &vadj) ;
   IVL_setPointerToList(gsub->adjIVL, vsub, vsize, vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      vadj[ii] = invmap[vadj[ii]] ;
   }
   IVqsortUp(vsize, vadj) ;
}
if ( nvbndsub > 0 ) {
   int   *ind = IVinit(nvtot, -1) ;
   for ( vsub = nvtxsub ; vsub < nvtotsub ; vsub++ ) {
      v = map[vsub] ;
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      for ( ii = 0, count = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( (wsub = invmap[w]) != -1 ) {
            ind[count++] = wsub ;
         }
      }
      IVqsortUp(count, ind) ;
      IVL_setList(gsub->adjIVL, vsub, count, ind) ;
   }
   IVfree(ind) ;
}
/*
   ---------------------------------
   fill vertex weights if applicable
   ---------------------------------
*/
if ( gsub->type % 2 == 1 ) {
   gsub->totvwght = 0 ;
   for ( vsub = 0 ; vsub < nvtotsub ; vsub++ ) {
      v = map[vsub] ;
      gsub->totvwght += g->vwghts[v] ;
      gsub->vwghts[vsub] = g->vwghts[v] ;
   }
} else {
   gsub->totvwght = gsub->nvtx ;
}
/*
   ----------------------------------------------------------------
   free the inverse map, create a new map[] vector the right size,
   copy the old map vector into the new and free the old map vector
   ----------------------------------------------------------------
*/
IVfree(invmap) ;
*pmap = IVinit(nvtotsub, -1) ;
IVcopy(nvtotsub, *pmap, map) ;
IVfree(map) ;

return(gsub) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return 1 if the graph is symmetric
   return 0 otherwise

   created -- 96oct31, cca
   ----------------------------------
*/
int
Graph_isSymmetric (
   Graph   *graph 
) {
int   ii, jj, nvtx, rc, v, vsize, w, wsize ;
int   *vadj, *wadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_isSymmetric(%p)"
           "\n bad input\n", graph) ;
   exit(-1) ;
}
/*
   ----------------------
   loop over the vertices
   ----------------------
*/
rc = 1 ;
nvtx = graph->nvtx ;
for ( v = 0 ; v < nvtx ; v++ ) {
   Graph_adjAndSize(graph, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      w = vadj[ii] ;
      Graph_adjAndSize(graph, w, &wsize, &wadj) ;
      for ( jj = 0 ; jj < wsize ; jj++ ) {
         if ( wadj[jj] == v ) {
            break ;
         }
      }
      if ( jj == wsize ) {
         fprintf(stdout, "\n edge (%d,%d) present, edge (%d,%d) is not",
                 v, w, w, v) ;
         rc = 0 ;
/*
         return(rc) ;
*/
      }
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
