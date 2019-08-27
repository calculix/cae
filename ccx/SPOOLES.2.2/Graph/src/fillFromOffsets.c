/*  fillFromOffsets.c  */

#include "../Graph.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- take an adjacency structure in the
              (offsets[neqns+1], adjncy[*]) form
              and load the Graph object

   g -- pointer to Graph object, must be initialized with nvtx = neqns
   neqns -- # of equations
   offsets -- offsets vector
   adjncy  -- big adjacency vector
      note, the adjacency for list v is found in
            adjncy[offsets[v]:offsets[v+1]-1]
      also note, offsets[] and adjncy[] must be zero based,
      if (offsets,adjncy) come from a harwell-boeing file, they use
      the fortran numbering, so each value must be decremented to
      conform with C's zero based numbering
   flag -- task flag
      flag = 0 --> just set the adjacency list for v to be that
                   found in adjncy[offsets[v]:offsets[v+1]-1]
      flag = 1 --> the input adjancency is just the upper triangle
                   (or strict upper triangle) as from a harwell-boeing
                   file. fill the Graph object with the full adjacency 
                   structure, including (v,v) edges

   created -- 96mar16, cca
   -------------------------------------------------------------------
*/
void
Graph_fillFromOffsets (
   Graph   *g,
   int     neqns,
   int     offsets[],
   int     adjncy[],
   int     flag
) {
IVL   *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if (  g == NULL 
   || neqns <= 0
   || offsets == NULL
   || adjncy == NULL
   || flag < 0
   || flag > 1 ) {
   fprintf(stderr, 
           "\n fatal error in Graph_fillFromOffsets(%p,%d,%p,%p,%d)"
           "\n bad input\n", g, neqns, offsets, adjncy, flag) ;
   exit(-1) ;
}
/*
   ---------------------------
   initialize the Graph object
   ---------------------------
*/
Graph_init1(g, 0, neqns, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
adjIVL = g->adjIVL ;
if ( flag == 0 ) {
   int   count, ii, nedge, v, w ;
   int   *list, *mark ;
/*
   ----------------------------------------------
   simple map, do not enforce symmetric structure
   ----------------------------------------------
*/
   list = IVinit(neqns, -1) ;
   mark = IVinit(neqns, -1) ;
   for ( v = 0, nedge = 0 ; v < neqns ; v++ ) {
      count = 0 ;
      for ( ii = offsets[v] ; ii < offsets[v+1] ; ii++ ) {
         w = adjncy[ii] ;
if ( v == neqns ) {
   fprintf(stdout, "\n hey there!! (v,w) = (%d,%d)", v, w) ;
}
         if ( 0 <= w && w < neqns && mark[w] != v ) {
            list[count++] = w ;
            mark[w] = v ;
         }
      }
      if ( mark[v] != v ) {
         list[count++] = v ;
         mark[v] = v ;
      }
      IVqsortUp(count, list) ;
      IVL_setList(adjIVL, v, count, list) ;
      nedge += count ;
   }
   g->totvwght = neqns ;
   g->totewght = g->nedges = nedge ;
/*
   ----------------------------
   now free the working storage
   ----------------------------
*/
   IVfree(list) ;
   IVfree(mark) ;
} else {
   int   ii, jj, u, v, vsize, w ;
   int   *head, *link, *list, *sizes, *vadj ;
   int   **p_adj ;
/*
   -------------------------------------------
   enforce symmetric structure and (v,v) edges
   make a first pass to check the input
   -------------------------------------------
*/
fprintf(stdout, "\n offsets") ;
IVfprintf(stdout, neqns+1, offsets) ;
   for ( v = 0 ; v < neqns ; v++ ) {
fprintf(stdout, "\n v = %d", v) ;
      for ( ii = offsets[v] ; ii < offsets[v+1] ; ii++ ) {
fprintf(stdout, "\n    w = %d", adjncy[ii]) ;
         if ( (w = adjncy[ii]) < v || neqns <= w ) {
            fprintf(stderr, 
               "\n fatal error in Graph_fillFromOffsets(%p,%d,%p,%p,%d)"
               "\n list %d, entry %d\n", g, neqns, offsets, adjncy,
               flag, v, w) ;
            exit(-1) ;
         }
      }
   }
   head  = IVinit(neqns, -1) ;
   link  = IVinit(neqns, -1) ;
   list  = IVinit(neqns, -1) ;
   sizes = IVinit(neqns, 0) ;
   p_adj = PIVinit(neqns) ;
   for ( v = 0 ; v < neqns ; v++ ) {
      vsize = 0 ;
/*
      -------------------------
      add edges to vertices < v
      -------------------------
*/
      while ( (u = head[v]) != -1 ) {
         head[v] = link[u] ;
         list[vsize++] = u ;
         if ( --sizes[u] > 0 ) {
            w = *(++p_adj[u]) ;
            link[u] = head[w] ;
            head[w] = u ;
         }
      }
/*
      -----------------
      add in edge (v,v)
      -----------------
*/
      list[vsize++] = v ;
      jj = vsize ;
/*
      -------------------------
      add edges to vertices > v
      -------------------------
*/
      for ( ii = offsets[v] ; ii < offsets[v+1] ; ii++ ) {
         if ( (w = adjncy[ii]) != v ) {
            list[vsize++] = w ;
         }
      }
/*
      ---------------------
      sort and set the list
      ---------------------
*/
      IVqsortUp(vsize, list) ;
      IVL_setList(adjIVL, v, vsize, list) ;
/*
      --------------------------------------------------
      link v to first vertex in its lists greater than v
      --------------------------------------------------
*/
      if ( jj < vsize ) {
         IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
         w        = vadj[jj]   ;
         link[v]  = head[w]    ;
         head[w]  = v          ;
         sizes[v] = vsize - jj ;
         p_adj[v] = &vadj[jj]  ;
      }
      g->nedges += vsize ;
   }
   g->totvwght = neqns     ;
   g->totewght = g->nedges ;
/*
   ----------------------------
   now free the working storage
   ----------------------------
*/
   IVfree(head)   ;
   IVfree(link)   ;
   IVfree(list)   ;
   IVfree(sizes)  ;
   PIVfree(p_adj) ;
}
return ; }

/*--------------------------------------------------------------------*/
