/*  makeGraphXbyX.c  */

#include "../BPG.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- return the X by X graph where (x1,x2) is an edge
              if there exists a y in Y such that (x1,y) and
              (x2,y) are edges in the bipartite graph.

   created -- 95dec06, cca
              to be used in the BKL object.
   -----------------------------------------------------------
*/
Graph *
BPG_makeGraphXbyX (
   BPG   *bpg
) {
Graph   *graph, *gXbyX ;
int     count, ii, jj, nX, x, xsize, y, ysize, z ;
int     *list, *mark, *xadj, *yadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL ) {
   fprintf(stdout, "\n fatal error in BPG_makeGraphXbyX(%p)"
           "\n bad input\n", bpg) ;
   exit(-1) ;
}
/*
   ----------------------
   check for quick return
   ----------------------
*/
if ( (graph = bpg->graph) == NULL || (nX = bpg->nX) <= 0 ) {
   return(NULL) ;
}
/*
   --------------------
   initialize the graph
   --------------------
*/
gXbyX = Graph_new() ;
Graph_init1(gXbyX, graph->type, nX, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
/*
   --------------
   fill the graph
   --------------
*/
mark = IVinit(nX, -1) ;
list = IVinit(nX, -1) ;
for ( x = 0 ; x < nX ; x++ ) {
   Graph_adjAndSize(graph, x, &xsize, &xadj) ;
   mark[x] = x ;
   for ( ii = 0, count = 0 ; ii < xsize ; ii++ ) {
      y = xadj[ii] ;
      Graph_adjAndSize(graph, y, &ysize, &yadj) ;
      for ( jj = 0 ; jj < ysize ; jj++ ) {
         z = yadj[jj] ;
         if ( mark[z] != x ) {
            mark[z] = x ;
            list[count++] = z ;
         }
      }
   }
   if ( count > 0 ) {
      IVqsortUp(count, list) ;
      IVL_setList(gXbyX->adjIVL, x, count, list) ;
   }
}
IVfree(list) ;
IVfree(mark) ;
/*
   ---------------------------------------
   set vertex weight vector if appropriate
   ---------------------------------------
*/
if ( graph->type % 2 == 1 ) {
   IVcopy(nX, gXbyX->vwghts, graph->vwghts) ;
}

return(gXbyX) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- return the Y by Y graph where (y1,y2) is an edge
              if there exists a x in X such that (x,y1) and
              (x,y2) are edges in the bipartite graph.

   created -- 95dec07, cca
   -----------------------------------------------------------
*/
Graph *
BPG_makeGraphYbyY (
   BPG   *bpg
) {
Graph   *graph, *gYbyY ;
int     count, ii, jj, nX, nY, x, xsize, y, ysize, z ;
int     *list, *mark, *xadj, *yadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL ) {
   fprintf(stdout, "\n fatal error in BPG_makeGraphXbyX(%p)"
           "\n bad input\n", bpg) ;
   exit(-1) ;
}
/*
   ----------------------
   check for quick return
   ----------------------
*/
if ( (graph = bpg->graph) == NULL || (nY = bpg->nY) <= 0 ) {
   return(NULL) ;
}
nX = bpg->nX ;
/*
   --------------------
   initialize the graph
   --------------------
*/
gYbyY = Graph_new() ;
Graph_init1(gYbyY, graph->type, nY, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
/*
   --------------
   fill the graph
   --------------
*/
mark = IVinit(nY, -1) ;
list = IVinit(nY, -1) ;
for ( y = 0 ; y < nY ; y++ ) {
   Graph_adjAndSize(graph, nX + y, &ysize, &yadj) ;
   mark[y] = y ;
   for ( ii = 0, count = 0 ; ii < ysize ; ii++ ) {
      x = yadj[ii] ;
      Graph_adjAndSize(graph, x, &xsize, &xadj) ;
      for ( jj = 0 ; jj < xsize ; jj++ ) {
         z = xadj[jj] ;
         if ( mark[z] != y ) {
            mark[z] = y ;
            list[count++] = z ;
         }
      }
   }
   if ( count > 0 ) {
      IVqsortUp(count, list) ;
      IVL_setList(gYbyY->adjIVL, nX + y, count, list) ;
   }
}
IVfree(list) ;
IVfree(mark) ;
/*
   ---------------------------------------
   set vertex weight vector if appropriate
   ---------------------------------------
*/
if ( graph->type % 2 == 1 ) {
   IVcopy(nY, gYbyY->vwghts, graph->vwghts + nX) ;
}

return(gYbyY) ; }

/*--------------------------------------------------------------------*/
