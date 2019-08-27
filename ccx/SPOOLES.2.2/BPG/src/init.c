/*  init.c  */

#include "../BPG.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   initialize a BPG object, set bpg->nX = nX, bpg->nY = nY,
   and bpg->g = g. convert g into bipartite form, i.e.,
      adj(v) := adj(v) \cap {nX, ..., nX + nY - 1} for v in [0, nX)
      adj(v) := adj(v) \cap {0, ..., nX - 1}       for v in [nX, nX+nY)

   created -- 95oct06, cca
   --------------------------------------------------------------------
*/
void
BPG_init (
   BPG     *bpg,
   int     nX,
   int     nY,
   Graph   *graph
) {
int   ierr, ii, jj, v, vsize, w ;
int   *vadj ;
IVL   *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || nX <= 0 || nY <= 0 || graph == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_init(%p,%d,%d,%p)"
           "\n bad input\n", bpg, nX, nY, graph) ;
   exit(-1) ;
}
/*
   --------------
   set the fields
   --------------
*/
BPG_clearData(bpg) ;
bpg->nX    =   nX  ;
bpg->nY    =   nY  ;
bpg->graph = graph ;
/*
   ------------------------------------
   work on the graph's adjacency lists,
   enforcing the bipartite property
   ------------------------------------
*/
adjIVL = graph->adjIVL ;
for ( v = 0 ; v < nX ; v++ ) {
   IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n old %d : ", v) ;
   IVfp80(stdout, vsize, vadj, 10, &ierr) ;
   fflush(stdout) ;
#endif
   ii = 0, jj = vsize - 1 ;
   while ( ii <= jj ) {
      w = vadj[ii] ;
      if ( nX <= w && w < nX + nY ) {
         ii++ ;
      } else {
         vadj[ii] = vadj[jj] ;
         vadj[jj] =    w     ;
         jj-- ;
      }
   }
   vsize = ii ;
#if MYDEBUG > 0
   fprintf(stdout, "\n new %d : ", v) ;
   IVfp80(stdout, vsize, vadj, 10, &ierr) ;
   fflush(stdout) ;
#endif
   IVL_setList(adjIVL, v, vsize, NULL) ;
}
for ( v = nX ; v < nX + nY ; v++ ) {
   IVL_listAndSize(adjIVL, v, &vsize, &vadj) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n old %d : ", v) ;
   IVfp80(stdout, vsize, vadj, 10, &ierr) ;
   fflush(stdout) ;
#endif
   ii = 0, jj = vsize - 1 ;
   while ( ii <= jj ) {
      w = vadj[ii] ;
      if ( 0 <= w && w < nX ) {
         ii++ ;
      } else {
         vadj[ii] = vadj[jj] ;
         vadj[jj] =    w     ;
         jj-- ;
      }
   }
   vsize = ii ;
   IVL_setList(adjIVL, v, vsize, NULL) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n new %d : ", v) ;
   IVfp80(stdout, vsize, vadj, 10, &ierr) ;
   fflush(stdout) ;
#endif
}

return ; }

/*--------------------------------------------------------------------*/
#define MYDEBUG 0
/*
   -----------------------------------------------------------
   initialize from a graph given a color vector and two target
   vectors. this method can be used to get the bipartite graph
   of the boundary of two components.

   graph  -- graph from which to extract the bipartite graph
   colors -- vector of colors for the vertices
   cX     -- color for the X vertices
   cY     -- color for the Y vertices
   cmap   -- map from vertices in g to vertices in bpg
   indX   -- vector to hold the global indices in X
   indY   -- vector to hold the global indices in Y

   created -- 95oct06, cca
   -----------------------------------------------------------
*/
void
BPG_initFromColoring (
   BPG     *bpg,
   Graph   *graph,
   int     colors[],
   int     cX,
   int     cY,
   int     cmap[],
   int     indX[],
   int     indY[]
) {
Graph   *bpg_g ;
int     count, ierr, ii, iX, iy, msize, nV, nX, nY, v, vsize, w, x, y ;
int     *ewghts, *list, *vadj, *vewghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || graph == NULL || colors == NULL 
   || cX < 0 || cY < 0 || cX == cY || cmap == NULL ) {
   fprintf(stderr, 
           "\n fatal error in BPG_initFromColoring(%p,%p,%p,%d,%d,%p)"
           "\n bad input\n", bpg, graph, colors, cX, cY, cmap) ;
   exit(-1) ;
}
BPG_clearData(bpg) ;
nV = graph->nvtx ;
IVfill(nV, cmap, -1) ;
#if MYDEBUG > 0
fprintf(stdout, 
        "\n\n just inside BPG_initFromColoring(%p,%p,%p,%d,%d,%p)",
        bpg, graph, colors, cX, cY, cmap) ;
fflush(stdout) ;
#endif
/*
   ------------------
   get the X vertices
   ------------------
*/
for ( v = 0, nX = 0 ; v < nV ; v++ ) {
   if ( colors[v] == cX ) {
      indX[nX] = v    ;
      cmap[v]  = nX++ ;
   }
}
/*
   ----------------------------------------------------
   now get Y = {v | v \in \bnd{X} and colors[v] == cY }
   ----------------------------------------------------
*/
nY = 0 ;
for ( iX = 0 ; iX < nX ; iX++ ) {
   v = indX[iX] ;
   Graph_adjAndSize(graph, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      if (  (w = vadj[ii]) < nV
         && colors[w] == cY && cmap[w] == -1 ) {
         indY[nY] = w ;
         cmap[w]  = nX + nY++ ;
      }
   }
}
bpg->nX = nX ;
bpg->nY = nY ;
#if MYDEBUG > 0
fprintf(stdout, "\n indX(%d) :", nX) ;
IVfp80(stdout, nX, indX, 15, &ierr) ;
fprintf(stdout, "\n indY(%d) :", nY) ;
IVfp80(stdout, nY, indY, 15, &ierr) ;
fflush(stdout) ;
#endif
#if MYDEBUG > 1
fprintf(stdout, "\n cmap(%d) :", g->nvtx) ;
IVfp80(stdout, graph->nvtx, cmap, 15, &ierr) ;
fflush(stdout) ;
#endif
if ( bpg->nX == 0 || bpg->nY == 0 ) {
   fprintf(stderr, "\n fatal error in BPG_initFromColoring"
           "\n nX = %d, nY = %d", nX, nY) ;
   fprintf(stderr, "\n colors") ;
   IVfp80(stderr, nV, colors, 80, &ierr) ;
   fprintf(stderr, "\n graph") ;
   Graph_writeForHumanEye(graph, stderr) ;
   exit(-1) ;
}
/*
   --------------------------------------
   initialize the bipartite graph's graph
   --------------------------------------
*/
bpg->graph = bpg_g = Graph_new() ;
Graph_init1(bpg_g, graph->type, 
            nX + nY, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
/*
   -------------------------------------
   set the vertex weights if appropriate
   -------------------------------------
*/
if ( graph->type % 2 == 1 ) {
   for ( x = 0 ; x < nX ; x++ ) {
      v = indX[x] ;
      bpg_g->vwghts[x] = graph->vwghts[v] ;
   }
   for ( iy = 0, y = nX ; iy < nY ; iy++, y++ ) {
      v = indY[iy] ;
      bpg_g->vwghts[y] = graph->vwghts[v] ;
   }
   bpg_g->totvwght = IVsum(nX + nY, bpg_g->vwghts) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n vertex weights for X") ;
   IVfp80(stdout, nX, bpg_g->vwghts, 80, &ierr) ;
   fprintf(stdout, "\n vertex weights for Y") ;
   IVfp80(stdout, nY, bpg_g->vwghts + nX, 80, &ierr) ;
   fprintf(stdout, "\n w(X) = %d", IVsum(nX, bpg_g->vwghts)) ;
   fprintf(stdout, "\n w(Y) = %d", IVsum(nY, bpg_g->vwghts + nX)) ;
   fflush(stdout) ;
#endif
}
/*
   -----------------------------------------------
   set the adjacency or adjacency and edge weights
   -----------------------------------------------
*/
if ( graph->type < 2 ) {
/*
   -------------------------------
   adjacency only, no edge weights
   -------------------------------
*/
   msize  = IVL_maxListSize(graph->adjIVL) ;
   list   = IVinit2(msize) ;
   for ( x = 0 ; x < nX ; x++ ) {
      v = indX[x] ;
      Graph_adjAndSize(graph, v, &vsize, &vadj) ;
      for ( ii = 0, count = 0 ; ii < vsize ; ii++ ) {
         if ( (w = vadj[ii]) < nV && colors[w] == cY ) {
            list[count++] = cmap[w] ;
         }
      }
      IVqsortUp(count, list) ;
      IVL_setList(bpg_g->adjIVL, x, count, list) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n setting x (%d) list %d : ", x, indX[x]) ;
      IVfp80(stdout, count, list, 20, &ierr) ;
      fflush(stdout) ;
#endif
   }
   for ( iy = 0, y = nX ; iy < nY ; iy++, y++ ) {
      v = indY[iy] ;
      Graph_adjAndSize(graph, v, &vsize, &vadj) ;
      for ( ii = 0, count = 0 ; ii < vsize ; ii++ ) {
         if ( (w = vadj[ii]) < nV && colors[w] == cX ) {
            list[count++] = cmap[w] ;
         }
      }
      IVqsortUp(count, list) ;
      IVL_setList(bpg_g->adjIVL, y, count, list) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n setting y (%d) list %d : ", y, indY[iy]) ;
      IVfp80(stdout, count, list, 20, &ierr) ;
      fflush(stdout) ;
#endif
   }
   IVfree(list) ;
} else {
/*
   --------------------------
   adjacency and edge weights
   --------------------------
*/
   msize  = IVL_maxListSize(graph->adjIVL) ;
   list   = IVinit2(msize) ;
   ewghts = IVinit2(msize) ;
   for ( x = 0 ; x < nX ; x++ ) {
      v = indX[x] ;
      Graph_adjAndEweights(graph, v, &vsize, &vadj, &vewghts) ;
      for ( ii = 0, count = 0 ; ii < vsize ; ii++ ) {
         if ( (w = vadj[ii]) < nV && colors[w] == cY ) {
            list[count]   = cmap[w] ;
            ewghts[count] = vewghts[ii] ;
            count++ ;
         }
      }
      IV2qsortUp(count, list, ewghts) ;
      IVL_setList(bpg_g->adjIVL,   x, count, list) ;
      IVL_setList(bpg_g->ewghtIVL, x, count, ewghts) ;
   }
   for ( iy = 0, y = nX ; iy < nY ; iy++, y++ ) {
      v = indY[iy] ;
      Graph_adjAndEweights(graph, v, &vsize, &vadj, &vewghts) ;
      for ( ii = 0, count = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( colors[w] == cX ) {
            list[count]   = cmap[w] ;
            ewghts[count] = vewghts[ii] ;
            count++ ;
         }
      }
      IV2qsortUp(count, list, ewghts) ;
      IVL_setList(bpg_g->adjIVL,   y, count, list) ;
      IVL_setList(bpg_g->ewghtIVL, y, count, ewghts) ;
   }
   IVfree(list)   ;
   IVfree(ewghts) ;
}
bpg_g->nedges = IVsum(nX + nY, graph->adjIVL->sizes) ;

return ; }

/*--------------------------------------------------------------------*/
