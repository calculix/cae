/*  compress.c  */

#include "../Graph.h"

#define MYDEBUG 0
#define MYTRACE 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   given a graph g and a fine-to-coarse map vector *mapIV,
   return a compressed graph with type coarseType.
   note, the compressed graph will have no trivial boundary
         even if the original graph did have a boundary.

   created -- 96mar02, cca
   ---------------------------------------------------------------
*/
Graph *
Graph_compress2 (
   Graph   *g,
   IV      *mapIV,
   int     coarseType
) {
/*
   -------------------------------------------
   check input and get dimensions and pointers
   -------------------------------------------
*/
if (  g == NULL || mapIV == NULL 
   || g->nvtx != IV_size(mapIV)
   || coarseType < 0 || 3 < coarseType ) {
   fprintf(stderr, "\n fatal error in Graph_compress2(%p,%p,%d)"
           "\n bad input\n", g, mapIV, coarseType) ;
   if ( g != NULL ) {
      Graph_writeStats(g, stderr) ;
   }
   if ( mapIV != NULL ) {
      IV_writeStats(mapIV, stderr) ;
   }
   exit(-1) ;
}
return(Graph_compress(g, IV_entries(mapIV), coarseType)) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   given a graph g and a fine-to-coarse map vector cmap[],
   return a compressed graph with type coarseType.
   note, the compressed graph will have no trivial boundary
         even if the original graph did have a boundary.

   created -- 95sep29, cca
   ---------------------------------------------------------------
*/
Graph *
Graph_compress (
   Graph   *g,
   int     cmap[],
   int     coarseType
) {
Graph   *g2 ;
int     ierr, ii, j, jj, jsize, J, Jsize, k, K, ncvtx, nvtx, wght ;
int     *head, *indices, *jind, *Jind, *jwghts, *Jwghts, 
        *link, *mark, *vwghts, *Vwghts ;
IVL     *adjIVL, *AdjIVL, *ewghtIVL, *EwghtIVL ;

#if MYTRACE > 0
fprintf(stdout, "\n just inside Graph_compress(%p,%p,%d)", 
        g, cmap, coarseType) ;
fflush(stdout) ;
#endif
/*
   -------------------------------------------
   check input and get dimensions and pointers
   -------------------------------------------
*/
if ( g == NULL || cmap == NULL || coarseType < 0 || 3 < coarseType ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n bad input\n", g, cmap, coarseType) ;
   exit(-1) ;
}
if ( (nvtx = g->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n nvtx = %d\n", g, cmap, coarseType, nvtx) ;
   exit(-1) ;
}
if ( (adjIVL = g->adjIVL) == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n adjIVL is NULL\n", g, cmap, coarseType) ;
   exit(-1) ;
}
if ( g->type % 2 == 1 && (vwghts = g->vwghts) == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n g->type = %d and g->vwghts is NULL\n", 
           g, cmap, coarseType, g->type) ;
   exit(-1) ;
}
if ( g->type >= 2 && (ewghtIVL = g->ewghtIVL) == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n g->type = %d and g->ewghtIVL is NULL\n", 
           g, cmap, coarseType, g->type) ;
   exit(-1) ;
}
if ( IVmin(nvtx, cmap, &j) < 0 ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n IVmin(cmap) = %d\n", 
           g, cmap, coarseType, IVmin(nvtx, cmap, &j)) ;
   exit(-1) ;
}
ncvtx = 1 + IVmax(nvtx, cmap, &j) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ncvtx = %d", ncvtx) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------
   initialize the coarse graph object
   ----------------------------------
*/
g2 = Graph_new() ;
Graph_init1(g2, coarseType, ncvtx, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
if ( (AdjIVL = g2->adjIVL) == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n AdjIVL is NULL\n", g, cmap, coarseType) ;
   exit(-1) ;
}
if ( g2->type % 2 == 1 && (Vwghts = g2->vwghts) == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n g2->type = %d and g2->vwghts is NULL\n", 
           g, cmap, coarseType, g2->type) ;
   exit(-1) ;
}
if ( g2->type >= 2 && (EwghtIVL = g2->ewghtIVL) == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_compress(%p,%p,%d)"
           "\n g2->type = %d and g2->ewghtIVL is NULL\n", 
           g, cmap, coarseType, g2->type) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n after initializing the coarse graph object") ;
Graph_writeStats(g2, stdout) ;
fflush(stdout) ;
#endif
/*
   --------------------------------
   set up the head and link vectors
   --------------------------------
*/
head = IVinit(ncvtx, -1) ;
link = IVinit(nvtx,  -1) ;
for ( j = 0 ; j < nvtx ; j++ ) {
   J       = cmap[j] ;
   link[j] = head[J] ;
   head[J] =    j    ;
}
/*
   ------------------------------------------------
   fill the adjacency structure of the coarse graph
   ------------------------------------------------
*/
indices = IVinit2(ncvtx) ;
mark    = IVinit(ncvtx, -1) ;
for ( J = 0 ; J < ncvtx ; J++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n checking out J = %d", J) ;
#endif
   Jsize = 0 ;
   for ( j = head[J] ; j != -1 ; j = link[j] ) {
      IVL_listAndSize(adjIVL, j, &jsize, &jind) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    adj j = %d : ", j) ;
      IVfp80(stdout, jsize, jind, 20, &ierr) ;
#endif
      for ( ii = 0 ; ii < jsize ; ii++ ) {
         if ( (k = jind[ii]) < nvtx ) {
            K = cmap[k] ;
#if MYDEBUG > 0
            fprintf(stdout, "\n    k = %d, K = %d", k, K) ;
#endif
            if ( mark[K] != J ) {
#if MYDEBUG > 0
               fprintf(stdout, ", added") ;
#endif
               mark[K] = J ;
               indices[Jsize++] = K ;
            }
         }
      }
   }
   if ( Jsize > 0 ) {
      IVqsortUp(Jsize, indices) ;
   }
   IVL_setList(AdjIVL, J, Jsize, indices) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n setting list %d :", J) ;
   IVfp80(stdout, Jsize, indices, 20, &ierr) ;
#endif
}
g2->nedges = AdjIVL->tsize ;
IVfree(indices) ;
IVfree(mark) ;

if ( coarseType % 2 == 1 ) {
/*
   -------------------------------------------
   fill the vertex weights of the coarse graph
   -------------------------------------------
*/
   for ( J = 0 ; J < ncvtx ; J++ ) {
      wght = 0 ;
      for ( j = head[J] ; j != -1 ; j = link[j] ) {
         if ( g->type % 2 == 1 ) {
            wght += vwghts[j] ;
         } else {
            wght++ ;
         }
      }
      Vwghts[J] = wght ;
    }
   g2->totvwght = IVsum(ncvtx, Vwghts) ;
#if MYDEBUG > 0
   { int ierr ;
      fprintf(stdout, "\n Vwghts(%d) : ", ncvtx) ;
      IVfp80(stdout, ncvtx, Vwghts, 80, &ierr) ;
      fflush(stdout) ;
   }
#endif
} else {
/*
   -------------------------------------------------
   coarse graph does not have vertex weights,
   set total vertex weight to the number of vertices
   -------------------------------------------------
*/
   g2->totvwght = ncvtx ;
}
if ( coarseType >= 2 ) {
/*
   -----------------------------------------
   fill the edge weights of the coarse graph
   -----------------------------------------
*/
   for ( J = 0 ; J < ncvtx ; J++ ) {
      IVL_listAndSize(AdjIVL, J, &Jsize, &Jind) ;
      IVL_setList(EwghtIVL, J, Jsize, NULL) ;
   }
   if ( g->type >= 2 ) {
/*
      ---------------------------
      fine graph had edge weights
      ---------------------------
*/
      for ( j = 0 ; j < nvtx ; j++ ) {
         J = cmap[j] ;
         IVL_listAndSize(adjIVL,   j, &jsize, &jind) ;
         IVL_listAndSize(ewghtIVL, j, &jsize, &jwghts) ;
         IVL_listAndSize(AdjIVL,   J, &Jsize, &Jind) ;
         IVL_listAndSize(EwghtIVL, J, &Jsize, &Jwghts) ;
         for ( ii = 0 ; ii < jsize ; ii++ ) {
            k = jind[ii] ;
            if ( k < nvtx ) {
               K = cmap[k] ;
               for ( jj = 0 ; jj < Jsize ; jj++ ) {
                  if ( Jind[jj] == K ) {
                     Jwghts[jj] += jwghts[ii] ;
                     break ;
                  }
               }
            }
         }
      }
   } else {
/*
      ------------------------------------
      fine graph did not have edge weights
      ------------------------------------
*/
      for ( j = 0 ; j < nvtx ; j++ ) {
         J = cmap[j] ;
         IVL_listAndSize(adjIVL,   j, &jsize, &jind) ;
         IVL_listAndSize(AdjIVL,   J, &Jsize, &Jind) ;
         IVL_listAndSize(EwghtIVL, J, &Jsize, &Jwghts) ;
         for ( ii = 0 ; ii < jsize ; ii++ ) {
            k = jind[ii] ;
            if ( k < nvtx ) {
               K = cmap[k] ;
               for ( jj = 0 ; jj < Jsize ; jj++ ) {
                  if ( Jind[jj] == K ) {
                     Jwghts[jj]++ ;
                     break ;
                  }
               }
            }
         }
      }
   }
   g2->totewght = IVL_sum(EwghtIVL) ;
} else {
/*
   --------------------------------------------
   coarse graph does not have edge weights,
   set total edge weight to the number of edges
   --------------------------------------------
*/
   g2->totewght = g2->nedges ;
}
/*
   ------------------------------
   free the head and link vectors
   ------------------------------
*/
IVfree(head) ;
IVfree(link) ;

return(g2) ; }

/*--------------------------------------------------------------------*/
