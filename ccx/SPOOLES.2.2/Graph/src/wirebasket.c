/*  wirebasket.c  */

#include "../Graph.h"

#define MYDEBUG 1

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   on input
      stages[v] = 0 --> vertex is in a domain
      stages[v] > 0 --> vertex is in the multisector
   this remains true on output, but the stage of a multisector
   vertex is equal to the number of different domains that
   are found within radius edges from itself
   
   created -- 97jul30, cca
   ---------------------------------------------------------------------
*/
void
Graph_wirebasketStages (
   Graph   *graph,
   IV      *stagesIV,
   int     radius
) {
int   count, idom, ierr, ii, last, ndom, now, nvtx, u, v, vsize, w ;
int   *dist, *dmark, *domids, *list, *stages, *vadj, *vmark ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || stagesIV == NULL || radius < 0 ) {
   fprintf(stderr, "\n fatal error in Graph_wirebasketStages(%p,%p,%d)"
           "\n bad input\n", graph, stagesIV, radius) ;
   exit(-1) ;
}
IV_sizeAndEntries(stagesIV, &nvtx, &stages) ;
if ( nvtx != graph->nvtx || stages == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_wirebasketStages(%p,%p,%d)"
           "\n stages->nvtx = %d, graph->nvtx = %d, stages = %p\n", 
           graph, stagesIV, nvtx, radius, graph->nvtx, stages) ;
   exit(-1) ;
}
/*
   -----------------------------------------------
   fill domids[] 
      domids[v] = -1 --> v is in the multisector
      domids[v] != -1 --> v is in domain domids[v]
   -----------------------------------------------
*/
domids = IVinit(nvtx, -1) ;
list   = IVinit(nvtx, -1) ;
for ( u = 0, ndom = 0 ; u < nvtx ; u++ ) {
#if MYDEBUG > 1
   fprintf(stdout, "\n vertex %d, stage %d, domid %d", 
           u, stages[u], domids[u]) ;
#endif
   if ( stages[u] == 0 && domids[u] == -1 ) {
#if MYDEBUG > 1
      fprintf(stdout, "\n vertex %d starts domain %d", u, ndom) ;
#endif
      list[now = last = 0] = u ;
      domids[u] = ndom ;
      while ( now <= last ) {
         v = list[now++] ;
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            if ( stages[w] == 0 && domids[w] == -1 ) {
#if MYDEBUG > 1
               fprintf(stdout, "\n    adding vertex %d", w) ;
#endif
               domids[w]    = ndom ;
               list[++last] = w ;
            }
         }
      }
      ndom++ ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n domids") ;
fprintf(stdout, "\n %d", nvtx) ;
IVfp80(stdout, nvtx, domids, 80, &ierr) ;
#endif
/*
   --------------------------------------------------------
   fill the stages of the multisector vertices based on the
   number of different domains that are within radius steps
   --------------------------------------------------------
*/
dmark = IVinit(ndom, -1) ;
vmark = IVinit(nvtx, -1) ;
dist  = IVinit(nvtx, -1) ;
for ( u = 0 ; u < nvtx ; u++ ) {
   if ( stages[u] != 0 ) {
#if MYDEBUG > 1
      fprintf(stdout, "\n\n checking out schur vertex %d", u) ;
#endif
      list[now = last = 0] = u ;
      vmark[u] = u ;
      dist[u]  = 0 ;
      count    = 0 ;
      while ( now <= last ) {
         v = list[now++] ;
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
#if MYDEBUG > 1
         fprintf(stdout, 
                 "\n   removing vertex %d from list, dist %d", 
                 v, dist[v]) ;
#endif
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
#if MYDEBUG > 1
               fprintf(stdout, 
                      "\n      adjacent vertex %d, vmark %d, domids %d",
                      w, vmark[w], domids[w]) ;
#endif
            if ( vmark[w] != u ) {
               vmark[w] = u ;
               if ( (idom = domids[w]) != -1 ) { 
                  if ( dmark[idom] != u ) {
#if MYDEBUG > 1
                     fprintf(stdout, ", marking domain %d", idom) ;
#endif
                     dmark[idom] = u ;
                     count++ ;
                  }
               } else if ( dist[v] < radius - 1 ) {
#if MYDEBUG > 1
                     fprintf(stdout, ", adding %d to list", w) ;
#endif
                  dist[w] = dist[v] + 1 ;
                  list[++last] = w ;
               }
            }
         }
      }
      stages[u] = count ;
#if MYDEBUG > 1
      fprintf(stdout, "\n   setting stage of %d to be %d", u, count) ;
#endif
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n stages") ;
fprintf(stdout, "\n %d", nvtx) ;
IVfp80(stdout, nvtx, stages, 80, &ierr) ;
#endif
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(domids) ;
IVfree(list)   ;
IVfree(dmark)  ;
IVfree(vmark)  ;
IVfree(dist)   ;

return ; }

/*--------------------------------------------------------------------*/
