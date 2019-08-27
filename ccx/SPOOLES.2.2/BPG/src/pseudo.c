/*  pseudo.c  */

#include "../BPG.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   return a pseudoperipheral node

   created -- 95oct07, cca
   ------------------------------
*/
int
BPG_pseudoperipheralnode (
   BPG   *bpg,
   int   seed
) {
int   last, mate, oldrad, rad, root, tag ;
int   *dist, *list, *mark ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_pseudoperipheralnode(%p,%d)"
           "\n bad input\n", bpg, seed) ;
   exit(-1) ;
}
if ( seed < 0 ) {
   seed = - seed ;
}
seed = seed % (bpg->nX + bpg->nY) ;
list = IVinit(bpg->nX + bpg->nY, -1) ;
dist = IVinit(bpg->nX + bpg->nY, -1) ;
mark = IVinit(bpg->nX + bpg->nY, -1) ;
oldrad = 0 ;
/*
   -------------------------------------------
   drop a level structure from the seed vertex
   -------------------------------------------
*/
tag = 1 ;
last = BPG_levelStructure(bpg, seed, list, dist, mark, tag) ;
mate = list[last] ;
rad  = dist[mate] ;
#if MYDEBUG > 0
fprintf(stdout, "\n BPG_pseudoperipheralnode") ;
fprintf(stdout, "\n node %d, mate %d, rad %d", seed, mate, rad) ;
fflush(stdout) ;
#endif
/*
   -----------------------------------------------
   loop while vertex with a larger radius is found
   -----------------------------------------------
*/
while ( oldrad < rad ) {
   oldrad = rad ;
   root   = mate ;
   tag++ ;
   last = BPG_levelStructure(bpg, root, list, dist, mark, tag) ;
   mate = list[last] ;
   rad  = dist[mate] ;
#if MYDEBUG > 0
   fprintf(stdout, "\n node %d, mate %d, rad %d", root, mate, rad) ;
   fflush(stdout) ;
#endif
}
#if MYDEBUG > 1
{ int root ;
for ( root = 0, tag++ ; root < bpg->nX ; root++, tag++ ) {
   last = BPG_levelStructure(bpg, root, list, dist, mark, tag) ;
   mate = list[last] ;
   rad  = dist[mate] ;
   fprintf(stdout, "\n node %d, mate %d, rad %d", root, mate, rad) ;
   fflush(stdout) ;
}
}
#endif
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(list) ;
IVfree(dist) ;
IVfree(mark) ;

return(root) ; }

/*--------------------------------------------------------------------*/
#define MYDEBUG 0
/*
   ----------------------------------------------------
   return value -- # of vertices in the level structure

   created -- 95oct07, cca
   ----------------------------------------------------
*/
int
BPG_levelStructure (
   BPG   *bpg,
   int   root,
   int   list[],
   int   dist[],
   int   mark[],
   int   tag
) {
int   ii, jj, last, now, u, usize, v, vsize, w ;
int   *uadj, *vadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || root < 0 || root >= bpg->nX + bpg->nY
   || list == NULL || dist == NULL || mark == NULL ) {
   fprintf(stderr, 
           "\n fatal error in BPG_levelStructure(%p,%d,%p,%p,%p,%d)"
           "\n bad input\n", bpg, root, list, dist, mark, tag) ;
   exit(-1) ;
}

#if MYDEBUG > 0
fprintf(stdout, "\n inside BPG_levelStructure(%p,%d,%p,%p,%p,%d)",
        bpg, root, list, dist, mark, tag) ;
fflush(stdout) ;
#endif

now = last = 0 ;
list[0] = root ;
dist[root] = 0 ;
mark[root] = tag ;
while ( now <= last ) {
   u = list[now++] ;
#if MYDEBUG > 0
   fprintf(stdout, "\n u = %d", u) ;
   fflush(stdout) ;
#endif
   Graph_adjAndSize(bpg->graph, u, &usize, &uadj) ;
   for ( ii = 0 ; ii < usize ; ii++ ) {
      v = uadj[ii] ;
      Graph_adjAndSize(bpg->graph, v, &vsize, &vadj) ;
      for ( jj = 0 ; jj < vsize ; jj++ ) {
         w = vadj[jj] ;
#if MYDEBUG > 0
         fprintf(stdout, "\n    w = %d", w) ;
#endif
         if ( mark[w] != tag ) {
#if MYDEBUG > 0
            fprintf(stdout, ", adding to list, dist = %d", dist[u] + 1);
            fflush(stdout) ;
#endif
            mark[w] = tag ;
            list[++last] = w ;
            dist[w] = dist[u] + 1 ;
         }
      }
   }
}

return(last) ; }

/*--------------------------------------------------------------------*/
