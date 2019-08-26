/*  equivMap.c  */

#include "../Graph.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   fill an IV object with an equivalence map
   if ( map[u] == map[v] ) then
      then u and v are adjacent to the same vertices
   endif
   NOTE : each empty list is mapped to a different value

   return value -- IV object that contains the equivalence map

   created -- 96feb23, cca
   -----------------------------------------------------------
*/
IV *
Graph_equivMap (
   Graph   *g
) {
int   ierr, ii, ismapped, jj, nvtx, nvtx2, u, usize, v, vsize, vsum ;
int   *chksum, *eqmap, *mark, *uadj, *vadj ;
IV    *eqmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL || (nvtx = g->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in Graph_equivMap(%p)"
           "\n bad input\n", g) ;
   exit(-1) ;
}
/*
   -----------------------------------------
   initialize map and set up working storage
   -----------------------------------------
*/
eqmapIV = IV_new() ;
IV_init(eqmapIV, nvtx, NULL) ;
eqmap = IV_entries(eqmapIV) ;
IVfill(nvtx, eqmap, -1) ;
mark   = IVinit(nvtx, -1) ;
chksum = IVinit(nvtx, -1) ;
/*
   -----------------------------------------
   loop over the vertices in ascending order
   -----------------------------------------
*/
nvtx2 = 0 ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( eqmap[v] == -1 ) {
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      if ( vsize == 0 ) {
/*
         ------------------------------
         list is empty, map to new list
         ------------------------------
*/
#if MYDEBUG > 1
         fprintf(stdout, 
                 "\n vertex %d, empty list, map to %d", v, eqmap[v]) ;
         fflush(stdout) ;
#endif
         eqmap[v] = nvtx2++ ;
      } else {
/*
         ----------------
         compute checksum
         ----------------
*/
         vsum = v ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            if ( vadj[ii] != v ) {
               vsum += vadj[ii] ;
            }
         }
         chksum[v] = vsum ;
#if MYDEBUG > 1
         fprintf(stdout, 
                 "\n vertex %d, checksum = %d", v, chksum[v]) ;
         fflush(stdout) ;
#endif
/*
         -------------------------------------------------------
         look at all adjacent vertices with numbers lower than v
         -------------------------------------------------------
*/
         ismapped = 0 ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            if ( (u = vadj[ii]) < v && chksum[u] == vsum ) {
/*
               ------------------------------
               u and v have the same checksum
               ------------------------------
*/
#if MYDEBUG > 1
               fprintf(stdout, 
                       "\n     checking vertex %d, checksum = %d", 
                       u, chksum[u]) ;
               fflush(stdout) ;
#endif
               Graph_adjAndSize(g, u, &usize, &uadj) ;
               if ( vsize == usize ) {
/*
                  ----------------------------------------------
                  the adjacency lists are the same size, compare
                  ----------------------------------------------
*/
                  if ( ismapped == 0 ) {
/*
                     ---------------------------
                     mark all vertices in adj(v)
                     ---------------------------
*/
                     mark[v] = v ;
                     for ( jj = 0 ; jj < vsize ; jj++ ) {
                         mark[vadj[jj]] = v ;
                     }
                     ismapped = 1 ;
                  }
/*
                  -------------------------------------------------
                  check to see if all vertices in adj(u) are marked
                  -------------------------------------------------
*/
                  for ( jj = 0 ; jj < usize ; jj++ ) {
                     if ( mark[uadj[jj]] != v ) {
#if MYDEBUG > 1
                        fprintf(stdout, 
       "\n   vertex %d is adjacent to %d but not %d", uadj[jj], u, v) ;
#endif
                        break ;
                     }
                  }
                  if ( jj == usize ) {
/*
                     ----------------------------------
                     lists are identical, set map for v 
                     to be the same as the map for u
                     ----------------------------------
*/
#if MYDEBUG > 1
                     fprintf(stdout, 
                    "\n    lists are identical, map[%d] = map[%d] = %d",
                    v, u, eqmap[u]) ;
                     fflush(stdout) ;
#endif
                     eqmap[v] = eqmap[u] ;
                     break ;
                  }
               }
            }
         }
         if ( eqmap[v] == -1 ) {
/*
            ------------------------------------------------
            v was not found to be indistinguishable from any
            adjacent vertex with lower number, set map value
            ------------------------------------------------
*/
#if MYDEBUG > 1
               fprintf(stdout, "\n    mapping %d to %d", v, nvtx2) ;
               fflush(stdout) ;
#endif
            eqmap[v] = nvtx2++ ;
         }
      }
   }
}
IVfree(mark) ;
IVfree(chksum) ;
#if MYDEBUG > 0
fprintf(stdout, "\n final map") ;
IVfp80(stdout, nvtx, eqmap, 80, &ierr) ;
fflush(stdout) ;
#endif

return(eqmapIV) ; }

/*--------------------------------------------------------------------*/
