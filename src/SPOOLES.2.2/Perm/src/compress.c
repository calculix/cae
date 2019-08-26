/*  compress.c  */

#include "../Perm.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   given a permutation and a vector to map vertices 
   into compressed vertices, create and return a 
   permutation object for the compressed vertices.

   created -- 96may02, cca
   ------------------------------------------------
*/
Perm *
Perm_compress (
   Perm   *perm,
   IV     *eqmapIV
) {
int    n, N, v, vcomp, vnew ;
int    *eqmap, *head, *link, *newToOld, *oldToNew, *vals ; 
Perm   *perm2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  perm == NULL 
   || (n = perm->size) <= 0
   || eqmapIV == NULL 
   || n != IV_size(eqmapIV)
   || (eqmap = IV_entries(eqmapIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_compress(%p,%p)"
           "\n bad input\n", perm, eqmapIV) ;
   if ( perm != NULL ) {
      Perm_writeStats(perm, stderr) ;
   }
   if ( eqmapIV != NULL ) {
      IV_writeStats(eqmapIV, stderr) ;
   }
   exit(-1) ;
}
n = perm->size ;
if ( (oldToNew = perm->oldToNew) == NULL ) {
   Perm_fillOldToNew(perm) ;
   oldToNew = perm->oldToNew ;
}
if ( (newToOld = perm->newToOld) == NULL ) {
   Perm_fillNewToOld(perm) ;
   newToOld = perm->newToOld ;
}
/*
   ---------------------------------
   create the new permutation object
   ---------------------------------
*/
N = 1 + IVmax(n, eqmap, &v) ;
perm2 = Perm_new() ;
Perm_initWithTypeAndSize(perm2, 3, N) ;
/*
   --------------------------------------------
   get the head/link structure for the vertices
   --------------------------------------------
*/
head = IVinit(N, -1) ;
link = IVinit(n, -1) ;
for ( v = 0 ; v < n ; v++ ) {
   vcomp = eqmap[v] ;
   link[v] = head[vcomp] ;
   head[vcomp] = v ;
}
/*
   ---------------------------
   get the two vectors to sort
   ---------------------------
*/
IVramp(N, perm2->newToOld, 0, 1) ;
vals = IVinit(N, -1) ;
for ( vcomp = 0 ; vcomp < N ; vcomp++ ) {
   v = head[vcomp] ;
   vnew = perm->oldToNew[v] ;
   for ( v = link[v] ; v != -1 ; v = link[v] ) {
      if ( vnew > perm->oldToNew[v] ) {
         vnew = perm->oldToNew[v] ;
      }
   }
   vals[vcomp] = vnew ;
}
IV2qsortUp(N, vals, perm2->newToOld) ;
for ( vcomp = 0 ; vcomp < N ; vcomp++ ) {
   perm2->oldToNew[perm2->newToOld[vcomp]] = vcomp ;
}
/*
   ---------------------
   free the working data
   ---------------------
*/
IVfree(head) ;
IVfree(link) ;
IVfree(vals) ;

return(perm2) ; }

/*--------------------------------------------------------------------*/
