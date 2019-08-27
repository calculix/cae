/*  compress.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   create and return an IV object that contains
   the map from vertices to fundamental chains.

   return value -- # of fundamental chains

   created -- 96jun23, cca
   -------------------------------------------
*/
IV *
Tree_fundChainMap (
   Tree   *tree
) {
int   nfc, u, v ;
int   *map ;
IV    *mapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n <= 0 ) {
   fprintf(stderr, "\n fatal error in Tree_fundChainMap(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
mapIV = IV_new() ;
IV_init(mapIV, tree->n, NULL) ;
map = IV_entries(mapIV) ;
for ( v = Tree_postOTfirst(tree), nfc = 0 ;
      v != -1 ;
      v = Tree_postOTnext(tree, v) ) {
   if ( (u = tree->fch[v]) == -1 || tree->sib[u] != -1 ) {
/*
      --------------------
      v starts a new chain
      --------------------
*/
      map[v] = nfc++ ;
   } else {
/*
      -----------------------------------------------
      v belongs in the same chain as its only child u
      -----------------------------------------------
*/
      map[v] = map[u] ;
   }
}
return(mapIV) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   compress a tree based on a map from old vertices to new vertices.
   the restriction on the map is that the set {u | map[u] = U} must
   be connected for all U.

   created  -- 95nov15, cca
   modified -- 96jan04, cca
      bug fixed, N computed incorrectly
   modified -- 96jun23, cca
      in calling sequence, int map[] converted to IV *mapIV 
   -----------------------------------------------------------------
*/
Tree *
Tree_compress (
   Tree   *tree,
   IV     *mapIV
) {
int    n, N, u, U, v, V ;
int    *head, *link, *map ;
Tree   *tree2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL 
   || (n = tree->n) <= 0 
   || mapIV == NULL 
   || n != IV_size(mapIV)
   || (map = IV_entries(mapIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_compress(%p,%p)"
           "\n bad input\n", tree, mapIV) ;
   exit(-1) ;
}
/*
   -----------------------
   initialize the new tree
   -----------------------
*/
N = 1 + IV_max(mapIV) ;
tree2 = Tree_new() ;
Tree_init1(tree2, N) ;
/*
   -----------------------------------------------------------
   get the head/link construct to map old nodes into new nodes
   -----------------------------------------------------------
*/
head = IVinit(N, -1) ;
link = IVinit(n, -1) ;
for ( v = 0 ; v < n ; v++ ) {
   if ( (V = map[v]) < 0 || V >= N ) {
      fprintf(stderr, "\n fatal error in Tree_compress(%p,%p)"
              "\n map[%d] = %d, N = %d\n", tree, map, v, V, N) ;
      exit(-1) ;
   }
   link[v] = head[V] ;
   head[V] =    v    ;
}
/*
   ---------------------
   fill the tree vectors
   ---------------------
*/
for ( U = 0 ; U < N ; U++ ) {
   for ( u = head[U] ; u != -1 ; u = link[u] ) {
      if ( (v = tree->par[u]) == -1 ) {
         tree2->par[U] = -1 ;
         break ;
      } else if ( (V = map[v]) != U ) {
         tree2->par[U] = V ;
         break ;
      }
   }
}
Tree_setFchSibRoot(tree2) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(head) ;
IVfree(link) ;
 
return(tree2) ; }

/*--------------------------------------------------------------------*/
