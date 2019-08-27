/*  util.c  */

#include "../DSTree.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 96mar10, cca
   ----------------------------------------------
*/
int
DSTree_sizeOf (
   DSTree   *dstree
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_sizeOf(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
nbytes = sizeof(struct _DSTree) ;
if ( dstree->tree != NULL ) {
   nbytes += Tree_sizeOf(dstree->tree) ;
}
if ( dstree->mapIV != NULL ) {
   nbytes += IV_sizeOf(dstree->mapIV) ;
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   renumber the fronts by a post-order traversal
   
   created -- 96apr13, cca
   ---------------------------------------------
*/
void
DSTree_renumberViaPostOT (
   DSTree * dstree 
) {
int    count, J, K, n, nvtx, v ;
int    *map, *oldmap, *temp ;
IV     *mapIV ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL 
   || (tree = dstree->tree) == NULL
   || (n = tree->n) <= 0
   || (mapIV  = dstree->mapIV) == NULL 
   || (nvtx   = IV_size(mapIV)) <= 0 
   || (oldmap = IV_entries(mapIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_renumberViaPostOT(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
/*
   ---------------------------------
   renumber the Tree object but
   preserve the post-order traversal
   ---------------------------------
*/
map   = IVinit(n, -1) ;
temp  = IVinit(n, -1) ;
count = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   map[J] = count++ ;
}
for ( J = 0 ; J < n ; J++ ) {
   if ( (K = tree->par[J]) != -1 ) {
      temp[map[J]] = map[K] ;
   } else {
      temp[map[J]] = -1 ;
   }
}
IVcopy(n, tree->par, temp) ;
for ( J = 0 ; J < n ; J++ ) {
   if ( (K = tree->fch[J]) != -1 ) {
      temp[map[J]] = map[K] ;
   } else {
      temp[map[J]] = -1 ;
   }
}
IVcopy(n, tree->fch, temp) ;
for ( J = 0 ; J < n ; J++ ) {
   if ( (K = tree->sib[J]) != -1 ) {
      temp[map[J]] = map[K] ;
   } else {
      temp[map[J]] = -1 ;
   }
}
IVcopy(n, tree->sib, temp) ;
if ( tree->root != -1 ) {
   tree->root = map[tree->root] ;
}
/*
   -----------------------------
   remap the vertex to front map
   -----------------------------
*/
for ( v = 0 ; v < nvtx ; v++ ) {
   J = oldmap[v] ;
   if ( 0 <= J && J < n ) {
      oldmap[v] = map[J] ;
   }
}
IVfree(map) ;
IVfree(temp) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- return the weight of the vertices in the domains

   created -- 97jun21, cca
   -----------------------------------------------------------
*/
int
DSTree_domainWeight (
   DSTree   *dstree,
   int      vwghts[]
) {
int    domwght, ireg, nvtx, v ;
int    *fch, *map ;
IV     *mapIV ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_domainWeight(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
tree  = DSTree_tree(dstree) ;
mapIV = DSTree_mapIV(dstree) ;
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
fch = tree->fch ;
if ( vwghts != NULL ) {
   for ( v = domwght = 0 ; v < nvtx ; v++ ) {
      ireg = map[v] ;
      if ( fch[ireg] == -1 ) {
         domwght += vwghts[v] ;
      }
   }
} else {
   for ( v = domwght = 0 ; v < nvtx ; v++ ) {
      ireg = map[v] ;
      if ( fch[ireg] == -1 ) {
         domwght++ ;
      }
   }
}
return(domwght) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- return the weight of the vertices in the separators

   created -- 97jun21, cca
   --------------------------------------------------------------
*/
int
DSTree_separatorWeight (
   DSTree   *dstree,
   int      vwghts[]
) {
int    ireg, nvtx, sepwght, v ;
int    *fch, *map ;
IV     *mapIV ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_separatorWeight(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
tree  = DSTree_tree(dstree) ;
mapIV = DSTree_mapIV(dstree) ;
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
fch = tree->fch ;
if ( vwghts != NULL ) {
   for ( v = sepwght = 0 ; v < nvtx ; v++ ) {
      ireg = map[v] ;
      if ( fch[ireg] != -1 ) {
         sepwght += vwghts[v] ;
      }
   }
} else {
   for ( v = sepwght = 0 ; v < nvtx ; v++ ) {
      ireg = map[v] ;
      if ( fch[ireg] != -1 ) {
         sepwght++ ;
      }
   }
}
return(sepwght) ; }

/*--------------------------------------------------------------------*/
