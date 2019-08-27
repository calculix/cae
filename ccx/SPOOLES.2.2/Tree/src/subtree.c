/*  subtree.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to initialize subtree with the subtree 
              of tree using nodes in nodeidsIV

   return values ---
      1 -- normal return
     -1 -- subtree is NULL
     -2 -- nodeidsIV is NULL
     -3 -- tree is NULL
     -4 -- nodeidsIV is invalid

   created -- 98oct15, cca
   -------------------------------------------------
*/
int
Tree_initFromSubtree (
   Tree   *subtree,
   IV     *nodeidsIV,
   Tree   *tree
) {
int   J, Jsub, K, Ksub, nnodeInSubtree, nnodeInTree ;
int   *localmap, *par, *parsub, *subtreeNodes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( subtree == NULL ) {
   fprintf(stderr, "\n\n error in Tree_initFromSubtree()"
           "\n subtree is NULL\n") ;
   return(-1) ;
}
if ( nodeidsIV == NULL ) {
   fprintf(stderr, "\n\n error in Tree_initFromSubtree()"
           "\n nodeidsIV is NULL\n") ;
   return(-2) ;
}
if ( tree == NULL ) {
   fprintf(stderr, "\n\n error in Tree_initFromSubtree()"
           "\n tree is NULL\n") ;
   return(-3) ;
}
nnodeInTree = Tree_nnodes(tree) ;
IV_sizeAndEntries(nodeidsIV, &nnodeInSubtree, &subtreeNodes) ;
if ( nnodeInSubtree < 0 || nnodeInSubtree >= nnodeInTree ) {
   fprintf(stderr, "\n\n error in Tree_initFromSubtree()"
           "\n nnodeInTree = %d, nnodeInSubtree = %d\n",
           nnodeInTree, nnodeInSubtree) ;
   return(-4) ;
}
for ( Jsub = 0 ; Jsub < nnodeInSubtree ; Jsub++ ) {
   J = subtreeNodes[Jsub] ;
   if ( J < 0 || J >= nnodeInTree ) {
      fprintf(stderr, "\n\n error in Tree_initFromSubtree()"
              "\n nnodeInTree = %d, subtreeNodes[%d] = %d\n",
              nnodeInTree, Jsub, subtreeNodes[Jsub]) ;
      return(-4) ;
   }
}
par = Tree_par(tree) ;
/*
   ----------------------------
   create a global-to-local map
   ----------------------------
*/
localmap = IVinit(nnodeInTree, -1) ;
for ( Jsub = 0 ; Jsub < nnodeInSubtree ; Jsub++ ) {
   localmap[subtreeNodes[Jsub]] = Jsub ;
}
/*
   ----------------------
   initialize the subtree
   ----------------------
*/
Tree_init1(subtree, nnodeInSubtree) ;
parsub = Tree_par(subtree) ;
for ( Jsub = 0 ; Jsub < nnodeInSubtree ; Jsub++ ) {
   J = subtreeNodes[Jsub] ;
   if ( (K = par[J]) != -1 && (Ksub = localmap[K]) != -1 ) {
      parsub[Jsub] = Ksub ;
   }
}
Tree_setFchSibRoot(subtree) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(localmap) ;

return(1) ; }

/*--------------------------------------------------------------------*/
