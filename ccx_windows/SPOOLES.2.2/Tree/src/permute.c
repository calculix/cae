/*  permute.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   return a permuted tree

   created -- 96jan04, cca
   -----------------------
*/
Tree *
Tree_permute (
   Tree   *tree,
   int    newToOld[],
   int    oldToNew[]
) {
int    n, u_old, v_new, v_old, w_old ;
Tree   *tree2 ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || (n = tree->n) <= 0 
   || newToOld == NULL || oldToNew == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_permute(%p,%p,%p)"
           "\n bad input\n", tree, newToOld, oldToNew) ;
   exit(-1) ;
}
/*
   -----------------
   create a new tree
   -----------------
*/
tree2 = Tree_new() ;
Tree_init1(tree2, n) ;
/*
   ---------------
   fill the fields
   ---------------
*/
for ( v_new = 0 ; v_new < n ; v_new++ ) {
   v_old = newToOld[v_new] ;
   if ( (w_old = tree->par[v_old]) != -1 ) {
      tree2->par[v_new] = oldToNew[w_old] ;
   }
   if ( (u_old = tree->fch[v_old]) != -1 ) {
      tree2->fch[v_new] = oldToNew[u_old] ;
   }
   if ( (u_old = tree->sib[v_old]) != -1 ) {
      tree2->sib[v_new] = oldToNew[u_old] ;
   }
}
if ( tree->root != -1 ) {
   tree2->root = oldToNew[tree->root] ;
}

return(tree2) ; }

/*--------------------------------------------------------------------*/
