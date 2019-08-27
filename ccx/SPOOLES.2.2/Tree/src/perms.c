/*  perms.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   fill the new-to-old permutation vector

   created -- 95nov15, cca
   --------------------------------------
*/
void
Tree_fillNewToOldPerm (
   Tree   *tree,
   int    newToOld[]
) {
int   i, v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 || newToOld == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_fillNewToOldPerm(%p,%p)"
           "\n bad input\n", tree, newToOld) ;
   exit(-1) ;
}
/*
   -----------------------------------------------
   post-order traversal to fill permutation vector
   -----------------------------------------------
*/
for ( v = Tree_postOTfirst(tree), i = 0 ;
      v != -1 ;
      v = Tree_postOTnext(tree, v) ) {
   newToOld[i++] = v ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   fill the old-to-new permutation vector

   created -- 95nov15, cca
   --------------------------------------
*/
void
Tree_fillOldToNewPerm (
   Tree   *tree,
   int    oldToNew[]
) {
int   i, v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 || oldToNew == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_fillOldToNewPerm(%p,%p)"
           "\n bad input\n", tree, oldToNew) ;
   exit(-1) ;
}
/*
   -----------------------------------------------
   post-order traversal to fill permutation vector
   -----------------------------------------------
*/
for ( v = Tree_postOTfirst(tree), i = 0 ;
      v != -1 ;
      v = Tree_postOTnext(tree, v) ) {
   oldToNew[v] = i++ ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   fill the new-to-old and old-to-new permutation vectors

   created -- 95nov15, cca
   ------------------------------------------------------
*/
void
Tree_fillBothPerms (
   Tree   *tree,
   int    newToOld[],
   int    oldToNew[]
) {
int   i, v ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || tree->n < 1 
   || newToOld == NULL || oldToNew == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_fillBothPerms(%p,%p,%p)"
           "\n bad input\n", tree, newToOld, oldToNew) ;
   exit(-1) ;
}
/*
   ------------------------------------------------
   post-order traversal to fill permutation vectors
   ------------------------------------------------
*/
for ( v = Tree_postOTfirst(tree), i = 0 ;
      v != -1 ;
      v = Tree_postOTnext(tree, v) ) {
   newToOld[i] =  v  ;
   oldToNew[v] = i++ ;
}
return ; }

/*--------------------------------------------------------------------*/
