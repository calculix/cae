/*  init.c  */

#include "../DSTree.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   initialize the object given the number of nodes

   created -- 96mar10, cca
   -----------------------------------------------
*/
void
DSTree_init1 (
   DSTree   *dstree,
   int      ndomsep,
   int      nvtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || ndomsep <= 0 ) {
   fprintf(stderr, "\n fatal error in DSTree_init1(%p,%d,%d)"
           "\n bad input\n", dstree, ndomsep, nvtx) ;
   exit(-1) ;
}
DSTree_clearData(dstree) ;
dstree->tree = Tree_new() ;
Tree_init1(dstree->tree, ndomsep) ;
dstree->mapIV = IV_new() ;
IV_init(dstree->mapIV, nvtx, NULL) ;
IV_fill(dstree->mapIV, -1) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   initialize the object given a Tree object and a map IV object

   created  -- 96mar10, cca
   -------------------------------------------------------------
*/
void
DSTree_init2 (
   DSTree   *dstree,
   Tree     *tree,
   IV       *mapIV
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || tree == NULL || tree->n < 1 
   || mapIV == NULL || IV_size(mapIV) < 1 ) {
   fprintf(stderr, "\n fatal error in DSTree_init2(%p,%p,%p)"
           "\n bad input\n", dstree, tree, mapIV) ;
   exit(-1) ;
}
DSTree_clearData(dstree) ;
dstree->tree  = tree  ;
dstree->mapIV = mapIV ;

return ; }

/*--------------------------------------------------------------------*/
