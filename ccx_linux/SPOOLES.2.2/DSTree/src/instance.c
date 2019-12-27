/*  instance.c  */

#include "../DSTree.h"
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   return a pointer to the internal Tree object

   created -- 97jun21, cca
   --------------------------------------------
*/
Tree *
DSTree_tree (
   DSTree   *dstree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_tree(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
return(dstree->tree) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   return a pointer to the map IV object

   created -- 97jun21, cca
   -------------------------------------
*/
IV *
DSTree_mapIV (
   DSTree   *dstree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_mapIV(%p)"
           "\n bad input\n", dstree) ;
   exit(-1) ;
}
return(dstree->mapIV) ; }

/*--------------------------------------------------------------------*/
