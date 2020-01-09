/*  justify.c  */

#include "../ETree.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   left-justify a tree by subtree size
   children are linked in ascending order of their subtree size

   created -- 96jan11, cca
   ------------------------------------------------------------
*/
void
ETree_leftJustify (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_leftJustify(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
Tree_leftJustify(etree->tree) ;
 
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   left-justify a etree by a metric
   children are linked in ascending order of their metric

   created -- 96jan11, cca
   ------------------------------------------------------
*/
void
ETree_leftJustifyI (
   ETree   *etree,
   IV      *metricIV
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || etree->nfront <= 0
   || etree->nvtx <= 0
   || metricIV == NULL 
   || IV_size(metricIV) != etree->nfront
   || IV_entries(metricIV) == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_leftJustifyI(%p,%p)"
           "\n bad input\n", etree, metricIV) ;
   exit(-1) ;
}
Tree_leftJustifyI(etree->tree, metricIV) ;
 
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   left-justify a etree by a metric
   children are linked in ascending order of their metric

   created -- 96jan11, cca
   ------------------------------------------------------
*/
void
ETree_leftJustifyD (
   ETree     *etree,
   DV        *metricDV
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || etree->nfront <= 0
   || etree->nvtx <= 0
   || metricDV == NULL 
   || DV_size(metricDV) != etree->nfront
   || DV_entries(metricDV) == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_leftJustifyD(%p,%p)"
           "\n bad input\n", etree, metricDV) ;
   exit(-1) ;
}
Tree_leftJustifyD(etree->tree, metricDV) ;
 
return ; }

/*--------------------------------------------------------------------*/
