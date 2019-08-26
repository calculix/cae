/*  metrics.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   create and return a subtree metric IV object
   input  : vmetricIV -- a metric defined on the vertices
   return : tmetricIV -- a metric defined on the subtrees
  
   created -- 96jun23, cca
   ------------------------------------------------------
*/
IV *
Tree_setSubtreeImetric (
   Tree   *tree,
   IV     *vmetricIV
) {
int   u, v ;
int   *tmetric, *vmetric ;
IV    *tmetricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || tree->n <= 0 
   || vmetricIV == NULL 
   || tree->n != IV_size(vmetricIV) 
   || (vmetric = IV_entries(vmetricIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setSubtreeImetric(%p,%p)"
           "\n bad input\n", tree, vmetricIV) ;
   exit(-1) ;
}
tmetricIV = IV_new() ;
IV_init(tmetricIV, tree->n, NULL) ;
tmetric = IV_entries(tmetricIV) ;
for ( v = Tree_postOTfirst(tree) ; 
      v != -1 ; 
      v = Tree_postOTnext(tree, v) ) {
   tmetric[v] = vmetric[v] ;
   for ( u = tree->fch[v] ; u != -1 ; u = tree->sib[u] ) {
      tmetric[v] += tmetric[u] ;
   }
}
return(tmetricIV) ; }
      
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   create and return a subtree metric DV object
   input  : vmetricDV -- a metric defined on the vertices
   return : tmetricDV -- a metric defined on the subtrees
  
   created -- 96jun23, cca
   ------------------------------------------------------
*/
DV *
Tree_setSubtreeDmetric (
   Tree   *tree,
   DV     *vmetricDV
) {
int      u, v ;
double   *tmetric, *vmetric ;
DV    *tmetricDV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || tree->n <= 0 
   || vmetricDV == NULL 
   || tree->n != DV_size(vmetricDV) 
   || (vmetric = DV_entries(vmetricDV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setSubtreeImetric(%p,%p)"
           "\n bad input\n", tree, vmetricDV) ;
   exit(-1) ;
}
tmetricDV = DV_new() ;
DV_init(tmetricDV, tree->n, NULL) ;
tmetric = DV_entries(tmetricDV) ;
for ( v = Tree_postOTfirst(tree) ; 
      v != -1 ; 
      v = Tree_postOTnext(tree, v) ) {
   tmetric[v] = vmetric[v] ;
   for ( u = tree->fch[v] ; u != -1 ; u = tree->sib[u] ) {
      tmetric[v] += tmetric[u] ;
   }
}
return(tmetricDV) ; }
      
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   create and return a depth metric IV object
   input  : vmetricIV -- a metric defined on the vertices
   output : dmetricIV -- a depth metric defined on the vertices
 
   dmetric[u] = vmetric[u] + dmetric[par[u]] if par[u] != -1
              = vmetric[u]                   if par[u] == -1

   created -- 96jun23, cca
   ------------------------------------------------------------
*/
IV *
Tree_setDepthImetric (
   Tree   *tree,
   IV     *vmetricIV
) {
int   u, v ;
int   *dmetric, *vmetric ;
IV    *dmetricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || tree->n < 1 
   || vmetricIV == NULL 
   || tree->n != IV_size(vmetricIV)
   || (vmetric = IV_entries(vmetricIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setDepthImetric(%p,%p)"
           "\n bad input\n", tree, vmetricIV) ;
   exit(-1) ;
}
dmetricIV = IV_new() ;
IV_init(dmetricIV, tree->n, NULL) ;
dmetric = IV_entries(dmetricIV) ;
for ( u = Tree_preOTfirst(tree) ; 
      u != -1 ; 
      u = Tree_preOTnext(tree, u) ) {
   dmetric[u] = vmetric[u] ;
   if ( (v = tree->par[u]) != -1 ) {
      dmetric[u] += dmetric[v] ;
   }
}
return(dmetricIV) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   create and return a depth metric DV object
   input  : vmetricDV -- a metric defined on the vertices
   output : dmetricDV -- a depth metric defined on the vertices
 
   dmetric[u] = vmetric[u] + dmetric[par[u]] if par[u] != -1
              = vmetric[u]                   if par[u] == -1

   created -- 96jun23, cca
   ------------------------------------------------------------
*/
DV *
Tree_setDepthDmetric (
   Tree   *tree,
   DV     *vmetricDV
) {
int      u, v ;
double   *dmetric, *vmetric ;
DV       *dmetricDV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || tree->n < 1 
   || vmetricDV == NULL 
   || tree->n != DV_size(vmetricDV)
   || (vmetric = DV_entries(vmetricDV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setDepthDmetric(%p,%p)"
           "\n bad input\n", tree, vmetricDV) ;
   exit(-1) ;
}
dmetricDV = DV_new() ;
DV_init(dmetricDV, tree->n, NULL) ;
dmetric = DV_entries(dmetricDV) ;
for ( u = Tree_preOTfirst(tree) ; 
      u != -1 ; 
      u = Tree_preOTnext(tree, u) ) {
   dmetric[u] = vmetric[u] ;
   if ( (v = tree->par[u]) != -1 ) {
      dmetric[u] += dmetric[v] ;
   }
}
return(dmetricDV) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   create and return a height metric IV object
   input  : vmetricIV -- a metric defined on the vertices
   output : dmetricIV -- a depth metric defined on the vertices
 
   hmetric[v] = vmetric[v] + max{p(u) = v} hmetric[u] if fch[v] != -1
              = vmetric[v]                            if fch[v] == -1

   created -- 96jun23, cca
   ------------------------------------------------------------------
*/
IV *
Tree_setHeightImetric (
   Tree   *tree,
   IV     *vmetricIV
) {
int   u, v, val ;
int   *hmetric, *vmetric ;
IV    *hmetricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || tree->n < 1 
   || vmetricIV == NULL 
   || tree->n != IV_size(vmetricIV)
   || (vmetric = IV_entries(vmetricIV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setHeightImetric(%p,%p)"
           "\n bad input\n", tree, vmetricIV) ;
   if ( tree != NULL ) {
      Tree_writeForHumanEye(tree, stderr) ;
   }
   if ( vmetricIV != NULL ) {
      IV_writeForHumanEye(vmetricIV, stderr) ;
   }
   exit(-1) ;
}
hmetricIV = IV_new() ; 
IV_init(hmetricIV, tree->n, NULL) ; 
hmetric = IV_entries(hmetricIV) ;
for ( v = Tree_postOTfirst(tree) ; 
      v != -1 ; 
      v = Tree_postOTnext(tree, v) ) {
   for ( u = tree->fch[v], val = 0 ; u != -1 ; u = tree->sib[u] ) {
      if ( val < hmetric[u] ) {
         val = hmetric[u] ;
      }
   }
   hmetric[v] = val + vmetric[v] ;
}
return(hmetricIV) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   create and return a height metric DV object
   input  : vmetricDV -- a metric defined on the vertices
   output : dmetricDV -- a depth metric defined on the vertices
 
   hmetric[v] = vmetric[v] + max{p(u) = v} hmetric[u] if fch[v] != -1
              = vmetric[v]                            if fch[v] == -1

   created -- 96jun23, cca
   ------------------------------------------------------------------
*/
DV *
Tree_setHeightDmetric (
   Tree   *tree,
   DV     *vmetricDV
) {
int      u, v, val ;
double   *hmetric, *vmetric ;
DV       *hmetricDV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || tree->n < 1 
   || vmetricDV == NULL 
   || tree->n != DV_size(vmetricDV)
   || (vmetric = DV_entries(vmetricDV)) == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setHeightDmetric(%p,%p)"
           "\n bad input\n", tree, vmetricDV) ;
   exit(-1) ;
}
hmetricDV = DV_new() ; 
DV_init(hmetricDV, tree->n, NULL) ; 
hmetric = DV_entries(hmetricDV) ;
for ( v = Tree_postOTfirst(tree) ; 
      v != -1 ; 
      v = Tree_postOTnext(tree, v) ) {
   for ( u = tree->fch[v], val = 0 ; u != -1 ; u = tree->sib[u] ) {
      if ( val < hmetric[u] ) {
         val = hmetric[u] ;
      }
   }
   hmetric[v] = val + vmetric[v] ;
}
return(hmetricDV) ; }

/*--------------------------------------------------------------------*/
