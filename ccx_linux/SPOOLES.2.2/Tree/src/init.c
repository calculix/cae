/*  init.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 95nov15, cca
   -----------------------
*/
void
Tree_init1 (
   Tree   *tree,
   int    size
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || size < 0 ) {
   fprintf(stderr, "\n fatal error in Tree_init1(%p,%d)"
           "\n bad input\n", tree, size) ;
   exit(-1) ;
}
/*
   -----------------------
   clear any previous data
   -----------------------
*/
Tree_clearData(tree) ;
/*
   -----------------------------------------------
   set size field and initialize the three vectors
   -----------------------------------------------
*/
tree->n   = size ;
if ( size > 0 ) {
  tree->par = IVinit(size, -1) ;
  tree->fch = IVinit(size, -1) ;
  tree->sib = IVinit(size, -1) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   initialize given a parent vector

   created -- 95nov15, cca
   --------------------------------
*/
void
Tree_init2 (
   Tree   *tree,
   int    size,
   int    par[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || size <= 0 || par == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_init2(%p,%d,%p)"
           "\n bad input\n", tree, size, par) ;
   exit(-1) ;
}
/*
   ----------------------------
   first use simple initializer
   ----------------------------
*/
Tree_init1(tree, size) ;
/*
   ------------------
   copy parent vector
   ------------------
*/
IVcopy(size, tree->par, par) ;
/*
   -------------------------
   set fch[], sib[] and root
   -------------------------
*/
Tree_setFchSibRoot(tree) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   initialize given the tree vectors

   created -- 95nov15, cca
   ---------------------------------
*/
void
Tree_init3 (
   Tree   *tree,
   int    size,
   int    par[],
   int    fch[],
   int    sib[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || size <= 0 
   || par == NULL || fch == NULL || sib == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_init3(%p,%d,%p,%p,%p)"
           "\n bad input\n", tree, size, par, fch, sib) ;
   exit(-1) ;
}
/*
   ----------------------------
   first use simple initializer
   ----------------------------
*/
Tree_init1(tree, size) ;
/*
   ----------------------
   copy the three vectors
   ----------------------
*/
IVcopy(size, tree->par, par) ;
IVcopy(size, tree->fch, fch) ;
IVcopy(size, tree->sib, sib) ;
/*
   --------
   set root
   --------
*/
Tree_setRoot(tree) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   set the fch[], sib[] and root fields

   created -- 95nov15, cca
   ------------------------------------
*/
void
Tree_setFchSibRoot (
   Tree   *tree
) {
int   n, root, u, v ;
int   *fch, *par, *sib ;
/*
   ---------------
   check the input
   ---------------
if (  tree == NULL || (n = tree->n) < 1 ) {
*/
if (  tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setFchSibRoot(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
if (  (n = tree->n) < 1 ) {
   return ;
}
par = tree->par ;
fch = tree->fch ;
sib = tree->sib ;
/*
   ---------------------
   initialize the fields
   ---------------------
*/
IVfill(n, tree->fch, -1) ;
IVfill(n, tree->sib, -1) ;
root = -1 ;
/*
   --------------
   set the fields
   --------------
for ( u = 0 ; u < n ; u++ ) {
*/
for ( u = n - 1 ; u >= 0 ; u-- ) {
   if ( (v = par[u]) != -1 ) {
      sib[u] = fch[v] ;
      fch[v] =    u   ;
   } else {
      sib[u] = root ;
      root   =   u  ;
   }
}
tree->root = root ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the root field

   created -- 95nov15, cca
   -----------------------
*/
void
Tree_setRoot (
   Tree   *tree
) {
int   n, root, u ;
int   *par, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if (  tree == NULL || (n = tree->n) < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_setRoot(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
n    = tree->n   ;
par  = tree->par ;
sib  = tree->sib ;
root = -1 ;
/*
   --------------
   set the fields
   --------------
*/
for ( u = 0 ; u < n ; u++ ) {
   if ( par[u] == -1 ) {
      sib[u] = root ;
      root   =   u  ;
   }
}
tree->root = root ;

return ; }

/*--------------------------------------------------------------------*/
