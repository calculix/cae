/*  util.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   return the first vertex in a post-order traversal
   
   created -- 95nov15, cca
   -------------------------------------------------
*/
int
Tree_postOTfirst (
   Tree   *tree
) {
int   v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_postOTfirst(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
/*
   ----------------------
   find the leftmost leaf
   ----------------------
*/
if ( (v = tree->root) != -1 ) {
   while ( tree->fch[v] != -1 ) {
      v = tree->fch[v] ;
   }
}
return(v) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   return the vertex that follows v in a post-order traversal
   ----------------------------------------------------------
*/
int
Tree_postOTnext (
   Tree   *tree,
   int    v
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 || v < 0 || v >= tree->n ) {
   fprintf(stderr, "\n fatal error in Tree_postOTnext(%p,%d)"
           "\n bad input\n", tree, v) ;
   exit(-1) ;
}
/*
   ---------------------------------------------------
   find leftmost leaf of sibling (if exists) or parent
   ---------------------------------------------------
*/
if ( tree->sib[v] != -1 ) {
   v = tree->sib[v] ;
   while ( tree->fch[v] != -1 ) {
      v = tree->fch[v] ;
   }
} else {
   v = tree->par[v] ;
}
return(v) ; }
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   return the first vertex in a pre-order traversal
   
   created -- 95nov15, cca
   ------------------------------------------------
*/
int
Tree_preOTfirst (
   Tree   *tree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_preOTfirst(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
return(tree->root) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   return the vertex that follows v in a pre-order traversal
   
   created -- 95nov15, cca
   ---------------------------------------------------------
*/
int
Tree_preOTnext (
   Tree   *tree,
   int    v
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 || v < 0 || v >= tree->n ) {
   fprintf(stderr, "\n fatal error in Tree_preOTnext(%p,%d)"
           "\n bad input\n", tree, v) ;
   exit(-1) ;
}
/*
   -------------------------------------
   find the next vertex in the traversal
   -------------------------------------
*/
if ( tree->fch[v] != -1 ) {
   v = tree->fch[v] ;
} else {
   while ( tree->sib[v] == -1 && tree->par[v] != -1 ) {
      v = tree->par[v] ;
   }
   v = tree->sib[v] ;
}
return(v) ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   return the number of leaves in the tree

   created -- 95nov15, cca
   ---------------------------------------
*/
int
Tree_nleaves (
   Tree   *tree
) {
int   nleaf, v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_nleaves(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
 
nleaf = 0 ;
for ( v = 0 ; v < tree->n ; v++ ) {
   if ( tree->fch[v] == -1 ) {
      nleaf++ ;
   }
}
return(nleaf) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   return the number of roots of the tree (forest)

   created -- 95nov15, cca
   -----------------------------------------------
*/
int
Tree_nroots (
   Tree   *tree
) {
int   nroot, v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_nroots(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
 
nroot = 0 ;
for ( v = 0 ; v < tree->n ;v++ ) {
   if ( tree->par[v] == -1 ) {
      nroot++ ;
   }
}
return(nroot) ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   return the number of children of a vertex

   created  -- 95nov15, cca
   modified -- 96jan07, cca
      v checked to be valid
   -----------------------------------------
*/
int
Tree_nchild (
   Tree   *tree,
   int    v
) {
int   nchild, w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_nchild(%p,%d)"
           "\n bad input\n", tree, v) ;
   exit(-1) ;
}
if ( v < 0 || v >= tree->n ) {
   fprintf(stderr, "\n fatal error in Tree_nchild(%p,%d)"
           "\n v = %d, size = %d\n", tree, v, v, tree->n) ;
   exit(-1) ;
}
nchild = 0 ;
for ( w = tree->fch[v] ; w != -1 ; w = tree->sib[w] ) {
   nchild++ ;
}
return(nchild) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   this method returns an IV object that holds 
   the number of children for the tree nodes.

   created -- 96dec18, cca
   -------------------------------------------
*/
IV *
Tree_nchildIV (
   Tree   *tree
) {
int   n, v, w ;
int   *nchild, *par ;
IV    *nchildIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || (n = tree->n) < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_nchildIV(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
nchildIV = IV_new() ;
IV_init(nchildIV, n, NULL) ;
IV_fill(nchildIV, 0) ;
par = tree->par ;
nchild = IV_entries(nchildIV) ;
for ( v = 0 ; v < n ; v++ ) {
   if ( (w = par[v]) != -1 ) {
      nchild[w]++ ;
   }
}

return(nchildIV) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   return the height of the tree

   created -- 96aug23, cca
   -----------------------------
*/
int
Tree_height (
   Tree   *tree
) {
int   u, v, vheight ;
int   *heights ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || tree->n < 1 ) {
   fprintf(stderr, "\n fatal error in Tree_height(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
heights = IVinit(tree->n, 1) ;
for ( v = Tree_postOTfirst(tree) ;
      v != -1 ;
      v = Tree_postOTnext(tree, v) ) {
   if ( (u = tree->fch[v]) == -1 ) {
      vheight = 1 ;
   } else {
      vheight = heights[u] ;
      for ( u = tree->sib[u] ; u != -1 ; u = tree->sib[u] ) {
         if ( vheight < heights[u] ) {
            vheight = heights[u] ;
         }
      }
      vheight++ ;
   }
   heights[v] = vheight ;
}
v = tree->root ;
vheight = heights[v] ;
for ( v = tree->sib[v] ; v != -1 ; v = tree->sib[v] ) {
   if ( vheight < heights[v] ) {
      vheight = heights[v] ;
   }
}
IVfree(heights) ;

return(vheight) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   return the maximum number of children of any node in the tree

   created -- 96sep05, cca
   -------------------------------------------------------------
*/
int
Tree_maxNchild (
   Tree   *tree
) {
int   maxnchild, n, nchild, u, v ;
int   *fch, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_maxNchild(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
n   = tree->n   ;
fch = tree->fch ;
sib = tree->sib ;
maxnchild = 0 ;
for ( v = 0 ; v < n ; v++ ) {
   for ( u = fch[v], nchild = 0 ; u != -1 ; u = sib[u] ) {
      nchild++ ;
   }
   if ( maxnchild < nchild ) {
      maxnchild = nchild ;
   }
}
return(maxnchild) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   return the number of bytes used by the object
   ---------------------------------------------
*/
int
Tree_sizeOf (
   Tree   *tree
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_sizeOf(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
 
nbytes = sizeof(struct _Tree) + 3*tree->n*sizeof(int) ;
 
return(nbytes) ; }

/*--------------------------------------------------------------------*/
