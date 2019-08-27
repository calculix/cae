/*  instance.c  */

#include "../Tree.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- return the number of nodes in the tree

   created -- 98jun12, cca
   -------------------------------------------------
*/
int
Tree_nnodes (
   Tree   *tree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_nnodes(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
return(tree->n) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- return the root of the tree

   created -- 98jun12, cca
   --------------------------------------
*/
int
Tree_root (
   Tree   *tree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_root(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
return(tree->root) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- return a pointer to the parent vector

   created -- 98jun12, cca
   ------------------------------------------------
*/
int *
Tree_par (
   Tree   *tree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_par(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
return(tree->par) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- return a pointer to the first child vector

   created -- 98jun12, cca
   -----------------------------------------------------
*/
int *
Tree_fch (
   Tree   *tree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_fch(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
return(tree->fch) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- return a pointer to the sibling vector

   created -- 98jun12, cca
   -------------------------------------------------
*/
int *
Tree_sib (
   Tree   *tree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_sib(%p)"
           "\n bad input\n", tree) ;
   exit(-1) ;
}
return(tree->sib) ; }

/*--------------------------------------------------------------------*/
