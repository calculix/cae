/*  basics.c  */

#include "../Tree.h"

#define MYTRACE 0
#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- create and return a new Tree object

   created -- 95nov15, cca
   -----------------------------------------------
*/
Tree *
Tree_new ( 
   void
) {
Tree   *tree ;

#if MYTRACE > 0
fprintf(stdout, "\n just inside Tree_new()") ;
fflush(stdout) ;
#endif

ALLOCATE(tree, struct _Tree, 1) ;

Tree_setDefaultFields(tree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Tree_new()") ;
fflush(stdout) ;
#endif

return(tree) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- set the default fields for the Tree object

   created -- 95nov15, cca
   ------------------------------------------------------
*/
void
Tree_setDefaultFields (
   Tree   *tree
) {

#if MYTRACE > 0
fprintf(stdout, "\n just inside Tree_setDefaultFields(%)", g) ;
fflush(stdout) ;
#endif

if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_setDefaultFields(%p)"
           "\n tree is NULL\n", tree) ;
   exit(-1) ;
}
tree->n    =   0  ;
tree->root =  -1  ;
tree->par  = NULL ;
tree->fch  = NULL ;
tree->sib  = NULL ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Tree_setDefaultFields(%)", tree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 95nov15, cca
   --------------------------------
*/
void
Tree_clearData (
   Tree   *tree
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside Tree_clearData(%)", tree) ;
fflush(stdout) ;
#endif

if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_clearData(%p)"
           "\n tree is NULL\n", tree) ;
   exit(-1) ;
}

if ( tree->par != NULL ) {
   IVfree(tree->par) ;
}
if ( tree->fch != NULL ) {
   IVfree(tree->fch) ;
}
if ( tree->sib != NULL ) {
   IVfree(tree->sib) ;
}
Tree_setDefaultFields(tree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Tree_clearData(%)", tree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- free the Tree object

   created -- 95nov15, cca
   --------------------------------
*/
void
Tree_free (
   Tree   *tree
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside Tree_free(%)", tree) ;
fflush(stdout) ;
#endif

if ( tree == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_free(%p)"
           "\n tree is NULL\n", tree) ;
   exit(-1) ;
}

Tree_clearData(tree) ;

FREE(tree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving Tree_free(%)", tree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
