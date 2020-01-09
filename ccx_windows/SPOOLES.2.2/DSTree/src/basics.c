/*  basics.c  */

#include "../DSTree.h"

#define MYTRACE 0
#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- create and return a new DSTree object

   created -- 96mar10, cca
   -----------------------------------------------
*/
DSTree *
DSTree_new ( 
   void
) {
DSTree   *dstree ;

#if MYTRACE > 0
fprintf(stdout, "\n just inside DSTree_new()") ;
fflush(stdout) ;
#endif

ALLOCATE(dstree, struct _DSTree, 1) ;

DSTree_setDefaultFields(dstree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving DSTree_new()") ;
fflush(stdout) ;
#endif

return(dstree) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- set the default fields for the DSTree object

   created -- 96mar10, cca
   ------------------------------------------------------
*/
void
DSTree_setDefaultFields (
   DSTree   *dstree
) {

#if MYTRACE > 0
fprintf(stdout, "\n just inside DSTree_setDefaultFields(%)", g) ;
fflush(stdout) ;
#endif

if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_setDefaultFields(%p)"
           "\n dstree is NULL\n", dstree) ;
   exit(-1) ;
}
dstree->tree  = NULL ;
dstree->mapIV = NULL ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving DSTree_setDefaultFields(%)", dstree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 96mar10, cca
   --------------------------------
*/
void
DSTree_clearData (
   DSTree   *dstree
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside DSTree_clearData(%)", dstree) ;
fflush(stdout) ;
#endif

if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_clearData(%p)"
           "\n dstree is NULL\n", dstree) ;
   exit(-1) ;
}

if ( dstree->tree != NULL ) {
   Tree_clearData(dstree->tree) ;
   Tree_free(dstree->tree) ;
}
if ( dstree->mapIV != NULL ) {
   IV_free(dstree->mapIV) ;
}
DSTree_setDefaultFields(dstree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving DSTree_clearData(%)", dstree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- free the DSTree object

   created -- 96mar10, cca
   --------------------------------
*/
void
DSTree_free (
   DSTree   *dstree
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside DSTree_free(%)", dstree) ;
fflush(stdout) ;
#endif

if ( dstree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_free(%p)"
           "\n dstree is NULL\n", dstree) ;
   exit(-1) ;
}

DSTree_clearData(dstree) ;
FREE(dstree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving DSTree_free(%)", dstree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
