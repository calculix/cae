/*  basics.c  */

#include "../ETree.h"

#define MYTRACE 0
#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- create and return a new ETree object

   created -- 95nov15, cca
   -----------------------------------------------
*/
ETree *
ETree_new ( 
   void
) {
ETree   *etree ;

#if MYTRACE > 0
fprintf(stdout, "\n just inside ETree_new()") ;
fflush(stdout) ;
#endif

ALLOCATE(etree, struct _ETree, 1) ;

ETree_setDefaultFields(etree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving ETree_new()") ;
fflush(stdout) ;
#endif

return(etree) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- set the default fields for the ETree object

   created -- 95nov15, cca
   ------------------------------------------------------
*/
void
ETree_setDefaultFields (
   ETree   *etree
) {

#if MYTRACE > 0
fprintf(stdout, "\n just inside ETree_setDefaultFields(%)", g) ;
fflush(stdout) ;
#endif

if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_setDefaultFields(%p)"
           "\n etree is NULL\n", etree) ;
   exit(-1) ;
 }
etree->nfront       =   0  ;
etree->nvtx         =   0  ;
etree->tree         = NULL ;
etree->nodwghtsIV   = NULL ;
etree->bndwghtsIV   = NULL ;
etree->vtxToFrontIV = NULL ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving ETree_setDefaultFields(%)", etree) ;
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
ETree_clearData (
   ETree   *etree
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside ETree_clearData(%)", etree) ;
fflush(stdout) ;
#endif

if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_clearData(%p)"
           "\n etree is NULL\n", etree) ;
   exit(-1) ;
}

if ( etree->tree != NULL ) {
   Tree_free(etree->tree) ;
}
if ( etree->nodwghtsIV != NULL ) {
   IV_free(etree->nodwghtsIV) ;
}
if ( etree->bndwghtsIV != NULL ) {
   IV_free(etree->bndwghtsIV) ;
}
if ( etree->vtxToFrontIV != NULL ) {
   IV_free(etree->vtxToFrontIV) ;
}
ETree_setDefaultFields(etree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving ETree_clearData(%)", etree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- free the ETree object

   created -- 95nov15, cca
   --------------------------------
*/
void
ETree_free (
   ETree   *etree
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside ETree_free(%)", etree) ;
fflush(stdout) ;
#endif

if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_free(%p)"
           "\n etree is NULL\n", etree) ;
   exit(-1) ;
}

ETree_clearData(etree) ;
FREE(etree) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving ETree_free(%)", etree) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
