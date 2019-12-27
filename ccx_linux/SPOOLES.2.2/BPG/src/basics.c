/*  basics.c  */

#include "../BPG.h"

#define MYTRACE 0
#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- create and return a new BPG object

   created -- 95oct06, cca
   -----------------------------------------------
*/
BPG *
BPG_new ( 
   void
) {
BPG   *bpg ;

#if MYTRACE > 0
fprintf(stdout, "\n just inside BPG_new()") ;
fflush(stdout) ;
#endif

ALLOCATE(bpg, struct _BPG, 1) ;

BPG_setDefaultFields(bpg) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving BPG_new()") ;
fflush(stdout) ;
#endif

return(bpg) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- set the default fields for the BPG object

   created -- 95oct06, cca
   ------------------------------------------------------
*/
void
BPG_setDefaultFields (
   BPG   *bpg
) {

#if MYTRACE > 0
fprintf(stdout, "\n just inside BPG_setDefaultFields(%)", bpg) ;
fflush(stdout) ;
#endif

if ( bpg == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_setDefaultFields(%p)"
           "\n bipartite graph is NULL\n", bpg) ;
   exit(-1) ;
}
bpg->nX    =  0   ;
bpg->nY    =  0   ;
bpg->graph = NULL ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving BPG_setDefaultFields(%)", bpg) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 95oct06, cca
   --------------------------------
*/
void
BPG_clearData (
   BPG   *bpg
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside BPG_clearData(%)", bpg) ;
fflush(stdout) ;
#endif

if ( bpg == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_clearData(%p)"
           "\n bipartite graph is NULL\n", bpg) ;
   exit(-1) ;
}

if ( bpg->graph != NULL ) {
   Graph_free(bpg->graph) ;
}
BPG_setDefaultFields(bpg) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving BPG_clearData(%)", bpg) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- free the BPG object

   created -- 95oct06, cca
   --------------------------------
*/
void
BPG_free (
   BPG   *bpg
) {
#if MYTRACE > 0
fprintf(stdout, "\n just inside BPG_free(%)", bpg) ;
fflush(stdout) ;
#endif

if ( bpg == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_free(%p)"
           "\n bipartite graph is NULL\n", bpg) ;
   exit(-1) ;
}

BPG_clearData(bpg) ;

FREE(bpg) ;

#if MYTRACE > 0
fprintf(stdout, "\n leaving BPG_free(%)", bpg) ;
fflush(stdout) ;
#endif

return ; }

/*--------------------------------------------------------------------*/
