/*  basics.C  */

#include "../BridgeMT.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor method

   created -- 98sep18, cca
   -----------------------
*/
BridgeMT *
BridgeMT_new ( 
   void 
) {
BridgeMT   *bridge ;

ALLOCATE(bridge, struct _BridgeMT, 1) ;

BridgeMT_setDefaultFields(bridge) ;

return(bridge) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   return value ---
      1 -- normal return
     -1 -- bridge is NULL

   created -- 98sep18, cca
   -----------------------
*/
int
BridgeMT_setDefaultFields ( 
   BridgeMT   *bridge
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n fatal error in BridgeMT_setDefaultFields(%p)"
           "\n bad input\n", bridge) ;
   return(-1) ;
}
/*
   ----------------
   graph statistics
   ----------------
*/
bridge->neqns  = 0 ;
bridge->nedges = 0 ;
bridge->Neqns  = 0 ;
bridge->Nedges = 0 ;
/*
   -------------------
   ordering parameters
   -------------------
*/
bridge->compressCutoff = 0.0 ;
bridge->maxdomainsize  = -1  ;
bridge->maxnzeros      = -1  ;
bridge->maxsize        = -1  ;
bridge->seed           = -1  ;
/*
   -------------------------------
   matrix/factorization parameters
   -------------------------------
*/
bridge->type         = SPOOLES_REAL ;
bridge->symmetryflag = SPOOLES_SYMMETRIC ;
bridge->sparsityflag = FRONTMTX_DENSE_FRONTS ;
bridge->pivotingflag = SPOOLES_NO_PIVOTING ;
bridge->tau          = 100.0 ;
bridge->droptol      = 1.e-3 ;
bridge->lookahead    =   0  ;
bridge->patchinfo    = NULL ;
/*
   ------------------------
   multithreaded parameters
   ------------------------
*/
bridge->nthread   =   0  ;
bridge->ownersIV  = NULL ;
bridge->solvemap  = NULL ;
bridge->cumopsDV  = NULL ;
/*
   ------------------------------------
   message info, statistics and timings
   ------------------------------------
*/
IVzero(6, bridge->stats) ;
DVzero(16, bridge->cpus) ;
bridge->msglvl  =    0   ;
bridge->msgFile = stdout ;
/*
   -------------------
   pointers to objects
   -------------------
*/
bridge->frontETree = NULL ;
bridge->symbfacIVL = NULL ;
bridge->mtxmanager = NULL ;
bridge->frontmtx   = NULL ;
bridge->oldToNewIV = NULL ;
bridge->newToOldIV = NULL ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   return value ---
      1 -- normal return
     -1 -- bridge is NULL

   created -- 98sep18, cca
   -----------------------
*/
int
BridgeMT_clearData ( 
   BridgeMT   *bridge
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n fatal error in BridgeMT_clearData(%p)"
           "\n bad input\n", bridge) ;
   return(-1) ;
}
if ( bridge->frontmtx != NULL ) {
   FrontMtx_free(bridge->frontmtx) ;
}
if ( bridge->frontETree != NULL ) {
   ETree_free(bridge->frontETree) ;
}
if ( bridge->symbfacIVL != NULL ) {
   IVL_free(bridge->symbfacIVL) ;
}
if ( bridge->mtxmanager != NULL ) {
   SubMtxManager_free(bridge->mtxmanager) ;
}
if ( bridge->oldToNewIV != NULL ) {
   IV_free(bridge->oldToNewIV) ;
}
if ( bridge->newToOldIV != NULL ) {
   IV_free(bridge->newToOldIV) ;
}
if ( bridge->ownersIV != NULL ) {
   IV_free(bridge->ownersIV) ;
}
if ( bridge->solvemap != NULL ) {
   SolveMap_free(bridge->solvemap) ;
}
if ( bridge->cumopsDV != NULL ) {
   DV_free(bridge->cumopsDV) ;
}
BridgeMT_setDefaultFields(bridge) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   return value ---
      1 -- normal return
     -1 -- bridge is NULL

   created -- 98sep18, cca
   -----------------------
*/
int
BridgeMT_free ( 
   BridgeMT   *bridge
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n fatal error in BridgeMT_free(%p)"
           "\n bad input\n", bridge) ;
   return(-1) ;
}
BridgeMT_clearData(bridge) ;
FREE(bridge) ;

return(1) ; }

/*--------------------------------------------------------------------*/
