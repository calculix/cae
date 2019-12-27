/*  basics.C  */

#include "../BridgeMPI.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor method

   created -- 98sep25, cca
   -----------------------
*/
BridgeMPI *
BridgeMPI_new ( 
   void 
) {
BridgeMPI   *bridge ;

ALLOCATE(bridge, struct _BridgeMPI, 1) ;

BridgeMPI_setDefaultFields(bridge) ;

return(bridge) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   return value ---
      1 -- normal return
     -1 -- bridge is NULL

   created -- 98sep25, cca
   -----------------------
*/
int
BridgeMPI_setDefaultFields ( 
   BridgeMPI   *bridge
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n fatal error in BridgeMPI_setDefaultFields(%p)"
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
bridge->type           = SPOOLES_REAL ;
bridge->symmetryflag   = SPOOLES_SYMMETRIC ;
bridge->sparsityflag   = FRONTMTX_DENSE_FRONTS ;
bridge->pivotingflag   = SPOOLES_NO_PIVOTING ;
bridge->tau            = 100.0 ;
bridge->droptol        = 1.e-3 ;
bridge->patchinfo      = NULL ;
bridge->lookahead      =   0  ;
/*
   ---------------------------
   MPI information and objects
   ---------------------------
*/
bridge->nproc          =   0  ;
bridge->myid           =  -1  ;
bridge->comm           = NULL ;
bridge->ownersIV       = NULL ;
bridge->solvemap       = NULL ;
bridge->cumopsDV       = NULL ;
bridge->vtxmapIV       = NULL ;
bridge->rowmapIV       = NULL ;
bridge->ownedColumnsIV = NULL ;
bridge->Aloc           = NULL ;
bridge->Xloc           = NULL ;
bridge->Yloc           = NULL ;
/*
   ------------------------------------
   message info, statistics and timings
   ------------------------------------
*/
IVzero(6, bridge->stats) ;
DVzero(22, bridge->cpus) ;
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

   created -- 98sep25, cca
   -----------------------
*/
int
BridgeMPI_clearData ( 
   BridgeMPI   *bridge
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n fatal error in BridgeMPI_clearData(%p)"
           "\n bad input\n", bridge) ;
   return(-1) ;
}
/*
   ------------
   free objects
   ------------
*/
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
/*
   ----------------
   free MPI objects
   ----------------
*/
if ( bridge->ownersIV != NULL ) {
   IV_free(bridge->ownersIV) ;
}
if ( bridge->solvemap != NULL ) {
   SolveMap_free(bridge->solvemap) ;
}
if ( bridge->cumopsDV != NULL ) {
   DV_free(bridge->cumopsDV) ;
}
if ( bridge->rowmapIV != bridge->vtxmapIV ) {
   IV_free(bridge->rowmapIV) ;
}
if ( bridge->vtxmapIV != NULL ) {
   IV_free(bridge->vtxmapIV) ;
}
if ( bridge->ownedColumnsIV != NULL ) { 
   IV_free(bridge->ownedColumnsIV) ; 
}
if ( bridge->Aloc != NULL ) { 
   InpMtx_free(bridge->Aloc) ; 
}
if ( bridge->Xloc != NULL ) { 
   DenseMtx_free(bridge->Xloc) ; 
}
if ( bridge->Yloc != NULL ) { 
   DenseMtx_free(bridge->Yloc) ; 
}
/*
   ------------------
   set default fields
   ------------------
*/
BridgeMPI_setDefaultFields(bridge) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL

   created -- 98sep25, cca
   -----------------------
*/
int
BridgeMPI_free ( 
   BridgeMPI   *bridge
) {
if ( bridge == NULL ) {
   fprintf(stderr, "\n fatal error in BridgeMPI_free(%p)"
           "\n bad input\n", bridge) ;
   return(-1) ;
}
BridgeMPI_clearData(bridge) ;
FREE(bridge) ;

return(1) ; }

/*--------------------------------------------------------------------*/
