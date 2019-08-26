/*  factorSetup.c  */

#include "../BridgeMPI.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- to construct the map from fronts to processors,
      and compute operations for each processor.

   maptype -- type of map for parallel factorization
      maptype = 1 --> wrap map
      maptype = 2 --> balanced map
      maptype = 3 --> subtree-subset map
      maptype = 4 --> domain decomposition map
   cutoff -- used when maptype = 4 as upper bound on
      relative domain size

   return value --
      1 -- success
     -1 -- bridge is NULL
     -2 -- front tree is NULL

   created -- 98sep25, cca
   ----------------------------------------------------------
*/
int
BridgeMPI_factorSetup (
   BridgeMPI   *bridge,
   int         maptype,
   double      cutoff
) {
double   t1, t2 ;
DV       *cumopsDV ;
ETree    *frontETree ;
FILE     *msgFile ;
int      msglvl, nproc ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t1) ;
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_factorSetup()"
           "\n bridge is NULL") ;
   return(-1) ;
}
if ( (frontETree = bridge->frontETree) == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_factorSetup()"
           "\n frontETree is NULL") ;
   return(-2) ;
}
nproc   = bridge->nproc   ;
msglvl  = bridge->msglvl  ;
msgFile = bridge->msgFile ;
/*
   -------------------------------------------
   allocate and initialize the cumopsDV object
   -------------------------------------------
*/
if ( (cumopsDV = bridge->cumopsDV) == NULL ) {
   cumopsDV = bridge->cumopsDV = DV_new() ;
}
DV_setSize(cumopsDV, nproc) ;
DV_zero(cumopsDV) ;
/*
   ----------------------------
   create the owners map object
   ----------------------------
*/
switch ( maptype ) {
case 1 :
   bridge->ownersIV = ETree_wrapMap(frontETree, bridge->type, 
                               bridge->symmetryflag, cumopsDV) ;
   break ;
case 2 :
   bridge->ownersIV = ETree_balancedMap(frontETree, bridge->type,
                               bridge->symmetryflag, cumopsDV) ;
   break ;
case 3 :
   bridge->ownersIV = ETree_subtreeSubsetMap(frontETree, bridge->type,
                               bridge->symmetryflag, cumopsDV) ;
   break ;
case 4 :
   bridge->ownersIV = ETree_ddMap(frontETree, bridge->type,
                               bridge->symmetryflag, cumopsDV, cutoff) ;
   break ;
default :
   bridge->ownersIV = ETree_ddMap(frontETree, bridge->type,
                         bridge->symmetryflag, cumopsDV, 1./(2*nproc)) ;
   break ;
}
MARKTIME(t2) ;
bridge->cpus[7] = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n parallel factor setup") ;
   fprintf(msgFile, "\n type = %d, symmetryflag = %d",
           bridge->type, bridge->symmetryflag) ;
   fprintf(msgFile, "\n total factor operations = %.0f",
           DV_sum(cumopsDV)) ;
   fprintf(msgFile, 
           "\n upper bound on speedup due to load balance = %.2f",
           DV_max(cumopsDV)/DV_sum(cumopsDV)) ;
   fprintf(msgFile, "\n operations distributions over threads") ;
   DV_writeForHumanEye(cumopsDV, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n owners map IV object") ;
   IV_writeForHumanEye(bridge->ownersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------
   create the vertex map object
   ----------------------------
*/
bridge->vtxmapIV = IV_new() ;
IV_init(bridge->vtxmapIV, bridge->neqns, NULL) ;
IVgather(bridge->neqns, IV_entries(bridge->vtxmapIV),
         IV_entries(bridge->ownersIV), 
         ETree_vtxToFront(bridge->frontETree)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n vertex map IV object") ;
   IV_writeForHumanEye(bridge->vtxmapIV, msgFile) ;
   fflush(msgFile) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
