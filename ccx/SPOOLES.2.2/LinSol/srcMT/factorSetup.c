/*  factorSetup.c  */

#include "../BridgeMT.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to construct the map from fronts to threads,
      and compute operations for each thread.

   nthread -- number of threads
   maptype -- type of map for parallel factorization
     maptype = 1 --> wrap map
     maptype = 2 --> balanced map
     maptype = 3 --> subtree-subset map
     maptype = 4 --> domain decomposition map
   cutoff -- used when maptype = 4 as upper bound on
     relative domain size
  default is maptype = 4 and cutoff = 1/(2*nthread)

   return value --
      1 -- success
     -1 -- bridge is NULL
     -2 -- nthread is invalid, must be > 0
     -3 -- front tree is NULL

   created -- 98sep24, cca
   -------------------------------------------------------
*/
int
BridgeMT_factorSetup (
   BridgeMT   *bridge,
   int        nthread,
   int        maptype,
   double     cutoff
) {
double   t1, t2 ;
DV       *cumopsDV ;
ETree    *frontETree ;
FILE     *msgFile ;
int      msglvl ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t1) ;
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_factorSetup()"
           "\n bridge is NULL") ;
   return(-1) ;
}
if ( nthread < 1 ) {
   fprintf(stderr, "\n error in BridgeMT_factorSetup()"
           "\n nthread = %d, is invalid", nthread) ;
   return(-2) ;
}
if ( (frontETree = bridge->frontETree) == NULL ) {
   fprintf(stderr, "\n error in BridgeMT_factorSetup()"
           "\n frontETree is NULL") ;
   return(-5) ;
}
bridge->nthread = nthread ;
/*
   -------------------------------------------
   allocate and initialize the cumopsDV object
   -------------------------------------------
*/
if ( (cumopsDV = bridge->cumopsDV) == NULL ) {
   cumopsDV = bridge->cumopsDV = DV_new() ;
}
DV_setSize(cumopsDV, nthread) ;
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
                       bridge->symmetryflag, cumopsDV, 1./(2*nthread)) ;
   break ;
}
MARKTIME(t2) ;
bridge->cpus[5] = t2 - t1 ;
msglvl  = bridge->msglvl ;
msgFile = bridge->msgFile ;
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
return(1) ; }

/*--------------------------------------------------------------------*/
