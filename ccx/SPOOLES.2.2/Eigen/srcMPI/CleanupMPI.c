/*  CleanupMPI.c  */

#include "../BridgeMPI.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_Cleanup = 0 ;
static double time_Cleanup = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to free the owned data structures

   return values --
      1 -- normal return
     -1 -- data is NULL

   created -- 98aug10, cca
   --------------------------------------------
*/
int
CleanupMPI (
   void   *data
) {
BridgeMPI   *bridge = (BridgeMPI *) data ;
#if MYDEBUG > 0
double   t1, t2 ;
MARKTIME(t1) ;
count_Cleanup++ ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) CleanupMPI()", count_Cleanup) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, "\n (%d) CleanupMPI()", count_Cleanup) ;
fflush(bridge->msgFile) ;
#endif

bridge->pencil->inpmtxA = NULL ;
bridge->pencil->inpmtxB = NULL ;
Pencil_free(bridge->pencil) ;
IVL_free(bridge->symbfacIVL) ;
FrontMtx_free(bridge->frontmtx) ;
ETree_free(bridge->frontETree) ;
SubMtxManager_free(bridge->mtxmanager) ;
IV_free(bridge->oldToNewIV) ;
IV_free(bridge->newToOldIV) ;
IV_free(bridge->vtxmapIV) ;
IV_free(bridge->ownersIV) ;
IV_free(bridge->myownedIV) ;
if ( bridge->rowmapIV != NULL ) {
   IV_free(bridge->rowmapIV) ;
}
if ( bridge->info != NULL ) {
   MatMul_cleanup(bridge->info) ;
   DenseMtx_free(bridge->Xloc) ;
   DenseMtx_free(bridge->Yloc) ;
}
SolveMap_free(bridge->solvemap) ;

#if MYDEBUG > 0
MARKTIME(t2) ;
time_Cleanup += t2 - t1 ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, ", %8.3f seconds, %8.3f total seconds, ",
           t2 - t1, time_Cleanup) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, ", %8.3f seconds, %8.3f total seconds, ",
        t2 - t1, time_Cleanup) ;
fflush(bridge->msgFile) ;
#endif

return(1) ; }

/*--------------------------------------------------------------------*/
