/*  CleanupMT.c  */

#include "../BridgeMT.h"

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
CleanupMT (
   void   *data
) {
BridgeMT   *bridge = (BridgeMT *) data ;
#if MYDEBUG > 0
double   t1, t2 ;
MARKTIME(t1) ;
count_Cleanup++ ;
fprintf(stdout, "\n (%d) CleanupMT()", count_Cleanup) ;
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
IV_free(bridge->ownersIV) ;
DenseMtx_free(bridge->X) ;
DenseMtx_free(bridge->Y) ;
SolveMap_free(bridge->solvemap) ;

#if MYDEBUG > 0
MARKTIME(t2) ;
time_Cleanup += t2 - t1 ;
fprintf(stdout, ", %8.3f seconds, %8.3f total seconds, ",
        t2 - t1, time_Cleanup) ;
#endif

return(1) ; }

/*--------------------------------------------------------------------*/
