/*  solveSetup.c  */

#include "../BridgeMT.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to setup for the parallel solve

   return value ---
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- frontmtx is NULL
     -3 -- frontmtx has not yet been postprocessed

   created -- 98sep24, cca
   -----------------------------------------------
*/
int
BridgeMT_solveSetup (
   BridgeMT   *bridge
) {
double     t1, t2 ;
FILE       *msgFile ;
FrontMtx   *frontmtx ;
int        msglvl ;
SolveMap   *solvemap ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t1) ;
if ( bridge == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMT_solveSetup()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( (frontmtx = bridge->frontmtx) == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMT_solveSetup()"
           "\n frontmtx is NULL\n") ;
   return(-2) ;
}
if ( ! FRONTMTX_IS_2D_MODE(frontmtx) ) {
   fprintf(stderr, "\n\n error in BridgeMT_solveSetup()"
           "\n frontmtx must be in 2-D mode\n") ;
   return(-2) ;
}
msglvl  = bridge->msglvl  ;
msgFile = bridge->msgFile ;
if ( (solvemap = bridge->solvemap) == NULL ) {
   solvemap = bridge->solvemap = SolveMap_new() ;
} else {
   SolveMap_clearData(solvemap) ;
}
if (  FRONTMTX_IS_NONSYMMETRIC(frontmtx)
   && FRONTMTX_IS_PIVOTING(frontmtx) ) {
   SolveMap_ddMap(solvemap, SPOOLES_NONSYMMETRIC, 
                  frontmtx->upperblockIVL, frontmtx->lowerblockIVL,
                  bridge->nthread, bridge->ownersIV, frontmtx->tree,
                  bridge->seed, msglvl, msgFile) ;
} else {
   SolveMap_ddMap(solvemap, SPOOLES_SYMMETRIC, 
                  frontmtx->upperblockIVL, NULL,
                  bridge->nthread, bridge->ownersIV, frontmtx->tree,
                  bridge->seed, msglvl, msgFile) ;
}
MARKTIME(t2) ;
bridge->cpus[11] = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n solve map created") ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n SolveMap") ;
   SolveMap_writeForHumanEye(solvemap, msgFile) ;
   fflush(msgFile) ;
}

return(1) ; }   

/*--------------------------------------------------------------------*/
