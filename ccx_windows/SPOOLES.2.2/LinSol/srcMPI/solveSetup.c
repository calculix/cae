/*  solveSetup.c  */

#include "../BridgeMPI.h"

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
BridgeMPI_solveSetup (
   BridgeMPI   *bridge
) {
double     t0, t1, t2 ;
FILE       *msgFile ;
FrontMtx   *frontmtx ;
int        firsttag, msglvl, myid, nproc ;
int        stats[4] ;
IV         *ownersIV ;
MPI_Comm   comm ;
SolveMap   *solvemap ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if ( bridge == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_solveSetup()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( (frontmtx = bridge->frontmtx) == NULL ) {
   fprintf(stderr, "\n\n error in BridgeMPI_solveSetup()"
           "\n frontmtx is NULL\n") ;
   return(-2) ;
}
if ( ! FRONTMTX_IS_2D_MODE(frontmtx) ) {
   fprintf(stderr, "\n\n error in BridgeMPI_solveSetup()"
           "\n frontmtx must be in 2-D mode\n") ;
   return(-2) ;
}
msglvl   = bridge->msglvl   ;
msgFile  = bridge->msgFile  ;
myid     = bridge->myid     ;
comm     = bridge->comm     ;
nproc    = bridge->nproc    ;
ownersIV = bridge->ownersIV ;

if ( (solvemap = bridge->solvemap) == NULL ) {
   solvemap = bridge->solvemap = SolveMap_new() ;
} else {
   SolveMap_clearData(solvemap) ;
}
if (  FRONTMTX_IS_NONSYMMETRIC(frontmtx)
   && FRONTMTX_IS_PIVOTING(frontmtx) ) {
   SolveMap_ddMap(solvemap, SPOOLES_NONSYMMETRIC, 
                  frontmtx->upperblockIVL, frontmtx->lowerblockIVL,
                  nproc, ownersIV, frontmtx->tree,
                  bridge->seed, msglvl, msgFile) ;
} else {
   SolveMap_ddMap(solvemap, SPOOLES_SYMMETRIC, 
                  frontmtx->upperblockIVL, NULL,
                  nproc, ownersIV, frontmtx->tree,
                  bridge->seed, msglvl, msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n solve map created") ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n SolveMap") ;
   SolveMap_writeForHumanEye(solvemap, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------
   split the front matrix
   ----------------------
*/
IVzero(4, stats) ;
firsttag = 0 ;
MARKTIME(t1) ;
FrontMtx_MPI_split(frontmtx, solvemap, stats, msglvl, msgFile,
                   firsttag, comm) ;
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n CPU %8.3f : split the matrix", t2 - t1) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n submatrix manager after split") ;
   SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n frontmtx after split") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   generate bridge->rowmapIV, the map from 
   rows of the factor to owning processors
   ---------------------------------------
*/
if ( FRONTMTX_IS_PIVOTING(bridge->frontmtx) ) {
/*
   -----------------------------------------
   factorization done with pivoting, row map
   is different than for the original matrix
   -----------------------------------------
*/
   if ( bridge->rowmapIV == NULL ) {
/*
      ------------------------------------------------
      this is the first solve since the factorization,
      create the rowmap IV object
      ------------------------------------------------
*/
      bridge->rowmapIV = FrontMtx_MPI_rowmapIV(frontmtx, ownersIV,
                                              msglvl, msgFile, comm) ;
   }
} else {
/*
   -------------------------------------------
   no pivoting, map is simply the original map
   -------------------------------------------
*/
   bridge->rowmapIV = bridge->vtxmapIV ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n row map IV object") ;
   IV_writeForHumanEye(bridge->rowmapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------
   create ownedColumnsIV, a vector of column 
   ids that are owned by this processor.
   -----------------------------------------
*/
bridge->ownedColumnsIV = FrontMtx_ownedColumnsIV(frontmtx, myid, 
                                            ownersIV, msglvl, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n ownedColumns IV object") ;
   IV_writeForHumanEye(bridge->ownedColumnsIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------
   create DenseMtx objects to hold the local 
   solution and right hand side matrices
   -----------------------------------------
*/
bridge->Xloc = DenseMtx_new() ;
bridge->Yloc = DenseMtx_new() ;

MARKTIME(t2) ;
bridge->cpus[14] += t2 - t0 ;

return(1) ; }

/*--------------------------------------------------------------------*/
