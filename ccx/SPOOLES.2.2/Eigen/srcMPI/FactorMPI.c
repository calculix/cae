/*  FactorMPI.c  */

#include "../BridgeMPI.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_Factor = 0 ;
static double time_Factor = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to compute the factorization of A - sigma * B

   note: all variables in the calling sequence are references
         to allow call from fortran.

   input parameters 

      data    -- pointer to bridge data object
      psigma  -- shift for the matrix pencil
      ppvttol -- pivot tolerance
         *ppvttol =  0.0 --> no pivoting used
         *ppvttol != 0.0 --> pivoting used, entries in factor are
                             bounded above by 1/pvttol in magnitude

   output parameters 

      *pinertia -- on return contains the number of negative eigenvalues
      *perror   -- on return contains an error code
          1 -- error found during factorization
          0 -- normal return
         -1 -- psigma is NULL
         -2 -- ppvttol is NULL
         -3 -- data is NULL
         -4 -- pinertia is NULL

   created -- 98aug10, cca & jcp
   ---------------------------------------------------------------------
*/
void
FactorMPI ( 
   double     *psigma, 
   double     *ppvttol, 
   void       *data,
   int        *pinertia,
   int        *perror
) {
BridgeMPI    *bridge = (BridgeMPI *) data ; 
Chv          *rootchv ;
ChvManager   *chvmanager ;
double       droptol=0.0, tau ;
double       cpus[20] ;
FILE         *msgFile ;
int          recvtemp[3], sendtemp[3], stats[20] ;
int          msglvl, nnegative, nzero, npositive, pivotingflag, tag ;
MPI_Comm     comm ;
int          nproc ;

#if MYDEBUG > 0
double   t1, t2 ;
count_Factor++ ;
MARKTIME(t1) ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, "\n (%d) FactorMPI()", count_Factor) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, "\n (%d) FactorMPI()", count_Factor) ;
fflush(bridge->msgFile) ;
#endif

nproc = bridge->nproc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( psigma == NULL ) {
   fprintf(stderr, "\n error in FactorMPI()"
           "\n psigma is NULL\n") ;
   *perror = -1 ; return ;
}
if ( ppvttol == NULL ) {
   fprintf(stderr, "\n error in FactorMPI()"
           "\n ppvttol is NULL\n") ;
   *perror = -2 ; return ;
}
if ( data == NULL ) {
   fprintf(stderr, "\n error in FactorMPI()"
           "\n data is NULL\n") ;
   *perror = -3 ; return ;
}
if ( pinertia == NULL ) {
   fprintf(stderr, "\n error in FactorMPI()"
           "\n pinertia is NULL\n") ;
   *perror = -4 ; return ;
}
if ( perror == NULL ) {
   fprintf(stderr, "\n error in FactorMPI()"
           "\n perror is NULL\n") ;
   return ;
}
comm    = bridge->comm    ;
msglvl  = bridge->msglvl  ;
msgFile = bridge->msgFile ;
/*
   ----------------------------------
   set the shift in the pencil object
   ----------------------------------
*/ 
bridge->pencil->sigma[0] = -(*psigma) ;
bridge->pencil->sigma[1] = 0.0 ;
/*
   -----------------------------------------
   if the matrices are in local coordinates
   (i.e., this is the first factorization 
    following a matrix-vector multiply) then
   map the matrix into global coordinates
   -----------------------------------------
*/
if ( bridge->coordFlag == LOCAL ) {
   if ( bridge->prbtype == 1 ) {
      MatMul_setGlobalIndices(bridge->info, bridge->B) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n matrix B in local coordinates") ;
         InpMtx_writeForHumanEye(bridge->B, msgFile) ;
         fflush(msgFile) ;
      }
   }
   if ( bridge->prbtype == 2 ) {
      MatMul_setGlobalIndices(bridge->info, bridge->A) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n matrix A in local coordinates") ;
         InpMtx_writeForHumanEye(bridge->A, msgFile) ;
         fflush(msgFile) ;
      }
   }
   bridge->coordFlag = GLOBAL ;
}
/*
   -----------------------------------------------------
   clear the front matrix and submatrix mananger objects
   -----------------------------------------------------
*/ 
FrontMtx_clearData(bridge->frontmtx);
SubMtxManager_clearData(bridge->mtxmanager);
SolveMap_clearData(bridge->solvemap) ;
if ( bridge->rowmapIV != NULL ) {
   IV_free(bridge->rowmapIV) ;
   bridge->rowmapIV = NULL ;
}
/*
   -----------------------------------------------------------
   set the pivot tolerance.
   NOTE: spooles's "tau" parameter is a bound on the magnitude 
   of the factor entries, and is the recipricol of that of the 
   pivot tolerance of the lanczos code
   -----------------------------------------------------------
*/ 
if ( *ppvttol == 0.0 ) {
   tau = 10.0 ;
   pivotingflag = SPOOLES_NO_PIVOTING ;
} else {
   tau = (1.0)/(*ppvttol) ;
   pivotingflag = SPOOLES_PIVOTING ;
}
/*
   ----------------------------------
   initialize the front matrix object
   ----------------------------------
*/ 
FrontMtx_init(bridge->frontmtx, bridge->frontETree, bridge->symbfacIVL,
              SPOOLES_REAL, SPOOLES_SYMMETRIC, FRONTMTX_DENSE_FRONTS,
              pivotingflag, NO_LOCK, bridge->myid, bridge->ownersIV, 
              bridge->mtxmanager, bridge->msglvl, bridge->msgFile) ;
/*
   -------------------------
   compute the factorization
   -------------------------
*/
tag = 0 ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, NO_LOCK, 0);
IVfill(20, stats, 0) ;
DVfill(20, cpus,  0.0) ;
rootchv = FrontMtx_MPI_factorPencil(bridge->frontmtx, bridge->pencil, 
                             tau, droptol, chvmanager, bridge->ownersIV,
                             0, perror, cpus, stats, bridge->msglvl, 
                             bridge->msgFile, tag, comm) ;
ChvManager_free(chvmanager);
tag += 3*FrontMtx_nfront(bridge->frontmtx) + 2 ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n numeric factorization") ;
   FrontMtx_writeForHumanEye(bridge->frontmtx, bridge->msgFile) ;
   fflush(bridge->msgFile) ;
}
/*
   ----------------------------
   if matrix is singular then
      set error flag and return
   ----------------------------
*/ 
if ( rootchv != NULL ) {
   fprintf(msgFile, "\n WHOA NELLY!, matrix is singular") ;
   fflush(msgFile) ;
   *perror = 1 ;
   return ;
}
/*
   ------------------------------------------------------------------
   post-process the factor matrix, convert from fronts to submatrices
   ------------------------------------------------------------------
*/ 
FrontMtx_MPI_postProcess(bridge->frontmtx, bridge->ownersIV, stats,
                         bridge->msglvl, bridge->msgFile, tag, comm);
tag += 5*bridge->nproc ;
/*
   -------------------
   compute the inertia
   -------------------
*/ 
FrontMtx_inertia(bridge->frontmtx, &nnegative, &nzero, &npositive) ;
sendtemp[0] = nnegative ;
sendtemp[1] = nzero     ;
sendtemp[2] = npositive ;
if ( bridge->msglvl > 2 && bridge->msgFile != NULL ) {
   fprintf(bridge->msgFile, "\n local inertia = < %d, %d, %d >",
           nnegative, nzero, npositive) ;
   fflush(bridge->msgFile) ;
}
MPI_Allreduce((void *) sendtemp, (void *) recvtemp, 3, MPI_INT, 
           MPI_SUM, comm) ;
nnegative = recvtemp[0] ;
nzero     = recvtemp[1] ;
npositive = recvtemp[2] ;
if ( bridge->msglvl > 2 && bridge->msgFile != NULL ) {
   fprintf(bridge->msgFile, "\n global inertia = < %d, %d, %d >",
           nnegative, nzero, npositive) ;
   fflush(bridge->msgFile) ;
}
*pinertia = nnegative;
/*
   ---------------------------
   create the solve map object
   ---------------------------
*/
SolveMap_ddMap(bridge->solvemap, SPOOLES_REAL,
               FrontMtx_upperBlockIVL(bridge->frontmtx),
               FrontMtx_lowerBlockIVL(bridge->frontmtx), nproc,
               bridge->ownersIV, FrontMtx_frontTree(bridge->frontmtx),
               bridge->seed, bridge->msglvl, bridge->msgFile) ;
/*
   -------------------------------
   redistribute the front matrices
   -------------------------------
*/
FrontMtx_MPI_split(bridge->frontmtx, bridge->solvemap, stats,
                   bridge->msglvl, bridge->msgFile, tag, comm) ;
if ( *ppvttol != 0.0 ) {
/*
   -------------------------------------------------------------
   pivoting for stability may have taken place. create rowmapIV, 
   the map from rows in the factorization to processes.
   -------------------------------------------------------------
*/
   bridge->rowmapIV = FrontMtx_MPI_rowmapIV(bridge->frontmtx,
                                       bridge->ownersIV, bridge->msglvl,
                                       bridge->msgFile, bridge->comm) ;
   if ( bridge->msglvl > 2 && bridge->msgFile != NULL ) {
      fprintf(bridge->msgFile, "\n\n bridge->rowmapIV") ;
      IV_writeForHumanEye(bridge->rowmapIV, bridge->msgFile) ;
      fflush(bridge->msgFile) ;
   }
} else {
   bridge->rowmapIV = NULL ;
}
/*
   ------------------------------------------------------------------
   set the error. (this is simple since when the spooles codes detect 
   a fatal error, they print out a message to stderr and exit.)
   ------------------------------------------------------------------
*/ 
*perror = 0 ;

#if MYDEBUG > 0
MARKTIME(t2) ;
time_Factor += t2 - t1 ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, ", %8.3f seconds, %8.3f total time",
           t2 - t1, time_Factor) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, ", %8.3f seconds, %8.3f total time",
        t2 - t1, time_Factor) ;
fflush(bridge->msgFile) ;
#endif
 
return; }

/*--------------------------------------------------------------------*/
