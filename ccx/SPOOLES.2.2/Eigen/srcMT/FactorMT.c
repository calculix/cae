/*  FactorMT.c  */

#include "../BridgeMT.h"

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
FactorMT ( 
   double   *psigma, 
   double   *ppvttol, 
   void     *data,
   int      *pinertia,
   int      *perror
) {
BridgeMT     *bridge = (BridgeMT *) data ; 
Chv          *rootchv ;
ChvManager   *chvmanager ;
double       droptol=0.0, tau ;
double       cpus[10] ;
int          stats[20] ;
int          nnegative, nzero, npositive, pivotingflag ;
int          nproc ;
#if MYDEBUG > 0
double   t1, t2 ;
MARKTIME(t1) ;
count_Factor++ ;
fprintf(stdout, "\n (%d) FactorMT()", count_Factor) ;
fflush(stdout) ;
#endif
/*
   ---------------
   check the input
   ---------------
*/
if ( psigma == NULL ) {
   fprintf(stderr, "\n error in FactorMT()"
           "\n psigma is NULL\n") ;
   *perror = -1 ; return ;
}
if ( ppvttol == NULL ) {
   fprintf(stderr, "\n error in FactorMT()"
           "\n ppvttol is NULL\n") ;
   *perror = -2 ; return ;
}
if ( data == NULL ) {
   fprintf(stderr, "\n error in FactorMT()"
           "\n data is NULL\n") ;
   *perror = -3 ; return ;
}
if ( pinertia == NULL ) {
   fprintf(stderr, "\n error in FactorMT()"
           "\n pinertia is NULL\n") ;
   *perror = -4 ; return ;
}
if ( perror == NULL ) {
   fprintf(stderr, "\n error in FactorMT()"
           "\n perror is NULL\n") ;
   return ;
}
nproc = bridge->nthread ;
/*
   ----------------------------------
   set the shift in the pencil object
   ----------------------------------
*/ 
bridge->pencil->sigma[0] = -(*psigma) ;
bridge->pencil->sigma[1] = 0.0 ;
/*
   ----------------------------------------------------------------
   clear the front matrix, submatrix mananger and solve map objects
   ----------------------------------------------------------------
*/ 
FrontMtx_clearData(bridge->frontmtx);
SubMtxManager_clearData(bridge->mtxmanager);
SolveMap_clearData(bridge->solvemap) ;
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
              pivotingflag, LOCK_IN_PROCESS, 0, NULL, 
              bridge->mtxmanager, bridge->msglvl, bridge->msgFile) ;
/*
   -------------------------
   compute the factorization
   -------------------------
*/
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, LOCK_IN_PROCESS, 1);
IVfill(20, stats, 0) ;
DVfill(10, cpus, 0.0) ;
rootchv = FrontMtx_MT_factorPencil(bridge->frontmtx, bridge->pencil, 
                 tau, droptol, chvmanager, bridge->ownersIV, 0,
                 perror, cpus, stats, bridge->msglvl, bridge->msgFile) ;
ChvManager_free(chvmanager);
/*
   ----------------------------
   if matrix is singular then
      set error flag and return
   ----------------------------
*/ 
if ( rootchv != NULL ) {
   *perror = 1 ;
   return ;
}
/*
   ------------------------------------------------------------------
   post-process the factor matrix, convert from fronts to submatrices
   ------------------------------------------------------------------
*/ 
FrontMtx_postProcess(bridge->frontmtx, bridge->msglvl, bridge->msgFile);
/*
   -------------------
   compute the inertia
   -------------------
*/ 
FrontMtx_inertia(bridge->frontmtx, &nnegative, &nzero, &npositive) ;
*pinertia = nnegative;
/*
   --------------------------
   set up the SolveMap object
   --------------------------
*/
SolveMap_ddMap(bridge->solvemap, SPOOLES_SYMMETRIC,
               FrontMtx_upperBlockIVL(bridge->frontmtx),
               FrontMtx_lowerBlockIVL(bridge->frontmtx), nproc,
               bridge->ownersIV, FrontMtx_frontTree(bridge->frontmtx),
               bridge->seed, bridge->msglvl, bridge->msgFile ) ;
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
fprintf(stdout, ", %8.3f seconds, %8.3f total time",
        t2 - t1, time_Factor) ;
fflush(stdout) ;
#endif
 
return; }

/*--------------------------------------------------------------------*/
