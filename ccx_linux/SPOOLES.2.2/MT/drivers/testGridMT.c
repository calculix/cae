/*  testGridMT.c  */

#include "../spoolesMT.h"
#include "../../FrontMtx.h"
#include "../../SolveMap.h"
#include "../../SymbFac.h"
#include "../../Drand.h"
#include "../../timings.h"
#include "../../misc.h"

/*--------------------------------------------------------------------*/
/*
void mkNDlinsys ( int n1, int n2, int n3, int maxzeros, int maxsize,
   int type, int symmetryflag, int nrhs, int seed, int msglvl,
   FILE *msgFile, ETree **pfrontETree, IVL **psymbfacIVL,
   InpMtx **pmtxA, DenseMtx **pmtxX, DenseMtx **pmtxB ) ;
*/
/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------------
   test the factor method for a grid matrix
   (1) construct a linear system for a nested dissection
       ordering on a regular grid
   (2) create a solution matrix object
   (3) multiply the solution with the matrix
       to get a right hand side matrix object
   (4) factor the matrix 
   (5) solve the system
 
   created -- 98may16, cca
   -----------------------------------------------------
*/
{
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxB, *mtxX, *mtxZ ;
FrontMtx        *frontmtx ;
InpMtx          *mtxA ;
SubMtxManager   *mtxmanager ;
double          cputotal, cutoff, droptol, factorops ;
double          cpus[11] ;
Drand           drand ;
double          nops, tau, t1, t2   ;
DV              *cumopsDV ;
ETree           *frontETree ;
FILE            *msgFile ;
int             error, lookahead, maptype, maxsize, maxzeros, msglvl, 
                neqns, n1, n2, n3, nrhs, nthread, nzf, pivotingflag, 
                seed, sparsityflag, symmetryflag,
                type ;
int             stats[16] ;
IV              *frontOwnersIV ;
IVL             *symbfacIVL ;
SolveMap        *solvemap ;

if ( argc != 20 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile n1 n2 n3 maxzeros maxsize" 
"\n         seed type symmetryflag sparsityflag pivotingflag "
"\n         tau droptol nrhs nthread maptype cutoff lookahead"
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    n1           -- number of grid points in the first direction"
"\n    n2           -- number of grid points in the second direction"
"\n    n3           -- number of grid points in the third direction"
"\n    maxzeros     -- max number of zeroes in a front"
"\n    maxsize      -- max number of internal nodes in a front"
"\n    seed         -- random number seed"
"\n    type         -- type of entries"
"\n       1 --> real entries"
"\n       2 --> complex entries"
"\n    symmetryflag -- symmetry flag"
"\n       0 --> symmetric"
"\n       1 --> hermitian"
"\n       2 --> nonsymmetric"
"\n    sparsityflag -- sparsity flag"
"\n       0 --> store dense fronts"
"\n       1 --> store sparse fronts, use droptol to drop entries"
"\n    pivotingflag -- pivoting flag"
"\n       0 --> do not pivot"
"\n       1 --> enable pivoting"
"\n    tau     -- upper bound on factor entries"
"\n               used only with pivoting"
"\n    droptol -- lower bound on factor entries"
"\n               used only with sparse fronts"
"\n    nrhs    -- # of right hand sides"
"\n    nthread -- number of threads"
"\n    maptype -- type of map from fronts to threads"
"\n       1 --> wrap map"
"\n       2 --> balanced map via a post-order traversal"
"\n       3 --> subtree-subset map"
"\n       4 --> domain decomposition map"
"\n    cutoff  -- cutoff used for domain decomposition map"
"\n       0 <= cutoff <= 1 used to define the multisector"
"\n lookahead -- lookahead for controlling the parallelism"
"\n", argv[0]) ;
   return(-1) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
n1            = atoi(argv[3]) ;
n2            = atoi(argv[4]) ;
n3            = atoi(argv[5]) ;
maxzeros      = atoi(argv[6]) ;
maxsize       = atoi(argv[7]) ;
seed          = atoi(argv[8]) ;
type          = atoi(argv[9]) ;
symmetryflag  = atoi(argv[10]) ;
sparsityflag  = atoi(argv[11]) ;
pivotingflag  = atoi(argv[12]) ;
tau           = atof(argv[13]) ;
droptol       = atof(argv[14]) ;
nrhs          = atoi(argv[15]) ;
nthread       = atoi(argv[16]) ;
maptype       = atoi(argv[17]) ;
cutoff        = atof(argv[18]) ;
lookahead     = atoi(argv[19]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n n1            -- %d" 
        "\n n2            -- %d" 
        "\n n3            -- %d" 
        "\n maxzeros      -- %d" 
        "\n maxsize       -- %d" 
        "\n seed          -- %d" 
        "\n type          -- %d" 
        "\n symmetryflag  -- %d" 
        "\n sparsityflag  -- %d" 
        "\n pivotingflag  -- %d" 
        "\n tau           -- %e" 
        "\n droptol       -- %e" 
        "\n nrhs          -- %d" 
        "\n nthread       -- %d" 
        "\n maptype       -- %d" 
        "\n cutoff        -- %f" 
        "\n lookahead     -- %d" 
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, maxzeros, maxsize, seed, 
        type, symmetryflag, sparsityflag, pivotingflag, tau, droptol, 
        nrhs, nthread, maptype, cutoff, lookahead);
fflush(msgFile) ;
neqns = n1 * n2 * n3 ;
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
Drand_setDefaultFields(&drand) ;
Drand_init(&drand) ;
Drand_setSeed(&drand, seed) ;
/*
Drand_setUniform(&drand, 0.0, 1.0) ;
*/
Drand_setNormal(&drand, 0.0, 1.0) ;
/*
   --------------------------
   generate the linear system
   --------------------------
*/
mkNDlinsys(n1, n2, n3, maxzeros, maxsize, type,
           symmetryflag, nrhs, seed, msglvl, msgFile,
           &frontETree, &symbfacIVL, &mtxA, &mtxX, &mtxB) ;
/*
   --------------------------------------------------
   initialize the cumulative operations metric object
   --------------------------------------------------
*/
cumopsDV = DV_new() ;
DV_init(cumopsDV, nthread, NULL) ;
DV_fill(cumopsDV, 0.0) ;
/*
   -------------------------------
   create the owners map IV object
   -------------------------------
*/
switch ( maptype ) {
case 1 :
   frontOwnersIV = ETree_wrapMap(frontETree, type, 
                                 symmetryflag, cumopsDV) ;
   break ;
case 2 :
   frontOwnersIV = ETree_balancedMap(frontETree, type,
                                 symmetryflag, cumopsDV) ;
   break ;
case 3 :
   frontOwnersIV = ETree_subtreeSubsetMap(frontETree, type,
                                 symmetryflag, cumopsDV) ;
   break ;
case 4 :
   frontOwnersIV = ETree_ddMap(frontETree, type,
                               symmetryflag, cumopsDV, cutoff) ;
   break ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n totalOps = %.0f", DV_sum(cumopsDV)) ;
   DVscale(DV_size(cumopsDV), DV_entries(cumopsDV),
            nthread/DV_sum(cumopsDV)) ;
   fprintf(msgFile, "\n\n cumopsDV") ;
   DV_writeForHumanEye(cumopsDV, msgFile) ;
}
DV_free(cumopsDV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n frontOwnersIV") ;
   IV_writeForHumanEye(frontOwnersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   initialize the FrontMtx object
   -------------------------------
*/
MARKTIME(t1) ;
frontmtx = FrontMtx_new() ;
mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(mtxmanager, LOCK_IN_PROCESS, 0) ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL,
              type, symmetryflag, sparsityflag, pivotingflag,
              LOCK_IN_PROCESS, 0, NULL, mtxmanager, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the front matrix",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n nentD  = %d, nentL = %d, nentU = %d",
        frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n front matrix initialized") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------
   factor the matrix
   -----------------
*/
nzf       = ETree_nFactorEntries(frontETree, symmetryflag) ;
factorops = ETree_nFactorOps(frontETree, type, symmetryflag) ;
fprintf(msgFile, 
        "\n %d factor entries, %.0f factor ops, %8.3f ratio",
        nzf, factorops, factorops/nzf) ;
IVzero(16, stats) ;
DVzero(11, cpus) ;
fprintf(msgFile, "\n\n starting the parallel factor") ;
fflush(msgFile) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, LOCK_IN_PROCESS, 1) ;
MARKTIME(t1) ;
rootchv = FrontMtx_MT_factorInpMtx(frontmtx, mtxA, tau, droptol, 
                                 chvmanager, frontOwnersIV, lookahead,
                                 &error, cpus, stats, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : factor matrix, %8.3f mflops",
        t2 - t1, 1.e-6*factorops/(t2-t1)) ;
if ( rootchv != NULL ) {
   fprintf(msgFile, "\n\n factorization did not complete") ;
   for ( chv = rootchv ; chv != NULL ; chv = chv->next ) {
      fprintf(stdout, "\n chv %d, nD = %d, nL = %d, nU = %d",
              chv->id, chv->nD, chv->nL, chv->nU) ;
   }
}
if ( error >= 0 ) {
   fprintf(msgFile, "\n\n fatal error encountered at front %d", error) ;
   exit(-1) ;
}
fprintf(msgFile,
        "\n %8d pivots, %8d pivot tests, %8d delayed rows and columns",
        stats[0], stats[1], stats[2]) ;
if ( frontmtx->rowadjIVL != NULL ) {
   fprintf(msgFile,
           "\n %d entries in rowadjIVL", frontmtx->rowadjIVL->tsize) ;
}
if ( frontmtx->coladjIVL != NULL ) {
   fprintf(msgFile,
           "\n %d entries in coladjIVL", frontmtx->coladjIVL->tsize) ;
}
if ( frontmtx->upperblockIVL != NULL ) {
   fprintf(msgFile,
           "\n %d fronts, %d entries in upperblockIVL",
           frontmtx->nfront, frontmtx->upperblockIVL->tsize) ;
}
if ( frontmtx->lowerblockIVL != NULL ) {
   fprintf(msgFile,
           ", %d entries in lowerblockIVL",
           frontmtx->lowerblockIVL->tsize) ;
}
fprintf(msgFile,
        "\n %d entries in D, %d entries in L, %d entries in U",
        stats[3], stats[4], stats[5]) ;
fprintf(msgFile, "\n %d locks", frontmtx->nlocks) ;
cputotal = cpus[10] ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    manage working storage  %8.3f %6.2f"
   "\n    initialize/load fronts  %8.3f %6.2f"
   "\n    update fronts           %8.3f %6.2f"
   "\n    aggregate insert        %8.3f %6.2f"
   "\n    aggregate remove/add    %8.3f %6.2f"
   "\n    assemble postponed data %8.3f %6.2f"
   "\n    factor fronts           %8.3f %6.2f"
   "\n    extract postponed data  %8.3f %6.2f"
   "\n    store factor entries    %8.3f %6.2f"
   "\n    miscellaneous           %8.3f %6.2f"
   "\n    total time              %8.3f",
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal, 
   cpus[5], 100.*cpus[5]/cputotal, 
   cpus[6], 100.*cpus[6]/cputotal, 
   cpus[7], 100.*cpus[7]/cputotal, 
   cpus[8], 100.*cpus[8]/cputotal, 
   cpus[9], 100.*cpus[9]/cputotal, cputotal) ;
}
SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
ChvManager_writeForHumanEye(chvmanager, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
}
/*
   ------------------------------
   post-process the factor matrix
   ------------------------------
*/
MARKTIME(t1) ;
FrontMtx_postProcess(frontmtx, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : post-process the matrix", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front factor matrix after post-processing") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
}
fprintf(msgFile, "\n\n after post-processing") ;
SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
/*
   ----------------
   solve the system
   ----------------
*/
neqns = mtxB->nrow ;
nrhs  = mtxB->ncol ;
mtxZ  = DenseMtx_new() ;
DenseMtx_init(mtxZ, type, 0, 0, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(mtxZ) ;
if ( type == SPOOLES_REAL ) {
   nops = frontmtx->nentD + 2*frontmtx->nentU ;
   if ( frontmtx->symmetryflag == 2 ) {
      nops += 2*frontmtx->nentL ;
   } else {
      nops += 2*frontmtx->nentU ;
   }
} if ( type == SPOOLES_COMPLEX ) {
   nops = 8*frontmtx->nentD + 8*frontmtx->nentU ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      nops += 8*frontmtx->nentL ;
   } else {
      nops += 8*frontmtx->nentU ;
   }
}
nops *= nrhs ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rhs") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(stdout) ;
}
DVzero(6, cpus) ;
MARKTIME(t1) ;
FrontMtx_solve(frontmtx, mtxZ, mtxB, mtxmanager, 
               cpus, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : solve the system, %.3f mflops",
        t2 - t1, 1.e-6*nops/(t2 - t1)) ;
cputotal = cpus[5] ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    set up solves               %8.3f %6.2f"
   "\n    load rhs and store solution %8.3f %6.2f"
   "\n    forward solve               %8.3f %6.2f"
   "\n    diagonal solve              %8.3f %6.2f"
   "\n    backward solve              %8.3f %6.2f"
   "\n    total time              %8.3f",
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal, cputotal) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n computed solution") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(stdout) ;
}
DenseMtx_sub(mtxZ, mtxX) ;
fprintf(msgFile, "\n\n serial maxabs error   = %12.4e",
        DenseMtx_maxabs(mtxZ)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n error") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
fprintf(msgFile, "\n\n after solve") ;
SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
/*
   --------------------------------
   now solve the system in parallel
   --------------------------------
*/
solvemap = SolveMap_new() ;
if (  FRONTMTX_IS_NONSYMMETRIC(frontmtx)
   && FRONTMTX_IS_PIVOTING(frontmtx) ) {
   SolveMap_ddMap(solvemap, SPOOLES_NONSYMMETRIC, 
                  FrontMtx_upperBlockIVL(frontmtx),
                  FrontMtx_lowerBlockIVL(frontmtx),
                  nthread, frontOwnersIV, frontmtx->tree,
                  seed, msglvl, msgFile) ;
} else {
   SolveMap_ddMap(solvemap, SPOOLES_SYMMETRIC, 
                  FrontMtx_upperBlockIVL(frontmtx), NULL,
                  nthread, frontOwnersIV, frontmtx->tree,
                  seed, msglvl, msgFile) ;
}
fprintf(msgFile, "\n solve map created") ;
fflush(msgFile) ;
DVzero(6, cpus) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rhs") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(stdout) ;
}
DenseMtx_zero(mtxZ) ;
MARKTIME(t1) ;
FrontMtx_MT_solve(frontmtx, mtxZ, mtxB, mtxmanager, solvemap,
                  cpus, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : solve the system, %.3f mflops",
        t2 - t1, 1.e-6*nops/(t2 - t1)) ;
cputotal = cpus[5] ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    set up solves               %8.3f %6.2f"
   "\n    load rhs and store solution %8.3f %6.2f"
   "\n    forward solve               %8.3f %6.2f"
   "\n    diagonal solve              %8.3f %6.2f"
   "\n    backward solve              %8.3f %6.2f"
   "\n    total time              %8.3f",
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal, cputotal) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n computed solution") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(stdout) ;
}
DenseMtx_sub(mtxZ, mtxX) ;
fprintf(msgFile, "\n\n parallel maxabs error   = %12.4e",
        DenseMtx_maxabs(mtxZ)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n error") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
fprintf(msgFile, "\n\n after solve") ;
SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
FrontMtx_free(frontmtx) ;
ETree_free(frontETree) ;
IV_free(frontOwnersIV) ;
IVL_free(symbfacIVL) ;
InpMtx_free(mtxA) ;
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxB) ;
DenseMtx_free(mtxZ) ;
ChvManager_free(chvmanager) ;
SubMtxManager_free(mtxmanager) ;
SolveMap_free(solvemap) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
