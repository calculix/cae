/*  testGrid.c  */

#include "../FrontMtx.h"
#include "../../Drand.h"
#include "../../SymbFac.h"
#include "../../timings.h"
#include "../../misc.h"

/*--------------------------------------------------------------------*/
void mkNDlinsys ( int n1, int n2, int n3, int maxzeros, int maxsize,
   int type, int symmetryflag, int nrhs, int seed, int msglvl,
   FILE *msgFile, ETree **pfrontETree, IVL **psymbfacIVL,
   InpMtx **pmtxA, DenseMtx **pmtxX, DenseMtx **pmtxB ) ;
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
double          cputotal, droptol, factorops ;
double          cpus[9] ;
Drand           drand ;
double          nops, tau, t1, t2   ;
ETree           *frontETree   ;
FILE            *msgFile ;
int             error, lockflag, maxsize, maxzeros, msglvl, neqns, 
                n1, n2, n3, nrhs, nzf, pivotingflag, 
                seed, sparsityflag, symmetryflag, type ;
int             stats[6] ;
IVL             *symbfacIVL ;

if ( argc != 17 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile n1 n2 n3 maxzeros maxsize" 
"\n         seed type symmetryflag sparsityflag "
"\n         pivotingflag tau droptol lockflag nrhs"
"\n    msglvl   -- message level"
"\n    msgFile  -- message file"
"\n    n1       -- number of grid points in the first direction"
"\n    n2       -- number of grid points in the second direction"
"\n    n3       -- number of grid points in the third direction"
"\n    maxzeros -- max number of zeroes in a front"
"\n    maxsize  -- max number of internal nodes in a front"
"\n    seed     -- random number seed"
"\n    type     -- type of entries"
"\n       1 --> real"
"\n       2 --> complex"
"\n    symmetryflag -- symmetry flag"
"\n       0 --> symmetric "
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
"\n    lockflag -- flag to specify lock status"
"\n       0 --> mutex lock is not allocated or initialized"
"\n       1 --> mutex lock is allocated and it can synchronize"
"\n             only threads in this process."
"\n       2 --> mutex lock is allocated and it can synchronize"
"\n             only threads in this and other processes."
"\n    nrhs     -- # of right hand sides"
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
lockflag      = atoi(argv[15]) ;
nrhs          = atoi(argv[16]) ;
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
        "\n lockflag      -- %d" 
        "\n nrhs          -- %d" 
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, maxzeros, maxsize,
        seed, type, symmetryflag, sparsityflag, pivotingflag, 
        tau, droptol, lockflag, nrhs) ;
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
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n mtxA") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
   fprintf(msgFile, "\n mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fprintf(msgFile, "\n mtxB") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(msgFile) ;
}
/*
fprintf(msgFile, "\n neqns = %d ;", n1*n2*n3) ;
fprintf(msgFile, "\n nrhs = %d ;", nrhs) ;
fprintf(msgFile, "\n A = zeros(neqns, neqns) ;") ;
fprintf(msgFile, "\n X = zeros(neqns, nrhs) ;") ;
fprintf(msgFile, "\n B = zeros(neqns, nrhs) ;") ;
InpMtx_writeForMatlab(mtxA, "A", msgFile) ;
DenseMtx_writeForMatlab(mtxX, "X", msgFile) ;
DenseMtx_writeForMatlab(mtxB, "B", msgFile) ;
{
int      *ivec1 = InpMtx_ivec1(mtxA) ;
int      *ivec2 = InpMtx_ivec2(mtxA) ;
double   *dvec = InpMtx_dvec(mtxA) ;
int      ichv, ii, col, offset, row, nent = InpMtx_nent(mtxA) ;
fprintf(msgFile, "\n coordType = %d", mtxA->coordType) ;
fprintf(msgFile, "\n start of matrix output file") ;
fprintf(msgFile, "\n %d %d %d", n1*n2*n3, n1*n2*n3, nent) ;
for ( ii = 0 ; ii < nent ; ii++ ) {
   ichv = ivec1[ii] ; 
   if ( (offset = ivec2[ii]) >= 0 ) {
      row = ichv, col = row + offset ;
   } else {
      col = ichv, row = col - offset ;
   }
   fprintf(msgFile, "\n %d %d %24.16e %24.16e",
           row, col, dvec[2*ii], dvec[2*ii+1]) ;
}
}
{
int      ii, jj ;
double   imag, real ;
fprintf(msgFile, "\n start of rhs output file") ;
fprintf(msgFile, "\n %d %d", n1*n2*n3, nrhs) ;
for ( ii = 0 ; ii < n1*n2*n3 ; ii++ ) {
   fprintf(msgFile, "\n %d ", ii) ;
   for ( jj = 0 ; jj < nrhs ; jj++ ) {
      DenseMtx_complexEntry(mtxB, ii, jj, &real, &imag) ;
      fprintf(msgFile, " %24.16e %24.16e", real, imag) ;
   }
}
}
*/
/*
   ------------------------------
   initialize the FrontMtx object
   ------------------------------
*/
MARKTIME(t1) ;
frontmtx   = FrontMtx_new() ;
mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(mtxmanager, lockflag, 0) ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL,
              type, symmetryflag, sparsityflag, pivotingflag,
              lockflag, 0, NULL, mtxmanager, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : initialize the front matrix",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile,
           "\n nendD  = %d, nentL = %d, nentU = %d",
           frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n front matrix initialized") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
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
IVzero(6, stats) ;
DVzero(9, cpus) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, lockflag, 1) ;
MARKTIME(t1) ;
rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, droptol, 
                                chvmanager, &error, cpus, 
                                stats, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : factor matrix, %8.3f mflops",
        t2 - t1, 1.e-6*factorops/(t2-t1)) ;
if ( rootchv != NULL ) {
   fprintf(msgFile, "\n\n factorization did not complete") ;
   for ( chv = rootchv ; chv != NULL ; chv = chv->next ) {
      fprintf(stdout, "\n chv %d, nD = %d, nL = %d, nU = %d",
              chv->id, chv->nD, chv->nL, chv->nU) ;
   }
}
if ( error >= 0 ) {
   fprintf(msgFile, "\n\n error encountered at front %d\n", error) ;
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
           ", %d entries in coladjIVL", frontmtx->coladjIVL->tsize) ;
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
cputotal = cpus[8] ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    initialize fronts       %8.3f %6.2f"
   "\n    load original entries   %8.3f %6.2f"
   "\n    update fronts           %8.3f %6.2f"
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
   cpus[7], 100.*cpus[7]/cputotal, cputotal) ;
}
SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
ChvManager_writeForHumanEye(chvmanager, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %% MATLAB file: front factor matrix") ;
   FrontMtx_writeForMatlab(frontmtx, "L", "D", "U", msgFile) ;
}
/*
   ------------------------------
   post-process the factor matrix
   ------------------------------
*/
MARKTIME(t1) ;
FrontMtx_postProcess(frontmtx, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : post-process the matrix", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front factor matrix after post-processing") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
}
fprintf(msgFile, "\n\n after post-processing") ;
SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
/*
  code to test out the IO methods.
  write the matrix to a file, free it,
  then read it back in.
  note: formatted files do not have much accuracy.
*/
/*
FrontMtx_writeToFile(frontmtx, "temp.frontmtxb") ;
FrontMtx_free(frontmtx) ;
frontmtx = FrontMtx_new() ;
FrontMtx_readFromFile(frontmtx, "temp.frontmtxb") ;
frontmtx->manager = mtxmanager ;
FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
*/
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
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      nops += 2*frontmtx->nentL ;
   } else {
      nops += 2*frontmtx->nentU ;
   }
} else if ( type == SPOOLES_COMPLEX ) {
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
fprintf(msgFile, "\n\n CPU %8.3f : solve the system, %.3f mflops",
        t2 - t1, 1.e-6*nops/(t2 - t1)) ;
cputotal = t2 - t1 ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    set up solves               %8.3f %6.2f"
   "\n    load rhs and store solution %8.3f %6.2f"
   "\n    forward solve               %8.3f %6.2f"
   "\n    diagonal solve              %8.3f %6.2f"
   "\n    backward solve              %8.3f %6.2f"
   "\n    total time                  %8.3f",
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal, cputotal) ;
}
/*
fprintf(msgFile, "\n Z = zeros(neqns, nrhs) ;") ;
DenseMtx_writeForMatlab(mtxZ, "Z", msgFile) ;
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n computed solution") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
DenseMtx_sub(mtxZ, mtxX) ;
fprintf(msgFile, "\n\n maxabs error = %12.4e", DenseMtx_maxabs(mtxZ)) ;
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
InpMtx_free(mtxA) ;
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxB) ;
DenseMtx_free(mtxZ) ;
FrontMtx_free(frontmtx) ;
ETree_free(frontETree) ;
IVL_free(symbfacIVL) ;
ChvManager_free(chvmanager) ;
SubMtxManager_free(mtxmanager) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
