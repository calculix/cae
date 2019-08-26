/*  testQRgridMT.c  */

#include "../spoolesMT.h"
#include "../../FrontMtx.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
void mkNDlinsysQR ( int n1, int n2, int n3, int type, int nrhs,
   int seed, int msglvl, FILE *msgFile, ETree **pfrontETree,
   IVL **psymbfacIVL, InpMtx **pmtxA, DenseMtx **pmtxX,
   DenseMtx **pmtxB) ;
/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   test the QR factor method for a FrontMtx object
   on an n1 x n2 x n3 grid

   (1) generate an overdetermined system AX = B
   (2) factor the matrix 
   (3) solve the systems

   created  -- 98may29, cca
   ---------------------------------------------------
*/
{
ChvManager      *chvmanager ;
DenseMtx        *mtxB, *mtxX, *mtxZ ;
double          cputotal, cutoff, factorops ;
double          cpus[10] ;
double          nops, t1, t2 ;
DV              *cumopsDV ;
ETree           *frontETree ;
FILE            *msgFile ;
FrontMtx        *frontmtx ;
InpMtx          *mtxA ;
int             maptype, msglvl, neqns, nrhs, nthread,
                n1, n2, n3, seed, type ;
IV              *frontOwnersIV ;
IVL             *symbfacIVL ;
SolveMap        *solvemap ;
SubMtxManager   *mtxmanager ;

if ( argc != 12 ) {
   fprintf(stdout,
      "\n\n usage : %s msglvl msgFile n1 n2 n3 seed nrhs type "
      "\n         nthread maptype cutoff"
      "\n    msglvl  -- message level"
      "\n    msgFile -- message file"
      "\n    n1      -- # of points in the first direction"
      "\n    n2      -- # of points in the second direction"
      "\n    n3      -- # of points in the third direction"
      "\n    seed    -- random number seed"
      "\n    nrhs    -- # of right hand sides"
      "\n    type    -- type of linear system"
      "\n      1 -- real"
      "\n      2 -- complex"
      "\n    nthread -- number of threads"
      "\n    maptype -- type of map from fronts to threads"
      "\n       1 --> wrap map"
      "\n       2 --> balanced map via a post-order traversal"
      "\n       3 --> subtree-subset map"
      "\n       4 --> domain decomposition map"
      "\n    cutoff  -- cutoff used for domain decomposition map"
      "\n       0 <= cutoff <= 1 used to define the multisector"
      "\n", argv[0]) ;
   return(0) ;
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
n1      = atoi(argv[3]) ;
n2      = atoi(argv[4]) ;
n3      = atoi(argv[5]) ;
seed    = atoi(argv[6]) ;
nrhs    = atoi(argv[7]) ;
type    = atoi(argv[8]) ;
nthread = atoi(argv[8]) ;
maptype = atoi(argv[8]) ;
cutoff  = atoi(argv[8]) ;
fprintf(msgFile,
        "\n %s "
        "\n msglvl  -- %d"
        "\n msgFile -- %s"
        "\n n1      -- %d"
        "\n n2      -- %d"
        "\n n3      -- %d"
        "\n seed    -- %d"
        "\n nrhs    -- %d"
        "\n type    -- %d"
        "\n nthread -- %d"
        "\n maptype -- %d"
        "\n cutoff  -- %f"
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, seed, nrhs, type,
        nthread, maptype, cutoff) ;
fflush(msgFile) ;
neqns = n1*n2*n3 ;
if ( type != SPOOLES_REAL && type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error, type must be real or complex") ;
   exit(-1) ;
}
/*
   ------------------------------------------
   generate the A X = B overdetermined system
   ------------------------------------------
*/
mkNDlinsysQR(n1, n2, n3, type, nrhs, seed, msglvl, msgFile,
             &frontETree, &symbfacIVL, &mtxA, &mtxX, &mtxB) ;
fprintf(msgFile, "\n\n %d entries in A", InpMtx_nent(mtxA)) ;
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
                                 SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 2 :
   frontOwnersIV = ETree_balancedMap(frontETree, type,
                                 SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 3 :
   frontOwnersIV = ETree_subtreeSubsetMap(frontETree, type,
                                 SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 4 :
   frontOwnersIV = ETree_ddMap(frontETree, type,
                                 SPOOLES_SYMMETRIC, cumopsDV, cutoff) ;
   break ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n totalOps = %.0f", DV_sum(cumopsDV)) ;
   DVscale(DV_size(cumopsDV), DV_entries(cumopsDV),
            nthread/DV_sum(cumopsDV)) ;
   fprintf(msgFile, "\n\n cumopsDV") ;
   DV_writeForHumanEye(cumopsDV, msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n frontOwnersIV") ;
   IV_writeForHumanEye(frontOwnersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------
   initialize the FrontMtx object
   ------------------------------
*/
MARKTIME(t1) ;
mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(mtxmanager, LOCK_IN_PROCESS, 0) ;
frontmtx = FrontMtx_new() ;
if ( type == SPOOLES_REAL ) {
   FrontMtx_init(frontmtx, frontETree, symbfacIVL, type,
                 SPOOLES_SYMMETRIC, FRONTMTX_DENSE_FRONTS,
                 SPOOLES_NO_PIVOTING, LOCK_IN_PROCESS,
                 0, NULL, mtxmanager, msglvl, msgFile) ;
} else if ( type == SPOOLES_COMPLEX ) {
   FrontMtx_init(frontmtx, frontETree, symbfacIVL, type,
                 SPOOLES_HERMITIAN, FRONTMTX_DENSE_FRONTS,
                 SPOOLES_NO_PIVOTING, LOCK_IN_PROCESS,
                 0, NULL, mtxmanager, msglvl, msgFile) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the front matrix",
        t2 - t1) ;
fprintf(msgFile, "\n\n %d entries in D, %d entries in U",
        frontmtx->nentD, frontmtx->nentU) ;
/*
   -----------------
   factor the matrix
   -----------------
*/
DVzero(10, cpus) ;
InpMtx_changeCoordType(mtxA, INPMTX_BY_ROWS) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, LOCK_IN_PROCESS, 1) ;
MARKTIME(t1) ;
FrontMtx_MT_QR_factor(frontmtx, mtxA, chvmanager, frontOwnersIV, 
                      cpus, &factorops, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n CPU %8.3f : FrontMtx_MT_QR_factor, %.0f ops, %.2f mflops",
        t2 - t1, factorops, 1.e-6*factorops/(t2-t1)) ;
cputotal = cpus[6] ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile, "\n"
   "\n                              CPU    %%"
   "\n    setup factorization  %8.3f %6.2f"
   "\n    setup fronts         %8.3f %6.2f"
   "\n    factor fronts        %8.3f %6.2f"
   "\n    store factor         %8.3f %6.2f"
   "\n    store update         %8.3f %6.2f"
   "\n    miscellaneous        %8.3f %6.2f"
   "\n    total time           %8.3f"
   "\n    wall clock time      %8.3f",
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal,
   cpus[5], 100.*cpus[5]/cputotal, 
   cpus[6], t2 - t1) ;
}
fprintf(msgFile, "\n\n ChvManager statistics") ;
ChvManager_writeForHumanEye(chvmanager, msgFile) ;
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
   ----------------
   solve the system
   ----------------
*/
mtxZ = DenseMtx_new() ;
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
DVzero(7, cpus) ;
MARKTIME(t1) ;
FrontMtx_QR_solve(frontmtx, mtxA, mtxZ, mtxB, mtxmanager,
                  cpus, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, 
     "\n\n CPU %8.3f : serial solve, %d rhs, %.0f ops, %.3f mflops",
     t2 - t1, nrhs, nops, 1.e-6*nops/(t2 - t1)) ;
cputotal = t2 - t1 ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile, "\n"
   "\n                                        CPU    %%"
   "\n    A^TB matrix-matrix multiply    %8.3f %6.2f"
   "\n    total solve time               %8.3f %6.2f"
   "\n       set up solves                 %8.3f %6.2f"
   "\n       load rhs and store solution   %8.3f %6.2f"
   "\n       forward solve                 %8.3f %6.2f"
   "\n       diagonal solve                %8.3f %6.2f"
   "\n       backward solve                %8.3f %6.2f"
   "\n    total QR solve time            %8.3f",
   cpus[6], 100.*cpus[6]/cputotal,
   cpus[5], 100.*cpus[5]/cputotal, 
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal,
   cputotal) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n serial solve computed solution") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
/*
   -----------------
   compute the error
   -----------------
*/
DenseMtx_sub(mtxZ, mtxX) ;
fprintf(msgFile, "\n\n serial solve: maxabs error = %12.4e",
DenseMtx_maxabs(mtxZ)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n error") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
fprintf(msgFile, "\n\n after serial solve") ;
SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
/*
   --------------------------------
   now solve the system in parallel
   --------------------------------
*/
solvemap = SolveMap_new() ;
SolveMap_ddMap(solvemap, SPOOLES_SYMMETRIC, 
               FrontMtx_upperBlockIVL(frontmtx), NULL,
               nthread, frontOwnersIV, frontmtx->tree,
               seed, msglvl, msgFile) ;
fprintf(msgFile, "\n solve map created") ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rhs") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(stdout) ;
}
DenseMtx_zero(mtxZ) ;
MARKTIME(t1) ;
FrontMtx_MT_QR_solve(frontmtx, mtxA, mtxZ, mtxB, mtxmanager,
                     solvemap, cpus, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, 
     "\n\n CPU %8.3f : parallel solve, %d rhs, %.0f ops, %.3f mflops",
     t2 - t1, nrhs, nops, 1.e-6*nops/(t2 - t1)) ;
cputotal = t2 - t1 ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile, "\n"
   "\n                                        CPU    %%"
   "\n    A^TB matrix-matrix multiply    %8.3f %6.2f"
   "\n    total solve time               %8.3f %6.2f"
   "\n       set up solves                 %8.3f %6.2f"
   "\n       load rhs and store solution   %8.3f %6.2f"
   "\n       forward solve                 %8.3f %6.2f"
   "\n       diagonal solve                %8.3f %6.2f"
   "\n       backward solve                %8.3f %6.2f"
   "\n    total QR solve time            %8.3f",
   cpus[6], 100.*cpus[6]/cputotal,
   cpus[5], 100.*cpus[5]/cputotal, 
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal,
   cputotal) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n computed solution from parallel solve") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
/*
   -----------------
   compute the error
   -----------------
*/
DenseMtx_sub(mtxZ, mtxX) ;
fprintf(msgFile, "\n\n parallel solve: maxabs error = %12.4e",
DenseMtx_maxabs(mtxZ)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n error") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
fprintf(msgFile, "\n\n after parallel solve") ;
SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
InpMtx_free(mtxA) ;
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxZ) ;
DenseMtx_free(mtxB) ;
FrontMtx_free(frontmtx) ;
IVL_free(symbfacIVL) ;
ETree_free(frontETree) ;
SubMtxManager_free(mtxmanager) ;
ChvManager_free(chvmanager) ;
DV_free(cumopsDV) ;
IV_free(frontOwnersIV) ;
SolveMap_free(solvemap) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
