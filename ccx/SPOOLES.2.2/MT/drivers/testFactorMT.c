/*  testFactorMT.c  */

#include "../spoolesMT.h"
#include "../../FrontMtx.h"
#include "../../Drand.h"
#include "../../SymbFac.h"
#include "../../timings.h"
#include "../../misc.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------------
   test the factor method for a grid matrix
   (1) read in an InpMtx object
   (2) read in an ETree object
   (3) create a solution matrix object
   (4) multiply the solution with the matrix
       to get a right hand side matrix object
   (5) factor the matrix 
   (6) solve the system

   created -- 98sep05, cca
   -----------------------------------------------------
*/
{
char            *etreeFileName, *mtxFileName ;
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxB, *mtxX, *mtxZ ;
double          one[2] = { 1.0, 0.0 }, cpus[11] ;
double          cputotal, cutoff, droptol, factorops ;
Drand           drand ;
double          nops, tau, t1, t2   ;
DV              *cumopsDV ;
ETree           *frontETree   ;
FILE            *msgFile ;
FrontMtx        *frontmtx ;
InpMtx          *mtxA ;
int             error, loc, lookahead, maptype, msglvl, neqns, nrhs, 
                nthread, nzf, pivotingflag, rc, seed, sparsityflag, 
                symmetryflag, type ;
int             stats[20] ;
IV              *frontOwnersIV, *newToOldIV, *oldToNewIV ;
IVL             *symbfacIVL ;
SolveMap        *solvemap ;
SubMtxManager   *mtxmanager ;

if ( argc != 16 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile mtxFile etreeFile"
"\n         seed symmetryflag sparsityflag "
"\n         pivotingflag tau droptol nrhs "
"\n         nthread maptype cutoff lookahead"
"\n    msglvl    -- message level"
"\n    msgFile   -- message file"
"\n    mtxFile   -- file to read in InpMtx matrix object"
"\n    etreeFile -- file to read in ETree front tree object"
"\n    seed      -- random number seed"
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
"\n    nrhs     -- # of right hand sides"
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
mtxFileName   = argv[3] ;
etreeFileName = argv[4] ;
seed          = atoi(argv[5]) ;
symmetryflag  = atoi(argv[6]) ;
sparsityflag  = atoi(argv[7]) ;
pivotingflag  = atoi(argv[8]) ;
tau           = atof(argv[9]) ;
droptol       = atof(argv[10]) ;
nrhs          = atoi(argv[11]) ;
nthread       = atoi(argv[12]) ;
maptype       = atoi(argv[13]) ;
cutoff        = atof(argv[14]) ;
lookahead     = atoi(argv[15]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n mtxFileName   -- %s"
        "\n etreeFileName -- %s"
        "\n seed          -- %d" 
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
        argv[0], msglvl, argv[2], mtxFileName, etreeFileName,
        seed, symmetryflag, sparsityflag, pivotingflag, 
        tau, droptol, nrhs, nthread, maptype, cutoff, lookahead) ;
fflush(msgFile) ;
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
   -------------------------
   read in the InpMtx object
   -------------------------
*/
if ( strcmp(mtxFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
mtxA = InpMtx_new() ;
MARKTIME(t1) ;
rc = InpMtx_readFromFile(mtxA, mtxFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in mtxA from file %s",
        t2 - t1, mtxFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from InpMtx_readFromFile(%p,%s)",
           rc, mtxA, mtxFileName) ;
   exit(-1) ;
}
type = mtxA->inputMode ;
neqns = 1 + IVmax(mtxA->nent, InpMtx_ivec1(mtxA), &loc) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after reading InpMtx object from file %s",
           mtxFileName) ;
   if ( msglvl == 2 ) {
      InpMtx_writeStats(mtxA, msgFile) ;
   } else {
      InpMtx_writeForHumanEye(mtxA, msgFile) ;
   }
   fflush(msgFile) ;
}
if ( INPMTX_IS_BY_ROWS(mtxA) ) {
   fprintf(msgFile, "\n matrix coordinate type is rows") ;
} else if ( INPMTX_IS_BY_COLUMNS(mtxA) ) {
   fprintf(msgFile, "\n matrix coordinate type is columns") ;
} else if ( INPMTX_IS_BY_CHEVRONS(mtxA) ) {
   fprintf(msgFile, "\n matrix coordinate type is chevrons") ;
} else {
   fprintf(msgFile, "\n\n, error, bad coordinate type") ;
   exit(-1) ;
}
if ( INPMTX_IS_RAW_DATA(mtxA) ) {
   fprintf(msgFile, "\n matrix storage mode is raw data\n") ;
} else if ( INPMTX_IS_SORTED(mtxA) ) {
   fprintf(msgFile, "\n matrix storage mode is sorted\n") ;
} else if ( INPMTX_IS_BY_VECTORS(mtxA) ) {
   fprintf(msgFile, "\n matrix storage mode is by vectors\n") ;
} else {
   fprintf(msgFile, "\n\n, error, bad storage mode") ;
   exit(-1) ;
}
/*
   --------------------------------------------------------
   generate the linear system
   1. generate solution matrix and fill with random numbers
   2. generate rhs matrix and fill with zeros
   3. compute matrix-matrix multiply
   --------------------------------------------------------
*/
MARKTIME(t1) ;
mtxX = DenseMtx_new() ;
DenseMtx_init(mtxX, type, 0, -1, neqns, nrhs, 1, neqns) ;
DenseMtx_fillRandomEntries(mtxX, &drand) ;
mtxB = DenseMtx_new() ;
DenseMtx_init(mtxB, type, 1, -1, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(mtxB) ;
switch ( symmetryflag ) {
case SPOOLES_SYMMETRIC : 
   InpMtx_sym_mmm(mtxA, mtxB, one, mtxX) ;
   break ;
case SPOOLES_HERMITIAN :
   InpMtx_herm_mmm(mtxA, mtxB, one, mtxX) ;
   break ;
case SPOOLES_NONSYMMETRIC :
   InpMtx_nonsym_mmm(mtxA, mtxB, one, mtxX) ;
   break ;
default :
   break ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : set up the solution and rhs ",
        t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n original mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fprintf(msgFile, "\n\n original mtxB") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
if ( strcmp(etreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
frontETree = ETree_new() ;
MARKTIME(t1) ;
rc = ETree_readFromFile(frontETree, etreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in frontETree from file %s",
        t2 - t1, etreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
           rc, frontETree, etreeFileName) ;
   exit(-1) ;
}
ETree_leftJustify(frontETree) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after reading ETree object from file %s",
           etreeFileName) ;
   if ( msglvl == 2 ) {
      ETree_writeStats(frontETree, msgFile) ;
   } else {
      ETree_writeForHumanEye(frontETree, msgFile) ;
   }
}
fflush(msgFile) ;
/*
   --------------------------------------------------
   get the permutations, permute the matrix and the 
   front tree, and compute the symbolic factorization
   --------------------------------------------------
*/
MARKTIME(t1) ;
oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : get permutations", t2 - t1) ;
MARKTIME(t1) ;
ETree_permuteVertices(frontETree, oldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute front tree", t2 - t1) ;
MARKTIME(t1) ;
InpMtx_permute(mtxA, IV_entries(oldToNewIV), IV_entries(oldToNewIV)) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute mtxA", t2 - t1) ;
MARKTIME(t1) ;
InpMtx_mapToUpperTriangle(mtxA) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : map to upper triangle", t2 - t1) ;
if ( ! INPMTX_IS_BY_CHEVRONS(mtxA) ) {
   MARKTIME(t1) ;
   InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : change coordinate type", t2 - t1) ;
}
if ( INPMTX_IS_RAW_DATA(mtxA) ) {
   MARKTIME(t1) ;
   InpMtx_changeStorageMode(mtxA, INPMTX_SORTED) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : sort entries ", t2 - t1) ;
}
if ( INPMTX_IS_SORTED(mtxA) ) {
   MARKTIME(t1) ;
   InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : convert to vectors ", t2 - t1) ;
}
MARKTIME(t1) ;
symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : symbolic factorization", t2 - t1) ;
MARKTIME(t1) ;
DenseMtx_permuteRows(mtxB, oldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute rhs", t2 - t1) ;
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
   ------------------------------
   initialize the FrontMtx object
   ------------------------------
*/
MARKTIME(t1) ;
frontmtx   = FrontMtx_new() ;
mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(mtxmanager, LOCK_IN_PROCESS, 0) ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL,
              type, symmetryflag, sparsityflag, pivotingflag,
              LOCK_IN_PROCESS, 0, NULL, mtxmanager, msglvl, msgFile) ;
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
IVzero(16, stats) ;
DVzero(11, cpus) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, LOCK_IN_PROCESS, 1) ;
MARKTIME(t1) ;
error = -1 ;
rootchv = FrontMtx_MT_factorInpMtx(frontmtx, mtxA, tau, droptol, 
                                chvmanager, frontOwnersIV, lookahead,
                                &error, cpus, stats, msglvl, msgFile) ;
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
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n computed solution") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
/*
   -------------------------------------------------------------
   permute the computed solution back into the original ordering
   -------------------------------------------------------------
*/
MARKTIME(t1) ;
DenseMtx_permuteRows(mtxZ, newToOldIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute solution", t2 - t1) ;
/*
   -----------------
   compute the error
   -----------------
*/
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
   ----------------------------
   solve the system in parallel
   ----------------------------
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
neqns = mtxB->nrow ;
nrhs  = mtxB->ncol ;
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
FrontMtx_MT_solve(frontmtx, mtxZ, mtxB, mtxmanager, solvemap,
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
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n computed solution") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(stdout) ;
}
/*
   -------------------------------------------------------------
   permute the computed solution back into the original ordering
   -------------------------------------------------------------
*/
MARKTIME(t1) ;
DenseMtx_permuteRows(mtxZ, newToOldIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute solution", t2 - t1) ;
/*
   -----------------
   compute the error
   -----------------
*/
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
IV_free(oldToNewIV) ;
IV_free(newToOldIV) ;
IV_free(frontOwnersIV) ;
InpMtx_free(mtxA) ;
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxB) ;
DenseMtx_free(mtxZ) ;
FrontMtx_free(frontmtx) ;
ETree_free(frontETree) ;
IVL_free(symbfacIVL) ;
ChvManager_free(chvmanager) ;
SubMtxManager_free(mtxmanager) ;
SolveMap_free(solvemap) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
