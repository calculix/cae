/*  testGridMPI.c  */

#include "../spoolesMPI.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
void mkNDlinsys ( int n1, int n2, int n3, int maxzeros, int maxsize,
   int type, int symmetryflag, int nrhs, int seed, int msglvl,
   FILE *msgFile, ETree **pfrontETree, IVL **psymbfacIVL,
   InpMtx **pmtxA, DenseMtx **pmtxX, DenseMtx **pmtxB ) ;
/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
*/
{
Chv             *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxB, *mtxkeep, *mtxX, *mtxZ ;
double          cputotal, cutoff, droptol, error, factorops, 
                solvecpu, solveops, tau, terror, t1, t2 ;
double          cpus[14], tcpus[14] ;
DV              *cumopsDV ;
ETree           *frontETree ;
FILE            *msgFile ;
FrontMtx        *frontmtx ;
InpMtx          *mtxA, *mtxAkeep ;
int             factorerror, lookahead, maptype, maxsize, maxzeros, 
                msglvl, myid, neqns, nproc, nrhs, nrow, nzf, n1, n2, 
                n3, pivotingflag, seed, sparsityflag, symmetryflag, 
                tag, type ;
int             *rowindX, *rowindZ ;
int             stats[17], tstats[17] ;
IV              *frontOwnersIV, *vtxmapIV ;
IVL             *symbfacIVL ;
SolveMap        *solvemap ;
SubMtxManager   *factormanager, *solvemanager ;
/*
   ---------------------------------------------------------------
   find out the identity of this process and the number of process
   ---------------------------------------------------------------
*/
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
fprintf(stdout, "\n process %d of %d, argc = %d", myid, nproc, argc) ;
fflush(stdout) ;
if ( argc != 19 ) {
   fprintf(stdout, 
      "\n\n usage : %s "
      "\n         msglvl msgFile n1 n2 n3 maxzeros maxsize seed "
      "\n         type symmetryflag sparsityflag pivotingflag "
      "\n         tau droptol lookahead nrhs maptype cutoff"
      "\n    msglvl       -- message level"
      "\n    msgFile      -- message file"
      "\n    n1           -- # of points in first dimension"
      "\n    n2           -- # of points in second dimension"
      "\n    n3           -- # of points in third dimension"
      "\n    maxzeros     -- max # of zeros in a front"
      "\n    maxsize      -- max # of internal vertices in a front"
      "\n    seed         -- random number seed"
      "\n    type         -- type of entries"
      "\n       1 --> real "
      "\n       2 --> complex "
      "\n    symmetryflag -- symmetry flag"
      "\n       0 --> symmetric "
      "\n       1 --> hermitian "
      "\n       2 --> nonsymmetric "
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
      "\n    lookahead -- parameter to schedule computation"
      "\n       0 --> no lookahead"
      "\n       nonzero --> look up the tree this many steps"
      "\n    nrhs -- number of right hand sides"
      "\n    maptype -- type of factorization map"
      "\n       1 --> wrap map"
      "\n       2 --> balanced map via post-order traversal"
      "\n       3 --> subtree-subset map"
      "\n       4 --> domain decomposition map"
      "\n    cutoff -- used when maptype = 4"
      "\n       0 < cutoff < 1 to define multisector"
      "\n       try cutoff = 1/nproc or 1/(2*nproc)"
      "\n", argv[0]) ;
   { int ii ;
   fprintf(stdout, "\n\n argc = %d", argc) ;
   for ( ii = 0 ; ii < argc ; ii++ ) {
      fprintf(stdout, "\n arg %d = %s", ii, argv[ii]) ;
   }
   }
   MPI_Finalize() ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else {
   int    length = strlen(argv[2]) + 1 + 4 ;
   char   *buffer = CVinit(length, '\0') ;
   sprintf(buffer, "%s.%d", argv[2], myid) ;
   if ( (msgFile = fopen(buffer, "w")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], buffer) ;
      return(-1) ;
   }
   CVfree(buffer) ;
}
n1           = atoi(argv[ 3]) ;
n2           = atoi(argv[ 4]) ;
n3           = atoi(argv[ 5]) ;
maxzeros     = atoi(argv[ 6]) ;
maxsize      = atoi(argv[ 7]) ;
seed         = atoi(argv[ 8]) ;
type         = atoi(argv[ 9]) ;
symmetryflag = atoi(argv[10]) ;
sparsityflag = atoi(argv[11]) ;
pivotingflag = atoi(argv[12]) ;
tau          = atof(argv[13]) ;
droptol      = atof(argv[14]) ;
lookahead    = atoi(argv[15]) ;
nrhs         = atoi(argv[16]) ;
maptype      = atoi(argv[17]) ;
cutoff       = atof(argv[18]) ;
fprintf(msgFile, "\n %s :"
        "\n   msglvl       = %d"
        "\n   msgFile      = %s"
        "\n   n1           = %d"
        "\n   n2           = %d"
        "\n   n3           = %d"
        "\n   maxzeros     = %d"
        "\n   maxsize      = %d"
        "\n   seed         = %d"
        "\n   type         = %d"
        "\n   symmetryflag = %d"
        "\n   sparsityflag = %d"
        "\n   pivotingflag = %d"
        "\n   tau          = %f"
        "\n   droptol      = %f"
        "\n   lookahead    = %d"
        "\n   nrhs         = %d"
        "\n   maptype      = %d"
        "\n   cutoff       = %f"
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, maxzeros, maxsize, seed,
        type, symmetryflag, sparsityflag, pivotingflag, tau, droptol,
        lookahead, nrhs, maptype, cutoff) ;
fflush(msgFile) ;
neqns = n1*n2*n3 ;
/*
   ---------------------------------------------------------------
   find out the identity of this process and the number of process
   ---------------------------------------------------------------
*/
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
if ( myid == 0 ) {
/*
   --------------------------
   generate the linear system
   --------------------------
*/
   mkNDlinsys(n1, n2, n3, maxzeros, maxsize, type, symmetryflag, nrhs,
              seed, msglvl, msgFile, &frontETree, &symbfacIVL, &mtxA,
              &mtxX, &mtxB) ;
/*
   -------------------------------
   free the symbolic factorization
   -------------------------------
*/
   IVL_free(symbfacIVL) ;
/*
   -------------------------------------------
   send the front tree to the other processors
   -------------------------------------------
*/
   frontETree = ETree_MPI_Bcast(frontETree, 0, 
                                msglvl, msgFile, MPI_COMM_WORLD) ;
} else {
/*
   ---------------------------------------
   receive the front tree from processor 0
   ---------------------------------------
*/
   frontETree = ETree_MPI_Bcast(NULL, 0, 
                                msglvl, msgFile, MPI_COMM_WORLD) ;
/*
   ------------------------------
   create the objects for X and B
   ------------------------------
*/
   mtxX = DenseMtx_new() ;
   mtxB = DenseMtx_new() ;
   DenseMtx_init(mtxX, type, -1, -1, 0, nrhs, 1, 0) ;
   DenseMtx_init(mtxB, type, -1, -1, 0, nrhs, 1, 0) ;
/*
   -----------------------
   create the object for A
   -----------------------
*/
   mtxA = InpMtx_new() ;
   InpMtx_init(mtxA, INPMTX_BY_CHEVRONS, type, 0, 0) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n front ETree object") ;
   ETree_writeForHumanEye(frontETree, msgFile) ;
   fprintf(msgFile, "\n\n mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fprintf(msgFile, "\n\n mtxB") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fprintf(msgFile, "\n\n mtxA") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
}
/*
   ---------------------------------
   generate the owners map IV object
   ---------------------------------
*/
cumopsDV = DV_new() ;
DV_init(cumopsDV, nproc, NULL) ;
DV_fill(cumopsDV, 0.0) ;
switch ( maptype ) {
case 1 :
   frontOwnersIV = ETree_wrapMap(frontETree, type, symmetryflag,
                                 cumopsDV) ;
   break ;
case 2 :
   frontOwnersIV = ETree_balancedMap(frontETree, type, symmetryflag,
                                     cumopsDV) ;
   break ;
case 3 :
   frontOwnersIV = ETree_subtreeSubsetMap(frontETree, type, 
                                          symmetryflag, cumopsDV) ;
   break ;
case 4 :
   frontOwnersIV = ETree_ddMap(frontETree, type, symmetryflag,
                               cumopsDV, cutoff) ;
   break ;
default :
   fprintf(stderr, "\n fatal error, maptype = %d", maptype) ;
   MPI_Finalize() ;
   exit(-1) ;
   break ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n owners map object") ;
   IV_writeForHumanEye(frontOwnersIV, msgFile) ;
}
fprintf(msgFile, "\n\n upper bound on load balance = %4.2f\n",
        DV_sum(cumopsDV)/DV_max(cumopsDV)) ;
DV_free(cumopsDV) ;
/*
   ---------------------
   create the vertex map
   ---------------------
*/
vtxmapIV = IV_new() ;
IV_init(vtxmapIV, neqns, NULL) ;
IVgather(neqns, IV_entries(vtxmapIV), IV_entries(frontOwnersIV), 
         ETree_vtxToFront(frontETree)) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n vertex map object") ;
   IV_writeForHumanEye(vtxmapIV, msgFile) ;
}
/*
   ------------------
   split the X matrix
   ------------------
*/
tag = 1 ;
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
mtxkeep = DenseMtx_MPI_splitByRows(mtxX, vtxmapIV, stats, msglvl,
                                   msgFile, tag, MPI_COMM_WORLD) ;
DenseMtx_free(mtxX) ;
mtxX = mtxkeep ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : split mtxX", t2 - t1) ;
fprintf(msgFile, 
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n mtxX object") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
}
/*
   ------------------
   split the B matrix
   ------------------
*/
tag = 1 ;
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
mtxkeep = DenseMtx_MPI_splitByRows(mtxB, vtxmapIV, stats, msglvl,
                                   msgFile, tag, MPI_COMM_WORLD) ;
DenseMtx_free(mtxB) ;
mtxB = mtxkeep ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : split mtxB", t2 - t1) ;
fprintf(msgFile, 
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n mtxB object") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
}
/*
   ------------------
   split the A matrix
   ------------------
*/
tag = 1 ;
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
mtxAkeep = InpMtx_MPI_split(mtxA, vtxmapIV, stats, msglvl,
                            msgFile, tag, MPI_COMM_WORLD) ;
InpMtx_free(mtxA) ;
mtxA = mtxAkeep ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : split mtxA", t2 - t1) ;
fprintf(msgFile, 
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n mtxA object") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
}
/*
   ----------------------------------
   compute the symbolic factorization
   ----------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
symbfacIVL = SymbFac_MPI_initFromInpMtx(frontETree, frontOwnersIV, mtxA,
                          stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : symbolic factorization", t2 - t1) ;
fprintf(msgFile, 
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n symbolic factorization object") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
}
/*
   ----------------------------------
   initialize the front matrix object
   ----------------------------------
*/
factormanager = SubMtxManager_new() ;
SubMtxManager_init(factormanager, 0, 0) ;
MARKTIME(t1) ;
frontmtx = FrontMtx_new() ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL, type, symmetryflag, 
              sparsityflag, pivotingflag, NO_LOCK, myid, frontOwnersIV, 
              factormanager, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the front matrix",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n front matrix initialized") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------
   factor the front matrix
   -----------------------
*/
nzf       = ETree_nFactorEntries(frontETree, symmetryflag) ;
factorops = ETree_nFactorOps(frontETree, type, symmetryflag) ;
IVzero(17, stats) ;
DVzero(14, cpus) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, NO_LOCK, 0) ;
MARKTIME(t1) ;
rootchv = FrontMtx_MPI_factorInpMtx(frontmtx, mtxA, tau, droptol,
                    chvmanager, frontOwnersIV, lookahead, &factorerror,
                    cpus, stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : factor matrix, %8.3f mflops",
        t2 - t1, 1.e-6*factorops/(t2-t1)) ;
if ( factorerror >= 0 ) {
   fprintf(msgFile, "\n\n error found in front %d", factorerror) ;
   exit(-1) ;
}
if ( rootchv != NULL ) {
   Chv   *chv ;
   fprintf(msgFile, "\n\n factorization did not complete") ;
   for ( chv = rootchv ; chv != NULL ; chv = chv->next ) {
      fprintf(stdout, "\n chv %d, nD = %d, nL = %d, nU = %d",
              chv->id, chv->nD, chv->nL, chv->nU) ;
      Chv_writeForHumanEye(chv, msgFile) ;
   }
}
DVzero(14, tcpus) ;
MPI_Reduce((void *) cpus, (void *) tcpus, 14, MPI_DOUBLE,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
cputotal = DVsum(14, tcpus) ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    initialize fronts         %8.3f %6.2f"
   "\n    load fronts               %8.3f %6.2f"
   "\n    update fronts             %8.3f %6.2f"
   "\n    aggregate insert          %8.3f %6.2f"
   "\n    aggregate remove/add      %8.3f %6.2f"
   "\n    assemble postponed data   %8.3f %6.2f"
   "\n    factor fronts             %8.3f %6.2f"
   "\n    extract postponed data    %8.3f %6.2f"
   "\n    store factor entries      %8.3f %6.2f"
   "\n    post initial receives     %8.3f %6.2f"
   "\n    check for recv'd messages %8.3f %6.2f"
   "\n    post initial sends        %8.3f %6.2f"
   "\n    check for sent messages   %8.3f %6.2f"
   "\n    total time                %8.3f",
   tcpus[0], 100.*tcpus[0]/cputotal,
   tcpus[1], 100.*tcpus[1]/cputotal,
   tcpus[2], 100.*tcpus[2]/cputotal,
   tcpus[3], 100.*tcpus[3]/cputotal,
   tcpus[4], 100.*tcpus[4]/cputotal,
   tcpus[5], 100.*tcpus[5]/cputotal,
   tcpus[6], 100.*tcpus[6]/cputotal,
   tcpus[7], 100.*tcpus[7]/cputotal,
   tcpus[8], 100.*tcpus[8]/cputotal,
   tcpus[9], 100.*tcpus[9]/cputotal, 
   tcpus[10], 100.*tcpus[10]/cputotal, 
   tcpus[11], 100.*tcpus[11]/cputotal, 
   tcpus[12], 100.*tcpus[12]/cputotal, cputotal) ;
}
fprintf(msgFile, 
     "\n\n Local statistics"
     "\n    %d pivots, %d pivot tests, %d delayed rows and columns"
     "\n    %d entries in D, %d entries in L, %d entries in U"
     "\n    %d active Chv objects, %d active bytes, %d requested bytes",
     stats[0], stats[1], stats[2], stats[3], stats[4], stats[5],
     stats[14], stats[15], stats[16]) ;
fprintf(msgFile, 
        "\n    local comm  : # sent   bytes sent   # recv   bytes recv"
        "\n    aggregate     %6d %12d   %6d %12d"
        "\n    postponed     %6d %12d   %6d %12d",
        stats[6],  stats[7],  stats[8],  stats[9], 
        stats[10], stats[11], stats[12], stats[13]) ;
IVzero(11, tstats) ;
MPI_Reduce((void *) stats, (void *) tstats, 17, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
     "\n\n Global statistics"
     "\n    %d pivots, %d pivot tests, %d delayed rows and columns"
     "\n    %d entries in D, %d entries in L, %d entries in U"
     "\n    %d active Chv objects, %d active bytes, %d requested bytes",
        tstats[0], tstats[1], tstats[2], 
        tstats[3], tstats[4], tstats[5],
        tstats[14], tstats[15], tstats[16]) ;
   fprintf(msgFile, 
        "\n    global comm  : # sent   bytes sent   # recv   bytes recv"
        "\n    aggregate      %6d %12d   %6d %12d"
        "\n    postponed      %6d %12d   %6d %12d",
        tstats[6],  tstats[7],  tstats[8],  tstats[9], 
        tstats[10], tstats[11], tstats[12], tstats[13]) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      solveops = 2*(tstats[3] + tstats[4] + tstats[5]) ;
   } else {
      solveops = 2*(tstats[3] + 2*tstats[5]) ;
   }
   if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
      solveops *= 4 ;
   }
   solveops *= nrhs ;
   fprintf(msgFile, "\n\n solve operations = %.0f", solveops) ;
}
SubMtxManager_writeForHumanEye(factormanager, msgFile) ;
ChvManager_writeForHumanEye(chvmanager, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------------
   redistribute the right hand side and solution if necessary
   ----------------------------------------------------------
*/
MARKTIME(t1) ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) { 
   IV   *colmapIV ;

   colmapIV = FrontMtx_MPI_colmapIV(frontmtx, frontOwnersIV,
                                    msglvl, msgFile, MPI_COMM_WORLD) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n column map IV object") ;
      IV_writeForHumanEye(colmapIV, msgFile) ;
      fflush(msgFile) ;
   }
   stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
   tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
   MARKTIME(t1) ;
   mtxkeep = DenseMtx_MPI_splitByRows(mtxX, colmapIV, stats, msglvl,
                                      msgFile, tag, MPI_COMM_WORLD) ;
   DenseMtx_free(mtxX) ;
   mtxX = mtxkeep ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n\n CPU %8.3f : solution matrix split ", t2 - t1) ;
   fprintf(msgFile, 
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
   fflush(msgFile) ;
   MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD) ;
   if ( myid == 0 ) {
      fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
              tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
      fflush(msgFile) ;
   }
   stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
   tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      IV   *rowmapIV ;

      rowmapIV = FrontMtx_MPI_rowmapIV(frontmtx, frontOwnersIV, msglvl, 
                                       msgFile, MPI_COMM_WORLD) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n row map IV object") ;
         IV_writeForHumanEye(rowmapIV, msgFile) ;
         fflush(msgFile) ;
      }
      MARKTIME(t1) ;
      mtxkeep = DenseMtx_MPI_splitByRows(mtxB, rowmapIV, stats, msglvl,
                                      msgFile, tag, MPI_COMM_WORLD) ;
      DenseMtx_free(mtxB) ;
      mtxB = mtxkeep ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %8.3f : rhs matrix split ", t2 - t1) ;
      IV_free(rowmapIV) ;
   } else {
      IVfill(4, stats, 0) ;
      MARKTIME(t1) ;
      mtxkeep = DenseMtx_MPI_splitByRows(mtxB, colmapIV, stats, msglvl,
                                      msgFile, tag, MPI_COMM_WORLD) ;
      DenseMtx_free(mtxB) ;
      mtxB = mtxkeep ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %8.3f : rhs matrix split ", t2 - t1) ;
   }
   fprintf(msgFile, 
        "\n\n Redistibute rhs and solution for this process"
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
   fflush(msgFile) ;
   MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD) ;
   if ( myid == 0 ) {
      fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
              tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n solution matrix") ;
      DenseMtx_writeForHumanEye(mtxX, msgFile) ;
      fprintf(msgFile, "\n\n right hand side matrix") ;
      DenseMtx_writeForHumanEye(mtxB, msgFile) ;
      fflush(msgFile) ;
   }
   IV_free(colmapIV) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : solution and rhs set up ", t2-t1);
/*
   -----------------------------------------------------------------
   post-process the matrix and convert to a submatrix representation
   -----------------------------------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
FrontMtx_MPI_postProcess(frontmtx, frontOwnersIV, 
                         stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : post-process the factor matrix",
        t2 - t1) ;
fprintf(msgFile, 
        "\n Post-process communication for this process"
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after post-processing, front factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   create the solve map object
   ---------------------------
*/
MARKTIME(t1) ;
solvemap = SolveMap_new() ;
if (  FRONTMTX_IS_NONSYMMETRIC(frontmtx)
   && FRONTMTX_IS_PIVOTING(frontmtx) ) {
   SolveMap_ddMap(solvemap, SPOOLES_NONSYMMETRIC, 
                  frontmtx->upperblockIVL, frontmtx->lowerblockIVL, 
                  nproc, frontOwnersIV, frontmtx->tree, 
                  seed, msglvl, msgFile);
} else {
   SolveMap_ddMap(solvemap, SPOOLES_SYMMETRIC, frontmtx->upperblockIVL, 
                  NULL, nproc, frontOwnersIV, frontmtx->tree, 
                  seed, msglvl, msgFile) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : solve map created", t2 - t1) ;
if ( msglvl > 2 ) {
   SolveMap_writeForHumanEye(solvemap, msgFile) ;
} else {
   SolveMap_writeStats(solvemap, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------------------
   split the front matrix to conform to the solve map
   --------------------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
FrontMtx_MPI_split(frontmtx, solvemap, 
                   stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : front matrix split", t2 - t1) ;
fprintf(msgFile, 
        "\n Split front matrix communication for this process"
        "\n                 # sent   bytes sent   # recv   bytes recv"
        "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d\n",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after splitting, front factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------
   solve the linear system
   -----------------------
*/
solvemanager = SubMtxManager_new() ;
SubMtxManager_init(solvemanager, 0, 0) ;
IVzero(8,  stats) ;
IVzero(8,  tstats) ;
DVzero(5, cpus) ;
mtxZ = DenseMtx_new() ;
DenseMtx_init(mtxZ, type, 0, 0, mtxX->nrow, mtxX->ncol, 1, mtxX->nrow) ;
DenseMtx_zero(mtxZ) ;
DenseMtx_rowIndices(mtxX, &nrow, &rowindX) ;
DenseMtx_rowIndices(mtxZ, &nrow, &rowindZ) ;
IVcopy(nrow, rowindZ, rowindX) ;
MARKTIME(t1) ;
FrontMtx_MPI_solve(frontmtx, mtxZ, mtxB, solvemanager, solvemap, 
                   cpus, stats, msglvl, msgFile, tag, MPI_COMM_WORLD);
MARKTIME(t2) ;
solvecpu = t2 - t1 ;
fprintf(msgFile, "\n CPU %8.3f : solve done ", solvecpu) ;
/*
fprintf(msgFile, "\n\n cpus") ;
DVfprintf(msgFile, 6, cpus) ;
*/
DVzero(6, tcpus) ;
MPI_Reduce((void *) cpus, (void *) tcpus, 6, MPI_DOUBLE,
           MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, ", %.3f mflops", 1.e-6*solveops/(t2-t1)) ;
/*
   fprintf(msgFile, "\n\n tcpus") ;
   DVfprintf(msgFile, 6, tcpus) ;
*/
   solvecpu = DVsum(6, tcpus) ;
   fprintf(msgFile, "\n\n Solve CPU breakdown:") ;
   if ( solvecpu > 0.0 ) {
      fprintf(msgFile,
      "\n    set up solves               %8.3f %6.2f"
      "\n    load rhs and store solution %8.3f %6.2f"
      "\n    forward solve               %8.3f %6.2f"
      "\n    diagonal solve              %8.3f %6.2f"
      "\n    backward solve              %8.3f %6.2f"
      "\n    miscellaneous               %8.3f %6.2f"
      "\n    total time                  %8.3f",
      tcpus[0], 100.*tcpus[0]/solvecpu,
      tcpus[1], 100.*tcpus[1]/solvecpu,
      tcpus[2], 100.*tcpus[2]/solvecpu,
      tcpus[3], 100.*tcpus[3]/solvecpu,
      tcpus[4], 100.*tcpus[4]/solvecpu, 
      tcpus[5], 100.*tcpus[5]/solvecpu, solvecpu) ;
   }
}
fprintf(msgFile, 
        "\n\n Solve communication for this processor"
        "\n                   aggregates              solution"
        "\n                   #    # bytes          #     #bytes"
        "\n   send   %10d %10d %10d %10d"
        "\n   recv   %10d %10d %10d %10d",
        stats[1], stats[3], stats[0], stats[2],
        stats[5], stats[7], stats[4], stats[6]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 8, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
           "\n\n Solve communication for all processors"
           "\n                   aggregates              solution"
           "\n                   #    # bytes          #     #bytes"
           "\n   send   %10d %10d %10d %10d"
           "\n   recv   %10d %10d %10d %10d",
           tstats[1], tstats[3], tstats[0], tstats[2],
           tstats[5], tstats[7], tstats[4], tstats[6]) ;
   fflush(msgFile) ;
}
SubMtxManager_writeForHumanEye(solvemanager, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n computed solution") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(msgFile) ;
}
DenseMtx_sub(mtxZ, mtxX) ;
error  = DenseMtx_maxabs(mtxZ) ;
fprintf(msgFile, "\n\n local error  = %12.4e", error) ;
MPI_Reduce((void *) &error, (void *) &terror, 1, MPI_DOUBLE,
           MPI_MAX, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global error = %12.4e", terror) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n error matrix") ;
   DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxB) ;
DenseMtx_free(mtxZ) ;
InpMtx_free(mtxA) ;
SubMtxManager_free(factormanager) ;
SubMtxManager_free(solvemanager) ;
ChvManager_free(chvmanager) ;
FrontMtx_free(frontmtx) ;
ETree_free(frontETree) ;
IVL_free(symbfacIVL) ;
IV_free(vtxmapIV) ;
IV_free(frontOwnersIV) ;
SolveMap_free(solvemap) ;

MPI_Barrier(MPI_COMM_WORLD) ;

fclose(msgFile) ;

MPI_Finalize() ;

return(1) ; }

/*--------------------------------------------------------------------*/
