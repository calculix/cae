/*  testWrapperMT.c  */

#include "../BridgeMT.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] ) {
/*
   -----------------------------------------------------------
   purpose -- main driver program to solve a linear system
      where the matrix and rhs are read in from files and
      the solution is written to a file.
      NOTE: multithreaded version

   created -- 98sep24, cca
   -----------------------------------------------------------
*/
BridgeMT   *bridge ;
char       *mtxFileName, *rhsFileName, *solFileName ;
double     nfactorops ; 
FILE       *msgFile ;
InpMtx     *mtxA ;
int        error, msglvl, neqns, nfent, nfind, nfront, nrhs, nrow, 
           nsolveops, nthread, permuteflag, rc, seed, symmetryflag, 
           type ;
DenseMtx   *mtxX, *mtxY ;
/*--------------------------------------------------------------------*/
/*
    --------------------
    get input parameters
    --------------------
*/
if ( argc != 11 ) {
   fprintf(stdout, 
         "\n\n usage : %s msglvl msgFile neqns type symmetryflag "
         "\n         mtxFile rhsFile solFile seed nthread\n"
         "\n   msglvl  -- message level"
         "\n      0 -- no output"
         "\n      1 -- timings and statistics"
         "\n      2 and greater -- lots of output"
         "\n   msgFile -- message file"
         "\n   neqns   -- # of equations"
         "\n   type    -- type of entries"
         "\n      1 -- real"
         "\n      2 -- complex"
         "\n   symmetryflag -- symmetry flag"
         "\n      0 -- symmetric"
         "\n      1 -- hermitian"
         "\n      2 -- nonsymmetric"
         "\n   neqns   -- # of equations"
         "\n   mtxFile -- input file for A matrix InpMtx object"
         "\n      must be *.inpmtxf or *.inpmtxb"
         "\n   rhsFile -- input file for Y DenseMtx object"
         "\n      must be *.densemtxf or *.densemtxb"
         "\n   solFile -- output file for X DenseMtx object"
         "\n      must be none, *.densemtxf or *.densemtxb"
         "\n   seed -- random number seed"
         "\n   nthread -- number of threads"
         "\n",
         argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "w")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
neqns        = atoi(argv[3]) ;
type         = atoi(argv[4]) ;
symmetryflag = atoi(argv[5]) ;
mtxFileName  = argv[6] ;
rhsFileName  = argv[7] ;
solFileName  = argv[8] ;
seed         = atoi(argv[9]) ;
nthread      = atoi(argv[10]) ;
fprintf(msgFile, 
        "\n\n %s input :"
        "\n msglvl       = %d"
        "\n msgFile      = %s"
        "\n neqns        = %d"
        "\n type         = %d"
        "\n symmetryflag = %d"
        "\n mtxFile      = %s"
        "\n rhsFile      = %s"
        "\n solFile      = %s"
        "\n nthread      = %d"
        "\n",
        argv[0], msglvl, argv[2], neqns, type, symmetryflag, 
        mtxFileName, rhsFileName, solFileName, nthread) ;
/*--------------------------------------------------------------------*/
/*
   ------------------
   read in the matrix
   ------------------
*/
mtxA = InpMtx_new() ;
rc = InpMtx_readFromFile(mtxA, mtxFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n fatal error reading mtxA from file %s, rc = %d",
           mtxFileName, rc) ;
   fflush(msgFile) ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n InpMtx object ") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   read in the right hand side matrix
   ----------------------------------
*/
mtxY = DenseMtx_new() ;
rc = DenseMtx_readFromFile(mtxY, rhsFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n fatal error reading mtxY from file %s, rc = %d",
           rhsFileName, rc) ;
   fflush(msgFile) ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n DenseMtx object for right hand side") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
   fflush(msgFile) ;
}
DenseMtx_dimensions(mtxY, &nrow, &nrhs) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   create and setup a BridgeMT object
   ----------------------------------
*/
bridge = BridgeMT_new() ;
BridgeMT_setMatrixParams(bridge, neqns, type, symmetryflag) ;
BridgeMT_setMessageInfo(bridge, msglvl, msgFile) ;
rc = BridgeMT_setup(bridge, mtxA) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error return %d from BridgeMT_setup()", rc) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n ----- SETUP -----\n") ;
fprintf(msgFile, 
        "\n    CPU %8.3f : time to construct Graph"
        "\n    CPU %8.3f : time to compress Graph"
        "\n    CPU %8.3f : time to order Graph"
        "\n    CPU %8.3f : time for symbolic factorization"
        "\n CPU %8.3f : total setup time\n",
        bridge->cpus[0],
        bridge->cpus[1],
        bridge->cpus[2],
        bridge->cpus[3],
        bridge->cpus[4]) ;
rc = BridgeMT_factorStats(bridge, type, symmetryflag, &nfront,
                          &nfind, &nfent, &nsolveops, &nfactorops) ;
if ( rc != 1 ) {
   fprintf(stderr,
           "\n error return %d from BridgeMT_factorStats()", rc) ;
   exit(-1) ;
}
fprintf(msgFile,
        "\n\n factor matrix statistics"
        "\n %d fronts, %d indices, %d entries"
        "\n %d solve operations, %12.4e factor operations",
        nfront, nfind, nfent, nsolveops, nfactorops) ;
fflush(msgFile) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------
   setup the parallel factorization
   --------------------------------
*/
rc = BridgeMT_factorSetup(bridge, nthread, 0, 0.0) ;
fprintf(msgFile, "\n\n ----- PARALLEL FACTOR SETUP -----\n") ;
fprintf(msgFile, 
        "\n    CPU %8.3f : time to setup parallel factorization",
        bridge->cpus[5]) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n total factor operations = %.0f",
           DV_sum(bridge->cumopsDV)) ;
   fprintf(msgFile, 
           "\n upper bound on speedup due to load balance = %.2f",
           DV_sum(bridge->cumopsDV)/DV_max(bridge->cumopsDV)) ;
   fprintf(msgFile, "\n operations distributions over threads") ;
   DV_writeForHumanEye(bridge->cumopsDV, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -----------------
   factor the matrix
   -----------------
*/
permuteflag  = 1 ;
rc = BridgeMT_factor(bridge, mtxA, permuteflag, &error) ;
if ( rc == 1 ) {
   fprintf(msgFile, "\n\n factorization completed successfully\n") ;
} else {
   fprintf(msgFile, 
           "\n return code from factorization = %d\n"
           "\n error code                     = %d\n",
           rc, error) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n ----- FACTORIZATION -----\n") ;
fprintf(msgFile, 
        "\n    CPU %8.3f : time to permute original matrix"
        "\n    CPU %8.3f : time to initialize factor matrix" 
        "\n    CPU %8.3f : time to compute factorization"
        "\n    CPU %8.3f : time to post-process factorization"
        "\n CPU %8.3f : total factorization time\n",
        bridge->cpus[6],
        bridge->cpus[7],
        bridge->cpus[8],
        bridge->cpus[9],
        bridge->cpus[10]) ;
fprintf(msgFile, "\n\n factorization statistics"
        "\n %d pivots, %d pivot tests, %d delayed vertices"
        "\n %d entries in D, %d entries in L, %d entries in U",
        bridge->stats[0], bridge->stats[1], bridge->stats[2], 
        bridge->stats[3], bridge->stats[4], bridge->stats[5]) ;
fprintf(msgFile, 
        "\n\n factorization: raw mflops %8.3f, overall mflops %8.3f",
        1.e-6*nfactorops/bridge->cpus[8],
        1.e-6*nfactorops/bridge->cpus[10]) ;
fflush(msgFile) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------
   setup the parallel solve
   ------------------------
*/
rc = BridgeMT_solveSetup(bridge) ;
fprintf(msgFile, "\n\n ----- PARALLEL SOLVE SETUP -----\n") ;
fprintf(msgFile, 
        "\n    CPU %8.3f : time to setup parallel solve",
        bridge->cpus[11]) ;
/*--------------------------------------------------------------------*/
/*
   ----------------
   solve the system
   ----------------
*/
mtxX = DenseMtx_new() ;
DenseMtx_init(mtxX, type, 0, 0, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(mtxX) ;
rc = BridgeMT_solve(bridge, permuteflag, mtxX, mtxY) ;
if (rc == 1) {
   fprintf(msgFile, "\n\n solve complete successfully\n") ;
} else {
   fprintf(msgFile, "\n" " return code from solve = %d\n", rc) ;
}
fprintf(msgFile, "\n\n ----- SOLVE -----\n") ;
fprintf(msgFile, 
        "\n    CPU %8.3f : time to permute rhs into new ordering"
        "\n    CPU %8.3f : time to solve linear system"
        "\n    CPU %8.3f : time to permute solution into old ordering"
        "\n CPU %8.3f : total solve time\n",
        bridge->cpus[12], bridge->cpus[13],
        bridge->cpus[14], bridge->cpus[15]) ;
fprintf(msgFile, 
        "\n\n solve: raw mflops %8.3f, overall mflops %8.3f",
        1.e-6*nsolveops/bridge->cpus[13],
        1.e-6*nsolveops/bridge->cpus[15]) ;
fflush(msgFile) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n solution matrix in original ordering") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
if ( strcmp(solFileName, "none") != 0 ) {
/*
   -----------------------------------
   write the solution matrix to a file
   -----------------------------------
*/
   rc = DenseMtx_writeToFile(mtxX, solFileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile,
              "\n fatal error writing mtxX to file %s, rc = %d",
              solFileName, rc) ;
      fflush(msgFile) ;
      exit(-1) ;
   }
}
/*--------------------------------------------------------------------*/
/*
   ---------------------
   free the working data
   ---------------------
*/
InpMtx_free(mtxA) ;
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxY) ;
BridgeMT_free(bridge) ;

/*--------------------------------------------------------------------*/

return(1) ; }

/*--------------------------------------------------------------------*/
