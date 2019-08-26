/*  testTransform.c  */

#include "../../ETree.h"
#include "../../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------
   read in an ETree object and testTransform it by 

   1 -- merge only children if possible
   2 -- merge all children if possible
   3 -- split large non-leaf fronts

   created -- 96jun27, cca
   -------------------------------------------------------
*/
{
char     *inETreeFileName, *inGraphFileName, *outETreeFileName ;
double   cpus[6], ops[6], t1, t2 ;
ETree    *etree0, *etree1, *etree2, *etree3, *etree4, *etree5 ;
FILE     *msgFile ;
Graph    *graph ;
int      maxsize, maxzeros, msglvl, rc, seed ;
int      nfronts[7], nfind[7], nzf[7] ;
IV       *nzerosIV ;
IVL      *symbfacIVL ;

if ( argc != 9 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile inGraphFile outETreeFile "
"\n         maxzeros maxsize seed"
      "\n    msglvl       -- message level"
      "\n    msgFile      -- message file"
      "\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
      "\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
      "\n    outETreeFile -- output file, must be *.etreef or *.etreeb"
      "\n    maxzeros     -- maximum number of zeros in a front"
      "\n    maxsize      -- maximum number of vertices in a front"
      "\n    seed         -- random number seed"
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
inETreeFileName  = argv[3] ;
inGraphFileName  = argv[4] ;
outETreeFileName = argv[5] ;
maxzeros         = atoi(argv[6]) ;
maxsize          = atoi(argv[7]) ;
seed             = atoi(argv[8]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl       -- %d" 
        "\n msgFile      -- %s" 
        "\n inETreeFile  -- %s" 
        "\n inGraphFile  -- %s" 
        "\n outETreeFile -- %s" 
        "\n maxzeros     -- %d" 
        "\n maxsize      -- %d" 
        "\n seed         -- %d" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName, inGraphFileName, 
        outETreeFileName, maxzeros, maxsize, seed) ;
fflush(msgFile) ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
if ( strcmp(inGraphFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
graph = Graph_new() ;
MARKTIME(t1) ;
rc = Graph_readFromFile(graph, inGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, graph, inGraphFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading Graph object from file %s",
        inGraphFileName) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
} else {
   Graph_writeStats(graph, msgFile) ;
}
fflush(msgFile) ;
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
if ( strcmp(inETreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
etree0 = ETree_new() ;
MARKTIME(t1) ;
rc = ETree_readFromFile(etree0, inETreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in etree from file %s",
        t2 - t1, inETreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
           rc, etree0, inETreeFileName) ;
   exit(-1) ;
}
nfronts[0] = ETree_nfront(etree0) ;
nfind[0]   = ETree_nFactorIndices(etree0) ;
nzf[0]     = ETree_nFactorEntries(etree0, SPOOLES_SYMMETRIC) ;
ops[0]     = ETree_nFactorOps(etree0, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n\n original  : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[0], nfind[0], nzf[0], ops[0]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n original front tree ") ;
   ETree_writeForHumanEye(etree0, msgFile) ;
   fflush(msgFile) ;
}
fflush(msgFile) ;
/*
   ----------------------------------
   get the fundamental supernode tree
   ----------------------------------
*/
nzerosIV = IV_new() ;
IV_init(nzerosIV, nfronts[0], NULL) ;
IV_fill(nzerosIV, 0) ;
MARKTIME(t1) ;
etree1 = ETree_mergeFrontsOne(etree0, 0, nzerosIV) ;
MARKTIME(t2) ;
cpus[1] = t2 - t1 ;
fprintf(msgFile, "\n CPU %9.5f : get fundamental supernode tree", 
        t2 - t1) ;
nfronts[1] = ETree_nfront(etree1) ;
nfind[1]   = ETree_nFactorIndices(etree1) ;
nzf[1]     = ETree_nFactorEntries(etree1, SPOOLES_SYMMETRIC) ;
ops[1]     = ETree_nFactorOps(etree1, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n merge one : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[1], nfind[1], nzf[1], ops[1]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after first merge") ;
   ETree_writeForHumanEye(etree1, msgFile) ;
}
fprintf(msgFile, "\n IV_sum(nzerosIV) = %d, IV_max(nzerosIV) = %d", 
        IV_sum(nzerosIV), IV_max(nzerosIV)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
ETree_writeToFile(etree1, 
                  "/local1/cleve/ARPA/matrices/R3D13824/nd1.etreef") ;
/*
   ---------------------------
   try to absorb only children
   ---------------------------
*/
MARKTIME(t1) ;
etree2 = ETree_mergeFrontsOne(etree1, maxzeros, nzerosIV) ;
MARKTIME(t2) ;
cpus[2] = t2 - t1 ;
fprintf(msgFile, "\n CPU %9.5f : merge only single child", t2 - t1) ;
nfronts[2] = ETree_nfront(etree2) ;
nfind[2]   = ETree_nFactorIndices(etree2) ;
nzf[2]     = ETree_nFactorEntries(etree2, SPOOLES_SYMMETRIC) ;
ops[2]     = ETree_nFactorOps(etree2, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n merge one : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[2], nfind[2], nzf[2], ops[2]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after first merge") ;
   ETree_writeForHumanEye(etree2, msgFile) ;
}
fprintf(msgFile, "\n IV_sum(nzerosIV) = %d, IV_max(nzerosIV) = %d", 
        IV_sum(nzerosIV), IV_max(nzerosIV)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
ETree_writeToFile(etree2, 
                  "/local1/cleve/ARPA/matrices/R3D13824/nd2.etreef") ;
/*
   --------------------------
   try to absorb all children
   --------------------------
*/
MARKTIME(t1) ;
etree3 = ETree_mergeFrontsAll(etree2, maxzeros, nzerosIV) ;
MARKTIME(t2) ;
cpus[3] = t2 - t1 ;
fprintf(msgFile, "\n CPU %9.5f : merge all children", t2 - t1) ;
nfronts[3] = ETree_nfront(etree3) ;
nfind[3]   = ETree_nFactorIndices(etree3) ;
nzf[3]     = ETree_nFactorEntries(etree3, SPOOLES_SYMMETRIC) ;
ops[3]     = ETree_nFactorOps(etree3, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n merge all : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
                 nfronts[3], nfind[3], nzf[3], ops[3]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after second merge") ;
   ETree_writeForHumanEye(etree2, msgFile) ;
}
fprintf(msgFile, "\n IV_sum(nzerosIV) = %d, IV_max(nzerosIV) = %d", 
        IV_sum(nzerosIV), IV_max(nzerosIV)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
ETree_writeToFile(etree3, 
                  "/local1/cleve/ARPA/matrices/R3D13824/nd3.etreef") ;
/*
   --------------------------
   try to absorb any children
   --------------------------
*/
MARKTIME(t1) ;
etree4 = ETree_mergeFrontsAny(etree3, maxzeros, nzerosIV) ;
MARKTIME(t2) ;
cpus[4] = t2 - t1 ;
fprintf(msgFile, "\n CPU %9.5f : merge any child", t2 - t1) ;
nfronts[4] = ETree_nfront(etree4) ;
nfind[4]   = ETree_nFactorIndices(etree4) ;
nzf[4]     = ETree_nFactorEntries(etree4, SPOOLES_SYMMETRIC) ;
ops[4]     = ETree_nFactorOps(etree4, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n merge any : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
                 nfronts[4], nfind[4], nzf[4], ops[4]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after second merge") ;
   ETree_writeForHumanEye(etree4, msgFile) ;
}
fprintf(msgFile, "\n IV_sum(nzerosIV) = %d, IV_max(nzerosIV) = %d", 
        IV_sum(nzerosIV), IV_max(nzerosIV)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
ETree_writeToFile(etree4, 
                  "/local1/cleve/ARPA/matrices/R3D13824/nd4.etreef") ;
/*
   --------------------
   split the front tree
   --------------------
*/
MARKTIME(t1) ;
etree5 = ETree_splitFronts(etree4, graph->vwghts, maxsize, 0) ;
MARKTIME(t2) ;
cpus[5] = t2 - t1 ;
fprintf(msgFile, "\n CPU %9.5f : split fronts", t2 - t1) ;
nfronts[5] = ETree_nfront(etree5) ;
nfind[5]   = ETree_nFactorIndices(etree5) ;
nzf[5]     = ETree_nFactorEntries(etree5, SPOOLES_SYMMETRIC) ;
ops[5]     = ETree_nFactorOps(etree5, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n split     : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[5], nfind[5], nzf[5], ops[5]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after split") ;
   ETree_writeForHumanEye(etree5, msgFile) ;
}
/*
   ----------------------------------------
   create the symbolic factorization object
   ----------------------------------------
*/
MARKTIME(t1) ;
symbfacIVL = SymbFac_initFromGraph(etree5, graph) ;
MARKTIME(t2) ;
cpus[6] = t2 - t1 ;
fprintf(msgFile, "\n CPU %9.5f : symbolic factorization", t2 - t1) ;
nfronts[6] = ETree_nfront(etree5) ;
nfind[6]   = ETree_nFactorIndices(etree5) ;
nzf[6]     = ETree_nFactorEntries(etree5, SPOOLES_SYMMETRIC) ;
ops[6]     = ETree_nFactorOps(etree5, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n final     : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[6], nfind[6], nzf[6], ops[6]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after symbolic factorization") ;
   ETree_writeForHumanEye(etree5, msgFile) ;
   fprintf(msgFile, "\n\n after symbolic factorization") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
}
fprintf(msgFile, "\n\n"
"\n                  CPU  #fronts  #indices  #entries       #ops"
"\n original  :          %8d  %8d  %8d %12.0f "
"\n fs tree   : %8.3f %8d  %8d  %8d %12.0f "
"\n merge one : %8.3f %8d  %8d  %8d %12.0f "
"\n merge all : %8.3f %8d  %8d  %8d %12.0f "
"\n merge any : %8.3f %8d  %8d  %8d %12.0f "
"\n split     : %8.3f %8d  %8d  %8d %12.0f "
"\n final     : %8.3f %8d  %8d  %8d %12.0f ",
         nfronts[0], nfind[0], nzf[0], ops[0],
cpus[1], nfronts[1], nfind[1], nzf[1], ops[1],
cpus[2], nfronts[2], nfind[2], nzf[2], ops[2],
cpus[3], nfronts[3], nfind[3], nzf[3], ops[3],
cpus[4], nfronts[4], nfind[4], nzf[4], ops[4],
cpus[5], nfronts[5], nfind[5], nzf[5], ops[5],
cpus[6], nfronts[6], nfind[6], nzf[6], ops[6]) ;
ETree_writeToFile(etree5, 
                  "/local1/cleve/ARPA/matrices/R3D13824/nd5.etreef") ;
/*
   --------------------------
   write out the ETree object
   --------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ETree_writeToFile(etree5, outETreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outETreeFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_writeToFile(%p,%s)",
           rc, etree5, outETreeFileName) ;
}
/*
   ----------------------
   free the ETree objects
   ----------------------
*/
/*
ETree_free(etree0) ;
ETree_free(etree1) ;
ETree_free(etree2) ;
ETree_free(etree3) ;
ETree_free(etree4) ;
ETree_free(etree5) ;
*/
Graph_free(graph) ;
IVL_free(symbfacIVL) ;
IV_free(nzerosIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
