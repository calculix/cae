/*  mkNDETreeNew.c  */

#include "../../ETree.h"
#include "../../EGraph.h"
#include "../../misc.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------------------
   make ETree objects for nested dissection on a regular grid

   1 -- vertex elimination tree
   2 -- fundamental supernode front tree
   3 -- merge only children if possible
   4 -- merge all children if possible
   5 -- split large non-leaf fronts

   created -- 98feb05, cca
   ------------------------------------------------------------
*/
{
char     *outETreeFileName ;
double   ops[6] ;
double   t1, t2 ;
EGraph   *egraph ;
ETree    *etree0, *etree1, *etree2, *etree3, *etree4, *etree5 ;
FILE     *msgFile ;
Graph    *graph ;
int      nfronts[6], nfind[6], nzf[6] ; 
int      maxsize, maxzeros, msglvl, n1, n2, n3, nvtx, rc, v ;
int      *newToOld, *oldToNew ;
IV       *nzerosIV ;

if ( argc != 9 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile n1 n2 n3 maxzeros maxsize outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    n1       -- number of points in the first direction"
      "\n    n2       -- number of points in the second direction"
      "\n    n3       -- number of points in the third direction"
      "\n    maxzeros -- number of points in the third direction"
      "\n    maxsize  -- maximum number of vertices in a front"
      "\n    outFile  -- output file, must be *.etreef or *.etreeb"
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
n1 = atoi(argv[3]) ;
n2 = atoi(argv[4]) ;
n3 = atoi(argv[5]) ;
maxzeros = atoi(argv[6]) ;
maxsize  = atoi(argv[7]) ;
outETreeFileName = argv[8] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n n1       -- %d" 
        "\n n2       -- %d" 
        "\n n3       -- %d" 
        "\n maxzeros -- %d" 
        "\n maxsize  -- %d" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, 
        maxzeros, maxsize, outETreeFileName) ;
fflush(msgFile) ;
/*
   ----------------------------
   create the grid graph object
   ----------------------------
*/
if ( n1 == 1 ) {
   egraph = EGraph_make9P(n2, n3, 1) ;
} else if ( n2 == 1 ) {
   egraph = EGraph_make9P(n1, n3, 1) ;
} else if ( n3 == 1 ) {
   egraph = EGraph_make9P(n1, n2, 1) ;
} else {
   egraph = EGraph_make27P(n1, n2, n3, 1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %d x %d x %d grid EGraph", n1, n2, n3) ;
   EGraph_writeForHumanEye(egraph, msgFile) ;
   fflush(msgFile) ;
}
graph = EGraph_mkAdjGraph(egraph) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %d x %d x %d grid Graph", n1, n2, n3) ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   get the nested dissection ordering
   ----------------------------------
*/
nvtx = n1*n2*n3 ;
newToOld = IVinit(nvtx, -1) ;
oldToNew = IVinit(nvtx, -1) ;
mkNDperm(n1, n2, n3, newToOld, 0, n1-1, 0, n2-1, 0, n3-1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   oldToNew[newToOld[v]] = v ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %d x %d x %d nd ordering", n1, n2, n3) ;
   IVfprintf(msgFile, nvtx, oldToNew) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------
   create the vertex elimination ETree object
   ------------------------------------------
*/
etree0 = ETree_new() ;
ETree_initFromGraphWithPerms(etree0, graph, newToOld, oldToNew) ;
nfronts[0] = ETree_nfront(etree0) ;
nfind[0]   = ETree_nFactorIndices(etree0) ;
nzf[0]     = ETree_nFactorEntries(etree0, SPOOLES_SYMMETRIC) ;
ops[0]     = ETree_nFactorOps(etree0, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile, 
        "\n vtx tree  : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[0], nfind[0], nzf[0], ops[0]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n vertex elimination tree") ;
   ETree_writeForHumanEye(etree0, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   create the fundamental supernode ETree object
   ---------------------------------------------
*/
nzerosIV = IV_new() ;
IV_init(nzerosIV, nvtx, NULL) ;
IV_fill(nzerosIV, 0) ;
etree1     = ETree_mergeFrontsOne(etree0, 0, nzerosIV) ;
nfronts[1] = ETree_nfront(etree1) ;
nfind[1]   = ETree_nFactorIndices(etree1) ;
nzf[1]     = ETree_nFactorEntries(etree1, SPOOLES_SYMMETRIC) ;
ops[1]     = ETree_nFactorOps(etree1, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile, 
        "\n fs tree   : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[1], nfind[1], nzf[1], ops[1]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n fundamental supernode front tree") ;
   ETree_writeForHumanEye(etree1, msgFile) ;
   fprintf(msgFile, "\n\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   try to absorb only children
   ---------------------------
*/
etree2 = ETree_mergeFrontsOne(etree1, maxzeros, nzerosIV) ;
nfronts[2] = ETree_nfront(etree2) ;
nfind[2]   = ETree_nFactorIndices(etree2) ;
nzf[2]     = ETree_nFactorEntries(etree2, SPOOLES_SYMMETRIC) ;
ops[2]     = ETree_nFactorOps(etree2, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile, 
        "\n merge one : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[2], nfind[2], nzf[2], ops[2]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after mergeOne") ;
   ETree_writeForHumanEye(etree2, msgFile) ;
   fprintf(msgFile, "\n\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------
   try to absorb all children
   --------------------------
*/
etree3 = ETree_mergeFrontsAll(etree2, maxzeros, nzerosIV) ;
nfronts[3] = ETree_nfront(etree3) ;
nfind[3]   = ETree_nFactorIndices(etree3) ;
nzf[3]     = ETree_nFactorEntries(etree3, SPOOLES_SYMMETRIC) ;
ops[3]     = ETree_nFactorOps(etree3, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile, 
        "\n merge all : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
                 nfronts[3], nfind[3], nzf[3], ops[3]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after mergeAll") ;
   ETree_writeForHumanEye(etree3, msgFile) ;
   fprintf(msgFile, "\n\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------
   try to absorb any other children
   --------------------------------
*/
etree4 = etree3 ;
/*
etree4 = ETree_mergeFrontsAny(etree3, maxzeros, nzerosIV) ;
nfronts[4] = ETree_nfront(etree4) ;
nfind[4]   = ETree_nFactorIndices(etree4) ;
nzf[4]     = ETree_nFactorEntries(etree4, SPOOLES_SYMMETRIC) ;
ops[4]     = ETree_nFactorOps(etree4, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile, 
        "\n merge any : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
                 nfronts[4], nfind[4], nzf[4], ops[4]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after mergeAny") ;
   ETree_writeForHumanEye(etree3, msgFile) ;
   fprintf(msgFile, "\n\n nzerosIV") ;
   IV_writeForHumanEye(nzerosIV, msgFile) ;
   fflush(msgFile) ;
}
*/
/*
   --------------------
   split the front tree
   --------------------
*/
etree5 = ETree_splitFronts(etree4, NULL, maxsize, 0) ;
nfronts[5] = ETree_nfront(etree5) ;
nfind[5]   = ETree_nFactorIndices(etree5) ;
nzf[5]     = ETree_nFactorEntries(etree5, SPOOLES_SYMMETRIC) ;
ops[5]     = ETree_nFactorOps(etree5, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile, 
        "\n split     : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[5], nfind[5], nzf[5], ops[5]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after split") ;
   ETree_writeForHumanEye(etree4, msgFile) ;
   fflush(msgFile) ;
}
fprintf(msgFile, "\n\n complex symmetric ops %.0f",
        ETree_nFactorOps(etree5, SPOOLES_COMPLEX, SPOOLES_SYMMETRIC)) ;
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
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from ETree_writeToFile(%p,%s)",
              rc, etree5, outETreeFileName) ;
   }
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree0) ;
ETree_free(etree1) ;
ETree_free(etree2) ;
ETree_free(etree3) ;
/*
ETree_free(etree4) ;
*/
ETree_free(etree5) ;
EGraph_free(egraph) ;
Graph_free(graph) ;
IVfree(newToOld) ;
IVfree(oldToNew) ;
IV_free(nzerosIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
