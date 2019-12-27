/*  testOrderViaMMD.c  */

#include "../misc.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------------
   order a graph via a multiple minimum degree ordering

   created -- 97nov08, cca
   ----------------------------------------------------
*/
{
char     *inGraphFileName, *outETreeFileName ;
double   t1, t2 ;
ETree    *frontETree ;
int      maxsize, maxzeros, msglvl, nvtx, rc, seed ;
Graph    *graph ;
FILE     *msgFile ;

if ( argc != 8 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile GraphFile "
"\n         maxsize maxzeros seed ETreeFile"
"\n    msglvl        -- message level"
"\n    msgFile       -- message file"
"\n    GraphFile     -- input graph file, must be *.graphf or *.graphb"
"\n    maxzeros      -- maximum # of zeros in a front"
"\n    maxsize       -- maximum # of internal columns in a front"
"\n    seed          -- random number seed"
"\n    ETreeFile     -- output ETree file, must be *.etreef or *.etreeb"
"\n",
argv[0]) ;
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
inGraphFileName  = argv[3] ;
maxsize          = atoi(argv[4]) ;
maxzeros         = atoi(argv[5]) ;
seed             = atoi(argv[6]) ;
outETreeFileName = argv[7] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n GraphFile     -- %s" 
        "\n maxsize       -- %d" 
        "\n maxzeros      -- %d" 
        "\n seed          -- %d" 
        "\n ETreeFile     -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, maxsize, 
        maxzeros, seed, outETreeFileName) ;
fflush(msgFile) ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
graph = Graph_new() ;
if ( strcmp(inGraphFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
rc = Graph_readFromFile(graph, inGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
nvtx = graph->nvtx ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, graph, inGraphFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading Graph object from file %s",
           inGraphFileName) ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   order the graph using minimum degree
   ------------------------------------
*/
MARKTIME(t1) ;
frontETree = orderViaMMD(graph, seed, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %9.2f : order the graph via MMD", t2 - t1) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(frontETree, msgFile) ;
   fflush(msgFile) ;
}
fprintf(msgFile,
        "\n\n %d fronts, %d indices, %d entries, %.0f operations",
        frontETree->nfront,
        ETree_nFactorIndices(frontETree),
        ETree_nFactorEntries(frontETree, SPOOLES_SYMMETRIC),
        ETree_nFactorOps(frontETree, SPOOLES_REAL, SPOOLES_SYMMETRIC)) ;
/*
   ------------------
   transform the tree
   ------------------
*/
MARKTIME(t1) ;
frontETree = ETree_transform(frontETree, graph->vwghts, 
                             maxzeros, maxsize, seed);
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : transform tree", t2 - t1) ;
   fflush(msgFile) ;
}
fprintf(msgFile,
        "\n\n %d fronts, %d indices, %d entries, %.0f operations",
        frontETree->nfront,
        ETree_nFactorIndices(frontETree),
        ETree_nFactorEntries(frontETree, SPOOLES_SYMMETRIC),
        ETree_nFactorOps(frontETree, SPOOLES_REAL, SPOOLES_SYMMETRIC)) ;
/*
   ------------------------------------
   optionally write the ETree to a file
   ------------------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   ETree_writeToFile(frontETree, outETreeFileName) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;
ETree_free(frontETree) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
