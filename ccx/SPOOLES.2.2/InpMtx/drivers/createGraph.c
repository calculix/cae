/*  createGraph.c  */

#include "../InpMtx.h"
#include "../../Graph.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------------
   read in a InpMtx object and create the Graph object

   created -- 97feb14, cca
   ----------------------------------------------------
*/
{
InpMtx   *inpmtx ;
FILE     *msgFile ;
Graph    *graph ;
int      count, msglvl, nvtx, rc ;
IVL      *adjIVL ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.inpmtxf or *.inpmtxb"
      "\n    outFile  -- output file, must be *.graphf or *.graphb"
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
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], argv[3], argv[4]) ;
fflush(msgFile) ;
/*
   --------------------------
   read in the InpMtx object
   --------------------------
*/
inpmtx = InpMtx_new() ;
if ( strcmp(argv[3], "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
rc = InpMtx_readFromFile(inpmtx, argv[3]) ;
fprintf(msgFile, "\n return value %d from InpMtx_readFromFile(%p,%s)",
        rc, inpmtx, argv[3]) ;
if ( rc != 1 ) {
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading InpMtx object from file %s",
           argv[3]) ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fflush(msgFile) ;
}
InpMtx_changeStorageMode(inpmtx, 3) ;
nvtx  = 1 + IV_max(&inpmtx->ivec1IV) ;
count = 1 + IV_max(&inpmtx->ivec2IV) ;
if ( nvtx < count ) {
   nvtx = count ;
}
/*
   ------------------------------------
   create the full adjacency IVL object
   ------------------------------------
*/
adjIVL = InpMtx_fullAdjacency(inpmtx) ;
/*
   ---------------------
   fill the Graph object
   ---------------------
*/
graph = Graph_new() ;
Graph_init2(graph, 0, nvtx, 0, adjIVL->tsize, nvtx, adjIVL->tsize, 
            adjIVL, NULL, NULL) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Graph object") ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------
   check that the graph is symmetric
   ---------------------------------
*/
if ( (rc = Graph_isSymmetric(graph)) == 1 ) {
   fprintf(msgFile, "\n\n graph is symmetric\n") ;
} else {
   fprintf(msgFile, "\n\n graph is not symmetric\n") ;
}
/*
   ---------------------------
   write out the Graph object
   ---------------------------
*/
if ( strcmp(argv[4], "none") != 0 ) {
   rc = Graph_writeToFile(graph, argv[4]) ;
   fprintf(msgFile, 
           "\n return value %d from Graph_writeToFile(%p,%s)",
           rc, graph, argv[4]) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;
InpMtx_free(inpmtx) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
