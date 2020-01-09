/*  testDomWeightStages.c  */

#include "../../Graph.h"
#include "../DSTree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   read in a DSTree object, read in a Graph file,
   read in a DV cutoffs file, get the stages IV object 
   based on domain weight and write it to a file.

   created -- 97jun12, cca
   ---------------------------------------------------
*/
{
char     *inCutoffDVfileName, *inDSTreeFileName, 
         *inGraphFileName, *outIVfileName ;
double   t1, t2 ;
DV       *cutoffDV ;
Graph    *graph ;
int      msglvl, rc ;
IV       *stagesIV ;
DSTree   *dstree ;
FILE     *msgFile ;

if ( argc != 7 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inDSTreeFile inGraphFile "
"\n         inCutoffDVfile outFile"
"\n    msglvl         -- message level"
"\n    msgFile        -- message file"
"\n    inDSTreeFile   -- input file, must be *.dstreef or *.dstreeb"
"\n    inGraphFile    -- input file, must be *.graphf or *.graphb"
"\n    inCutoffDVfile -- input file, must be *.dvf or *.dvb"
"\n    outFile        -- output file, must be *.ivf or *.ivb"
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
inDSTreeFileName   = argv[3] ;
inGraphFileName    = argv[4] ;
inCutoffDVfileName = argv[5] ;
outIVfileName      = argv[6] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl             -- %d" 
        "\n msgFile            -- %s" 
        "\n inDSTreeFileName   -- %s" 
        "\n inGraphFileName    -- %s" 
        "\n inCutoffDVfileName -- %s" 
        "\n outFile            -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inDSTreeFileName, 
        inGraphFileName, inCutoffDVfileName, outIVfileName) ;
fflush(msgFile) ;
/*
   -------------------------
   read in the DSTree object
   -------------------------
*/
if ( strcmp(inDSTreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
dstree = DSTree_new() ;
MARKTIME(t1) ;
rc = DSTree_readFromFile(dstree, inDSTreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in dstree from file %s",
        t2 - t1, inDSTreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from DSTree_readFromFile(%p,%s)",
           rc, dstree, inDSTreeFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading DSTree object from file %s",
        inDSTreeFileName) ;
if ( msglvl > 2 ) {
   DSTree_writeForHumanEye(dstree, msgFile) ;
} else {
   DSTree_writeStats(dstree, msgFile) ;
}
fflush(msgFile) ;
/*
   -------------------------
   read in the Graph object
   -------------------------
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
   -----------------------------
   read in the cutoffs DV object
   -----------------------------
*/
if ( strcmp(inCutoffDVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
cutoffDV = DV_new() ;
MARKTIME(t1) ;
rc = DV_readFromFile(cutoffDV, inCutoffDVfileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inCutoffDVfileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from DV_readFromFile(%p,%s)",
           rc, cutoffDV, inCutoffDVfileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading DV object from file %s",
        inCutoffDVfileName) ;
if ( msglvl > 0 ) {
   DV_writeForHumanEye(cutoffDV, msgFile) ;
} else {
   DV_writeStats(cutoffDV, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------
   get the stages vector
   ---------------------
*/
stagesIV = DSTree_stagesViaDomainWeight(dstree, 
                                        graph->vwghts, cutoffDV) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(stagesIV, msgFile) ;
} else {
   IV_writeStats(stagesIV, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------------
   write out the DSTree object
   ---------------------------
*/
if ( stagesIV != NULL && strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(stagesIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write dstree to file %s",
           t2 - t1, outIVfileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from IV_writeToFile(%p,%s)",
              rc, stagesIV, outIVfileName) ;
   }
}
/*
   ----------------------
   free the DSTree object
   ----------------------
*/
DSTree_free(dstree) ;
if ( stagesIV != NULL ) {
   IV_free(stagesIV) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
