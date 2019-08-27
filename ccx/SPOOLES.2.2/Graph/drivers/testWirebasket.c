/*  testWirebasket.c  */

#include "../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   read in a Graph and a stages id IV object,
   replace the stages IV object with wirebasket stages

   created -- 97jul30, cca
   ---------------------------------------------------
*/
{
char     *inCompidsFileName, *inGraphFileName, *outStagesIVfileName ;
double   t1, t2 ;
Graph    *graph ;
int      msglvl, nvtx, radius, rc, v ;
int      *compids, *stages ;
IV       *compidsIV, *stagesIV ;
FILE     *msgFile ;

if ( argc != 7 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile inStagesFile "
      "\n         outStagesFile radius"
      "\n    msglvl        -- message level"
      "\n    msgFile       -- message file"
      "\n    inGraphFile   -- input file, must be *.graphf or *.graphb"
      "\n    inStagesFile  -- output file, must be *.ivf or *.ivb"
      "\n    outStagesFile -- output file, must be *.ivf or *.ivb"
      "\n    radius        -- radius to set the stage "
      "\n                     of a separator vertex"
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
inGraphFileName     = argv[3] ;
inCompidsFileName  = argv[4] ;
outStagesIVfileName = argv[5] ;
radius              = atoi(argv[6]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inGraphFile   -- %s" 
        "\n inStagesFile  -- %s" 
        "\n outStagesFile -- %s" 
        "\n radius        -- %d" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, inCompidsFileName,
        outStagesIVfileName, radius) ;
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
   ---------------------
   read in the IV object
   ---------------------
*/
if ( strcmp(inCompidsFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
compidsIV = IV_new() ;
MARKTIME(t1) ;
rc = IV_readFromFile(compidsIV, inCompidsFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in compidsIV from file %s",
        t2 - t1, inCompidsFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
           rc, compidsIV, inCompidsFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IV object from file %s",
        inCompidsFileName) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(compidsIV, msgFile) ;
} else {
   IV_writeStats(compidsIV, msgFile) ;
}
fflush(msgFile) ;
IV_sizeAndEntries(compidsIV, &nvtx, &compids) ;
/*
   ----------------------------
   convert to the stages vector
   ----------------------------
*/
stagesIV = IV_new() ;
IV_init(stagesIV, nvtx, NULL) ;
stages = IV_entries(stagesIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] == 0 ) {
      stages[v] = 1 ;
   } else {
      stages[v] = 0 ;
   }
}
/*
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] == 0 ) {
      stages[v] = 0 ;
   } else {
      stages[v] = 1 ;
   }
}
*/
/*
   -------------------------
   get the wirebasket stages
   -------------------------
*/
Graph_wirebasketStages(graph, stagesIV, radius) ;
IV_sizeAndEntries(stagesIV, &nvtx, &stages) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( stages[v] == 2 ) {
      stages[v] = 1 ;
   } else if ( stages[v] > 2 ) {
      stages[v] = 2 ;
   }
}
fprintf(msgFile, "\n\n new stages IV object") ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(stagesIV, msgFile) ;
} else {
   IV_writeStats(stagesIV, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------------
   write out the stages object
   ---------------------------
*/
if ( strcmp(outStagesIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   IV_writeToFile(stagesIV, outStagesIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write stagesIV to file %s",
           t2 - t1, outStagesIVfileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from IV_writeToFile(%p,%s)",
              rc, stagesIV, outStagesIVfileName) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;
IV_free(stagesIV) ;
IV_free(compidsIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
