/*  testIO.c  */

#include "../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   test Graph_readFromFile and Graph_writeToFile,
   useful for translating between formatted *.graphf
   and binary *.graphb files.

   created -- 96fmar02, cca
   -------------------------------------------------
*/
{
char     *inGraphFileName, *outGraphFileName ;
double   t1, t2 ;
int      msglvl, rc ;
Graph    *graph ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.graphf or *.graphb"
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
inGraphFileName  = argv[3] ;
outGraphFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, outGraphFileName) ;
fflush(msgFile) ;
/*
   ----------------------
   read in the Graph object
   ----------------------
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
   write out the Graph object
   ------------------------
*/
if ( strcmp(outGraphFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = Graph_writeToFile(graph, outGraphFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write graph to file %s",
           t2 - t1, outGraphFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_writeToFile(%p,%s)",
           rc, graph, outGraphFileName) ;
}
/*
   ---------------------
   free the Graph object
   ---------------------
*/
Graph_free(graph) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
