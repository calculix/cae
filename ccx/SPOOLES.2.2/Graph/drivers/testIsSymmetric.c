/*  testIsSymmetric.c  */

#include "../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------------------
   check to see whether a graph is symmetric

   created -- 96oct31, cca
   -----------------------------------------
*/
{
char     *inGraphFileName ;
double   t1, t2 ;
int      msglvl, rc ;
Graph    *graph ;
FILE     *msgFile ;

if ( argc != 4 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile "
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.graphf or *.graphb"
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
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName) ;
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
   -------------------------------------------
   check to see whether the graph is symmetric
   -------------------------------------------
*/
rc = Graph_isSymmetric(graph) ;
if ( rc == 1 ) {
   fprintf(msgFile, "\n graph is symmetric") ;
} else {
   fprintf(msgFile, "\n graph is not symmetric") ;
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
