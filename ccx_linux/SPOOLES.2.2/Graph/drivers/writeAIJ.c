/*  writeAIJ.c  */

#include "../Graph.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   read in a Graph and write out an AIJ file
   
   nrow  ncol  nent
   row col entry
   ...

   created -- 98dec05, cca
   -------------------------------------------------
*/
{
char     *inGraphFileName, *outFileName ;
double   t1, t2 ;
Drand    *drand ;
int      ii, msglvl, nedges, nvtx, rc, v, vsize, w ;
int      *vadj ;
Graph    *graph ;
FILE     *msgFile, *outFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFileName"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.graphf or *.graphb"
      "\n    outFile  -- output file"
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
outFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, outFileName) ;
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
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new();
Drand_setSeed(drand, 47) ;
Drand_setUniform(drand, -1., 1.) ;
/*
   --------------------------
   write out the AIJ entries
   --------------------------
*/
if ( strcmp(outFileName, "stdout") == 0 ) {
   outFile = stdout ;
} else if ( (outFile = fopen(outFileName, "w")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n", argv[0], outFile) ;
   return(-1) ;
}
nvtx   = graph->nvtx ;
nedges = graph->nedges ;
fprintf(outFile, "\n %d %d %d", nvtx, nvtx, nvtx + (nedges-nvtx)/2) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   Graph_adjAndSize(graph, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      w = vadj[ii] ;
      if ( v <= w ) {
         fprintf(outFile,  
                 "\n %8d %8d %20.12e", v, w, Drand_value(drand)) ;
      }
   }
}
fclose(outFile) ;
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
