/*  expandGraph.c  */

#include "../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------
   1) read in graph
   2) read in the expand map
   4) get expanded graph
   5) optionally write out expanded graph map

   created -- 95mar02, cca
   ------------------------------------------
*/
{
char     *inGraphFileName, *inIVfileName, *outGraphFileName ;
double   t1, t2 ;
int      msglvl, rc ;
Graph    *gc, *gf ;
FILE     *msgFile ;
IV       *mapIV ;

if ( argc != 6 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile "
      "\n         inMapFile outGraphFile"
      "\n    msglvl       -- message level"
      "\n    msgFile      -- message file"
      "\n    inGraphFile  -- input file for graph"
      "\n                    must be *.graphf or *.graphb"
      "\n    inMapFile    -- output file for map vector"
      "\n                    must be *.ivf or *.ivb"
      "\n    outGraphFile -- output file for compressed graph"
      "\n                    must be *.graphf or *.graphb"
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
inIVfileName     = argv[4] ;
outGraphFileName = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl       -- %d" 
        "\n msgFile      -- %s" 
        "\n inGraphFile  -- %s" 
        "\n inMapFile    -- %s" 
        "\n outGraphFile -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, 
        inIVfileName, outGraphFileName) ;
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
MARKTIME(t1) ;
gc = Graph_new() ;
rc = Graph_readFromFile(gc, inGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
        rc, gc, inGraphFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading Graph object from file %s",
           inGraphFileName) ;
   Graph_writeForHumanEye(gc, msgFile) ;
   fflush(msgFile) ;
} else {
   Graph_writeStats(gc, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------
   read in the IV object
   ---------------------
*/
if ( strcmp(inIVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
mapIV = IV_new() ;
rc = IV_readFromFile(mapIV, inIVfileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in IV from file %s",
        t2 - t1, inIVfileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
        rc, mapIV, inIVfileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading IV object from file %s",
           inIVfileName) ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
} else {
   IV_writeStats(mapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   create the expanded graph
   -------------------------
*/
MARKTIME(t1) ;
gf = Graph_expand2(gc, mapIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : create expanded graph", t2 - t1) ;
fprintf(msgFile, "\n\n expanded Graph object") ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(gf, msgFile) ;
   fflush(msgFile) ;
} else {
   Graph_writeStats(gf, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------
   write out the expanded Graph object
   -----------------------------------
*/
if ( strcmp(outGraphFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = Graph_writeToFile(gf, outGraphFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write expanded graph to file %s",
           t2 - t1, outGraphFileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from Graph_writeToFile(%p,%s)",
              rc, gf, outGraphFileName) ;
   }
}
/*
   ----------------
   free the objects
   ----------------
*/
Graph_free(gf) ;
Graph_free(gc) ;
IV_free(mapIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
