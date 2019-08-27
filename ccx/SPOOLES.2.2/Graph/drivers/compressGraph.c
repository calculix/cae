/*  compressGraph.c  */

#include "../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   1) read in graph
   2) get the equivalence map
   3) optionally write out equivalence map
   4) get compressed graph
   5) optionally write out compressed graph map

   created -- 95mar02, cca
   -------------------------------------------------
*/
{
char     *inGraphFileName, *outIVfileName, *outGraphFileName ;
double   t1, t2 ;
int      coarseType, msglvl, rc ;
Graph    *gc, *gf ;
FILE     *msgFile ;
IV       *mapIV ;

if ( argc != 7 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile "
      "\n         coarseType outMapFile outGraphFile"
      "\n    msglvl       -- message level"
      "\n    msgFile      -- message file"
      "\n    inGraphFile  -- input file for graph"
      "\n                    must be *.graphf or *.graphb"
      "\n    coarseType   -- type for compressed graph"
      "\n    outMapFile   -- output file for map vector"
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
coarseType       = atoi(argv[4]) ;
outIVfileName    = argv[5] ;
outGraphFileName = argv[6] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl       -- %d" 
        "\n msgFile      -- %s" 
        "\n inGraphFile  -- %s" 
        "\n coarseType   -- %d" 
        "\n outMapFile   -- %s" 
        "\n outGraphFile -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, 
        coarseType, outIVfileName, outGraphFileName) ;
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
gf = Graph_new() ;
rc = Graph_readFromFile(gf, inGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
        rc, gf, inGraphFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading Graph object from file %s",
           inGraphFileName) ;
   Graph_writeForHumanEye(gf, msgFile) ;
   fflush(msgFile) ;
} else {
   Graph_writeStats(gf, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------
   get the equivalence map
   -----------------------
*/
MARKTIME(t1) ;
mapIV = Graph_equivMap(gf) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : create equivalence map", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n equivalence map IV object") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the map IV object
   ---------------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(mapIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write map IV to file %s",
           t2 - t1, outIVfileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, "\n return value %d from IV_writeToFile(%p,%s)",
              rc, mapIV, outIVfileName) ;
   }
}
/*
   ------------------------
   get the compressed graph
   ------------------------
*/
MARKTIME(t1) ;
gc = Graph_compress(gf, IV_entries(mapIV), coarseType) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : compress the graph", t2 - t1) ;
fprintf(msgFile, "\n\n compressed graph") ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(gc, msgFile) ;
   fflush(msgFile) ;
} else {
   Graph_writeStats(gc, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------
   write out the Graph object
   --------------------------
*/
if ( strcmp(outGraphFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = Graph_writeToFile(gc, outGraphFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write compressed graph to file %s",
           t2 - t1, outGraphFileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from Graph_writeToFile(%p,%s)",
              rc, gc, outGraphFileName) ;
      exit(-1) ;
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
