/*  testSmooth4.c  */

#include "../GPart.h"
#include "../../DSTree.h"
#include "../../MSMD.h"
#include "../../BKL.h"
#include "../../ETree.h"
#include "../../Perm.h"
#include "../../IV.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------
   (1) read in a Graph, 
   (2) read in a 2-set partition
   (3) smooth the separator
   (4) optionally write out new partition

   created -- 96oct21, cca
   ------------------------------------------
*/
{
char        *inGraphFileName, *inIVfileName, 
            *msgFileName, *outIVfileName ;
double      alpha, cost, t1, t2 ;
FILE        *msgFile ;
GPart       *gpart ;
Graph       *graph ;
int         ierr, msglvl, option, nvtx, rc ;
IV          *tagsIV ;

if ( argc != 8 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile inIVfile "
      "\n         option alpha outIVfile"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inGraphFile -- input file, must be *.graphf or *.graphb"
      "\n    inIVFile    -- input file, must be *.ivf or *.ivb"
      "\n    option      -- smoothing option"
      "\n       1 --> bipartite 2-layers"
      "\n       2 --> non-bipartite 2-layers"
      "\n       k > 2 --> k/2 layers on each side"
      "\n    alpha       -- cost function parameter"
      "\n    outIVFile   -- output file, must be *.ivf or *.ivb"
      "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
msgFileName = argv[2] ;
if ( strcmp(msgFileName, "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(msgFileName, "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], msgFileName) ;
   return(-1) ;
}
inGraphFileName = argv[3] ;
inIVfileName    = argv[4] ;
option          = atoi(argv[5]) ;
alpha           = atof(argv[6]) ;
outIVfileName   = argv[7] ;
fprintf(msgFile, 
        "\n %s : "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inGraphFile -- %s" 
        "\n inIVfile    -- %s" 
        "\n option      -- %d" 
        "\n alpha       -- %f" 
        "\n outIVfile   -- %s" 
        "\n", argv[0], msglvl, 
        msgFileName, inGraphFileName, inIVfileName, 
        option, alpha, outIVfileName) ;
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
graph = Graph_new() ;
Graph_setDefaultFields(graph) ;
if ( (rc = Graph_readFromFile(graph, inGraphFileName)) != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
        rc, graph, inGraphFileName) ;
   exit(-1) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s", 
        t2 - t1, inGraphFileName) ;
nvtx = graph->nvtx ;
if ( msglvl < 4 ) {
   Graph_writeStats(graph, msgFile) ;
   fflush(msgFile) ;
} else {
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   read in the IV object
   ------------------------
*/
if ( strcmp(inIVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
tagsIV = IV_new() ;
IV_setDefaultFields(tagsIV) ;
if ( (rc = IV_readFromFile(tagsIV, inIVfileName)) != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
        rc, tagsIV, inIVfileName) ;
   exit(-1) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in tagsIV from file %s", 
        t2 - t1, inIVfileName) ;
nvtx = IV_size(tagsIV) ;
if ( msglvl < 4 ) {
   IV_writeStats(tagsIV, msgFile) ;
   fflush(msgFile) ;
} else {
   IV_writeForHumanEye(tagsIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------
   create the GPart object
   -----------------------
*/
MARKTIME(t1) ;
gpart = GPart_new() ;
GPart_init(gpart, graph) ;
GPart_setMessageInfo(gpart, msglvl, msgFile) ;
IV_copy(&gpart->compidsIV, tagsIV) ;
GPart_setCweights(gpart) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : GPart initialized", t2 - t1) ;
fprintf(msgFile, "\n cweights : ") ;
IVfp80(msgFile, 1 + gpart->ncomp, IV_entries(&gpart->cweightsIV),
       30, &ierr) ;
fflush(msgFile) ;
if ( msglvl > 0 ) {
   int   wB, wS, wW ;
   int   *cweights = IV_entries(&gpart->cweightsIV) ;

   fprintf(msgFile, "\n initial component weights :") ;
   IVfp80(msgFile, 1+gpart->ncomp, cweights, 20, &ierr) ;
   wS = cweights[0] ;
   wB = cweights[1] ;
   wW = cweights[2] ;
   if ( wB < wW ) {
      wB = wW ;
      wW = cweights[1] ;
   }
   cost = wS*(1. + (alpha*wB)/wW) ;
   fprintf(msgFile, 
   "\n initial |S| = %d , balance = %6.3f , cpu = %8.3f , cost = %8.1f\n",
           wS, ((double) wB)/wW, t2 - t1, cost) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   check that we have a valid vertex separator
   -------------------------------------------
*/
if ( 1 != GPart_validVtxSep(gpart) ) {
   fprintf(msgFile, "\n\n WHOA : separator is not valid"
           "\n\n WHOA : separator is not valid\n\n") ;
   fflush(msgFile) ;
}
/*
   --------------------
   smooth the separator
   --------------------
*/
MARKTIME(t1) ;
if ( option <= 2 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, 
              "\n calling GPart_smoothBy2layers with option %d",
              option) ;
      fflush(stdout) ;
   }
   GPart_smoothBy2layers(gpart, option, alpha) ;
} else {
   if ( msglvl > 1 ) {
      fprintf(msgFile, 
              "\n calling GPart_smoothBisector with option %d",
              option/2) ;
      fflush(stdout) ;
   }
   GPart_smoothBisector(gpart, option/2, alpha) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %9.5f : improve bisector", t2 - t1) ;
if ( msglvl > 0 ) {
   int   wB, wS, wW ;
   int   *cweights = IV_entries(&gpart->cweightsIV) ;

   fprintf(msgFile, "\n final component weights :") ;
   IVfp80(msgFile, 1+gpart->ncomp, cweights, 20, &ierr) ;
   wS = cweights[0] ;
   wB = cweights[1] ;
   wW = cweights[2] ;
   if ( wB < wW ) {
      wB = wW ;
      wW = cweights[1] ;
   }
   cost = wS*(1. + (alpha*wB)/wW) ;
   fprintf(msgFile, 
   "\n final |S| = %d , balance = %6.3f , cpu = %8.3f , cost = %8.1f\n",
           wS, ((double) wB)/wW, t2 - t1, cost) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   check that we have a valid vertex separator
   -------------------------------------------
*/
if ( 1 != GPart_validVtxSep(gpart) ) {
   fprintf(msgFile, "\n\n WHOA : separator is not valid"
           "\n\n WHOA : separator is not valid\n\n") ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------
   optionally write the tags file out to disk
   ------------------------------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   IV_writeToFile(&gpart->compidsIV, outIVfileName) ;
}
/*
   ------------------------
   free the data structures
   ------------------------
*/
GPart_free(gpart) ;
Graph_free(graph) ;
IV_free(tagsIV) ;

fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
