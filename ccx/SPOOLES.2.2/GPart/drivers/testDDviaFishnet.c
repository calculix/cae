/*  testDDviaFishnet.c  */

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
   ---------------------------------------------------------------
   read in a Graph and generate a domain decomposition via fishnet

   created -- 96oct21, cca
   ---------------------------------------------------------------
*/
{
char        *inGraphFileName, *msgFileName, *outIVfileName ;
double      t1, t2 ;
FILE        *msgFile ;
float       freeze ;
GPart       *gpart ;
Graph       *graph ;
int         maxweight, minweight, msglvl, nvtx, rc, seed ;

if ( argc != 9 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile "
      "\n           freeze minweight maxweight seed outIVfile"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inGraphFile -- input file, must be *.graphf or *.graphb"
      "\n    freeze      -- freeze factor, try 4"
      "\n    minweight   -- minimum domain weight"
      "\n    maxweight   -- maximum domain weight"
      "\n    seed        -- random number seed"
      "\n    outIVfile   -- output file for vertex component ids,"
      "\n                   must be *.ivf or *.ivb"
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
freeze          = atof(argv[4]) ;
minweight       = atoi(argv[5]) ;
maxweight       = atoi(argv[6]) ;
seed            = atoi(argv[7]) ;
outIVfileName   = argv[8] ;
fprintf(msgFile, 
        "\n %s : "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inGraphFile -- %s" 
        "\n freeze      -- %f" 
        "\n minweight   -- %d" 
        "\n maxweight   -- %d" 
        "\n seed        -- %d" 
        "\n outIVfile   -- %s" 
        "\n", argv[0], msglvl, msgFileName, inGraphFileName, freeze, 
        minweight, maxweight, seed, outIVfileName) ;
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
   -----------------------
   create the GPart object
   -----------------------
*/
MARKTIME(t1) ;
gpart = GPart_new() ;
GPart_init(gpart, graph) ;
GPart_setMessageInfo(gpart, msglvl, msgFile) ;
MARKTIME(t2) ;
/*
   ---------------------------------------------
   generate the domain decomposition via Fishnet
   ---------------------------------------------
*/
MARKTIME(t1) ;
GPart_DDviaFishnet(gpart, freeze, minweight, maxweight, seed) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : generate domain decomposition",
        t2 - t1) ;
/*
   -------------------------------------------
   check that we have a valid vertex separator
   -------------------------------------------
*/
if ( 1 != GPart_validVtxSep(gpart) ) {
   fprintf(msgFile, "\n\n WHOA : multisector is not valid\n\n") ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   get the boundary weights IV object
   ----------------------------------
*/
if ( msglvl > 1 ) {
   int   icomp, ncomp ;
   int   *bnd, *cweights ;
   IV    *bndIV ;

   bndIV = GPart_bndWeightsIV(gpart) ;
   bnd = IV_entries(bndIV) ;
   ncomp = gpart->ncomp ;
   cweights = IV_entries(&gpart->cweightsIV) ;
   fprintf(stdout, 
           "\n %% component id, component weight, boundary weight"
           "\n data = [ ...") ;
   for ( icomp = 0 ; icomp <= ncomp ; icomp++ ) {
      fprintf(stdout, "\n %8d %8d %8d", 
              icomp, cweights[icomp], bnd[icomp]) ;
   }
   fprintf(stdout, " ] ") ;
   IV_free(bndIV) ;
}
/*
   ------------------------------------------
   optionally write out the compids IV object
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

fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
