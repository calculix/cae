/*  testTwoSetViaBKL.c  */

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
   -----------------------------------------------
   read in a Graph, 
   read in a components IV object
   find the two set partition via BKL

   created -- 96oct21, cca
   -----------------------------------------------
*/
{
char        *inGraphFileName, *inIVfileName, *msgFileName,
            *outIVfileName ;
double      alpha, t1, t2 ;
double      cpus[3] ;
FILE        *msgFile ;
GPart       *gpart ;
Graph       *graph ;
int         ierr, msglvl, nvtx, rc, seed ;
IV          *compidsIV ;

if ( argc != 8 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile inIVfileName "
      "\n           seed alpha outIVfileName"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inGraphFile -- input file, must be *.graphf or *.graphb"
      "\n    inIVfile    -- input file, must be *.ivf or *.ivb"
      "\n    seed        -- random number seed"
      "\n    alpha       -- partition cost parameter"
      "\n    outIVfile   -- output file, must be *.ivf or *.ivb"
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
seed            = atoi(argv[5]) ;
alpha           = atof(argv[6]) ;
outIVfileName   = argv[7] ;
fprintf(msgFile, 
        "\n %s : "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inGraphFile -- %s" 
        "\n inIVfile    -- %s" 
        "\n seed        -- %d" 
        "\n alpha       -- %f" 
        "\n outIVfile   -- %s" 
        "\n", argv[0], msglvl, msgFileName, inGraphFileName, 
        inIVfileName, seed, alpha, outIVfileName) ;
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
   ---------------------
   read in the IV object
   ---------------------
*/
if ( strcmp(inIVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
compidsIV = IV_new() ;
IV_setDefaultFields(compidsIV) ;
if ( (rc = IV_readFromFile(compidsIV, inIVfileName)) != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
        rc, compidsIV, inIVfileName) ;
   exit(-1) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in compidsIV from file %s", 
        t2 - t1, inIVfileName) ;
if ( msglvl < 4 ) {
   IV_writeStats(compidsIV, msgFile) ;
   fflush(msgFile) ;
} else {
   IV_writeForHumanEye(compidsIV, msgFile) ;
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
IV_copy(&gpart->compidsIV, compidsIV) ;
GPart_setMessageInfo(gpart, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : GPart initialized", t2 - t1) ;
/*
   --------------------------------
   find an initial bisector via BKL
   --------------------------------
*/
MARKTIME(t1) ;
GPart_TwoSetViaBKL(gpart, alpha, seed, cpus) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : generate bisector",
        t2 - t1) ;
fprintf(msgFile, "\n initial component weights :") ;
IVfp80(msgFile, 1+gpart->ncomp, IV_entries(&gpart->cweightsIV), 
       30, &ierr) ;
fflush(msgFile) ;
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
   ------------------------------------------------
   optionally write the compids IV object to a file
   ------------------------------------------------
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
IV_free(compidsIV) ;

fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
