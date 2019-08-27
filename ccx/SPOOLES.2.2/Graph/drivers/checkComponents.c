/*  checkComponents.c  */

#include "../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------
   read in a Graph and check for connected components

   created -- 96fmar02, cca
   --------------------------------------------------
*/
{
char     *inGraphFileName ;
double   t1, t2 ;
int      icomp, msglvl, ncomp, rc, v ;
int      *counts, *weights ;
IV       *mapIV ;
Graph    *g ;
FILE     *msgFile ;

if ( argc != 4 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inGraphFile -- input file, must be *.graphf or *.graphb"
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
inGraphFileName = argv[3] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName) ;
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
g = Graph_new() ;
MARKTIME(t1) ;
rc = Graph_readFromFile(g, inGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, g, inGraphFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading Graph object from file %s",
        inGraphFileName) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(g, msgFile) ;
} else {
   Graph_writeStats(g, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------------------------
   get the map from vertices to components
   ---------------------------------------
*/
MARKTIME(t1) ;
mapIV = Graph_componentMap(g) ;
MARKTIME(t2) ;
ncomp = 1 + IV_max(mapIV) ;
fprintf(msgFile, "\n CPU %9.5f : find %d components",
        t2 - t1, ncomp) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n component map") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
}
/*
   ------------------------------------
   get the statistics on the components
   ------------------------------------
*/
counts  = IVinit(ncomp, 0) ;
weights = IVinit(ncomp, 0) ;
MARKTIME(t1) ;
Graph_componentStats(g, IV_entries(mapIV), counts, weights) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : compute component statistics",
        t2 - t1) ;
fprintf(msgFile, "\n\n component   # vertices  weight") ;
for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
   fprintf(msgFile, "\n %7d    %7d    %7d", 
           icomp, counts[icomp], weights[icomp]) ;
}
/*
   ----------------
   free the storage
   ----------------
*/
Graph_free(g)   ;
IV_free(mapIV)  ;
IVfree(counts)  ;
IVfree(weights) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
