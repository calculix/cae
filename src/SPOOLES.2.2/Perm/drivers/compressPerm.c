/*  compressPerm.c  */

#include "../Perm.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------------
   read in a permutation and an equivalence map vector.
   created the permutation for the compressed graph.

   created -- 96may02, cca
   ----------------------------------------------------
*/
{
char     *inIVfileName, *inPermFileName, *outPermFileName ;
double   t1, t2 ;
int      msglvl, rc ;
IV       *eqmapIV ;
Perm     *perm, *perm2 ;
FILE     *msgFile ;

if ( argc != 6 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inPermFile inMapFile outPermFile"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inPermFile  -- input file, must be *.permf or *.permb"
      "\n    inMapFile   -- input file, must be *.ivf or *.ivb"
      "\n    outPermFile -- output file, must be *.permf or *.permb"
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
inPermFileName  = argv[3] ;
inIVfileName    = argv[4] ;
outPermFileName = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inPermFile  -- %s" 
        "\n inMapFile   -- %s" 
        "\n outPermFile -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inPermFileName, inIVfileName,
        outPermFileName) ;
fflush(msgFile) ;
/*
   -----------------------
   read in the Perm object
   -----------------------
*/
if ( strcmp(inPermFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
perm = Perm_new() ;
MARKTIME(t1) ;
rc = Perm_readFromFile(perm, inPermFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in perm from file %s",
        t2 - t1, inPermFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Perm_readFromFile(%p,%s)",
           rc, perm, inPermFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading Perm object from file %s",
        inPermFileName) ;
if ( msglvl > 2 ) {
   Perm_writeForHumanEye(perm, msgFile) ;
} else {
   Perm_writeStats(perm, msgFile) ;
}
fflush(msgFile) ;
/*
   -------------------------
   read in the IV map object
   -------------------------
*/
if ( strcmp(inIVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
eqmapIV = IV_new() ;
MARKTIME(t1) ;
rc = IV_readFromFile(eqmapIV, inIVfileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in eqmapIV from file %s",
        t2 - t1, inIVfileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
           rc, eqmapIV, inIVfileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IV object from file %s",
        inIVfileName) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(eqmapIV, msgFile) ;
} else {
   IV_writeStats(eqmapIV, msgFile) ;
}
fflush(msgFile) ;
/*
   -------------------------------------
   get the compressed permutation object
   -------------------------------------
*/
perm2 = Perm_compress(perm, eqmapIV) ;
fprintf(msgFile, "\n\n compressed Perm object") ;
if ( msglvl > 2 ) {
   Perm_writeForHumanEye(perm2, msgFile) ;
} else {
   Perm_writeStats(perm2, msgFile) ;
}
/*
   -------------------------
   write out the Perm object
   -------------------------
*/
if ( strcmp(outPermFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = Perm_writeToFile(perm2, outPermFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write perm to file %s",
           t2 - t1, outPermFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Perm_writeToFile(%p,%s)",
           rc, perm2, outPermFileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
Perm_free(perm) ;
IV_free(eqmapIV) ;
Perm_free(perm2) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
