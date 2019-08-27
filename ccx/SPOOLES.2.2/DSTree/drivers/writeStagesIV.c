/*  writeStagesIV.c  */

#include "../DSTree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------
   read in a DSTree object, create a stagesIV object
   and write it to a file.

   created -- 97feb01, cca
   --------------------------------------------------
*/
{
char     *inDSTreeFileName, *outIVfileName, *type ;
double   t1, t2 ;
int      msglvl, rc ;
IV       *stagesIV ;
DSTree   *dstree ;
FILE     *msgFile ;

if ( argc != 6 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile type outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.dstreef or *.dstreeb"
      "\n    type     -- type of stages vector"
      "\n    outFile  -- output file, must be *.ivf or *.ivb"
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
inDSTreeFileName = argv[3] ;
type             = argv[4] ;
outIVfileName    = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inDSTreeFileName, outIVfileName) ;
fflush(msgFile) ;
/*
   -------------------------
   read in the DSTree object
   -------------------------
*/
if ( strcmp(inDSTreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
dstree = DSTree_new() ;
MARKTIME(t1) ;
rc = DSTree_readFromFile(dstree, inDSTreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in dstree from file %s",
        t2 - t1, inDSTreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from DSTree_readFromFile(%p,%s)",
           rc, dstree, inDSTreeFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading DSTree object from file %s",
        inDSTreeFileName) ;
if ( msglvl > 2 ) {
   DSTree_writeForHumanEye(dstree, msgFile) ;
} else {
   DSTree_writeStats(dstree, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------
   get the stages vector
   ---------------------
*/
if ( strcmp(type, "ND") == 0 ) {
   stagesIV = DSTree_NDstages(dstree) ;
} else if ( strcmp(type, "ND2") == 0 ) {
   stagesIV = DSTree_ND2stages(dstree) ;
} else if ( strcmp(type, "MS2") == 0 ) {
   stagesIV = DSTree_MS2stages(dstree) ;
} else if ( strcmp(type, "MS3") == 0 ) {
   stagesIV = DSTree_MS3stages(dstree) ;
} else {
   stagesIV = NULL ;
}
/*
   ---------------------------
   write out the DSTree object
   ---------------------------
*/
if ( stagesIV != NULL && strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(stagesIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write dstree to file %s",
           t2 - t1, outIVfileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from IV_writeToFile(%p,%s)",
              rc, stagesIV, outIVfileName) ;
   }
}
/*
   ----------------------
   free the DSTree object
   ----------------------
*/
DSTree_free(dstree) ;
if ( stagesIV != NULL ) {
   IV_free(stagesIV) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
