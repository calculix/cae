/*  testIO.c  */

#include "../DSTree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------
   test DSTree_readFromFile and DSTree_writeToFile,
   useful for translating between formatted *.dstreef
   and binary *.dstreeb files.

   created -- 97feb01, cca
   --------------------------------------------------
*/
{
char     *inDSTreeFileName, *outDSTreeFileName ;
double   t1, t2 ;
int      msglvl, rc ;
DSTree   *dstree ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.dstreef or *.dstreeb"
      "\n    outFile  -- output file, must be *.dstreef or *.dstreeb"
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
inDSTreeFileName  = argv[3] ;
outDSTreeFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inDSTreeFileName, outDSTreeFileName) ;
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
   ---------------------------
   write out the DSTree object
   ---------------------------
*/
if ( strcmp(outDSTreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = DSTree_writeToFile(dstree, outDSTreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write dstree to file %s",
           t2 - t1, outDSTreeFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from DSTree_writeToFile(%p,%s)",
           rc, dstree, outDSTreeFileName) ;
}
/*
   ----------------------
   free the DSTree object
   ----------------------
*/
DSTree_free(dstree) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
