/*  testIO.c  */

#include "../ETree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   test ETree_readFromFile and ETree_writeToFile,
   useful for translating between formatted *.etreef
   and binary *.etreeb files.

   created -- 96may02, cca
   -------------------------------------------------
*/
{
char     *inETreeFileName, *outETreeFileName ;
double   t1, t2 ;
int      msglvl, rc ;
ETree    *etree ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.etreef or *.etreeb"
      "\n    outFile  -- output file, must be *.etreef or *.etreeb"
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
inETreeFileName  = argv[3] ;
outETreeFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName, outETreeFileName) ;
fflush(msgFile) ;
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
if ( strcmp(inETreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
etree = ETree_new() ;
MARKTIME(t1) ;
rc = ETree_readFromFile(etree, inETreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in etree from file %s",
        t2 - t1, inETreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
           rc, etree, inETreeFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading ETree object from file %s",
        inETreeFileName) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
fflush(msgFile) ;
{
IV   *oldToNewIV = ETree_oldToNewVtxPerm(etree) ;
IV_writeForHumanEye(oldToNewIV, msgFile) ;
}
/*
   --------------------------
   write out the ETree object
   --------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ETree_writeToFile(etree, outETreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outETreeFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_writeToFile(%p,%s)",
           rc, etree, outETreeFileName) ;
}
/*
   ---------------------
   free the ETree object
   ---------------------
*/
ETree_free(etree) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
