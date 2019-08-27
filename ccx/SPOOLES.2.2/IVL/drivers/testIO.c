/*  testIO.c  */

#include "../IVL.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   test IVL_readFromFile and IVL_writeToFile,
   useful for translating between formatted *.ivlf
   and binary *.ivlb files.

   created -- 98sep05, cca
   -------------------------------------------------
*/
{
char     *inIVLfileName, *outIVLfileName ;
double   t1, t2 ;
int      msglvl, rc ;
IVL    *ivl ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.ivlf or *.ivlb"
      "\n    outFile  -- output file, must be *.ivlf or *.ivlb"
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
inIVLfileName  = argv[3] ;
outIVLfileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inIVLfileName, outIVLfileName) ;
fflush(msgFile) ;
/*
   ----------------------
   read in the IVL object
   ----------------------
*/
if ( strcmp(inIVLfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
ivl = IVL_new() ;
IVL_init1(ivl, IVL_CHUNKED, 0) ;
MARKTIME(t1) ;
rc = IVL_readFromFile(ivl, inIVLfileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in ivl from file %s",
        t2 - t1, inIVLfileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IVL_readFromFile(%p,%s)",
           rc, ivl, inIVLfileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IVL object from file %s",
        inIVLfileName) ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(ivl, msgFile) ;
} else {
   IVL_writeStats(ivl, msgFile) ;
}
fflush(msgFile) ;
/*
   ------------------------
   write out the IVL object
   ------------------------
*/
if ( strcmp(outIVLfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IVL_writeToFile(ivl, outIVLfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write ivl to file %s",
           t2 - t1, outIVLfileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IVL_writeToFile(%p,%s)",
           rc, ivl, outIVLfileName) ;
}
/*
   -------------------
   free the IVL object
   -------------------
*/
IVL_free(ivl) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
