/*  testIO.c  */

#include "../Perm.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   test Perm_readFromFile and Perm_writeToFile,
   useful for translating between formatted *.permf
   and binary *.permb files.

   created -- 96may02, cca
   -------------------------------------------------
*/
{
char     *inPermFileName, *outPermFileName ;
double   t1, t2 ;
int      msglvl, rc ;
Perm     *perm ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.permf or *.permb"
      "\n    outFile  -- output file, must be *.permf or *.permb"
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
outPermFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inPermFileName, outPermFileName) ;
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
 -  ------------------------
   write out the Perm object
  - ------------------------
*/
if ( strcmp(outPermFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = Perm_writeToFile(perm, outPermFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write perm to file %s",
           t2 - t1, outPermFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Perm_writeToFile(%p,%s)",
           rc, perm, outPermFileName) ;
}
/*
   --------------------
   free the Perm object
   --------------------
*/
Perm_free(perm) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
