/*  testIO.c  */

#include "../IV.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------
   test IV_readFromFile and IV_writeToFile,
   useful for translating between formatted *.ivf
   and binary *.ivb files.

   created -- 97dec13, cca
   ----------------------------------------------
*/
{
char     *inIVFileName, *outIVFileName ;
double   t1, t2 ;
int      msglvl, rc ;
IV       *ivobj ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.ivf or *.ivb"
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
inIVFileName  = argv[3] ;
outIVFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inIVFileName, outIVFileName) ;
fflush(msgFile) ;
/*
   ---------------------
   read in the IV object
   ---------------------
*/
if ( strcmp(inIVFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
ivobj = IV_new() ;
MARKTIME(t1) ;
rc = IV_readFromFile(ivobj, inIVFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in iv from file %s",
        t2 - t1, inIVFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
           rc, ivobj, inIVFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IV object from file %s",
        inIVFileName) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(ivobj, msgFile) ;
} else {
   IV_writeStats(ivobj, msgFile) ;
}
fflush(msgFile) ;
/*
   -----------------------
   write out the IV object
   -----------------------
*/
if ( strcmp(outIVFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(ivobj, outIVFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write iv to file %s",
           t2 - t1, outIVFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_writeToFile(%p,%s)",
           rc, ivobj, outIVFileName) ;
}
/*
   ------------------
   free the IV object
   ------------------
*/
IV_free(ivobj) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
