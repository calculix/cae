/*  testIO.c  */

#include "../ZV.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------
   test ZV_readFromFile and ZV_writeToFile,
   useful for translating between formatted *.zvf
   and binary *.zvb files.

   created -- 98jun02, cca
   ----------------------------------------------
*/
{
char     *inZVFileName, *outZVFileName ;
double   t1, t2 ;
int      msglvl, rc ;
ZV       *zvobj ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.zvf or *.zvb"
      "\n    outFile  -- output file, must be *.zvf or *.zvb"
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
inZVFileName  = argv[3] ;
outZVFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inZVFileName, outZVFileName) ;
fflush(msgFile) ;
/*
   ---------------------
   read in the ZV object
   ---------------------
*/
if ( strcmp(inZVFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
zvobj = ZV_new() ;
MARKTIME(t1) ;
rc = ZV_readFromFile(zvobj, inZVFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in iv from file %s",
        t2 - t1, inZVFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ZV_readFromFile(%p,%s)",
           rc, zvobj, inZVFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading ZV object from file %s",
        inZVFileName) ;
if ( msglvl > 2 ) {
   ZV_writeForHumanEye(zvobj, msgFile) ;
} else {
   ZV_writeStats(zvobj, msgFile) ;
}
fflush(msgFile) ;
/*
   -----------------------
   write out the ZV object
   -----------------------
*/
if ( strcmp(outZVFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ZV_writeToFile(zvobj, outZVFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write iv to file %s",
           t2 - t1, outZVFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ZV_writeToFile(%p,%s)",
           rc, zvobj, outZVFileName) ;
}
/*
   ------------------
   free the ZV object
   ------------------
*/
ZV_free(zvobj) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
