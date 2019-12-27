/*  testIO.c  */

#include "../DV.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------
   test DV_readFromFile and DV_writeToFile,
   useful for translating between formatted *.dvf
   and binary *.dvb files.

   created -- 98jun02, cca
   ----------------------------------------------
*/
{
char     *inDVFileName, *outDVFileName ;
double   t1, t2 ;
int      msglvl, rc ;
DV       *dvobj ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.dvf or *.dvb"
      "\n    outFile  -- output file, must be *.dvf or *.dvb"
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
inDVFileName  = argv[3] ;
outDVFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inDVFileName, outDVFileName) ;
fflush(msgFile) ;
/*
   ---------------------
   read in the DV object
   ---------------------
*/
if ( strcmp(inDVFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
dvobj = DV_new() ;
MARKTIME(t1) ;
rc = DV_readFromFile(dvobj, inDVFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in iv from file %s",
        t2 - t1, inDVFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from DV_readFromFile(%p,%s)",
           rc, dvobj, inDVFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading DV object from file %s",
        inDVFileName) ;
if ( msglvl > 2 ) {
   DV_writeForHumanEye(dvobj, msgFile) ;
} else {
   DV_writeStats(dvobj, msgFile) ;
}
fflush(msgFile) ;
/*
   -----------------------
   write out the DV object
   -----------------------
*/
if ( strcmp(outDVFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = DV_writeToFile(dvobj, outDVFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write iv to file %s",
           t2 - t1, outDVFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from DV_writeToFile(%p,%s)",
           rc, dvobj, outDVFileName) ;
}
/*
   ------------------
   free the DV object
   ------------------
*/
DV_free(dvobj) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
