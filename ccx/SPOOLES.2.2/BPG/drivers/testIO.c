/*  testIO.c  */

#include "../BPG.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   test BPG_readFromFile and BPG_writeToFile,
   useful for translating between formatted *.bpgf
   and binary *.bpgb files.

   created -- 96mar08, cca
   -------------------------------------------------
*/
{
char     *inBPGFileName, *outBPGFileName ;
double   t1, t2 ;
int      msglvl, rc ;
BPG      *bpg ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.bpgf or *.bpgb"
      "\n    outFile  -- output file, must be *.bpgf or *.bpgb"
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
inBPGFileName  = argv[3] ;
outBPGFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inBPGFileName, outBPGFileName) ;
fflush(msgFile) ;
/*
   ----------------------
   read in the BPG object
   ----------------------
*/
if ( strcmp(inBPGFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
bpg = BPG_new() ;
MARKTIME(t1) ;
rc = BPG_readFromFile(bpg, inBPGFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inBPGFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from BPG_readFromFile(%p,%s)",
           rc, bpg, inBPGFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading BPG object from file %s",
        inBPGFileName) ;
if ( msglvl > 2 ) {
   BPG_writeForHumanEye(bpg, msgFile) ;
} else {
   BPG_writeStats(bpg, msgFile) ;
}
fflush(msgFile) ;
/*
   ------------------------
   write out the BPG object
   ------------------------
*/
if ( strcmp(outBPGFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = BPG_writeToFile(bpg, outBPGFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write graph to file %s",
           t2 - t1, outBPGFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from BPG_writeToFile(%p,%s)",
           rc, bpg, outBPGFileName) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
