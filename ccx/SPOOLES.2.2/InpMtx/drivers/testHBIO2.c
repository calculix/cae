/*  testHBIO2.c  */

#include "../../misc.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   read in a Harwell-Boeing matrix, 
   convert to a InpMtx object and write to a file.

   created -- 98sep11, cca
   ---------------------------------------------------
*/
{
char      *inFileName, *outFileName ;
double    t1, t2 ;
int       msglvl, rc ;
InpMtx    *inpmtx ;
FILE      *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be Harwell-Boeing format"
      "\n    outFile  -- output file, must be *.inpmtxf or *.inpmtxb"
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
inFileName  = argv[3] ;
outFileName = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inFileName, outFileName) ;
fflush(msgFile) ;
/*
   ---------------------------------------------
   read in the Harwell-Boeing matrix information
   ---------------------------------------------
*/
if ( strcmp(inFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
inpmtx = InpMtx_new() ;
fprintf(msgFile, "\n 1. inpmtx = %p", inpmtx) ;
rc = InpMtx_readFromHBfile(inpmtx, inFileName) ;
fprintf(msgFile, "\n 2. inpmtx = %p", inpmtx) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in matrix", t2 - t1) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n error reading in matrix") ;
   return(1) ;
}
if ( msglvl > 2 ) {
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------
   write out the InpMtx object
   ----------------------------
*/
if ( strcmp(outFileName, "none") != 0 ) {
   rc = InpMtx_writeToFile(inpmtx, outFileName) ;
   fprintf(msgFile, 
           "\n return value %d from InpMtx_writeToFile(%p,%s)",
           rc, inpmtx, outFileName) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
InpMtx_free(inpmtx) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
