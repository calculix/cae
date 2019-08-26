/*  testIO.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   test InpMtx_readFromFile and InpMtx_writeToFile,
   useful for translating between formatted *.inpmtxf
   and binary *.inpmtxb files.

   created -- 95dec17, cca
   ---------------------------------------------------
*/
{
int      msglvl, rc ;
InpMtx   *inpmtx ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : testIO msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.inpmtxf or *.inpmtxb"
      "\n    outFile  -- output file, must be *.inpmtxf or *.inpmtxb"
      "\n") ;
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
fprintf(msgFile, 
        "\n testIO "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        msglvl, argv[2], argv[3], argv[4]) ;
fflush(msgFile) ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
inpmtx = InpMtx_new() ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n after setting default fields") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------
   read in the InpMtx object
   --------------------------
*/
if ( strcmp(argv[3], "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
rc = InpMtx_readFromFile(inpmtx, argv[3]) ;
fprintf(msgFile, "\n return value %d from InpMtx_readFromFile(%p,%s)",
        rc, inpmtx, argv[3]) ;
if ( rc != 1 ) {
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading InpMtx object from file %s",
           argv[3]) ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   change the storage mode to vectors
   ----------------------------------
*/
InpMtx_changeCoordType(inpmtx, INPMTX_BY_ROWS) ;
InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after changing storage mode") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the InpMtx object
   ---------------------------
*/
if ( strcmp(argv[4], "none") != 0 ) {
   rc = InpMtx_writeToFile(inpmtx, argv[4]) ;
   fprintf(msgFile, 
           "\n return value %d from InpMtx_writeToFile(%p,%s)",
           rc, inpmtx, argv[4]) ;
}
/*
   ---------------
   free the object
   ---------------
*/
InpMtx_free(inpmtx) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
