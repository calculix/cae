/*  readAIJ.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------
   read in Dave Young format
   construct a InpMtx object and
   write it out to a file

   created -- 97oct17, cca
   -------------------------------
*/
{
char     *inFileName, *outFileName ;
double   ent ;
InpMtx   *inpmtx ;
FILE     *inputFile, *msgFile ;
int      icol, irow, jcol, jrow, 
         msglvl, ncol, nentInRow, nrow, rc ;
int      *ivec1, *ivec2 ;

if ( argc != 5 ) {
   fprintf(stdout, 
  "\n\n usage : readAIJ2 msglvl msgFile inputFile outFile "
   "\n    msglvl    -- message level"
   "\n    msgFile   -- message file"
   "\n    inputFile -- input file for a(i,j) entries"
   "\n       the first line must be \"nrow ncol \""
   "\n          next lines are \"irow nentInRow\""
   "\n                         \"jcol entry \""
   "\n                         ..."
   "\n                         \"jcol entry \""
   "\n          ..." "\n       endif"
   "\n    outFile -- output file, must be *.inpmtxf or *.inpmtxb"
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
inFileName  = argv[3] ;
outFileName = argv[4] ;
fprintf(msgFile, 
        "\n readAIJ "
        "\n msglvl    -- %d" 
        "\n msgFile   -- %s" 
        "\n inputFile -- %s" 
        "\n outFile   -- %s" 
        "\n",
        msglvl, argv[2], inFileName, outFileName ) ;
fflush(msgFile) ;
/*
   ----------------------------
   open the input file and read
   #rows #columns 
   ----------------------------
*/
if ( (inputFile = fopen(inFileName, "r")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], inFileName) ;
   return(-1) ;
}
rc = fscanf(inputFile, "%d %d", &nrow, &ncol) ;
if ( rc != 2 ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n %d of 2 fields read on first line of file %s",
           argv[0], rc, inFileName) ;
   return(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n read in nrow = %d, ncol = %d", nrow, ncol) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------
   initialize the object
   set coordType = INPMTX_BY_ROWS --> row coordinates
   --------------------------------------------------
*/
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, INPMTX_BY_ROWS, SPOOLES_REAL, 10*nrow, 0) ;
/*
   -------------------------------------------------
   read in the entries and load them into the object
   -------------------------------------------------
*/
for ( irow = 0 ; irow < nrow ; irow++ ) {
   rc = fscanf(inputFile, "%d %d", &jrow, &nentInRow) ;
   if ( rc != 2 ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n %d of 2 fields read on entry %d of file %s",
              argv[0], rc, irow, inFileName) ;
      return(-1) ;
   }
   for ( icol = 0 ; icol < nentInRow ; icol++ ) {
      rc = fscanf(inputFile, "%d %lf", &jcol, &ent) ;
      if ( rc != 2 ) {
         fprintf(stderr, "\n fatal error in %s"
                 "\n %d of 2 fields read on entry %d of file %s",
                 argv[0], rc, irow, inFileName) ;
         return(-1) ;
      }
      InpMtx_inputRealEntry(inpmtx, jrow-1, jcol-1, ent) ;
   }
}
/*
   -----------------------------
   sort and compress the entries
   -----------------------------
*/
InpMtx_changeStorageMode(inpmtx, 3) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n sorted, compressed and vector form") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the InpMtx object
   ---------------------------
*/
if ( strcmp(outFileName, "none") != 0 ) {
   rc = InpMtx_writeToFile(inpmtx, outFileName) ;
   fprintf(msgFile, 
           "\n return value %d from InpMtx_writeToFile(%p,%s)",
           rc, inpmtx, outFileName) ;
}
/*
   ---------------------
   free the working data
   ---------------------
*/
InpMtx_free(inpmtx) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
