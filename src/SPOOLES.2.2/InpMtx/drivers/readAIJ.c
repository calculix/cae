/*  readAIJ.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   read in (i, j, a(i,j)) triples, 
   construct a InpMtx object and
   write it out to a file

   created -- 97oct17, cca
   ---------------------------------------------------
*/
{
char     *inFileName, *outFileName ;
InpMtx   *inpmtx ;
FILE      *inputFile, *msgFile ;
int       dataType, flag, ient, msglvl, 
          ncol, nent, nrow, rc ;
int       *ivec1, *ivec2 ;

if ( argc != 7 ) {
   fprintf(stdout, 
   "\n\n usage : readAIJ msglvl msgFile dataType inputFile outFile flag"
   "\n    msglvl    -- message level"
   "\n    msgFile   -- message file"
   "\n    dataType  -- 0 for indices only, 1 for double, 2 for complex"
   "\n    inputFile -- input file for a(i,j) entries"
   "\n       the first line must be \"nrow ncol nentries\""
   "\n       if dataType == 0 then"
   "\n          next lines are \"irow jcol\""
   "\n       else if dataType == 1 then"
   "\n          next lines are \"irow jcol entry\""
   "\n       else if dataType == 2 then"
   "\n          next lines are \"irow jcol realEntry imagEntry\""
   "\n       endif"
   "\n    outFile -- output file, must be *.inpmtxf or *.inpmtxb"
   "\n    flag    -- flag for 0-based or 1-based addressing"
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
dataType    = atoi(argv[3]) ;
inFileName  = argv[4] ;
outFileName = argv[5] ;
flag        = atoi(argv[6]) ;
fprintf(msgFile, 
        "\n readAIJ "
        "\n msglvl    -- %d" 
        "\n msgFile   -- %s" 
        "\n dataType  -- %d" 
        "\n inputFile -- %s" 
        "\n outFile   -- %s" 
        "\n flag      -- %d" 
        "\n",
        msglvl, argv[2], dataType, inFileName, outFileName, flag) ;
fflush(msgFile) ;
/*
   ----------------------------
   open the input file and read
   #rows #columns #entries
   ----------------------------
*/
if ( (inputFile = fopen(inFileName, "r")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], inFileName) ;
   return(-1) ;
}
rc = fscanf(inputFile, "%d %d %d", &nrow, &ncol, &nent) ;
if ( rc != 3 ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n %d of 3 fields read on first line of file %s",
           argv[0], rc, inFileName) ;
   return(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n read in nrow = %d, ncol = %d, nent = %d",
           nrow, ncol, nent) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------
   initialize the object
   set coordType = INPMTX_BY_ROWS --> row coordinates
   set inputMode = dataType 
   --------------------------------------------------
*/
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, INPMTX_BY_ROWS, dataType, nent, 0) ;
/*
   -------------------------------------------------
   read in the entries and load them into the object
   -------------------------------------------------
*/
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
if ( INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   for ( ient = 0 ; ient < nent ; ient++ ) {
      rc = fscanf(inputFile, "%d %d", ivec1 + ient, ivec2 + ient) ;
      if ( rc != 2 ) {
         fprintf(stderr, "\n fatal error in %s"
                 "\n %d of 2 fields read on entry %d of file %s",
                 argv[0], rc, ient, inFileName) ;
         return(-1) ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n entry %d, row %d, column %d",
                 ient, ivec1[ient], ivec2[ient]) ;
         fflush(msgFile) ;
      }
   }
} else if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double   *dvec = InpMtx_dvec(inpmtx) ;
   for ( ient = 0 ; ient < nent ; ient++ ) {
      rc = fscanf(inputFile, "%d %d %le", 
                  ivec1 + ient, ivec2 + ient, dvec + ient) ;
      if ( rc != 3 ) {
         fprintf(stderr, "\n fatal error in %s"
                 "\n %d of 3 fields read on entry %d of file %s",
                 argv[0], rc, ient, argv[3]) ;
         return(-1) ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n entry %d, row %d, column %d, value %e",
                 ient, ivec1[ient], ivec2[ient], dvec[ient]) ;
         fflush(msgFile) ;
      }
   }
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   double   *dvec = InpMtx_dvec(inpmtx) ;
   for ( ient = 0 ; ient < nent ; ient++ ) {
      rc = fscanf(inputFile, "%d %d %le %le", 
                  ivec1 + ient, ivec2 + ient, 
                  dvec + 2*ient, dvec + 2*ient+1) ;
      if ( rc != 4 ) {
         fprintf(stderr, "\n fatal error in %s"
                 "\n %d of 4 fields read on entry %d of file %s",
                 argv[0], rc, ient, argv[3]) ;
         return(-1) ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
              "\n entry %d, row %d, column %d, value %12.4e + %12.4e*i",
              ient, ivec1[ient], ivec2[ient], 
              dvec[2*ient], dvec[2*ient+1]) ;
         fflush(msgFile) ;
      }
   }
}
inpmtx->nent = nent ;
if ( flag == 1 ) {
/*
   --------------------------------------------------
   indices were in FORTRAN mode, decrement for C mode
   --------------------------------------------------
*/
   for ( ient = 0 ; ient < nent ; ient++ ) {
      ivec1[ient]-- ; ivec2[ient]-- ;
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
