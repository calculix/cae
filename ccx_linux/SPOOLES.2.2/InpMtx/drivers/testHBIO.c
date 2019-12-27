/*  testHBIO.c  */

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
char      *inFileName, *outFileName, *type ;
double    t1, t2 ;
double    *values ;
int       ierr, ii, iiend, iistart, inputMode, jcol, msglvl, ncol, 
          nnonzeros, nrhs, nrow, rc ;
int       *colptr, *colind, *rowind ;
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
readHB_info(inFileName, &nrow, &ncol, &nnonzeros, &type, &nrhs) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in information", t2 - t1) ;
fprintf(msgFile, 
        "\n matrix is %d x %d with %d entries, type = %s, nrhs = %d",
        nrow, ncol, nnonzeros, type, nrhs) ;
fflush(msgFile) ;
switch ( type[0] ) {
case 'P' :
   inputMode = INPMTX_INDICES_ONLY ;
   break ;
case 'R' :
   inputMode = SPOOLES_REAL ;
   break ;
case 'C' :
   inputMode = SPOOLES_COMPLEX ;
   break ;
default :
   fprintf(stderr, "\n fatal error in %s, type = %s"
           "\n first character must be 'P', 'R' or 'C'",
           argv[0], type) ;
   exit(-1) ;
   break ;
}
FREE(type) ;
/*
   -----------------------------
   initialize the InpMtx object
   -----------------------------
*/
MARKTIME(t1) ;
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, INPMTX_BY_COLUMNS, inputMode, nnonzeros, 0) ;
colptr = IVinit(ncol+1, -1) ;
colind = InpMtx_ivec1(inpmtx)   ;
rowind = InpMtx_ivec2(inpmtx)   ;
values = InpMtx_dvec(inpmtx)    ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize InpMtx object", t2 - t1) ;
/*
   -------------------------------
   read in the indices and entries
   -------------------------------
*/
MARKTIME(t1) ;
readHB_mat_double(inFileName, colptr, rowind, values) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in matrix entries", t2 - t1) ;
/*
   --------------------------------------------
   decrement the column offsets and row indices
   --------------------------------------------
*/
MARKTIME(t1) ;
for ( jcol = 0 ; jcol <= ncol ; jcol++ ) {
   colptr[jcol]-- ;
}
for ( ii = 0 ; ii < nnonzeros ; ii++ ) {
   rowind[ii]-- ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : decrement indices", t2 - t1) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n colptr") ;
   IVfp80(msgFile, ncol+1, colptr, 80, &ierr) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n rowind") ;
   IVfp80(msgFile, nnonzeros, rowind, 80, &ierr) ;
}
if ( values != NULL && msglvl > 3 ) {
   fprintf(msgFile, "\n\n values") ;
   DVfprintf(msgFile, nnonzeros, values) ;
}
/*
   -------------------------------------------
   fill the ivec1[] vector with column indices
   -------------------------------------------
*/
MARKTIME(t1) ;
for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
   iistart = colptr[jcol] ;
   iiend   = colptr[jcol+1] - 1 ;
   for ( ii = iistart ; ii <= iiend ; ii++ ) {
      colind[ii] = jcol ;
   }
}
InpMtx_setNent(inpmtx, nnonzeros) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : fill column indices", t2 - t1) ;
/*
   ---------------------------
   sort and convert to vectors
   ---------------------------
*/
MARKTIME(t1) ;
InpMtx_changeStorageMode(inpmtx, INPMTX_SORTED) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : sort entries", t2 - t1) ;
MARKTIME(t1) ;
InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : convert to vectors", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n InpMtx object ") ;
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
IVfree(colptr) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
