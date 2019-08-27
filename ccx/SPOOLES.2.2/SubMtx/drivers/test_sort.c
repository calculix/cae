/*  test_sort.c  */

#include "../SubMtx.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------
   test the sort methods.

   created -- 98apr15, cca
   -----------------------
*/
{
SubMtx     *mtxA ;
double   t1, t2 ;
Drand    *drand ;
FILE     *msgFile ;
int      mode, msglvl, ncolA, nentA, nrowA, seed, type ;
int      *colind, *ivtemp, *rowind ;

if ( argc != 9 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile type mode "
           "\n         nrowA ncolA nentA seed"
           "\n    msglvl  -- message level"
           "\n    msgFile -- message file"
           "\n    type    -- type of entries"
           "\n       1 -- real"
           "\n       2 -- complex"
           "\n    mode    -- type of matrix A"
           "\n       0 -- dense stored by rows"
           "\n       1 -- dense stored by columns"
           "\n       2 -- sparse stored by rows"
           "\n       3 -- sparse stored by columns"
           "\n    nrowA   -- # of rows in matrix A"
           "\n    ncolA   -- # of columns in matrix A"
           "\n    nentA   -- # of entries in matrix A"
           "\n    seed    -- random number seed"
           "\n", argv[0]) ;
   return(0) ;
}
if ( (msglvl = atoi(argv[1])) < 0 ) {
   fprintf(stderr, "\n message level must be positive\n") ;
   exit(-1) ;
}
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n unable to open file %s\n", argv[2]) ;
   return(-1) ;
}
type  = atoi(argv[3]) ;
mode  = atoi(argv[4]) ;
nrowA = atoi(argv[5]) ;
ncolA = atoi(argv[6]) ;
nentA = atoi(argv[7]) ;
seed  = atoi(argv[8]) ;
fprintf(msgFile, "\n %% %s:"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% mode    = %d"
        "\n %% nrowA   = %d"
        "\n %% ncolA   = %d"
        "\n %% nentA   = %d"
        "\n %% seed    = %d",
        argv[0], msglvl, argv[2], type, mode, 
        nrowA, ncolA, nentA, seed) ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if ( nrowA <= 0 
   || ncolA <= 0 
   || nentA > nrowA*ncolA ) {
   fprintf(stderr, "\n invalid input\n") ;
   exit(-1) ;
}
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
Drand_init(drand) ;
Drand_setSeed(drand, seed) ;
Drand_setNormal(drand, 0.0, 1.0) ;
/*
   -----------------------------------
   initialize the A matrix SubMtx object
   -----------------------------------
*/
mtxA  = SubMtx_new() ;
SubMtx_initRandom(mtxA, type, mode, 0, 0, nrowA, ncolA, nentA, seed) ;
SubMtx_rowIndices(mtxA, &nrowA, &rowind) ;
ivtemp = IVinit(nrowA + 5, -1) ;
IVramp(nrowA + 5, ivtemp, 0, 1) ;
IVshuffle(nrowA + 5, ivtemp, ++seed) ;
IVcopy(nrowA, rowind, ivtemp) ;
IVfree(ivtemp) ;
SubMtx_columnIndices(mtxA, &ncolA, &colind) ;
ivtemp = IVinit(ncolA + 5, -1) ;
IVramp(ncolA + 5, ivtemp, 0, 1) ;
IVshuffle(ncolA + 5, ivtemp, ++seed) ;
IVcopy(ncolA, colind, ivtemp) ;
IVfree(ivtemp) ;
SubMtx_writeToFile(mtxA, "temp.submtxb") ;
SubMtx_writeToFile(mtxA, "temp.submtxf") ;
/*
   --------------------
   print out the matrix
   --------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% A SubMtx object") ;
   fprintf(msgFile, "\n A = zeros(%d,%d) ;", nrowA+5, ncolA+5) ;
   SubMtx_writeForMatlab(mtxA, "A", msgFile) ;
   fflush(msgFile) ;
}
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_SPARSE_ROWS :
/*
   --------------------------------
   sort the rows in ascending order
   --------------------------------
*/
   SubMtx_sortRowsUp(mtxA) ;
   break ;
}
switch ( mode ) {
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_SPARSE_COLUMNS :
/*
   -----------------------------------
   sort the columns in ascending order
   -----------------------------------
*/
   SubMtx_sortColumnsUp(mtxA) ;
   break ;
}
/*
   --------------------------
   print out the matrix again
   --------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% B SubMtx object") ;
   fprintf(msgFile, "\n B = zeros(%d,%d) ;", nrowA+5, ncolA+5) ;
   SubMtx_writeForMatlab(mtxA, "B", msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------
   check with matlab
   -----------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n maxabs = max(max(abs(A - B))) ") ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
SubMtx_free(mtxA) ;
Drand_free(drand) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
