/*  test_maxabs.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------------
   test the Chv_maxabsInDiagonal(), Chv_maxabsInRow() 
   and Chv_maxabsInColumn() methods.
   the program's output is a matlab file
   to check correctness of the code.

   created -- 98apr22, cca
   ----------------------------------------------------
*/
{
Chv     *chv ;
double   imag, maxval, real, t1, t2 ;
double   *colmaxes, *entries, *rowmaxes ;
Drand    *drand ;
FILE     *msgFile ;
int      icol, ii, irow, jcol, jrow, msglvl, ncol, nD, nent, 
         nL, nrow, nU, rc, seed, symflag, tag, type ;
int      *colind, *colmark, *rowind, *rowmark ;

if ( argc != 8 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile nD nU type symflag seed "
           "\n    msglvl  -- message level"
           "\n    msgFile -- message file"
           "\n    nD      -- # of rows and columns in the (1,1) block"
           "\n    nU      -- # of columns in the (1,2) block"
           "\n    type    -- entries type"
           "\n       1 --> real"
           "\n       2 --> complex"
           "\n    symflag -- symmetry flag"
           "\n       0 --> symmetric"
           "\n       1 --> hermitian"
           "\n       2 --> nonsymmetric "
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
nD      = atoi(argv[3]) ;
nU      = atoi(argv[4]) ;
type    = atoi(argv[5]) ;
symflag = atoi(argv[6]) ;
seed    = atoi(argv[7]) ;
fprintf(msgFile, "\n %% testChv:"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% nD      = %d"
        "\n %% nU      = %d"
        "\n %% type    = %d"
        "\n %% symflag = %d"
        "\n %% seed    = %d",
        msglvl, argv[2], nD, nU, type, symflag, seed) ;
nL   = nU ;
nrow = nD + nL ;
ncol = nD + nU ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if (  nD <= 0 || nL < 0 || nU < 0 
   || symflag < 0 || symflag > 3 ) {
   fprintf(stderr, "\n invalid input"
      "\n nD = %d, nL = %d, nU = %d, symflag = %d\n",
           nD, nL, nU, symflag) ;
   exit(-1) ;
}
fprintf(msgFile,
        "\n nD = %d ;"
        "\n nL = %d ;"
        "\n nU = %d ;"
        "\n nrow = nD + nL ;"
        "\n ncol = nD + nU ;",
        nD, nL, nU) ;
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
   ----------------------------
   initialize the Chv object
   ----------------------------
*/
MARKTIME(t1) ;
chv = Chv_new() ;
Chv_init(chv, 0, nD, nL, nU, type, symflag) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize chv object",
        t2 - t1) ;
fflush(msgFile) ;
Chv_columnIndices(chv, &ncol, &colind) ;
IVramp(ncol, colind, 0, 1) ;
if ( CHV_IS_NONSYMMETRIC(chv) ) {
   Chv_rowIndices(chv, &nrow, &rowind) ;
   IVramp(nrow, rowind, 0, 1) ;
}
/*
   ------------------------------------
   load the entries with random entries
   ------------------------------------
*/
nent    = Chv_nent(chv) ;
entries = Chv_entries(chv) ;
if ( CHV_IS_REAL(chv) ) {
   Drand_fillDvector(drand, nent, entries) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   Drand_fillDvector(drand, 2*nent, entries) ;
}
if ( CHV_IS_HERMITIAN(chv) ) {
/*
   ---------------------------------------------------------
   hermitian example, set imaginary part of diagonal to zero
   ---------------------------------------------------------
*/
   for ( irow = 0 ; irow < nD ; irow++ ) {
      Chv_complexEntry(chv, irow, irow, &real, &imag) ;
      Chv_setComplexEntry(chv, irow, irow, real, 0.0) ;
   }
}
fprintf(msgFile, "\n %% matrix entries") ;
Chv_writeForMatlab(chv, "a", msgFile) ;
/*
   -----------------------------
   find the row and column maxes
   -----------------------------
*/
rowmaxes = DVinit(nrow, 0.0) ;
colmaxes = DVinit(nrow, 0.0) ;
for ( irow = 0 ; irow < nD ; irow++ ) {
   jcol = Chv_maxabsInRow(chv, irow, &maxval) ;
   rowmaxes[irow] = maxval ;
}
for ( jcol = 0 ; jcol < nD ; jcol++ ) {
   irow = Chv_maxabsInColumn(chv, jcol, &maxval) ;
   colmaxes[jcol] = maxval ;
}
fprintf(msgFile, "\n\n rowmaxes = [ ...") ;
for ( irow = 0 ; irow < nD ; irow++ ) {
   fprintf(msgFile, "\n %20.12e", rowmaxes[irow]) ;
}
fprintf(msgFile, " ] ;") ;
fprintf(msgFile, "\n\n colmaxes = [ ...") ;
for ( jcol = 0 ; jcol < nD ; jcol++ ) {
   fprintf(msgFile, "\n %20.12e", colmaxes[jcol]) ;
}
fprintf(msgFile, " ] ;") ;
/*
   -----------------
   check with matlab
   -----------------
*/
fprintf(msgFile,
        "\n\n for irow = 1:nD"
        "\n   maxval = max(abs(a(irow,:))) ;"
        "\n   rmaxes(irow) = maxval ;"
        "\nend"
        "\nrowerror = norm(rmaxes' - rowmaxes) "
        "\nfor jcol = 1:nD"
        "\n   maxval = max(abs(a(:,jcol))) ;"
        "\n   cmaxes(jcol) = maxval ;"
        "\nend"
        "\ncolerror = norm(cmaxes' - colmaxes) ") ;
/*
   -----------------------------------------------------
   find the row and column maxes of just the (1,1) block
   -----------------------------------------------------
*/
rowmark = IVinit(nD, -1) ;
colmark = IVinit(nD, -1) ;
tag     = -1 ;
for ( irow = 0 ; irow < nD ; irow++ ) {
   jcol = Chv_maxabsInRow11(chv, irow, colmark, tag, &maxval) ;
   rowmaxes[irow] = maxval ;
}
for ( jcol = 0 ; jcol < nD ; jcol++ ) {
   irow = Chv_maxabsInColumn11(chv, jcol, rowmark, tag, &maxval) ;
   colmaxes[jcol] = maxval ;
}
fprintf(msgFile, "\n\n rowmaxes = [ ...") ;
for ( irow = 0 ; irow < nD ; irow++ ) {
   fprintf(msgFile, "\n %20.12e", rowmaxes[irow]) ;
}
fprintf(msgFile, " ] ;") ;
fprintf(msgFile, "\n\n colmaxes = [ ...") ;
for ( jcol = 0 ; jcol < nD ; jcol++ ) {
   fprintf(msgFile, "\n %20.12e", colmaxes[jcol]) ;
}
fprintf(msgFile, " ] ;") ;
/*
   -----------------
   check with matlab
   -----------------
*/
fprintf(msgFile,
        "\n\n for irow = 1:nD"
        "\n   maxval = max(abs(a(irow,1:nD))) ;"
        "\n   rmaxes(irow) = maxval ;"
        "\nend"
        "\nrow11error = norm(rmaxes' - rowmaxes) "
        "\nfor jcol = 1:nD"
        "\n   maxval = max(abs(a(1:nD,jcol))) ;"
        "\n   cmaxes(jcol) = maxval ;"
        "\nend"
        "\ncol11error = norm(cmaxes' - colmaxes) ") ;
/*
   ---------------------------------------------
   find the diagonal max of just the (1,1) block
   ---------------------------------------------
*/
jcol = Chv_maxabsInDiagonal11(chv, colmark, tag, &maxval) ;
fprintf(msgFile, "\n\n maxval = %20.12e ;", maxval) ;
/*
   -----------------
   check with matlab
   -----------------
*/
fprintf(msgFile,
        "\n\n maxabs = abs(a(1,1)) ;"
        "\nfor irow = 2:nD"
        "\n   val = abs(a(irow,irow)) ;"
        "\n   if maxabs < val"
        "\n      maxabs = val ;"
        "\n   end "
        "\nend"
        "\ndiag11error = abs(maxabs - maxval)") ;
fprintf(msgFile,
        "\n [ rowerror colerror row11error col11error diag11error]") ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Chv_free(chv) ;
Drand_free(drand) ;
DVfree(rowmaxes) ;
DVfree(colmaxes) ;
IVfree(rowmark) ;
IVfree(colmark) ;
           
fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
