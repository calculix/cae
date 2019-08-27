/*  test_r1upd.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------
   test the Chv_r1upd() method.
   the program's output is a matlab file
   to check correctness of the code.

   created -- 98apr30, cca
   -------------------------------------
*/
{
Chv     *chv ;
double   imag, real, t1, t2 ;
double   *entries ;
Drand    *drand ;
FILE     *msgFile ;
int      ii, irow, jcol, msglvl, ncol, nD, nent, nL, nrow, nU, 
         rc, seed, symflag, tag, type ;
int      *colind, *rowind ;

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
           "\n       0 --> hermitian"
           "\n       1 --> symmetric"
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
nL = nU ;
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
if ( symflag <= 2 && nL != nU ) {
   fprintf(stderr, "\n invalid input"
      "\n symflag = %d, nL = %d, nU = %d", symflag, nL, nU) ;
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
   for ( irow = 0 ; irow < nD ; irow++ ) {
      Chv_complexEntry(chv, irow, irow, &real, &imag) ;
      Chv_setComplexEntry(chv, irow, irow, real, 0.0) ;
   }
}
fprintf(msgFile, "\n %% matrix entries") ;
Chv_writeForMatlab(chv, "a", msgFile) ;
/*
   ---------------------------------------
   write out matlab code for rank-1 update
   ---------------------------------------
*/
fprintf(msgFile,
        "\n nD = %d ;"
        "\n nL = %d ;"
        "\n nU = %d ;"
        "\n nrow = nD + nL ;"
        "\n ncol = nD + nU ;"
        "\n b = a ; "
        "\n d = a(1,1) ;"
        "\n l = a(2:nrow,1) / d ; "
        "\n u = a(1,2:ncol) ; "
        "\n b(2:nrow,2:ncol) = a(2:nrow,2:ncol) - l * u ; "
        "\n u = u / d ; "
        "\n b(1,1) = d ; "
        "\n b(1,2:ncol) = u ; "
        "\n b(2:nrow,1) = l ; ",
        nD, nL, nU) ;
if ( nL > 0 && nU > 0 ) {
   fprintf(msgFile, "\n b(nD+1:nrow,nD+1:ncol) = 0.0 ;") ;
}
/*
   -------------------------
   perform the rank-1 update
   -------------------------
*/
rc = Chv_r1upd(chv) ;
/*
fprintf(msgFile, "\n raw entries vector") ;
DVfprintf(msgFile, 2*nent, entries) ;
*/
fprintf(msgFile, "\n %% matrix entries after update") ;
Chv_writeForMatlab(chv, "c", msgFile) ;
fprintf(msgFile, "\n maxerr = max(max(abs(c-b)))") ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Chv_free(chv) ;
Drand_free(drand) ;
           
fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
