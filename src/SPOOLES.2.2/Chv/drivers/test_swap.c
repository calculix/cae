/*  test_swap.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------
   test the Chv_swapRows(), Chv_swapColumns() and
   Chv_swapRowsAndColumns() methods.
   the program's output is a matlab file
   to check correctness of the code.

   created -- 98apr18, cca
   ------------------------------------------------
*/
{
Chv      *chv ;
double   imag, real, t1, t2 ;
double   *entries ;
Drand    *drand ;
FILE     *msgFile ;
int      icol, ii, irow, jcol, jrow, msglvl, ncol, nD, nent, 
         nL, nrow, nU, rc, seed, symflag, tag, type ;
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
           "\n       0 --> symmetric"
           "\n       1 --> hermitian"
           "\n       2 --> nonsymmetric"
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
   || !(type == SPOOLES_REAL || type == SPOOLES_COMPLEX) 
   || !(symflag == SPOOLES_SYMMETRIC 
       || symflag == SPOOLES_HERMITIAN 
       || symflag == SPOOLES_NONSYMMETRIC) ) {
   fprintf(stderr, "\n invalid input"
      "\n nD = %d, nL = %d, nU = %d, type = %d, symflag = %d\n",
           nD, nL, nU, type, symflag) ;
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
   --------------------------
   initialize the Chv object
   --------------------------
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
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n raw entries vector") ;
   if ( CHV_IS_REAL(chv) ) {
      DVfprintf(msgFile, nent, entries) ;
   } else if ( CHV_IS_COMPLEX(chv) ) {
      DVfprintf(msgFile, 2*nent, entries) ;
   }
   fflush(msgFile) ;
}
if ( CHV_IS_COMPLEX(chv) && CHV_IS_HERMITIAN(chv) ) {
/*
   ---------------------------------------------------------
   hermitian example, set imaginary part of diagonal to zero
   ---------------------------------------------------------
*/
   for ( irow = 0 ; irow < nD ; irow++ ) {
      Chv_complexEntry(chv, irow, irow, &real, &imag) ;
fprintf(msgFile, "\n %% setting a_{%d,%d} = %20.12e ;",
        irow, irow, real) ;
      Chv_setComplexEntry(chv, irow, irow, real, 0.0) ;
   }
}
fprintf(msgFile, "\n %% matrix entries") ;
Chv_writeForMatlab(chv, "a", msgFile) ;
if ( CHV_IS_NONSYMMETRIC(chv) ) {
/*
   ---------------------------
   choose the two rows to swap
   ---------------------------
*/
   Drand_setUniform(drand, 0.0, nD) ;
   irow = Drand_value(drand) ;
   jrow = Drand_value(drand) ;
   while ( jrow == irow ) {
      jrow = Drand_value(drand) ;
   }
   if ( irow < 0 || irow >= nD || jrow < 0 || jrow >= nD ) {
      fprintf(stderr, "\n fatal error, irow = %d, jrow = %d\n",
              irow, jrow) ;
      exit(-1) ;
   }
   fprintf(msgFile, "\n %% swapping rows %d and %d", irow+1, jrow+1) ;
   fprintf(msgFile, 
           "\n irow = %d ;"
           "\n jrow = %d ;"
           "\n b = a ;"
           "\n b(irow,:) = a(jrow,:) ;"
           "\n b(jrow,:) = a(irow,:) ;",
           irow+1, jrow+1) ;
   fflush(msgFile) ;
   Chv_swapRows(chv, irow, jrow) ;
   fprintf(msgFile, "\n %% matrix entries") ;
   Chv_writeForMatlab(chv, "c", msgFile) ;
   fprintf(msgFile, "\n maxerrrowswap1 = norm(c - a)") ;
/*
   ------------------------------
   choose the two columns to swap
   ------------------------------
*/
   Drand_setUniform(drand, 0.0, nD) ;
   icol = Drand_value(drand) ;
   jcol = Drand_value(drand) ;
   while ( jcol == icol ) {
      jcol = Drand_value(drand) ;
   }
   if ( icol < 0 || icol >= nD || jcol < 0 || jcol >= nD ) {
      fprintf(stderr, "\n fatal error, icol = %d, jcol = %d\n",
              icol, jrow) ;
      exit(-1) ;
   }
   fprintf(msgFile, 
           "\n %% swapping columns %d and %d", icol+1, jcol+1) ;
   fprintf(msgFile, 
           "\n icol = %d ;"
           "\n jcol = %d ;"
           "\n c = b ;"
           "\n c(:,icol) = b(:,jcol) ;"
           "\n c(:,jcol) = b(:,icol) ;",
           icol+1, jcol+1) ;
   fflush(msgFile) ;
   Chv_swapColumns(chv, icol, jcol) ;
   fprintf(msgFile, "\n %% matrix entries") ;
   Chv_writeForMatlab(chv, "d", msgFile) ;
   fprintf(msgFile, "\n maxerrcolswap1 = norm(d - a)") ;
   IVramp(ncol, colind, 0, 1) ;
   IVramp(nrow, rowind, 0, 1) ;
   Chv_writeForMatlab(chv, "e", msgFile) ;
   fprintf(msgFile, "\n maxerrswap = norm(e - c)") ;
} else {
/*
   ---------------------------------------
   choose the two rows and columns to swap
   ---------------------------------------
*/
   Drand_setUniform(drand, 0.0, nD) ;
   irow = Drand_value(drand) ;
   jrow = Drand_value(drand) ;
   while ( jrow == irow ) {
      jrow = Drand_value(drand) ;
   }
   if ( irow < 0 || irow >= nD || jrow < 0 || jrow >= nD ) {
      fprintf(stderr, "\n fatal error, irow = %d, jrow = %d\n",
              irow, jrow) ;
      exit(-1) ;
   }
   fprintf(msgFile, 
          "\n %% swapping rows and columns %d and %d", irow+1, jrow+1) ;
   fprintf(msgFile, 
           "\n irow = %d ;"
           "\n jrow = %d ;"
           "\n b = a ;"
           "\n b(irow,:) = a(jrow,:) ;"
           "\n b(jrow,:) = a(irow,:) ;"
           "\n c = b ;"
           "\n c(:,irow) = b(:,jrow) ;"
           "\n c(:,jrow) = b(:,irow) ;",
           irow+1, jrow+1) ;
   fflush(msgFile) ;
   Chv_swapRowsAndColumns(chv, irow, jrow) ;
   fprintf(msgFile, "\n %% matrix entries") ;
   Chv_writeForMatlab(chv, "d", msgFile) ;
   fprintf(msgFile, "\n maxerrsymswap1 = norm(d - a)") ;
   IVramp(ncol, colind, 0, 1) ;
   Chv_writeForMatlab(chv, "e", msgFile) ;
   fprintf(msgFile, "\n maxerrsymswap2 = norm(e - c)") ;
}
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
