/*  test_findPivot.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   test the Chv_findPivot(), swap and update methods.
   the program's output is a matlab file
   to check correctness of the code.

   created -- 98jan24, cca
   ---------------------------------------------------
*/
{
Chv     *chv ;
double   imag, real, tau, t1, t2 ;
double   *entries ;
Drand    *drand ;
DV       *workDV ;
FILE     *msgFile ;
int      icol, ii, ipvt, irow, jcol, jpvt, jrow, msglvl, ncol, nD, 
         ndelay, nent, nL, nrow, ntest, nU, rc, pivotsize, seed, 
         symflag, tag, temp, type ;
int      *colind, *rowind ;

if ( argc != 9 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile nD nU type symflag seed tau "
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
           "\n    tau     -- bound on magnitudes of factor entries"
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
tau     = atof(argv[8]) ;
fprintf(msgFile, "\n %% testChv:"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% nD      = %d"
        "\n %% nU      = %d"
        "\n %% type    = %d"
        "\n %% symflag = %d"
        "\n %% seed    = %d"
        "\n %% tau     = %12.4e",
        msglvl, argv[2], nD, nU, type, symflag, seed, tau) ;
nL   = nU ;
nrow = nD + nL ;
ncol = nD + nU ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if (  nD <= 0 || nU < 0 
   || (symflag != SPOOLES_SYMMETRIC
   &&  symflag !=  SPOOLES_HERMITIAN
   &&  symflag !=  SPOOLES_NONSYMMETRIC) ) {
   fprintf(stderr, "\n invalid input"
      "\n nD = %d, nL = %d, nU = %d, symflag = %d\n",
           nD, nL, nU, symflag) ;
   exit(-1) ;
}
if (  (symflag ==  SPOOLES_SYMMETRIC || symflag ==  SPOOLES_HERMITIAN) 
   && nL != nU ) {
   fprintf(stderr, "\n invalid input"
      "\n symflag = %d, nL = %d, nU = %d", symflag, nL, nU) ;
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
   DVfprintf(msgFile, 2*nent, entries) ;
   fflush(msgFile) ;
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
   ------------
   find a pivot 
   ------------
*/
workDV = DV_new() ;
ndelay = 0 ;
ntest  = 0 ;
pivotsize = Chv_findPivot(chv, workDV, tau, ndelay, 
                           &irow, &jcol, &ntest) ;
fprintf(msgFile, "\n\n %% pivotsize = %d", pivotsize) ;
ipvt = irow ;
jpvt = jcol ;
if (  (symflag == SPOOLES_SYMMETRIC || symflag == SPOOLES_HERMITIAN)
   && irow > jcol ) {
   temp = irow ;
   irow = jcol ;
   jcol = temp ;
}
fprintf(msgFile, "\n\n irow = %d ; \n jcol = %d ;", irow+1, jcol+1) ;
/*
   -------------------------
   swap the rows and columns
   -------------------------
*/
if ( pivotsize == 0 ) {
   exit(0) ;
} else if ( pivotsize == 1 ) {
   fprintf(msgFile, 
           "\n b = a ;"
           "\n xtemp = b(irow,:) ;"
           "\n b(irow,:) = b(1,:) ;"
           "\n b(1,:) = xtemp ;"
           "\n xtemp = b(:,jcol) ;"
           "\n b(:,jcol) = b(:,1) ;"
           "\n b(:,1) = xtemp ;") ;
   if ( CHV_IS_SYMMETRIC(chv) || symflag == CHV_IS_HERMITIAN(chv) ) {
      Chv_swapRowsAndColumns(chv, 0, irow) ;
   } else {
      Chv_swapRows(chv, 0, irow) ;
      Chv_swapColumns(chv, 0, jcol) ;
   }
} else if ( pivotsize == 2 ) {
   if ( symflag < 2 ) {
      fprintf(msgFile, 
              "\n b = a ;"
              "\n xtemp = b(irow,:) ;"
              "\n b(irow,:) = b(1,:) ;"
              "\n b(1,:) = xtemp ;"
              "\n xtemp = b(:,irow) ;"
              "\n b(:,irow) = b(:,1) ;"
              "\n b(:,1) = xtemp ;"
              "\n xtemp = b(jcol,:) ;"
              "\n b(jcol,:) = b(2,:) ;"
              "\n b(2,:) = xtemp ;"
              "\n xtemp = b(:,jcol) ;"
              "\n b(:,jcol) = b(:,2) ;"
              "\n b(:,2) = xtemp ;") ;
      Chv_swapRowsAndColumns(chv, 0, irow) ;
      Chv_swapRowsAndColumns(chv, 1, jcol) ;
   } else {
      fprintf(stderr, "\n fatal error, symflag = %d, pvtsize = %d",
              symflag, pivotsize) ;
      exit(-1) ;
   }
}
/*
   -----------------------------------------
   check that the swap was executed properly
   -----------------------------------------
*/
fprintf(msgFile, "\n %% matrix entries") ;
Chv_writeForMatlab(chv, "c", msgFile) ;
fprintf(msgFile, "\n maxerrswap = norm(c - a)") ;
/*
   ---------------------------
   ramp the indices once again
   ---------------------------
*/
IVramp(ncol, colind, 0, 1) ;
if ( CHV_IS_NONSYMMETRIC(chv) ) {
   Chv_rowIndices(chv, &nrow, &rowind) ;
   IVramp(nrow, rowind, 0, 1) ;
}
/*
   -----------------------------------
   perform the rank-1 or rank-2 update
   -----------------------------------
*/
fprintf(msgFile, "\n\n ckeep = b ;") ;
fprintf(msgFile, "\n\n c = b ;") ;
if ( pivotsize == 1 ) {
   rc = Chv_r1upd(chv) ;
   fprintf(msgFile, 
           "\n\n d = c(1,1) ;"
           "\n l = c(2:nrow,1)/d ;"
           "\n u = c(1,2:ncol) ;") ;
   if ( nD > 1 ) {
      fprintf(msgFile, 
           "\n c(2:nrow,2:ncol) = c(2:nrow,2:ncol) - l*u ;") ;
   }
   fprintf(msgFile, 
           "\n u = u / d ;"
           "\n c(1:1,1:1) = d ; "
           "\n c(1:1,2:ncol) = u ; "
           "\n c(2:ncol,1:1) = l ; ") ;
   fprintf(msgFile, "\n c(nD+1:nrow,nD+1:ncol) = 0 ;") ;
} else {
   rc = Chv_r2upd(chv) ;
   fprintf(msgFile, 
           "\n\n d = c(1:2,1:2) ;"
           "\n l = c(3:nrow,1:2) / d ;"
           "\n u = c(1:2,3:ncol) ;") ;
   if ( nD > 2 ) {
      fprintf(msgFile, 
              "\n c(3:nrow,3:ncol) = c(3:nrow,3:ncol) - l*u ;") ;
   }
   fprintf(msgFile, 
           "\n u = d \\ u ; "
           "\n c(1:2,1:2) = d ; "
           "\n c(1:2,3:ncol) = u ; "
           "\n c(3:ncol,1:2) = l ; ") ;
   if ( nU > 0 ) {
      fprintf(msgFile, 
           "\n c(nD+1:nrow,nD+1:ncol) = 0 ;") ;
   }
}
fprintf(msgFile, "\n %% matrix entries after update") ;
Chv_writeForMatlab(chv, "f", msgFile) ;
fprintf(msgFile, "\n maxerrupd = norm(f - c)") ;
/*
   ------------------------------------------------------
   check out the maximum magnitude of elements in l and u
   ------------------------------------------------------
*/
fprintf(msgFile, "\n ipvt = %d", ipvt + 1) ;
fprintf(msgFile, "\n jpvt = %d", jpvt + 1) ;
fprintf(msgFile, "\n pivotsize = %d", pivotsize) ;
fprintf(msgFile, "\n tau = %12.4e", tau) ;
if ( symflag < 2 ) {
   fprintf(msgFile, "\n ubound = max(max(abs(u))) ") ;
} else {
   fprintf(msgFile, 
           "\n lbound = max(max(abs(l))) "
           "\n ubound = max(max(abs(u))) ") ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Chv_free(chv) ;
Drand_free(drand) ;
DV_free(workDV) ;
           
fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
