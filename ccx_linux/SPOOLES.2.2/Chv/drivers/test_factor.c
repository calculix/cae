/*  test_factor.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------
   test the Chv_factor() method.
   the program's output is a matlab file
   to check correctness of the code.

   created -- 98apr22, cca
   -------------------------------------
*/
{
Chv     *chv ;
double   imag, real, tau, t1, t2 ;
double   *entries ;
Drand    *drand ;
DV       *workDV ;
FILE     *msgFile ;
int      ii, ipivot, irow, jcol, msglvl, ncol, nD, ndelay, 
         nelim, nent, nL, nrow, npivot, ntest, nU, pivotflag, 
         rc, seed, symflag, tag, type ;
int      *colind, *pivotsizes, *rowind ;
IV       *pivotsizesIV ;

if ( argc != 10 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile nD nU type symflag pivotflag seed tau "
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
"\n       2 --> nonsymmetric"
"\n    pivotflag -- pivoting flag"
"\n       0 --> no pivoting"
"\n       1 --> pivoting"
"\n    tau  -- bound on magnitude of factor entries"
"\n    seed -- random number seed"
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
nD        = atoi(argv[3]) ;
nU        = atoi(argv[4]) ;
type      = atoi(argv[5]) ;
symflag   = atoi(argv[6]) ;
pivotflag = atoi(argv[7]) ;
seed      = atoi(argv[8]) ;
tau       = atof(argv[9]) ;
fprintf(msgFile, "\n %% testChv:"
        "\n %% msglvl    = %d"
        "\n %% msgFile   = %s"
        "\n %% nD        = %d"
        "\n %% nU        = %d"
        "\n %% type      = %d"
        "\n %% symflag   = %d"
        "\n %% pivotflag = %d"
        "\n %% seed      = %d"
        "\n %% tau       = %f",
        msglvl, argv[2], nD, nU, type, symflag, pivotflag, seed, tau) ;
nL = nU ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if (  nD <= 0 || nL < 0 || nU < 0 
   || symflag < 0 || symflag > 2 ) {
   fprintf(stderr, "\n invalid input"
      "\n nD = %d, nL = %d, nU = %d, symflag = %d\n",
           nD, nL, nU, symflag) ;
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
Chv_writeForMatlab(chv, "A", msgFile) ;
if ( pivotflag == 1 ) {
   pivotsizesIV = IV_new() ;
} else {
   pivotsizesIV = NULL ;
}
workDV = DV_new() ;
/*
   -----------------
   factor the matrix
   -----------------
*/
ndelay = ntest = 0 ;
if ( pivotflag == SPOOLES_PIVOTING ) {
   nelim = Chv_factorWithPivoting(chv, ndelay, pivotflag, pivotsizesIV, 
                                  workDV, tau, &ntest) ;
} else {
   nelim = Chv_factorWithNoPivoting(chv, NULL) ;
}
fprintf(msgFile, "\n nD = %d ;\n nelim = %d", nD, nelim) ;
/*
   ---------------------
   write out the factors
   ---------------------
*/
Chv_rowIndices(chv, &nrow, &rowind) ;
Chv_columnIndices(chv, &ncol, &colind) ;
fprintf(msgFile, 
        "\n\n L = eye(%d,%d); "
        "\n D = zeros(%d,%d); "
        "\n T = zeros(%d,%d); "
        "\n U = eye(%d,%d); ", 
        ncol, ncol, ncol, ncol, ncol, ncol, ncol, ncol) ;
if ( pivotflag == 0 ) {
   if ( CHV_IS_REAL(chv) ) {
      for ( irow = 0 ; irow < nD ; irow++ ) {
         Chv_realEntry(chv, irow, irow, &real) ;
         fprintf(msgFile, "\n D(%d,%d) = %20.12e ;",
                 rowind[irow]+1, colind[irow]+1, real) ;
         for ( jcol = irow + 1 ; jcol < nD + nU ; jcol++ ) {
            Chv_realEntry(chv, irow, jcol, &real) ;
            fprintf(msgFile, "\n U(%d,%d) = %20.12e ;",
                    rowind[irow]+1, colind[jcol]+1, real) ;
            Chv_realEntry(chv, jcol, irow, &real) ;
            fprintf(msgFile, "\n L(%d,%d) = %20.12e ;",
                    rowind[jcol]+1, colind[irow]+1, real) ;
         }
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( irow = 0 ; irow < nD ; irow++ ) {
         Chv_complexEntry(chv, irow, irow, &real, &imag) ;
         fprintf(msgFile, "\n D(%d,%d) = %20.12e + %20.12e*i ;",
                 rowind[irow]+1, colind[irow]+1, real, imag) ;
         for ( jcol = irow + 1 ; jcol < nD + nU ; jcol++ ) {
            Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
            fprintf(msgFile, "\n U(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[irow]+1, colind[jcol]+1, real, imag) ;
            Chv_complexEntry(chv, jcol, irow, &real, &imag) ;
            fprintf(msgFile, "\n L(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[jcol]+1, colind[irow]+1, real, imag) ;
         }
      }
   }
} else {
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      fprintf(msgFile, "\n colind(%d) = %d;", 
              irow + 1, 1+colind[irow]) ;
   }
   if ( CHV_IS_NONSYMMETRIC(chv) ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         fprintf(msgFile, "\n rowind(%d) = %d;", 
                 irow + 1, 1+rowind[irow]) ;
      }
      IV_setSize(pivotsizesIV, nelim) ;
      IV_fill(pivotsizesIV, 1) ;
   } else {
      fprintf(msgFile, "\n rowind = colind ;") ;
   }
   fprintf(msgFile, "\n A = A(rowind,colind) ;") ;
   IVramp(nrow, rowind, 0, 1) ;
   IVramp(ncol, colind, 0, 1) ;
   IV_sizeAndEntries(pivotsizesIV, &npivot, &pivotsizes) ;
   fprintf(msgFile, "\n npivot = %d ;", npivot) ;
   if ( CHV_IS_REAL(chv) ) {
      for ( ipivot = irow = 0 ; ipivot < npivot ; ipivot++ ) {
         if ( pivotsizes[ipivot] == 1 ) {
            Chv_realEntry(chv, irow, irow, &real) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e ;",
                    rowind[irow]+1, colind[irow]+1, real) ;
            for ( jcol = irow + 1 ; jcol < nD + nU ; jcol++ ) {
               Chv_realEntry(chv, irow, jcol, &real) ;
               fprintf(msgFile, "\n U(%d,%d) = %20.12e ;",
                       rowind[irow]+1, colind[jcol]+1, real) ;
               Chv_realEntry(chv, jcol, irow, &real) ;
               fprintf(msgFile, "\n L(%d,%d) = %20.12e ;",
                       rowind[jcol]+1, colind[irow]+1, real) ;
            }
            irow += 1 ;
         } else {
            Chv_realEntry(chv, irow, irow, &real) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e ;",
                    rowind[irow]+1, colind[irow]+1, real) ;
            Chv_realEntry(chv, irow, irow+1, &real) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e ;",
                    rowind[irow]+1, colind[irow+1]+1, real) ;
            Chv_realEntry(chv, irow+1, irow, &real) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e ;",
                    rowind[irow+1]+1, colind[irow]+1, real) ;
            Chv_realEntry(chv, irow+1, irow+1, &real) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e ;",
                    rowind[irow+1]+1, colind[irow+1]+1, real) ;
            for ( jcol = irow + 2 ; jcol < nD + nU ; jcol++ ) {
               Chv_realEntry(chv, irow, jcol, &real) ;
               fprintf(msgFile, "\n U(%d,%d) = %20.12e ;",
                       rowind[irow]+1, colind[jcol]+1, real) ;
               Chv_realEntry(chv, jcol, irow, &real) ;
               fprintf(msgFile, "\n L(%d,%d) = %20.12e ;",
                       rowind[jcol]+1, colind[irow]+1, real) ;
               Chv_realEntry(chv, irow+1, jcol, &real) ;
               fprintf(msgFile, "\n U(%d,%d) = %20.12e ;",
                       rowind[irow+1]+1, colind[jcol]+1, real) ;
               Chv_realEntry(chv, jcol, irow+1, &real) ;
               fprintf(msgFile, "\n L(%d,%d) = %20.12e ;",
                       rowind[jcol]+1, colind[irow+1]+1, real) ;
            }
            irow += 2 ;
         }
      }
      for ( irow = nelim ; irow < nD ; irow++ ) {
         Chv_realEntry(chv, irow, irow, &real) ;
         fprintf(msgFile, "\n T(%d,%d) = %20.12e ;",
                 rowind[irow]+1, colind[irow]+1, real) ;
         for ( jcol = irow + 1 ; jcol < ncol ; jcol++ ) {
            Chv_realEntry(chv, irow, jcol, &real) ;
            fprintf(msgFile, "\n T(%d,%d) = %20.12e ;",
                    rowind[irow]+1, colind[jcol]+1, real) ;
            Chv_realEntry(chv, jcol, irow, &real) ;
            fprintf(msgFile, "\n T(%d,%d) = %20.12e ;",
                    rowind[jcol]+1, colind[irow]+1, real) ;
         }
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      for ( ipivot = irow = 0 ; ipivot < npivot ; ipivot++ ) {
         if ( pivotsizes[ipivot] == 1 ) {
            Chv_complexEntry(chv, irow, irow, &real, &imag) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[irow]+1, colind[irow]+1, real, imag) ;
            for ( jcol = irow + 1 ; jcol < nD + nU ; jcol++ ) {
               Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
               fprintf(msgFile, "\n U(%d,%d) = %20.12e + %20.12e*i ;",
                       rowind[irow]+1, colind[jcol]+1, real, imag) ;
               Chv_complexEntry(chv, jcol, irow, &real, &imag) ;
               fprintf(msgFile, "\n L(%d,%d) = %20.12e + %20.12e*i ;",
                       rowind[jcol]+1, colind[irow]+1, real, imag) ;
            }
            irow += 1 ;
         } else {
            Chv_complexEntry(chv, irow, irow, &real, &imag) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[irow]+1, colind[irow]+1, real, imag) ;
            Chv_complexEntry(chv, irow, irow+1, &real, &imag) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[irow]+1, colind[irow+1]+1, real, imag) ;
            Chv_complexEntry(chv, irow+1, irow, &real, &imag) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[irow+1]+1, colind[irow]+1, real, imag) ;
            Chv_complexEntry(chv, irow+1, irow+1, &real, &imag) ;
            fprintf(msgFile, "\n D(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[irow+1]+1, colind[irow+1]+1, real, imag) ;
            for ( jcol = irow + 2 ; jcol < nD + nU ; jcol++ ) {
               Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
               fprintf(msgFile, "\n U(%d,%d) = %20.12e + %20.12e*i ;",
                       rowind[irow]+1, colind[jcol]+1, real, imag) ;
               Chv_complexEntry(chv, jcol, irow, &real, &imag) ;
               fprintf(msgFile, "\n L(%d,%d) = %20.12e + %20.12e*i ;",
                       rowind[jcol]+1, colind[irow]+1, real, imag) ;
               Chv_complexEntry(chv, irow+1, jcol, &real, &imag) ;
               fprintf(msgFile, "\n U(%d,%d) = %20.12e + %20.12e*i ;",
                       rowind[irow+1]+1, colind[jcol]+1, real, imag) ;
               Chv_complexEntry(chv, jcol, irow+1, &real, &imag) ;
               fprintf(msgFile, "\n L(%d,%d) = %20.12e + %20.12e*i ;",
                       rowind[jcol]+1, colind[irow+1]+1, real, imag) ;
            }
            irow += 2 ;
         }
      }
      for ( irow = nelim ; irow < nD ; irow++ ) {
         Chv_complexEntry(chv, irow, irow, &real, &imag) ;
         fprintf(msgFile, "\n T(%d,%d) = %20.12e + %20.12e*i ;",
                 rowind[irow]+1, colind[irow]+1, real, imag) ;
         for ( jcol = irow + 1 ; jcol < ncol ; jcol++ ) {
            Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
            fprintf(msgFile, "\n T(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[irow]+1, colind[jcol]+1, real, imag) ;
            Chv_complexEntry(chv, jcol, irow, &real, &imag) ;
            fprintf(msgFile, "\n T(%d,%d) = %20.12e + %20.12e*i ;",
                    rowind[jcol]+1, colind[irow]+1, real, imag) ;
         }
      }
   }
}
fprintf(msgFile, "\n B = A ;") ;
fprintf(msgFile, 
        "\n B = A - T - L(:,1:%d) * D(1:%d,1:%d) * U(1:%d,:) ; ",
        nelim, nelim, nelim, nelim) ;
fprintf(msgFile, "\n B(%d:%d,%d:%d) = 0.0 ; ", nD+1, ncol, nD+1, ncol) ;
fprintf(msgFile, 
        "\n maxabsB = max(max(abs(B)))"
        "\n maxabsL = max(max(abs(L - eye(%d,%d))))"
        "\n maxabsU = max(max(abs(U - eye(%d,%d))))"
        "\n [ maxabsB maxabsL maxabsU ]",
        ncol, ncol, ncol, ncol) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Chv_free(chv) ;
DV_free(workDV) ;
if ( pivotsizesIV != NULL ) {
   IV_free(pivotsizesIV) ;
}
Drand_free(drand) ;
           
fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
