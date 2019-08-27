/*  test_addChevron.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------
   test the Chv_addChevron() method.

   created -- 98apr18, cca
   ---------------------------------------
*/
{
Chv     *chv ;
double   alpha[2] ;
double   imag, real, t1, t2 ;
double   *chvent, *entries ;
Drand    *drand ;
FILE     *msgFile ;
int      chvsize, count, ichv, ierr, ii, iloc, irow, jcol,
         lastcol, msglvl, ncol, nD, nent, nL, nrow, nU, 
         off, seed, symflag, type, upper ;
int      *chvind, *colind, *keys, *rowind, *temp ;

if ( argc != 10 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile nD nU type symflag seed "
           "\n         alphareal alphaimag"
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
           "\n    alpha   -- scaling parameter"
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
nD       = atoi(argv[3]) ;
nU       = atoi(argv[4]) ;
type     = atoi(argv[5]) ;
symflag  = atoi(argv[6]) ;
seed     = atoi(argv[7]) ;
alpha[0] = atof(argv[8]) ;
alpha[1] = atof(argv[9]) ;
if (  nD <= 0 || nU < 0 || symflag < 0 || symflag > 2 ) {
   fprintf(stderr, "\n invalid input"
           "\n nD = %d, nU = %d, symflag = %d\n", nD, nU, symflag) ;
   exit(-1) ;
}
fprintf(msgFile, "\n alpha = %12.4e + %12.4e*i ;", alpha[0], alpha[1]) ;
nL = nU ;
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
Drand_init(drand) ;
Drand_setSeed(drand, seed) ;
Drand_setUniform(drand, -1.0, 1.0) ;
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
temp = IVinit(2*(nD+nU), -1) ;
IVramp(2*(nD+nU), temp, 0, 1) ;
IVshuffle(2*(nD+nU), temp, ++seed) ;
IVcopy(ncol, colind, temp) ;
IVqsortUp(ncol, colind) ;
if ( CHV_IS_NONSYMMETRIC(chv) ) {
   Chv_rowIndices(chv, &nrow, &rowind) ;
   IVcopy(nrow, rowind, colind) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n %% column indices") ;
   IVfprintf(msgFile, ncol, colind) ;
}
lastcol = colind[ncol-1] ;
nent = Chv_nent(chv) ;
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

if ( msglvl > 1 ) {
   fprintf(msgFile, "\n a = zeros(%d,%d) ;", lastcol+1, lastcol+1) ;
   Chv_writeForMatlab(chv, "a", msgFile) ;
}
/*
   --------------------------------------------------
   fill a chevron with random numbers and indices
   that are a subset of a front's, as in the assembly
   of original matrix entries.
   --------------------------------------------------
*/
Drand_setUniform(drand, 0, nD) ;
iloc = (int) Drand_value(drand) ;
ichv = colind[iloc] ;
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   upper = nD - iloc + nU ;
} else {
   upper = 2*(nD - iloc) - 1 + nL + nU ;
}
Drand_setUniform(drand, 1, upper) ;
chvsize = (int) Drand_value(drand) ;
fprintf(msgFile, "\n %% iloc = %d, ichv = %d, chvsize = %d", 
        iloc, ichv, chvsize) ;
chvind  = IVinit(chvsize, -1) ;
chvent  = DVinit(2*chvsize, 0.0) ;
Drand_setNormal(drand, 0.0, 1.0) ;
if ( CHV_IS_REAL(chv) ) {
   Drand_fillDvector(drand, chvsize, chvent) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   Drand_fillDvector(drand, 2*chvsize, chvent) ;
}
keys    = IVinit(upper+1, -1) ;
keys[0] = 0 ;
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   for ( ii = iloc + 1, count = 1 ; ii < nD + nU ; ii++ ) {
      keys[count++] = colind[ii] - ichv ;
   }
} else {
   for ( ii = iloc + 1, count = 1 ; ii < nD + nU ; ii++ ) {
      keys[count++] =   colind[ii] - ichv ;
      keys[count++] = - colind[ii] + ichv ;
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% iloc = %d, ichv = %d", iloc, ichv) ;
   fprintf(msgFile, "\n %% upper = %d", upper) ;
   fprintf(msgFile, "\n %% chvsize = %d", chvsize) ;
   fprintf(msgFile, "\n %% initial keys") ;
   IVfprintf(msgFile, count, keys) ;
}
   IVshuffle(count, keys, ++seed) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% shuffled keys") ;
   IVfp80(msgFile, count, keys, 80, &ierr) ;
}
IVcopy(chvsize, chvind, keys) ;
if ( CHV_IS_REAL(chv) ) {
   IVDVqsortUp(chvsize, chvind, chvent) ;
} else if ( CHV_IS_COMPLEX(chv) ) {
   IVZVqsortUp(chvsize, chvind, chvent) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% chvind") ;
   IVfprintf(msgFile, chvsize, chvind) ;
}
if ( CHV_IS_HERMITIAN(chv) ) {
   for ( ii = 0 ; ii < chvsize ; ii++ ) {
      if ( chvind[ii] == 0 ) {
         chvent[2*ii+1] = 0.0 ;
      }
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n b = zeros(%d,%d) ;", lastcol+1, lastcol+1) ;
   if ( CHV_IS_REAL(chv) ) {
      if ( CHV_IS_SYMMETRIC(chv) ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            off = chvind[ii] ;
            fprintf(msgFile, "\n b(%d,%d) = %20.12e ;",
                    colind[iloc]+1, colind[iloc]+off+1, chvent[ii]) ;
            fprintf(msgFile, "\n b(%d,%d) = %20.12e ;",
                    colind[iloc]+off+1, colind[iloc]+1, chvent[ii]) ;
         }
      } else {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            off = chvind[ii] ;
            if ( off > 0 ) {
               fprintf(msgFile, "\n b(%d,%d) = %20.12e ;",
                       colind[iloc]+1, colind[iloc]+off+1, chvent[ii]) ;
            } else {
               fprintf(msgFile, "\n b(%d,%d) = %20.12e ;",
                       colind[iloc]-off+1, colind[iloc]+1, chvent[ii]) ;
            }
         }
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            off = chvind[ii] ;
            fprintf(msgFile, "\n b(%d,%d) = %20.12e + %20.12e*i;",
                    colind[iloc]+1, colind[iloc]+off+1,
                    chvent[2*ii], chvent[2*ii+1]) ;
            if ( CHV_IS_HERMITIAN(chv) ) {
               fprintf(msgFile, "\n b(%d,%d) = %20.12e + %20.12e*i;",
                       colind[iloc]+off+1, colind[iloc]+1, 
                       chvent[2*ii], -chvent[2*ii+1]) ;
            } else {
               fprintf(msgFile, "\n b(%d,%d) = %20.12e + %20.12e*i;",
                       colind[iloc]+off+1, colind[iloc]+1, 
                       chvent[2*ii], chvent[2*ii+1]) ;
            }
         }
      } else {
         for ( ii = 0 ; ii < chvsize ; ii++ ) {
            off = chvind[ii] ;
            if ( off > 0 ) {
               fprintf(msgFile, "\n b(%d,%d) = %20.12e + %20.12e*i;",
                       colind[iloc]+1, colind[iloc]+off+1,
                       chvent[2*ii], chvent[2*ii+1]) ;
            } else {
               fprintf(msgFile, "\n b(%d,%d) = %20.12e + %20.12e*i;",
                       colind[iloc]-off+1, colind[iloc]+1, 
                       chvent[2*ii], chvent[2*ii+1]) ;
            }
         }
      }
   }
}
/*
   ------------------------------------
   add the chevron into the Chv object
   ------------------------------------
*/
Chv_addChevron(chv, alpha, ichv, chvsize, chvind, chvent) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n %% after adding the chevron") ;
   fprintf(msgFile, "\n c = zeros(%d,%d) ;", lastcol+1, lastcol+1) ;
   Chv_writeForMatlab(chv, "c", msgFile) ;
}
/*
   -----------------
   compute the error
   -----------------
*/
fprintf(msgFile, "\n max(max(abs(c - (a + alpha*b))))") ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Chv_free(chv) ;
Drand_free(drand) ;
IVfree(temp) ;
IVfree(chvind) ;
DVfree(chvent) ;
IVfree(keys) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
