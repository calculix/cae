/*  test_copyEntriesToVector.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------
   test the copyEntriesToVector routine

   created -- 98may01, cca,
   ------------------------------------
*/
{
Chv      *chvJ, *chvI ;
double   imag, real, t1, t2 ;
double   *dvec, *entries ;
Drand    *drand ;
FILE     *msgFile ;
int      count, first, ierr, ii, iilast, ipivot, irow, jcol, jj, 
         jjlast, maxnent, mm, msglvl, ncol, nD, nent, nentD, nentL, 
         nentL11, nentL21, nentU, nentU11, nentU12, nL, npivot, nrow,
         nU, pivotingflag, seed, storeflag, symflag, total, type ;
int      *colind, *pivotsizes, *rowind ;

if ( argc != 10 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile nD nU type symflag "
"\n         pivotingflag storeflag seed"
"\n    msglvl    -- message level"
"\n    msgFile   -- message file"
"\n    nD        -- # of rows and columns in the (1,1) block"
"\n    nU        -- # of columns in the (1,2) block"
"\n    type      -- entries type"
"\n        1 --> real"
"\n        2 --> complex"
"\n    symflag   -- symmetry flag"
"\n        0 --> symmetric"
"\n        1 --> nonsymmetric"
"\n    pivotingflag -- pivoting flag"
"\n        if symflag = 1 and pivotingflag = 1 then"
"\n           construct pivotsizes[] vector"
"\n        endif"
"\n    storeflag -- flag to denote how to store entries"
"\n        0 --> store by rows"
"\n        1 --> store by columns"
"\n    seed      -- random number seed"
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
nD           = atoi(argv[3]) ;
nU           = atoi(argv[4]) ;
type         = atoi(argv[5]) ;
symflag      = atoi(argv[6]) ;
pivotingflag = atoi(argv[7]) ;
storeflag    = atoi(argv[8]) ;
seed         = atoi(argv[9]) ;
if ( msglvl > 0 ) {
   switch ( storeflag ) {
   case 0  : fprintf(msgFile, "\n\n %% STORE BY ROWS") ; break ;
   case 1  : fprintf(msgFile, "\n\n %% STORE BY COLUMNS") ; break ;
   default : 
      fprintf(stderr, "\n bad value %d for storeflag", storeflag) ;
      break ;
   }
}
nL = nU ;
if ( symflag == SPOOLES_NONSYMMETRIC ) {
   pivotingflag = 0 ;
}
/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
drand = Drand_new() ;
Drand_init(drand) ;
Drand_setNormal(drand, 0.0, 1.0) ;
Drand_setSeed(drand, seed) ;
/*
   --------------------------
   initialize the chvJ object
   --------------------------
*/
MARKTIME(t1) ;
chvJ = Chv_new() ;
Chv_init(chvJ, 0, nD, nL, nU, type, symflag) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix objects",
        t2 - t1) ;
nent = Chv_nent(chvJ) ;
entries = Chv_entries(chvJ) ;
if ( CHV_IS_REAL(chvJ) ) {
   Drand_fillDvector(drand, nent, entries) ;
} else if ( CHV_IS_COMPLEX(chvJ) ) {
   Drand_fillDvector(drand, 2*nent, entries) ;
}
Chv_columnIndices(chvJ, &ncol, &colind) ;
IVramp(ncol, colind, 0, 1) ;
if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
   Chv_rowIndices(chvJ, &nrow, &rowind) ;
   IVramp(nrow, rowind, 0, 1) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% chevron a") ;
   Chv_writeForMatlab(chvJ, "a", msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------
   initialize the chvI object
   --------------------------
*/
MARKTIME(t1) ;
chvI = Chv_new() ;
Chv_init(chvI, 0, nD, nL, nU, type, symflag) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix objects",
        t2 - t1) ;
Chv_zero(chvI) ;
Chv_columnIndices(chvI, &ncol, &colind) ;
IVramp(ncol, colind, 0, 1) ;
if ( CHV_IS_NONSYMMETRIC(chvI) ) {
   Chv_rowIndices(chvI, &nrow, &rowind) ;
   IVramp(nrow, rowind, 0, 1) ;
}
if ( symflag == 0 && pivotingflag == 1 ) {
/*
   ------------------------------
   create the pivotsizes[] vector
   ------------------------------
*/
   Drand_setUniform(drand, 1, 2.999) ;
   pivotsizes = IVinit(nD, 0) ;
   Drand_fillIvector(drand, nD, pivotsizes) ;
/*
   fprintf(msgFile, "\n initial pivotsizes[] : ") ;
   IVfp80(msgFile, nD, pivotsizes, 80, &ierr) ;
*/
   for ( npivot = count = 0 ; npivot < nD ; npivot++ ) {
      count += pivotsizes[npivot] ;
      if ( count > nD ) {
         pivotsizes[npivot]-- ;
         count-- ;
      } 
      if ( count == nD ) {
         break ;
      }
   }
   npivot++ ;
/*
   fprintf(msgFile, "\n final pivotsizes[] : ") ;
   IVfp80(msgFile, npivot, pivotsizes, 80, &ierr) ;
*/
} else {
   npivot = 0 ;
   pivotsizes = NULL ;
}
/*
   --------------------------------------------------
   first test: copy lower, diagonal and upper entries
   --------------------------------------------------
*/
if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
   nentL = Chv_countEntries(chvJ, npivot, pivotsizes, CHV_STRICT_LOWER);
} else {
   nentL = 0 ;
}
nentD = Chv_countEntries(chvJ, npivot, pivotsizes, CHV_DIAGONAL) ;
nentU = Chv_countEntries(chvJ, npivot, pivotsizes, CHV_STRICT_UPPER) ;
maxnent = nentL ;
if ( maxnent < nentD ) { maxnent = nentD ; }
if ( maxnent < nentU ) { maxnent = nentU ; }
if ( CHV_IS_REAL(chvJ) ) {
   dvec = DVinit(maxnent, 0.0) ;
} else if ( CHV_IS_COMPLEX(chvJ) ) {
   dvec = DVinit(2*maxnent, 0.0) ;
}
if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
/*
   --------------------------------------
   copy the entries in the lower triangle,
   then move into the chvI object
   --------------------------------------
*/
   nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                                  dvec, CHV_STRICT_LOWER, storeflag) ;
   if ( nent != nentL ) {
      fprintf(stderr, "\n error: nentL = %d, nent = %d", nentL, nent) ;
      exit(-1) ;
   }
   if ( storeflag == 0 ) {
      for ( irow = 0, mm = 0 ; irow < nrow ; irow++ ) {
         jjlast = (irow < nD) ? irow - 1 : nD - 1 ;
         for ( jj = 0 ; jj <= jjlast ; jj++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jj, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ;
               imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jj, real, imag) ;
            }
         }
      }
   } else {
      for ( jcol = 0, mm = 0 ; jcol < nD ; jcol++ ) {
         for ( irow = jcol + 1 ; irow < nrow ; irow++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
/*
fprintf(msgFile, "\n %% mm = %d, a(%d,%d) = %20.12e + %20.12e*i",
        mm, irow, jcol, real, imag) ;
*/
               Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
            }
         }
      }
   }
}
/*
   ---------------------------------------
   copy the entries in the diagonal matrix
   then move into the chvI object
   ---------------------------------------
*/
nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                               dvec, CHV_DIAGONAL, storeflag) ;
if ( nent != nentD ) {
   fprintf(stderr, "\n error: nentD = %d, nent = %d", nentD, nent) ;
   exit(-1) ;
}
if ( pivotsizes == NULL ) {
   for ( jcol = 0, mm = 0 ; jcol < nD ; jcol++, mm++ ) {
      if ( CHV_IS_REAL(chvJ) ) {
         real = dvec[mm] ; 
         Chv_setRealEntry(chvI, jcol, jcol, real) ;
      } else if ( CHV_IS_COMPLEX(chvJ) ) {
         real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
         Chv_setComplexEntry(chvI, jcol, jcol, real, imag) ;
      }
   }
} else {
   for ( ipivot = irow = mm = 0 ; ipivot < npivot ; ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         if ( CHV_IS_REAL(chvJ) ) {
            real = dvec[mm] ; 
            Chv_setRealEntry(chvI, irow, irow, real) ;
         } else if ( CHV_IS_COMPLEX(chvJ) ) {
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, irow, real, imag) ;
         }
         mm++ ; irow++ ;
      } else {
         if ( CHV_IS_REAL(chvJ) ) {
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow, irow, real) ;
            mm++ ; 
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow, irow+1, real) ;
            mm++ ; 
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow+1, irow+1, real) ;
            mm++ ; 
            irow += 2 ;
         } else if ( CHV_IS_COMPLEX(chvJ) ) {
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, irow, real, imag) ;
            mm++ ; 
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, irow+1, real, imag) ;
            mm++ ; 
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow+1, irow+1, real, imag) ;
            mm++ ; 
            irow += 2 ;
         }
      }
   }
}
/*
   --------------------------------------
   copy the entries in the upper triangle,
   then move into the chvI object
   --------------------------------------
*/
nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                               dvec, CHV_STRICT_UPPER, storeflag) ;
if ( nent != nentU ) {
   fprintf(stderr, "\n error: nentU = %d, nent = %d", nentU, nent) ;
   exit(-1) ;
}
if ( storeflag == 1 ) {
   if ( pivotsizes == NULL ) {
      for ( jcol = mm = 0 ; jcol < ncol ; jcol++ ) {
         iilast = (jcol < nD) ? jcol - 1 : nD - 1 ;
         for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ; 
               Chv_setRealEntry(chvI, ii, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, ii, jcol, real, imag) ;
            }
         }
      }
   } else {
      for ( ipivot = jcol = mm = 0 ; ipivot < npivot ; ipivot++ ) {
         iilast = jcol - 1 ;
         for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, ii, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, ii, jcol, real, imag) ;
            }
         }
         jcol++ ;
         if ( pivotsizes[ipivot] == 2 ) {
            for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, ii, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, ii, jcol, real, imag) ;
               }
            }
            jcol++ ;
         }
      }
      for ( jcol = nD ; jcol < ncol ; jcol++ ) {
         for ( irow = 0 ; irow < nD ; irow++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
            }
         }
      }
   }
} else {
   if ( pivotsizes == NULL ) {
      for ( irow = mm = 0 ; irow < nD ; irow++ ) {
         for ( jcol = irow + 1 ; jcol < ncol ; jcol++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
            }
         }
      }
   } else {
      for ( ipivot = irow = mm = 0 ; ipivot < npivot ; ipivot++ ) {
         if ( pivotsizes[ipivot] == 1 ) {
            for ( jcol = irow + 1 ; jcol < ncol ; jcol++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, irow, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
               }
            }
            irow++ ;
         } else {
            for ( jcol = irow + 2 ; jcol < ncol ; jcol++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, irow, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
               }
            }
            for ( jcol = irow + 2 ; jcol < ncol ; jcol++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, irow+1, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, irow+1, jcol, real, imag) ;
               }
            }
            irow += 2 ;
         }
      }
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% chevron b") ;
   Chv_writeForMatlab(chvI, "b", msgFile) ;
   fprintf(msgFile, 
           "\n\n emtx1 = abs(a - b) ; enorm1 = max(max(emtx1))") ;
   fflush(msgFile) ;
}
DVfree(dvec) ;
/*
   -----------------------------------------------------
   second test: copy lower (1,1), lower (2,1), diagonal,
                upper(1,1) and upper(1,2) blocks
   -----------------------------------------------------
*/
if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
   nentL11 = Chv_countEntries(chvJ, npivot, pivotsizes, 
                              CHV_STRICT_LOWER_11) ;
   nentL21 = Chv_countEntries(chvJ, npivot, pivotsizes, 
                              CHV_LOWER_21) ;
} else {
   nentL11 = 0 ;
   nentL21 = 0 ;
}
nentD   = Chv_countEntries(chvJ, npivot, pivotsizes, CHV_DIAGONAL) ;
nentU11 = Chv_countEntries(chvJ, npivot, pivotsizes, 
                           CHV_STRICT_UPPER_11) ;
nentU12 = Chv_countEntries(chvJ, npivot, pivotsizes, 
                           CHV_UPPER_12) ;
maxnent = nentL11 ;
if ( maxnent < nentL21 ) { maxnent = nentL21 ; }
if ( maxnent < nentD   ) { maxnent = nentD   ; }
if ( maxnent < nentU11 ) { maxnent = nentU11 ; }
if ( maxnent < nentU12 ) { maxnent = nentU12 ; }
fprintf(msgFile, 
        "\n %% nentL11 = %d, nentL21 = %d"
        "\n %% nentD = %d, nentU11 = %d, nentU12 = %d",
        nentL11, nentL21, nentD, nentU11, nentU12) ;
if ( CHV_IS_REAL(chvJ) ) {
   dvec = DVinit(maxnent, 0.0) ;
} else if ( CHV_IS_COMPLEX(chvJ) ) {
   dvec = DVinit(2*maxnent, 0.0) ;
}
Chv_zero(chvI) ;
if ( CHV_IS_NONSYMMETRIC(chvJ) ) {
/*
   ------------------------------------------
   copy the entries in the lower (1,1) block,
   then move into the chvI object
   ------------------------------------------
*/
   nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                                 dvec, CHV_STRICT_LOWER_11, storeflag) ;
   if ( nent != nentL11 ) {
      fprintf(stderr, "\n error: nentL = %d, nent = %d", nentL, nent) ;
      exit(-1) ;
   }
   if ( storeflag == 0 ) {
      for ( irow = 0, mm = 0 ; irow < nD ; irow++ ) {
         jjlast = (irow < nD) ? irow - 1 : nD - 1 ;
         for ( jj = 0 ; jj <= jjlast ; jj++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jj, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jj, real, imag) ;
            }
         }
      }
   } else {
      for ( jcol = 0, mm = 0 ; jcol < nD ; jcol++ ) {
         for ( irow = jcol + 1 ; irow < nD ; irow++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
            }
         }
      }
   }
/*
   ------------------------------------------
   copy the entries in the lower (2,1) block,
   then move into the chvI object
   ------------------------------------------
*/
   nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                                  dvec, CHV_LOWER_21, storeflag);
   if ( nent != nentL21 ) {
      fprintf(stderr, "\n error: nentL21 = %d, nent = %d", 
              nentL21, nent) ;
      exit(-1) ;
   }
   if ( storeflag == 0 ) {
      for ( irow = nD, mm = 0 ; irow < nrow ; irow++ ) {
         for ( jcol = 0 ; jcol < nD ; jcol++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
            }
         }
      }
   } else {
      for ( jcol = 0, mm = 0 ; jcol < nD ; jcol++ ) {
         for ( irow = nD ; irow < nrow ; irow++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
            }
         }
      }
   }
}
/*
   ---------------------------------------
   copy the entries in the diagonal matrix
   then move into the chvI object
   ---------------------------------------
*/
nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                               dvec, CHV_DIAGONAL, storeflag) ;
if ( nent != nentD ) {
   fprintf(stderr, "\n error: nentD = %d, nent = %d", nentD, nent) ;
   exit(-1) ;
}
if ( pivotsizes == NULL ) {
   for ( jcol = 0, mm = 0 ; jcol < nD ; jcol++, mm++ ) {
      if ( CHV_IS_REAL(chvJ) ) {
         real = dvec[mm] ;
         Chv_setRealEntry(chvI, jcol, jcol, real) ;
      } else if ( CHV_IS_COMPLEX(chvJ) ) {
         real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
         Chv_setComplexEntry(chvI, jcol, jcol, real, imag) ;
      }
   }
} else {
   for ( ipivot = irow = mm = 0 ; ipivot < npivot ; ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         if ( CHV_IS_REAL(chvJ) ) {
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow, irow, real) ;
         } else if ( CHV_IS_COMPLEX(chvJ) ) {
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, irow, real, imag) ;
         }
         mm++ ; irow++ ;
      } else {
         if ( CHV_IS_REAL(chvJ) ) {
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow, irow, real) ;
            mm++ ; 
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow, irow+1, real) ;
            mm++ ; 
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow+1, irow+1, real) ;
            mm++ ; 
         } else if ( CHV_IS_COMPLEX(chvJ) ) {
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, irow, real, imag) ;
            mm++ ; 
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, irow+1, real, imag) ;
            mm++ ; 
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow+1, irow+1, real, imag) ;
            mm++ ; 
         }
         irow += 2 ;
      }
   }
}
/*
   -----------------------------------------
   copy the entries in the upper (1,1) block
   then move into the chvI object
   -----------------------------------------
*/
nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                               dvec, CHV_STRICT_UPPER_11, storeflag) ;
if ( nent != nentU11 ) {
   fprintf(stderr, "\n error: nentU11 = %d, nent = %d", nentU11, nent) ;
   exit(-1) ;
}
if ( storeflag == 1 ) {
   if ( pivotsizes == NULL ) {
      for ( jcol = mm = 0 ; jcol < nD ; jcol++ ) {
         iilast = (jcol < nD) ? jcol - 1 : nD - 1 ;
         for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, ii, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, ii, jcol, real, imag) ;
            }
         }
      }
   } else {
      for ( ipivot = jcol = mm = 0 ; ipivot < npivot ; ipivot++ ) {
         iilast = jcol - 1 ;
         for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, ii, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, ii, jcol, real, imag) ;
            }
         }
         jcol++ ;
         if ( pivotsizes[ipivot] == 2 ) {
            for ( ii = 0 ; ii <= iilast ; ii++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, ii, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, ii, jcol, real, imag) ;
               }
            }
            jcol++ ;
         }
      }
   }
} else {
   if ( pivotsizes == NULL ) {
      for ( irow = mm = 0 ; irow < nD ; irow++ ) {
         for ( jcol = irow + 1 ; jcol < nD ; jcol++, mm++ ) {
            if ( CHV_IS_REAL(chvJ) ) {
               real = dvec[mm] ;
               Chv_setRealEntry(chvI, irow, jcol, real) ;
            } else if ( CHV_IS_COMPLEX(chvJ) ) {
               real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
               Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
            }
         }
      }
   } else {
      for ( ipivot = irow = mm = 0 ; ipivot < npivot ; ipivot++ ) {
         if ( pivotsizes[ipivot] == 1 ) {
            for ( jcol = irow + 1 ; jcol < nD ; jcol++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, irow, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
               }
            }
            irow++ ;
         } else {
            for ( jcol = irow + 2 ; jcol < nD ; jcol++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, irow, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
               }
            }
            for ( jcol = irow + 2 ; jcol < nD ; jcol++, mm++ ) {
               if ( CHV_IS_REAL(chvJ) ) {
                  real = dvec[mm] ;
                  Chv_setRealEntry(chvI, irow+1, jcol, real) ;
               } else if ( CHV_IS_COMPLEX(chvJ) ) {
                  real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
                  Chv_setComplexEntry(chvI, irow+1, jcol, real, imag) ;
               }
            }
            irow += 2 ;
         }
      }
   }
}
/*
   -----------------------------------------
   copy the entries in the upper (1,2) block
   then move into the chvI object
   -----------------------------------------
*/
nent = Chv_copyEntriesToVector(chvJ, npivot, pivotsizes, maxnent, 
                               dvec, CHV_UPPER_12, storeflag) ;
if ( nent != nentU12 ) {
   fprintf(stderr, "\n error: nentU12 = %d, nent = %d", nentU12, nent) ;
   exit(-1) ;
}
if ( storeflag == 1 ) {
   for ( jcol = nD, mm = 0 ; jcol < ncol ; jcol++ ) {
      for ( irow = 0 ; irow < nD ; irow++, mm++ ) {
         if ( CHV_IS_REAL(chvJ) ) {
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow, jcol, real) ;
         } else if ( CHV_IS_COMPLEX(chvJ) ) {
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
         }
      }
   }
} else {
   for ( irow = mm = 0 ; irow < nD ; irow++ ) {
      for ( jcol = nD ; jcol < ncol ; jcol++, mm++ ) {
         if ( CHV_IS_REAL(chvJ) ) {
            real = dvec[mm] ;
            Chv_setRealEntry(chvI, irow, jcol, real) ;
         } else if ( CHV_IS_COMPLEX(chvJ) ) {
            real = dvec[2*mm] ; imag = dvec[2*mm+1] ;
            Chv_setComplexEntry(chvI, irow, jcol, real, imag) ;
         }
      }
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% chevron b") ;
   Chv_writeForMatlab(chvI, "b", msgFile) ;
   fprintf(msgFile, 
           "\n\n emtx2 = abs(a - b) ; enorm2 = max(max(emtx2))") ;
   fprintf(msgFile, "\n\n [ enorm1 enorm2]") ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( pivotsizes != NULL ) {
   IVfree(pivotsizes) ;
}
Chv_free(chvJ) ;
Chv_free(chvI) ;
Drand_free(drand) ;
DVfree(dvec) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
