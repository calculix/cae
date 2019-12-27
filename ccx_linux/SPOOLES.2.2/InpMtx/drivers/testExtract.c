/*  testExtract.c  */

#include "../InpMtx.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------------------------
   generate a random matrix and test the submatrix extraction method.
   the output is a matlab file to test correctness.

   created -- 98oct15, cca
 --------------------------------------------------------------------
*/
{
int        coordType, dataType, ii, msglvl, ncolA, ncol1, ncol2, 
           nitem, nrowA, nrow1, nrow2, rc, seed, symmetryflag ;
int        *colids, *rowids, *temp ;
InpMtx     *A ;
IV         *cols1IV, *cols2IV, *rows1IV, *rows2IV ;
FILE       *msgFile ;

if ( argc != 14 ) {
   fprintf(stdout, 
   "\n\n %% usage : %s msglvl msgFile dataType symmetryflag coordType "
   "\n %%         nrowA ncolA nitem nrow1 nrow2 ncol1 ncol2 seed"
   "\n %%    msglvl   -- message level"
   "\n %%    msgFile  -- message file"
   "\n %%    dataType -- type of matrix entries"
   "\n %%       1 -- real"
   "\n %%       2 -- complex"
   "\n %%    symmetryflag -- type of matrix symmetry"
   "\n %%       0 -- symmetric"
   "\n %%       1 -- hermitian"
   "\n %%       2 -- nonsymmetric"
   "\n %%    coordType -- coordinate Type"
   "\n %%       1 -- by rows"
   "\n %%       2 -- by columns"
   "\n %%       3 -- by chevrons"
   "\n %%    nrowA    -- number of rows in A"
   "\n %%    ncolA    -- number of columns in A"
   "\n %%    nitem    -- number of items to be loaded into A"
   "\n %%    nrow1    -- number of rows in the first block"
   "\n %%    nrow2    -- number of rows in the second block"
   "\n %%    ncol1    -- number of columns in the first block"
   "\n %%    ncol2    -- number of columns in the second block"
   "\n %%    seed     -- random number seed"
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
dataType      = atoi(argv[3]) ;
symmetryflag  = atoi(argv[4]) ;
coordType     = atoi(argv[5]) ;
nrowA         = atoi(argv[6]) ;
ncolA         = atoi(argv[7]) ;
nitem         = atoi(argv[8]) ;
nrow1         = atoi(argv[9]) ;
nrow2         = atoi(argv[10]) ;
ncol1         = atoi(argv[11]) ;
ncol2         = atoi(argv[12]) ;
seed          = atoi(argv[13]) ;
fprintf(msgFile, 
        "\n %% %s "
        "\n %% msglvl        -- %d" 
        "\n %% msgFile       -- %s" 
        "\n %% dataType      -- %d" 
        "\n %% symmetryflag  -- %d" 
        "\n %% coordType     -- %d" 
        "\n %% nrowA         -- %d" 
        "\n %% ncolA         -- %d" 
        "\n %% nitem         -- %d" 
        "\n %% nrow1         -- %d" 
        "\n %% nrow2         -- %d" 
        "\n %% ncol1         -- %d" 
        "\n %% ncol2         -- %d" 
        "\n %% seed          -- %d"
        "\n",
        argv[0], msglvl, argv[2], dataType, symmetryflag, coordType,
        nrowA, ncolA, nitem, nrow1, nrow2, ncol1, ncol2, seed) ;
fflush(msgFile) ;
if ( dataType != SPOOLES_REAL && dataType != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n invalid value %d for dataType\n", dataType) ;
   exit(-1) ;
}
if (  symmetryflag != SPOOLES_SYMMETRIC
   && symmetryflag != SPOOLES_HERMITIAN
   && symmetryflag != SPOOLES_NONSYMMETRIC ) {
   fprintf(stderr, 
           "\n invalid value %d for symmetryflag\n", symmetryflag) ;
   exit(-1) ;
}
if ( symmetryflag == SPOOLES_HERMITIAN && dataType != SPOOLES_COMPLEX ){
   fprintf(stderr, 
           "\n symmetryflag is hermitian, data type is not complex\n") ;
   exit(-1) ;
}
if (  coordType != INPMTX_BY_ROWS 
   && coordType != INPMTX_BY_COLUMNS 
   && coordType != INPMTX_BY_CHEVRONS ) {
   fprintf(stderr, 
           "\n invalid value %d for coordType\n", coordType) ;
   exit(-1) ;
}
if ( nrowA <= 0 || ncolA <= 0 || nitem <= 0 ) {
   fprintf(stderr, 
           "\n invalid value: nrow = %d, ncol = %d, nitem = %d",
           nrowA, ncolA, nitem) ;
   exit(-1) ;
}
fprintf(msgFile, "\n %% symmetryflag %d, nrowA %d, ncolA %d\n",
        symmetryflag, nrowA, ncolA) ;
if (   (   symmetryflag == SPOOLES_SYMMETRIC
        || symmetryflag == SPOOLES_HERMITIAN)
    && nrowA != ncolA ) {
   fprintf(stderr, 
           "\n symmetric matrix, nrowA %d, ncolA %d\n", nrowA, ncolA) ;
   exit(-1) ;
}
if ( nrow1 < 0 || nrow2 < 0 || (nrow1 + nrow2 != nrowA) ) {
   fprintf(stderr, 
           "\n invalid value: nrow = %d, nrow1 = %d, nrow2 = %d",
           nrowA, nrow1, nrow2) ;
   exit(-1) ;
}
if ( ncol1 < 0 || ncol2 < 0 || (ncol1 + ncol2 != ncolA) ) {
   fprintf(stderr, 
           "\n invalid value: ncol = %d, ncol1 = %d, ncol2 = %d",
           ncolA, ncol1, ncol2) ;
   exit(-1) ;
}
/*
   ----------------------------
   initialize the matrix object
   ----------------------------
*/
A = InpMtx_new() ;
rc = InpMtx_randomMatrix(A, dataType, coordType, INPMTX_BY_VECTORS,
                   nrowA, ncolA, symmetryflag, 0, nitem, seed) ;
/*
   -------------------------------------------
   write the assembled matrix to a matlab file
   -------------------------------------------
*/
fprintf(msgFile, "\n A = zeros(%d,%d) ;", nrowA, ncolA) ;
InpMtx_writeForMatlab(A, "A", msgFile) ;
if ( symmetryflag == SPOOLES_SYMMETRIC ) {
   fprintf(msgFile, 
          "\n A = diag(diag(A)) + triu(A,1) + transpose(triu(A,1)) ;") ;
} else if ( symmetryflag == SPOOLES_HERMITIAN ) {
   fprintf(msgFile, 
         "\n A = diag(diag(A)) + triu(A,1) + ctranspose(triu(A,1)) ;") ;
} 
/*
   ----------------------------------
   generate row and column id vectors
   ----------------------------------
*/
temp = IVinit(nrowA, -1) ;
IVramp(nrowA, temp, 0, 1) ;
IVshuffle(nrowA, temp, ++seed) ;
if ( nrow1 > 0 ) {
   rows1IV = IV_new() ;
   IV_init(rows1IV, nrow1, NULL) ;
   IVqsortUp(nrow1, temp) ;
   IVcopy(nrow1, IV_entries(rows1IV), temp) ;
} else {
   rows1IV = NULL ;
}
if ( nrow2 > 0 ) {
   rows2IV = IV_new() ;
   IV_init(rows2IV, nrow2, NULL) ;
   IVqsortUp(nrow2, temp + nrow1) ;
   IVcopy(nrow2, IV_entries(rows2IV), temp + nrow1) ;
} else {
   rows2IV = NULL ;
}
IVfree(temp) ;
if ( symmetryflag == SPOOLES_NONSYMMETRIC ) {
   temp = IVinit(ncolA, -1) ;
   IVramp(ncolA, temp, 0, 1) ;
   IVshuffle(ncolA, temp, ++seed) ;
   if ( ncol1 > 0 ) {
      cols1IV = IV_new() ;
      IV_init(cols1IV, ncol1, NULL) ;
      IVqsortUp(ncol1, temp) ;
      IVcopy(ncol1, IV_entries(cols1IV), temp) ;
   } else {
      cols1IV = NULL ;
   }
   if ( ncol2 > 0 ) {
      cols2IV = IV_new() ;
      IV_init(cols2IV, ncol2, NULL) ;
      IVqsortUp(ncol2, temp + ncol1) ;
      IVcopy(ncol2, IV_entries(cols2IV), temp + ncol1) ;
   } else {
      cols2IV = NULL ;
   }
   IVfree(temp) ;
} else {
   if ( ncol1 > 0 ) {
      cols1IV = IV_new() ;
      IV_init(cols1IV, ncol1, NULL) ;
      IV_copy(cols1IV, rows1IV) ;
   } else {
      cols1IV = NULL ;
   }
   if ( ncol2 > 0 ) {
      cols2IV = IV_new() ;
      IV_init(cols2IV, ncol2, NULL) ;
      IV_copy(cols2IV, rows2IV) ;
   } else {
      cols2IV = NULL ;
   }
}
if ( nrow1 > 0 ) {
   fprintf(msgFile, "\n rows1 = zeros(%d,1) ;", nrow1) ;
   IV_writeForMatlab(rows1IV, "rows1", msgFile) ;
}
if ( nrow2 > 0 ) {
   fprintf(msgFile, "\n rows2 = zeros(%d,1) ;", nrow2) ;
   IV_writeForMatlab(rows2IV, "rows2", msgFile) ;
}
if ( ncol1 > 0 ) {
   fprintf(msgFile, "\n cols1 = zeros(%d,1) ;", ncol1) ;
   IV_writeForMatlab(cols1IV, "cols1", msgFile) ;
}
if ( ncol2 > 0 ) {
   fprintf(msgFile, "\n cols2 = zeros(%d,1) ;", ncol2) ;
   IV_writeForMatlab(cols2IV, "cols2", msgFile) ;
}
/*
   -----------------------
   extract the submatrices
   -----------------------
*/
if ( nrow1 > 0 ) {
   if ( ncol1 > 0 ) {
      InpMtx   *A11 = InpMtx_new() ;
      rc = InpMtx_initFromSubmatrix(A11, A, rows1IV, cols1IV,
                                   symmetryflag, msglvl, msgFile) ;
      if ( rc != 1 ) {
         fprintf(stderr, "\n error return %d for A11\n", rc) ;
         exit(-1) ;
      }
      fprintf(msgFile, "\n A11 = zeros(%d,%d) ;", nrow1, ncol1) ;
      InpMtx_writeForMatlab(A11, "A11", msgFile) ;
      if ( symmetryflag == SPOOLES_SYMMETRIC ) {
         fprintf(msgFile, 
  "\n A11 = diag(diag(A11)) + triu(A11,1) + transpose(triu(A11,1)) ;") ;
      } else if ( symmetryflag == SPOOLES_HERMITIAN ) {
         fprintf(msgFile, 
 "\n A11 = diag(diag(A11)) + triu(A11,1) + ctranspose(triu(A11,1)) ;") ;
      } 
      fprintf(msgFile, 
              "\n err11 = max(max(abs(A11 - A(rows1,cols1)))) ;") ;
      InpMtx_free(A11) ;
   } else {
      fprintf(msgFile, "\n err11 = 0 ;") ;
   }
   if ( ncol2 > 0 ) {
      InpMtx   *A12 = InpMtx_new() ;
      rc = InpMtx_initFromSubmatrix(A12, A, rows1IV, cols2IV,
                                   symmetryflag, msglvl, msgFile) ;
      if ( rc != 1 ) {
         fprintf(stderr, "\n error return %d for A12\n", rc) ;
         exit(-1) ;
      }
      fprintf(msgFile, "\n A12 = zeros(%d,%d) ;", nrow1, ncol2) ;
      InpMtx_writeForMatlab(A12, "A12", msgFile) ;
      fprintf(msgFile, 
              "\n err12 = max(max(abs(A12 - A(rows1,cols2)))) ;") ;
      InpMtx_free(A12) ;
   } else {
      fprintf(msgFile, "\n err12 = 0 ;") ;
   }
} else {
   fprintf(msgFile, "\n err11 = 0 ;") ;
   fprintf(msgFile, "\n err12 = 0 ;") ;
}
if ( nrow2 > 0 ) {
   if ( ncol1 > 0 ) {
      InpMtx   *A21 = InpMtx_new() ;
      rc = InpMtx_initFromSubmatrix(A21, A, rows2IV, cols1IV,
                                   symmetryflag, msglvl, msgFile) ;
      if ( rc != 1 ) {
         fprintf(stderr, "\n error return %d for A21\n", rc) ;
         exit(-1) ;
      }
      fprintf(msgFile, "\n A21 = zeros(%d,%d) ;", nrow2, ncol1) ;
      InpMtx_writeForMatlab(A21, "A21", msgFile) ;
      fprintf(msgFile, 
              "\n err21 = max(max(abs(A21 - A(rows2,cols1)))) ;") ;
      InpMtx_free(A21) ;
   } else {
      fprintf(msgFile, "\n err21 = 0 ;") ;
   }
   if ( ncol2 > 0 ) {
      InpMtx   *A22 = InpMtx_new() ;
      rc = InpMtx_initFromSubmatrix(A22, A, rows2IV, cols2IV,
                                   symmetryflag, msglvl, msgFile) ;
      if ( rc != 1 ) {
         fprintf(stderr, "\n error return %d for A22\n", rc) ;
         exit(-1) ;
      }
      fprintf(msgFile, "\n A22 = zeros(%d,%d) ;", nrow2, ncol2) ;
      InpMtx_writeForMatlab(A22, "A22", msgFile) ;
      if ( symmetryflag == SPOOLES_SYMMETRIC ) {
         fprintf(msgFile, 
  "\n A22 = diag(diag(A22)) + triu(A22,1) + transpose(triu(A22,1)) ;") ;
      } else if ( symmetryflag == SPOOLES_HERMITIAN ) {
         fprintf(msgFile, 
 "\n A22 = diag(diag(A22)) + triu(A22,1) + ctranspose(triu(A22,1)) ;") ;
      } 
      fprintf(msgFile, 
              "\n err22 = max(max(abs(A22 - A(rows2,cols2)))) ;") ;
      InpMtx_free(A22) ;
   } else {
      fprintf(msgFile, "\n err22 = 0 ;") ;
   }
} else {
   fprintf(msgFile, "\n err21 = 0 ;") ;
   fprintf(msgFile, "\n err22 = 0 ;") ;
}
fprintf(msgFile, "\n error = [ err11 err12 ; err21 err22 ] ") ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
InpMtx_free(A) ;
if ( rows1IV != NULL ) { IV_free(rows1IV) ; }
if ( rows2IV != NULL ) { IV_free(rows2IV) ; }
if ( cols1IV != NULL ) { IV_free(cols1IV) ; }
if ( cols2IV != NULL ) { IV_free(cols2IV) ; }

fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
