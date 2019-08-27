/*  test_factorupd.c  */

#include "../Chv.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------
   test the Chv_update{H,S,N}() methods.
   T := T - U^T * D * U
   T := T - U^H * D * U
   T := T - L   * D * U

   created -- 98apr23, cca
   -------------------------------------
*/
{
Chv     *chvT ;
SubMtx     *mtxD, *mtxL, *mtxU ;
double   imag, ops, real, t1, t2 ;
Drand    *drand ;
DV       *tempDV ;
FILE     *msgFile ;
int      irow, msglvl, ncolT, nDT, ncolU, nentT, nentU, nrowD, 
         nrowL, nrowT, offset, seed, size, sparsityflag, symflag, type ;
int      *colindT, *colindU, *ivec, *rowindL, *rowindT ;

if ( argc != 13 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile type symflag sparsityflag"
           "\n         ncolT ncolU nrowD nentU offset seed"
           "\n    msglvl  -- message level"
           "\n    msgFile -- message file"
           "\n    type    -- entries type"
           "\n       1 -- real"
           "\n       2 -- complex"
           "\n    symflag -- type of matrix U"
           "\n       0 -- symmetric"
           "\n       1 -- hermitian"
           "\n       2 -- nonsymmetric"
           "\n    sparsityflag -- dense or sparse"
           "\n       0 -- dense"
           "\n       1 -- sparse"
           "\n    ncolT   -- # of rows and columns in matrix T"
           "\n    nDT     -- # of internal rows and columns in matrix T"
           "\n    ncolU   -- # of rows and columns in matrix U"
           "\n    nrowD   -- # of rows and columns in matrix D"
           "\n    nentU   -- # of entries in matrix U"
           "\n    offset  -- distance between D_I and T"
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
type         = atoi(argv[3]) ;
symflag      = atoi(argv[4]) ;
sparsityflag = atoi(argv[5]) ;
ncolT        = atoi(argv[6]) ;
nDT          = atoi(argv[7]) ;
ncolU        = atoi(argv[8]) ;
nrowD        = atoi(argv[9]) ;
nentU        = atoi(argv[10]) ;
offset       = atoi(argv[11]) ;
seed         = atoi(argv[12]) ;
fprintf(msgFile, "\n %% %s:"
        "\n %% msglvl       = %d"
        "\n %% msgFile      = %s"
        "\n %% type         = %d"
        "\n %% symflag      = %d"
        "\n %% sparsityflag = %d"
        "\n %% ncolT        = %d"
        "\n %% nDT          = %d"
        "\n %% ncolU        = %d"
        "\n %% nrowD        = %d"
        "\n %% nentU        = %d"
        "\n %% offset       = %d"
        "\n %% seed         = %d",
        argv[0], msglvl, argv[2], type, symflag, sparsityflag, 
        ncolT, nDT, ncolU, nrowD, nentU, offset, seed) ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if (  (type != SPOOLES_REAL 
       && type != SPOOLES_COMPLEX) 
   || (symflag != SPOOLES_SYMMETRIC 
       && symflag != SPOOLES_HERMITIAN 
       && symflag != SPOOLES_NONSYMMETRIC) 
   || (sparsityflag < 0 || sparsityflag > 1)
   || ncolT <= 0 || ncolU > (ncolT + offset) || nrowD <= 0 ) {
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
Drand_setSeed(drand, ++seed) ;
Drand_setNormal(drand, 0.0, 1.0) ;
/*
   -----------------------
   get a vector of indices
   -----------------------
*/
size = nrowD + offset + ncolT ;
ivec = IVinit(size, -1) ;
IVramp(size, ivec, 0, 1) ;
/*
   ----------------------------
   initialize the T Chv object
   ----------------------------
*/
fprintf(msgFile, "\n\n %% symflag = %d", symflag) ;
MARKTIME(t1) ;
chvT = Chv_new() ;
Chv_init(chvT, 0, nDT, ncolT - nDT, ncolT - nDT, type, symflag) ;
nentT = Chv_nent(chvT) ;
if ( CHV_IS_REAL(chvT) ) {
   Drand_fillDvector(drand, nentT, Chv_entries(chvT)) ;
} else if ( CHV_IS_COMPLEX(chvT) ) {
   Drand_fillDvector(drand, 2*nentT, Chv_entries(chvT)) ;
}
Chv_columnIndices(chvT, &ncolT, &colindT) ;
IVcopy(ncolT, colindT, ivec + nrowD + offset) ;
if ( CHV_IS_NONSYMMETRIC(chvT) ) {
   Chv_rowIndices(chvT, &nrowT, &rowindT) ;
   IVcopy(nrowT, rowindT, colindT) ;
}
IVfree(ivec) ;
if ( CHV_IS_HERMITIAN(chvT) ) {
   fprintf(msgFile, "\n\n %% hermitian\n") ;
/*
   ---------------------------------------------------------
   hermitian example, set imaginary part of diagonal to zero
   ---------------------------------------------------------
*/
   for ( irow = 0 ; irow < nDT ; irow++ ) {
      Chv_complexEntry(chvT, irow, irow, &real, &imag) ;
      Chv_setComplexEntry(chvT, irow, irow, real, 0.0) ;
   }
}
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize chvT Chv object",
        t2 - t1) ;
fprintf(msgFile, "\n T = zeros(%d,%d); ", size, size) ;
Chv_writeForMatlab(chvT, "T", msgFile) ;
/*
   ---------------------------
   initialize the D Mtx object
   ---------------------------
*/
MARKTIME(t1) ;
mtxD = SubMtx_new() ;
if ( CHV_IS_REAL(chvT) ) {
   if ( CHV_IS_SYMMETRIC(chvT) ) {
      SubMtx_initRandom(mtxD, SPOOLES_REAL, SUBMTX_BLOCK_DIAGONAL_SYM,
                      0, 0, nrowD, nrowD, nrowD*nrowD, ++seed) ;
   } else {
      SubMtx_initRandom(mtxD, SPOOLES_REAL, SUBMTX_DIAGONAL, 
                      0, 0, nrowD, nrowD, nrowD*nrowD, ++seed) ;
   }
} else if ( CHV_IS_COMPLEX(chvT) ) {
   if ( CHV_IS_HERMITIAN(chvT) ) {
      SubMtx_initRandom(mtxD,SPOOLES_COMPLEX,SUBMTX_BLOCK_DIAGONAL_HERM,
                      0, 0, nrowD, nrowD, nrowD*nrowD, ++seed) ;
   } else if ( CHV_IS_SYMMETRIC(chvT) ) {
      SubMtx_initRandom(mtxD,SPOOLES_COMPLEX, SUBMTX_BLOCK_DIAGONAL_SYM,
                      0, 0, nrowD, nrowD, nrowD*nrowD, ++seed) ;
   } else {
      SubMtx_initRandom(mtxD, SPOOLES_COMPLEX, SUBMTX_DIAGONAL, 
                      0, 0, nrowD, nrowD, nrowD*nrowD, ++seed) ;
   }
}
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize D SubMtx object",
        t2 - t1) ;
fprintf(msgFile, "\n D = zeros(%d,%d) ;", nrowD, nrowD) ;
SubMtx_writeForMatlab(mtxD, "D", msgFile) ;
/*
   ----------------------------
   initialize the U SubMtx object
   ----------------------------
*/
MARKTIME(t1) ;
mtxU = SubMtx_new() ;
if ( CHV_IS_REAL(chvT) ) {
   if ( sparsityflag == 0 ) {
      SubMtx_initRandom(mtxU, SPOOLES_REAL, SUBMTX_DENSE_COLUMNS, 
                      0, 0, nrowD, ncolU, nentU, ++seed) ;
   } else {
      SubMtx_initRandom(mtxU, SPOOLES_REAL, SUBMTX_SPARSE_COLUMNS, 
                      0, 0, nrowD, ncolU, nentU, ++seed) ;
   }
} else if ( CHV_IS_COMPLEX(chvT) ) {
   if ( sparsityflag == 0 ) {
      SubMtx_initRandom(mtxU, SPOOLES_COMPLEX, SUBMTX_DENSE_COLUMNS, 
                      0, 0, nrowD, ncolU, nentU, ++seed) ;
   } else {
      SubMtx_initRandom(mtxU, SPOOLES_COMPLEX, SUBMTX_SPARSE_COLUMNS, 
                      0, 0, nrowD, ncolU, nentU, ++seed) ;
   }
}
ivec = IVinit(offset + ncolT, -1) ;
IVramp(offset + ncolT, ivec, nrowD, 1) ;
IVshuffle(offset + ncolT, ivec, ++seed) ;
SubMtx_columnIndices(mtxU, &ncolU, &colindU) ;
IVcopy(ncolU, colindU, ivec) ;
IVqsortUp(ncolU, colindU) ;
IVfree(ivec) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize U SubMtx object",
        t2 - t1) ;
fprintf(msgFile, "\n U = zeros(%d,%d) ;", nrowD, size) ;
SubMtx_writeForMatlab(mtxU, "U", msgFile) ;
if ( CHV_IS_NONSYMMETRIC(chvT) ) {
/*
   ----------------------------
   initialize the L SubMtx object
   ----------------------------
*/
   MARKTIME(t1) ;
   mtxL = SubMtx_new() ;
   if ( CHV_IS_REAL(chvT) ) {
      if ( sparsityflag == 0 ) {
         SubMtx_initRandom(mtxL, SPOOLES_REAL, SUBMTX_DENSE_ROWS,
                         0, 0, ncolU, nrowD, nentU, ++seed) ;
      } else {
         SubMtx_initRandom(mtxL, SPOOLES_REAL, SUBMTX_SPARSE_ROWS,
                         0, 0, ncolU, nrowD, nentU, ++seed) ;
      }
   } else if ( CHV_IS_COMPLEX(chvT) ) {
      if ( sparsityflag == 0 ) {
         SubMtx_initRandom(mtxL, SPOOLES_COMPLEX, SUBMTX_DENSE_ROWS,
                         0, 0, ncolU, nrowD, nentU, ++seed) ;
      } else {
         SubMtx_initRandom(mtxL, SPOOLES_COMPLEX, SUBMTX_SPARSE_ROWS,
                         0, 0, ncolU, nrowD, nentU, ++seed) ;
      }
   }
   SubMtx_rowIndices(mtxL, &nrowL, &rowindL) ;
   IVcopy(nrowL, rowindL, colindU) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n %% CPU : %.3f to initialize L SubMtx object",
           t2 - t1) ;
   fprintf(msgFile, "\n L = zeros(%d,%d) ;", size, nrowD) ;
   SubMtx_writeForMatlab(mtxL, "L", msgFile) ;
} else {
   mtxL = NULL ;
}
/*
   --------------------------------
   compute the matrix-matrix update
   --------------------------------
*/
tempDV = DV_new() ;
ops = 8*nrowD*nrowD*ncolU ;
if ( CHV_IS_SYMMETRIC(chvT) ) {
   Chv_updateS(chvT, mtxD, mtxU, tempDV) ;
} else if ( CHV_IS_HERMITIAN(chvT) ) {
   Chv_updateH(chvT, mtxD, mtxU, tempDV) ;
} else if ( CHV_IS_NONSYMMETRIC(chvT) ) {
   Chv_updateN(chvT, mtxL, mtxD, mtxU, tempDV) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to compute m-m, %.3f mflops",
        t2 - t1, ops*1.e-6/(t2 - t1)) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% Z Chv object") ;
   fprintf(msgFile, "\n Z = zeros(%d,%d); ", size, size) ;
   Chv_writeForMatlab(chvT, "Z", msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------
   check with matlab
   -----------------
*/
if ( msglvl > 1 ) {
   if ( CHV_IS_HERMITIAN(chvT) ) {
      fprintf(msgFile, "\n\n B  =  ctranspose(U) * D * U ;") ;
   } else if ( CHV_IS_SYMMETRIC(chvT) ) {
      fprintf(msgFile, "\n\n B  =  transpose(U) * D * U ;") ;
   } else {
      fprintf(msgFile, "\n\n B  =  L * D * U ;") ;
   }
   fprintf(msgFile, 
           "\n\n for irow = 1:%d"
           "\n      for jcol = 1:%d"
           "\n         if T(irow,jcol) ~= 0.0"
           "\n            T(irow,jcol) = T(irow,jcol) - B(irow,jcol) ;"
           "\n         end"
           "\n      end"
           "\n   end"
           "\n emtx   = abs(Z - T) ;",
           size, size) ;
   fprintf(msgFile, "\n\n maxabs = max(max(emtx)) ") ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( mtxL != NULL ) {
   SubMtx_free(mtxL) ;
}
Chv_free(chvT) ;
SubMtx_free(mtxD) ;
SubMtx_free(mtxU) ;
DV_free(tempDV) ;
Drand_free(drand) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
