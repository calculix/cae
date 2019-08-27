/*  test_solveupd.c  */

#include "../SubMtx.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------
   test the SubMtx_solveupd() method.

   created -- 98may02, cca
   ----------------------------------
*/
{
SubMtx   *mtxA, *mtxX, *mtxY ;
double   ops, t1, t2 ;
double   *entX, *entY ;
Drand    *drand ;
FILE     *msgFile ;
int      inc1, inc2, irowA, jcolX, mode, msglvl, ncolA, nentA, nrowA, 
         ncolX, nrowX, ncolY, nrowY, seed, type ;
int      *colindA, *ivec, *rowindA ;

if ( argc != 12 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile type mode"
           "\n         nrowY ncolY nrowA ncolA nentA nrowX seed"
           "\n    msglvl  -- message level"
           "\n    msgFile -- message file"
           "\n    type    -- type of matrix A"
           "\n       1 -- real"
           "\n       2 -- complex"
           "\n    mode    -- mode of matrix A"
           "\n       0 -- dense stored by rows"
           "\n       1 -- dense stored by columns"
           "\n       2 -- sparse stored by rows"
           "\n       3 -- sparse stored by columns"
           "\n    nrowY   -- # of rows in vector Y"
           "\n    ncolY   -- # of columns in vector Y"
           "\n    nrowA   -- # of rows in matrix A"
           "\n    ncolA   -- # of columns in matrix A"
           "\n    nentA   -- # of entries in matrix A"
           "\n    nrowX   -- # of rows in matrix X"
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
nrowY = atoi(argv[5]) ;
ncolY = atoi(argv[6]) ;
nrowA = atoi(argv[7]) ;
ncolA = atoi(argv[8]) ;
nentA = atoi(argv[9]) ;
nrowX = atoi(argv[10]) ;
seed  = atoi(argv[11]) ;
fprintf(msgFile, "\n %% %s:"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% mode    = %d"
        "\n %% nrowY   = %d"
        "\n %% ncolY   = %d"
        "\n %% nrowA   = %d"
        "\n %% ncolA   = %d"
        "\n %% nentA   = %d"
        "\n %% nrowX   = %d"
        "\n %% seed    = %d",
        argv[0], msglvl, argv[2], type, mode, nrowY, ncolY, 
        nrowA, ncolA, nentA, nrowX, seed) ;
ncolX = ncolY ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if (  nrowA <= 0 || nrowA > nrowY
   || ncolA <= 0 || ncolA > nrowX
   || nentA > nrowA*ncolA 
   || nrowX <= 0 ) {
   fprintf(stderr, "\n invalid input\n") ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n invalid type %d\n", type) ;
   exit(-1) ;
}
switch ( mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS :
   break ;
default :
   fprintf(stderr, "\n invalid mode %d\n", mode) ;
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
   ------------------------------
   initialize the X SubMtx object
   ------------------------------
*/
MARKTIME(t1) ;
mtxX = SubMtx_new() ;
SubMtx_init(mtxX, type, SUBMTX_DENSE_COLUMNS, 0, 0, 
            nrowX, ncolX, nrowX*ncolX) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
if ( SUBMTX_IS_REAL(mtxX) ) {
   Drand_fillDvector(drand, nrowX*ncolX, entX) ;
} else if ( SUBMTX_IS_COMPLEX(mtxX) ) {
   Drand_fillDvector(drand, 2*nrowX*ncolX, entX) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize X SubMtx object",
        t2 - t1) ;
/*
   ------------------------------
   initialize the Y SubMtx object
   ------------------------------
*/
MARKTIME(t1) ;
mtxY = SubMtx_new() ;
SubMtx_init(mtxY, type, SUBMTX_DENSE_COLUMNS, 0, 0, 
            nrowY, ncolY, nrowY*ncolY) ;
SubMtx_denseInfo(mtxY, &nrowY, &ncolY, &inc1, &inc2, &entY) ;
if ( SUBMTX_IS_REAL(mtxX) ) {
   Drand_fillDvector(drand, nrowY*ncolY, entY) ;
   DVzero(nrowY*ncolY, entY) ;
} else if ( SUBMTX_IS_COMPLEX(mtxX) ) {
   Drand_fillDvector(drand, 2*nrowY*ncolY, entY) ;
   DVzero(2*nrowY*ncolY, entY) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize Y SubMtx object",
        t2 - t1) ;
/*
   -------------------------------------
   initialize the A matrix SubMtx object
   -------------------------------------
*/
mtxA  = SubMtx_new() ;
SubMtx_initRandom(mtxA, type, mode, 0, 0, nrowA, ncolA, nentA, seed) ;
/*
   -------------------------
   load the row indices of A
   -------------------------
*/
SubMtx_rowIndices(mtxA, &nrowA, &rowindA) ;
ivec = IVinit(nrowY, -1) ;
IVramp(nrowY, ivec, 0, 1) ;
IVshuffle(nrowY, ivec, seed+1) ;
IVcopy(nrowA, rowindA, ivec) ;
IVqsortUp(nrowA, rowindA) ;
IVfree(ivec) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% row indices of A") ;
   IVfprintf(msgFile, nrowA, rowindA) ;
   fflush(msgFile) ;
}
/*
   ----------------------------
   load the column indices of A
   ----------------------------
*/
SubMtx_columnIndices(mtxA, &ncolA, &colindA) ;
ivec = IVinit(nrowX, -1) ;
IVramp(nrowX, ivec, 0, 1) ;
IVshuffle(nrowX, ivec, seed+2) ;
IVcopy(ncolA, colindA, ivec) ;
IVqsortUp(ncolA, colindA) ;
IVfree(ivec) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n %% column indices of A") ;
   IVfprintf(msgFile, ncolA, colindA) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   compute the matrix-matrix multiply
   ----------------------------------
*/
if ( type == SPOOLES_REAL ) {
   double   *colX, *pYij, *rowA ;
   double   sum ;
   DV       *colDV, *rowDV ;
   int      ii ;

   ops = 2*nrowA*ncolA*ncolX ;
   colDV = DV_new() ;
   DV_init(colDV, nrowX, NULL) ;
   colX = DV_entries(colDV) ;
   rowDV = DV_new() ;
   DV_init(rowDV, ncolA, NULL) ;
   rowA = DV_entries(rowDV) ;
   for ( jcolX = 0 ; jcolX < ncolX ; jcolX++ ) {
      SubMtx_fillColumnDV(mtxX, jcolX, colDV) ;
      for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
         SubMtx_fillRowDV(mtxA, irowA, rowDV) ;
         if ( ncolA == nrowX ) {
            for ( ii = 0, sum = 0.0 ; ii < ncolA ; ii++ ) {
               sum += rowA[ii] * colX[ii] ;
            }
         } else {
            for ( ii = 0, sum = 0.0 ; ii < ncolA ; ii++ ) {
               sum += rowA[ii] * colX[colindA[ii]] ;
            }
         }
         if ( nrowA == nrowY ) {
            SubMtx_locationOfRealEntry(mtxY, irowA, jcolX, &pYij) ;
         } else {
            SubMtx_locationOfRealEntry(mtxY, rowindA[irowA], jcolX,
                                       &pYij) ;
         }
         *pYij = sum ;
      }
   }
   DV_free(colDV) ;
   DV_free(rowDV) ;
} else if ( type == SPOOLES_COMPLEX ) {
   double   *colX, *pYIij, *pYRij, *rowA ;
   double   idot, rdot ;
   ZV       *colZV, *rowZV ;

   ops = 8*nrowA*ncolA*ncolX ;
   colZV = ZV_new() ;
   ZV_init(colZV, nrowX, NULL) ;
   colX = ZV_entries(colZV) ;
   rowZV = ZV_new() ;
   ZV_init(rowZV, ncolA, NULL) ;
   rowA = ZV_entries(rowZV) ;
   for ( jcolX = 0 ; jcolX < ncolX ; jcolX++ ) {
      SubMtx_fillColumnZV(mtxX, jcolX, colZV) ;
      for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
         SubMtx_fillRowZV(mtxA, irowA, rowZV) ;
         if ( ncolA == nrowX ) {
            ZVdotU(ncolA, colX, rowA, &rdot, &idot) ;
         } else {
            ZVdotiU(ncolA, colX, colindA, rowA, &rdot, &idot) ;
         }
         if ( nrowA == nrowY ) {
            SubMtx_locationOfComplexEntry(mtxY, 
                                 irowA, jcolX, &pYRij, &pYIij) ;
         } else {
            SubMtx_locationOfComplexEntry(mtxY, rowindA[irowA], jcolX,
                                 &pYRij, &pYIij) ;
         }
         *pYRij = rdot ;
         *pYIij = idot ;
      }
   }
   ZV_free(colZV) ;
   ZV_free(rowZV) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to compute m-m, %.3f mflops",
        t2 - t1, ops*1.e-6/(t2 - t1)) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% Z SubMtx object") ;
   fprintf(msgFile, "\n Z = zeros(%d,%d) ;", nrowY, ncolY) ;
   SubMtx_writeForMatlab(mtxY, "Z", msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------
   print out the matrices
   ----------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% Y SubMtx object") ;
   fprintf(msgFile, "\n Y = zeros(%d,%d) ;", nrowY, ncolY) ;
   SubMtx_writeForMatlab(mtxY, "Y", msgFile) ;
   fprintf(msgFile, "\n\n %% A SubMtx object") ;
   fprintf(msgFile, "\n A = zeros(%d,%d) ;", nrowY, nrowX) ;
   SubMtx_writeForMatlab(mtxA, "A", msgFile) ;
   fprintf(msgFile, "\n\n %% X SubMtx object") ;
   fprintf(msgFile, "\n X = zeros(%d,%d) ;", nrowX, ncolY) ;
   SubMtx_writeForMatlab(mtxX, "X", msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------
   check with matlab
   -----------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile,
           "\n\n emtx   = abs(Y - A*X) ;"
           "\n\n maxabs = max(max(emtx)) ") ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   compute the update Y := Y - A*X
   (Y should now be zero)
   -------------------------------
*/
SubMtx_solveupd(mtxY, mtxA, mtxX) ;
/*
   ----------------------
   print out the Y matrix
   ----------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% Z SubMtx object") ;
   fprintf(msgFile, "\n Z = zeros(%d,%d) ;", nrowY, ncolY) ;
   SubMtx_writeForMatlab(mtxY, "Z", msgFile) ;
   fflush(msgFile) ;
}
fprintf(msgFile, "\n RES %4d %4d %4d %4d %4d %4d %4d %4d %12.4e", 
        type, mode, nrowY, ncolY, nrowA, ncolA, nrowX, ncolX, 
        SubMtx_maxabs(mtxY)) ;
if ( msglvl > 1 ) {
   fprintf(msgFile,
           "\n\n emtx   = abs(Z) ;"
           "\n\n maxerr = max(max(emtx)) ") ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
SubMtx_free(mtxA) ;
SubMtx_free(mtxX) ;
SubMtx_free(mtxY) ;
Drand_free(drand) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
