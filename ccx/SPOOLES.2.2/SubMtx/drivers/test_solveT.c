/*  test_solveT.c  */

#include "../SubMtx.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------
   test the SubMtx_solveT() method.

   created -- 98may01, cca
   --------------------------------
*/
{
SubMtx   *mtxA, *mtxB, *mtxX ;
double   idot, rdot, t1, t2 ;
double   *entB, *entX ;
Drand    *drand ;
FILE     *msgFile ;
int      inc1, inc2, mode, msglvl, ncolA, nentA, nrowA, 
         ncolB, nrowB, ncolX, nrowX, seed, type ;

if ( argc != 9 ) {
   fprintf(stdout, 
       "\n\n usage : %s msglvl msgFile type mode nrowA nentA ncolB seed"
       "\n    msglvl  -- message level"
       "\n    msgFile -- message file"
       "\n    type    -- type of matrix A"
       "\n       1 -- real"
       "\n       2 -- complex"
       "\n    mode    -- mode of matrix A"
       "\n       2 -- sparse stored by rows"
       "\n       3 -- sparse stored by columns"
       "\n       5 -- sparse stored by subrows"
       "\n       6 -- sparse stored by subcolumns"
       "\n    nrowA   -- # of rows in matrix A"
       "\n    nentA   -- # of entries in matrix A"
       "\n    ncolB   -- # of columns in matrix B"
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
nrowA = atoi(argv[5]) ;
nentA = atoi(argv[6]) ;
ncolB = atoi(argv[7]) ;
seed  = atoi(argv[8]) ;
fprintf(msgFile, "\n %% %s:"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% mode    = %d"
        "\n %% nrowA   = %d"
        "\n %% nentA   = %d"
        "\n %% ncolB   = %d"
        "\n %% seed    = %d",
        argv[0], msglvl, argv[2], type, mode, 
        nrowA, nentA, ncolB, seed) ;
ncolA = nrowA ;
nrowB = nrowA ;
nrowX = nrowA ;
ncolX = ncolB ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if ( nrowA <= 0 || nentA <= 0 || ncolB <= 0 ) {
   fprintf(stderr, "\n invalid input\n") ;
   exit(-1) ;
}
switch ( type ) {
case SPOOLES_REAL :
   switch ( mode ) {
   case SUBMTX_DENSE_SUBROWS :
   case SUBMTX_SPARSE_ROWS :
   case SUBMTX_DENSE_SUBCOLUMNS :
   case SUBMTX_SPARSE_COLUMNS :
      break ;
   default :
      fprintf(stderr, "\n invalid mode %d\n", mode) ;
      exit(-1) ;
   }
   break ;
case SPOOLES_COMPLEX :
   switch ( mode ) {
   case SUBMTX_DENSE_SUBROWS :
   case SUBMTX_SPARSE_ROWS :
   case SUBMTX_DENSE_SUBCOLUMNS :
   case SUBMTX_SPARSE_COLUMNS :
      break ;
   default :
      fprintf(stderr, "\n invalid mode %d\n", mode) ;
      exit(-1) ;
   }
   break ;
default :
   fprintf(stderr, "\n invalid type %d\n", type) ;
   exit(-1) ;
   break ;
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
SubMtx_initRandom(mtxX, type, SUBMTX_DENSE_COLUMNS, 0, 0, 
                  nrowX, ncolX, nrowX*ncolX, ++seed) ;
SubMtx_denseInfo(mtxX, &nrowX, &ncolX, &inc1, &inc2, &entX) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize X SubMtx object",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% X SubMtx object") ;
   fprintf(msgFile, "\n X = zeros(%d,%d) ;", nrowX, ncolX) ;
   SubMtx_writeForMatlab(mtxX, "X", msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------
   initialize the B SubMtx object
   ------------------------------
*/
MARKTIME(t1) ;
mtxB = SubMtx_new() ;
SubMtx_init(mtxB, type,
            SUBMTX_DENSE_COLUMNS, 0, 0, nrowB, ncolB, nrowB*ncolB) ;
SubMtx_denseInfo(mtxB, &nrowB, &ncolB, &inc1, &inc2, &entB) ;
if ( SUBMTX_IS_REAL(mtxX) ) {
    DVcopy(nrowB*ncolB, entB, entX) ;
} else if ( SUBMTX_IS_COMPLEX(mtxX) ) {
   ZVcopy(nrowB*ncolB, entB, entX) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize B SubMtx object",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% B SubMtx object") ;
   fprintf(msgFile, "\n B = zeros(%d,%d) ;", nrowB, ncolB) ;
   SubMtx_writeForMatlab(mtxB, "B", msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------
   initialize the A matrix SubMtx object
   -------------------------------------
*/
seed++ ;
mtxA = SubMtx_new() ;
switch ( mode ) {
case SUBMTX_DENSE_SUBROWS :
case SUBMTX_SPARSE_ROWS :
   SubMtx_initRandomLowerTriangle(mtxA, type, mode, 0, 0, 
                                  nrowA, ncolA, nentA, seed, 1) ;
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
case SUBMTX_SPARSE_COLUMNS :
   SubMtx_initRandomUpperTriangle(mtxA, type, mode, 0, 0, 
                                  nrowA, ncolA, nentA, seed, 1) ;
   break ;
default :
   fprintf(stderr, "\n fatal error in test_solve"
           "\n invalid mode = %d", mode) ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% A SubMtx object") ;
   fprintf(msgFile, "\n A = zeros(%d,%d) ;", nrowA, ncolA) ;
   SubMtx_writeForMatlab(mtxA, "A", msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------------
   compute B = (I + A^T) * X (for lower and upper triangular)
   ----------------------------------------------------------
*/
if ( SUBMTX_IS_REAL(mtxA) ) {
   DV       *colDV, *rowDV ;
   double   *colX, *rowA, *pBij, *pXij ;
   int      irowA, jcolX ;

   colDV = DV_new() ;
   DV_init(colDV, nrowA, NULL) ;
   colX = DV_entries(colDV) ;
   rowDV = DV_new() ;
   DV_init(rowDV, nrowA, NULL) ;
   rowA = DV_entries(rowDV) ;
   for ( jcolX = 0 ; jcolX < ncolB ; jcolX++ ) {
      SubMtx_fillColumnDV(mtxX, jcolX, colDV) ;
      for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
         SubMtx_fillColumnDV(mtxA, irowA, rowDV) ;
         SubMtx_locationOfRealEntry(mtxX, irowA, jcolX, &pXij) ;
         SubMtx_locationOfRealEntry(mtxB, irowA, jcolX, &pBij) ;
         *pBij = *pXij + DVdot(nrowA, rowA, colX) ;
      }
   }
   DV_free(colDV) ;
   DV_free(rowDV) ;
} else if ( SUBMTX_IS_COMPLEX(mtxA) ) {
   ZV       *colZV, *rowZV ;
   double   *colX, *rowA, *pBIij, *pBRij, *pXIij, *pXRij ;
   int      irowA, jcolX ;

   colZV = ZV_new() ;
   ZV_init(colZV, nrowA, NULL) ;
   colX = ZV_entries(colZV) ;
   rowZV = ZV_new() ;
   ZV_init(rowZV, nrowA, NULL) ;
   rowA = ZV_entries(rowZV) ;
   for ( jcolX = 0 ; jcolX < ncolB ; jcolX++ ) {
      SubMtx_fillColumnZV(mtxX, jcolX, colZV) ;
      for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
         SubMtx_fillColumnZV(mtxA, irowA, rowZV) ;
         SubMtx_locationOfComplexEntry(mtxX, 
                                       irowA, jcolX, &pXRij, &pXIij) ;
         SubMtx_locationOfComplexEntry(mtxB, 
                                       irowA, jcolX, &pBRij, &pBIij) ;
         ZVdotU(nrowA, rowA, colX, &rdot, &idot) ;
         *pBRij = *pXRij + rdot ;
         *pBIij = *pXIij + idot ;
      }
   }
   ZV_free(colZV) ;
   ZV_free(rowZV) ;
}
/*
   ----------------------
   print out the matrices
   ----------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% X SubMtx object") ;
   fprintf(msgFile, "\n X = zeros(%d,%d) ;", nrowX, ncolX) ;
   SubMtx_writeForMatlab(mtxX, "X", msgFile) ;
   fprintf(msgFile, "\n\n %% A SubMtx object") ;
   fprintf(msgFile, "\n A = zeros(%d,%d) ;", nrowA, ncolA) ;
   SubMtx_writeForMatlab(mtxA, "A", msgFile) ;
   fprintf(msgFile, "\n\n %% B SubMtx object") ;
   fprintf(msgFile, "\n B = zeros(%d,%d) ;", nrowB, ncolB) ;
   SubMtx_writeForMatlab(mtxB, "B", msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------
   check with matlab
   -----------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile,
           "\n\n emtx   = abs(B - X - transpose(A)*X) ;"
           "\n\n condA = cond(eye(%d,%d) + transpose(A))"
           "\n\n maxabsZ = max(max(abs(emtx))) ", nrowA, nrowA) ;
   fflush(msgFile) ;
}
/*
   --------------------------------
   compute the solve (I + A^T)Y = B
   --------------------------------
*/
SubMtx_solveT(mtxA, mtxB) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% Y SubMtx object") ;
   fprintf(msgFile, "\n Y = zeros(%d,%d) ;", nrowB, ncolB) ;
   SubMtx_writeForMatlab(mtxB, "Y", msgFile) ;
   fprintf(msgFile,
           "\n\n %% solerror   = abs(Y - X) ;"
           "\n\n solerror   = abs(Y - X) ;"
           "\n\n maxabserror = max(max(solerror)) ") ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
SubMtx_free(mtxA) ;
SubMtx_free(mtxX) ;
SubMtx_free(mtxB) ;
Drand_free(drand) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
