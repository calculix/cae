/*  test_scalevec.c  */

#include "../SubMtx.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------
   test the SubMtx_scale{1,2,3}vec() methods.

   created -- 98may02, cca
   ------------------------------------------
*/
{
SubMtx   *mtxA ;
double   t1, t2 ;
double   *x0, *x1, *x2, *y0, *y1, *y2 ;
Drand    *drand ;
DV       *xdv0, *xdv1, *xdv2, *ydv0, *ydv1, *ydv2 ;
ZV       *xzv0, *xzv1, *xzv2, *yzv0, *yzv1, *yzv2 ;
FILE     *msgFile ;
int      mode, msglvl, nrowA, seed, type ;

if ( argc != 7 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile type nrowA seed"
           "\n    msglvl  -- message level"
           "\n    msgFile -- message file"
           "\n    type -- type of matrix A"
           "\n       1 -- real"
           "\n       2 -- complex"
           "\n    mode -- mode of matrix A"
           "\n       7 -- diagonal"
           "\n       8 -- block diagonal symmetric"
           "\n       9 -- block diagonal hermitian (complex only)"
           "\n    nrowA -- # of rows in matrix A"
           "\n    seed  -- random number seed"
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
seed  = atoi(argv[6]) ;
fprintf(msgFile, "\n %% %s:"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% mode    = %d"
        "\n %% nrowA   = %d"
        "\n %% seed    = %d",
        argv[0], msglvl, argv[2], type, mode, nrowA, seed) ;
/*
   -----------------------------
   check for errors in the input
   -----------------------------
*/
if ( nrowA <= 0 ) {
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
Drand_setSeed(drand, seed) ;
Drand_setNormal(drand, 0.0, 1.0) ;
/*
   ----------------------------
   initialize the X ZV objects
   ----------------------------
*/
MARKTIME(t1) ;
if ( type == SPOOLES_REAL ) {
   xdv0 = DV_new() ;
   DV_init(xdv0, nrowA, NULL) ;
   x0 = DV_entries(xdv0) ;
   Drand_fillDvector(drand, nrowA, x0) ;
   xdv1 = DV_new() ;
   DV_init(xdv1, nrowA, NULL) ;
   x1 = DV_entries(xdv1) ;
   Drand_fillDvector(drand, nrowA, x1) ;
   xdv2 = DV_new() ;
   DV_init(xdv2, nrowA, NULL) ;
   x2 = DV_entries(xdv2) ;
   Drand_fillDvector(drand, nrowA, x2) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n %% CPU : %.3f to initialize X ZV objects",
           t2 - t1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n %% X DV objects") ;
      fprintf(msgFile, "\n X0 = zeros(%d,1) ;", nrowA) ;
      DV_writeForMatlab(xdv0, "X0", msgFile) ;
      fprintf(msgFile, "\n X1 = zeros(%d,1) ;", nrowA) ;
      DV_writeForMatlab(xdv1, "X1", msgFile) ;
      fprintf(msgFile, "\n X2 = zeros(%d,1) ;", nrowA) ;
      DV_writeForMatlab(xdv2, "X2", msgFile) ;
      fflush(msgFile) ;
   }
} else if ( type == SPOOLES_COMPLEX ) {
   xzv0 = ZV_new() ;
   ZV_init(xzv0, nrowA, NULL) ;
   x0 = ZV_entries(xzv0) ;
   Drand_fillDvector(drand, 2*nrowA, x0) ;
   xzv1 = ZV_new() ;
   ZV_init(xzv1, nrowA, NULL) ;
   x1 = ZV_entries(xzv1) ;
   Drand_fillDvector(drand, 2*nrowA, x1) ;
   xzv2 = ZV_new() ;
   ZV_init(xzv2, nrowA, NULL) ;
   x2 = ZV_entries(xzv2) ;
   Drand_fillDvector(drand, 2*nrowA, x2) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n %% CPU : %.3f to initialize X ZV objects",
           t2 - t1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n %% X ZV objects") ;
      fprintf(msgFile, "\n X0 = zeros(%d,1) ;", nrowA) ;
      ZV_writeForMatlab(xzv0, "X0", msgFile) ;
      fprintf(msgFile, "\n X1 = zeros(%d,1) ;", nrowA) ;
      ZV_writeForMatlab(xzv1, "X1", msgFile) ;
      fprintf(msgFile, "\n X2 = zeros(%d,1) ;", nrowA) ;
      ZV_writeForMatlab(xzv2, "X2", msgFile) ;
      fflush(msgFile) ;
   }
}
/*
   ---------------------------------
   initialize the Y DV or ZV objects
   ---------------------------------
*/
MARKTIME(t1) ;
if ( type == SPOOLES_REAL ) {
   ydv0 = DV_new() ;
   DV_init(ydv0, nrowA, NULL) ;
   y0 = DV_entries(ydv0) ;
   ydv1 = DV_new() ;
   DV_init(ydv1, nrowA, NULL) ;
   y1 = DV_entries(ydv1) ;
   ydv2 = DV_new() ;
   DV_init(ydv2, nrowA, NULL) ;
   y2 = DV_entries(ydv2) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n %% CPU : %.3f to initialize Y DV objects",
           t2 - t1) ;
} else if ( type == SPOOLES_COMPLEX ) {
   yzv0 = ZV_new() ;
   ZV_init(yzv0, nrowA, NULL) ;
   y0 = ZV_entries(yzv0) ;
   yzv1 = ZV_new() ;
   ZV_init(yzv1, nrowA, NULL) ;
   y1 = ZV_entries(yzv1) ;
   yzv2 = ZV_new() ;
   ZV_init(yzv2, nrowA, NULL) ;
   y2 = ZV_entries(yzv2) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n %% CPU : %.3f to initialize Y ZV objects",
           t2 - t1) ;
}
/*
   -----------------------------------
   initialize the A matrix SubMtx object
   -----------------------------------
*/
seed++ ;
mtxA = SubMtx_new() ;
switch ( mode ) {
case SUBMTX_DIAGONAL :
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   SubMtx_initRandom(mtxA, type, mode, 0, 0, nrowA, nrowA, 0, seed) ;
   break ;
default :
   fprintf(stderr, "\n fatal error in test_solve"
           "\n invalid mode = %d", mode) ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% A SubMtx object") ;
   fprintf(msgFile, "\n A = zeros(%d,%d) ;", nrowA, nrowA) ;
   SubMtx_writeForMatlab(mtxA, "A", msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------
   compute Y0 = A * X0
   -------------------
*/
if ( type == SPOOLES_REAL ) {
   DVzero(nrowA, y0) ;
} else if ( type == SPOOLES_COMPLEX ) {
   DVzero(2*nrowA, y0) ;
}
SubMtx_scale1vec(mtxA, y0, x0) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n Z0 = zeros(%d,1) ;", nrowA) ;
   if ( type == SPOOLES_REAL ) {
      DV_writeForMatlab(ydv0, "Z0", msgFile) ;
   } else if ( type == SPOOLES_COMPLEX ) {
      ZV_writeForMatlab(yzv0, "Z0", msgFile) ;
   }
   fprintf(msgFile, "\n err0 = Z0 - A*X0 ;") ;
   fprintf(msgFile, "\n error0 = max(abs(err0))") ;
   fflush(msgFile) ;
}
if ( type == SPOOLES_REAL ) {
   DVzero(nrowA, y0) ;
   DVzero(nrowA, y1) ;
} else if ( type == SPOOLES_COMPLEX ) {
   DVzero(2*nrowA, y0) ;
   DVzero(2*nrowA, y1) ;
}
SubMtx_scale2vec(mtxA, y0, y1, x0, x1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n Z0 = zeros(%d,1) ;", nrowA) ;
   fprintf(msgFile, "\n\n Z1 = zeros(%d,1) ;", nrowA) ;
   if ( type == SPOOLES_REAL ) {
      DV_writeForMatlab(ydv0, "Z0", msgFile) ;
      DV_writeForMatlab(ydv1, "Z1", msgFile) ;
   } else if ( type == SPOOLES_COMPLEX ) {
      ZV_writeForMatlab(yzv0, "Z0", msgFile) ;
      ZV_writeForMatlab(yzv1, "Z1", msgFile) ;
   }
   fprintf(msgFile, "\n err1 = [Z0 Z1] - A*[X0 X1] ;") ;
   fprintf(msgFile, "\n error1 = max(abs(err1))") ;
   fflush(msgFile) ;
}
if ( type == SPOOLES_REAL ) {
   DVzero(nrowA, y0) ;
   DVzero(nrowA, y1) ;
   DVzero(nrowA, y2) ;
} else if ( type == SPOOLES_COMPLEX ) {
   DVzero(2*nrowA, y0) ;
   DVzero(2*nrowA, y1) ;
   DVzero(2*nrowA, y2) ;
}
SubMtx_scale3vec(mtxA, y0, y1, y2, x0, x1, x2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n Z0 = zeros(%d,1) ;", nrowA) ;
   fprintf(msgFile, "\n\n Z1 = zeros(%d,1) ;", nrowA) ;
   fprintf(msgFile, "\n\n Z2 = zeros(%d,1) ;", nrowA) ;
   if ( type == SPOOLES_REAL ) {
      DV_writeForMatlab(ydv0, "Z0", msgFile) ;
      DV_writeForMatlab(ydv1, "Z1", msgFile) ;
      DV_writeForMatlab(ydv2, "Z2", msgFile) ;
   } else if ( type == SPOOLES_COMPLEX ) {
      ZV_writeForMatlab(yzv0, "Z0", msgFile) ;
      ZV_writeForMatlab(yzv1, "Z1", msgFile) ;
      ZV_writeForMatlab(yzv2, "Z2", msgFile) ;
   }
   fprintf(msgFile, "\n err2 = [Z0 Z1 Z2] - A*[X0 X1 X2] ;") ;
   fprintf(msgFile, "\n error3 = max(abs(err2))") ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
SubMtx_free(mtxA) ;
if ( type == SPOOLES_REAL ) {
   DV_free(xdv0) ;
   DV_free(xdv1) ;
   DV_free(xdv2) ;
   DV_free(ydv0) ;
   DV_free(ydv1) ;
   DV_free(ydv2) ;
} else if ( type == SPOOLES_COMPLEX ) {
   ZV_free(xzv0) ;
   ZV_free(xzv1) ;
   ZV_free(xzv2) ;
   ZV_free(yzv0) ;
   ZV_free(yzv1) ;
   ZV_free(yzv2) ;
}
Drand_free(drand) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
