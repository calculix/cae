/*  testQR.c  */

#include "../A2.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------------------------
   test the A2_QRreduce(), A2_QRcomputeQ() and A2_applyQT() methods

   created -- 98dec10, cca
   ----------------------------------------------------------------
*/
{
A2       *A, *Q, *R, *X, *Y ;
double   nops, t1, t2 ;
DV       workDV ;
FILE     *msgFile ;
int      inc1, inc2, irow, jcol,
         msglvl, nrow, ncol, ncolX, seed, type ;

if ( argc != 10 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile type nrow ncol inc1 inc2 seed ncolX"
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    type    -- entries type"
"\n      1 -- real"
"\n      2 -- complex"
"\n    nrow    -- # of rows in A"
"\n    ncol    -- # of columns in A"
"\n    inc1    -- row increment, must be ncol"
"\n    inc2    -- column increment, must be 1"
"\n    seed    -- random number seed"
"\n    ncolX   -- # of columns in X"
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
type = atoi(argv[3]) ;
nrow = atoi(argv[4]) ;
ncol = atoi(argv[5]) ;
inc1 = atoi(argv[6]) ;
inc2 = atoi(argv[7]) ;
if (   type < 1 || type > 2 || nrow < 0 || ncol < 0 
    || inc1 < 1 || inc2 < 1 ) {
   fprintf(stderr, 
       "\n fatal error, type %d, nrow %d, ncol %d, inc1 %d, inc2 %d",
       type, nrow, ncol, inc1, inc2) ;
   exit(-1) ;
}
seed   = atoi(argv[8]) ;
ncolX  = atoi(argv[9]) ;
fprintf(msgFile, "\n\n %% %s :"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% nrow    = %d"
        "\n %% ncol    = %d"
        "\n %% inc1    = %d"
        "\n %% inc2    = %d"
        "\n %% seed    = %d"
        "\n %% ncolX   = %d"
        "\n",
        argv[0], msglvl, argv[2], type, 
        nrow, ncol, inc1, inc2, seed, ncolX) ;
if ( inc1 != 1 && inc2 != 1 ) {
   fprintf(stderr, "\n inc1 = %d, inc2 = %d\n", inc1, inc2) ;
   exit(-1) ;
}
/*
   -----------------------------
   initialize the matrix objects
   -----------------------------
*/
MARKTIME(t1) ;
A = A2_new() ;
A2_init(A, type, nrow, ncol, inc1, inc2, NULL) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix object",
        t2 - t1) ;
MARKTIME(t1) ;
A2_fillRandomUniform(A, -1, 1, seed++) ;
MARKTIME(t2) ;
fprintf(msgFile, 
      "\n %% CPU : %.3f to fill matrix with random numbers", t2 - t1) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n matrix A") ;
   A2_writeForHumanEye(A, msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n %% matrix A") ;
   A2_writeForMatlab(A, "A", msgFile) ;
}
/*
   --------------------
   compute the R matrix
   --------------------
*/
DV_setDefaultFields(&workDV) ;
/*
   ------------------
   use rank-1 updates
   ------------------
*/
MARKTIME(t1) ;
nops = A2_QRreduce(A, &workDV, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n %% CPU : rank-1: %.3f to compute R, %.0f ops, %.3f mflops",
        t2 - t1, nops, 1.e-6*nops/(t2-t1)) ;
/*
   ----------------------
   write out the R matrix
   ----------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n %% matrix R") ;
   R = A2_new() ;
   A2_subA2(R, A, 0, ncol - 1, 0, ncol - 1) ;
   A2_writeForMatlab(R, "R", msgFile) ;
   fprintf(msgFile, 
           "\n for jj = 1:%d"
           "\n    for ii = jj+1:%d"
           "\n       R(ii,jj) = 0.0 ;"
           "\n    end"
           "\n end", ncol, ncol) ;
   fflush(msgFile) ;
   A2_free(R) ;
}
/*
   -------------------------------
   check the error |A^H*A - R^H*R|
   -------------------------------
*/
if ( msglvl > 1 ) {
   if ( type == SPOOLES_REAL ) {
      fprintf(msgFile, 
              "\n emtx1 = transpose(A)*A - transpose(R)*R ;"
              "\n error = max(max(abs(emtx1))) " ) ;
   } else if ( type == SPOOLES_COMPLEX ) {
      fprintf(msgFile, 
              "\n emtx1 = ctranspose(A)*A - ctranspose(R)*R ;"
              "\n error = max(max(abs(emtx1))) " ) ;
   }
   DV_clearData(&workDV) ;
}
if ( inc1 == 1 ) {
/*
   ---------
   compute Q
   ---------
*/
   Q = A2_new() ;
   A2_init(Q, type, nrow, ncol, inc1, inc2, NULL) ;
   A2_computeQ(Q, A, &workDV, msglvl, msgFile) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %% matrix Q") ;
      A2_writeForMatlab(Q, "Q", msgFile) ;
      fprintf(msgFile, 
              "\n emtx2 = A - Q * R ;"
              "\n error = max(max(abs(emtx2))) " ) ;
   }
   A2_free(Q) ;
/*
   ---------------------------------------------------
   create a matrix X with the same number of rows as A
   ---------------------------------------------------
*/
   MARKTIME(t1) ;
   X = A2_new() ;
   A2_init(X, type, nrow, ncolX, inc1, inc2, NULL) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix object",
           t2 - t1) ;
   MARKTIME(t1) ;
   A2_fillRandomUniform(X, -1, 1, seed++) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
      "\n %% CPU : %.3f to fill matrix with random numbers", t2 - t1) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n matrix X") ;
      A2_writeForHumanEye(X, msgFile) ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %% matrix X") ;
      A2_writeForMatlab(X, "X", msgFile) ;
   }
/*
   -------------------
   compute Y = Q^T * X
   -------------------
*/
   MARKTIME(t1) ;
   Y = A2_new() ;
   A2_init(Y, type, nrow, ncolX, inc1, inc2, NULL) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix object",
           t2 - t1) ;
   A2_applyQT(Y, A, X, &workDV, msglvl, msgFile) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %% matrix Y") ;
      A2_writeForMatlab(Y, "Y", msgFile) ;
      fprintf(msgFile, "\n [Qexact,Rexact] = qr(A) ;") ;
      if ( A2_IS_REAL(A) ) {
         fprintf(msgFile, "\n emtx3 = Y - transpose(Qexact) * X ;") ;
      } else {
         fprintf(msgFile, "\n emtx3 = Y - ctranspose(Qexact) * X ;") ;
      }
      fprintf(msgFile, "\n error = max(max(abs(emtx3))) " ) ;
   }
   A2_free(Q) ;
   A2_free(X) ;
   A2_free(Y) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_clearData(&workDV) ;
A2_free(A) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
