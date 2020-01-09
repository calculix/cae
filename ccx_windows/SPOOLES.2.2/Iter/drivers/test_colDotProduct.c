/*  test_colDotProduct.c  */

#include "../Iter.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------
   test the DenseMtx_colDotProduct routine.
   value = x^H * y, where x and y are i-th and j-th columns of a 
   DenseMtx A respectively.  (In real application, y could be from
   another DenseMtx.)

   when msglvl > 1, the output of this program
   can be fed into Matlab to check for errors

   created -- 98dec03, jwu
   -----------------------------------------------
*/
{
DenseMtx   *A ;
double     t1, t2, value[2] ;
Drand      *drand ;
FILE       *msgFile ;
int        inc1, inc2, icol, jcol, msglvl, nrow, ncol, rc, seed, type ;

if ( argc != 11 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile type nrow ncol inc1 inc2 "
"\n         icol, jcol, seed "
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    type    -- entries type"
"\n      1 -- real"
"\n      2 -- complex"
"\n    nrow    -- # of rows "
"\n    ncol    -- # of columns "
"\n    inc1    -- row increment "
"\n    inc2    -- column increment "
"\n    icol    -- vector x: i-th column of A "
"\n    jcol    -- vector y: j-th column of A "
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
icol = atoi(argv[8]) ;
jcol = atoi(argv[9]) ;
seed   = atoi(argv[10]) ;
fprintf(msgFile, "\n\n %% %s :"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% nrow    = %d"
        "\n %% ncol    = %d"
        "\n %% inc1    = %d"
        "\n %% inc2    = %d"
        "\n %% icol    = %d"
        "\n %% jcol    = %d"
        "\n %% seed    = %d"
        "\n",
        argv[0], msglvl, argv[2], type, nrow, ncol, inc1, inc2, 
        icol, jcol, seed) ;
/*
   ----------------------------
   initialize the matrix object
   ----------------------------
*/
MARKTIME(t1) ;
A = DenseMtx_new() ;
DenseMtx_init(A, type, 0, 0, nrow, ncol, inc1, inc2) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix object",
        t2 - t1) ;
MARKTIME(t1) ;
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
seed++ ;
Drand_setUniform(drand, -1.0, 1.0) ;
DenseMtx_fillRandomEntries(A, drand) ;
MARKTIME(t2) ;
fprintf(msgFile, 
      "\n %% CPU : %.3f to fill matrix with random numbers", t2 - t1) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n matrix A") ;
   DenseMtx_writeForHumanEye(A, msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n %% matrix A") ;
   fprintf(msgFile, "\n nrow = %d ;", nrow) ;
   fprintf(msgFile, "\n ncol = %d ;", ncol) ;
   fprintf(msgFile, "\n");
   DenseMtx_writeForMatlab(A, "A", msgFile) ;
}
/*
   ----------------
   compute the inner product
   ----------------
*/
DenseMtx_colDotProduct(A, icol, A, jcol, value);

if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% Dot product A(:,%d)^H * A(:,%d):",icol+1, jcol+1) ;
   fprintf(msgFile, "\n real = %e", value[0]) ;
   if (DENSEMTX_IS_COMPLEX(A)) fprintf(msgFile, "\n imag = %e", value[1]);
   fprintf(msgFile, "\n");
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DenseMtx_free(A) ;
Drand_free(drand) ;

return(1) ; }

/*--------------------------------------------------------------------*/
