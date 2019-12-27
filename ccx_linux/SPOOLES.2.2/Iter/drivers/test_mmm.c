/*  test_DenseMtx_mmm.c  */

#include "../Iter.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------
   test the DenseMtx_mmm routine.

   C = alpha*A*B + beta*C, where A, B and C  are DenseMtx. 
       alpha and beta are scalars.

   when msglvl > 1, the output of this program
   can be fed into Matlab to check for errors

   created -- 98dec14, ycp
   -------------------------------------------------------
*/
{
DenseMtx   *mtxA, *mtxB, *mtxC;
double     t1, t2, value[2] = {1.0, 1.0} ;
Drand      *drand ;
FILE       *msgFile ;
int        i, j, k, msglvl, nrow, nk, ncol, cnrow, cncol, seed, type ;
int        ainc1, ainc2, binc1, binc2, cinc1, cinc2;
double     alpha[2], beta[2], one[2] = {1.0, 0.0}, rvalue;
char       A_opt[1]=" ", B_opt[1]=" ";

if ( argc != 20 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile type nrow nk ncol ainc1 ainc2 binc1 "
"\n         binc2 cinc1 cinc2 A_opt B_opt ralpha ialpha rbeta ibeta seed "
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    type    -- entries type"
"\n      1 -- real"
"\n      2 -- complex"
"\n    nrow    -- # of rows of mtxA "
"\n    nk      -- # of columns of mtxA "
"\n    ncol    -- # of columns of mtxB "
"\n    ainc1   -- A row increment "
"\n    ainc2   -- A column increment "
"\n    binc1   -- B row increment "
"\n    binc2   -- B column increment "
"\n    binc1   -- C row increment "
"\n    binc2   -- C column increment "
"\n    A_opt   -- A option "
"\n    B_opt   -- B option "
"\n    ralpha  -- real(alpha)"
"\n    ialpha  -- imag(alpha)"
"\n    rbeta   -- real(beta)"
"\n    ibeta   -- imag(beta)"
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
nk   = atoi(argv[5]) ;
ncol = atoi(argv[6]) ;
ainc1= atoi(argv[7]) ;
ainc2= atoi(argv[8]) ;
binc1= atoi(argv[9]) ;
binc2= atoi(argv[10]) ;
cinc1= atoi(argv[11]) ;
cinc2= atoi(argv[12]) ;
if (  type < 1 ||  type > 2 ||  nrow < 0 ||  ncol < 0 ||
     ainc1 < 1 || ainc2 < 1 || binc1 < 1 || binc2 < 1  ) {
   fprintf(stderr, 
       "\n fatal error, type %d, nrow %d, ncol %d, ainc1 %d, ainc2 %d"
       ", binc1 %d, binc2 %d", type, nrow, ncol, ainc1, ainc2, binc1, binc2) ;
   exit(-1) ;
}
A_opt[0] = *argv[13] ;
B_opt[0] = *argv[14] ;
alpha[0]= atof (argv[15]);
alpha[1]= atof (argv[16]);
beta[0] = atof (argv[17]);
beta[1] = atof (argv[18]);
seed    = atoi (argv[19]) ;
fprintf(msgFile, "\n\n %% %s :"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% nrow    = %d"
        "\n %% nk      = %d"
        "\n %% ncol    = %d"
        "\n %% ainc1   = %d"
        "\n %% ainc2   = %d"
        "\n %% binc1   = %d"
        "\n %% binc2   = %d"
        "\n %% cinc1   = %d"
        "\n %% cinc2   = %d"
        "\n %% a_opt   = %c"
        "\n %% b_opt   = %c"
        "\n %% ralpha  = %e"
        "\n %% ialpha  = %e"
        "\n %% rbeta   = %e"
        "\n %% ibeta   = %e"
        "\n %% seed    = %d"
        "\n",
        argv[0], msglvl, argv[2], type, nrow, nk, ncol, ainc1, ainc2, 
        binc1, binc2, cinc1, cinc2, A_opt[0], B_opt[0], alpha[0], 
        alpha[1], beta[0], beta[1], seed) ;
/*
   ----------------------------
   initialize the matrix object
   ----------------------------
*/
MARKTIME(t1) ;
mtxA = DenseMtx_new() ;
DenseMtx_init(mtxA, type, 0, 0, nrow, nk, ainc1, ainc2) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix object",
        t2 - t1) ;
MARKTIME(t1) ;
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
seed++ ;
Drand_setUniform(drand, -1.0, 1.0) ;
DenseMtx_fillRandomEntries(mtxA, drand) ;
MARKTIME(t2) ;
fprintf(msgFile, 
      "\n %% CPU : %.3f to fill matrix A with random numbers", t2 - t1) ;
MARKTIME(t1) ;
mtxB = DenseMtx_new() ;
DenseMtx_init(mtxB, type, 0, 0, nk, ncol, binc1, binc2) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix object",
        t2 - t1) ;
MARKTIME(t1) ;
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
seed++ ;
Drand_setUniform(drand, -1.0, 1.0) ;
DenseMtx_fillRandomEntries(mtxB, drand) ;
MARKTIME(t2) ;
fprintf(msgFile,
      "\n %% CPU : %.3f to fill matrix B with random numbers", t2 - t1) ;

cnrow = nrow;
cncol = ncol;
MARKTIME(t1) ;
mtxC = DenseMtx_new() ;
if ( A_opt[0] == 't' || A_opt[0] == 'T' || 
     A_opt[0] == 'c' || A_opt[0] == 'C') {
  cnrow = nk; 
}
if ( B_opt[0] == 't' || B_opt[0] == 'T' ||
     B_opt[0] == 'c' || B_opt[0] == 'C') {
  cncol = nk; 
}
if ( cinc1 == 1 && cinc2 == nrow ){ /* stored by column */
  cinc1 = 1;
  cinc2 = cnrow;
} else { /* stored by row */
  cinc1 = cncol;
  cinc2 = 1; 
}
DenseMtx_init(mtxC, type, 0, 0, cnrow, cncol, cinc1, cinc2) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n %% CPU : %.3f to initialize matrix object",
        t2 - t1) ;
MARKTIME(t1) ;
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
seed++ ;
Drand_setUniform(drand, -1.0, 1.0) ;
DenseMtx_fillRandomEntries(mtxC, drand) ;
MARKTIME(t2) ;
fprintf(msgFile,
      "\n %% CPU : %.3f to fill matrix C with random numbers", t2 - t1) ;

if ( msglvl > 3 ) {
   fprintf(msgFile, "\n matrix A") ;
   DenseMtx_writeForHumanEye(mtxA, msgFile) ;
   fprintf(msgFile, "\n matrix B") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fprintf(msgFile, "\n matrix C") ;
   DenseMtx_writeForHumanEye(mtxC, msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %% beta  = (%f, %f)", beta[0], beta[1]) ;
   fprintf(msgFile, "\n %% alpha = (%f, %f)\n", alpha[0], alpha[1]) ;
   fprintf(msgFile, "\n %% matrix A") ;
   fprintf(msgFile, "\n nrow = %d ;", nrow) ;
   fprintf(msgFile, "\n ncol = %d ;", nk) ;
   DenseMtx_writeForMatlab(mtxA, "A", msgFile) ;
   fprintf(msgFile, "\n");
   fprintf(msgFile, "\n %% matrix B") ;
   fprintf(msgFile, "\n nrow = %d ;", nk) ;
   fprintf(msgFile, "\n ncol = %d ;", ncol) ;
   DenseMtx_writeForMatlab(mtxB, "B", msgFile) ;
   fprintf(msgFile, "\n");
   fprintf(msgFile, "\n %% matrix C") ;
   fprintf(msgFile, "\n nrow = %d ;", cnrow) ;
   fprintf(msgFile, "\n ncol = %d ;", cncol) ;
   DenseMtx_writeForMatlab(mtxC, "C", msgFile) ;
}
/*
   --------------------------
   performs the matrix-matrix operations
   C = alpha*(A)*(B) + beta*C
   --------------------------
*/
   DenseMtx_mmm(A_opt, B_opt, &beta, mtxC, &alpha, mtxA, mtxB);

if ( msglvl > 1 ) {
   fprintf(msgFile, "\n");
   fprintf(msgFile, "\n %% *** Output matrix C ***") ;
   fprintf(msgFile, "\n nrow = %d ;", cnrow) ;
   fprintf(msgFile, "\n ncol = %d ;", cncol) ;
   DenseMtx_writeForMatlab(mtxC, "C", msgFile) ;
   fprintf(msgFile, "\n");
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DenseMtx_free(mtxA) ;
DenseMtx_free(mtxB) ;
DenseMtx_free(mtxC) ;
Drand_free(drand)   ;

return(1) ; }

/*--------------------------------------------------------------------*/
