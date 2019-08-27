/*  test_norms.c  */

#include "../A2.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------
   simple test program

   created -- 98apr15, cca
   -----------------------
*/
{
A2       *A ;
double   t1, t2, value ;
FILE     *msgFile ;
int      inc1, inc2, irow, jcol,
         msglvl, nrow, ncol, seed, type ;

if ( argc != 9 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile type nrow ncol inc1 inc2 seed "
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    type    -- entries type"
"\n      1 -- real"
"\n      2 -- complex"
"\n    nrow    -- # of rows "
"\n    ncol    -- # of columns "
"\n    inc1    -- row increment "
"\n    inc2    -- column increment "
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
seed = atoi(argv[7]) ;
fprintf(msgFile, "\n\n %% %s :"
        "\n %% msglvl  = %d"
        "\n %% msgFile = %s"
        "\n %% type    = %d"
        "\n %% nrow    = %d"
        "\n %% ncol    = %d"
        "\n %% inc1    = %d"
        "\n %% inc2    = %d"
        "\n %% seed    = %d"
        "\n",
        argv[0], msglvl, argv[2], type, nrow, ncol, inc1, inc2, seed) ;
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
A2_fillRandomUniform(A, -1, 1, seed) ;
seed++ ;
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
   -------------
   get the norms
   -------------
*/
value = A2_maxabs(A) ;
fprintf(msgFile, "\n error_maxabs = abs(%20.12e - max(max(abs(A))))",
        value) ;
value = A2_frobNorm(A) ;
fprintf(msgFile, "\n error_frob = abs(%20.12e - norm(A, 'fro'))",
        value) ;
value  = A2_oneNorm(A) ;
fprintf(msgFile, "\n error_one = abs(%20.12e - norm(A, 1))",
        value) ;
value  = A2_infinityNorm(A) ;
fprintf(msgFile, "\n error_inf = abs(%20.12e - norm(A, inf))",
        value) ;
for ( irow = 0 ; irow < nrow ; irow++ ) {
   value = A2_infinityNormOfRow(A, irow) ;
   fprintf(msgFile, 
    "\n error_infNormsOfRows(%d) = abs(%20.12e - norm(A(%d,:), inf)) ;",
    irow+1, value, irow+1) ;
   value = A2_oneNormOfRow(A, irow) ;
   fprintf(msgFile, 
    "\n error_oneNormsOfRows(%d) = abs(%20.12e - norm(A(%d,:), 1)) ;",
    irow+1, value, irow+1) ;
   value = A2_twoNormOfRow(A, irow) ;
   fprintf(msgFile, 
    "\n error_twoNormsOfRows(%d) = abs(%20.12e - norm(A(%d,:), 2)) ;",
    irow+1, value, irow+1) ;
}
for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
   value = A2_infinityNormOfColumn(A, jcol) ;
   fprintf(msgFile, 
 "\n error_infNormsOfColumns(%d) = abs(%20.12e - norm(A(:,%d), inf)) ;",
    jcol+1, value, jcol+1) ;
   value = A2_oneNormOfColumn(A, jcol) ;
   fprintf(msgFile, 
   "\n error_oneNormsOfColumns(%d) = abs(%20.12e - norm(A(:,%d), 1)) ;",
      jcol+1, value, jcol+1) ;
   value = A2_twoNormOfColumn(A, jcol) ;
   fprintf(msgFile, 
   "\n error_twoNormsOfColumns(%d) = abs(%20.12e - norm(A(:,%d), 2)) ;",
      jcol+1, value, jcol+1) ;
}
fprintf(msgFile, 
"\n error_in_row_norms = [ max(error_infNormsOfRows) "
"\n                        max(error_oneNormsOfRows) "
"\n                        max(error_twoNormsOfRows) ]"
"\n error_in_column_norms = [ max(error_infNormsOfColumns) "
"\n                           max(error_oneNormsOfColumns) "
"\n                           max(error_twoNormsOfColumns) ]") ;
fprintf(msgFile, "\n") ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
A2_free(A) ;

return(1) ; }

/*--------------------------------------------------------------------*/
