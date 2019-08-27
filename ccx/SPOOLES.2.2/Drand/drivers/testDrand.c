/*  testDrand.c  */

#include "../Drand.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------
   test the Drand random number generator object
   ---------------------------------------------
*/
{
double   ddot, dmean, param1, param2 ;
double   *dvec ;
Drand    drand ;
FILE     *msgFile ;
int      distribution, ierr, imean, msglvl, n, seed1, seed2 ;
int      *ivec ;

if ( argc != 9 ) {
   fprintf(stderr, 
"\n\n usage : testDrand msglvl msgFile "
"\n         distribution param1 param2 seed1 seed2 n"
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    distribution -- 1 for uniform(param1,param2)"
"\n                 -- 2 for normal(param1,param2)"
"\n    param1       -- first parameter"
"\n    param2       -- second parameter"
"\n    seed1        -- first random number seed"
"\n    seed2        -- second random number seed"
"\n    n            -- length of the vector"
"\n"
) ;
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
distribution = atoi(argv[3]) ;
if ( distribution < 1 || distribution > 2 ) {
   fprintf(stderr, "\n fatal error in testDrand"
           "\n distribution must be 1 (uniform) or 2 (normal)") ;
   exit(-1) ;
}
param1 = atof(argv[4]) ;
param2 = atof(argv[5]) ;
seed1  = atoi(argv[6]) ;
seed2  = atoi(argv[7]) ;
n      = atoi(argv[8]) ;

Drand_init(&drand) ;
Drand_setSeeds(&drand, seed1, seed2) ;
switch ( distribution ) {
case 1 : 
   fprintf(msgFile, "\n uniform in [%f,%f]", param1, param2) ;
   Drand_setUniform(&drand, param1, param2) ; 
   break ;
case 2 : 
   fprintf(msgFile, "\n normal(%f,%f)", param1, param2) ;
   Drand_setNormal(&drand, param1, param2) ; 
   break ;
}
/*
   ---------------------------------------------
   fill the integer and double precision vectors
   ---------------------------------------------
*/
dvec = DVinit(n, 0.0) ;
Drand_fillDvector(&drand, n, dvec) ;
dmean = DVsum(n, dvec)/n ;
ddot  = DVdot(n, dvec, dvec) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n dvec mean = %.4f, variance = %.4f",
           dmean, sqrt(fabs(ddot - n*dmean)/n)) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n dvec") ;
   DVfprintf(msgFile, n, dvec) ;
}
DVqsortUp(n, dvec) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n sorted dvec") ;
   DVfprintf(msgFile, n, dvec) ;
}
ivec = IVinit(n, 0) ;
Drand_fillIvector(&drand, n, ivec) ;
imean = IVsum(n, ivec)/n ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n ivec") ;
   IVfp80(msgFile, n, ivec, 80, &ierr) ;
}
IVqsortUp(n, ivec) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n sorted ivec") ;
   IVfp80(msgFile, n, ivec, 80, &ierr) ;
}

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
