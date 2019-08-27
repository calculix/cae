/*  fpgrids.c  */

#include "../misc.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to print a vector on a 2-d grid

   input --

      n1   -- number of points in the first direction
      n2   -- number of points in the second direction
      ivec -- integer vector to be printed in %4d format
      fp   -- file pointer

   created -- 95nov16, cca
   -----------------------------------------------------
*/
void
fp2DGrid ( 
   int    n1, 
   int    n2, 
   int    ivec[], 
   FILE   *fp 
) {
int   i, j ;

if ( n1 <= 0 || n2 <= 0 || ivec == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in fp2DGrid(%d,%d,%p,%p)"
           "\n bad input\n", n1, n2, ivec, fp) ;
   exit(-1) ; 
}
for ( j = n2 - 1 ; j >= 0 ; j-- ) {
   fprintf(fp, "\n") ;
   for ( i = 0 ; i < n1 ; i++ ) {
      fprintf(fp, "%4d", ivec[i + j*n1]) ; 
   } 
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to print a vector on a 3-d grid

   input --

      n1   -- number of points in the first direction
      n2   -- number of points in the second direction
      n3   -- number of points in the third direction
      ivec -- integer vector to be printed in %4d format
      fp   -- file pointer

   created -- 95nov16, cca
   -----------------------------------------------------
*/
void
fp3DGrid ( 
   int    n1, 
   int    n2, 
   int    n3, 
   int    ivec[], 
   FILE   *fp 
) {
int   k ;

if ( n1 <= 0 || n2 <= 0 || n3 <= 0 || ivec == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in fp3DGrid(%d,%d,%d,%p,%p)"
           "\n bad input\n", n1, n2, n3, ivec, fp) ;
   exit(-1) ; 
}
for ( k = 0 ; k < n3 ; k++ ) {
   fprintf(fp, "\n") ;
   fp2DGrid(n1, n2, ivec + k*n1*n2, fp) ; 
}
return ; }

/*--------------------------------------------------------------------*/
