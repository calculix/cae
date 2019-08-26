/*  PIV.c  */

#include "../Utilities.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to free a pointer to int vector
              must have been created by PIVinit

   created -- 95sep22, cca
   --------------------------------------------
*/
void
PIVfree ( 
   int **p_ivec 
) {
if ( p_ivec != NULL ) {
   FREE(p_ivec) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to allocate and initialize to NULL
              a vector of pointer to int

   created -- 95sep22, cca
   ---------------------------------------------
*/
int **
PIVinit ( 
   int size 
) {
int   **p_ivec = NULL ;
if ( size > 0 ) {
   int   i ;
   ALLOCATE(p_ivec, int *, size) ;
   for ( i = 0 ; i < size ; i++ ) {
      p_ivec[i] = NULL ; 
   }
}
return(p_ivec) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- to set up a pointer vector

   created -- 95sep22, cca
   -------------------------------------
*/
void
PIVsetup ( 
   int   length, 
   int   sizes[], 
   int   ivec[], 
   int   *p_ivec[] 
) {
if ( length > 0 ) {
   if ( sizes == NULL || ivec == NULL || p_ivec == NULL ) {
      fprintf(stderr, "\n fatal error in PIVsetup, invalid data"
              "\n length = %d, sizes = %p, ivec = %p, p_ivec = %p\n",
              length, sizes, ivec, p_ivec) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         if ( sizes[j] > 0 ) {
            p_ivec[j] = ivec ;
            ivec += sizes[j] ;
         } else {
            p_ivec[j] = NULL ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- to copy a pointer vector

   created -- 95sep22, cca
   -----------------------------------
*/
void
PIVcopy ( 
   int   length, 
   int   *p_ivec1[], 
   int   *p_ivec2[] 
) {
if ( length > 0 ) {
   if ( p_ivec1 == NULL || p_ivec2 == NULL ) {
      fprintf(stdout, "\n fatal error in PIVcopy, invalid data"
              "\n length = %d, p_ivec1 = %p, p_ivec2 = %p\n",
              length, p_ivec1, p_ivec2) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         p_ivec1[j] = p_ivec2[j] ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
