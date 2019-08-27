/*  PFV.c  */

#include "../Utilities.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to free a pointer to float vector
              must have been created by PFVinit

   created -- 95sep22, cca
   --------------------------------------------
*/
void
PFVfree ( 
   float **p_fvec 
) {
if ( p_fvec != NULL ) {
   FREE(p_fvec) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to allocate and initialize to NULL
              a vector of pointer to float

   created -- 95sep22, cca
   ---------------------------------------------
*/
float **
PFVinit ( 
   int size 
) {
float   **p_fvec = NULL ;
if ( size > 0 ) {
   int   i ;
   ALLOCATE(p_fvec, float *, size) ;
   for ( i = 0 ; i < size ; i++ ) {
      p_fvec[i] = NULL ; 
   }
}
return(p_fvec) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- to set up a pointer vector

   created -- 95sep22, cca
   -------------------------------------
*/
void
PFVsetup ( 
   int     length, 
   int     sizes[], 
   float   fvec[], 
   float   *p_fvec[] 
) {
if ( length > 0 ) {
   if ( sizes == NULL || fvec == NULL || p_fvec == NULL ) {
      fprintf(stderr, "\n fatal error in PFVsetup, invalid data"
              "\n length = %d, sizes = %p, fvec = %p, p_fvec = %p\n",
              length, sizes, fvec, p_fvec) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         if ( sizes[j] > 0 ) {
            p_fvec[j] = fvec ;
            fvec += sizes[j] ;
         } else {
            p_fvec[j] = NULL ;
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
PFVcopy ( 
   int     length, 
   float   *p_fvec1[], 
   float   *p_fvec2[] 
) {
if ( length > 0 ) {
   if ( p_fvec1 == NULL || p_fvec2 == NULL ) {
      fprintf(stdout, "\n fatal error in PFVcopy, invalid data"
              "\n length = %d, p_fvec1 = %p, p_fvec2 = %p\n",
              length, p_fvec1, p_fvec2) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         p_fvec1[j] = p_fvec2[j] ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
