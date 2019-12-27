/*  PCV.c  */

#include "../Utilities.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to free a pointer to char vector
              must have been created by PCVinit

   created -- 95sep22, cca
   --------------------------------------------
*/
void
PCVfree ( 
   char **p_cvec 
) {
if ( p_cvec != NULL ) {
   FREE(p_cvec) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to allocate and initialize to NULL
              a vector of pointer to char

   created -- 95sep22, cca
   ---------------------------------------------
*/
char **
PCVinit ( 
   int size 
) {
char   **p_cvec = NULL ;
if ( size > 0 ) {
   int   i ;
   ALLOCATE(p_cvec, char *, size) ;
   for ( i = 0 ; i < size ; i++ ) {
      p_cvec[i] = NULL ; 
   }
}
return(p_cvec) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- to set up a pointer vector

   created -- 95sep22, cca
   -------------------------------------
*/
void
PCVsetup ( 
   int    length, 
   int    sizes[], 
   char   cvec[], 
   char   *p_cvec[] 
) {
if ( length > 0 ) {
   if ( sizes == NULL || cvec == NULL || p_cvec == NULL ) {
      fprintf(stderr, "\n fatal error in PCVsetup, invalid data"
              "\n length = %d, sizes = %p, cvec = %p, p_cvec = %p\n",
              length, sizes, cvec, p_cvec) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         if ( sizes[j] > 0 ) {
            p_cvec[j] = cvec ;
            cvec += sizes[j] ;
         } else {
            p_cvec[j] = NULL ;
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
PCVcopy ( 
   int    length, 
   char   *p_cvec1[], 
   char   *p_cvec2[] 
) {
if ( length > 0 ) {
   if ( p_cvec1 == NULL || p_cvec2 == NULL ) {
      fprintf(stdout, "\n fatal error in PCVcopy, invalid data"
              "\n length = %d, p_cvec1 = %p, p_cvec2 = %p\n",
              length, p_cvec1, p_cvec2) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         p_cvec1[j] = p_cvec2[j] ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
