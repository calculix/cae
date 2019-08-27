/*  PDV.c  */

#include "../Utilities.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to free a pointer to double vector
              must have been created by PDVinit

   created -- 95sep22, cca
   ---------------------------------------------
*/
void
PDVfree ( 
   double **p_dvec 
) {
if ( p_dvec != NULL ) {
   FREE(p_dvec) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to allocate and initialize to NULL
              a vector of pointer to double

   created -- 95sep22, cca
   ---------------------------------------------
*/
double **
PDVinit ( 
   int size 
) {
double   **p_dvec = NULL ;
if ( size > 0 ) {
   int   i ;
   ALLOCATE(p_dvec, double *, size) ;
   for ( i = 0 ; i < size ; i++ ) {
      p_dvec[i] = NULL ; 
   }
}
return(p_dvec) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- to set up a pointer vector

   created -- 95sep22, cca
   -------------------------------------
*/
void
PDVsetup ( 
   int      length, 
   int      sizes[], 
   double   dvec[], 
   double   *p_dvec[] 
) {
if ( length > 0 ) {
   if ( sizes == NULL || dvec == NULL || p_dvec == NULL ) {
      fprintf(stderr, "\n fatal error in PDVsetup, invalid data"
              "\n length = %d, sizes = %p, dvec = %p, p_dvec = %p\n",
              length, sizes, dvec, p_dvec) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         if ( sizes[j] > 0 ) {
            p_dvec[j] = dvec ;
            dvec += sizes[j] ;
         } else {
            p_dvec[j] = NULL ;
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
PDVcopy ( 
   int      length, 
   double   *p_dvec1[], 
   double   *p_dvec2[] 
) {
if ( length > 0 ) {
   if ( p_dvec1 == NULL || p_dvec2 == NULL ) {
      fprintf(stdout, "\n fatal error in PDVcopy, invalid data"
              "\n length = %d, p_dvec1 = %p, p_dvec2 = %p\n",
              length, p_dvec1, p_dvec2) ;
      exit(-1) ;
   } else {
      int   j ;
      for ( j = 0 ; j < length ; j++ ) {
         p_dvec1[j] = p_dvec2[j] ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
