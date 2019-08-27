/*  util.c  */

#include "../Drand.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   return a random double precision value

   created -- 96may26, cca
   --------------------------------------
*/
double
Drand_value (
   Drand   *drand
) {
double   sum, t ;
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_value(%p)"
           "\n bad input\n", drand) ;
   exit(-1) ;
}
/*
   --------------------
   switch over the mode
   --------------------
*/
if ( drand->mode == 1 ) {
/*
   ------------
   uniform mode
   ------------
*/
   drand->seed1 = fmod(40014*drand->seed1, drand->base1) ;
   drand->seed2 = fmod(40692*drand->seed2, drand->base2) ;
   t = drand->seed1 - drand->seed2 ;
   if ( t <= 0 ) {
      t = t + (drand->base1 - 1) ;
   }
   t = drand->lower + (t/drand->base1)*(drand->upper - drand->lower) ;
} else {
/*
   -----------
   normal mode
   -----------
*/
   drand->seed1 = fmod(40014*drand->seed1, drand->base1) ;
   drand->seed2 = fmod(40692*drand->seed2, drand->base2) ;
   t = drand->seed1 - drand->seed2 ;
   if ( t <= 0 ) {
      t = t + (drand->base1 - 1) ;
   }
   t = t / drand->base1 ;
   sum = t ;
   drand->seed1 = fmod(40014*drand->seed1, drand->base1) ;
   drand->seed2 = fmod(40692*drand->seed2, drand->base2) ;
   t = drand->seed1 - drand->seed2 ;
   if ( t <= 0 ) {
      t = t + (drand->base1 - 1) ;
   }
   t = t / drand->base1 ;
   sum += t ;
   drand->seed1 = fmod(40014*drand->seed1, drand->base1) ;
   drand->seed2 = fmod(40692*drand->seed2, drand->base2) ;
   t = drand->seed1 - drand->seed2 ;
   if ( t <= 0 ) {
      t = t + (drand->base1 - 1) ;
   }
   t = t / drand->base1 ;
   sum += t ;
   t = drand->mean + drand->sigma*(2.*sum - 3.) ;
}

return(t) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   fill a double precision complex vector with random numbers

   created -- 98jun02, cca
   ----------------------------------------------------------
*/
void
Drand_fillZvector (
   Drand    *drand,
   int      size,
   double   dvec[] 
) {
int   i ;
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL || size < 0 || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_fillZvector(%p,%d,%p)"
           "\n bad input\n", drand, size, dvec) ;
   exit(-1) ;
}
/*
   ---------------
   fill the vector
   ---------------
*/
for ( i = 0 ; i < 2*size ; i++ ) {
   dvec[i] = Drand_value(drand) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   fill a double precision vector with random numbers

   created -- 96may26, cca
   --------------------------------------------------
*/
void
Drand_fillDvector (
   Drand    *drand,
   int      size,
   double   dvec[] 
) {
int   i ;
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL || size < 0 || dvec == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_fillDvector(%p,%d,%p)"
           "\n bad input\n", drand, size, dvec) ;
   exit(-1) ;
}
/*
   ---------------
   fill the vector
   ---------------
*/
for ( i = 0 ; i < size ; i++ ) {
   dvec[i] = Drand_value(drand) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   fill a integer vector with random numbers

   created -- 96may26, cca
   -----------------------------------------
*/
void
Drand_fillIvector (
   Drand    *drand,
   int      size,
   int      ivec[] 
) {
int   i ;
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL || size < 0 || ivec == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_fillIvector(%p,%d,%p)"
           "\n bad input\n", drand, size, ivec) ;
   exit(-1) ;
}
/*
   ---------------
   fill the vector
   ---------------
*/
for ( i = 0 ; i < size ; i++ ) {
   ivec[i] = (int) Drand_value(drand) ;
}

return ; }

/*--------------------------------------------------------------------*/
