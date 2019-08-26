/*  init.c  */

#include "../Drand.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   initialize the Drand object

   created -- 96may26, cca
   ---------------------------
*/
void
Drand_init (
   Drand   *drand
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_init(%p)"
           "\n bad input\n", drand) ;
   exit(-1) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
Drand_setDefaultFields(drand) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   set the seeds using one input seed

   created -- 96may26, cca
   ----------------------------------
*/
void
Drand_setSeed (
   Drand   *drand,
   int     u
) {
if ( drand == NULL || u <= 0 || u >= drand->base1 ) {
   fprintf(stderr, "\n fatal error in Drand_setSeed(%p,%d)"
           "\n first seed must in in (0,%.0f)", 
           drand, u, drand->base1) ;
   exit(-1) ;
}
drand->seed1 = u ;
drand->seed2 = fmod(2718.*u, drand->base2) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   set the seeds using two input seeds

   created -- 96may26, cca
   -----------------------------------
*/
void
Drand_setSeeds (
   Drand   *drand,
   int     u,
   int     v
) {
if (  drand == NULL
   || u <= 0 || u >= drand->base1 
   || v <= 0 || v >= drand->base2 ) {
   fprintf(stderr, "\n fatal error in Drand_setSeeds(%p,%d,%d)"
           "\n first seed must in in (0,%.0f)"
           "\n second seed must in in (0,%.0f)\n",
           drand, u, v, drand->base1, drand->base2) ;
   exit(-1) ;
}
drand->seed1 = u ;
drand->seed2 = v ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   set the mode to be uniform in [lower, upper]

   created -- 96may26, cca
   --------------------------------------------
*/
void
Drand_setUniform (
   Drand    *drand,
   double   lower,
   double   upper
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL || lower > upper ) {
   fprintf(stderr, "\n fatal error in Drand_setUniform(%p,%f,%f)"
           "\n bad input\n", drand, lower, upper) ;
   exit(-1) ;
}
drand->mode  =   1   ;
drand->lower = lower ;
drand->upper = upper ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   set the mode to be normal(mean, sigma)

   created -- 96may26, cca
   --------------------------------------
*/
void
Drand_setNormal (
   Drand    *drand,
   double   mean,
   double   sigma
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL || sigma <= 0 ) {
   fprintf(stderr, "\n fatal error in Drand_setNormal(%p,%f,%f)"
           "\n bad input\n", drand, mean, sigma) ;
   exit(-1) ;
}
drand->mode  =   2   ;
drand->mean  = mean  ;
drand->sigma = sigma ;

return ; }

/*--------------------------------------------------------------------*/
