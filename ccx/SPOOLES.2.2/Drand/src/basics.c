/*  basics.c  */

#include "../Drand.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 96may26, cca
   -----------------------
*/
Drand *
Drand_new ( 
   void 
) {
Drand   *drand ;

ALLOCATE(drand, struct _Drand, 1) ;
Drand_setDefaultFields(drand) ;

return(drand) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 96may26, cca
   -----------------------
*/
void
Drand_setDefaultFields (
   Drand   *drand
) {
if ( drand == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_setDefaultFields(%p)"
           "\n bad input", drand) ;
   exit(-1) ;
}
drand->seed1 =  123456789.0 ;
drand->seed2 =  987654321.0 ;
drand->base1 = 2147483563.0 ;
drand->base2 = 2147483399.0 ;
drand->lower = 0.0 ;
drand->upper = 1.0 ;
drand->mean  = 0.0 ;
drand->sigma = 1.0 ;
drand->mode  =  1  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 96may26, cca
   --------------------------------------------------
*/
void
Drand_clearData ( 
   Drand   *drand 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( drand == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_clearData(%p)"
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
   ------------------------------------------
   destructor, free's the object and its data

   created -- 96may26, cca
   ------------------------------------------
*/
Drand *
Drand_free ( 
   Drand   *drand 
) {
if ( drand == NULL ) {
   fprintf(stderr, "\n fatal error in Drand_free(%p)"
           "\n bad input\n", drand) ;
   exit(-1) ;
}
Drand_clearData(drand) ;
FREE(drand) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
