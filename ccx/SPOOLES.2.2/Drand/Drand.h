/*  Drand.h  */

#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   Double precision random number generator object

   there are two modes --- 
   uniform on [lower,upper] and normal(mean, sigma)

   created -- 96may26, cca
   ------------------------------------------------
*/
typedef 
struct _Drand {
   double   seed1 ;
   double   seed2 ;
   double   base1 ;
   double   base2 ;
   double   lower ;
   double   upper ;
   double   mean  ;
   double   sigma ;
   int      mode  ;
} Drand ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor

   created -- 96may26, cca
   -----------------------
*/
Drand *
Drand_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields

   created -- 96may26, cca
   -----------------------
*/
void
Drand_setDefaultFields (
   Drand   *drand
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 96may26, cca
   --------------------------------------------------
*/
void
Drand_clearData ( 
   Drand   *drand 
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 96may26, cca
   ------------------------------------------
*/
Drand *
Drand_free ( 
   Drand   *drand 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------
   initialize the Drand object

   created -- 96may26, cca
   ---------------------------
*/
void
Drand_init (
   Drand   *drand
) ;
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
) ;
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
) ;
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------
   return a random double precision value

   created -- 96may26, cca
   --------------------------------------
*/
double
Drand_value (
   Drand   *drand
) ;
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
) ;
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
) ;
/*
   -----------------------------------------
   fill a integer vector with random numbers

   created -- 96may26, cca
   -----------------------------------------
*/
void
Drand_fillIvector (
   Drand   *drand,
   int     size,
   int     ivec[] 
) ;
/*--------------------------------------------------------------------*/
