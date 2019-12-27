/*  PDV.h  */

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
   double   **p_dvec 
) ;
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
   int   size 
) ;
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
