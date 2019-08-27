/*  PFV.h  */

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to free a pointer to float vector
              must have been created by PFVinit

   created -- 95sep22, cca
   ---------------------------------------------
*/
void
PFVfree ( 
   float   **p_fvec 
) ;
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
PFVsetup ( 
   int     length, 
   int     sizes[], 
   float   fvec[], 
   float   *p_fvec[] 
) ;
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
) ;
/*--------------------------------------------------------------------*/
