/*  PIV.h  */

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
) ;
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
) ;
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
