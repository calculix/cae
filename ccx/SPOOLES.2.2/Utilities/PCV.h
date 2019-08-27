/*  PCV.h  */

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
) ;
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
) ;
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
