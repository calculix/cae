/*  init.c  */

#include "../Perm.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   initializer

   created -- 96jan05, cca
   -----------------------
*/
void
Perm_initWithTypeAndSize ( 
   Perm   *perm,
   int    isPresent, 
   int    size 
) {
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Perm_clearData(perm) ;
/*
   ---------------
   check the input
   ---------------
*/
if ( isPresent < 1 || isPresent > 3 || size <= 0 ) {
   fprintf(stderr, 
           "\n\n fatal error in Perm_initWithTypeAndSize(%p,%d,%d)"
           "\n isPresent = %d, must be 1, 2 or 3"
           "\n size = %d, must be positive", 
           perm, isPresent, size, isPresent, size) ;
   exit(-1) ;
} 
perm->isPresent = isPresent ;
perm->size      = size      ;
/*
   ----------------------------------------
   allocate storage for permutation vectors
   ----------------------------------------
*/
switch ( isPresent ) {
case 1 :
   perm->newToOld = IVinit(size, -1) ;
   break ;
case 2 :
   perm->oldToNew = IVinit(size, -1) ;
   break ;
case 3 :
   perm->oldToNew = IVinit(size, -1) ;
   perm->newToOld = IVinit(size, -1) ;
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
