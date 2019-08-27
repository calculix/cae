/*  util.c  */

#include "../Perm.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   return the storage taken by this object

   created -- 96jan05, cca
   ---------------------------------------
*/
int
Perm_sizeOf ( 
   Perm   *perm
) {
int   bytes = sizeof(struct _Perm) ;

if ( perm == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_sizeOf(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}

if ( perm->newToOld != NULL ) {
   bytes += perm->size * sizeof(int) ;
}
if ( perm->oldToNew != NULL ) {
   bytes += perm->size * sizeof(int) ;
}
return(bytes) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   check that the permutation object does house a permutation

   return value --
      1 if a true permutation
      0 otherwise
   ----------------------------------------------------------
*/
int
Perm_checkPerm (
   Perm   *perm
) {
int   inew, iold, rc, size ;
int   *counts, *newToOld, *oldToNew ;
/*
   ---------------
   check the input
   ---------------
*/
if (  perm == NULL 
   || perm->isPresent < 1 || perm->isPresent > 3
   || (size = perm->size) <= 0 ) {
   fprintf(stderr, "\n fatal error in Perm_checkPerm(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}
rc = 1 ;
counts = IVinit(size, 0) ;
if ( (newToOld = perm->newToOld) != NULL ) {
   for ( inew = 0 ; inew < size ; inew++ ) {
      if ( 0 <= (iold = newToOld[inew]) && iold < size ) {
         counts[iold]++ ;
      } else {
         IVfree(counts) ;
         return(0) ;
      }
   }
   for ( iold = 0 ; iold < size ; iold++ ) {
      if ( counts[iold] != 1 ) {
         IVfree(counts) ;
         return(0) ;
      }
   }
}
if ( (oldToNew = perm->oldToNew) != NULL ) {
   IVzero(size, counts) ;
   for ( iold = 0 ; iold < size ; iold++ ) {
      if ( 0 <= (inew = oldToNew[iold]) && inew < size ) {
         counts[inew]++ ;
      } else {
         IVfree(counts) ;
         return(0) ;
      }
   }
   for ( inew = 0 ; inew < size ; inew++ ) {
      if ( counts[inew] != 1 ) {
         IVfree(counts) ;
         return(0) ;
      }
   }
}
IVfree(counts) ;

return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   if the old-to-new vector is not present, 
   create it and fill its entries

   created -- 96mar16, cca
   ----------------------------------------
*/
void
Perm_fillOldToNew (
   Perm   *perm
) {
int   size ;
/*
   ---------------
   check the input
   ---------------
*/
if (  perm == NULL 
   || perm->isPresent < 1 || perm->isPresent > 3
   || (size = perm->size) <= 0 ) {
   fprintf(stderr, "\n fatal error in Perm_fillOldToNew(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}
if ( perm->isPresent == 1 ) {
   int   inew ;
   int   *newToOld = perm->newToOld ;
   int   *oldToNew = perm->oldToNew = IVinit(size, -1) ;

   for ( inew = 0 ; inew < size ; inew++ ) {
      oldToNew[newToOld[inew]] = inew ;
   }
   perm->isPresent = 3 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   if the new-to-old vector is not present, 
   create it and fill its entries

   created -- 96mar16, cca
   ----------------------------------------
*/
void
Perm_fillNewToOld (
   Perm   *perm
) {
int   size ;
/*
   ---------------
   check the input
   ---------------
*/
if (  perm == NULL 
   || perm->isPresent < 1 || perm->isPresent > 3
   || (size = perm->size) <= 0 ) {
   fprintf(stderr, "\n fatal error in Perm_fillNewToOld(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}
if ( perm->isPresent == 2 ) {
   int   iold ;
   int   *newToOld = perm->newToOld = IVinit(size, -1) ;
   int   *oldToNew = perm->oldToNew ;

   for ( iold = 0 ; iold < size ; iold++ ) {
      newToOld[oldToNew[iold]] = iold ;
   }
   perm->isPresent = 3 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   if the old-to-new vector is present, 
   release it and free its entries

   created -- 96mar16, cca
   ------------------------------------
*/
void
Perm_releaseOldToNew (
   Perm   *perm
) {
int   size ;
/*
   ---------------
   check the input
   ---------------
*/
if (  perm == NULL 
   || perm->isPresent < 1 || perm->isPresent > 3
   || (size = perm->size) <= 0 ) {
   fprintf(stderr, "\n fatal error in Perm_fillOldToNew(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}
switch ( perm->isPresent ) {
case 1 :
   break ;
case 2 :
   IVfree(perm->oldToNew) ;
   perm->isPresent = 0 ;
   break ;
case 3 :
   IVfree(perm->oldToNew) ;
   perm->isPresent = 1 ;
   break ;
default :
   break ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   if the new-to-old vector is present, 
   release it and free its entries

   created -- 96mar16, cca
   ------------------------------------
*/
void
Perm_releaseNewToOld (
   Perm   *perm
) {
int   size ;
/*
   ---------------
   check the input
   ---------------
*/
if (  perm == NULL 
   || perm->isPresent < 1 || perm->isPresent > 3
   || (size = perm->size) <= 0 ) {
   fprintf(stderr, "\n fatal error in Perm_fillOldToNew(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}
switch ( perm->isPresent ) {
case 1 :
   IVfree(perm->newToOld) ;
   perm->isPresent = 0 ;
   break ;
case 2 :
   break ;
case 3 :
   IVfree(perm->newToOld) ;
   perm->isPresent = 1 ;
   break ;
default :
   break ;
}

return ; }

/*--------------------------------------------------------------------*/
