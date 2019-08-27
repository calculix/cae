/*  instance.c  */

#include "../IV.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   return value is 0 if the entries are not owned by the object
   otherwise the return value is the number of entries at the base
   storage of the vector

   created  -- 96jun22, cca
   modified -- 96aug28, cca
   ---------------------------------------------------------------
*/
int
IV_owned (
   IV   *iv
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_owned(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
return(iv->owned) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   return the vector size

   created -- 95oct06, cca
   -----------------------
*/
int
IV_size (
   IV   *iv
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_size(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
return(iv->size) ; }
/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   return the maximum size of the vector

   created -- 95dec08, cca
   -------------------------------------
*/
int
IV_maxsize (
   IV   *iv
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_maxsize(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
return(iv->maxsize) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   return the loc'th entry of a vector.
   note: if loc is out of range then -1 is returned

   created -- 96jun29, cca
   ------------------------------------------------
*/
int 
IV_entry (
   IV    *iv,
   int   loc
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_entries(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
if ( loc < 0 || loc >= iv->size ) {
   return(-1) ;
} else {
   return(iv->vec[loc]) ; 
}
}

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return a pointer to the object's entries array

   created -- 95oct06, cca
   ----------------------------------------------
*/
int *
IV_entries (
   IV   *iv
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_entries(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
return(iv->vec) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   fill *psize with the vector's size
   and *pentries with the address of the vector

   created -- 95oct06, cca
   --------------------------------------------
*/
void
IV_sizeAndEntries (
   IV    *iv,
   int   *psize,
   int   **pentries
) {
if ( iv == NULL || psize == NULL || pentries == NULL ) {
   fprintf(stderr, "\n fatal error in IV_sizeAndEntries(%p,%p,%p)"
           "\n bad input\n", iv, psize, pentries) ;
   exit(-1) ;
}
*psize    = iv->size ;
*pentries = iv->vec  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   set and entry in the vector

   created -- 96jul14, cca
   ---------------------------
*/
void
IV_setEntry ( 
   IV    *iv,
   int   loc,
   int   value
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || loc < 0 ) {
   fprintf(stderr, "\n fatal error in IV_setEntry(%p,%d,%d)"
           "\n bad input\n", iv, loc, value) ;
   exit(-1) ;
}
if ( loc >= iv->maxsize ) {
   int newmaxsize = (int) 1.25*iv->maxsize ;
   if ( newmaxsize < 10 ) {
      newmaxsize = 10 ;
   }
   if ( loc >= newmaxsize ) {
      newmaxsize = loc + 1 ;
   }
   IV_setMaxsize(iv, newmaxsize) ;
}
if ( loc >= iv->size ) {
   iv->size = loc + 1 ;
}
iv->vec[loc] = value ;

return ; }
   
/*--------------------------------------------------------------------*/
