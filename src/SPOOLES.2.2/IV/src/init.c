/*  init.C  */

#include "../IV.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   simplest initialization method

   if entries != NULL
      the object does not own the entries,
      it just points to the entries base address
   else if size > 0
      the object will own the entries, 
      it allocates a vector of size int's.
   else 
      nothing happens
   endif

   created -- 96aug28, cca
   ---------------------------------------------
*/
void
IV_init (
   IV    *iv,
   int   size,
   int   *entries 
) {
if ( iv == NULL || size < 0 ) {
   fprintf(stderr, "\n fatal error in IV_init(%p,%d,%p)"
           "\n bad input\n", iv, size, entries) ;
   exit(-1) ;
}
/*
   --------------
   clear any data
   --------------
*/
IV_clearData(iv) ;
/*
   -----------------------------
   set the size and maximum size
   -----------------------------
*/
iv->maxsize = iv->size = size ;
/*
   -------------------------
   set vector and owner flag
   -------------------------
*/
if ( entries != NULL ) {
   iv->owned = 0 ;
   iv->vec   = entries ; 
} else if ( size > 0 ) {
   iv->owned = 1 ;
   iv->vec   = IVinit(size, -1) ;
}
/*
fprintf(stdout, 
        "\n %% leaving IV_init, iv %p, size %d, maxsize %d, entries %p",
        iv, iv->size, iv->maxsize, iv->vec) ;
fflush(stdout) ;
*/

return ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------
   basic initializion method
 
   created -- 95oct06, cca
   -------------------------
*/
void
IV_init1 ( 
   IV    *iv,
   int   size
) {
IV_init(iv, size, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------
   total initializion method
 
   created -- 95oct06, cca
   -------------------------
*/
void
IV_init2 ( 
   IV    *iv,
   int   size, 
   int   maxsize, 
   int   owned, 
   int   *vec 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_init2(%p,%d,%d,%d,%p)"
           "\n bad input\n", iv, size, maxsize, owned, vec) ;
   exit(-1) ;
}
if ( size < 0 || maxsize < size ) {
   fprintf(stderr, "\n fatal error in IV_init2(%p,%d,%d,%d,%p)"
           "\n size = %d, maxsize = %d \n", 
           iv, size, maxsize, owned, vec, size, maxsize) ;
   exit(-1) ;
}
if ( owned < 0 || 1 < owned ) {
   fprintf(stderr, "\n fatal error in IV_init2(%p,%d,%d,%d,%p)"
           "\n owned = %d\n", iv, size, maxsize, owned, vec, owned) ;
   exit(-1) ;
}
if ( owned == 1 && vec == NULL ) {
   fprintf(stderr, "\n fatal error in IV_init2(%p,%d,%d,%d,%p)"
           "\n owned = %d and vec = %p", 
           iv, size, maxsize, owned, vec, owned, vec) ;
   exit(-1) ;
} 
/*
   --------------
   clear any data
   --------------
*/
IV_clearData(iv) ;

if ( vec == NULL ) {
/*
   ----------------------------------------------
   no entries input, use the simplest initializer
   ----------------------------------------------
*/
   IV_init(iv, size, NULL) ;
} else {
/*
   ---------------------------------
   entries are input, set the fields
   ---------------------------------
*/
   iv->size    = size    ;
   iv->maxsize = maxsize ;
   iv->owned   = owned   ;
   iv->vec     = vec     ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   set the maximum size of the vector

   created -- 96dec08, cca
   ----------------------------------
*/
void
IV_setMaxsize (
   IV    *iv,
   int   newmaxsize
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || newmaxsize < 0 ) {
   fprintf(stderr, "\n fatal error in IV_setMaxsize(%p,%d)"
           "\n bad input\n", iv, newmaxsize) ;
   exit(-1) ;
}
if ( iv->maxsize > 0 && iv->owned == 0 ) {
   fprintf(stderr, "\n fatal error in IV_setMaxsize(%p,%d)"
           "\n iv->maxsize = %d, iv->owned = %d\n", 
           iv, newmaxsize, iv->maxsize, iv->owned) ;
   exit(-1) ;
}
if ( iv->maxsize != newmaxsize ) {
/*
   -----------------------------------
   allocate new storage for the vector
   -----------------------------------
*/
   int   *vec = IVinit(newmaxsize, -1) ;
   if ( iv->size > 0 ) {
/*
      ---------------------------------
      copy old entries into new entries
      ---------------------------------
*/
      if ( iv->vec == NULL ) {
         fprintf(stderr, "\n fatal error in IV_setMaxsize(%p,%d)"
                 "\n iv->size = %d, iv->vec is NULL\n", 
                 iv, newmaxsize, iv->size) ;
         exit(-1) ;
      }
      if ( iv->size <= newmaxsize ) {
/*
         -----------------------------------------
         new maximum size is greater than old size
         -----------------------------------------
*/
         IVcopy(iv->size, vec, iv->vec) ;
      } else {
/*
         -----------------------
         note, data is truncated
         -----------------------
*/
         IVcopy(newmaxsize, vec, iv->vec) ;
         iv->size = newmaxsize ;
      }
   }
   if ( iv->vec != NULL ) {
/*
      ----------------
      free old entries
      ----------------
*/
      IVfree(iv->vec) ;
   }
/*
   ----------
   set fields
   ----------
*/
   iv->maxsize = newmaxsize ;
   iv->owned   = 1 ;
   iv->vec     = vec ;
}
/*
fprintf(stdout, 
        "\n %% leaving IV_setMaxsize, iv %p, size %d, maxsize %d, entries %p",
        iv, iv->size, iv->maxsize, iv->vec) ;
fflush(stdout) ;
*/
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   set the size of the vector

   created -- 96dec08, cca
   --------------------------
*/
void
IV_setSize (
   IV    *iv,
   int   newsize
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || newsize < 0 ) {
   fprintf(stderr, "\n fatal error in IV_setSize(%p,%d)"
           "\n bad input\n", iv, newsize) ;
   exit(-1) ;
}
if ( 0 < iv->maxsize && iv->maxsize < newsize && iv->owned == 0 ) {
   fprintf(stderr, "\n fatal error in IV_setSize(%p,%d)"
           "\n iv->maxsize = %d, newsize = %d, iv->owned = %d\n", 
           iv, newsize, iv->maxsize, newsize, iv->owned) ;
   exit(-1) ;
}
if ( iv->maxsize < newsize ) {
/*
   -------------------------------------------------------------
   new size requested is more than maxsize, set new maximum size 
   -------------------------------------------------------------
*/
   IV_setMaxsize(iv, newsize) ;
}
iv->size = newsize ;
/*
fprintf(stdout, 
        "\n %% leaving IV_setSize, iv %p, size %d, maxsize %d, entries %p",
        iv, iv->size, iv->maxsize, iv->vec) ;
fflush(stdout) ;
*/

return ; }

/*--------------------------------------------------------------------*/
