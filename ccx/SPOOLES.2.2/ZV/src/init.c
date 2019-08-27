/*  init.C  */

#include "../ZV.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   simplest initialization method

   data is cleared
   if entries != NULL
      the object does not own the entries,
      it just points to the entries base address
   else if size > 0
      the object will own the entries, 
      it allocates a vector of 2*size doubles's.
   else 
      nothing happens
   endif

   created -- 98jan22, cca
   ---------------------------------------------
*/
void
ZV_init (
   ZV       *zv,
   int      size,
   double   *entries 
) {
if ( zv == NULL || size < 0 ) {
   fprintf(stderr, "\n fatal error in ZV_init(%p,%d,%p)"
           "\n bad input\n", zv, size, entries) ;
   exit(-1) ;
}
/*
   --------------
   clear any data
   --------------
*/
ZV_clearData(zv) ;
/*
   -----------------------------
   set the size and maximum size
   -----------------------------
*/
zv->maxsize = zv->size = size ;
/*
   -------------------------
   set vector and owner flag
   -------------------------
*/
if ( entries != NULL ) {
   zv->owned = 0 ;
   zv->vec   = entries ; 
} else if ( size > 0 ) {
   zv->owned = 1 ;
   zv->vec   = DVinit(2*size, 0.0) ;
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------
   basic initializion method
 
   created -- 98jan22, cca
   -------------------------
*/
void
ZV_init1 ( 
   ZV    *zv,
   int   size
) {
ZV_init(zv, size, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------
   total initializion method
 
   created -- 98jan22, cca
   -------------------------
*/
void
ZV_init2 ( 
   ZV       *zv,
   int      size, 
   int      maxsize, 
   int      owned, 
   double   *vec 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_init2(%p,%d,%d,%d,%p)"
           "\n bad input\n", zv, size, maxsize, owned, vec) ;
   exit(-1) ;
}
if ( size < 0 || maxsize < size ) {
   fprintf(stderr, "\n fatal error in ZV_init2(%p,%d,%d,%d,%p)"
           "\n size = %d, maxsize = %d \n", 
           zv, size, maxsize, owned, vec, size, maxsize) ;
   exit(-1) ;
}
if ( owned < 0 || 1 < owned ) {
   fprintf(stderr, "\n fatal error in ZV_init2(%p,%d,%d,%d,%p)"
           "\n owned = %d\n", zv, size, maxsize, owned, vec, owned) ;
   exit(-1) ;
}
if ( owned == 1 && vec == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_init2(%p,%d,%d,%d,%p)"
           "\n owned = %d and vec = %p", 
           zv, size, maxsize, owned, vec, owned, vec) ;
   exit(-1) ;
} 
/*
   --------------
   clear any data
   --------------
*/
ZV_clearData(zv) ;

if ( vec == NULL ) {
/*
   ----------------------------------------------
   no entries input, use the simplest initializer
   ----------------------------------------------
*/
   ZV_init(zv, size, NULL) ;
} else {
/*
   ---------------------------------
   entries are input, set the fields
   ---------------------------------
*/
   zv->size    = size    ;
   zv->maxsize = maxsize ;
   zv->owned   = owned   ;
   zv->vec     = vec     ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   set the maximum size of the vector

   created -- 98jan22, cca
   ----------------------------------
*/
void
ZV_setMaxsize (
   ZV    *zv,
   int   newmaxsize
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || newmaxsize < 0 ) {
   fprintf(stderr, "\n fatal error in ZV_setMaxsize(%p,%d)"
           "\n bad input\n", zv, newmaxsize) ;
   exit(-1) ;
}
if ( zv->maxsize > 0 && zv->owned == 0 ) {
   fprintf(stderr, "\n fatal error in ZV_setMaxsize(%p,%d)"
           "\n zv->maxsize = %d, zv->owned = %d\n", 
           zv, newmaxsize, zv->maxsize, zv->owned) ;
   exit(-1) ;
}
if ( zv->maxsize != newmaxsize ) {
/*
   -----------------------------------
   allocate new storage for the vector
   -----------------------------------
*/
   double   *vec = DVinit(2*newmaxsize, 0.0) ;
   if ( zv->size > 0 ) {
/*
      ---------------------------------
      copy old entries into new entries
      ---------------------------------
*/
      if ( zv->vec == NULL ) {
         fprintf(stderr, "\n fatal error in ZV_setMaxsize(%p,%d)"
                 "\n zv->size = %d, zv->vec is NULL\n", 
                 zv, newmaxsize, zv->size) ;
         exit(-1) ;
      }
      if ( zv->size <= newmaxsize ) {
         DVcopy(2*zv->size, vec, zv->vec) ;
      } else {
/*
         -----------------------
         note, data is truncated
         -----------------------
*/
         DVcopy(2*newmaxsize, vec, zv->vec) ;
         zv->size = newmaxsize ;
      }
   }
   if ( zv->vec != NULL ) {
/*
      ----------------
      free old entries
      ----------------
*/
      DVfree(zv->vec) ;
   }
/*
   ----------
   set fields
   ----------
*/
   zv->maxsize = newmaxsize ;
   zv->owned   = 1 ;
   zv->vec     = vec ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   set the size of the vector

   created -- 98jan22, cca
   --------------------------
*/
void
ZV_setSize (
   ZV    *zv,
   int   newsize
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || newsize < 0 ) {
   fprintf(stderr, "\n fatal error in ZV_setSize(%p,%d)"
           "\n bad input\n", zv, newsize) ;
   exit(-1) ;
}
if ( 0 < zv->maxsize && zv->maxsize < newsize && zv->owned == 0 ) {
   fprintf(stderr, "\n fatal error in ZV_setSize(%p,%d)"
           "\n zv->maxsize = %d, newsize = %d, zv->owned = %d\n", 
           zv, newsize, zv->maxsize, newsize, zv->owned) ;
   exit(-1) ;
}
if ( zv->maxsize < newsize ) {
/*
   -------------------------------------------------------------
   new size requested is more than maxsize, set new maximum size 
   -------------------------------------------------------------
*/
   ZV_setMaxsize(zv, newsize) ;
}
zv->size = newsize ;

return ; }

/*--------------------------------------------------------------------*/
