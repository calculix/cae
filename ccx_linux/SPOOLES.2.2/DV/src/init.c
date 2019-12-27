/*  init.C  */

#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   simplest initialization method

   if entries != NULL
      the object does not own the entries,
      it just points to the entries base address
   else if size > 0
      the object will own the entries, 
      it allocates a vector of size doubles's.
   else 
      nothing happens
   endif

   created -- 96aug28, cca
   ---------------------------------------------
*/
void
DV_init (
   DV       *dv,
   int      size,
   double   *entries 
) {
if ( dv == NULL || size < 0 ) {
   fprintf(stderr, "\n fatal error in DV_init(%p,%d,%p)"
           "\n bad input\n", dv, size, entries) ;
   exit(-1) ;
}
/*
   --------------
   clear any data
   --------------
*/
DV_clearData(dv) ;
/*
   -----------------------------
   set the size and maximum size
   -----------------------------
*/
dv->maxsize = dv->size = size ;
/*
   -------------------------
   set vector and owner flag
   -------------------------
*/
if ( entries != NULL ) {
   dv->owned = 0 ;
   dv->vec   = entries ; 
} else if ( size > 0 ) {
   dv->owned = 1 ;
/*
   dv->vec   = DVinit(size, 0.0) ;
*/
   dv->vec   = DVinit2(size) ;
}
return ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------
   basic initializion method
 
   created -- 95oct06, cca
   -------------------------
*/
void
DV_init1 ( 
   DV    *dv,
   int   size
) {
DV_init(dv, size, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------
   total initializion method
 
   created -- 95oct06, cca
   -------------------------
*/
void
DV_init2 ( 
   DV       *dv,
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
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_init2(%p,%d,%d,%d,%p)"
           "\n bad input\n", dv, size, maxsize, owned, vec) ;
   exit(-1) ;
}
if ( size < 0 || maxsize < size ) {
   fprintf(stderr, "\n fatal error in DV_init2(%p,%d,%d,%d,%p)"
           "\n size = %d, maxsize = %d \n", 
           dv, size, maxsize, owned, vec, size, maxsize) ;
   exit(-1) ;
}
if ( owned < 0 || 1 < owned ) {
   fprintf(stderr, "\n fatal error in DV_init2(%p,%d,%d,%d,%p)"
           "\n owned = %d\n", dv, size, maxsize, owned, vec, owned) ;
   exit(-1) ;
}
if ( owned == 1 && vec == NULL ) {
   fprintf(stderr, "\n fatal error in DV_init2(%p,%d,%d,%d,%p)"
           "\n owned = %d and vec = %p", 
           dv, size, maxsize, owned, vec, owned, vec) ;
   exit(-1) ;
} 
/*
   --------------
   clear any data
   --------------
*/
DV_clearData(dv) ;

if ( vec == NULL ) {
/*
   ----------------------------------------------
   no entries input, use the simplest initializer
   ----------------------------------------------
*/
   DV_init(dv, size, NULL) ;
} else {
/*
   ---------------------------------
   entries are input, set the fields
   ---------------------------------
*/
   dv->size    = size    ;
   dv->maxsize = maxsize ;
   dv->owned   = owned   ;
   dv->vec     = vec     ;
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
DV_setMaxsize (
   DV    *dv,
   int   newmaxsize
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || newmaxsize < 0 ) {
   fprintf(stderr, "\n fatal error in DV_setMaxsize(%p,%d)"
           "\n bad input\n", dv, newmaxsize) ;
   exit(-1) ;
}
if ( dv->maxsize > 0 && dv->owned == 0 ) {
   fprintf(stderr, "\n fatal error in DV_setMaxsize(%p,%d)"
           "\n dv->maxsize = %d, dv->owned = %d\n", 
           dv, newmaxsize, dv->maxsize, dv->owned) ;
   exit(-1) ;
}
if ( dv->maxsize != newmaxsize ) {
/*
   -----------------------------------
   allocate new storage for the vector
   -----------------------------------
*/
   double   *vec ;
/*
   vec = DVinit(newmaxsize, 0.0) ;
*/
   vec = DVinit2(newmaxsize) ;
   if ( dv->size > 0 ) {
/*
      ---------------------------------
      copy old entries into new entries
      ---------------------------------
*/
      if ( dv->vec == NULL ) {
         fprintf(stderr, "\n fatal error in DV_setMaxsize(%p,%d)"
                 "\n dv->size = %d, dv->vec is NULL\n", 
                 dv, newmaxsize, dv->size) ;
         exit(-1) ;
      }
      if ( dv->size <= newmaxsize ) {
         DVcopy(dv->size, vec, dv->vec) ;
      } else {
/*
         -----------------------
         note, data is truncated
         -----------------------
*/
         DVcopy(newmaxsize, vec, dv->vec) ;
         dv->size = newmaxsize ;
      }
   }
   if ( dv->vec != NULL ) {
/*
      ----------------
      free old entries
      ----------------
*/
      DVfree(dv->vec) ;
   }
/*
   ----------
   set fields
   ----------
*/
   dv->maxsize = newmaxsize ;
   dv->owned   = 1 ;
   dv->vec     = vec ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   set the size of the vector

   created -- 96dec08, cca
   --------------------------
*/
void
DV_setSize (
   DV    *dv,
   int   newsize
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || newsize < 0 ) {
   fprintf(stderr, "\n fatal error in DV_setSize(%p,%d)"
           "\n bad input\n", dv, newsize) ;
   exit(-1) ;
}
if ( 0 < dv->maxsize && dv->maxsize < newsize && dv->owned == 0 ) {
   fprintf(stderr, "\n fatal error in DV_setSize(%p,%d)"
           "\n dv->maxsize = %d, newsize = %d, dv->owned = %d\n", 
           dv, newsize, dv->maxsize, newsize, dv->owned) ;
   exit(-1) ;
}
if ( dv->maxsize < newsize ) {
/*
   -------------------------------------------------------------
   new size requested is more than maxsize, set new maximum size 
   -------------------------------------------------------------
*/
   DV_setMaxsize(dv, newsize) ;
}
dv->size = newsize ;

return ; }

/*--------------------------------------------------------------------*/
