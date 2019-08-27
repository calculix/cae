/*  instance.c  */

#include "../ZV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   return 1 if the entries are owned by the object
   return 0 otherwise

   created -- 98jan22, cca
   -----------------------------------------------
*/
int
ZV_owned (
   ZV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_owned(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->owned) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   return the vector size

   created -- 98jan22, cca
   -----------------------
*/
int
ZV_maxsize (
   ZV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_maxsize(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->maxsize) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   return the vector size

   created -- 98jan22, cca
   -----------------------
*/
int
ZV_size (
   ZV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_size(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->size) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   return the loc'th entry of a vector.
 
   created -- 98jan22, cca
   -------------------------------------------------
*/
void
ZV_entry (
   ZV       *dv,
   int      loc,
   double   *pReal,
   double   *pImag
) {
if ( dv == NULL || pReal == NULL || pImag == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_entry(%p,%d,%p,%p)"
           "\n bad input\n", dv, loc, pReal, pImag) ;
   exit(-1) ;
}
if ( loc < 0 || loc >= dv->size || dv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_entry(%p,%d,%p,%p)"
           "\n bad state: size = %d, vec = %p\n",
           dv, loc, pReal, pImag, dv->size, dv->vec) ;
   exit(-1) ;
}
*pReal = dv->vec[2*loc] ;
*pImag = dv->vec[2*loc+1] ;

return ; }
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   return pointers to the loc'th entry of a vector.
 
   created -- 98jan22, cca
   -------------------------------------------------
*/
void
ZV_pointersToEntry (
   ZV       *dv,
   int      loc,
   double   **ppReal,
   double   **ppImag
) {
if ( dv == NULL || ppReal == NULL || ppImag == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_pointersToEntry(%p,%d,%p,%p)"
           "\n bad input\n", dv, loc, ppReal, ppImag) ;
   exit(-1) ;
}
if ( loc < 0 || loc >= dv->size || dv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_pointersToEntry(%p,%d,%p,%p)"
           "\n bad state: size = %d, vec = %p\n",
           dv, loc, ppReal, ppImag, dv->size, dv->vec) ;
   exit(-1) ;
}
*ppReal = &dv->vec[2*loc] ;
*ppImag = &dv->vec[2*loc+1] ;

return ; }
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return a pointer to the object's entries array

   created -- 98jan22, cca
   ----------------------------------------------
*/
double *
ZV_entries (
   ZV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_entries(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->vec) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   fill *psize with the vector's size
   and *pentries with the address of the vector

   created -- 98jan22, cca
   --------------------------------------------
*/
void
ZV_sizeAndEntries (
   ZV       *dv,
   int      *psize,
   double   **pentries
) {
if ( dv == NULL || psize == NULL || pentries == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_sizeAndEntries(%p,%p,%p)"
           "\n bad input\n", dv, psize, pentries) ;
   exit(-1) ;
}
*psize    = dv->size ;
*pentries = dv->vec  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   set and entry in the vector
 
   created -- 98jan22, cca
   ---------------------------
*/
void
ZV_setEntry (
   ZV       *dv,
   int      loc,
   double   real,
   double   imag
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || loc < 0 ) {
   fprintf(stderr, "\n fatal error in ZV_setEntry(%p,%d,%f,%f)"
           "\n bad input\n", dv, loc, real, imag) ;
   exit(-1) ;
}
if ( loc >= dv->maxsize ) {
   int newmaxsize = (int) 1.25*dv->maxsize ;
   if ( newmaxsize < 10 ) {
      newmaxsize = 10 ;
   }
   if ( loc >= newmaxsize ) {
      newmaxsize = loc + 1 ;
   }
   ZV_setMaxsize(dv, newmaxsize) ;
}
if ( loc >= dv->size ) {
   dv->size = loc + 1 ;
}
dv->vec[2*loc]   = real ;
dv->vec[2*loc+1] = imag ;
 
return ; }

/*--------------------------------------------------------------------*/
