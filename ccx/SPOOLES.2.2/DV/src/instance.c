/*  instance.c  */

#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   return 1 if the entries are owned by the object
   return 0 otherwise

   created -- 96jun22, cca
   -----------------------------------------------
*/
int
DV_owned (
   DV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_owned(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->owned) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   return the vector size

   created -- 95oct06, cca
   -----------------------
*/
int
DV_maxsize (
   DV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_maxsize(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->maxsize) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   return the vector size

   created -- 95oct06, cca
   -----------------------
*/
int
DV_size (
   DV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_size(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->size) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   return the loc'th entry of a vector.
   note: if loc is out of range then 0.0 is returned
 
   created -- 96jun29, cca
   -------------------------------------------------
*/
double 
DV_entry (
   DV    *dv,
   int   loc
) {
if ( dv == NULL || dv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in DV_entry(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
if ( loc < 0 || loc >= dv->size ) {
   return(0.0) ;
} else {
   return(dv->vec[loc]) ; 
}
}
 
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return a pointer to the object's entries array

   created -- 95oct06, cca
   ----------------------------------------------
*/
double *
DV_entries (
   DV   *dv
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_entries(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
return(dv->vec) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   fill *psize with the vector's size
   and *pentries with the address of the vector

   created -- 95oct06, cca
   --------------------------------------------
*/
void
DV_sizeAndEntries (
   DV       *dv,
   int      *psize,
   double   **pentries
) {
if ( dv == NULL || psize == NULL || pentries == NULL ) {
   fprintf(stderr, "\n fatal error in DV_sizeAndEntries(%p,%p,%p)"
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
 
   created -- 96jul14, cca
   ---------------------------
*/
void
DV_setEntry (
   DV       *dv,
   int      loc,
   double   value
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || loc < 0 ) {
   fprintf(stderr, "\n fatal error in DV_setEntry(%p,%d,%f)"
           "\n bad input\n", dv, loc, value) ;
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
   DV_setMaxsize(dv, newmaxsize) ;
}
if ( loc >= dv->size ) {
   dv->size = loc + 1 ;
}
dv->vec[loc] = value ;
 
return ; }

/*--------------------------------------------------------------------*/
