/*  basics.C  */

#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   shift the base of the entries and adjust the dimensions
   note: this is a dangerous operation because the dv->vec
   does not point to the base of the entries any longer,
   and thus if the object owns its entries and it is called
   to resize them or to free them, malloc and free will choke.
 
   USE WITH CAUTION!
 
   created  -- 96aug25, cca
   modified -- 96aug28, cca
      structure changed
   -----------------------------------------------------------
*/
void
DV_shiftBase (
   DV    *dv,
   int    offset
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_shiftBase(%p,%d)"
           "\n bad input\n", dv, offset) ;
   exit(-1) ;
}
dv->vec     += offset ;
dv->maxsize -= offset ;
dv->size    -= offset ;
 
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   push an entry onto the list

   created -- 95oct06, cca
   --------------------------------------
*/
void
DV_push (
   DV       *dv,
   double   val
) {
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_push(%p,%f)"
           "\n bad input\n", dv, val) ;
   exit(-1) ;
}
if ( dv->size == dv->maxsize ) {
   if ( dv->maxsize > 0 ) {
      DV_setMaxsize(dv, 2*dv->maxsize) ;
   } else {
      DV_setMaxsize(dv, 10) ;
   }
}
dv->vec[dv->size++] = val ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   minimum and maximum entries and sum

   created -- 95oct06, cca
   -----------------------------------
*/
double 
DV_min ( 
   DV   *dv
) {
int   i ;

if ( dv == NULL || dv->size <= 0 || dv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in DV_min(%p), size = %d, vec = %p",
           dv, dv->size, dv->vec) ;
   exit(-1) ;
}
return(DVmin(dv->size, dv->vec, &i)) ; }

double 
DV_max ( 
   DV   *dv
) {
int   i ;

if ( dv == NULL || dv->size <= 0 || dv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in DV_max(%p), size = %d, vec = %p",
           dv, dv->size, dv->vec) ;
   exit(-1) ;
}
return(DVmax(dv->size, dv->vec, &i)) ; }

double 
DV_sum ( 
   DV   *dv
) {
if ( dv == NULL || dv->size <= 0 || dv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in DV_sum(%p), size = %d, vec = %p",
           dv, dv->size, dv->vec) ;
   exit(-1) ;
}
return(DVsum(dv->size, dv->vec)) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   sort each index list into ascending or descending order

   created -- 95oct06, cca
   -------------------------------------------------------
*/
void
DV_sortUp ( 
   DV   *dv
) {
if ( dv == NULL || dv->size <= 0 || dv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DV_sortUp(%p), size = %d, vec = %p",
           dv, dv->size, dv->vec) ;
   exit(-1) ;
}
DVqsortUp(dv->size, dv->vec) ;

return ; }

void
DV_sortDown ( 
   DV   *dv
) {
if ( dv == NULL || dv->size <= 0 || dv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DV_sortDown(%p), size = %d, vec = %p",
           dv, dv->size, dv->vec) ;
   exit(-1) ;
}
DVqsortDown(dv->size, dv->vec) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   ramp the entries

   created -- 95oct06, cca
   -----------------------
*/
void
DV_ramp (
   DV       *dv,
   double   base,
   double   incr
) {
if ( dv == NULL || dv->size <= 0 || dv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DV_ramp(%p,%f,%f), size = %d, vec = %p",
           dv, base, incr, dv->size, dv->vec) ;
   exit(-1) ;
}
DVramp(dv->size, dv->vec, base, incr) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   shuffle the list

   created -- 95oct06, cca
   -----------------------
*/
void
DV_shuffle (
   DV    *dv, 
   int   seed
) {
if ( dv == NULL || dv->size <= 0 || dv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DV_shuffle(%p,%d), size = %d, vec = %p",
           dv, seed, dv->size, dv->vec) ;
   exit(-1) ;
}
DVshuffle(dv->size, dv->vec, seed) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
DV_sizeOf ( 
   DV   *dv
) {
int   nbytes ;

if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_sizeOf(%p)"
           "\n bad input \n", dv) ;
   exit(-1) ;
}
nbytes = sizeof(struct _DV) ;
if ( dv->owned == 1 ) {
   nbytes += + dv->maxsize*sizeof(double) ; 
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   iterator :
   return the pointer to the head of the vector

   created -- 95oct06, cca
   --------------------------------------------
*/
double *
DV_first (
   DV   *dv
) {
double   *pd ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_first(%p)"
           "\n bad input", dv) ;
   exit(-1) ;
}
if ( dv->size == 0 ) {
   pd = NULL ;
} else {
   pd = dv->vec ;
}
return(pd) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   iterator :
   return the pointer to the next location in the vector

   created -- 95oct06, cca
   -----------------------------------------------------
*/
double *
DV_next (
   DV       *dv,
   double   *pd
) {
int   offset ;
/*
   ---------------
   check the input
   ---------------
*/
if ( pd == NULL ) {
   fprintf(stderr, "\n fatal error in DV_next(%p,%p)"
           "\n bad input", dv, pd) ;
   fflush(stderr) ;
   exit(-1) ;
}
/*
   ---------------
   check the input
   ---------------
*/
if ( (offset = pd - dv->vec) < 0 || offset >= dv->size ) {
/*
   -----------------------------
   error, offset is out of range
   -----------------------------
*/
   fprintf(stderr, "\n fatal error in DV_next(%p,%p)"
           "\n offset = %d, must be in [0,%d)",
           dv, pd, offset, dv->size) ;
   fflush(stderr) ;
   exit(-1) ;
} else if ( offset == dv->size - 1 ) {
/*
   ----------------------------
   end of the list, return NULL
   ----------------------------
*/
   pd = NULL ;
} else {
/*
   ----------------------------------------
   middle of the list, return next location
   ----------------------------------------
*/
   pd++ ;
}
return(pd) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   fill a vector with a value

   created -- 96jun22, cca
   --------------------------
*/
void
DV_fill (
   DV       *dv,
   double   value
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_fill(%p,%f)"
           "\n bad input\n", dv, value) ;
   exit(-1) ;
}
if ( dv->size > 0 ) {
   DVfill(dv->size, dv->vec, value) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   fill a vector with zeros

   created -- 98jun02, cca
   --------------------------
*/
void
DV_zero (
   DV   *dv
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL ) {
   fprintf(stderr, "\n fatal error in DV_zero(%p)"
           "\n bad input\n", dv) ;
   exit(-1) ;
}
if ( dv->size > 0 ) {
   DVzero(dv->size, dv->vec) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   copy entries from dv2 into dv1.
   note: this is a "mapped" copy, 
   dv1 and dv2 need not be the same size.
 
   created -- 96aug31, cca
   --------------------------------------
*/
void
DV_copy (
   DV   *dv1,
   DV   *dv2
) {
int      ii, size ;
double   *vec1, *vec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv1 == NULL || dv2 == NULL ) {
   fprintf(stderr, "\n fatal error in DV_copy(%p,%p)"
           "\n bad input\n", dv1, dv2) ;
   exit(-1) ;
}
size = dv1->size ;
if ( size > dv2->size ) {
   size = dv2->size ;
}
vec1 = dv1->vec ;
vec2 = dv2->vec ;
for ( ii = 0 ; ii < size ; ii++ ) {
   vec1[ii] = vec2[ii] ;
}
return ; }
 
/*--------------------------------------------------------------------*/
