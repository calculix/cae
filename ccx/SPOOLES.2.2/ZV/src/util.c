/*  basics.C  */

#include "../ZV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   shift the base of the entries and adjust the dimensions
   note: this is a dangerous operation because the zv->vec
   does not point to the base of the entries any longer,
   and thus if the object owns its entries and it is called
   to resize them or to free them, malloc and free will choke.
 
   USE WITH CAUTION!
 
   created  -- 98jan22, cca
   -----------------------------------------------------------
*/
void
ZV_shiftBase (
   ZV    *zv,
   int    offset
) {
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_shiftBase(%p,%d)"
           "\n bad input\n", zv, offset) ;
   exit(-1) ;
}
zv->vec     += 2*offset ;
zv->maxsize -= offset ;
zv->size    -= offset ;
 
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   push an entry onto the list

   created -- 95oct06, cca
   --------------------------------------
*/
void
ZV_push (
   ZV       *zv,
   double   real,
   double   imag
) {
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_push(%p,%f,%f)"
           "\n bad input\n", zv, real, imag) ;
   exit(-1) ;
}
if ( zv->size == zv->maxsize ) {
   if ( zv->maxsize > 0 ) {
      ZV_setMaxsize(zv, 2*zv->maxsize) ;
   } else {
      ZV_setMaxsize(zv, 10) ;
   }
}
zv->vec[2*zv->size]   = real ;
zv->vec[2*zv->size+1] = real ;
zv->size++ ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   minimum and maximum entries

   created -- 95oct06, cca
   ---------------------------
*/
double 
ZV_minabs ( 
   ZV   *zv
) {
if ( zv == NULL || zv->size <= 0 || zv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ZV_minabs(%p), size = %d, vec = %p",
           zv, zv->size, zv->vec) ;
   exit(-1) ;
}
return(ZVminabs(zv->size, zv->vec)) ; }

double 
ZV_maxabs ( 
   ZV   *zv
) {
if ( zv == NULL || zv->size <= 0 || zv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ZV_maxabs(%p), size = %d, vec = %p",
           zv, zv->size, zv->vec) ;
   exit(-1) ;
}
return(ZVmaxabs(zv->size, zv->vec)) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
ZV_sizeOf ( 
   ZV   *zv
) {
int   nbytes ;

if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_sizeOf(%p)"
           "\n bad input \n", zv) ;
   exit(-1) ;
}
nbytes = sizeof(struct _ZV) ;
if ( zv->owned == 1 ) {
   nbytes += 2*zv->maxsize*sizeof(double) ; 
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   fill a vector with a value

   created -- 96jun22, cca
   --------------------------
*/
void
ZV_fill (
   ZV       *zv,
   double   real,
   double   imag
) {
double   *vec ;
int      ii, jj, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_fill(%p,%f,%f)"
           "\n bad input\n", zv, real, imag) ;
   exit(-1) ;
}
if ( (size = zv->size) > 0 ) {
   if ( (vec = zv->vec) == NULL ) {
      fprintf(stderr, "\n fatal error in ZV_fill(%p,%f,%f)"
              "\n vec = NULL\n", zv, real, imag) ;
      exit(-1) ;
   }
   for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
      vec[jj]   = real ;
      vec[jj+1] = imag ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------
   fill a vector with zeros

   created -- 98jun02, cca
   ------------------------
*/
void
ZV_zero (
   ZV   *zv
) {
double   *vec ;
int      ii, jj, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_zero(%p)"
           "\n bad input\n", zv) ;
   exit(-1) ;
}
if ( (size = zv->size) > 0 ) {
   if ( (vec = zv->vec) == NULL ) {
      fprintf(stderr, "\n fatal error in ZV_zero(%p)"
              "\n vec = NULL\n", zv) ;
      exit(-1) ;
   }
   for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
      vec[jj]   = 0.0 ;
      vec[jj+1] = 0.0 ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   copy entries from zv2 into zv1.
   note: this is a "mapped" copy, 
   zv1 and zv2 need not be the same size.
 
   created -- 96aug31, cca
   --------------------------------------
*/
void
ZV_copy (
   ZV   *zv1,
   ZV   *zv2
) {
int      ii, jj, size ;
double   *vec1, *vec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv1 == NULL || zv2 == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_copy(%p,%p)"
           "\n bad input\n", zv1, zv2) ;
   exit(-1) ;
}
size = zv1->size ;
if ( size > zv2->size ) {
   size = zv2->size ;
}
vec1 = zv1->vec ;
vec2 = zv2->vec ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   vec1[jj]   = vec2[jj] ;
   vec1[jj+1] = vec2[jj+1] ;
}
return ; }
 
/*--------------------------------------------------------------------*/
