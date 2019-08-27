/*  sortutil.c  */

#include "../Utilities.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   sort the entries in ivec[] into ascending order,
   then compress out the duplicates

   return value -- # of compressed entries

   created -- 97dec18, cca
   ------------------------------------------------
*/
int
IVsortUpAndCompress (
   int   n,
   int   ivec[]
) {
int   first, ierr, ii, key ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n < 0 || ivec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in IVsortAndCompress(%d,%p)"
           "\n bad input, n = %d, ivec = %p",
           n, ivec, n, ivec) ;
   exit(-1) ;
}
if ( n == 0 ) {
   return(0) ;
}
/*
   ----------------------------------
   sort the vector in ascending order
   ----------------------------------
*/
IVqsortUp(n, ivec) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ivec[] after sort up") ;
IVfp80(stdout, n, ivec, 80, &ierr) ;
#endif
/*
   --------------------
   purge the duplicates
   --------------------
*/
first = 1 ;
key   = ivec[0] ;
#if MYDEBUG > 0
fprintf(stdout, "\n first = %d, key = %d, ivec[%d] = %d",
        first, key, 0, ivec[0]) ;
#endif
for ( ii = 1 ; ii < n ; ii++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n first = %d, key = %d, ivec[%d] = %d",
           first, key, 0, ivec[0]) ;
#endif
   if ( key != ivec[ii] ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n setting ivec[%d] = %d",
              first, ivec[ii]) ;
#endif
      ivec[first++] = key = ivec[ii] ;
   }
}
return(first) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   sort the entries in ivec[] into ascending order,
   using dvec[] as a companion vector.
   then compress out the duplicates integer keys,
   summing the double values

   return value -- # of compressed entries

   created -- 97dec18, cca
   ------------------------------------------------
*/
int
IVDVsortUpAndCompress (
   int      n,
   int      ivec[],
   double   dvec[]
) {
int   first, ierr, ii, key ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n < 0 || ivec == NULL || dvec == NULL) {
   fprintf(stderr, 
           "\n fatal error in IVDVsortAndCompress(%d,%p,%p)"
           "\n bad input, n = %d, ivec = %p, dvec = %p",
           n, ivec, dvec, n, ivec, dvec) ;
   exit(-1) ;
}
if ( n == 0 ) {
   return(0) ;
}
/*
   ---------------------------------
   sort ivec[] in ascending order
   with dvec[] as a companion vector
   ---------------------------------
*/
IVDVqsortUp(n, ivec, dvec) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ivec[] after sort up") ;
IVfp80(stdout, n, ivec, 80, &ierr) ;
fprintf(stdout, "\n dvec[] after sort up") ;
DVfprintf(stdout, n, dvec) ;
#endif
first = 1 ;
key   = ivec[0] ;
#if MYDEBUG > 0
fprintf(stdout, "\n first = %d, key = %d, ivec[%d] = %d",
        first, key, 0, ivec[0]) ;
#endif
for ( ii = 1 ; ii < n ; ii++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n first = %d, key = %d, ivec[%d] = %d",
           first, key, 0, ivec[0]) ;
#endif
   if ( key == ivec[ii] ) {
      dvec[first-1] += dvec[ii] ;
   } else if ( key != ivec[ii] ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n setting ivec[%d] = %d",
              first, ivec[ii]) ;
#endif
      ivec[first] = key = ivec[ii] ;
      dvec[first] = dvec[ii] ;
      first++ ;
   }
}
return(first) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   sort the entries in ivec[] into ascending order,
   using zvec[] as a complex companion vector.
   then compress out the duplicates integer keys,
   summing the complex values

   return value -- # of compressed entries

   created -- 98jan28, cca
   ------------------------------------------------
*/
int
IVZVsortUpAndCompress (
   int      n,
   int      ivec[],
   double   zvec[]
) {
int   first, ierr, ii, key ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n < 0 || ivec == NULL || zvec == NULL) {
   fprintf(stderr, 
           "\n fatal error in IVZVsortAndCompress(%d,%p,%p)"
           "\n bad input, n = %d, ivec = %p, zvec = %p",
           n, ivec, zvec, n, ivec, zvec) ;
   exit(-1) ;
}
if ( n == 0 ) {
   return(0) ;
}
/*
   ---------------------------------
   sort ivec[] in ascending order
   with zvec[] as a companion vector
   ---------------------------------
*/
IVZVqsortUp(n, ivec, zvec) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ivec[] after sort up") ;
IVfp80(stdout, n, ivec, 80, &ierr) ;
fprintf(stdout, "\n zvec[] after sort up") ;
DVfprintf(stdout, 2*n, zvec) ;
#endif
first = 1 ;
key   = ivec[0] ;
#if MYDEBUG > 0
fprintf(stdout, "\n first = %d, key = %d, ivec[%d] = %d",
        first, key, 0, ivec[0]) ;
#endif
for ( ii = 1 ; ii < n ; ii++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n first = %d, key = %d, ivec[%d] = %d",
           first, key, 0, ivec[0]) ;
#endif
   if ( key == ivec[ii] ) {
      zvec[2*(first-1)] += zvec[2*ii] ;
      zvec[2*(first-1)+1] += zvec[2*ii+1] ;
   } else if ( key != ivec[ii] ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n setting ivec[%d] = %d",
              first, ivec[ii]) ;
#endif
      ivec[first] = key = ivec[ii] ;
      zvec[2*first]   = zvec[2*ii] ;
      zvec[2*first+1] = zvec[2*ii+1] ;
      first++ ;
   }
}
return(first) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   sort the entries in ivec1[] and ivec2[] into lexicographic order,
   then compress out the duplicates

   return value -- # of compressed entries

   created -- 97dec18, cca
   -----------------------------------------------------------------
*/
int
IV2sortUpAndCompress (
   int   n,
   int   ivec1[],
   int   ivec2[]
) {
int   first, ierr, ii, key, length, start ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n < 0 || ivec1 == NULL || ivec2 == NULL) {
   fprintf(stderr, 
           "\n fatal error in IV2sortAndCompress(%d,%p,%p)"
           "\n bad input, n = %d, ivec1 = %p, ivec2 = %p",
           n, ivec1, ivec2, n, ivec1, ivec2) ;
   exit(-1) ;
}
if ( n == 0 ) {
   return(0) ;
}
/*
   ------------------------------------------------------------------
   sort ivec1[] in ascending order with ivec2[] as a companion vector
   ------------------------------------------------------------------
*/
IV2qsortUp(n, ivec1, ivec2) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ivec1[] after sort up") ;
IVfp80(stdout, n, ivec1, 80, &ierr) ;
fprintf(stdout, "\n ivec2[] after sort up") ;
IVfp80(stdout, n, ivec2, 80, &ierr) ;
#endif
first = start = 0 ;
key   = ivec1[0] ;
for ( ii = 1 ; ii < n ; ii++ ) {
   if ( key != ivec1[ii] ) {
      length = IVsortUpAndCompress(ii - start, ivec2 + start) ;
      IVfill(length, ivec1 + first, key) ;
      IVcopy(length, ivec2 + first, ivec2 + start) ;
      start = ii ;
      first += length ;
      key = ivec1[ii] ;
   }
}
length = IVsortUpAndCompress(n - start, ivec2 + start) ;
IVfill(length, ivec1 + first, key) ;
IVcopy(length, ivec2 + first, ivec2 + start) ;
first += length ;

return(first) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   sort the entries in ivec1[] and ivec2[] into lexicographic order,
   using dvec[] as a companion vector. then compress out the duplicates
   and sum the same values of dvec[]

   return value -- # of compressed entries

   created -- 97dec18, cca
   --------------------------------------------------------------------
*/
int
IV2DVsortUpAndCompress (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   dvec[]
) {
int   first, ierr, ii, key, length, start ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n < 0 || ivec1 == NULL || ivec2 == NULL || dvec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in IV2DVsortAndCompress(%d,%p,%p,%p)"
           "\n bad input, n = %d, ivec1 = %p, ivec2 = %p, dvec = %p",
           n, ivec1, ivec2, dvec, n, ivec1, ivec2, dvec) ;
   exit(-1) ;
}
if ( n == 0 ) {
   return(0) ;
}
/*
   ---------------------------------------
   sort ivec1[] in ascending order with 
   ivec2[] and dvec[] as companion vectors
   ---------------------------------------
*/
IV2DVqsortUp(n, ivec1, ivec2, dvec) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ivec1[] after sort up") ;
IVfp80(stdout, n, ivec1, 80, &ierr) ;
fprintf(stdout, "\n ivec2[] after sort up") ;
IVfp80(stdout, n, ivec2, 80, &ierr) ;
fprintf(stdout, "\n dvec[] after sort up") ;
DVfprintf(stdout, n, dvec) ;
#endif
first = start = 0 ;
key   = ivec1[0] ;
for ( ii = 1 ; ii < n ; ii++ ) {
   if ( key != ivec1[ii] ) {
      length = IVDVsortUpAndCompress(ii - start, 
                                     ivec2 + start, dvec + start) ;
      IVfill(length, ivec1 + first, key) ;
      IVcopy(length, ivec2 + first, ivec2 + start) ;
      DVcopy(length, dvec  + first, dvec  + start) ;
      start = ii ;
      first += length ;
      key = ivec1[ii] ;
   }
}
length = IVDVsortUpAndCompress(n - start, ivec2 + start, dvec + start) ;
IVfill(length, ivec1 + first, key) ;
IVcopy(length, ivec2 + first, ivec2 + start) ;
DVcopy(length, dvec  + first, dvec  + start) ;
first += length ;

return(first) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   sort the entries in ivec1[] and ivec2[] into lexicographic order,
   using zvec[] as a companion vector. then compress out the duplicates
   and sum the same values of zvec[]

   return value -- # of compressed entries

   created -- 98jan28, cca
   --------------------------------------------------------------------
*/
int
IV2ZVsortUpAndCompress (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   zvec[]
) {
int   first, ierr, ii, key, length, start ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n < 0 || ivec1 == NULL || ivec2 == NULL || zvec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in IV2ZVsortAndCompress(%d,%p,%p,%p)"
           "\n bad input, n = %d, ivec1 = %p, ivec2 = %p, zvec = %p",
           n, ivec1, ivec2, zvec, n, ivec1, ivec2, zvec) ;
   exit(-1) ;
}
if ( n == 0 ) {
   return(0) ;
}
/*
   ---------------------------------------
   sort ivec1[] in ascending order with 
   ivec2[] and zvec[] as companion vectors
   ---------------------------------------
*/
IV2ZVqsortUp(n, ivec1, ivec2, zvec) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n after IV2ZVqsortUp") ;
for ( ii = 0 ; ii < n ; ii++ ) {
   fprintf(stdout, "\n < %12d, %12d, %12.4e, %12.4e >",
           ivec1[ii], ivec2[ii], zvec[2*ii], zvec[2*ii+1]) ;
}
#endif
first = start = 0 ;
key   = ivec1[0] ;
for ( ii = 1 ; ii < n ; ii++ ) {
   if ( key != ivec1[ii] ) {
      length = IVZVsortUpAndCompress(ii - start, 
                                     ivec2 + start, zvec + 2*start) ;
      IVfill(length,   ivec1 + first,   key) ;
      IVcopy(length,   ivec2 + first,   ivec2 + start) ;
      DVcopy(2*length, zvec  + 2*first, zvec  + 2*start) ;
      start = ii ;
      first += length ;
      key = ivec1[ii] ;
   }
}
length = IVZVsortUpAndCompress(n - start, ivec2 + start, 
                               zvec + 2*start) ;
IVfill(length,   ivec1 + first,   key) ;
IVcopy(length,   ivec2 + first,   ivec2 + start) ;
DVcopy(2*length, zvec  + 2*first, zvec  + 2*start) ;
first += length ;

return(first) ; }

/*--------------------------------------------------------------------*/
