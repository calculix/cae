/*  IVsort.h  */

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   sort an array of integers into ascending order using insert sort

   created -- 95sep28, cca
   ----------------------------------------------------------------
*/
void
IVisortUp (
   int   n,
   int   ivec[]
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   sort an array of integers into ascending order using insert sort

   created -- 95sep29, cca
   ----------------------------------------------------------------
*/
void
IVisortDown (
   int   n,
   int   ivec[]
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   sort an array of integers into ascending order using quick sort

   created -- 95sep28, cca
   ---------------------------------------------------------------
*/
void
IVqsortUp (
   int   n,
   int   ivec[]
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   sort an array of integers into descending order using quick sort

   created -- 95sep29, cca
   ---------------------------------------------------------------
*/
void
IVqsortDown (
   int   n,
   int   ivec[]
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   return 1 if elements in array are in ascending order
   return 0 otherwise

   created -- 95sep28, cca
   ----------------------------------------------------
*/
int
IVisascending (
   int   n,
   int   ivec[]
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   return 1 if elements in array are in descending order
   return 0 otherwise

   created -- 95sep29, cca
   -----------------------------------------------------
*/
int
IVisdescending (
   int   n,
   int   ivec[]
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   use insert sort to sort a pair of arrays.
   on output the first is in ascending order.
   ivec1[] is the array of keys
   ivec2[] is the companion array

   example,
      ivec1[] = 5 8 3 2 9
      ivec2[] = 1 2 3 4 5
   becomes
      ivec1[] = 2 3 5 8 9
      ivec2[] = 4 3 1 2 5

   created -- 95sep28, cca
   ------------------------------------------
*/
void
IV2isortUp (
   int   n,
   int   ivec1[],
   int   ivec2[]
) ;
/*
   -------------------------------------------
   use insert sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec1[] is the array of keys
   ivec2[] is the companion array

   example,
      ivec1[] = 5 8 3 2 9
      ivec2[] = 1 2 3 4 5
   becomes
      ivec1[] = 9 8 5 3 2
      ivec2[] = 5 2 1 3 4

   created -- 95sep29, cca
   -------------------------------------------
*/
void
IV2isortDown (
   int   n,
   int   ivec1[],
   int   ivec2[]
) ;
/*
   ------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in ascending order.
   ivec1[] is the array of keys
   ivec2[] is the companion array

   example,
      ivec1[] = 5 8 3 2 9
      ivec2[] = 1 2 3 4 5
   becomes
      ivec1[] = 2 3 5 8 9
      ivec2[] = 4 3 1 2 5

   created -- 95sep29, cca
   ------------------------------------------
*/
void
IV2qsortUp (
   int   n,
   int   ivec1[],
   int   ivec2[]
) ;
/*
   -------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec1[] is the array of keys
   ivec2[] is the companion array

   example,
      ivec1[] = 5 8 3 2 9
      ivec2[] = 1 2 3 4 5
   becomes
      ivec1[] = 9 8 5 3 2
      ivec2[] = 5 2 1 3 4

   created -- 95sep29, cca
   -------------------------------------------
*/
void
IV2qsortDown (
   int   n,
   int   ivec1[],
   int   ivec2[]
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   use insert sort to sort a pair of arrays.
   on output the first is in ascending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec[] =   5   8   3   2   9
      dvec[] = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec[] =   2   3   5   8   9
      dvec[] = 4.0 3.0 1.0 2.0 5.0

   created -- 95sep28, cca
   ------------------------------------------
*/
void
IVDVisortUp (
   int      n,
   int      ivec[],
   double   dvec[]
) ;
/*
   -------------------------------------------
   use insert sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec[] =   5   8   3   2   9
      dvec[] = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec[] =   9   8   5   3   2
      dvec[] = 5.0 2.0 1.0 3.0 4.0

   created -- 95sep29, cca
   -------------------------------------------
*/
void
IVDVisortDown (
   int      n,
   int      ivec[],
   double   dvec[]
) ;
/*
   ------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in ascending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec[] =   5   8   3   2   9
      dvec[] = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec[] =   2   3   5   8   9
      dvec[] = 4.0 3.0 1.0 2.0 5.0

   created -- 95sep29, cca
   ------------------------------------------
*/
void
IVDVqsortUp (
   int      n,
   int      ivec[],
   double   dvec[]
) ;
/*
   -------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec[] =   5   8   3   2   9
      dvec[] = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec[] =   9   8   5   3   2
      dvec[] = 5.0 2.0 1.0 3.0 4.0

   created -- 95sep29, cca
   -------------------------------------------
*/
void
IVDVqsortDown (
   int      n,
   int      ivec[],
   double   dvec[]
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   use insert sort to sort a trio of arrays.
   on output the first is in ascending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0

   created -- 95sep28, cca
   ------------------------------------------
*/
void
IV2DVisortUp (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   dvec[]
) ;
/*
   -------------------------------------------
   use insert sort to sort a trio or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0

   created -- 95sep29, cca
   -------------------------------------------
*/
void
IV2DVisortDown (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   dvec[]
) ;
/*
   ------------------------------------------
   use quick sort to sort a trio or arrays.
   on output the first is in ascending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0

   created -- 95sep29, cca
   ------------------------------------------
*/
void
IV2DVqsortUp (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   dvec[]
) ;
/*
   -------------------------------------------
   use quick sort to sort a trio or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   dvec[] is the companion array

   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0

   created -- 95sep29, cca
   -------------------------------------------
*/
void
IV2DVqsortDown (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   dvec[]
) ;
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
) ;
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
) ;
/*

-----------------------------------------------------------------
   sort the entries in ivec1[] and ivec2[] into lexicographic
order,
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
