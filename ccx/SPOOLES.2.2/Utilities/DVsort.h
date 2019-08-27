/*  DVsort.h  */

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   The insert and quick sort methods in this file are based on

   Jon L. Bentley and M. Douglas McIlroy,
   "Engineering a sort function",
   Software -- Practice and Experience, vol 23(11), 1249-1265,
   November 1993.

   This quick sort method uses 
      1) a median of three medians to find a split value, and
      2) split-end partitioning to handle many elements equal
         to the split value.
   -----------------------------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   sort an array of doubles into ascending order using insert sort

   created -- 95nov25, cca
   ----------------------------------------------------------------
*/
void
DVisortUp (
   int      n,
   double   dvec[]
) ;
/*
   -----------------------------------------------------------------
   sort an array of doubles into descending order using insert sort

   created -- 95sep29, cca
   -----------------------------------------------------------------
*/
void
DVisortDown (
   int      n,
   double   dvec[]
) ;
/*
   ------------------------------------------
   use insert sort to sort a pair or arrays.
   on output the first is in ascending order.
   dvec1[] is the array of keys
   dvec2[] is the companion array

   example,
      dvec1[] = 5 8 3 2 9
      dvec2[] = 1 2 3 4 5
   becomes
      dvec1[] = 2 3 5 8 9
      dvec2[] = 4 3 1 2 5

   created -- 95nov25, cca
   ------------------------------------------
*/
void
DV2isortUp (
   int      n,
   double   dvec1[],
   double   dvec2[]
) ;
/*
   -------------------------------------------
   use insert sort to sort a pair or arrays.
   on output the first is in descending order.
   dvec1[] is the array of keys
   dvec2[] is the companion array

   example,
      dvec1[] = 5 8 3 2 9
      dvec2[] = 1 2 3 4 5
   becomes
      dvec1[] = 9 8 5 3 2
      dvec2[] = 5 2 1 3 4

   created -- 95sep29, cca
   -------------------------------------------
*/
void
DV2isortDown (
   int      n,
   double   dvec1[],
   double   dvec2[]
) ;
/*
   ---------------------------------------------------------------
   sort an array of doubles into ascending order using quick sort

   created -- 95nov25, cca
   ---------------------------------------------------------------
*/
void
DVqsortUp (
   int      n,
   double   dvec[]
) ;
/*
   ----------------------------------------------------------------
   sort an array of doubles into descending order using quick sort

   created -- 95sep29, cca
   ----------------------------------------------------------------
*/
void
DVqsortDown (
   int      n,
   double   dvec[]
) ;
/*
   ------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in ascending order.
   dvec1[] is the array of keys
   dvec2[] is the companion array

   example,
      dvec1[] = 5 8 3 2 9
      dvec2[] = 1 2 3 4 5
   becomes
      dvec1[] = 2 3 5 8 9
      dvec2[] = 4 3 1 2 5

   created -- 95sep29, cca
   ------------------------------------------
*/
void
DV2qsortUp (
   int      n,
   double   dvec1[],
   double   dvec2[]
) ;
/*
   -------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in descending order.
   dvec1[] is the array of keys
   dvec2[] is the companion array

   example,
      dvec1[] = 5 8 3 2 9
      dvec2[] = 1 2 3 4 5
   becomes
      dvec1[] = 9 8 5 3 2
      dvec2[] = 5 2 1 3 4

   created -- 95sep29, cca
   -------------------------------------------
*/
void
DV2qsortDown (
   int      n,
   double   dvec1[],
   double   dvec2[]
) ;
/*
   ----------------------------------------------------
   return 1 if elements in array are in ascending order
   return 0 otherwise

   created -- 95nov25, cca
   ----------------------------------------------------
*/
int
DVisascending (
   int      n,
   double   dvec[]
) ;
/*
   -----------------------------------------------------
   return 1 if elements in array are in descending order
   return 0 otherwise

   created -- 95sep29, cca
   -----------------------------------------------------
*/
int
DVisdescending (
   int      n,
   double   dvec[]
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   use insert sort to sort a pair of arrays.
   on output the first is in ascending order.
   dvec[] is the array of keys
   ivec[] is the companion array
 
   example,
      ivec[] =   4   2   7   5   6
      dvec[] = 5.0 8.0 3.0 2.0 9.0
   becomes
      ivec[] =   5   7   4   2   6
      dvec[] = 2.0 3.0 5.0 8.0 9.0
 
   created -- 97apr03, cca
   ------------------------------------------
*/
void
DVIVisortUp (
   int      n,
   double   dvec[],
   int      ivec[]
) ;
/*
   -------------------------------------------
   use insert sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   dvec[] is the companion array
 
   example,
      ivec[] =   4   2   7   5   6
      dvec[] = 5.0 8.0 3.0 2.0 9.0
   becomes
      ivec[] =   5   7   4   2   6
      dvec[] = 2.0 3.0 5.0 8.0 9.0
 
   created -- 95sep29, cca
   -------------------------------------------
*/
void
DVIVisortDown (
   int      n,
   double   dvec[],
   int      ivec[]
) ;
/*
   ------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in ascending order.
   ivec[] is the array of keys
   dvec[] is the companion array
 
   example,
      ivec[] =   4   2   7   5   6
      dvec[] = 5.0 8.0 3.0 2.0 9.0
   becomes
      ivec[] =   5   7   4   2   6
      dvec[] = 2.0 3.0 5.0 8.0 9.0
 
   created -- 97apr03, cca
   ------------------------------------------
*/
void
DVIVqsortUp (
   int      n,
   double   dvec[],
   int      ivec[]
) ;
/*
   -------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   dvec[] is the companion array
 
   example,
      ivec[] =   4   2   7   5   6
      dvec[] = 5.0 8.0 3.0 2.0 9.0
   becomes
      ivec[] =   5   7   4   2   6
      dvec[] = 2.0 3.0 5.0 8.0 9.0
 
   created -- 97apr03, cca
   -------------------------------------------
*/
void
DVIVqsortDown (
   int      n,
   double   dvec[],
   int      ivec[]
) ;
/*--------------------------------------------------------------------*/
