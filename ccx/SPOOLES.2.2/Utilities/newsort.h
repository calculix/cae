/*  newsort.h.c  */

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   return 1 if elements in array are in ascending order
   return 0 otherwise

   created -- 98jan28, cca
   ----------------------------------------------------
*/
int
IVisascending (
   int   n,
   int   ivec[]
) ;
/*
   -----------------------------------------------------
   return 1 if elements in array are in descending order
   return 0 otherwise

   created -- 98jan28, cca
   -----------------------------------------------------
*/
int
IVisdescending (
   int   n,
   int   ivec[]
) ;
/*
   ----------------------------------------------------
   return 1 if elements in array are in ascending order
   return 0 otherwise

   created -- 98jan28, cca
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

   created -- 98jan28, cca
   -----------------------------------------------------
*/
int
DVisdescending (
   int      n,
   double   dvec[]
) ;
/*====   INSERT SORT METHODS    ======================================*/
/*
   ----------------------------------------------------------------
   sort an array of integers into ascending order using insert sort

   created -- 98jan28, cca
   ----------------------------------------------------------------
*/
void
IVisortUp (
   int   n,
   int   ivec[]
) ;
/*
   ------------------------------------------
   use insert sort to sort a pair or arrays.
   on output the first is in ascending order.
   ivec1[] is the array of keys
   ivec2[] is the companion array
 
   example,
      ivec1[] = 5 8 3 2 9
      ivec2[] = 1 2 3 4 5
   becomes
      ivec1[] = 2 3 5 8 9
      ivec2[] = 4 3 1 2 5
 
   created -- 98jan28, cca
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
 
   created -- 98jan28, cca
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
 
   created -- 98jan28, cca
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
 
   created -- 98jan28, cca
   -------------------------------------------
*/
void
IVDVisortDown (
   int      n,
   int      ivec[],
   double   dvec[]
) ;
/*
   -----------------------------------------------------------------
   sort an array of integers into descending order using insert sort

   created -- 98jan28, cca
   -----------------------------------------------------------------
*/
void
IVisortDown (
   int   n,
   int   ivec[]
) ;
/*
   ------------------------------------------
   use insert sort to sort a trio of arrays.
   on output the first is in ascending order.
   ivec1[] is the array of keys
   ivec2[] is a companion array
   dvec[]  is a companion array
 
   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0
 
   created -- 98jan28, cca
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
   ivec1[] is the array of keys
   ivec2[] is a companion array
   dvec[]  is a companion array
 
   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0
 
   created -- 98jan28, cca
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
   -------------------------------------------------------------
   use insert sort to sort a pair of arrays.
   on output the first is in ascending order.
   ivec[] is the array of keys
   zvec[] is the companion array of complex
 
   example,
      ivec[] =     5         8         3         2         9
      zvec[] = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec[] =     2         3         5         8         9
      zvec[] = (4.0,4.5) (3.0,3.5) (1.0,1.5) (2.0,2.5) (5.0,5.5)
 
   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
IVZVisortUp (
   int      n,
   int      ivec[],
   double   zvec[]
) ;
/*
   -------------------------------------------------------------
   use insert sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   zvec[] is the companion array
 
   example,
      ivec[] =     5         8         3         2         9
      zvec[] = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec[] =     9         8         5         3         2
      zvec[] = (5.0,5.5) (2.0,2.5) (1.0,1.5) (3.0,3.5) (4.0,4.5)
 
   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
IVZVisortDown (
   int      n,
   int      ivec[],
   double   zvec[]
) ;
/*
   --------------------------------------------------------------
   use insert sort to sort a trio of arrays.
   on output the first is in ascending order.
   ivec1[] is the array of keys
   ivec2[] is a companion array
   zvec[]  is a companion array
 
   example,
      ivec1[] =     5         8         3         2         9
      ivec2[] =     6         7         8         9         5
      zvec[]  = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec1[] =     2         3         5         8         9
      ivec2[] =     9         8         6         7         5
      zvec[]  = (4.0,4.5) (3.0,3.5) (1.0,1.5) (2.0,2.5) (5.0,5.5)
 
 
   created -- 98jan28, cca
   --------------------------------------------------------------
*/
void
IV2ZVisortUp (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   zvec[]
) ;
/*
  --------------------------------------------------------------
   use insert sort to sort a trio or arrays.
   on output the first is in descending order.
   ivec1[] is the array of keys
   ivec2[] is a companion array
   zvec[]  is a companion array
 
   example,
      ivec1[] =     5         8         3         2         9
      ivec2[] =     6         7         8         9         5
      zvec[]  = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec1[] =     9         8         5         3         2
      ivec2[] =     5         7         6         8         9
      zvec[]  = (5.0,5.5) (2.0,2.5) (1.0,1.5) (4.0,4.5) (4.0,4.5)
 
   created -- 98jan28, cca
   --------------------------------------------------------------
*/
void
IV2ZVisortDown (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   zvec[]
) ;
/*
   ----------------------------------------------------------------
  sort an array of doubles into ascending order using insert sort
 
   created -- 98jan28, cca
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

   created -- 98jan28, cca
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
 
   created -- 98jan28, cca
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
 
   created -- 98jan28, cca
   -------------------------------------------
*/
void
DV2isortDown (
   int      n,
   double   dvec1[],
   double   dvec2[]
) ;
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
 
   created -- 98jan28, cca
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
 
   created -- 98jan28, cca
   -------------------------------------------
*/
void
DVIVisortDown (
   int      n,
   double   dvec[],
   int      ivec[]
) ;
/*====   QUICKSORT METHODS    ========================================*/
/*
   ---------------------------------------------------------------
   sort an array of integers into ascending order using quick sort

   created -- 98jan28, cca
   ---------------------------------------------------------------
*/
void
IVqsortUp (
   int   n,
   int   ivec[]
) ;
/*
   ----------------------------------------------------------------
   sort an array of integers into descending order using quick sort

   created -- 98jan28, cca
   ----------------------------------------------------------------
*/
void
IVqsortDown (
   int   n,
   int   ivec[]
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
 
   created -- 98jan28, cca
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
 
   created -- 98jan28, cca
   -------------------------------------------
*/
void
IV2qsortDown (
   int   n,
   int   ivec1[],
   int   ivec2[]
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

   created -- 98jan28, cca
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

   created -- 98jan28, cca
   -------------------------------------------
*/
void
IVDVqsortDown (
   int      n,
   int      ivec[],
   double   dvec[]
) ;
/*
   ------------------------------------------
   use quick sort to sort a trio or arrays.
   on output the first is in ascending order.
   ivec1[] is the array of keys
   ivec2[] is a companion array
   dvec[]  is a companion array

   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0

   created -- 98jan28, cca
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
   ivec1[] is the array of keys
   ivec2[] is a companion array
   dvec[]  is a companion array

   example,
      ivec1[] =   5   8   3   2   9
      ivec2[] =   6   7   8   9   5
      dvec[]  = 1.0 2.0 3.0 4.0 5.0
   becomes
      ivec1[] =   2   3   5   8   9
      ivec2[] =   9   8   6   7   5
      dvec[]  = 4.0 3.0 1.0 2.0 5.0

   created -- 98jan28, cca
   -------------------------------------------
*/
void
IV2DVqsortDown (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   dvec[]
) ;
/*
   -------------------------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in ascending order.
   ivec[] is the array of keys
   zvec[] is the companion array

   example,
      ivec[] =     5         8         3         2         9
      zvec[] = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec[] =     2         3         5         8         9
      zvec[] = (4.0,4.5) (3.0,3.5) (1.0,1.5) (2.0,2.5) (5.0,5.5)

   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
IVZVqsortUp (
   int      n,
   int      ivec[],
   double   zvec[]
) ;
/*
   -------------------------------------------------------------
   use quick sort to sort a pair or arrays.
   on output the first is in descending order.
   ivec[] is the array of keys
   zvec[] is the companion array

   example,
      ivec[] =     5         8         3         2         9
      zvec[] = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec[] =     9         8         5         3         2
      zvec[] = (5.0,5.5) (2.0,2.5) (1.0,1.5) (3.0,3.5) (4.0,4.5)

   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
IVZVqsortDown (
   int      n,
   int      ivec[],
   double   zvec[]
) ;
/*
   --------------------------------------------------------------
   use quick sort to sort a trio or arrays.
   on output the first is in ascending order.
   ivec1[] is the array of keys
   ivec2[] is a companion array
   zvec[]  is a companion array

   example,
      ivec1[] =     5         8         3         2         9
      ivec2[] =     6         7         8         9         5
      zvec[]  = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec1[] =     2         3         5         8         9
      ivec2[] =     9         8         6         7         5
      zvec[]  = (4.0,4.5) (3.0,3.5) (1.0,1.5) (2.0,2.5) (5.0,5.5)

   created -- 98jan28, cca
   --------------------------------------------------------------
*/
void
IV2ZVqsortUp (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   zvec[]
) ;
/*
   --------------------------------------------------------------
   use quick sort to sort a trio or arrays.
   on output the first is in descending order.
   ivec1[] is the array of keys
   ivec2[] is a companion array
   zvec[]  is a companion array

   example,
      ivec1[] =     5         8         3         2         9
      ivec2[] =     6         7         8         9         5
      zvec[]  = (1.0,1.5) (2.0,2.5) (3.0,3.5) (4.0,4.5) (5.0,5.5)
   becomes
      ivec1[] =     9         8         5         3         2
      ivec2[] =     5         7         6         8         9
      zvec[]  = (5.0,5.5) (2.0,2.5) (1.0,1.5) (4.0,4.5) (4.0,4.5)

   created -- 98jan28, cca
   --------------------------------------------------------------
*/
void
IV2ZVqsortDown (
   int      n,
   int      ivec1[],
   int      ivec2[],
   double   zvec[]
) ;
/*
   ---------------------------------------------------------------
   sort an array of doubles into ascending order using quick sort

   created -- 98jan28, cca
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

   created -- 98jan28, cca
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

   created -- 98jan28, cca
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

   created -- 98jan28, cca
   -------------------------------------------
*/
void
DV2qsortDown (
   int      n,
   double   dvec1[],
   double   dvec2[]
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

   created -- 98jan28, cca
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

   created -- 98jan28, cca
   -------------------------------------------
*/
void
DVIVqsortDown (
   int      n,
   double   dvec[],
   int      ivec[]
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
