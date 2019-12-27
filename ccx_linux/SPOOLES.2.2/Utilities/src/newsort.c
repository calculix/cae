/*  IVsort.c  */

#include "../Utilities.h"

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
   ---------------------------------------
   macros to make swapping elements easier
   ---------------------------------------
*/
#define I_SWAP(a,b) { \
      itemp     = ivec[(a)] ; \
      ivec[(a)] = ivec[(b)] ; \
      ivec[(b)] = itemp     ; }
#define I2_SWAP(a,b) { \
      itemp      = ivec1[(a)] ; \
      ivec1[(a)] = ivec1[(b)] ; \
      ivec1[(b)] = itemp      ; \
      itemp      = ivec2[(a)] ; \
      ivec2[(a)] = ivec2[(b)] ; \
      ivec2[(b)] = itemp      ; }
#define ID_SWAP(a,b) { \
      itemp     = ivec[(a)] ; \
      ivec[(a)] = ivec[(b)] ; \
      ivec[(b)] = itemp     ; \
      dtemp     = dvec[(a)] ; \
      dvec[(a)] = dvec[(b)] ; \
      dvec[(b)] = dtemp     ; }
#define I2D_SWAP(a,b) { \
      itemp      = ivec1[(a)] ; \
      ivec1[(a)] = ivec1[(b)] ; \
      ivec1[(b)] = itemp      ; \
      itemp      = ivec2[(a)] ; \
      ivec2[(a)] = ivec2[(b)] ; \
      ivec2[(b)] = itemp      ; \
      dtemp      = dvec[(a)]  ; \
      dvec[(a)]  = dvec[(b)]  ; \
      dvec[(b)]  = dtemp      ; }
#define IZ_SWAP(a,b) { \
      itemp         = ivec[(a)]     ; \
      ivec[(a)]     = ivec[(b)]     ; \
      ivec[(b)]     = itemp         ; \
      dtemp         = zvec[2*(a)]   ; \
      zvec[2*(a)]   = zvec[2*(b)]   ; \
      zvec[2*(b)]   = dtemp         ; \
      dtemp         = zvec[2*(a)+1] ; \
      zvec[2*(a)+1] = zvec[2*(b)+1] ; \
      zvec[2*(b)+1] = dtemp         ; }
#define I2Z_SWAP(a,b) { \
      itemp         = ivec1[(a)]    ; \
      ivec1[(a)]    = ivec1[(b)]    ; \
      ivec1[(b)]    = itemp         ; \
      itemp         = ivec2[(a)]    ; \
      ivec2[(a)]    = ivec2[(b)]    ; \
      ivec2[(b)]    = itemp         ; \
      dtemp         = zvec[2*(a)]   ; \
      zvec[2*(a)]   = zvec[2*(b)]   ; \
      zvec[2*(b)]   = dtemp         ; \
      dtemp         = zvec[2*(a)+1] ; \
      zvec[2*(a)+1] = zvec[2*(b)+1] ; \
      zvec[2*(b)+1] = dtemp         ; \
}
#define D_SWAP(a,b) { \
      dtemp     = dvec[(a)] ; \
      dvec[(a)] = dvec[(b)] ; \
      dvec[(b)] = dtemp     ; }
#define D2_SWAP(a,b) { \
      dtemp      = dvec1[(a)] ; \
      dvec1[(a)] = dvec1[(b)] ; \
      dvec1[(b)] = dtemp      ; \
      dtemp      = dvec2[(a)] ; \
      dvec2[(a)] = dvec2[(b)] ; \
      dvec2[(b)] = dtemp      ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   static methods to choose partitions
   -----------------------------------
*/
/*
   ----------------------------------------------------
   static function, return the median of three integers
   ----------------------------------------------------
*/
static int
Imedian3 ( 
   int   i,
   int   j,
   int   k,
   int   a[]
) {
if ( a[i] < a[j] ) {    
   /* a[i] < a[j] */
   if ( a[j] < a[k] ) { 
      /* a[i] < a[j] < a[k] */
      return(j) ;
   } else if ( a[i] < a[k] ) {
      /* a[i] < a[k] <= a[j] */
      return(k) ;
   } else {
      /* a[k] <= a[i] < a[j] */
      return(i) ;
   }
} else {
   /* a[j] <= a[i] */
   if ( a[i] < a[k] ) {
      /*  a[j] <= a[i] < a[k] */
      return(i) ;
   } else if ( a[j] < a[k] ) {
      /*  a[j] < a[k] <= a[i] */
      return(k) ;
   } else {
      /*  a[k] <= a[j] <= a[i] */
      return(j) ;
   }
}
}
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   static function, 
   returns an approximation to the median value of a vector
   if n < 7 then
      returns a[n/2]
   else if n < 40 then
      returns median(a[0], a[n/2], a[n-1])
   else 
      returns median( median(a[0], a[s], a[2s])
                      median(a[n/2-s], a[n/2], a[n/2+s])
                      median(a[n-1-2s], a[n-1-s], a[n-1]) )
      where s = n / 8
   endif

   created -- 98jan28, cca
   ---------------------------------------------------------
*/
static int
Icentervalue (
   int   n,
   int   a[]
) {
int   i, j, k, s ;

j = n / 2 ;
if ( n > 7 ) {
   i = 0 ;
   k = n - 1 ;
   if ( n >= 40 ) {
      s = n / 8 ;
      i = Imedian3(i, i+s, i+s+s, a) ;
      j = Imedian3(j-s, j, j+s, a) ;
      k = Imedian3(k-s-s, k-s, k, a) ;
   }
   j = Imedian3(i, j, k, a) ;
}
return(a[j]) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   static function, return the median of three doubles
   ----------------------------------------------------
*/
static int
Dmedian3 ( 
   int      i,
   int      j,
   int      k,
   double   a[]
) {
if ( a[i] < a[j] ) {    
   /* a[i] < a[j] */
   if ( a[j] < a[k] ) { 
      /* a[i] < a[j] < a[k] */
      return(j) ;
   } else if ( a[i] < a[k] ) {
      /* a[i] < a[k] <= a[j] */
      return(k) ;
   } else {
      /* a[k] <= a[i] < a[j] */
      return(i) ;
   }
} else {
   /* a[j] <= a[i] */
   if ( a[i] < a[k] ) {
      /*  a[j] <= a[i] < a[k] */
      return(i) ;
   } else if ( a[j] < a[k] ) {
      /*  a[j] < a[k] <= a[i] */
      return(k) ;
   } else {
      /*  a[k] <= a[j] <= a[i] */
      return(j) ;
   }
}
}
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   static function, 
   returns an approximation to the median value of a vector
   if n < 7 then
      returns a[n/2]
   else if n < 40 then
      returns median(a[0], a[n/2], a[n-1])
   else 
      returns median( median(a[0], a[s], a[2s])
                      median(a[n/2-s], a[n/2], a[n/2+s])
                      median(a[n-1-2s], a[n-1-s], a[n-1]) )
      where s = n / 8
   endif

   created -- 98jan28, cca
   ---------------------------------------------------------
*/
static double
Dcentervalue (
   int      n,
   double   a[]
) {
int   i, j, k, s ;

j = n / 2 ;
if ( n > 7 ) {
   i = 0 ;
   k = n - 1 ;
   if ( n >= 40 ) {
      s = n / 8 ;
      i = Dmedian3(i, i+s, i+s+s, a) ;
      j = Dmedian3(j-s, j, j+s, a) ;
      k = Dmedian3(k-s-s, k-s, k, a) ;
   }
   j = Dmedian3(i, j, k, a) ;
}
return(a[j]) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   methods to check if the int vector 
   is in ascending or descending order
   -----------------------------------
*/
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
) {
if ( n <= 0 ) {
   return(0) ;
} else if ( n == 1 ) {
   return(1) ;
} else {
   int   i ;
   for ( i = 1 ; i < n ; i++ ) {
      if ( ivec[i-1] > ivec[i] ) {
         return(0) ;
      }
   }
   return(1) ;
}
}
/*--------------------------------------------------------------------*/
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
) {
if ( n <= 0 ) {
   return(0) ;
} else if ( n == 1 ) {
   return(1) ;
} else {
   int   i ;
   for ( i = 1 ; i < n ; i++ ) {
      if ( ivec[i-1] < ivec[i] ) {
         return(0) ;
      }
   }
   return(1) ;
}
}
/*--------------------------------------------------------------------*/
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
) {
if ( n <= 0 ) {
   return(0) ;
} else if ( n == 1 ) {
   return(1) ;
} else {
   int   i ;
   for ( i = 1 ; i < n ; i++ ) {
      if ( dvec[i-1] > dvec[i] ) {
         return(0) ;
      }
   }
   return(1) ;
}
}
/*--------------------------------------------------------------------*/
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
) {
if ( n <= 0 ) {
   return(0) ;
} else if ( n == 1 ) {
   return(1) ;
} else {
   int   i ;
   for ( i = 1 ; i < n ; i++ ) {
      if ( dvec[i-1] < dvec[i] ) {
         return(0) ;
      }
   }
   return(1) ;
}
}
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
) {
int   i, itemp, j ;

for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec[j-1] > ivec[j] ; j-- ) {
      I_SWAP(j-1,j) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
int   i, itemp, j ;

for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec[j-1] < ivec[j] ; j-- ) {
      I_SWAP(j-1,j) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
int   i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec1[j-1] > ivec1[j] ; j-- ) {
      I2_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
int   i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec1[j-1] < ivec1[j] ; j-- ) {
      I2_SWAP(j-1,j)
   }
}
return ; }

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
 
   created -- 98jan28, cca
   ------------------------------------------
*/
void
IVDVisortUp (
   int      n,
   int      ivec[],
   double   dvec[]
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec[j-1] > ivec[j] ; j-- ) {
      ID_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec[j-1] < ivec[j] ; j-- ) {
      ID_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec1[j-1] > ivec1[j] ; j-- ) {
      I2D_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec1[j-1] < ivec1[j] ; j-- ) {
      I2D_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec[j-1] > ivec[j] ; j-- ) {
      IZ_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec[j-1] < ivec[j] ; j-- ) {
      IZ_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec1[j-1] > ivec1[j] ; j-- ) {
      I2Z_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && ivec1[j-1] < ivec1[j] ; j-- ) {
      I2Z_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && dvec[j-1] > dvec[j] ; j-- ) {
      D_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && dvec[j-1] < dvec[j] ; j-- ) {
      D_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && dvec1[j-1] > dvec1[j] ; j-- ) {
      D2_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && dvec1[j-1] < dvec1[j] ; j-- ) {
      D2_SWAP(j-1,j)
   }
}
return ; }

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
 
   created -- 98jan28, cca
   ------------------------------------------
*/
void
DVIVisortUp (
   int      n,
   double   dvec[],
   int      ivec[]
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && dvec[j-1] > dvec[j] ; j-- ) {
      ID_SWAP(j-1,j)
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      i, itemp, j ;
 
for ( i = 1 ; i < n ; i++ ) {
   for ( j = i ; j > 0 && dvec[j-1] < dvec[j] ; j-- ) {
      ID_SWAP(j-1,j)
   }
}
return ; }

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
) {
int   a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IVisortUp(n, ivec) ;
} else {
   v = Icentervalue(n, ivec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec[b] <= v ) {
         if ( ivec[b] == v ) {
            I_SWAP(a,b) ;
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec[c] >= v ) {
         if ( ivec[c] == v ) {
            I_SWAP(c,d) ;
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I_SWAP(b,c) ;
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I_SWAP(l,h) ;
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I_SWAP(l,h) ;
      l++ ;
      h++ ;
   }
   IVqsortUp(b - a, ivec) ;
   IVqsortUp(d - c, ivec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
int   a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IVisortDown(n, ivec) ;
} else {
   v = Icentervalue(n, ivec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec[b] >= v ) {
         if ( ivec[b] == v ) {
            I_SWAP(a,b) ;
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec[c] <= v ) {
         if ( ivec[c] == v ) {
            I_SWAP(c,d) ;
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I_SWAP(b,c) ;
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I_SWAP(l,h) ;
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I_SWAP(l,h) ;
      l++ ;
      h++ ;
   }
   IVqsortDown(b - a, ivec) ;
   IVqsortDown(d - c, ivec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
int   a, b, c, d, itemp, l, h, s, v ;
 
if ( n <= 10 ) {
   IV2isortUp(n, ivec1, ivec2) ;
} else {
   v = Icentervalue(n, ivec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec1[b] <= v ) {
         if ( ivec1[b] == v ) {
            I2_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec1[c] >= v ) {
         if ( ivec1[c] == v ) {
            I2_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I2_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I2_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I2_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IV2qsortUp(b - a, ivec1, ivec2) ;
   IV2qsortUp(d - c, ivec1 + n - (d - c), ivec2 + n - (d - c)) ;
}
 
return ; }
 
/*--------------------------------------------------------------------*/
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
) {
int   a, b, c, d, itemp, l, h, s, v ;
 
if ( n <= 10 ) {
   IV2isortDown(n, ivec1, ivec2) ;
} else {
   v = Icentervalue(n, ivec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec1[b] >= v ) {
         if ( ivec1[b] == v ) {
            I2_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec1[c] <= v ) {
         if ( ivec1[c] == v ) {
            I2_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I2_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I2_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I2_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IV2qsortDown(b - a, ivec1, ivec2) ;
   IV2qsortDown(d - c, ivec1 + n - (d - c), ivec2 + n - (d - c)) ;
}
 
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IVDVisortUp(n, ivec, dvec) ;
} else {
   v = Icentervalue(n, ivec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec[b] <= v ) {
         if ( ivec[b] == v ) {
            ID_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec[c] >= v ) {
         if ( ivec[c] == v ) {
            ID_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      ID_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IVDVqsortUp(b - a, ivec, dvec) ;
   IVDVqsortUp(d - c, ivec + n - (d - c), dvec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IVDVisortDown(n, ivec, dvec) ;
} else {
   v = Icentervalue(n, ivec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec[b] >= v ) {
         if ( ivec[b] == v ) {
            ID_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec[c] <= v ) {
         if ( ivec[c] == v ) {
            ID_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      ID_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IVDVqsortDown(b - a, ivec, dvec) ;
   IVDVqsortDown(d - c, ivec + n - (d - c), dvec + n - (d - c)) ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IV2DVisortUp(n, ivec1, ivec2, dvec) ;
} else {
   v = Icentervalue(n, ivec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec1[b] <= v ) {
         if ( ivec1[b] == v ) {
            I2D_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec1[c] >= v ) {
         if ( ivec1[c] == v ) {
            I2D_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I2D_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I2D_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I2D_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IV2DVqsortUp(b - a, ivec1, ivec2, dvec) ;
   IV2DVqsortUp(d - c, ivec1 + n - (d - c), ivec2 + n - (d - c),
                dvec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IV2DVisortDown(n, ivec1, ivec2, dvec) ;
} else {
   v = Icentervalue(n, ivec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec1[b] >= v ) {
         if ( ivec1[b] == v ) {
            I2D_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec1[c] <= v ) {
         if ( ivec1[c] == v ) {
            I2D_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I2D_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I2D_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I2D_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IV2DVqsortDown(b - a, ivec1, ivec2, dvec) ;
   IV2DVqsortDown(d - c, ivec1 + n - (d - c), ivec2 + n - (d - c),
                  dvec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IVZVisortUp(n, ivec, zvec) ;
} else {
   v = Icentervalue(n, ivec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec[b] <= v ) {
         if ( ivec[b] == v ) {
            IZ_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec[c] >= v ) {
         if ( ivec[c] == v ) {
            IZ_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      IZ_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      IZ_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      IZ_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IVZVqsortUp(b - a, ivec, zvec) ;
   IVZVqsortUp(d - c, ivec + n - (d - c), zvec + 2*(n - (d - c))) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IVZVisortDown(n, ivec, zvec) ;
} else {
   v = Icentervalue(n, ivec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec[b] >= v ) {
         if ( ivec[b] == v ) {
            IZ_SWAP(a,b)
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec[c] <= v ) {
         if ( ivec[c] == v ) {
            IZ_SWAP(c,d)
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      IZ_SWAP(b,c)
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      IZ_SWAP(l,h)
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      IZ_SWAP(l,h)
      l++ ;
      h++ ;
   }
   IVZVqsortDown(b - a, ivec, zvec) ;
   IVZVqsortDown(d - c, ivec + n - (d - c), zvec + 2*(n - (d - c))) ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 10 ) {
   IV2ZVisortUp(n, ivec1, ivec2, zvec) ;
} else {
   v = Icentervalue(n, ivec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec1[b] <= v ) {
         if ( ivec1[b] == v ) {
            I2Z_SWAP(a,b) 
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec1[c] >= v ) {
         if ( ivec1[c] == v ) {
            I2Z_SWAP(c,d) 
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I2Z_SWAP(b,c) 
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I2Z_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I2Z_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   IV2ZVqsortUp(b - a, ivec1, ivec2, zvec) ;
   IV2ZVqsortUp(d - c, ivec1 + n - (d - c), ivec2 + n - (d - c),
                zvec + 2*(n - (d - c))) ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp ;
int      a, b, c, d, itemp, l, h, s, v ;

if ( n <= 0 ) {
   IV2ZVisortDown(n, ivec1, ivec2, zvec) ;
} else {
   v = Icentervalue(n, ivec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && ivec1[b] >= v ) {
         if ( ivec1[b] == v ) {
            I2Z_SWAP(a,b) 
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && ivec1[c] <= v ) {
         if ( ivec1[c] == v ) {
            I2Z_SWAP(c,d) 
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      I2Z_SWAP(b,c) 
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      I2Z_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      I2Z_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   IV2ZVqsortDown(b - a, ivec1, ivec2, zvec) ;
   IV2ZVqsortDown(d - c, ivec1 + n - (d - c), ivec2 + n - (d - c),
                  zvec + 2*(n - (d - c))) ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp, v ;
int      a, b, c, d, l, h, s ;

if ( n <= 10 ) {
   DVisortUp(n, dvec) ;
} else {
   v = Dcentervalue(n, dvec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && dvec[b] <= v ) {
         if ( dvec[b] == v ) {
            D_SWAP(a,b) 
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && dvec[c] >= v ) {
         if ( dvec[c] == v ) {
            D_SWAP(c,d) 
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      D_SWAP(b,c) 
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      D_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      D_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   DVqsortUp(b - a, dvec) ;
   DVqsortUp(d - c, dvec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp, v ;
int      a, b, c, d, l, h, s ;

if ( n <= 10 ) {
   DVisortDown(n, dvec) ;
} else {
   v = Dcentervalue(n, dvec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && dvec[b] >= v ) {
         if ( dvec[b] == v ) {
            D_SWAP(a,b) 
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && dvec[c] <= v ) {
         if ( dvec[c] == v ) {
            D_SWAP(c,d) 
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      D_SWAP(b,c) 
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      D_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      D_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   DVqsortDown(b - a, dvec) ;
   DVqsortDown(d - c, dvec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp, v ;
int      a, b, c, d, l, h, s ;

if ( n <= 10 ) {
   DV2isortUp(n, dvec1, dvec2) ;
} else {
   v = Dcentervalue(n, dvec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && dvec1[b] <= v ) {
         if ( dvec1[b] == v ) {
            D2_SWAP(a,b) 
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && dvec1[c] >= v ) {
         if ( dvec1[c] == v ) {
            D2_SWAP(c,d) 
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      D2_SWAP(b,c) 
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      D2_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      D2_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   DV2qsortUp(b - a, dvec1, dvec2) ;
   DV2qsortUp(d - c, dvec1 + n - (d - c), dvec2 + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp, v ;
int      a, b, c, d, l, h, s ;

if ( n <= 10 ) {
   DV2isortDown(n, dvec1, dvec2) ;
} else {
   v = Dcentervalue(n, dvec1) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && dvec1[b] >= v ) {
         if ( dvec1[b] == v ) {
            D2_SWAP(a,b) 
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && dvec1[c] <= v ) {
         if ( dvec1[c] == v ) {
            D2_SWAP(c,d) 
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      D2_SWAP(b,c) 
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      D2_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      D2_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   DV2qsortDown(b - a, dvec1, dvec2) ;
   DV2qsortDown(d - c, dvec1 + n - (d - c), dvec2 + n - (d - c)) ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp, v ;
int      a, b, c, d, itemp, l, h, s ;

if ( n <= 10 ) {
   DVIVisortUp(n, dvec, ivec) ;
} else {
   v = Dcentervalue(n, dvec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && dvec[b] <= v ) {
         if ( dvec[b] == v ) {
            ID_SWAP(a,b) ;
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && dvec[c] >= v ) {
         if ( dvec[c] == v ) {
            ID_SWAP(c,d) ;
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      ID_SWAP(b,c) ;
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h) ;
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h) ;
      l++ ;
      h++ ;
   }
   DVIVqsortUp(b - a, dvec, ivec) ;
   DVIVqsortUp(d - c, dvec + n - (d - c), ivec + n - (d - c)) ;
}

return ; }

/*--------------------------------------------------------------------*/
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
) {
double   dtemp, v ;
int      a, b, c, d, itemp, l, h, s ;

if ( n <= 10 ) {
   DVIVisortDown(n, dvec, ivec) ;
} else {
   v = Dcentervalue(n, dvec) ;
   a = b = 0 ;
   c = d = n - 1 ;
   for ( ; ; ) {
      while ( b <= c && dvec[b] >= v ) {
         if ( dvec[b] == v ) {
            ID_SWAP(a,b) 
            a++ ;
         }
         b++ ;
      }
      while ( c >= b && dvec[c] <= v ) {
         if ( dvec[c] == v ) {
            ID_SWAP(c,d) 
            d-- ;
         }
         c-- ;
      }
      if ( b > c ) {
         break ;
      }
      ID_SWAP(b,c) 
      b++ ;
      c-- ;
   }
   s = (a <= b - a) ? a : b - a ;
   for ( l = 0, h = b - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   s = ((d - c) <= (n - 1 - d)) ? (d - c) : ( n - 1 - d) ;
   for ( l = b, h = n - s ; s != 0 ; s-- ) {
      ID_SWAP(l,h) 
      l++ ;
      h++ ;
   }
   DVIVqsortDown(b - a, dvec, ivec) ;
   DVIVqsortDown(d - c, dvec + n - (d - c), ivec + n - (d - c)) ;
}
return ; }

/*--------------------------------------------------------------------*/
