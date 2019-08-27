/*  util.c  */

#include "../IV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   shift the base of the entries and adjust the dimensions
   note: this is a dangerous operation because the iv->vec
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
IV_shiftBase (
   IV    *iv,
   int    offset
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_shiftBase(%p,%d)"
           "\n bad input\n", iv, offset) ;
   exit(-1) ;
}
iv->vec     += offset ;
iv->size    -= offset ;
iv->maxsize -= offset ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   push an entry onto the list

   created -- 95oct06, cca
   modified -- 95aug28, cca
      structure has changed
   modified -- 96dec08, cca
      maxsize added to structure
   -----------------------------
*/
void
IV_push (
   IV    *iv,
   int   val
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_push(%p,%d)"
           "\n bad input\n", iv, val) ;
   exit(-1) ;
}
/*
fprintf(stdout, "\n %% iv %p, size %d, maxsize %d, entries %p",
        iv, iv->size, iv->maxsize, iv->vec) ;
fflush(stdout) ;
*/
if ( iv->size == iv->maxsize ) {
   if ( iv->maxsize == 0 ) {
      IV_setMaxsize(iv, 10) ;
   } else {
      IV_setMaxsize(iv, 2*iv->maxsize) ;
   }
}
iv->vec[iv->size++] = val ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   minimum and maximum entries

   created -- 95oct06, cca
   ---------------------------
*/
int 
IV_min ( 
   IV   *iv
) {
int   i ;

if ( iv == NULL || iv->size <= 0 || iv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in IV_min(%p), size = %d, vec = %p",
           iv, iv->size, iv->vec) ;
   exit(-1) ;
}
return(IVmin(iv->size, iv->vec, &i)) ; }

int 
IV_max ( 
   IV   *iv
) {
int   i ;

if ( iv == NULL || iv->size <= 0 || iv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in IV_max(%p), size = %d, vec = %p",
           iv, iv->size, iv->vec) ;
   exit(-1) ;
}
return(IVmax(iv->size, iv->vec, &i)) ; }

int 
IV_sum ( 
   IV   *iv
) {
int   i ;

if ( iv == NULL || iv->size <= 0 || iv->vec == NULL ) {
   fprintf(stderr, "\n fatal error in IV_sum(%p), size = %d, vec = %p",
           iv, iv->size, iv->vec) ;
   exit(-1) ;
}
return(IVsum(iv->size, iv->vec)) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   sort each index list into ascending or descending order

   created -- 95oct06, cca
   -------------------------------------------------------
*/
void
IV_sortUp ( 
   IV   *iv
) {
if ( iv == NULL || iv->size <= 0 || iv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in IV_sortUp(%p), size = %d, vec = %p",
           iv, iv->size, iv->vec) ;
   exit(-1) ;
}
IVqsortUp(iv->size, iv->vec) ;

return ; }

void
IV_sortDown ( 
   IV   *iv
) {
if ( iv == NULL || iv->size <= 0 || iv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in IV_sortDown(%p), size = %d, vec = %p",
           iv, iv->size, iv->vec) ;
   exit(-1) ;
}
IVqsortDown(iv->size, iv->vec) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   ramp the entries

   created -- 95oct06, cca
   -----------------------
*/
void
IV_ramp (
   IV    *iv,
   int   base,
   int   incr
) {
if ( iv == NULL || iv->size <= 0 || iv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in IV_ramp(%p,%d,%d), size = %d, vec = %p",
           iv, base, incr, iv->size, iv->vec) ;
   exit(-1) ;
}
IVramp(iv->size, iv->vec, base, incr) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   shuffle the list

   created -- 95oct06, cca
   -----------------------
*/
void
IV_shuffle (
   IV    *iv, 
   int   seed
) {
if ( iv == NULL || iv->size <= 0 || iv->vec == NULL ) {
   fprintf(stderr, 
           "\n fatal error in IV_shuffle(%p,%d), size = %d, vec = %p",
           iv, seed, iv->size, iv->vec) ;
   exit(-1) ;
}
IVshuffle(iv->size, iv->vec, seed) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
IV_sizeOf ( 
   IV   *iv
) {
int   nbytes ;

if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_sizeOf(%p)"
           "\n bad input \n", iv) ;
   exit(-1) ;
}
nbytes = sizeof(struct _IV) ;
if ( iv->owned == 1 ) {
   nbytes += iv->maxsize * sizeof(int) ;
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------- 
   keep entries that are tagged, move others to end and adjust size

   created -- 95oct06, cca
   ---------------------------------------------------------------- 
*/
void
IV_filterKeep (
   IV    *iv,
   int   tags[],
   int   keepTag
) {
int   i, j, w, size ;
int   *vec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || tags == NULL ) {
   fprintf(stderr, "\n fatal error in IV_filterKeep(%p,%p,%d)"
           "\n bad input", iv, tags, keepTag) ;
   exit(-1) ;
}
size = iv->size ;
vec = iv->vec ;
/*
   --------------------------------------------
   move untagged entries to the end of the list
   --------------------------------------------
*/
for ( i = 0, j = size ; i < j ;     ) {
   w = vec[i] ;
   if ( tags[w] != keepTag ) {
      vec[i]   = vec[j-1] ;
      vec[j-1] = w ;
      j-- ;
   } else {
      i++ ;
   }
}
/*
   -------------------
   reset the list size
   -------------------
*/
iv->size = j ;
 
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------- 
   move purge entries that are tagged to end and adjust size

   created -- 95oct06, cca
   --------------------------------------------------------- 
*/
void
IV_filterPurge (
   IV    *iv,
   int   tags[],
   int   purgeTag
) {
int   i, j, size, w ;
int   *vec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || tags == NULL ) {
   fprintf(stderr, "\n fatal error in IV_filterPurge(%p,%p,%d)"
           "\n bad input", iv, tags, purgeTag) ;
   exit(-1) ;
}
size = iv->size ;
vec = iv->vec ;
/*
   --------------------------------------------
   move untagged entries to the end of the list
   --------------------------------------------
*/
for ( i = 0, j = size ; i < j ;   ) {
   w = vec[i] ;
   if ( tags[w] == purgeTag ) {
      vec[i]   = vec[j-1] ;
      vec[j-1] = w ;
      j-- ;
   } else {
      i++ ;
   }
}
/*
   -------------------
   reset the list size
   -------------------
*/
iv->size = j ;
 
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   iterator :
   return the pointer to the head of the vector

   created -- 95oct06, cca
   --------------------------------------------
*/
int *
IV_first (
   IV   *iv
) {
int   *pi ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_first(%p)"
           "\n bad input", iv) ;
   exit(-1) ;
}
if ( iv->size == 0 ) {
   pi = NULL ;
} else {
   pi = iv->vec ;
}
return(pi) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   iterator :
   return the pointer to the next location in the vector

   created -- 95oct06, cca
   -----------------------------------------------------
*/
int *
IV_next (
   IV    *iv,
   int   *pi
) {
int   offset ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || pi == NULL ) {
   fprintf(stderr, "\n fatal error in IV_next(%p,%p)"
           "\n bad input", iv, pi) ;
   fflush(stderr) ;
   exit(-1) ;
}
/*
   ---------------
   check the input
   ---------------
*/
if ( (offset = pi - iv->vec) < 0 || offset >= iv->size ) {
/*
   -----------------------------
   error, offset is out of range
   -----------------------------
*/
   fprintf(stderr, "\n fatal error in IV_next(%p,%p)"
           "\n offset = %d, must be in [0,%d)",
           iv, pi, offset, iv->size) ;
   fflush(stderr) ;
   exit(-1) ;
} else if ( offset == iv->size - 1 ) {
/*
   ----------------------------
   end of the list, return NULL
   ----------------------------
*/
   pi = NULL ;
} else {
/*
   ----------------------------------------
   middle of the list, return next location
   ----------------------------------------
*/
   pi++ ;
}
return(pi) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   fill a vector with a value

   created -- 96jun22, cca
   --------------------------
*/
void
IV_fill (
   IV    *iv,
   int   value
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_fill(%p,%d)"
           "\n bad input\n", iv, value) ;
   exit(-1) ;
}
if ( iv->size > 0 ) {
   IVfill(iv->size, iv->vec, value) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   copy entries from iv2 into iv1.
   note: this is a "mapped" copy, 
   iv1 and iv2 need not be the same size.

   created -- 96aug31, cca
   --------------------------------------
*/
void
IV_copy (
   IV   *iv1,
   IV   *iv2
) {
int   ii, size ;
int   *vec1, *vec2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv1 == NULL || iv2 == NULL ) {
   fprintf(stderr, "\n fatal error in IV_copy(%p,%p)"
           "\n bad input\n", iv1, iv2) ;
   exit(-1) ;
}
size = iv1->size ;
if ( size > iv2->size ) {
   size = iv2->size ;
}
vec1 = iv1->vec ;
vec2 = iv2->vec ;
for ( ii = 0 ; ii < size ; ii++ ) {
   vec1[ii] = vec2[ii] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   increment the loc'th location of the vector by one
   and return the new value

   created -- 96dec18, cca
   --------------------------------------------------
*/
int
IV_increment (
   IV    *iv,
   int   loc
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || loc < 0 || loc >= iv->size ) {
   fprintf(stderr, "\n fatal error in IV_increment(%p,%d)"
           "\n bad input\n", iv, loc) ;
   if ( iv != NULL ) {
      fprintf(stderr, "\n loc = %d, size = %d", loc, iv->size) ;
   }
   exit(-1) ;
}
return(++iv->vec[loc]) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   decrement the loc'th location of the vector by one
   and return the new value

   created -- 96dec18, cca
   --------------------------------------------------
*/
int
IV_decrement (
   IV    *iv,
   int   loc
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || loc < 0 || loc >= iv->size ) {
   fprintf(stderr, "\n fatal error in IV_decrement(%p,%d)"
           "\n bad input\n", iv, loc) ;
   if ( iv != NULL ) {
      fprintf(stderr, "\n loc = %d, size = %d", loc, iv->size) ;
   }
   exit(-1) ;
}
return(--iv->vec[loc]) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   return the first location in the vector that contains value.
   if value is not present, -1 is returned. cost is linear in 
   the size of the vector

   created -- 96jan15, cca
   ------------------------------------------------------------
*/
int
IV_findValue (
   IV   *iv,
   int  value
) {
int   ii, n ;
int   *vec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_findValue(%p,%d)"
           "\n bad input\n", iv, value) ;
   exit(-1) ;
}
if ( (n = iv->size) <= 0 || (vec = iv->vec) == NULL ) {
   return(-1) ;
}
for ( ii = 0 ; ii < n ; ii++ ) {
   if ( vec[ii] == value ) {
      return(ii) ;
   }
} 
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   return a location in the vector that contains value.
   the entries in the vector are assumed to be in ascending order.
   if value is not present, -1 is returned. cost is logarithmic in 
   the size of the vector

   created -- 96jan15, cca
   ---------------------------------------------------------------
*/
int
IV_findValueAscending (
   IV   *iv,
   int  value
) {
int   left, mid, n, right ;
int   *vec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_findValueAscending(%p,%d)"
           "\n bad input\n", iv, value) ;
   exit(-1) ;
}
if ( (n = iv->size) <= 0 || (vec = iv->vec) == NULL ) {
   return(-1) ;
}
left  = 0 ;
right = n - 1 ;
if ( vec[left] == value ) {
   return(left) ;
} else if ( vec[right] == value ) {
   return(right) ;
} else {
   while ( left < right - 1 ) {
      mid = (left + right)/2 ;
      if ( vec[mid] == value ) {
         return(mid) ;
      } else if ( vec[mid] < value ) {
         left = mid ;
      } else {
         right = mid ;
      }
   }
}
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   return a location in the vector that contains value.
   the entries in the vector are assumed to be in descending order.
   if value is not present, -1 is returned. cost is logarithmic in 
   the size of the vector

   created -- 96jan15, cca
   ----------------------------------------------------------------
*/
int
IV_findValueDescending (
   IV   *iv,
   int  value
) {
int   left, mid, n, right ;
int   *vec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_findValueDescending(%p,%d)"
           "\n bad input\n", iv, value) ;
   exit(-1) ;
}
if ( (n = iv->size) <= 0 || (vec = iv->vec) == NULL ) {
   return(-1) ;
}
left  = 0 ;
right = n - 1 ;
if ( vec[left] == value ) {
   return(left) ;
} else if ( vec[right] == value ) {
   return(right) ;
} else {
   while ( left < right - 1 ) {
      mid = (left + right)/2 ;
      if ( vec[mid] == value ) {
         return(mid) ;
      } else if ( vec[mid] > value ) {
         left = mid ;
      } else {
         right = mid ;
      }
   }
}
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- return invlistIV, an IV object
              that contains the inverse map,
              i.e., invlist[list[ii]] = ii.
              other entries of invlist[] are -1.
              all entries in listIV must be nonnegative
 
   created -- 98aug12, cca
   ----------------------------------------------------
*/
IV *
IV_inverseMap (
   IV   *listIV
) {
int   ii, maxval, n ;
int   *invlist, *list ;
IV    *invlistIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( listIV == NULL ) {
   fprintf(stderr, "\n fatal error in IV_inverseMap()"
           "\n bad input\n") ;
   exit(-1) ;
}
IV_sizeAndEntries(listIV, &n, &list) ;
if ( n <= 0 && list == NULL ) {
   fprintf(stderr, "\n fatal error in IV_inverseMap()"
           "\n size = %d, list = %p\n", n, list) ;
   exit(-1) ;
}
for ( ii = 0, maxval = -1 ; ii < n ; ii++ ) {
   if ( list[ii] < 0 ) {
      fprintf(stderr, "\n fatal error in IV_inverseMap()"
              "\n list[%d] = %d, must be positive\n", ii, list[ii]) ;
      exit(-1) ;
   }
   if ( maxval < list[ii] ) {
      maxval = list[ii] ;
   }
}
invlistIV = IV_new() ;
IV_init(invlistIV, 1 + maxval, NULL) ;
IV_fill(invlistIV, -1) ;
invlist = IV_entries(invlistIV) ;
for ( ii = 0 ; ii < n ; ii++ ) {
   if ( invlist[list[ii]] != -1 ) {
      fprintf(stderr, "\n fatal error in IV_inverseMap()"
              "\n repeated list value %d\n", list[ii]) ;
      exit(-1) ;
   }
   invlist[list[ii]] = ii ;
}
return(invlistIV) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- return an IV object that contains the locations 
              of all instances of target in listIV. this is 
              used when listIV is a map from [0,n) to a finite
              set (like processors) and we want to find all
              entries that are mapped to a specific processor.

   created -- 98aug12, cca
   -----------------------------------------------------------
*/
IV *
IV_targetEntries (
   IV   *listIV,
   int  target
) {
int   count, ii, n ;
int   *entries, *list ;
IV    *entriesIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( listIV == NULL ) {
   fprintf(stderr, "\n fatal error in IV_targetEntries()"
           "\n bad input\n") ;
   exit(-1) ;
}
IV_sizeAndEntries(listIV, &n, &list) ;
if ( n <= 0 && list == NULL ) {
   fprintf(stderr, "\n fatal error in IV_targetEntries()"
           "\n size = %d, list = %p\n", n, list) ;
   exit(-1) ;
}
for ( ii = count = 0 ; ii < n ; ii++ ) {
   if ( list[ii] == target ) {
      count++ ;
   }
}
entriesIV = IV_new() ;
if ( count > 0 ) {
   IV_init(entriesIV, count, NULL) ;
   entries = IV_entries(entriesIV) ;
   for ( ii = count = 0 ; ii < n ; ii++ ) {
      if ( list[ii] == target ) {
         entries[count++] = ii ;
      }
   }
}
return(entriesIV) ; }

/*--------------------------------------------------------------------*/
