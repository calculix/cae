/*  basics.c  */

#include "../IIheap.h"

/*--------------------------------------------------------------------*/
static void IIheap_siftUp ( IIheap *heap, int loc ) ;
static void IIheap_siftDown ( IIheap *heap, int loc ) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   create and return a new instance of the IIheap object

   created -- 95sep30, cca
   -----------------------------------------------------
*/
IIheap *
IIheap_new (
   void
) {
IIheap   *heap ;

ALLOCATE(heap, struct _IIheap, 1) ;

IIheap_setDefaultFields(heap) ;

return(heap) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   set the default fields for an IIheap object

   created -- 95sep30, cca
   -------------------------------------------
*/
void
IIheap_setDefaultFields (
   IIheap   *heap
) {
if ( heap == NULL ) {
   fprintf(stderr, "\n fatal error in IIheap_setDefaultFields(%p)"
           "\n heap is NULL\n", heap) ;
   exit(-1) ;
}
heap->size    =   0  ;
heap->maxsize =   0  ;
heap->heapLoc = NULL ;
heap->keys    = NULL ;
heap->values  = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 95sep30, cca
   -----------------------
*/
void
IIheap_clearData (
   IIheap   *heap
) {
if ( heap == NULL ) {
   fprintf(stderr, "\n fatal error in IIheap_clearData(%p)"
           "\n heap is NULL\n", heap) ;
   exit(-1) ;
}
if ( heap->heapLoc != NULL ) {
   IVfree(heap->heapLoc) ;
}
if ( heap->keys != NULL ) {
   IVfree(heap->keys) ;
}
if ( heap->values != NULL ) {
   IVfree(heap->values) ;
}
IIheap_setDefaultFields(heap) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   free the IIheap object

   created -- 95sep30, cca
   -----------------------
*/
void
IIheap_free (
   IIheap   *heap
) {
if ( heap == NULL ) {
   fprintf(stderr, "\n fatal error in IIheap_free(%p)"
           "\n heap is NULL\n", heap) ;
   exit(-1) ;
}
IIheap_clearData(heap) ;
FREE(heap) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   initializer, 
   set heap maximum size to maxsize
   and allocate the arrays

   created -- 95sep30, cca
   --------------------------------
*/
void
IIheap_init ( 
   IIheap   *heap,
   int      maxsize 
) {
if ( heap == NULL || maxsize <= 0 ) {
   fprintf(stderr, "\n\n error in IIheap_init(%p,%d)"
           "\n heap is NULL or size = %d is nonpositive\n",
           heap, maxsize, maxsize) ;
   exit(-1) ;
}
IIheap_clearData(heap) ;
heap->maxsize = maxsize ;
heap->heapLoc = IVinit(maxsize, -1) ;
heap->keys    = IVinit(maxsize, -1) ;
heap->values  = IVinit(maxsize, -1) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   fill pkey with the key and pvalue with the value 
   of the minimum (key, value) pair in the heap

   created -- 95sep30, cca
   ------------------------------------------------
*/
void
IIheap_root ( 
   IIheap   *heap,
   int      *pkey, 
   int      *pvalue 
) {
if ( heap == NULL || pkey == NULL || pvalue == NULL ) {
   fprintf(stderr, "\n\n error in IIheap_root(%p,%p,%p)"
           "\n heap is NULL or pid is NULL or pkey is NULL\n",
           heap, pkey, pvalue) ;
   exit(-1) ;
}

if ( heap->size > 0 ) {
   *pkey  = heap->keys[0]   ;
   *pvalue = heap->values[0] ;
} else {
   *pkey = *pvalue = -1 ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   insert the (key, value) pair into the heap

   created -- 95sep30, cca
   ------------------------------------------
*/
void
IIheap_insert ( 
   IIheap   *heap,
   int      key, 
   int      value 
) {
int   loc ;

if ( heap == NULL || key < 0 || key >= heap->maxsize ) {
   fprintf(stderr, "\n error in IIheap_insert(%p,%d,%d)"
           "\n heap is NULL or pair (%d,%d) is out of bounds\n",
           heap, key, value, key, value) ;
   exit(-1) ;
}
if ( heap->heapLoc[key] != -1 ) {
   fprintf(stderr, "\n error in IIheap_insert(%p,%d,%d)"
           "\n object (%d,%d) is already in heap\n",
           heap, key, value, key, value) ;
   exit(-1) ;
}
if ( heap->size == heap->maxsize ) {
   fprintf(stderr, "\n error in IIheap_insert(%p,%d,%d)"
           "\n heap size exceeded\n", heap, key, value) ;
   exit(-1) ;
}
/*
   -----------------------------
   insert the object and sift up
   -----------------------------
*/
heap->heapLoc[key] = loc = heap->size++ ;
heap->keys[loc]    = key   ;
heap->values[loc]  = value ;
IIheap_siftUp(heap, loc) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   remove (key,*) from the heap

   created -- 95sep30, cca
   ----------------------------
*/
void
IIheap_remove ( 
   IIheap   *heap,
   int      key 
) {
int   last, loc, newkey, newVal, oldVal ;
/*
   ---------------
   check the input
   ---------------
*/
if ( heap == NULL || key < 0 || key >= heap->maxsize ) {
   fprintf(stderr, "\n error in IIheap_remove(%p,%d)"
           "\n heap is NULL or object (%d) is out of bounds\n",
           heap, key, key) ;
   exit(-1) ;
}
if ( (loc = heap->heapLoc[key]) == -1 ) {
   fprintf(stderr, "\n error in IIheap_remove(%p,%d)"
           "\n object %d not in heap", heap, key, key) ;
   exit(-1) ;
}
/*
   -----------------------------------
   save the old value at this location
   and update the three heap vectors
   -----------------------------------
*/
if ( loc == (last = --heap->size) ) {
/*
   -----------------------------------------------------------------
   the element occupied the last position in the heap, simple return
   -----------------------------------------------------------------
*/
   heap->heapLoc[key] = -1 ;
   heap->keys[loc]    = -1 ;
   heap->values[loc]  =  0 ;
} else {
/*
   ---------------------------------------------
   move the last element to the current location
   ---------------------------------------------
*/
   heap->heapLoc[key]    = -1 ;
   newkey                = heap->keys[last] ;
   heap->heapLoc[newkey] = loc ;
   heap->keys[loc]       = newkey ;
   heap->keys[last]      = -1 ;
   oldVal                = heap->values[loc] ;
   heap->values[loc]     = newVal = heap->values[last] ;
   heap->values[last]    =  0 ;
   if ( oldVal > newVal ) {
      IIheap_siftUp(heap, loc) ;
   } else if ( oldVal < newVal ) {
      IIheap_siftDown(heap, loc) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to write the IIheap object to a file

   created -- 95sep30, cca
   -----------------------------------------------
*/
void
IIheap_print ( 
   IIheap   *heap,
   FILE     *fp 
) {
int   ierr, j ;

if ( heap == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in IIheap_print(%p,%p)"
           "\n heap is NULL or file pointer is NULL", heap, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n\n IIheap : present size %d, max size %d",
        heap->size, heap->maxsize) ;
if ( heap->size > 0 ) {
   fprintf(fp, "\n contents of heap : (location id value)") ;
   for ( j = 0 ; j < heap->size ; j++ ) {
      fprintf(fp, "\n %8d %8d %8d", j, 
              heap->keys[j], heap->values[j]) ;
   }
   fprintf(fp, "\n locations of ids") ;
   IVfp80(fp, heap->maxsize, heap->heapLoc, 80, &ierr) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95sep30, cca
   ----------------------------------------------
*/
int
IIheap_sizeOf ( 
   IIheap   *heap
) {
return(sizeof(IIheap) + 3*heap->maxsize*sizeof(int)) ; }

/*====================================================================*/
/*
   -----------------------------------------------
   purpose -- to sift a heap element down the tree

   created -- 95sep30, cca
   -----------------------------------------------
*/
static void
IIheap_siftDown ( 
   IIheap   *heap,
   int      loc 
) {
int   desc, itemp, left, right, size ;
int   *heapLoc, *keys, *values ;
/*
   ---------------
   check the input
   ---------------
*/
if ( heap == NULL || loc < 0 || loc >= heap->size ) {
   fprintf(stderr, "\n fatal error in IIheap_siftDown(%p,%d)"
           "\n heap is NULL or loc = %d out of range\n",
           heap, loc, loc) ;
   exit(-1) ;
}
size    = heap->size    ;
heapLoc = heap->heapLoc ;
keys    = heap->keys  ;
values  = heap->values ;

while ( 1 ) {
   left  = 2 * loc + 1 ;
   right = left + 1 ;
   if ( left >= size ) {
      break ;
   } else if ( right >= size ) {
      desc = left ;
   } else {
      if ( values[left] <= values[right] ) {
         desc = left ;
      } else {
         desc = right ;
      }
   }
/*
   if ( values[desc] < values[loc] ) {
*/
   if ( values[desc] <= values[loc] ) {
      itemp        = values[desc] ;
      values[desc] = values[loc]  ;
      values[loc]  = itemp        ;
      itemp        = keys[desc]   ;
      keys[desc]   = keys[loc]    ;
      keys[loc]    = itemp        ;
      heapLoc[keys[loc]]  = loc   ;
      heapLoc[keys[desc]] = desc  ;
      loc = desc ;
   } else {
      break ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to sift a heap element up the tree

   created -- 95sep30, cca
   ---------------------------------------------
*/
static void
IIheap_siftUp ( 
   IIheap   *heap,
   int      loc 
) {
int   itemp, par ;
int   *heapLoc, *keys, *values ;
/*
   ---------------
   check the input
   ---------------
*/
if ( heap == NULL || loc < 0 || loc >= heap->size ) {
   fprintf(stderr, "\n fatal error in IIheap_siftUp(%p,%d)"
           "\n heap is NULL or loc = %d out of range\n",
           heap, loc, loc) ;
   exit(-1) ;
}
heapLoc = heap->heapLoc ;
keys    = heap->keys  ;
values = heap->values ;

while ( loc != 0 ) {
   par = (loc - 1)/2 ;
/*
   if ( values[par] > values[loc] ) {
*/
   if ( values[par] >= values[loc] ) {
      itemp       = values[par] ;
      values[par] = values[loc] ;
      values[loc] = itemp       ;
      itemp       = keys[par]   ;
      keys[par]   = keys[loc]   ;
      keys[loc]   = itemp       ;
      heapLoc[keys[loc]] = loc  ;
      heapLoc[keys[par]] = par  ;
      loc = par ;
   } else {
      break ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
