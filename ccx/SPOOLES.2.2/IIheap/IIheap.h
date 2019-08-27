/*  IIheap.h  */

#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   IIheap -- heap with integer ids, integer values and a map from
             each integer value to a heap location

   size    -- present size of the heap
   maxsize -- maximum size of the heap
   heapLoc -- heap location of each id, size maxsize
   keys    -- object key of each location in the heap, size maxsize
   values  -- object value of each location in the heap, size maxsize

   created -- 95sep30, cca
   ------------------------------------------------------------------
*/
typedef struct _IIheap   IIheap ;
struct _IIheap {
   int    size     ;
   int    maxsize  ;
   int    *heapLoc ;
   int    *keys    ;
   int    *values  ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
--- methods found in basics.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   create and return a new instance of the IIheap object
 
   created -- 95sep30, cca
   -----------------------------------------------------
*/
IIheap *
IIheap_new (
   void
) ;
/*
   -------------------------------------------
   set the default fields for an IIheap object
 
   created -- 95sep30, cca
   -------------------------------------------
*/
void
IIheap_setDefaultFields (
   IIheap   *heap
) ;
/*
   -----------------------
   clear the data fields
 
   created -- 95sep30, cca
   -----------------------
*/
void
IIheap_clearData (
   IIheap   *heap
) ;
/*
   -----------------------
   free the IIheap object
 
   created -- 95sep30, cca
   -----------------------
*/
void
IIheap_free (
   IIheap   *heap
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
/*
   ----------------------------------------------
   return the number of bytes taken by the object
 
   created -- 95sep30, cca
   ----------------------------------------------
*/
int
IIheap_sizeOf (
   IIheap   *heap
) ;
/*--------------------------------------------------------------------*/
