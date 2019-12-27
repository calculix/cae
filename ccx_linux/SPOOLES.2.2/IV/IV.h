/*  IV.h  */

#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   IV -- integer vector object
 
   size    -- size of the vector
   maxsize -- maximum size of the vector
   owned   -- owner flag
      when == 1, storage pointed to by entries
      has been allocated here and can be free'd.
      when == 0, storage pointed to by entries
      has not been allocated here and cannot be free'd.
   vec -- pointer to base address
   ---------------------------------------------------------
*/
typedef struct _IV   IV ;
struct _IV {
   int   size    ;
   int   maxsize ;
   int   owned   ;
   int   *vec    ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   constructor method

   created -- 95oct06, cca
   -----------------------
*/
IV *
IV_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields

   created -- 95oct06, cca
   -----------------------
*/
void
IV_setDefaultFields ( 
   IV   *iv
) ;
/*
   -----------------------
   clear the data fields

   created -- 95oct06, cca
   -----------------------
*/
void
IV_clearData ( 
   IV   *iv
) ;
/*
   -----------------------
   destructor

   created -- 95oct06, cca
   -----------------------
*/
void
IV_free ( 
   IV   *iv
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------
   simplest initialization method
 
   if entries != NULL
      the object does not own the entries,
      it just points to the entries base address
   else if size > 0
      the object will own the entries,
      it allocates a vector of size int's.
   else
      nothing happens
   endif
 
   created -- 96aug28, cca
   ---------------------------------------------
*/
void
IV_init (
   IV    *iv,
   int   size,
   int   *entries 
) ;
/*
   -------------------------
   basic initializion method
 
   created -- 95oct06, cca
   -------------------------
*/
void
IV_init1 ( 
   IV    *iv,
   int   Maxsize
) ;
/*
   -------------------------
   total initializion method
 
   created -- 95oct06, cca
   -------------------------
*/
void
IV_init2 ( 
   IV    *iv,
   int   size, 
   int   maxsize, 
   int   owned, 
   int   *vec 
) ;
/*
   ----------------------------------
   set the maximum size of the vector
 
   created -- 96dec08, cca
   ----------------------------------
*/
void
IV_setMaxsize (
   IV    *iv,
   int   newmaxsize
) ;
/*
   --------------------------
   set the size of the vector
 
   created -- 96dec08, cca
   --------------------------
*/
void
IV_setSize (
   IV    *iv,
   int   newsize
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------
   return value is 0 if the entries are not owned by the object
   otherwise the return value is the number of entries at the base
   storage of the vector
 
   created  -- 96jun22, cca
   modified -- 96aug28, cca
   ---------------------------------------------------------------
*/
int
IV_owned (
   IV   *iv
) ;
/*
   ----------------------------------
   return the vector size and entries

   created -- 95oct06, cca
   ----------------------------------
*/
int
IV_size (
   IV   *iv
) ;
/*
   -------------------------------------
   return the maximum size of the vector
 
   created -- 95dec08, cca
   -------------------------------------
*/
int
IV_maxsize (
   IV   *iv
) ;
/*
   ------------------------------------------------
   return the loc'th entry of a vector.
   note: if loc is out of range then -1 is returned

   created -- 96jun29, cca
   ------------------------------------------------
*/
int 
IV_entry (
   IV    *iv,
   int   loc
) ;
/*
   ----------------------------------------------
   return a pointer to the object's entries array
 
   created -- 95oct06, cca
   ----------------------------------------------
*/
int *
IV_entries (
   IV   *iv
) ;
/*
   --------------------------------------------
   fill *psize with the vector's size
   and *pentries with the address of the vector
 
   created -- 95oct06, cca
   --------------------------------------------
*/
void
IV_sizeAndEntries (
   IV    *iv,
   int   *psize,
   int   **pentries
) ;
/*
   --------------------------
   set the size of the vector
 
   created -- 96jun22, cca
   --------------------------
*/
void
IV_setSize ( 
   IV    *iv,
   int   Size
) ;
/*
   ---------------------------
   set and entry in the vector
 
   created -- 96jul14, cca
   ---------------------------
*/
void
IV_setEntry ( 
   IV    *iv,
   int   location,
   int   value
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   shift the base of the entries and adjust the dimensions
   note: this is a dangerous operation because the iv->vec
   does not point to the base of the entries any longer,
   and thus if the object owns its entries and it is called
   to resize them or to free them, malloc and free will choke.
 
   USE WITH CAUTION!
 
   created -- 96aug25, cca
   modified -- 96aug28, cca
      structure changed
   -----------------------------------------------------------
*/
void
IV_shiftBase (
   IV    *iv,
   int    offset
) ;
/*
   ---------------------------
   push an entry onto the list

   created -- 95oct06, cca
   modified -- 95aug28, cca
      structure has changed
   ---------------------------
*/
void
IV_push (
   IV    *iv,
   int   val
) ;
/*
   ---------------------------
   minimum and maximum entries

   created -- 95oct06, cca
   ---------------------------
*/
int 
IV_min ( 
   IV   *iv
) ;
int 
IV_max ( 
   IV   *iv
) ;
int 
IV_sum ( 
   IV   *iv
) ;
/*
   -------------------------------------------------------
   sort each index list into ascending or descending order

   created -- 95oct06, cca
   -------------------------------------------------------
*/
void
IV_sortUp ( 
   IV   *iv
) ;
void
IV_sortDown ( 
   IV   *iv
) ;
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
) ;
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
) ;
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
IV_sizeOf ( 
   IV   *iv
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
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
) ;
/*
   ---------------------------------------------------
   purpose -- return invlistIV, an IV object
              that contains the inverse map,
              i.e., invlist[list[ii]] = ii.
              other entries of invlist[] are -1.
              all entris in listIV must be nonnegative
 
   created -- 98aug12, cca
   ---------------------------------------------------
*/
IV *
IV_inverseMap (
   IV   *listIV
) ;
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   purpose -- to read in an IV object from a file

   input --

      fn -- filename, must be *.ivb or *.ivf

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
IV_readFromFile ( 
   IV    *iv, 
   char   *fn 
) ;
/*
   -----------------------------------------------------
   purpose -- to read an IV object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   -----------------------------------------------------
*/
int
IV_readFromFormattedFile ( 
   IV    *iv, 
   FILE   *fp 
) ;
/*
   ---------------------------------------------------
   purpose -- to read an IV object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95oct06, cca
   ---------------------------------------------------
*/
int
IV_readFromBinaryFile ( 
   IV    *iv, 
   FILE   *fp 
) ;
/*
   -------------------------------------------
   purpose -- to write an IV object to a file

   input --

      fn -- filename
        *.ivb -- binary
        *.ivf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -------------------------------------------
*/
int
IV_writeToFile ( 
   IV    *iv, 
   char   *fn 
) ;
/*
   -----------------------------------------------------
   purpose -- to write an IV object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -----------------------------------------------------
*/
int
IV_writeToFormattedFile ( 
   IV    *iv, 
   FILE   *fp 
) ;
/*
   --------------------------------------------------
   purpose -- to write an IV object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   --------------------------------------------------
*/
int
IV_writeToBinaryFile ( 
   IV    *iv, 
   FILE   *fp 
) ;
/*
   -------------------------------------------------
   purpose -- to write an IV object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -------------------------------------------------
*/
int
IV_writeForHumanEye ( 
   IV    *iv, 
   FILE   *fp 
) ;
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the IV object

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   ---------------------------------------------------------
*/
int
IV_writeStats ( 
   IV    *iv, 
   FILE   *fp 
) ;
/*
   -------------------------------------------------------------------
   purpose -- to write out an integer vector with eighty column lines

   input --
 
      fp     -- file pointer, must be formatted and write access
      column -- present column
      pierr  -- pointer to int to hold return value, 
                should be 1 if any print was successful,
                if fprintf() failed, then ierr = -1
  
   return value -- present column
 
   created -- 96jun22, cca
   -------------------------------------------------------------------
*/
int
IV_fp80 (
   IV     *iv,
   FILE   *fp,
   int    column,
   int    *pierr
) ;
/*
   ---------------------------------------------------
   purpose -- to write the IV object for a matlab file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98oct22, cca
   ---------------------------------------------------
*/
int
IV_writeForMatlab ( 
   IV     *iv, 
   char   *name,
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
