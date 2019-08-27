/*  DV.h  */

#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   DV -- double vector object
 
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
typedef struct _DV   DV ;
struct _DV {
   int      size    ;
   int      maxsize ;
   int      owned   ;
   double   *vec    ;
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

   created -- 96jun23, cca
   -----------------------
*/
DV *
DV_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields

   created -- 96jun23, cca
   -----------------------
*/
void
DV_setDefaultFields ( 
   DV   *dv
) ;
/*
   -----------------------
   clear the data fields

   created -- 96jun23, cca
   -----------------------
*/
void
DV_clearData ( 
   DV   *dv
) ;
/*
   -----------------------
   destructor

   created -- 96jun23, cca
   -----------------------
*/
void
DV_free ( 
   DV   *dv
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
DV_init (
   DV       *iv,
   int      size,
   double   *entries 
) ;
/*
   -------------------------
   basic initializion method
 
   created -- 96jun23, cca
   -------------------------
*/
void
DV_init1 ( 
   DV    *dv,
   int   Maxsize
) ;
/*
   -------------------------
   total initializion method
 
   created -- 96jun23, cca
   -------------------------
*/
void
DV_init2 ( 
   DV       *dv,
   int      Size, 
   int      Maxsize, 
   int      owned, 
   double   *Dvec 
) ;
/*
   ----------------------------------
   set the maximum size of the vector
 
   created -- 96dec08, cca
   ----------------------------------
*/
void
DV_setMaxsize (
   DV    *dv,
   int   newmaxsize
) ;
/*
   --------------------------
   set the size of the vector
 
   created -- 96dec08, cca
   --------------------------
*/
void
DV_setSize (
   DV    *dv,
   int   newsize
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   return 1 if the entries are owned by the object
   return 0 otherwise

   created -- 96jun22, cca
   -----------------------------------------------
*/
int
DV_owned (
   DV   *dv
) ;
/*
   -----------------------
   return the vector size

   created -- 95oct06, cca
   -----------------------
*/
int
DV_maxsize (
   DV   *dv
) ;
/*
   -----------------------
   return the vector size

   created -- 95oct06, cca
   -----------------------
*/
int
DV_size (
   DV   *dv
) ;
/*
   -------------------------------------------------
   return the loc'th entry of a vector.
   note: if loc is out of range then 0.0 is returned
 
   created -- 96jun29, cca
   -------------------------------------------------
*/
double 
DV_entry (
   DV    *dv,
   int   loc
) ;
/*
   ----------------------------------------------
   return a pointer to the object's entries array

   created -- 95oct06, cca
   ----------------------------------------------
*/
double *
DV_entries (
   DV   *dv
) ;
/*
   --------------------------------------------
   fill *psize with the vector's size
   and *pentries with the address of the vector

   created -- 95oct06, cca
   --------------------------------------------
*/
void
DV_sizeAndEntries (
   DV       *dv,
   int      *psize,
   double   **pentries
) ;
/*
   ---------------------------
   set and entry in the vector
 
   created -- 96jul14, cca
   ---------------------------
*/
void
DV_setEntry (
   DV       *dv,
   int      loc,
   double   value
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
   note: this is a dangerous operation because the dv->vec
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
DV_shiftBase (
   DV    *dv,
   int    offset
) ;
/*
   -----------------------
   resize the vector
 
   created -- 96jun22, cca
   -----------------------
*/
void
DV_resize (
   DV    *dv,
   int   Maxsize
) ;
/*
   -----------------------
   set the size to be zero

   created -- 96jun23, cca
   -----------------------
*/
void
DV_clear (
   DV   *dv
) ;
/*
   ---------------------------
   push an entry onto the list

   created -- 96jun23, cca
   ---------------------------
*/
void
DV_push (
   DV       *dv,
   double   val
) ;
/*
   ------------------------------------
   minimum and maximum entries and sum

   created -- 96jun23, cca
   ------------------------------------
*/
double 
DV_min ( 
   DV   *dv
) ;
double 
DV_max ( 
   DV   *dv
) ;
double 
DV_sum ( 
   DV   *dv
) ;
/*
   -------------------------------------------------------
   sort each index list into ascending or descending order

   created -- 96jun23, cca
   -------------------------------------------------------
*/
void
DV_sortUp ( 
   DV   *dv
) ;
void
DV_sortDown ( 
   DV   *dv
) ;
/*
   -----------------------
   ramp the entries

   created -- 96jun23, cca
   -----------------------
*/
void
DV_ramp (
   DV       *dv,
   double   base,
   double   incr
) ;
/*
   -----------------------
   shuffle the list

   created -- 96jun23, cca
   -----------------------
*/
void
DV_shuffle (
   DV    *dv, 
   int   seed
) ;
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 96jun23, cca
   ----------------------------------------------
*/
int
DV_sizeOf ( 
   DV   *dv
) ;
/*
   --------------------------------------------
   iterator :
   return the pointer to the head of the vector

   created -- 96jun23, cca
   --------------------------------------------
*/
double *
DV_first (
   DV   *dv
) ;
/*
   -----------------------------------------------------
   iterator :
   return the pointer to the next location in the vector

   created -- 96jun23, cca
   -----------------------------------------------------
*/
double *
DV_next (
   DV       *dv,
   double   *pi
) ;
/*
   --------------------------
   fill a vector with a value
 
   created -- 96jun22, cca
   --------------------------
*/
void
DV_fill (
   DV       *dv,
   double   value
) ;
/*
   --------------------------
   fill a vector with zeros
 
   created -- 98jun02, cca
   --------------------------
*/
void
DV_zero (
   DV   *dv
) ;
/*
   --------------------------------------
   copy entries from dv2 into dv1.
   note: this is a "mapped" copy, 
   dv1 and dv2 need not be the same size.
 
   created -- 96aug31, cca
   --------------------------------------
*/
void
DV_copy (
   DV   *dv1,
   DV   *dv2
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in profile.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   to fill xDV and yDV with a log10 profile of the magnitudes of
   the entries in the DV object. tausmall and tau big provide
   cutoffs within which to examine the entries. pnzero, pnsmall 
   and pnbig are addresses to hold the number of entries zero,
   smaller than tausmall and larger than taubig, respectively.
 
   created -- 97feb14, cca
   ------------------------------------------------------------------
*/
void
DV_log10profile (
   DV      *dv,
   int      npts,
   DV       *xDV,
   DV       *yDV,
   double   tausmall,
   double   taubig,
   int      *pnzero,
   int      *pnsmall,
   int      *pnbig
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   purpose -- to read in an DV object from a file

   input --

      fn -- filename, must be *.dvb or *.dvf

   return value -- 1 if success, 0 if failure

   created -- 96jun23, cca
   ----------------------------------------------
*/
int
DV_readFromFile ( 
   DV    *dv, 
   char   *fn 
) ;
/*
   -----------------------------------------------------
   purpose -- to read an DV object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 96jun23, cca
   -----------------------------------------------------
*/
int
DV_readFromFormattedFile ( 
   DV    *dv, 
   FILE   *fp 
) ;
/*
   ---------------------------------------------------
   purpose -- to read an DV object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 96jun23, cca
   ---------------------------------------------------
*/
int
DV_readFromBinaryFile ( 
   DV    *dv, 
   FILE   *fp 
) ;
/*
   -------------------------------------------
   purpose -- to write an DV object to a file

   input --

      fn -- filename
        *.dvb -- binary
        *.dvf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   -------------------------------------------
*/
int
DV_writeToFile ( 
   DV    *dv, 
   char   *fn 
) ;
/*
   -----------------------------------------------------
   purpose -- to write an DV object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   -----------------------------------------------------
*/
int
DV_writeToFormattedFile ( 
   DV    *dv, 
   FILE   *fp 
) ;
/*
   --------------------------------------------------
   purpose -- to write an DV object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   --------------------------------------------------
*/
int
DV_writeToBinaryFile ( 
   DV    *dv, 
   FILE   *fp 
) ;
/*
   -------------------------------------------------
   purpose -- to write an DV object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   -------------------------------------------------
*/
int
DV_writeForHumanEye ( 
   DV    *dv, 
   FILE   *fp 
) ;
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the DV object

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   ---------------------------------------------------------
*/
int
DV_writeStats ( 
   DV    *dv, 
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
DV_fp80 (
   DV     *dv,
   FILE   *fp,
   int    column,
   int    *pierr
) ;
/*
   ---------------------------------------------------
   purpose -- to write the DV object for a matlab file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98feb07, cca
   ---------------------------------------------------
*/
int
DV_writeForMatlab ( 
   DV     *dv, 
   char   *name,
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
