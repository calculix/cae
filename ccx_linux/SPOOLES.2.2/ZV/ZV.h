/*  ZV.h  */

#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   ZV -- double complex vector object
 
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
typedef struct _ZV   ZV ;
struct _ZV {
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
 
   created -- 98jan22, cca
   -----------------------
*/
ZV *
ZV_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98jan22, cca
   -----------------------
*/
void
ZV_setDefaultFields ( 
   ZV   *zv
) ;
/*
   -----------------------
   clear the data fields
 
   created -- 98jan22, cca
   -----------------------
*/
void
ZV_clearData (
   ZV   *zv
) ;
/*
   -----------------------
   destructor
 
   created -- 98jan22, cca
   -----------------------
*/
void
ZV_free (
   ZV   *zv
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------
   simplest initialization method
 
   data is cleared
   if entries != NULL
      the object does not own the entries,
      it just points to the entries base address
   else if size > 0
      the object will own the entries, 
      it allocates a vector of 2*size doubles's.
   else 
      nothing happens
   endif
 
   created -- 98jan22, cca
   ---------------------------------------------
*/
void
ZV_init (
   ZV       *zv,
   int      size,
   double   *entries 
) ;
/*
   -------------------------
   basic initializion method
 
   created -- 98jan22, cca
   -------------------------
*/
void
ZV_init1 (
   ZV    *zv,
   int   size
) ;
/*
   -------------------------
   total initializion method
 
   created -- 98jan22, cca
   -------------------------
*/
void
ZV_init2 (
   ZV       *zv,
   int      size,
   int      maxsize,
   int      owned,
   double   *vec
) ;
/*
   ----------------------------------
   set the maximum size of the vector
 
   created -- 98jan22, cca
   ----------------------------------
*/
void
ZV_setMaxsize (
   ZV    *zv,
   int   newmaxsize
) ;
/*
   --------------------------
   set the size of the vector
 
   created -- 98jan22, cca
   --------------------------
*/
void
ZV_setSize (
   ZV    *zv,
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
 
   created -- 98jan22, cca
   -----------------------------------------------
*/
int
ZV_owned (
   ZV   *dv
) ;
/*
   -----------------------
   return the vector size
 
   created -- 98jan22, cca
   -----------------------
*/
int
ZV_maxsize (
   ZV   *dv
) ;
/*
   -----------------------
   return the vector size
 
   created -- 98jan22, cca
   -----------------------
*/
int
ZV_size (
   ZV   *dv
) ;
/*
   -------------------------------------------------
   return the loc'th entry of a vector.
 
   created -- 98jan22, cca
   -------------------------------------------------
*/
void
ZV_entry (
   ZV       *dv,
   int      loc,
   double   *pReal,
   double   *pImag
) ;
/*
   -------------------------------------------------
   return pointers to the loc'th entry of a vector.
 
   created -- 98jan22, cca
   -------------------------------------------------
*/
void
ZV_pointersToEntry (
   ZV       *dv,
   int      loc,
   double   **ppReal,
   double   **ppImag
) ;
/*
   ----------------------------------------------
   return a pointer to the object's entries array
 
   created -- 98jan22, cca
   ----------------------------------------------
*/
double *
ZV_entries (
   ZV   *dv
) ;
/*
   --------------------------------------------
   fill *psize with the vector's size
   and *pentries with the address of the vector
 
   created -- 98jan22, cca
   --------------------------------------------
*/
void
ZV_sizeAndEntries (
   ZV       *dv,
   int      *psize,
   double   **pentries
) ;
/*
   ---------------------------
   set and entry in the vector
 
   created -- 98jan22, cca
   ---------------------------
*/
void
ZV_setEntry (
   ZV       *dv,
   int      loc,
   double   real,
   double   imag
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   shift the base of the entries and adjust the dimensions
   note: this is a dangerous operation because the zv->vec
   does not point to the base of the entries any longer,
   and thus if the object owns its entries and it is called
   to resize them or to free them, malloc and free will choke.
 
   USE WITH CAUTION!
 
   created  -- 98jan22, cca
   -----------------------------------------------------------
*/
void
ZV_shiftBase (
   ZV    *zv,
   int    offset
) ;
/*
   --------------------------------------
   push an entry onto the list
 
   created -- 95oct06, cca
   --------------------------------------
*/
void
ZV_push (
   ZV       *zv,
   double   real,
   double   imag
) ;
/*
   ---------------------------
   minimum and maximum entries
 
   created -- 95oct06, cca
   ---------------------------
*/
double
ZV_minabs (
   ZV   *zv
) ;
/*
   ----------------------------------------------
   return the number of bytes taken by the object
 
   created -- 95oct06, cca
   ----------------------------------------------
*/
int
ZV_sizeOf (
   ZV   *zv
) ;
/*
   --------------------------
   fill a vector with a value
 
   created -- 96jun22, cca
   --------------------------
*/
void
ZV_fill (
   ZV       *zv,
   double   real,
   double   imag
) ;
/*
   ------------------------
   fill a vector with zeros
 
   created -- 98jun02, cca
   ------------------------
*/
void
ZV_zero (
   ZV   *zv
) ;
/*
   --------------------------------------
   copy entries from zv2 into zv1.
   note: this is a "mapped" copy,
   zv1 and zv2 need not be the same size.
 
   created -- 96aug31, cca
   --------------------------------------
*/
void
ZV_copy (
   ZV   *zv1,
   ZV   *zv2
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
   the entries in the ZV object. tausmall and tau big provide
   cutoffs within which to examine the entries. pnzero, pnsmall
   and pnbig are addresses to hold the number of entries zero,
   smaller than tausmall and larger than taubig, respectively.
 
   created -- 97feb14, cca
   ------------------------------------------------------------------
*/
void
ZV_log10profile (
   ZV      *zv,
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
   purpose -- to read in an ZV object from a file
 
   input --
 
      fn -- filename, must be *.zvb or *.zvf
 
   return value -- 1 if success, 0 if failure
 
   created -- 98jan22, cca
   ----------------------------------------------
*/
int
ZV_readFromFile (
   ZV     *zv,
   char   *fn
) ;
/*
   -----------------------------------------------------
   purpose -- to read an ZV object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 98jan22, cca
   -----------------------------------------------------
*/
int
ZV_readFromFormattedFile (
   ZV     *zv,
   FILE   *fp
) ;
/*
   ---------------------------------------------------
   purpose -- to read an ZV object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 98jan22, cca
   ---------------------------------------------------
*/
int
ZV_readFromBinaryFile (
   ZV    *zv,
   FILE   *fp
) ;
/*
   -------------------------------------------
   purpose -- to write an ZV object to a file
 
   input --
 
      fn -- filename
        *.zvb -- binary
        *.zvf -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan22, cca
   -------------------------------------------
*/
int
ZV_writeToFile (
   ZV    *zv,
   char   *fn
) ;
/*
   -----------------------------------------------------
   purpose -- to write an ZV object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan22, cca
   -----------------------------------------------------
*/
int
ZV_writeToFormattedFile (
   ZV    *zv,
   FILE   *fp
) ;
/*
   --------------------------------------------------
   purpose -- to write an ZV object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan22, cca
   --------------------------------------------------
*/
int
ZV_writeToBinaryFile (
   ZV    *zv,
   FILE   *fp
) ;
/*
   -------------------------------------------------
   purpose -- to write an ZV object for a human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan22, cca
   -------------------------------------------------
*/
int
ZV_writeForHumanEye (
   ZV    *zv,
   FILE   *fp
) ;
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the ZV object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan22, cca
   ---------------------------------------------------------
*/
int
ZV_writeStats (
   ZV    *zv,
   FILE   *fp
) ;
/*
   --------------------------------------------------
   purpose -- write the vector entries out for matlab
  
   created -- 98apr15, cca
   --------------------------------------------------
*/
void
ZV_writeForMatlab (
   ZV     *zv,
   char   *vecname,
   FILE   *fp
) ;
/*--------------------------------------------------------------------*/
