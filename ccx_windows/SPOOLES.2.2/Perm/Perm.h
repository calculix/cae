/*  Perm.h  */

#include "../cfiles.h"
#include "../IV.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   the Perm object contains one or two permutation vectors

   isPresent -- flag to tell which vectors are present
     0 --> neither are present
     1 --> newToOld is present
     2 --> oldToNew is present
     3 --> both are present
   size -- size of the vectors
   newToOld -- pointer to new-to-old permutation vector
   oldToNew -- pointer to old-to-new permutation vector

   created -- 96mar16
   ---------------------------------------------------------
*/
typedef struct _Perm   Perm ;
struct _Perm {
   int   isPresent ;
   int   size      ;
   int   *newToOld ;
   int   *oldToNew ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method founds in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor

   created -- 96jan05, cca
   -----------------------
*/
Perm *
Perm_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields

   created -- 96jan05, cca
   -----------------------
*/
void
Perm_setDefaultFields (
   Perm   *perm
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 96jan05, cca
   --------------------------------------------------
*/
void
Perm_clearData ( 
   Perm   *perm 
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 96jan05, cca
   ------------------------------------------
*/
Perm *
Perm_free ( 
   Perm   *perm 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method founds in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   initializer

   created -- 96jan05, cca
   -----------------------
*/
void
Perm_initWithTypeAndSize ( 
   Perm   *perm,
   int    isPresent, 
   int    size 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method founds in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------
   return the storage taken by this object

   created -- 96jan05, cca
   ---------------------------------------
*/
int
Perm_sizeOf ( 
   Perm   *perm
) ;
/*
   ----------------------------------------------------------
   check that the permutation object does house a permutation

   return value --
      1 if a true permutation
      0 otherwise
   ----------------------------------------------------------
*/
int
Perm_checkPerm (
   Perm   *perm
) ;
/*
   ----------------------------------------
   if the old-to-new vector is not present,
   create it and fill its entries
 
   created -- 96mar16, cca
   ----------------------------------------
*/
void
Perm_fillOldToNew (
   Perm   *perm
) ;
/*
   ----------------------------------------
   if the new-to-old vector is not present,
   create it and fill its entries
 
   created -- 96mar16, cca
   ----------------------------------------
*/
void
Perm_fillNewToOld (
   Perm   *perm
) ;
/*
   ------------------------------------
   if the old-to-new vector is present,
   release it and free its entries
 
   created -- 96mar16, cca
   ------------------------------------
*/
void
Perm_releaseOldToNew (
   Perm   *perm
) ;
/*
   ------------------------------------
   if the new-to-old vector is present,
   release it and free its entries
 
   created -- 96mar16, cca
   ------------------------------------
*/
void
Perm_releaseNewToOld (
   Perm   *perm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method founds in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   purpose -- to read in an Perm object from a file

   input --

      fn -- filename, must be *.permb or *.permf

   return value -- 1 if success, 0 if failure

   created -- 96jan05, cca
   -----------------------------------------------
*/
int
Perm_readFromFile ( 
   Perm    *perm, 
   char   *fn 
) ;
/*
   ------------------------------------------------------
   purpose -- to read an Perm object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 96jan05, cca
   ------------------------------------------------------
*/
int
Perm_readFromFormattedFile ( 
   Perm    *perm, 
   FILE   *fp 
) ;
/*
   ---------------------------------------------------
   purpose -- to read an Perm object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 96jan05, cca
   ---------------------------------------------------
*/
int
Perm_readFromBinaryFile ( 
   Perm    *perm, 
   FILE   *fp 
) ;
/*
   -------------------------------------------
   purpose -- to write an Perm object to a file

   input --

      fn -- filename
        *.permb -- binary
        *.permf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   -------------------------------------------
*/
int
Perm_writeToFile ( 
   Perm    *perm, 
   char   *fn 
) ;
/*
   -----------------------------------------------------
   purpose -- to write an Perm object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   -----------------------------------------------------
*/
int
Perm_writeToFormattedFile ( 
   Perm    *perm, 
   FILE   *fp 
) ;
/*
   --------------------------------------------------
   purpose -- to write an Perm object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   --------------------------------------------------
*/
int
Perm_writeToBinaryFile ( 
   Perm    *perm, 
   FILE   *fp 
) ;
/*
   -------------------------------------------------
   purpose -- to write an Perm object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   -------------------------------------------------
*/
int
Perm_writeForHumanEye ( 
   Perm    *perm, 
   FILE   *fp 
) ;
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the Perm object

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   ---------------------------------------------------------
*/
int
Perm_writeStats ( 
   Perm    *perm, 
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method founds in compress.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   given a permutation and a vector to map vertices 
   into compressed vertices, create and return a 
   permutation object for the compressed vertices.
 
   created -- 96may02, cca
   ------------------------------------------------
*/
Perm *
Perm_compress (
   Perm   *perm,
   IV     *eqmapIV
) ;
/*--------------------------------------------------------------------*/
