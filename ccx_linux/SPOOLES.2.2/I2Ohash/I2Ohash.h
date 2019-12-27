/*  I2Ohash.h  */

#include "../Utilities.h"

/*--------------------------------------------------------------------*/
typedef struct _I2Ohash   I2Ohash ;
struct _I2Ohash {
   int    nlist     ;
   int    grow      ;
   int    nitem     ;
   I2OP   *baseI2OP ;
   I2OP   *freeI2OP ;
   I2OP   **heads   ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c  ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   create and return a new instance of the I2Ohash object
 
   created -- 98jan29, cca
   -----------------------------------------------------
*/
I2Ohash *
I2Ohash_new (
   void
) ;
/*
   -------------------------------------------
   set the default fields for an I2Ohash object
 
   created -- 98jan29, cca
   -------------------------------------------
*/
void
I2Ohash_setDefaultFields (
   I2Ohash   *hashtbl
) ;
/*
   -----------------------
   clear the data fields
 
   created -- 98jan29, cca
   -----------------------
*/
void
I2Ohash_clearData (
   I2Ohash   *hashtbl
) ;
/*
   -----------------------
   free the I2Ohash object
 
   created -- 98jan29, cca
   -----------------------
*/
void
I2Ohash_free (
   I2Ohash   *hashtbl
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c  -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------
   initializer, 
   (1) set the number of lists to nlist
       and allocate the heads[] vector
   (2) initialize nobj I2OP objects
   (2) set hashtbl item growth factor to grow
 
   created -- 98jan29, cca
   ------------------------------------------
*/
void
I2Ohash_init ( 
   I2Ohash   *hashtbl,
   int      nlist,
   int      nobj,
   int      grow 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c  -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------
   insert the (key1, key2, value) triple into the hash table
 
   created -- 98jan29, cca
   ---------------------------------------------------------
*/
void
I2Ohash_insert (
   I2Ohash   *hashtbl,
   int       key1,
   int       key2,
   void      *value
) ;
/*
   ------------------------------------------------
   locate the first (key1,key2,*) in the hash table
 
   return value --
      0 -- no (key1,key2,*) value triple
      1 -- (key1,key2,value) value triple found,
           value put into *pvalue
 
   created -- 98jan29, cca
   ------------------------------------------------
*/
int
I2Ohash_locate (
   I2Ohash   *hashtbl,
   int       key1,
   int       key2,
   void      **pvalue
) ;
/*
   --------------------------------------------------
   remove the first (key1,key2,*) from the hash table
 
   return value --
      0 -- no (key1,key2,*) value triple
      1 -- (key1,key2,value) value triple found,
           value put into *pvalue
 
   created -- 98jan29, cca
   --------------------------------------------------
*/
int
I2Ohash_remove (
   I2Ohash   *hashtbl,
   int       key1,
   int       key2,
   void      **pvalue
) ;
/*
   ---------------------------------------------------------------
   return a measure of the nonuniform distribution of the entries.
   a value of 1.0 is best.
 
   created -- 98jan29, cca
   ---------------------------------------------------------------
*/
double
I2Ohash_measure (
   I2Ohash   *hashtbl
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c  -------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- to write the I2Ohash object to a file
 
   created -- 98jan29, cca
   -----------------------------------------------
*/
void
I2Ohash_writeForHumanEye ( 
   I2Ohash   *hashtbl,
   FILE     *fp 
) ;
/*--------------------------------------------------------------------*/
