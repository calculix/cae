/*  ChvList.h  */

#include "../Chv.h"
#include "../Lock.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   this object handles a list of lists of Chv objects

   nlist  -- # of lists
   heads  -- heads[ilist] contains a pointer 
             to the first Chv object in list ilist
   counts -- when not-NULL, counts[ilist] contains the remaining number 
             of objects to be added to list ilist before it is complete
   lock   -- mutex object, can be NULL
   flags  -- when not NULL, a vector to specify when a list needs 
             to be locked before adding an object to it.
      flags[ilist] = 'N' --> no need to lock
      flags[ilist] = 'Y' --> must lock
   nlocks -- number of times the list was locked
   --------------------------------------------------------------------
*/
typedef struct _ChvList   ChvList ;
struct _ChvList {
   int       nlist   ;
   Chv       **heads ;
   int       *counts ;
   Lock      *lock   ;
   char      *flags  ;
   int       nlocks  ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98may02, cca
   -----------------------
*/
ChvList *
ChvList_new (
   void
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may02, cca
   -----------------------
*/
void
ChvList_setDefaultFields (
   ChvList   *chvlist
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may02, cca
   --------------------------------------------------
*/
void
ChvList_clearData (
   ChvList   *chvlist
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may02, cca
   ------------------------------------------
*/
void
ChvList_free (
   ChvList   *chvlist
) ;

/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- basic initializer
 
   nlist  -- number of lists to be held by this object
   counts -- vector that contains number of items expected
             for each list.
      counts == NULL --> unknown number of items expected
      counts != NULL --> known number of items expected
   lockflag -- flag to specify lock status
      lockflag = 0 --> mutex lock is not allocated or initialized
      lockflag = 1 --> mutex lock is allocated and it can synchronize
                       only threads in this process.
      lockflag = 2 --> mutex lock is allocated and it can synchronize
                       threads in this and other processes.
   flags -- vector to specify whether to lock individual lists
      flags == NULL --> none or all lists must be locked,
                        use lockflag to determine
      flags[ilist] = 'N' --> no need to lock list ilist
      flags[ilist] = 'Y' --> must lock list ilist
 
   created -- 98may02, cca
   ------------------------------------------------------------------
*/
void
ChvList_init (
   ChvList   *chvlist,
   int       nlist,
   int       counts[],
   int       lockflag,
   char      flags[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
ChvList_writeForHumanEye (
   ChvList   *chvlist,
   FILE       *fp
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------
   return 1 if list ilist is not empty
   return 0 if list ilist is empty
 
   created -- 98may02, cca
   -----------------------------------
*/
int
ChvList_isListNonempty (
   ChvList   *chvlist,
   int       ilist
) ;
/*
   ---------------------------------------------------------
   return 1 if the count for list ilist is zero
   return 0 if the count for list ilist is greater than zero
 
   created -- 98may02, cca
   ---------------------------------------------------------
*/
int
ChvList_isCountZero (
   ChvList   *chvlist,
   int       ilist
) ;
/*
   ---------------------------------
   if chv is not NULL then
      add chv to list ilist
   endif
   decrement the count of list ilist
 
   created -- 98may02, cca
   ---------------------------------
*/
void
ChvList_addObjectToList (
   ChvList   *chvlist,
   Chv       *chv,
   int       ilist
) ;
/*
   ------------------------------------
   return pointer to head of list ilist
   and set head to NULL
 
   created -- 98may02, cca
   ------------------------------------
*/
Chv *
ChvList_getList (
   ChvList   *chvlist,
   int       ilist
) ;
/*--------------------------------------------------------------------*/
