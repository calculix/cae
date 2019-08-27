/*  init.c  */

#include "../SubMtxList.h"

/*--------------------------------------------------------------------*/
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
                       only threads in this and other processes.
   flags -- vector to specify whether to lock individual lists
      flags == NULL --> none or all lists must be locked,
                        use lockflag to determine
      flags[ilist] = 'N' --> no need to lock list ilist
      flags[ilist] = 'Y' --> must lock list ilist

   created -- 98may02, cca
   ------------------------------------------------------------------
*/
void
SubMtxList_init (
   SubMtxList   *list,
   int          nlist,
   int          counts[],
   int          lockflag,
   char         flags[]
) {
int   ilist ;
/*
   ---------------
   check the input
   ---------------
*/
if ( list == NULL || nlist <= 0 || lockflag < 0 || lockflag > 2 ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxList_init(%p,%d,%p,%d,%p)"
           "\n bad input\n", list, nlist, counts, lockflag, flags) ;
   exit(-1) ;
}
/*
   --------------
   clear all data
   --------------
*/
SubMtxList_clearData(list) ;
/*
   -------------------------------------------------------
   set the number of lists and allocate the heads[] vector
   -------------------------------------------------------
*/
list->nlist = nlist ;
ALLOCATE(list->heads, struct _SubMtx *, nlist) ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   list->heads[ilist] = NULL ;
}
if ( counts != NULL ) {
/*
   -------------------------------------
   allocate and fill the counts[] vector
   -------------------------------------
*/
   list->counts = IVinit(nlist, 0) ;
   IVcopy(nlist, list->counts, counts) ;
}
if ( lockflag > 0 ) {
/*
   -----------------
   allocate the lock
   -----------------
*/
   list->lock = Lock_new() ;
   Lock_init(list->lock, lockflag) ;
}
if ( flags != NULL ) {
/*
   ------------------------------------
   allocate and fill the flags[] vector
   ------------------------------------
*/
   list->flags = CVinit(nlist, 'N') ;
   CVcopy(nlist, list->flags, flags) ;
}
return ; }

/*--------------------------------------------------------------------*/
