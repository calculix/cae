/*  init.c  */

#include "../ChvList.h"

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
ChvList_init (
   ChvList   *chvlist,
   int       nlist,
   int       counts[],
   int       lockflag,
   char      flags[]
) {
int   ilist ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvlist == NULL || nlist <= 0 || lockflag < 0 || lockflag > 1 ) {
   fprintf(stderr, "\n fatal error in ChvList_init(%p,%d,%p,%d,%p)"
           "\n bad input\n", chvlist, nlist, counts, lockflag, flags) ;
   exit(-1) ;
}
/*
   --------------
   clear all data
   --------------
*/
ChvList_clearData(chvlist) ;
/*
   -------------------------------------------------------
   set the number of lists and allocate the heads[] vector
   -------------------------------------------------------
*/
chvlist->nlist = nlist ;
ALLOCATE(chvlist->heads, struct _Chv *, nlist) ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   chvlist->heads[ilist] = NULL ;
}
if ( counts != NULL ) {
/*
   -------------------------------------
   allocate and fill the counts[] vector
   -------------------------------------
*/
   chvlist->counts = IVinit(nlist, 0) ;
   IVcopy(nlist, chvlist->counts, counts) ;
}
if ( lockflag > 0 ) {
/*
   -----------------
   allocate the lock
   -----------------
*/
   chvlist->lock = Lock_new() ;
   Lock_init(chvlist->lock, lockflag) ;
}
if ( flags != NULL ) {
/*
   ------------------------------------
   allocate and fill the flags[] vector
   ------------------------------------
*/
   chvlist->flags = CVinit(nlist, 'N') ;
   CVcopy(nlist, chvlist->flags, flags) ;
}
return ; }

/*--------------------------------------------------------------------*/
