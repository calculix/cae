/*  util.c  */

#include "../ChvList.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
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
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chvlist == NULL || ilist < 0 || ilist >= chvlist->nlist ) {
   fprintf(stderr, "\n fatal error in ChvList_isListNonempty(%p,%d)"
           "\n bad input\n", chvlist, ilist) ;
}
return(chvlist->heads[ilist] != NULL) ; }

/*--------------------------------------------------------------------*/
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
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chvlist == NULL || ilist < 0 || ilist >= chvlist->nlist ) {
   fprintf(stderr, "\n fatal error in ChvList_isCountZero(%p,%d)"
           "\n bad input\n", chvlist, ilist) ;
}
if ( chvlist->counts == NULL ) {
   return(1) ;
} else {
   return(chvlist->counts[ilist] == 0) ;
}
}

/*--------------------------------------------------------------------*/
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
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chvlist == NULL || ilist < 0 || ilist >= chvlist->nlist ) {
   fprintf(stderr, 
           "\n fatal error in ChvList_addObjectToList(%p,%p,%d)"
           "\n bad input\n", chvlist, chv, ilist) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n ChvList %p : adding chv %p to list %d",
        chvlist, chv, ilist) ;
fflush(stdout) ;
#endif
if ( chvlist->lock != NULL
   && (chvlist->flags == NULL || chvlist->flags[ilist] == 'Y' ) ) {
/*
   --------------------------------------------------
   we must lock the list to (possibly) add the object
   and decrement the list's count
   --------------------------------------------------
*/
   Lock_lock(chvlist->lock) ;
   if ( chv != NULL ) {
      chv->next = chvlist->heads[ilist] ;
      chvlist->heads[ilist] = chv ;
   }
   if ( chvlist->counts != NULL ) {
      chvlist->counts[ilist]-- ;
   }
   chvlist->nlocks++ ;
   Lock_unlock(chvlist->lock) ;
} else {
/*
   ---------------------------------------------
   no need to lock the list, just (possibly)
   add the object and decrement the list's count
   ---------------------------------------------
*/
   if ( chv != NULL ) {
      chv->next = chvlist->heads[ilist] ;
      chvlist->heads[ilist] = chv ;
   }
   if ( chvlist->counts != NULL ) {
      chvlist->counts[ilist]-- ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n ChvList %p : heads[%d] = %p, counts[%d] = %d",
        chvlist, ilist, chvlist->heads[ilist],
        ilist, chvlist->counts[ilist]) ;
fflush(stdout) ;
#endif
return ; }

/*--------------------------------------------------------------------*/
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
) {
Chv   *chv ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvlist == NULL || ilist < 0 || ilist >= chvlist->nlist ) {
   fprintf(stderr, 
           "\n fatal error in ChvList_getList(%p,%d)"
           "\n bad input\n", chvlist, ilist) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n ChvList %p : get list %d", chvlist, ilist) ;
fflush(stdout) ;
#endif
if ( chvlist->heads[ilist] != NULL ) {
   if ( chvlist->lock == NULL
     || (chvlist->flags != NULL && chvlist->flags[ilist] == 'N')
     || (chvlist->counts != NULL && chvlist->counts[ilist] == 0) ) {
/*
      ------------------------
      no need to lock the list
      ------------------------
*/
      chv = chvlist->heads[ilist] ;
      chvlist->heads[ilist] = NULL ;
   } else {
/*
      ----------------------------------------------------
      we must lock the list to return the head of the list
      ----------------------------------------------------
*/
      Lock_lock(chvlist->lock) ;
      chv = chvlist->heads[ilist] ;
      chvlist->heads[ilist] = NULL ;
      chvlist->nlocks++ ;
      Lock_unlock(chvlist->lock) ;
   }
} else {
   chv = NULL ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n ChvList %p : chv %p", chvlist, chv) ;
fflush(stdout) ;
#endif
return(chv) ; }

/*--------------------------------------------------------------------*/
