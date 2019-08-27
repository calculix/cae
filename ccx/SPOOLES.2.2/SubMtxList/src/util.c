/*  util.c  */

#include "../SubMtxList.h"

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
SubMtxList_isListNonempty (
   SubMtxList   *list,
   int          ilist
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( list == NULL || ilist < 0 || ilist >= list->nlist ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxList_isListNonempty(%p,%d)"
           "\n bad input\n", list, ilist) ;
}
return(list->heads[ilist] != NULL) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   return 1 if the count for list ilist is zero
   return 0 if the count for list ilist is greater than zero

   created -- 98may02, cca
   ---------------------------------------------------------
*/
int
SubMtxList_isCountZero (
   SubMtxList   *list,
   int          ilist
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( list == NULL || ilist < 0 || ilist >= list->nlist ) {
   fprintf(stderr, "\n fatal error in SubMtxList_isCountZero(%p,%d)"
           "\n bad input\n", list, ilist) ;
}
if ( list->counts == NULL ) {
   return(1) ;
} else {
   return(list->counts[ilist] == 0) ;
}
}

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   add an object to the list and decrement the count.
   note, an object can be NULL.

   created -- 98may02, cca
   --------------------------------------------------
*/
void
SubMtxList_addObjectToList (
   SubMtxList   *list,
   SubMtx       *mtx,
   int          ilist
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( list == NULL || ilist < 0 || ilist >= list->nlist ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxList_addObjectToList(%p,%p,%d)"
           "\n bad input\n", list, mtx, ilist) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n SubMtxList %p : adding mtx %p to list %d",
        list, mtx, ilist) ;
fflush(stdout) ;
#endif
if ( list->lock != NULL
   && (list->flags == NULL || list->flags[ilist] == 'Y' ) ) {
/*
   --------------------------------------------------
   we must lock the list to (possibly) add the object
   and decrement the list's count
   --------------------------------------------------
*/
   Lock_lock(list->lock) ;
   if ( mtx != NULL ) {
      mtx->next = list->heads[ilist] ;
      list->heads[ilist] = mtx ;
   }
   if ( list->counts != NULL ) {
      list->counts[ilist]-- ;
   }
   list->nlocks++ ;
   Lock_unlock(list->lock) ;
} else {
/*
   ---------------------------------------------
   no need to lock the list, just (possibly)
   add the object and decrement the list's count
   ---------------------------------------------
*/
   if ( mtx != NULL ) {
      mtx->next = list->heads[ilist] ;
      list->heads[ilist] = mtx ;
   }
   if ( list->counts != NULL ) {
      list->counts[ilist]-- ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n SubMtxList %p : heads[%d] = %p, counts[%d] = %d",
        list, ilist, list->heads[ilist],
        ilist, list->counts[ilist]) ;
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
SubMtx *
SubMtxList_getList (
   SubMtxList   *list,
   int          ilist
) {
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if ( list == NULL || ilist < 0 || ilist >= list->nlist ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxList_getList(%p,%d)"
           "\n bad input\n", list, ilist) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n SubMtxList %p : get list %d", list, ilist) ;
fflush(stdout) ;
#endif
if ( list->heads[ilist] != NULL ) {
   if ( list->lock == NULL
     || (list->flags != NULL && list->flags[ilist] == 'N')
     || (list->counts != NULL && list->counts[ilist] == 0) ) {
/*
      ------------------------
      no need to lock the list
      ------------------------
*/
      mtx = list->heads[ilist] ;
      list->heads[ilist] = NULL ;
   } else {
/*
      ----------------------------------------------------
      we must lock the list to return the head of the list
      ----------------------------------------------------
*/
      Lock_lock(list->lock) ;
      mtx = list->heads[ilist] ;
      list->heads[ilist] = NULL ;
      list->nlocks++ ;
      Lock_unlock(list->lock) ;
   }
} else {
   mtx = NULL ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n SubMtxList %p : mtx %p", list, mtx) ;
fflush(stdout) ;
#endif
return(mtx) ; }

/*--------------------------------------------------------------------*/
