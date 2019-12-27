/*  basics.c  */

#include "../ChvList.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may02, cca
   -----------------------
*/
ChvList *
ChvList_new ( 
   void 
) {
ChvList   *chvlist ;

ALLOCATE(chvlist, struct _ChvList, 1) ;
ChvList_setDefaultFields(chvlist) ;

return(chvlist) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may02, cca
   -----------------------
*/
void
ChvList_setDefaultFields (
   ChvList   *chvlist
) {
if ( chvlist == NULL ) {
   fprintf(stderr, "\n fatal error in ChvList_setDefaultFields(%p)"
           "\n bad input", chvlist) ;
   exit(-1) ;
}
chvlist->nlist  =   0  ;
chvlist->heads  = NULL ;
chvlist->counts = NULL ;
chvlist->lock   = NULL ;
chvlist->flags  = NULL ;
chvlist->nlocks =   0  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may02, cca
   --------------------------------------------------
*/
void
ChvList_clearData ( 
   ChvList   *chvlist 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chvlist == NULL ) {
   fprintf(stderr, "\n fatal error in ChvList_clearData(%p)"
           "\n bad input\n", chvlist) ;
   exit(-1) ;
}
/*
   -------------
   free the data
   -------------
*/
if ( chvlist->heads != NULL ) {
   FREE(chvlist->heads) ;
}
if ( chvlist->counts != NULL ) {
   IVfree(chvlist->counts) ;
}
if ( chvlist->flags != NULL ) {
   CVfree(chvlist->flags) ;
}
if ( chvlist->lock != NULL ) {
/*
   -------------------------
   destroy and free the lock
   -------------------------
*/
   Lock_free(chvlist->lock) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
ChvList_setDefaultFields(chvlist) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may02, cca
   ------------------------------------------
*/
void
ChvList_free ( 
   ChvList   *chvlist 
) {
if ( chvlist == NULL ) {
   fprintf(stderr, "\n fatal error in ChvList_free(%p)"
           "\n bad input\n", chvlist) ;
   exit(-1) ;
}
ChvList_clearData(chvlist) ;
FREE(chvlist) ;

return ; }

/*--------------------------------------------------------------------*/
