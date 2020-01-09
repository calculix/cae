/*  basics.c  */

#include "../SubMtxList.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may02, cca
   -----------------------
*/
SubMtxList *
SubMtxList_new ( 
   void 
) {
SubMtxList   *list ;

ALLOCATE(list, struct _SubMtxList, 1) ;
SubMtxList_setDefaultFields(list) ;

return(list) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may02, cca
   -----------------------
*/
void
SubMtxList_setDefaultFields (
   SubMtxList   *list
) {
if ( list == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxList_setDefaultFields(%p)"
           "\n bad input", list) ;
   exit(-1) ;
}
list->nlist  =   0  ;
list->heads  = NULL ;
list->counts = NULL ;
list->lock   = NULL ;
list->flags  = NULL ;
list->nlocks =   0  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may02, cca
   --------------------------------------------------
*/
void
SubMtxList_clearData ( 
   SubMtxList   *list 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( list == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtxList_clearData(%p)"
           "\n bad input\n", list) ;
   exit(-1) ;
}
/*
   -------------
   free the data
   -------------
*/
if ( list->heads != NULL ) {
   FREE(list->heads) ;
}
if ( list->counts != NULL ) {
   IVfree(list->counts) ;
}
if ( list->flags != NULL ) {
   CVfree(list->flags) ;
}
if ( list->lock != NULL ) {
/*
   -------------------------
   destroy and free the lock
   -------------------------
*/
   Lock_free(list->lock) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
SubMtxList_setDefaultFields(list) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may02, cca
   ------------------------------------------
*/
void
SubMtxList_free ( 
   SubMtxList   *list 
) {
if ( list == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtxList_free(%p)"
           "\n bad input\n", list) ;
   exit(-1) ;
}
SubMtxList_clearData(list) ;
FREE(list) ;

return ; }

/*--------------------------------------------------------------------*/
