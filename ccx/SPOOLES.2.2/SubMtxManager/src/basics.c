/*  basics.c  */

#include "../SubMtxManager.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may02, cca
   -----------------------
*/
SubMtxManager *
SubMtxManager_new ( 
   void 
) {
SubMtxManager   *manager ;

ALLOCATE(manager, struct _SubMtxManager, 1) ;
SubMtxManager_setDefaultFields(manager) ;

return(manager) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may02, cca
   -----------------------
*/
void
SubMtxManager_setDefaultFields (
   SubMtxManager   *manager
) {
if ( manager == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxManager_setDefaultFields(%p)"
           "\n bad input", manager) ;
   exit(-1) ;
}
manager->head            = NULL ;
manager->lock            = NULL ;
manager->mode            =   0  ;
manager->nactive         =   0  ;
manager->nbytesactive    =   0  ;
manager->nbytesrequested =   0  ;
manager->nbytesalloc     =   0  ;
manager->nrequests       =   0  ;
manager->nreleases       =   0  ;
manager->nlocks          =   0  ;
manager->nunlocks        =   0  ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98may02, cca
   --------------------------------------------------
*/
void
SubMtxManager_clearData ( 
   SubMtxManager   *manager 
) {
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtxManager_clearData(%p)"
           "\n bad input\n", manager) ;
   exit(-1) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
while ( (mtx = manager->head) != NULL ) {
   manager->head = mtx->next ;
   SubMtx_free(mtx) ;
}
if ( manager->lock != NULL ) {
/*
   -------------------------
   destroy and free the lock
   -------------------------
*/
   Lock_free(manager->lock) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
SubMtxManager_setDefaultFields(manager) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may02, cca
   ------------------------------------------
*/
void
SubMtxManager_free ( 
   SubMtxManager   *manager 
) {
if ( manager == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtxManager_free(%p)"
           "\n bad input\n", manager) ;
   exit(-1) ;
}
SubMtxManager_clearData(manager) ;
FREE(manager) ;

return ; }

/*--------------------------------------------------------------------*/
