/*  basics.c  */

#include "../ChvManager.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98may02, cca
   -----------------------
*/
ChvManager *
ChvManager_new ( 
   void 
) {
ChvManager   *manager ;

ALLOCATE(manager, struct _ChvManager, 1) ;
ChvManager_setDefaultFields(manager) ;

return(manager) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may02, cca
   -----------------------
*/
void
ChvManager_setDefaultFields (
   ChvManager   *manager
) {
if ( manager == NULL ) {
   fprintf(stderr, "\n fatal error in ChvManager_setDefaultFields(%p)"
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
ChvManager_clearData ( 
   ChvManager   *manager 
) {
Chv   *chv ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL ) {
   fprintf(stderr, "\n fatal error in ChvManager_clearData(%p)"
           "\n bad input\n", manager) ;
   exit(-1) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
while ( (chv = manager->head) != NULL ) {
   manager->head = chv->next ;
   Chv_free(chv) ;
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
ChvManager_setDefaultFields(manager) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98may02, cca
   ------------------------------------------
*/
void
ChvManager_free ( 
   ChvManager   *manager 
) {
if ( manager == NULL ) {
   fprintf(stderr, "\n fatal error in ChvManager_free(%p)"
           "\n bad input\n", manager) ;
   exit(-1) ;
}
ChvManager_clearData(manager) ;
FREE(manager) ;

return ; }

/*--------------------------------------------------------------------*/
