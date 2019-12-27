/*  basics.c  */

#include "../Lock.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 97aug22, cca
   -----------------------
*/
Lock *
Lock_new ( 
   void 
) {
Lock   *lock ;

ALLOCATE(lock, struct _Lock, 1) ;
Lock_setDefaultFields(lock) ;

return(lock) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 97aug22, cca
   -----------------------
*/
void
Lock_setDefaultFields (
   Lock   *lock
) {
if ( lock == NULL ) {
   fprintf(stderr, "\n fatal error in Lock_setDefaultFields(%p)"
           "\n bad input", lock) ;
   exit(-1) ;
}
lock->nlocks   =   0  ;
lock->nunlocks =   0  ;
lock->mutex    = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 97aug22, cca
   --------------------------------------------------
*/
void
Lock_clearData ( 
   Lock   *lock 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( lock == NULL ) {
   fprintf(stderr, "\n fatal error in Lock_clearData(%p)"
           "\n bad input\n", lock) ;
   exit(-1) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( lock->mutex != NULL ) {
/*
   -------------------------
   destroy and free the lock
   -------------------------
*/
#if THREAD_TYPE == TT_SOLARIS
   mutex_destroy(lock->mutex) ;
#endif
#if THREAD_TYPE == TT_POSIX
   pthread_mutex_destroy(lock->mutex) ;
#endif
   FREE(lock->mutex) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
Lock_setDefaultFields(lock) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 97aug22, cca
   ------------------------------------------
*/
void
Lock_free ( 
   Lock   *lock 
) {
if ( lock == NULL ) {
   fprintf(stderr, "\n fatal error in Lock_free(%p)"
           "\n bad input\n", lock) ;
   exit(-1) ;
}
Lock_clearData(lock) ;
FREE(lock) ;

return ; }

/*--------------------------------------------------------------------*/
