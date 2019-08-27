/*  basics.c  */

#include "../Perm.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 96jan05, cca
   -----------------------
*/
Perm *
Perm_new ( 
   void 
) {
Perm   *perm ;
ALLOCATE(perm, struct _Perm, 1) ;
Perm_setDefaultFields(perm) ;

return(perm) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 96jan05, cca
   -----------------------
*/
void
Perm_setDefaultFields (
   Perm   *perm
) {
if ( perm == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_setDefaultFields(%p)"
           "\n bad input", perm) ;
   exit(-1) ;
}
perm->isPresent =   0  ;
perm->size      =   0  ;
perm->newToOld  = NULL ;
perm->oldToNew  = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 96jan05, cca
   --------------------------------------------------
*/
void
Perm_clearData ( 
   Perm   *perm 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_clearData(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}
if ( perm->newToOld != NULL ) {
   IVfree(perm->newToOld) ;
}
if ( perm->oldToNew != NULL ) {
   IVfree(perm->oldToNew) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
Perm_setDefaultFields(perm) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 96jan05, cca
   ------------------------------------------
*/
Perm *
Perm_free ( 
   Perm   *perm 
) {
if ( perm == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_free(%p)"
           "\n bad input\n", perm) ;
   exit(-1) ;
}
Perm_clearData(perm) ;
FREE(perm) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
