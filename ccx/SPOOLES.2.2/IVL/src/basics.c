/*  basics.c  */

#include "../IVL.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 95sep22, cca
   -----------------------
*/
IVL *
IVL_new ( 
   void 
) {
IVL   *ivl ;

ALLOCATE(ivl, struct _IVL, 1) ;
IVL_setDefaultFields(ivl) ;

return(ivl) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 95sep22, cca
   -----------------------
*/
void
IVL_setDefaultFields (
   IVL   *ivl
) {
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_setDefaultFields(%p)"
           "\n bad input", ivl) ;
   exit(-1) ;
}
ivl->type     = IVL_NOTYPE ;
ivl->maxnlist = 0          ;
ivl->nlist    = 0          ;
ivl->tsize    = 0          ;
ivl->sizes    = NULL       ;
ivl->p_vec    = NULL       ;
ivl->incr     = IVL_INCR   ;
ivl->chunk   = NULL       ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 95sep22, cca
   --------------------------------------------------
*/
void
IVL_clearData ( 
   IVL   *ivl 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_clearData(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
/*
   ----------------------------------------------------
   switch over the storage type to free list entries.
   action is taken when type is IVL_SOLO or IVL_CHUNKED
   ----------------------------------------------------
*/
switch ( ivl->type ) {
case IVL_SOLO : {
   int   ilist ;
   for ( ilist = 0 ; ilist < ivl->nlist ; ilist++ ) {
      if ( ivl->p_vec[ilist] != NULL ) {
         IVfree(ivl->p_vec[ilist]) ;
         ivl->p_vec[ilist] = NULL ;
         ivl->tsize -= ivl->sizes[ilist] ;
      }
   }
   } break ;
case IVL_CHUNKED : {
   Ichunk   *chunk ;
   while ( (chunk = ivl->chunk) != NULL ) {
      ivl->chunk = chunk->next ;
      if ( chunk->base != NULL ) {
         IVfree(chunk->base) ;
         chunk->base = NULL ;
      }
      FREE(chunk) ;
   }
   } break ;
case IVL_NOTYPE  :
case IVL_UNKNOWN :
   break ;
default :
   fprintf(stderr, "\n fatal error in IVL_clearData(%p)"
           "\n invalid type = %d\n", ivl, ivl->type) ;
   exit(-1) ;
}
/*
   -----------------------------------------------
   free storage for the sizes[] and p_vec[] arrays
   -----------------------------------------------
*/
if ( ivl->sizes != NULL ) {
   IVfree(ivl->sizes) ;
   ivl->sizes = NULL ;
}
if ( ivl->p_vec != NULL ) {
   PIVfree(ivl->p_vec) ;
   ivl->p_vec = NULL ;
}
ivl->nlist = ivl->maxnlist = 0 ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
IVL_setDefaultFields(ivl) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 95sep22, cca
   ------------------------------------------
*/
IVL *
IVL_free ( 
   IVL   *ivl 
) {
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_free(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
IVL_clearData(ivl) ;
FREE(ivl) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
