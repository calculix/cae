/*  init.c  */

#include "../IVL.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   initialize given the type and maximum number of lists.
   used for type IVL_CHUNKED, IVL_SOLO or IVL_UNKNOWN

   created -- 95sep22, cca
   ------------------------------------------------------
*/
void
IVL_init1 ( 
   IVL   *ivl, 
   int   type, 
   int   maxnlist
) {
/*
   -------------------
   check for bad input
   -------------------
*/
if ( ivl == NULL 
  || (type != IVL_CHUNKED && type != IVL_SOLO && type != IVL_UNKNOWN)
  || maxnlist < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_init1(%p,%d,%d)"
           "\n bad input", ivl, type, maxnlist) ;
   exit(-1) ;
}
/*
   --------------------------
   clear the data, if present
   --------------------------
*/
IVL_clearData(ivl) ;
/*
   ----------
   initialize
   ----------
*/
ivl->type     = type  ;
ivl->maxnlist = maxnlist ;
ivl->nlist    = maxnlist ;
if ( maxnlist > 0 ) {
   ivl->sizes = IVinit(maxnlist, 0) ;
   ivl->p_vec = PIVinit(maxnlist) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   initialize given the type, number of lists and their total size
   only used when type is IVL_CHUNKED.

   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
void
IVL_init2 ( 
   IVL   *ivl, 
   int   type, 
   int   maxnlist, 
   int   tsize 
) {
/*
   -------------------
   check for bad input
   -------------------
*/
if ( ivl == NULL || type != IVL_CHUNKED || maxnlist < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_init2(%p,%d,%d,%d)"
           "\n bad input", ivl, type, maxnlist, tsize) ;
   exit(-1) ;
}
/*
   ---------------------------------
   initialize via IVL_init1() method
   ---------------------------------
*/
IVL_init1(ivl, type, maxnlist) ;
/*
   ----------------------------------
   create chunk to hold tsize entries
   ----------------------------------
*/
if ( tsize > 0 ) {
   ALLOCATE(ivl->chunk, struct _Ichunk, 1) ;
   ivl->chunk->size  = tsize ;
   ivl->chunk->inuse = 0 ;
   ivl->chunk->base  = IVinit(tsize, -1) ;
   ivl->chunk->next  = NULL ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   initialize from a vector of list sizes.
   used with IVL_SOLO or IVL_CHUNKED.

   created -- 95sep22, cca
   --------------------------------------
*/
void
IVL_init3 ( 
   IVL   *ivl, 
   int   type, 
   int   maxnlist, 
   int   sizes[] 
) {
int   ilist ;
/*
   -------------------
   check for bad input
   -------------------
*/
if ( ivl == NULL || (type != IVL_CHUNKED && type != IVL_SOLO)
  || maxnlist < 0 || sizes == NULL ) {
   fprintf(stderr,
           "\n fatal error in IVL_init3(%p,%d,%d,%p)"
           "\n bad input", ivl, type, maxnlist, sizes) ;
   exit(-1) ;
}
switch ( type ) {
case IVL_SOLO :
/*
   ---------------------------------
   initialize via IVL_init1() method
   ---------------------------------
*/
   IVL_init1(ivl, type, maxnlist) ;
   break ;
case IVL_CHUNKED :
/*
   ---------------------------------
   initialize via IVL_init2() method
   ---------------------------------
*/
   IVL_init2(ivl, type, maxnlist, IVsum(maxnlist, sizes)) ;
   break ;
}
/*
   -------------------------
   set the size of each list
   -------------------------
*/
for ( ilist = 0 ; ilist < maxnlist ; ilist++ ) {
   IVL_setList(ivl, ilist, sizes[ilist], NULL) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   this method resizes the maximum number of lists,
   replacing the old sizes[] and p_vec[] vectors
   as necessary. the nlist value is NOT reset.

   created -- 96dec05, cca
   ------------------------------------------------
*/
void
IVL_setMaxnlist (
   IVL   *ivl,
   int   newmaxnlist
) {
if ( ivl == NULL || newmaxnlist < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_setMaxnlist(%p,%d)"
           "\n bad input\n", ivl, newmaxnlist) ;
   exit(-1) ;
}
if ( newmaxnlist != ivl->maxnlist ) {
   int   *ivec, **pivec ;
/*
   --------------------------------------------
   allocate (and fill) the new sizes[] array
   --------------------------------------------
*/
   ivec = IVinit(newmaxnlist, 0) ;
   if ( ivl->sizes != NULL ) {
      if ( ivl->nlist > newmaxnlist ) {
         IVcopy(newmaxnlist, ivec, ivl->sizes) ;
      } else if ( ivl->nlist > 0 ) {
         IVcopy(ivl->nlist, ivec, ivl->sizes) ;
      }
      IVfree(ivl->sizes) ;
   }
   ivl->sizes = ivec ;
/*
   --------------------------------------------
   allocate (and fill) the larger p_vec[] array
   --------------------------------------------
*/
   pivec = PIVinit(newmaxnlist) ;
   if ( ivl->p_vec != NULL ) {
      if ( ivl->nlist > newmaxnlist ) {
         PIVcopy(newmaxnlist, pivec, ivl->p_vec) ;
      } else if ( ivl->nlist > 0 ) {
         PIVcopy(ivl->nlist, pivec, ivl->p_vec) ;
      }
      PIVfree(ivl->p_vec) ;
   }
   ivl->p_vec = pivec ;
/*
   -----------------------------------
   set the new maximum number of lists
   -----------------------------------
*/
   ivl->maxnlist = newmaxnlist ;
   if ( ivl->nlist > newmaxnlist ) {
      ivl->nlist = newmaxnlist ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   this method resizes the number of lists,
   replacing the old sizes[] and p_vec[] vectors
   as necessary. 

   created -- 96dec05, cca
   ------------------------------------------------
*/
void
IVL_setNlist (
   IVL   *ivl,
   int   newnlist
) {
if ( ivl == NULL || newnlist < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_setNlist(%p,%d)"
           "\n bad input\n", ivl, newnlist) ;
   exit(-1) ;
}
if ( newnlist > ivl->maxnlist ) {
/*
   ------------------------------------
   increase the maximum number of lists
   ------------------------------------
*/
   IVL_setMaxnlist(ivl, newnlist) ;
}
/*
   -------------------
   set the nlist field
   -------------------
*/
ivl->nlist = newnlist ;

return ; }

/*--------------------------------------------------------------------*/
