/*  listmanip.c  */

#include "../IVL.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   fills size of and pointer to the entries of a list

   set *psize = size of list ilist
   set *pivec = base address of list ilist

   use as follows :

   IVL_listAndSize(ivl, ilist, &isize, &ivec) ;
   for ( i = 0 ; i < isize ; i++ ) {
      do something with ivec[i] ;
   }

   created -- 95sep22, cca
   ----------------------------------------------------------------
*/
void
IVL_listAndSize (
   IVL   *ivl,
   int   ilist,
   int   *psize,
   int   **pivec
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || ilist < 0 || ilist >= ivl->nlist
     || psize == NULL || pivec == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_listAndSize(%p,%d,%p,%p)"
           "\n bad input\n", ivl, ilist, psize, pivec) ;
   if ( ivl != NULL ) {
      fprintf(stderr, "\n ilist = %d, nlist = %d",
              ilist, ivl->nlist) ;
      IVL_writeForHumanEye(ivl, stderr) ;
   }
   exit(-1) ;
}
/*
   --------------------------
   set the two pointer fields
   --------------------------
*/
*psize = ivl->sizes[ilist] ;
*pivec = ivl->p_vec[ilist] ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   returns a pointer to the first element in list ilist

   to be used as a general iterator, e.g.,

   for ( pi = IVL_firstInList(ivl, ilist) ;
         pi != NULL ;
         pi = IVL_nextInList(ivl, ilist, pi) ) {
      do something ;
   }

   created -- 95sep27, cca
   ----------------------------------------------------
*/
int *
IVL_firstInList (
   IVL   *ivl,
   int   ilist
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_firstInList(%p,%d)"
           "\n bad input, ivl is NULL\n", ivl, ilist) ;
   exit(-1) ;
}
if ( ilist < 0 || ilist >= ivl->nlist ) {
   fprintf(stderr, "\n fatal error in IVL_firstInList(%p,%d)"
           "\n bad input, ilist = %d, must be in [0,%d) \n", 
           ivl, ilist, ilist, ivl->nlist) ;
   exit(-1) ;
}
if ( ivl->sizes[ilist] == 0 ) {
   return(NULL) ;
} else if ( ivl->p_vec[ilist] != NULL ) {
   return(ivl->p_vec[ilist]) ;
} else {
   fprintf(stderr, "\n fatal error in IVL_firstInList(%p,%d)"
           "\n size > 0 but list is NULL\n", ivl, ilist) ;
   exit(-1) ;
}
}
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   returns a pointer to the next element in list ilist

   to be used as a general iterator, e.g.,

   for ( pi = IVL_firstInList(ivl, ilist) ;
         pi != NULL ;
         pi = IVL_nextInList(ivl, ilist, pi) ) {
      do something ;
   }

   created -- 95sep27, cca
   ----------------------------------------------------
*/
int *
IVL_nextInList (
   IVL   *ivl,
   int   ilist,
   int   *pi
) {
int   offset ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_nextInList(%p,%d,%p)"
           "\n bad input, ivl is NULL\n", ivl, ilist, pi) ;
   exit(-1) ;
}
if ( ilist < 0 || ilist >= ivl->nlist ) {
   fprintf(stderr, "\n fatal error in IVL_nextInList(%p,%d,%p)"
           "\n bad input, ilist = %d, must be in [0,%d) \n", 
           ivl, ilist, pi, ilist, ivl->nlist) ;
   exit(-1) ;
}
if (  (pi == NULL) 
   || ((offset = pi - ivl->p_vec[ilist]) < 0)
   || offset >= ivl->sizes[ilist] ) {
   fprintf(stderr, "\n fatal error in IVL_nextInList(%p,%d,%p)"
           "\n bad pointer\n", ivl, ilist, pi) ;
   exit(-1) ;
} else if ( offset == ivl->sizes[ilist] - 1 ) {
   return(NULL) ;
} else {
   return(pi+1) ;
}
}
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to set or reset a list.

   ilist -- list id to set or reset
   isize -- size of the list
      if the present size of list ilist is smaller than isize,
         the old list is free'd (if ivl->type = IVL_SOLO) 
         or lost (if ivl->type = IVL_CHUNKED)
         or un-set (if ivl->type = IVL_UNKNOWN)
         and new storage is allocated (for IVL_SOLO and IVL_CHUNKED)
   ivec  -- list vector
      if ivl->type is IVL_UNKNOWN then
         if ivec != NULL then
            we set the ivl list pointer to be ivec
         endif
      else if ivec != NULL
         we copy ivec[] into ivl's storage for the list
      endif

   created   -- 95sep27, cca
   last mods -- 95oct06, cca
      type = IVL_UNKNOWN, p_vec[ilist] set to ivec 
      only when ivec is not NULL
      bug fixed, ivl->sizes[ilist] = isize ;
   -----------------------------------------------------------------
*/
void
IVL_setList ( 
   IVL   *ivl, 
   int   ilist, 
   int   isize,
   int   ivec[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_setList(%p,%d,%d,%p)"
           "\n bad input, ivl is NULL\n", ivl, ilist, isize, ivec) ;
   exit(-1) ;
}
if ( ilist < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_setList(%p,%d,%d,%p)"
           "\n bad input, ilist < 0", ivl, ilist, isize, ivec) ;
   exit(-1) ;
}
if ( ilist >= ivl->maxnlist ) {
   int   newmaxnlist = (int) 1.25*ivl->maxnlist ;
   if ( newmaxnlist < 10 ) {
      newmaxnlist = 10 ;
   } 
   if ( ilist >= newmaxnlist ) {
      newmaxnlist = ilist + 1 ;
   } 
   IVL_setMaxnlist(ivl, newmaxnlist) ;
}
if ( ilist >= ivl->nlist ) {
   ivl->nlist = ilist + 1 ;
}
if ( isize == 0 ) {
/*
   ------------------------------------------
   new list is empty, free storage if present
   ------------------------------------------
*/
   if ( ivl->type == IVL_SOLO ) {
      if ( ivl->p_vec[ilist] != NULL ) {
         IVfree(ivl->p_vec[ilist]) ;
      }
   }
   ivl->tsize -= ivl->sizes[ilist] ;
   ivl->sizes[ilist] =   0  ;
   ivl->p_vec[ilist] = NULL ;
} else if ( ivl->type == IVL_UNKNOWN ) {
/*
   -------------------------------
   simply set the size and pointer
   -------------------------------
*/
   ivl->tsize += isize - ivl->sizes[ilist] ;
   ivl->sizes[ilist] = isize ;
   if ( ivec != NULL ) {
      ivl->p_vec[ilist] = ivec ;
   }
} else {
/*
   --------------------------------------------------
   the list entries will be copied into ivl's storage
   --------------------------------------------------
*/
   if ( ivl->sizes[ilist] < isize ) {
/*
      --------------------------------------------
      old size (might be zero) is not large enough
      switch over the storage types
      --------------------------------------------
*/
      switch ( ivl->type ) {
      case IVL_SOLO :
         if ( ivl->p_vec[ilist] != NULL ) {
/*
            ---------------------------------------
            free old storage before allocating more
            and decrement the total list size
            ---------------------------------------
*/
            IVfree(ivl->p_vec[ilist]) ;
         }
         ivl->p_vec[ilist] = IVinit(isize, -1) ;
         break ;
      case IVL_CHUNKED : {
         Ichunk   *chunk ;
         if (  (chunk = ivl->chunk) == NULL
            || (chunk->size - chunk->inuse) < isize ) {
/*
            -------------------------------
            allocate a new chunk of storage
            -------------------------------
*/
            ALLOCATE(chunk, struct _Ichunk, 1) ;
            if ( isize < ivl->incr ) {
               chunk->size = ivl->incr ;
            } else {
               chunk->size = isize ;
            }
            chunk->inuse = 0 ;
            chunk->base  = IVinit(chunk->size, -1) ;
            chunk->next  = ivl->chunk ;
            ivl->chunk   = chunk ;
         }
/*
         --------------------------
         set pointer for list ilist
         --------------------------
*/
         ivl->p_vec[ilist] = chunk->base + chunk->inuse ;
         chunk->inuse += isize ;
         } break ;
      default :
         fprintf(stderr, "\n fatal error in IVL_setList(%p,%d,%d,%p)"
               "\n you are trying to grow a list but type = %d"
               "\n type must be IVL_CHUNKED = 1 or IVL_SOLO = 2\n",
               ivl, ilist, isize, ivec, ivl->type) ;
         exit(-1) ;
      }
   }
   ivl->tsize += isize - ivl->sizes[ilist] ;
   ivl->sizes[ilist] = isize ;
   if ( ivec != NULL ) {
/*
      --------------------------------
      copy the list into ivl's storage
      --------------------------------
*/
      IVcopy(isize, ivl->p_vec[ilist], ivec) ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   set a pointer to a list but don't allocate storage for the list.
   this method was needed when we form a subgraph with a boundary.
   lists for the interior vertices point into the parent graph,
   but lists for the boundary vertices must be allocated and owned.
   used only for type = IVL_CHUNKED. at some point in the future we
   should rethink the storage semantics for the IVL object.

   created -- 95nov11, cca
   ----------------------------------------------------------------
*/
void
IVL_setPointerToList (
   IVL   *ivl, 
   int   ilist, 
   int   isize,
   int   ivec[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_setPointerToList(%p,%d,%d,%p)"
           "\n bad input, ivl is NULL\n", ivl, ilist, isize, ivec) ;
   exit(-1) ;
}
if ( ivl->type != IVL_CHUNKED ) {
   fprintf(stderr, "\n fatal error in IVL_setPointerToList(%p,%d,%d,%p)"
           "\n this method is only used with type IVL_CHUNKED\n", 
           ivl, ilist, isize, ivec) ;
   exit(-1) ;
}
if ( ilist < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_setPointerToList(%p,%d,%d,%p)"
           "\n bad input, ilist < 0", ivl, ilist, isize, ivec) ;
   exit(-1) ;
}
if ( ilist >= ivl->maxnlist ) {
   int   newmaxnlist = (int) 1.25*ivl->maxnlist ;
   if ( newmaxnlist < 10 ) {
      newmaxnlist = 10 ;
   } 
   if ( ilist >= newmaxnlist ) {
      newmaxnlist = ilist + 1 ;
   } 
   IVL_setMaxnlist(ivl, newmaxnlist) ;
}
if ( ilist >= ivl->nlist ) {
   ivl->nlist = ilist + 1 ;
}
if ( ivl->type == IVL_SOLO && ivl->p_vec[ilist] != NULL ) {
   IVfree(ivl->p_vec[ilist]) ;
}
ivl->tsize += isize - ivl->sizes[ilist] ;
ivl->sizes[ilist] = isize ;
ivl->p_vec[ilist] = ivec  ;

return ; }

/*--------------------------------------------------------------------*/
