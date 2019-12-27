/*  util.c  */

#include "../SubMtxManager.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   return a pointer to a SubMtx object that has 
   been initialized with the input parameters

   created -- 98may02, cca
   -----------------------------------------------
*/
SubMtx *
SubMtxManager_newObjectOfSizeNbytes (
   SubMtxManager   *manager,
   int             nbytesNeeded
) {
SubMtx   *mtx, *prev ;
int    nbytesAvailable, newinstance ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || nbytesNeeded <= 0 ) {
   fprintf(stderr, 
      "\n fatal error in SubMtxMananger_newObjectOfSizeNbytes(%p,%d)"
      "\n bad input\n", manager, nbytesNeeded) ;
   exit(-1) ;
}
#if MYDEBUG > 1
   fprintf(stdout, 
           "\n\n new mtx request, nbytes needed = %d", nbytesNeeded) ;
   fflush(stdout) ;
#endif
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------
   lock the lock, get exclusive access to the list
   -----------------------------------------------
*/
   Lock_lock(manager->lock) ;
#if MYDEBUG > 1
   fprintf(stdout, "\n manager: lock is locked") ;
   fflush(stdout) ;
#endif
   manager->nlocks++ ;
#if MYDEBUG > 1
   fprintf(stdout, "\n %d locks so far", manager->nlocks) ;
   fflush(stdout) ;
#endif
}
/*
   ---------------------------------------------------------
   find a SubMtx object with the required number of bytes
   ---------------------------------------------------------
*/
for ( mtx = manager->head, prev = NULL ; 
      mtx != NULL ; 
      mtx = mtx->next ) {
   nbytesAvailable = SubMtx_nbytesInWorkspace(mtx) ;
#if MYDEBUG > 1
   fprintf(stdout, "\n free mtx %p, nbytes = %d",
           mtx, nbytesAvailable) ;
   fflush(stdout) ;
#endif
   if ( nbytesNeeded <= nbytesAvailable ) {
      break ;
   }
   prev = mtx ;
}
if ( mtx != NULL ) {
/*
   ---------------------------------------
   suitable object found, remove from list
   ---------------------------------------
*/
#if MYDEBUG > 1
   fprintf(stdout, "\n mtx = %p, %d nbytes available",
           mtx, nbytesAvailable) ;
   fflush(stdout) ;
#endif
   if ( prev == NULL ) {
      manager->head = mtx->next ;
   } else {
      prev->next = mtx->next ;
   }
   newinstance = 0 ;
} else {
/*
   ------------------------------------------------------------------
   no suitable object found, create new instance and allocate storage
   ------------------------------------------------------------------
*/
   mtx = SubMtx_new() ;
#if MYDEBUG > 1
   fprintf(stdout, 
           "\n no suitable object found, new mtx = %p, bytes = %d", 
           mtx, nbytesNeeded) ;
   fflush(stdout) ;
#endif
   newinstance = 1 ;
   DV_setSize(&mtx->wrkDV, nbytesNeeded/sizeof(double)) ;
}
if ( newinstance == 1 ) {
   manager->nbytesalloc  += SubMtx_nbytesInWorkspace(mtx) ;
}
manager->nactive++ ;
manager->nbytesactive    += SubMtx_nbytesInWorkspace(mtx) ;
manager->nbytesrequested += nbytesNeeded ;
#if MYDEBUG > 1
fprintf(stdout, "\n %d bytes active, %d bytes requested",
        manager->nbytesactive, manager->nbytesrequested) ;
#endif
manager->nrequests++ ;
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------------
   unlock the lock, release exclusive access to the list
   -----------------------------------------------------
*/
   manager->nunlocks++ ;
#if MYDEBUG > 1
   fprintf(stdout, "\n manager: unlocking, %d unlocks so far",
           manager->nunlocks) ;
   fflush(stdout) ;
#endif
   Lock_unlock(manager->lock) ;
}
return(mtx) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   release a SubMtx instance

   created -- 98may02, cca
   ----------------------------
*/
void
SubMtxManager_releaseObject (
   SubMtxManager   *manager,
   SubMtx          *mtx1
) {
SubMtx   *mtx2, *prev ;
int    nbytes1, nbytes2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || mtx1 == NULL ) {
   fprintf(stderr, 
       "\n fatal error in SubMtxManager_releaseObject(%p,%p)"
       "\n bad input\n", manager, mtx1) ;
   exit(-1) ;
}
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------
   lock the lock, get exclusive access to the list
   -----------------------------------------------
*/
   Lock_lock(manager->lock) ;
   manager->nlocks++ ;
#if MYDEBUG > 1
   fprintf(stdout, "\n\n manager : locking in releaseObject, %d locks",
           manager->nlocks) ;
   fflush(stdout) ;
#endif
}
manager->nreleases++ ;
manager->nbytesactive -= SubMtx_nbytesInWorkspace(mtx1) ;
manager->nactive-- ;
if ( manager->mode == 0 ) {
/*
   ---------------
   release storage
   ---------------
*/
   SubMtx_free(mtx1) ;
} else {
/*
   --------------------------------------------------------
   find a place in the list where the SubMtx objects are 
   sorted in ascending order of the size of their workspace
   --------------------------------------------------------
*/
   nbytes1 = SubMtx_nbytesInWorkspace(mtx1) ;
#if MYDEBUG > 1
   fprintf(stdout, "\n\n trying to release mtx %p with %d bytes",
           mtx1, nbytes1) ;
#endif
   for ( mtx2 = manager->head, prev = NULL ; 
         mtx2 != NULL ; 
         mtx2 = mtx2->next ) {
      nbytes2 = SubMtx_nbytesInWorkspace(mtx2) ;
#if MYDEBUG > 1
      fprintf(stdout, "\n list mtx %p with %d bytes", mtx2, nbytes2) ;
#endif
      if ( nbytes2 >= nbytes1 ) {
         break ;
      }
      prev = mtx2 ;
   }
   if ( prev == NULL ) {
      manager->head = mtx1 ;
#if MYDEBUG > 1
      fprintf(stdout, "\n manager->head = %p", mtx1) ;
#endif
   } else {
      prev->next = mtx1 ;
#if MYDEBUG > 1
      fprintf(stdout, "\n %p->next = %p", prev, mtx1) ;
#endif
   }
   mtx1->next = mtx2 ;
#if MYDEBUG > 1
   fprintf(stdout, "\n %p->next = %p", mtx1, mtx2) ;
#endif
}
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------------
   unlock the lock, release exclusive access to the list
   -----------------------------------------------------
*/
   manager->nunlocks++ ;
#if MYDEBUG > 1
   fprintf(stdout, 
           "\n manager : unlocking in releaseObject, %d unlocks",
           manager->nunlocks) ;
   fflush(stdout) ;
#endif
   Lock_unlock(manager->lock) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   release a list of SubMtx objects

   created -- 98may02, cca
   -----------------------------------
*/
void
SubMtxManager_releaseListOfObjects (
   SubMtxManager   *manager,
   SubMtx          *head
) {
SubMtx   *mtx1, *mtx2, *prev ;
int    nbytes1, nbytes2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || head == NULL ) {
   fprintf(stderr, 
       "\n fatal error in SubMtxManager_releaseListOfObjects(%p,%p)"
       "\n bad input\n", manager, head) ;
   exit(-1) ;
}
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------
   lock the lock, get exclusive access to the list
   -----------------------------------------------
*/
   Lock_lock(manager->lock) ;
   manager->nlocks++ ;
#if MYDEBUG > 1
   fprintf(stdout, 
           "\n\n manager : locking in releaseListOfObjects, %d locks",
           manager->nlocks) ;
   fflush(stdout) ;
#endif
}
if ( manager->mode == 0 ) {
/*
   ---------------
   release storage
   ---------------
*/
   while ( (mtx1 = head) != NULL ) {
      head = head->next ;
      manager->nbytesactive -= SubMtx_nbytesInWorkspace(mtx1) ;
      manager->nactive-- ;
      manager->nreleases++ ;
      SubMtx_free(mtx1) ;
   }
} else {
/*
   -------------------
   recycle the objects
   -------------------
*/
   while ( head != NULL ) {
      mtx1 = head ;
      head = mtx1->next ;
/*
      --------------------------------------------------------
      find a place in the list where the SubMtx objects are 
      sorted in ascending order of the size of their workspace
      --------------------------------------------------------
*/
      nbytes1 = SubMtx_nbytesInWorkspace(mtx1) ;
#if MYDEBUG > 1
      fprintf(stdout, "\n\n trying to release mtx %p with %d bytes",
              mtx1, nbytes1) ;
#endif
      for ( mtx2 = manager->head, prev = NULL ; 
            mtx2 != NULL ; 
            mtx2 = mtx2->next ) {
         nbytes2 = SubMtx_nbytesInWorkspace(mtx2) ;
#if MYDEBUG > 1
         fprintf(stdout, 
                 "\n list mtx %p with %d bytes", mtx2, nbytes2) ;
#endif
         if ( nbytes2 >= nbytes1 ) {
            break ;
         }
         prev = mtx2 ;
      }
      if ( prev == NULL ) {
         manager->head = mtx1 ;
#if MYDEBUG > 1
         fprintf(stdout, "\n manager->head = %p", mtx1) ;
#endif
      } else {
         prev->next = mtx1 ;
#if MYDEBUG > 1
         fprintf(stdout, "\n %p->next = %p", prev, mtx1) ;
#endif
      }
      mtx1->next = mtx2 ;
#if MYDEBUG > 1
      fprintf(stdout, "\n %p->next = %p", mtx1, mtx2) ;
#endif
      manager->nbytesactive -= SubMtx_nbytesInWorkspace(mtx1) ;
      manager->nactive-- ;
      manager->nreleases++ ;
#if MYDEBUG > 1
      fprintf(stdout, "\n # releases = %d", manager->nreleases) ;
      for ( mtx1 = manager->head ; mtx1 != NULL ; mtx1 = mtx1->next ) {
         fprintf(stdout, "\n mtx (%d,%d), nbytes %d",
                 mtx1->rowid, mtx1->colid, 
                 SubMtx_nbytesInWorkspace(mtx1)) ;
      }
#endif
   }
}
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------------
   unlock the lock, release exclusive access to the list
   -----------------------------------------------------
*/
   manager->nunlocks++ ;
#if MYDEBUG > 1
   fprintf(stdout, 
           "\n manager : unlocking in releaseListOfObjects, %d unlocks",
           manager->nunlocks) ;
   fflush(stdout) ;
#endif
   Lock_unlock(manager->lock) ;
}
return ; }

/*--------------------------------------------------------------------*/
