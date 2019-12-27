/*  util.c  */

#include "../ChvManager.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return a pointer to a Chv object that has 
   been initialized with the input parameters

   created -- 98may02, cca
   ------------------------------------------
*/
Chv *
ChvManager_newObjectOfSizeNbytes (
   ChvManager   *manager,
   int          nbytesNeeded
) {
Chv   *chv, *prev ;
int    nbytesAvailable, newinstance ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || nbytesNeeded <= 0 ) {
   fprintf(stderr, 
           "\n fatal error in ChvMananger_newObjectOfSizeNbytes(%p,%d)"
           "\n bad input\n", manager, nbytesNeeded) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n %d bytes needed", nbytesNeeded) ;
fflush(stdout) ;
#endif
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------
   lock the lock, get exclusive access to the list
   -----------------------------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n manager: locking") ;
   fflush(stdout) ;
#endif
   Lock_lock(manager->lock) ;
   manager->nlocks++ ;
#if MYDEBUG > 0
   fprintf(stdout, ", %d locks so far", manager->nlocks) ;
   fflush(stdout) ;
#endif
}
/*
   ----------------------------------------------------
   find a Chv object with the required number of bytes
   ----------------------------------------------------
*/
for ( chv = manager->head, prev = NULL ; 
      chv != NULL ; 
      chv = chv->next ) {
   nbytesAvailable = Chv_nbytesInWorkspace(chv) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n free chev %p, nbytes = %d",
           chv, nbytesAvailable) ;
   fflush(stdout) ;
#endif
   if ( nbytesNeeded <= nbytesAvailable ) {
      break ;
   }
   prev = chv ;
}
if ( chv != NULL ) {
/*
   ---------------------------------------
   suitable object found, remove from list
   ---------------------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n chv = %p, %d nbytes available",
           chv, nbytesAvailable) ;
   fflush(stdout) ;
#endif
   if ( prev == NULL ) {
      manager->head = chv->next ;
   } else {
      prev->next = chv->next ;
   }
   newinstance = 0 ;
} else {
/*
   ------------------------------------------------------------------
   no suitable object found, create new instance and allocate storage
   ------------------------------------------------------------------
*/
   chv = Chv_new() ;
#if MYDEBUG > 0
   fprintf(stdout, "\n new postponed chv = %p", chv) ;
   fflush(stdout) ;
#endif
   newinstance = 1 ;
   DV_setSize(&chv->wrkDV, nbytesNeeded/sizeof(double)) ;
}
/*
   -------------------------------
   increment the statistics fields
   -------------------------------
*/
if ( newinstance == 1 ) {
   manager->nbytesalloc += Chv_nbytesInWorkspace(chv) ;
}
manager->nactive++ ;
manager->nbytesactive += Chv_nbytesInWorkspace(chv) ;
manager->nbytesrequested += nbytesNeeded ;
#if MYDEBUG > 0
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
   Lock_unlock(manager->lock) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n manager: unlocking, %d unlocks so far",
           manager->nunlocks) ;
   fflush(stdout) ;
#endif
}
return(chv) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   release a Chv instance

   created -- 98may02, cca
   -----------------------
*/
void
ChvManager_releaseObject (
   ChvManager   *manager,
   Chv          *chv1
) {
Chv   *chv2, *prev ;
int    nbytes1, nbytes2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || chv1 == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ChvMananger_releaseObject(%p,%p)"
           "\n bad input\n", manager, chv1) ;
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
}
manager->nreleases++ ;
manager->nbytesactive -= Chv_nbytesInWorkspace(chv1) ;
manager->nactive-- ;
if ( manager->mode == 0 ) {
/*
   -------------------
   release the storage
   -------------------
*/
   Chv_free(chv1) ;
} else {
/*
   --------------------------------------------------------
   find a place in the list where the Chv objects are 
   sorted in ascending order of the size of their workspace
   --------------------------------------------------------
*/
   nbytes1 = Chv_nbytesInWorkspace(chv1) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n trying to release chevron %p with %d bytes",
           chv1, nbytes1) ;
#endif
   for ( chv2 = manager->head, prev = NULL ; 
         chv2 != NULL ; 
         chv2 = chv2->next ) {
      nbytes2 = Chv_nbytesInWorkspace(chv2) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n list chv %p with %d bytes", chv2, nbytes2) ;
#endif
      if ( nbytes2 > nbytes1 ) {
         break ;
      }
      prev = chv2 ;
   }
   if ( prev == NULL ) {
      manager->head = chv1 ;
#if MYDEBUG > 0
      fprintf(stdout, "\n manager->head = %p", chv1) ;
#endif
   } else {
      prev->next = chv1 ;
#if MYDEBUG > 0
      fprintf(stdout, "\n %p->next = %p", prev, chv1) ;
#endif
   }
   chv1->next = chv2 ;
#if MYDEBUG > 0
   fprintf(stdout, "\n %p->next = %p", chv1, chv2) ;
#endif
}
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------------
   unlock the lock, release exclusive access to the list
   -----------------------------------------------------
*/
   manager->nunlocks++ ;
   Lock_unlock(manager->lock) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   release a list of Chv objects

   created -- 98may02, cca
   ------------------------------
*/
void
ChvManager_releaseListOfObjects (
   ChvManager   *manager,
   Chv          *head
) {
Chv   *chv1, *chv2, *prev ;
int    nbytes1, nbytes2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || head == NULL ) {
   fprintf(stderr, 
       "\n fatal error in ChvManager_releaseListOfObjects(%p,%p)"
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
}
if ( manager->mode == 0 ) {
/*
   ---------------
   release storage
   ---------------
*/
   while ( (chv1 = head) != NULL ) {
      head = head->next ;
      manager->nbytesactive -= Chv_nbytesInWorkspace(chv1) ;
      manager->nactive-- ;
      manager->nreleases++ ;
      Chv_free(chv1) ;
   }
} else {
/*
   -------------------
   recycle the objects
   -------------------
*/
   while ( head != NULL ) {
      chv1 = head ;
      head = chv1->next ;
/*
      --------------------------------------------------------
      find a place in the list where the Chv objects are 
      sorted in ascending order of the size of their workspace
      --------------------------------------------------------
*/
      nbytes1 = Chv_nbytesInWorkspace(chv1) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n\n trying to release chevron %p with %d bytes",
              chv1, nbytes1) ;
#endif
      for ( chv2 = manager->head, prev = NULL ; 
            chv2 != NULL ; 
            chv2 = chv2->next ) {
         nbytes2 = Chv_nbytesInWorkspace(chv2) ;
#if MYDEBUG > 0
         fprintf(stdout, 
                 "\n list chevron %p with %d bytes", chv2, nbytes2) ;
#endif
         if ( nbytes2 > nbytes1 ) {
            break ;
         }
         prev = chv2 ;
      }
      if ( prev == NULL ) {
         manager->head = chv1 ;
#if MYDEBUG > 0
         fprintf(stdout, "\n manager->head = %p", chv1) ;
#endif
      } else {
         prev->next = chv1 ;
#if MYDEBUG > 0
         fprintf(stdout, "\n %p->next = %p", prev, chv1) ;
#endif
      }
      chv1->next = chv2 ;
#if MYDEBUG > 0
      fprintf(stdout, "\n %p->next = %p", chv1, chv2) ;
#endif
      manager->nbytesactive -= Chv_nbytesInWorkspace(chv1) ;
      manager->nactive-- ;
      manager->nreleases++ ;
   }
}
if ( manager->lock != NULL ) {
/*
   -----------------------------------------------------
   unlock the lock, release exclusive access to the list
   -----------------------------------------------------
*/
   manager->nunlocks++ ;
   Lock_unlock(manager->lock) ;
}
return ; }

/*--------------------------------------------------------------------*/
