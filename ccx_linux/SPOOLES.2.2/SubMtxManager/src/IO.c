/*  IO.c  */

#include "../SubMtxManager.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form

   created -- 98may02, cca
   ----------------------------------------
*/
void
SubMtxManager_writeForHumanEye (
   SubMtxManager   *manager,
   FILE            *fp
) {
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxManager_writeForHumanEye(%p,%p)"
           "\n bad input\n", manager, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n\n SubMtxManager object at address %p"
        "\n     %d active objects, %d bytes active"
        "\n     %d total bytes requested, %d total bytes allocated " 
        "\n     %d requests, %d releases, %d locks, %d unlocks", 
        manager, manager->nactive, manager->nbytesactive,
        manager->nbytesrequested, manager->nbytesalloc,
        manager->nrequests, manager->nreleases,
        manager->nlocks, manager->nunlocks) ;
/*
for ( mtx = manager->head ; mtx != NULL ; mtx = mtx->next ) {
   fprintf(fp, "\n mtx (%d,%d), nbytes %d",
           mtx->rowid, mtx->colid, SubMtx_nbytesInWorkspace(mtx)) ;
}
*/
return ; }

/*--------------------------------------------------------------------*/
