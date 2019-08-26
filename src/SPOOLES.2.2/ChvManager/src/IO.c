/*  IO.c  */

#include "../ChvManager.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form

   created -- 98may02, cca
   ----------------------------------------
*/
void
ChvManager_writeForHumanEye (
   ChvManager   *manager,
   FILE         *fp
) {
Chv   *chv ;
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ChvManager_writeForHumanEye(%p,%p)"
           "\n bad input\n", manager, fp) ;
   exit(-1) ;
}
fprintf(fp, 
        "\n\n ChvManager object at address %p"
        "\n     %d active objects, %d bytes active"
        "\n     %d total bytes requested, %d total bytes allocated " 
        "\n     %d requests, %d releases, %d locks, %d unlocks", 
        manager, manager->nactive, manager->nbytesactive,
        manager->nbytesrequested, manager->nbytesalloc,
        manager->nrequests, manager->nreleases,
        manager->nlocks, manager->nunlocks) ;
/*
for ( chv = manager->head ; chv != NULL ; chv = chv->next ) {
   fprintf(fp, "\n chv %d, nbytes %d",
           chv->id, Chv_nbytesInWorkspace(chv)) ;
}
*/
return ; }

/*--------------------------------------------------------------------*/
