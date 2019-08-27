/*  IO.c  */

#include "../ChvList.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form

   created -- 98may02, cca
   ----------------------------------------
*/
void
ChvList_writeForHumanEye (
   ChvList   *chvlist,
   FILE      *fp
) {
Chv   *chv ;
int   ilist ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvlist == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ChvList_writeForHumanEye(%p,%p)"
           "\n bad input\n", chvlist, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n ChvList object at address %p"
            "\n    %d lists, %d locks",
        chvlist, chvlist->nlist, chvlist->nlocks) ;
if ( chvlist->lock != NULL ) {
   fprintf(fp, "\n    lock located at %p", chvlist->lock) ;
} else {
   fprintf(fp, "\n    lock is NULL") ;
}
for ( ilist = 0 ; ilist < chvlist->nlist ; ilist++ ) {
   fprintf(fp, "\n %6d", ilist) ;
   if ( chvlist->counts != NULL ) {
      fprintf(fp, " %6d", chvlist->counts[ilist]) ;
   } else {
      fprintf(fp, " %6d", 0) ;
   }
   if ( chvlist->flags != NULL ) {
      fprintf(fp, " %6c", chvlist->flags[ilist]) ;
   } else {
      fprintf(fp, " %6c", 'N') ;
   }
   if ( chvlist->heads != NULL && chvlist->heads[ilist] != NULL ) {
      fprintf(fp, " %10p", chvlist->heads[ilist]) ;
   } else {
      fprintf(fp, "      NULL") ;
   }
}
/*
for ( chv = chvlist->head ; chv != NULL ; chv = chv->next ) {
   fprintf(fp, "\n chv %d, nbytes %d",
           chv->id, Chv_nbytesInWorkspace(chv)) ;
}
*/
return ; }

/*--------------------------------------------------------------------*/
