/*  IO.c  */

#include "../SubMtxList.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form

   created -- 98may02, cca
   ----------------------------------------
*/
void
SubMtxList_writeForHumanEye (
   SubMtxList   *list,
   FILE         *fp
) {
SubMtx   *mtx ;
int    ilist ;
/*
   ---------------
   check the input
   ---------------
*/
if ( list == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtxList_writeForHumanEye(%p,%p)"
           "\n bad input\n", list, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n SubMtxList object at address %p"
            "\n %d lists, %d locks" 
            "\n heads %p, counts %p, flags %p",
            list, list->nlist, list->nlocks,
            list->heads, list->counts, list->flags) ;
for ( ilist = 0 ; ilist < list->nlist ; ilist++ ) {
   fprintf(fp, "\n list %d : ", ilist) ;
   if ( list->heads == NULL ) {
      fprintf(fp, " head NULL") ;
   } else {
      fprintf(fp, " head %p", list->heads[ilist]) ;
   }
   if ( list->counts == NULL ) {
      fprintf(fp, ", counts NULL") ;
   } else {
      fprintf(fp, ", counts %d", list->counts[ilist]) ;
   }
   if ( list->flags == NULL ) {
      fprintf(fp, ", flags NULL") ;
   } else {
      fprintf(fp, ", flags %c", list->flags[ilist]) ;
   }
}
for ( ilist = 0 ; ilist < list->nlist ; ilist++ ) {
   if ( (mtx = list->heads[ilist]) != NULL ) {
      fprintf(fp, "\n list %d :", ilist) ;
      while ( mtx != NULL ) {
         fprintf(fp, "\n    mtx (%d,%d), nbytes %d",
                 mtx->rowid, mtx->colid, 
                 SubMtx_nbytesInWorkspace(mtx)) ;
         mtx = mtx->next ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
