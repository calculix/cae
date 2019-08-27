/*  update.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   accumulate updates to front J, store them in the Chv object

   created -- 98may04, cca
   ------------------------------------------------------------
*/
void
FrontMtx_update (
   FrontMtx   *frontmtx,
   Chv        *frontJ,
   IP         *heads[],
   char       status[],
   DV         *tempDV,
   int        msglvl,
   FILE       *msgFile
) {
SubMtx   *mtxD, *mtxL, *mtxU ;
int      I, J, nfront ;
IP       *first, *ip, *last, *nextip ;
/*
   --------------------------------------
   loop over the fronts that may update J
   --------------------------------------
*/
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n inside FrontMtx_update(%d)", frontJ->id) ;
   fflush(stdout) ;
}
J = frontJ->id ;
nfront = frontmtx->nfront ;
for ( ip = heads[J], heads[J] = first = last = NULL ; 
      ip != NULL ; 
      ip = nextip ) {
   nextip = ip->next ;
   I = ip->val ;
   if ( status == NULL || status[I] == 'F' ) {
      mtxD = FrontMtx_diagMtx(frontmtx, I) ;
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n   update from I %d, mtxD = %p", I, mtxD) ;
         fflush(stdout) ;
      }
      if ( mtxD != NULL ) {
/*
         -------------------------------------------------
         front I did have some rows and columns eliminated
         -------------------------------------------------
*/
         mtxU = FrontMtx_upperMtx(frontmtx, I, nfront) ;
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n   mtxU = %p", mtxU) ;
            fflush(stdout) ;
         }
         if ( mtxU != NULL ) {
/*
            ------------------------------
            the U_{I,bnd{I}} matrix exists
            ------------------------------
*/
            if ( FRONTMTX_IS_SYMMETRIC(frontmtx) ) {
               Chv_updateS(frontJ, mtxD, mtxU, tempDV) ;
            } else if ( FRONTMTX_IS_HERMITIAN(frontmtx) ) {
               Chv_updateH(frontJ, mtxD, mtxU, tempDV) ;
/*
fprintf(msgFile, "\n after update from front %d", mtxD->rowid) ;
Chv_writeForHumanEye(frontJ, msgFile) ;
*/
            } else if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
               mtxL = FrontMtx_lowerMtx(frontmtx, nfront, I) ;
               if ( msglvl > 3 ) {
                  fprintf(msgFile, "\n   mtxL = %p", mtxL) ;
                  fflush(stdout) ;
               }
               if ( mtxL != NULL ) {
/*
                  ------------------------------
                  the L_{bnd{I},I} matrix exists
                  ------------------------------
*/
                  Chv_updateN(frontJ, mtxL, mtxD, mtxU, tempDV) ;
               }
            }
         }
      }
      if ( last == NULL ) {
         last = ip ;
      }
      ip->next = first ;
      first    = ip    ;
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n   update from I %d is finished", I) ;
         fflush(stdout) ;
      }
   } else {
      ip->next = heads[J] ;
      heads[J] = ip ;
   }
}
/*
if ( frontJ->id == frontmtx->nfront - 1 ) {
   fprintf(msgFile, "\n\n last front after updates made") ;
   Chv_writeForHumanEye(frontJ, msgFile) ;
   fflush(msgFile) ;
}
*/
if ( last != NULL ) {
/*
   ------------------------------------
   link the IP objects to the free list
   ------------------------------------
*/
   last->next = heads[nfront] ;
   heads[nfront] = first ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n leaving FrontMtx_update(%d)", frontJ->id) ;
   fflush(stdout) ;
}
return ; }

/*--------------------------------------------------------------------*/
