/*  permute.c  */

#include "../Pencil.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------
   permute the matrix pencil
 
   created -- 98may02, cca
   -------------------------
*/
void
Pencil_permute (
   Pencil   *pencil,
   IV       *rowOldToNewIV,
   IV       *colOldToNewIV
) {
if ( pencil->inpmtxA != NULL ) {
   InpMtx_permute(pencil->inpmtxA,
                  IV_entries(rowOldToNewIV),
                  IV_entries(colOldToNewIV)) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_permute(pencil->inpmtxB,
                  IV_entries(rowOldToNewIV),
                  IV_entries(colOldToNewIV)) ;
}
return ; }
 
/*--------------------------------------------------------------------*/
