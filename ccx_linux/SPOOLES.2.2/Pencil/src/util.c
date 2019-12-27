/*  util.c  */

#include "../Pencil.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   sort and compress the pencil's entries

   created -- 98may02, cca
   --------------------------------------
*/
void
Pencil_sortAndCompress (
   Pencil   *pencil 
) {

if ( pencil->inpmtxA != NULL ) {
   InpMtx_sortAndCompress(pencil->inpmtxA) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_sortAndCompress(pencil->inpmtxB) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   convert the storage to vectors

   created -- 98may02, cca
   ------------------------------
*/
void
Pencil_convertToVectors (
   Pencil   *pencil 
) {

if ( pencil->inpmtxA != NULL ) {
   InpMtx_convertToVectors(pencil->inpmtxA) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_convertToVectors(pencil->inpmtxB) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   map entries to the lower triangle,
   used after a permutation of a symmetric matrix

   created -- 98may02, cca
   ----------------------------------------------
*/
void
Pencil_mapToLowerTriangle (
   Pencil   *pencil 
) {

if ( pencil->inpmtxA != NULL ) {
   InpMtx_mapToLowerTriangle(pencil->inpmtxA) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_mapToLowerTriangle(pencil->inpmtxB) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   map entries to the upper triangle,
   used after a permutation of a symmetric matrix

   created -- 98may02, cca
   ----------------------------------------------
*/
void
Pencil_mapToUpperTriangle (
   Pencil   *pencil 
) {

if ( pencil->inpmtxA != NULL ) {
   InpMtx_mapToUpperTriangle(pencil->inpmtxA) ;
}
if ( pencil->inpmtxB != NULL ) {
   InpMtx_mapToUpperTriangle(pencil->inpmtxB) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to return the full, symmetric adjacency IVL object
              for the graph of (A + B) + sigma * (A + B)^T
 
   created -- 98may02, cca
   -------------------------------------------------------------
*/
IVL *
Pencil_fullAdjacency (
   Pencil  *pencil
) {
IVL   *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( pencil == NULL ) {
   fprintf(stderr, "\n fatal error in Pencil_fullAdjacency(%p)"
           "\n NULL input\n\n", pencil) ;
   exit(-1) ;
}
if ( pencil->sigma[0] == 0.0 && pencil->sigma[1] == 0.0 ) {
   if ( pencil->inpmtxA == NULL ) {
      fprintf(stderr, "\n fatal error in Pencil_fullAdjacency(%p)"
              "\n pencil is identity \n\n", pencil) ;
      exit(-1) ;
   } else {
      adjIVL = InpMtx_fullAdjacency(pencil->inpmtxA) ;
   }
} else {
   if ( pencil->inpmtxB == NULL ) {
      adjIVL = InpMtx_fullAdjacency(pencil->inpmtxA) ;
   } else if ( pencil->inpmtxA == NULL ) {
      adjIVL = InpMtx_fullAdjacency(pencil->inpmtxB) ;
   } else {
      adjIVL = InpMtx_fullAdjacency2(pencil->inpmtxA, pencil->inpmtxB);
   }
}
return(adjIVL) ; }

/*--------------------------------------------------------------------*/
