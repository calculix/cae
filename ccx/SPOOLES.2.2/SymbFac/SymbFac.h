/*  SymbFac.h  */

#include "../ETree.h"
#include "../Graph.h"
#include "../InpMtx.h"
#include "../Pencil.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- using one graph, create and return an IVL object 
              that contains a symbolic factorization.

   created -- 97aug29, cca
   -----------------------------------------------------------
*/
IVL *
SymbFac_initFromGraph (
   ETree   *etree,
   Graph   *graph
) ;
/*
   -------------------------------------------------------
   purpose -- create and return an IVL object that 
              contains a symbolic factorization.
   note: we assume that both the ETree and DInpMtx objects
         are in their final permutation

   created -- 97mar15, cca
   -------------------------------------------------------
*/
IVL *
SymbFac_initFromInpMtx (
   ETree    *etree,
   InpMtx   *inpmtx
) ;
/*
   -----------------------------------------------------
   compute the symbolic factorization of a matrix pencil
   note: we assume that both the ETree and Pencil 
         objects are in their final permutation
 
   created -- 98may04, cca
   -----------------------------------------------------
*/
IVL *
SymbFac_initFromPencil (
   ETree    *etree,
   Pencil   *pencil
) ;
/*--------------------------------------------------------------------*/
