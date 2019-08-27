/*  permute.c  */

#include "../ETree.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   fill the new-to-old permutation vector for the fronts

   created -- 96jun23, cca
   -----------------------------------------------------
*/
IV *
ETree_newToOldFrontPerm (
   ETree  *etree
) {
int   nfront, nvtx ;
IV    *newToOldIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_newToOldFrontPerm(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
newToOldIV = IV_new() ;
IV_init(newToOldIV, nfront, NULL) ;
/*
   ------------------------------------------------------
   fill the permutation vector by calling the Tree method
   ------------------------------------------------------
*/
Tree_fillNewToOldPerm(etree->tree, IV_entries(newToOldIV)) ;

return(newToOldIV) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   fill the old-to-new permutation vector for the fronts

   created -- 96jun23, cca
   -----------------------------------------------------
*/
IV *
ETree_oldToNewFrontPerm (
   ETree  *etree
) {
int   nfront, nvtx ;
IV    *oldToNewIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_oldToNewFrontPerm(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
oldToNewIV = IV_new() ;
IV_init(oldToNewIV, nfront, NULL) ;
/*
   ------------------------------------------------------
   fill the permutation vector by calling the Tree method
   ------------------------------------------------------
*/
Tree_fillOldToNewPerm(etree->tree, IV_entries(oldToNewIV)) ;

return(oldToNewIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   fill the new-to-old permutation vector for the vertices

   created -- 96oct05, cca
   -------------------------------------------------------
*/
IV *
ETree_newToOldVtxPerm (
   ETree  *etree
) {
int   count, front, nfront, nvtx, v ;
int   *head, *link, *newToOld, *vtxToFront ;
IV    *newToOldVtxIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_newToOldVtxPerm(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   ----------------------------------------------------
   allocate the old-to-new vertex permutation IV object
   ----------------------------------------------------
*/
newToOldVtxIV = IV_new() ;
IV_init(newToOldVtxIV, nvtx, NULL) ;
newToOld = IV_entries(newToOldVtxIV) ;
/*
   -------------------------------------------------------
   get the head/link structure for the vertices and fronts
   -------------------------------------------------------
*/
head = IVinit(nfront, -1) ;
link = IVinit(nvtx,   -1) ;
for ( v = nvtx - 1 ; v >= 0 ; v-- ) {
   front = vtxToFront[v] ;
   link[v] = head[front] ; 
   head[front] = v ; 
}
/*
   ----------------------------------------------
   loop over the fronts in a post-order traversal
   ----------------------------------------------
*/
count = 0 ;
for ( front = Tree_postOTfirst(etree->tree) ;
      front != -1 ;
      front = Tree_postOTnext(etree->tree, front) ) {
   for ( v = head[front] ; v != -1 ; v = link[v] ) {
      newToOld[count++] = v ;
   }
}
IVfree(head) ;
IVfree(link) ;

return(newToOldVtxIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   fill the old-to-new permutation vector for the vertices

   created -- 96jun23, cca
   -------------------------------------------------------
*/
IV *
ETree_oldToNewVtxPerm (
   ETree  *etree
) {
int   count, front, nfront, nvtx, v ;
int   *head, *link, *oldToNew, *vtxToFront ;
IV    *oldToNewVtxIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_oldToNewVtxPerm(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   ----------------------------------------------------
   allocate the old-to-new vertex permutation IV object
   ----------------------------------------------------
*/
oldToNewVtxIV = IV_new() ;
IV_init(oldToNewVtxIV, nvtx, NULL) ;
oldToNew = IV_entries(oldToNewVtxIV) ;
/*
   -------------------------------------------------------
   get the head/link structure for the vertices and fronts
   -------------------------------------------------------
*/
head = IVinit(nfront, -1) ;
link = IVinit(nvtx,   -1) ;
for ( v = nvtx - 1 ; v >= 0 ; v-- ) {
   front = vtxToFront[v] ;
   link[v] = head[front] ; 
   head[front] = v ; 
}
/*
   ----------------------------------------------
   loop over the fronts in a post-order traversal
   ----------------------------------------------
*/
count = 0 ;
for ( front = Tree_postOTfirst(etree->tree) ;
      front != -1 ;
      front = Tree_postOTnext(etree->tree, front) ) {
   for ( v = head[front] ; v != -1 ; v = link[v] ) {
      oldToNew[v] = count++ ;
   }
}
IVfree(head) ;
IVfree(link) ;

return(oldToNewVtxIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- permute the vertices, 
              overwrite entries in the vertex-to-front map

   created -- 96oct03, cca
   -------------------------------------------------------
*/
void
ETree_permuteVertices (
   ETree   *etree,
   IV      *vtxOldToNewIV
) {
int   nvtx, v ;
int   *oldToNew, *temp, *vtxToFront ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || vtxOldToNewIV == NULL
   || (nvtx = etree->nvtx) <= 0 || nvtx != IV_size(vtxOldToNewIV) ) {
   fprintf(stderr, "\n fatal error in ETree_permuteVertices(%p,%p)"
           "\n bad input\n", etree, vtxOldToNewIV) ;
   if ( etree != NULL && vtxOldToNewIV != NULL ) {
      fprintf(stderr, 
              "\n etree->nvtx = %d, IV_size(vtxOldToNewIV) = %d\n",
              etree->nvtx, IV_size(vtxOldToNewIV)) ;
   }
   exit(-1) ;
}
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
oldToNew   = IV_entries(vtxOldToNewIV) ;
temp = IVinit(nvtx, -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   temp[oldToNew[v]] = vtxToFront[v] ;
}
IVcopy(nvtx, vtxToFront, temp) ;
IVfree(temp) ;

return ; }
   
/*--------------------------------------------------------------------*/
