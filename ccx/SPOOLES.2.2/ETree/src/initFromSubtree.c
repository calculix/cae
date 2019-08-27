/*  initFromSubtree.c  */
 
#include "../ETree.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to initialize subtree with the subtree 
              of the front tree using nodes in nodeidsIV.
              vtxIV is filled with the vertices in the subtree
 
   return values ---
      1 -- normal return
     -1 -- subtree is NULL
     -2 -- nodeidsIV is NULL
     -3 -- etree is NULL
     -4 -- nodeidsIV is invalid
     -5 -- vtxIV is NULL
 
   created -- 98oct15, cca
   -----------------------------------------------------------
*/
int
ETree_initFromSubtree (
   ETree   *subtree,
   IV      *nodeidsIV,
   ETree   *etree,
   IV      *vtxIV
) {
int   J, Jsub, nfrontInETree, nfrontInSubtree, 
      nvtxInETree, nvtxInSubtree, v, vSub ;
int   *bndwghts, *bndwghtsSub, *localmap, *nodwghts, *nodwghtsSub,
      *subtreeNodes, *vtxInSubtree, *vtxToFront, *vtxToFrontSub ;
/*
   ---------------
   check the input
   ---------------
*/
if ( subtree == NULL ) {
   fprintf(stderr, "\n\n error in ETree_initFromSubtree()"
           "\n subtree is NULL\n") ;
   return(-1) ;
}
if ( nodeidsIV == NULL ) {
   fprintf(stderr, "\n\n error in ETree_initFromSubtree()"
           "\n nodeidsIV is NULL\n") ;
   return(-2) ;
}
if ( etree == NULL ) {
   fprintf(stderr, "\n\n error in ETree_initFromSubtree()"
           "\n etree is NULL\n") ;
   return(-3) ;
}
nfrontInETree = ETree_nfront(etree) ;
IV_sizeAndEntries(nodeidsIV, &nfrontInSubtree, &subtreeNodes) ;
if ( nfrontInSubtree < 0 || nfrontInSubtree >= nfrontInETree ) {
   fprintf(stderr, "\n\n error in ETree_initFromSubtree()"
           "\n nfrontInETree = %d, nfrontInSubtree = %d\n",
           nfrontInETree, nfrontInSubtree) ;
   return(-4) ;
}
for ( Jsub = 0 ; Jsub < nfrontInSubtree ; Jsub++ ) {
   J = subtreeNodes[Jsub] ;
   if ( J < 0 || J >= nfrontInETree ) {
      fprintf(stderr, "\n\n error in ETree_initFromSubtree()"
              "\n nfrontInETree = %d, subtreeNodes[%d] = %d\n",
              nfrontInETree, Jsub, subtreeNodes[Jsub]) ;
      return(-4) ;
   }
}
if ( vtxIV == NULL ) {
   fprintf(stderr, "\n\n error in ETree_initFromSubtree()"
           "\n vtxIV is NULL\n") ;
   return(-5) ;
}
nvtxInETree = ETree_nvtx(etree) ;
vtxToFront  = ETree_vtxToFront(etree) ;
/*
   ----------------------------
   create a global-to-local map
   ----------------------------
*/
localmap = IVinit(nfrontInETree, -1) ;
for ( Jsub = 0 ; Jsub < nfrontInSubtree ; Jsub++ ) {
   J = subtreeNodes[Jsub] ;
   localmap[J] = Jsub ;
}
/*
   ---------------------------------------------
   compute the number of vertices in the subtree
   ---------------------------------------------
*/
nvtxInSubtree = 0 ;
for ( v = 0 ; v < nvtxInETree ; v++ ) {
   J = vtxToFront[v] ;
   if ( (Jsub = localmap[J]) != -1 ) {
      nvtxInSubtree++ ;
   }
}
/*
   ----------------------
   initialize the subtree
   ----------------------
*/
ETree_init1(subtree, nfrontInSubtree, nvtxInSubtree) ;
/*
   -----------------------------
   initialize the subtree's tree
   -----------------------------
*/
Tree_initFromSubtree(subtree->tree, nodeidsIV, etree->tree) ;
/*
   -----------------------------------
   set the nodwght and bndwght vectors
   -----------------------------------
*/
nodwghts    = ETree_nodwghts(etree) ;
bndwghts    = ETree_bndwghts(etree) ;
nodwghtsSub = ETree_nodwghts(subtree) ;
bndwghtsSub = ETree_bndwghts(subtree) ;
for ( Jsub = 0 ; Jsub < nfrontInSubtree ; Jsub++ ) {
   J = subtreeNodes[Jsub] ;
   nodwghtsSub[Jsub] = nodwghts[J] ;
   bndwghtsSub[Jsub] = bndwghts[J] ;
}
/*
   -------------------------------------
   set the subtree's vtxToFront[] vector
   and fill vtxIV with the vertices
   -------------------------------------
*/
IV_init(vtxIV, nvtxInSubtree, NULL) ;
vtxInSubtree  = IV_entries(vtxIV) ;
vtxToFrontSub = ETree_vtxToFront(subtree) ;
for ( v = vSub = 0 ; v < nvtxInETree ; v++ ) {
   J = vtxToFront[v] ;
   if ( (Jsub = localmap[J]) != -1 ) {
      vtxInSubtree[vSub] = v ;
      vtxToFrontSub[vSub] = Jsub ;
      vSub++ ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(localmap) ;

return(1) ; }

/*--------------------------------------------------------------------*/
