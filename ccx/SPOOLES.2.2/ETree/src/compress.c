/*  compress.c  */

#include "../ETree.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- 
   to create and return an IV object that contains the map 
   from old to new fronts that are fundamental chains.

   created  -- 96jun23, cca
   -------------------------------------------------------
*/
IV *
ETree_fundChainMap (
   ETree   *etree
) {
int   nfront, nvtx ;
IV    *frontmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->tree == NULL 
   || (nfront = etree->nfront) <= 0 || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_fundChainMap(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
/*
   -------------------------------------
   call the Tree object's method to get 
   the map from old fronts to new fronts
   -------------------------------------
*/
frontmapIV = Tree_fundChainMap(etree->tree) ;

return(frontmapIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- 
   to create and return an IV object that contains the map 
   from old to new fronts that are fundamental supernodes

   created  -- 96jun23, cca
   -------------------------------------------------------
*/
IV *
ETree_fundSupernodeMap (
   ETree   *etree
) {
int   child, front, nfront, nfs, nvtx ;
int   *bndwghts, *fch, *frontmap, *nodwghts, *par, *sib ;
IV    *frontmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->tree == NULL 
   || (nfront = etree->nfront) <= 0 || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_fundSupernodeMap(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
par      = etree->tree->par ;
fch      = etree->tree->fch ;
sib      = etree->tree->sib ;
nodwghts = IV_entries(etree->nodwghtsIV) ;
bndwghts = IV_entries(etree->bndwghtsIV) ;
/*
   ------------------------------------------
   fill the map from old fronts to new fronts
   ------------------------------------------
*/
frontmapIV = IV_new() ;
IV_init(frontmapIV, nfront, NULL) ;
frontmap = IV_entries(frontmapIV) ;
nfs = 0 ;
front = etree->tree->root ;
while ( front != -1 ) {
   while ( fch[front] != -1 ) {
      front = fch[front] ;
   }
   frontmap[front] = nfs++ ;
   while ( sib[front] == -1 && par[front] != -1 ) {
      front = par[front] ;
      child = fch[front] ;
      if (   sib[child] != -1 
          || (nodwghts[front] + bndwghts[front] != bndwghts[child]) ) {
         frontmap[front] = nfs++ ;
      } else {
         frontmap[front] = frontmap[child] ;
      }
   }
   front = sib[front] ;
}

return(frontmapIV) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   compress an ETree object given a map from old to new nodes.
   note, a new node must be a connected set of the old nodes.

   return value -- pointer to new ETree object

   created -- 96jun23, cca.
   -----------------------------------------------------------
*/
ETree *
ETree_compress (
   ETree   *etree,
   IV      *frontmapIV
) {
ETree   *etree2 ;
int     nfront, newfront, newnfront, newparfront, nvtx, oldfront,
        parfront, v ;
int     *bndwghts, *frontmap, *newbndwghts, *newnodwghts, *newpar,
        *newvtxToFront, *nodwghts, *par, *vtxToFront ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 || frontmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_compress(%p,%p)"
           "\n bad input\n", etree, frontmapIV) ;
   exit(-1) ;
}
/*
   --------------------------------
   pull out pointers and dimensions
   --------------------------------
*/
nfront     = etree->nfront ;
nvtx       = etree->nvtx   ;
par        = etree->tree->par ;
nodwghts   = IV_entries(etree->nodwghtsIV) ;
bndwghts   = IV_entries(etree->bndwghtsIV) ;
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
newnfront  = 1 + IV_max(frontmapIV) ;
frontmap   = IV_entries(frontmapIV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n newnfront = %d", newnfront) ;
#endif
/*
   -------------------------------
   initialize the new ETree object
   -------------------------------
*/
etree2 = ETree_new() ;
ETree_init1(etree2, newnfront, nvtx) ;
newpar        = etree2->tree->par ;
newnodwghts   = IV_entries(etree2->nodwghtsIV) ;
newbndwghts   = IV_entries(etree2->bndwghtsIV) ;
newvtxToFront = IV_entries(etree2->vtxToFrontIV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n newnodwghts") ;
IV_writeForHumanEye(etree2->nodwghtsIV, stdout) ;
#endif
/*
   ------------------------
   fill the new tree fields
   ------------------------
*/
for ( oldfront = 0 ; oldfront < nfront ; oldfront++ ) {
   newfront = frontmap[oldfront] ;
   parfront = par[oldfront] ;
#if MYDEBUG > 0
   fprintf(stdout, 
        "\n oldfront = %d, nodwght = %d, parfront = %d, newfront = %d",
        oldfront, nodwghts[oldfront], parfront, newfront) ;
   fflush(stdout) ;
#endif
   newnodwghts[newfront] += nodwghts[oldfront] ;
#if MYDEBUG > 0
   fprintf(stdout, 
        "\n nodwghts[%d] = %d, newnodwghts[%d] = %d",
        oldfront, nodwghts[oldfront],
        newfront, newnodwghts[newfront]) ;
   fflush(stdout) ;
#endif
   if (  parfront != -1
      && (newparfront = frontmap[parfront]) != newfront ) {
      newpar[newfront]      = newparfront        ;
      newbndwghts[newfront] = bndwghts[oldfront] ;
#if MYDEBUG > 0
      fprintf(stdout, "\n newparfront = %d", newparfront) ;
      fprintf(stdout, 
              "\n setting newpar[%d] = %d, newbndwghts[%d] = %d",
              newfront, newpar[newfront],
              newfront, newbndwghts[newfront]) ;
      fflush(stdout) ;
#endif
   }
}
Tree_setFchSibRoot(etree2->tree) ;
/*
   ---------------------------------------
   set the map from vertices to new fronts
   ---------------------------------------
*/
for ( v = 0 ; v < nvtx ; v++ ) {
   newvtxToFront[v] = frontmap[vtxToFront[v]] ;
}

return(etree2) ; }

/*--------------------------------------------------------------------*/
