/*  split.c  */

#include "../GPart.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   split the graph partition object into pieces

   created -- 95nov29, cca
   --------------------------------------------
*/
void
GPart_split (
   GPart   *gpart
) {
FILE    *msgFile ;
GPart   *gpartchild ;
Graph   *g, *gchild ;
int     domwght, icomp, ierr, msglvl, ncomp, nvtot, nvtx, sepwght ;
int     *compids, *cweights, *map ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || (g = gpart->g) == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_split(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
if ( gpart->fch != NULL ) {
   fprintf(stderr, "\n fatal error in GPart_split(%p)"
           "\n child(ren) exist, already split\n", gpart) ;
   exit(-1) ;
}
msgFile = gpart->msgFile ;
msglvl  = gpart->msglvl  ;
/*
   ------------------------------
   count the number of subgraphs
   and fill the cweights[] vector
   ------------------------------
*/
nvtx = g->nvtx ;
GPart_setCweights(gpart) ;
ncomp    = gpart->ncomp ;
cweights = IV_entries(&gpart->cweightsIV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n\n inside GPart_split, %d components, cweights : ", 
           ncomp) ;
   IV_fp80(&gpart->cweightsIV, msgFile, 25, &ierr) ;
}
if ( ncomp == 1 ) {
   return ;
}
/*
   -----------------------------------------
   compute the weight of the components and
   count the number of nontrivial components
   -----------------------------------------
*/
sepwght = cweights[0] ;
domwght = 0 ;
for ( icomp = 1 ; icomp <= ncomp ; icomp++ ) {
   domwght += cweights[icomp] ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n separator weight = %d, weight of components = %d",
           sepwght, domwght) ;
}
/*
   ------------------------------------------------------
   for each component
      create its subgraph with boundary
      create a GPart object to contain the subgraph
         and set as the child of the present GPart object
   end for
   ------------------------------------------------------
*/
compids = IV_entries(&gpart->compidsIV) ;
for ( icomp = 1 ; icomp <= ncomp ; icomp++ ) {
   gpartchild = GPart_new() ;
   gchild = Graph_subGraph(g, icomp, compids, &map) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n component %d", icomp) ;
      fprintf(msgFile, "\n map to parent") ;
      IVfp80(msgFile, gchild->nvtx + gchild->nvbnd, map, 80, &ierr) ;
      Graph_writeForHumanEye(gchild, msgFile) ;
      fflush(msgFile) ;
   }
   GPart_init(gpartchild, gchild) ;
   nvtot = gpartchild->nvtx + gpartchild->nvbnd ;
   IV_init2(&gpartchild->vtxMapIV, nvtot, nvtot, 1, map) ;
   gpartchild->par = gpart      ;
   gpartchild->sib = gpart->fch ;
   gpart->fch      = gpartchild ;
   gpartchild->msglvl  = gpart->msglvl  ;
   gpartchild->msgFile = gpart->msgFile ;
}

return ; }

/*--------------------------------------------------------------------*/
