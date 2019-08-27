/*  semi.c  */

#include "../ETree.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- 

   to find the optimal domain/schur complement partition
   for a semi-implict factorization.

   the gain of a subtree sbt(J) is equal to

   |L_{bnd{J},sbt{J}}| - |A_{bnd{J},sbt{J}}|
      - alpha *|L_{sbt{J},sbt{J}}| 

   when alpha = 0 we minimize active storage
   when alpha = 1 we minimize solve operations

   *ptotalgain is filled with the total gain

   the return value is compidsIV,
      compids[J] = 0 --> J is in the schur complement
      compids[J] != 0 --> J is in domain compids[J]

   created -- 98jun20, cca
   -----------------------------------------------------
*/
IV *
ETree_optPart (
   ETree    *etree,
   Graph    *graph,
   IVL      *symbfacIVL,
   double   alpha, 
   int      *ptotalgain,
   int      msglvl,
   FILE     *msgFile
) {
int    ii, I, J, K, nfront, nvtx, sizeJ, v, vsize, w, wght, wghtJ ;
int    *adjJ, *colCountsA, *colCountsL, *colSbtCountsA, *colSbtCountsL, 
       *diagCountsL, *diagSbtCountsL, *fch, *gain, *nodwghts, 
       *rowCountsA, *rowCountsL, *sib, *vadj, *vtxToFront, *vwghts ;
IV     *compidsIV, *gainIV ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || graph == NULL || symbfacIVL == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in ETree_optPart()"
           "\n bad input\n") ;
   exit(-1) ;
}
nfront = etree->nfront ;
nodwghts = ETree_nodwghts(etree) ;
vtxToFront = ETree_vtxToFront(etree) ;
tree = etree->tree ;
fch  = tree->fch  ;
sib  = tree->sib  ;
if ( (nvtx = graph->nvtx) != etree->nvtx ) {
   fprintf(stderr, "\n fatal error in ETree_optPart()"
           "\n etree->nvtx = %d, graph->nvtx = %d\n",
           etree->nvtx, graph->nvtx) ;
   exit(-1) ;
}
vwghts = graph->vwghts ;
/*
   ---------------------------------------
   compute the # of entries in the 
   offdiagonal front rows and columns of A
   ---------------------------------------
*/
rowCountsA = IVinit(nfront, 0) ;
colCountsA = IVinit(nfront, 0) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   J = vtxToFront[v] ;
   Graph_adjAndSize(graph, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      w = vadj[ii] ;
      if ( (K = vtxToFront[w]) > J ) {
         wght = (vwghts == NULL) ? 1 : vwghts[v]*vwghts[w] ;
         colCountsA[J] += wght ;
         rowCountsA[K] += wght ;
      }
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n rowCountsA") ;
   IVfprintf(msgFile, nfront, rowCountsA) ;
   fprintf(msgFile, "\n\n colCountsA") ;
   IVfprintf(msgFile, nfront, colCountsA) ;
}
/*
   ---------------------------------------------------------
   compute colSbtCountsA[J] = # entries in A_{bnd{J},sbt{J}}
   ---------------------------------------------------------
*/
colSbtCountsA = IVinit(nfront, 0) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   colSbtCountsA[J] = colCountsA[J] - rowCountsA[J] ;
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      colSbtCountsA[J] += colSbtCountsA[I] ;
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n colSbtCountsA") ;
   IVfprintf(msgFile, nfront, colSbtCountsA) ;
}
/*
   ---------------------------------------
   compute the # of entries in the 
   offdiagonal front rows and columns of L
   ---------------------------------------
*/
rowCountsL  = IVinit(nfront, 0) ;
colCountsL  = IVinit(nfront, 0) ;
diagCountsL = IVinit(nfront, 0) ;
for ( J = 0 ; J < nfront ; J++ ) {
   IVL_listAndSize(symbfacIVL, J, &sizeJ, &adjJ) ;
   wghtJ = nodwghts[J] ;
   diagCountsL[J] = (wghtJ*(wghtJ+1))/2 ;
   for ( ii = 0 ; ii < sizeJ ; ii++ ) {
      v = adjJ[ii] ;
      if ( (K = vtxToFront[J]) > J ) {
         wght = (vwghts == NULL) ? 1 : wghtJ*vwghts[v] ;
         colCountsL[J] += wght ;
         rowCountsL[K] += wght ;
      }
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n rowCountsL") ;
   IVfprintf(msgFile, nfront, rowCountsL) ;
   fprintf(msgFile, "\n\n colCountsL") ;
   IVfprintf(msgFile, nfront, colCountsL) ;
   fprintf(msgFile, "\n\n diagCountsL") ;
   IVfprintf(msgFile, nfront, diagCountsL) ;
}
/*
   ---------------------------------------------------------
   compute colSbtCountsL[J]  = # entries in L_{bnd{J},sbt{J}}
   compute diagSbtCountsL[J] = # entries in L_{sbt{J},sbt{J}}
   ---------------------------------------------------------
*/
colSbtCountsL  = IVinit(nfront, 0) ;
diagSbtCountsL = IVinit(nfront, 0) ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   colSbtCountsL[J]  = colCountsL[J] - rowCountsL[J] ;
   diagSbtCountsL[J] = rowCountsL[J] ;
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      colSbtCountsL[J]  +=  colSbtCountsL[I] ;
      diagSbtCountsL[J] += diagSbtCountsL[I] ;
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n colSbtCountsL") ;
   IVfprintf(msgFile, nfront, colSbtCountsL) ;
   fprintf(msgFile, "\n\n diagSbtCountsL") ;
   IVfprintf(msgFile, nfront, diagSbtCountsL) ;
}
/*
   -----------------------
   compute the gain vector
   -----------------------
*/
gainIV = IV_new() ;
IV_init(gainIV, nfront, NULL) ;
gain = IV_entries(gainIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   gain[J] = colSbtCountsL[J] - colSbtCountsA[J] 
             - alpha*diagCountsL[J] ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n gain") ;
   IVfprintf(msgFile, nfront, gain) ;
}
/*
   ----------------------------------------------
   get the component ids of the optimal partition
   ----------------------------------------------
*/
compidsIV = Tree_maximizeGainIV(tree, gainIV, ptotalgain, 
                                msglvl, msgFile) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n total gain = %d", *ptotalgain) ;
   fprintf(msgFile, "\n\n compidsIV") ;
   IV_writeForHumanEye(compidsIV, msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(colCountsA) ;
IVfree(rowCountsA) ;
IVfree(colSbtCountsA) ;
IVfree(colCountsL) ;
IVfree(rowCountsL) ;
IVfree(diagCountsL) ;
IVfree(colSbtCountsL) ;
IVfree(diagSbtCountsL) ;
IV_free(gainIV) ;

return(compidsIV) ; }

/*--------------------------------------------------------------------*/
