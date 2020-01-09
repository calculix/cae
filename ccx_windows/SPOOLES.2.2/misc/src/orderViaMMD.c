/*  orderViaMMD.c  */

#include "../misc.h"
#include "../../timings.h"
/*
   ---------------------------------------------------------------------
   COMPRESS_FRACTION --- 
     if # coarse graph vertices < COMPRESS_FRACTION * # of vertices then
        use the compressed graph
     endif
   ---------------------------------------------------------------------
*/
#define COMPRESS_FRACTION 0.75

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- return an ETree object for 
              a multiple minimum degree ordering

   graph -- graph to order
   seed    -- random number seed
   msglvl  -- message level, 0 --> no output, 1 --> timings
   msgFile -- message file

   created -- 97nov08, cca
   --------------------------------------------------------
*/
ETree *
orderViaMMD (
   Graph   *graph,
   int     seed,
   int     msglvl,
   FILE    *msgFile
) {
double   t1, t2 ;
ETree    *etree ;
int      nvtx, Nvtx ;
IV       *eqmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in orderViaMMD(%p,%d,%d,%p)"
           "\n bad input\n", 
           graph, seed, msglvl, msgFile) ;
   exit(-1) ;
}
/*
   ------------------------------
   compress the graph if worth it
   ------------------------------
*/
nvtx = graph->nvtx ;
MARKTIME(t1) ;
eqmapIV = Graph_equivMap(graph) ;
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : get equivalence map", t2 - t1) ;
   fflush(msgFile) ;
}
Nvtx = 1 + IV_max(eqmapIV) ;
if ( Nvtx <= COMPRESS_FRACTION * nvtx ) {
   MARKTIME(t1) ;
   graph = Graph_compress2(graph, eqmapIV, 1) ;
   MARKTIME(t2) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : compress graph", t2 - t1) ;
      fflush(msgFile) ;
   }
} else {
   IV_free(eqmapIV) ;
   eqmapIV = NULL ;
}
MARKTIME(t1) ;
IVL_sortUp(graph->adjIVL) ;
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : sort adjacency", t2 - t1) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   order the vertices and extract the front tree
   ---------------------------------------------
*/
{
MSMDinfo   *info ;
MSMD       *msmd ;

info = MSMDinfo_new() ;
info->seed         = seed    ;
info->compressFlag = 2       ;
info->msglvl       = msglvl  ;
info->msgFile      = msgFile ;
msmd = MSMD_new() ;
MSMD_order(msmd, graph, NULL, info) ;
etree = MSMD_frontETree(msmd) ;
if ( msglvl > 1 ) {
   MSMDinfo_print(info, msgFile) ;
}
MSMDinfo_free(info) ;
MSMD_free(msmd) ;
}
/*
   -------------------------------------------------
   expand the front tree if the graph was compressed
   -------------------------------------------------
*/
if ( eqmapIV != NULL ) {
   ETree *etree2 = ETree_expand(etree, eqmapIV) ;
   ETree_free(etree) ;
   etree = etree2 ;
   Graph_free(graph) ;
   IV_free(eqmapIV) ;
} else {
   MARKTIME(t1) ;
   IVL_sortUp(graph->adjIVL) ;
   MARKTIME(t2) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : sort adjacency", t2 - t1) ;
      fflush(msgFile) ;
   }
}
return(etree) ; }

/*--------------------------------------------------------------------*/
