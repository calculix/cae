/*  orderViaMS.c  */

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
   ------------------------------------------------------------------
   purpose -- return an ETree object for a multisection ordering

   graph -- graph to order
   maxdomainsize -- used to control the incomplete nested dissection 
     process. any subgraph whose weight is less than maxdomainsize 
     is not split further.
   seed    -- random number seed
   msglvl  -- message level, 0 --> no output, 1 --> timings
   msgFile -- message file

   created -- 97nov06, cca
   ------------------------------------------------------------------
*/
ETree *
orderViaMS (
   Graph   *graph,
   int     maxdomainsize,
   int     seed,
   int     msglvl,
   FILE    *msgFile
) {
double   t1, t2 ;
DSTree   *dstree ;
ETree    *etree ;
int      nvtx, Nvtx ;
IV       *eqmapIV, *stagesIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || maxdomainsize <= 0 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in orderViaMS(%p,%d,%d,%d,%p)"
           "\n bad input\n", 
           graph, maxdomainsize, seed, msglvl, msgFile) ;
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
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n CPU %8.3f : get equivalence map", t2 - t1) ;
   fflush(msgFile) ;
}
Nvtx = 1 + IV_max(eqmapIV) ;
if ( Nvtx <= COMPRESS_FRACTION * nvtx ) {
   MARKTIME(t1) ;
   graph = Graph_compress2(graph, eqmapIV, 1) ;
   MARKTIME(t2) ;
   if ( msglvl > 0 ) {
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
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n CPU %8.3f : sort adjacency", t2 - t1) ;
   fflush(msgFile) ;
}
/*
   -----------------------------
   get the domain separator tree
   -----------------------------
*/
{
GPart       *gpart ;
DDsepInfo   *info ;

info = DDsepInfo_new() ;
info->seed = seed ;
info->maxcompweight = maxdomainsize ; 
gpart = GPart_new() ;
GPart_init(gpart, graph) ;
GPart_setMessageInfo(gpart, msglvl, msgFile) ;
dstree = GPart_RBviaDDsep(gpart, info) ;
DSTree_renumberViaPostOT(dstree) ;
if ( msglvl > 0 ) {
   DDsepInfo_writeCpuTimes(info, msgFile) ;
}
DDsepInfo_free(info) ;
GPart_free(gpart) ;
}
/*
   ---------------------
   get the stages vector
   ---------------------
*/
stagesIV = DSTree_MS2stages(dstree) ;
DSTree_free(dstree) ;
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
MSMD_order(msmd, graph, IV_entries(stagesIV), info) ;
etree = MSMD_frontETree(msmd) ;
if ( msglvl > 0 ) {
   MSMDinfo_print(info, msgFile) ;
}
MSMDinfo_free(info) ;
MSMD_free(msmd) ;
IV_free(stagesIV) ;
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
   if ( msglvl > 0 ) {
      fprintf(msgFile, "\n CPU %8.3f : sort adjacency", t2 - t1) ;
      fflush(msgFile) ;
   }
}
return(etree) ; }

/*--------------------------------------------------------------------*/
