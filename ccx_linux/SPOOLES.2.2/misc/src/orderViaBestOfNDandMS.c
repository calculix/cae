/*  orderViaBestOfNDandMS.c  */

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
   purpose -- return an ETree object for the better of a 
              nested dissection and multisection orderings

   graph -- graph to order
   maxdomainsize -- used to control the incomplete nested dissection 
     process. any subgraph whose weight is less than maxdomainsize 
     is not split further.
   maxzeros -- maximum number of zero entries allowed in a front
   maxsize  -- maximum number of internal columns in a front
   seed     -- random number seed
   msglvl  -- message level
      1 -- timings and statistics
      2 -- more timings and statistics
      3 -- lots of output
   msgFile -- message file

   created -- 98aug26, cca
   ------------------------------------------------------------------
*/
ETree *
orderViaBestOfNDandMS (
   Graph    *graph,
   int      maxdomainsize,
   int      maxzeros,
   int      maxsize,
   int      seed,
   int      msglvl,
   FILE     *msgFile
) {
double   compressCPU, dstreeCPU, eqmapCPU, miscCPU, msCPU, msnzf, 
         msops, ndCPU, ndops, ndnzf, nzf, ops, totalCPU, transformCPU, 
         t0, t1, t2, t3 ;
DSTree   *dstree ;
ETree    *etree, *etree2, *etreeMS, *etreeND ;
int      nfront, msnfront, msnind, ndnfront, ndnind, nind, nvtx, Nvtx ;
IV       *eqmapIV, *stagesIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL ) {
   fprintf(stderr, "\n fatal error in orderViaBestOfNDandMS()"
           "\n graph is NULL\n") ;
   exit(-1) ;
}
if ( maxdomainsize <= 0 ) {
   fprintf(stderr, "\n fatal error in orderViaBestOfNDandMS()"
           "\n maxdomainsize %d\n", maxdomainsize) ;
   exit(-1) ;
}
if ( maxzeros < 0 ) {
   fprintf(stderr, "\n fatal error in orderViaBestOfNDandMS()"
           "\n maxzeros %d\n", maxzeros) ;
   exit(-1) ;
}
if ( maxsize <= 0 ) {
   fprintf(stderr, "\n fatal error in orderViaBestOfNDandMS()"
           "\n maxsize %d\n", maxsize) ;
   exit(-1) ;
}
if ( (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in orderViaBestOfNDandMS()"
           "\n msglvl %d, msgFile %p\n", msglvl, msgFile) ;
   exit(-1) ;
}
MARKTIME(t0) ;
/*
   ------------------------------
   compress the graph if worth it
   ------------------------------
*/
nvtx = graph->nvtx ;
MARKTIME(t1) ;
eqmapIV = Graph_equivMap(graph) ;
MARKTIME(t2) ;
eqmapCPU = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : get equivalence map", t2 - t1) ;
   fflush(msgFile) ;
}
Nvtx = 1 + IV_max(eqmapIV) ;
if ( Nvtx <= COMPRESS_FRACTION * nvtx ) {
   MARKTIME(t1) ;
   graph = Graph_compress2(graph, eqmapIV, 1) ;
   MARKTIME(t2) ;
   compressCPU = t2 - t1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : compress graph", t2 - t1) ;
      fflush(msgFile) ;
   }
} else {
   compressCPU = 0.0 ;
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
   -----------------------------
   get the domain separator tree
   -----------------------------
*/
MARKTIME(t1) ;
{
GPart       *gpart ;
DDsepInfo   *info ;

info = DDsepInfo_new() ;
info->seed = seed ;
info->maxcompweight = maxdomainsize ; 
info->alpha         = 0.1 ; 
gpart = GPart_new() ;
GPart_init(gpart, graph) ;
GPart_setMessageInfo(gpart, msglvl, msgFile) ;
dstree = GPart_RBviaDDsep(gpart, info) ;
DSTree_renumberViaPostOT(dstree) ;
if ( msglvl > 1 ) {
   DDsepInfo_writeCpuTimes(info, msgFile) ;
}
DDsepInfo_free(info) ;
GPart_free(gpart) ;
}
MARKTIME(t2) ;
dstreeCPU = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n CPU %8.3f : construct domain/separator tree", t2 - t1) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   get the stages vector for nested dissection
   -------------------------------------------
*/
MARKTIME(t1) ;
stagesIV = DSTree_NDstages(dstree) ;
MARKTIME(t2) ;
ndCPU = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : get stages for ND", t2 - t1) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   order the vertices and extract the front tree
   ---------------------------------------------
*/
MARKTIME(t2) ;
{
MSMDinfo   *info ;
MSMD       *msmd ;

MARKTIME(t1) ;
info = MSMDinfo_new() ;
info->seed         = seed    ;
info->compressFlag = 2       ;
info->msglvl       = msglvl  ;
info->msgFile      = msgFile ;
msmd = MSMD_new() ;
MSMD_order(msmd, graph, IV_entries(stagesIV), info) ;
etreeND = MSMD_frontETree(msmd) ;
MARKTIME(t2) ;
ndCPU += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : get tree for ND", t2 - t1) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n Nested Dissection information") ;
   MSMDinfo_print(info, msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Nested Dissection tree") ;
   ETree_writeForHumanEye(etreeND, msgFile) ;
}
MARKTIME(t1) ;
MSMDinfo_free(info) ;
MSMD_free(msmd) ;
IV_free(stagesIV) ;
MARKTIME(t2) ;
ndCPU += t2 - t1 ;
}
/*
   --------------------------------------
   get the stages vector for multisection
   --------------------------------------
*/
MARKTIME(t1) ;
stagesIV = DSTree_MS2stages(dstree) ;
MARKTIME(t2) ;
msCPU = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : get stages for MS", t2 - t1) ;
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

MARKTIME(t1) ;
info = MSMDinfo_new() ;
info->seed         = seed    ;
info->compressFlag = 2       ;
info->msglvl       = msglvl  ;
info->msgFile      = msgFile ;
msmd = MSMD_new() ;
MSMD_order(msmd, graph, IV_entries(stagesIV), info) ;
etreeMS = MSMD_frontETree(msmd) ;
MARKTIME(t2) ;
msCPU += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : get tree for ND", t2 - t1) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Multisection information") ;
   MSMDinfo_print(info, msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n Multisection tree") ;
   ETree_writeForHumanEye(etreeMS, msgFile) ;
}
MARKTIME(t1) ;
MSMDinfo_free(info) ;
MSMD_free(msmd) ;
IV_free(stagesIV) ;
MARKTIME(t2) ;
msCPU += t2 - t1 ;
}
/*
   --------------------------------------
   keep the better of the two front trees
   --------------------------------------
*/
ndnfront = ETree_nfront(etreeND) ;
msnfront = ETree_nfront(etreeMS) ;
ndnind   = ETree_nFactorIndices(etreeND) ;
msnind   = ETree_nFactorIndices(etreeMS) ;
ndnzf    = ETree_nFactorEntries(etreeND, SPOOLES_SYMMETRIC) ;
msnzf    = ETree_nFactorEntries(etreeMS, SPOOLES_SYMMETRIC) ;
ndops    = ETree_nFactorOps(etreeND, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
msops    = ETree_nFactorOps(etreeMS, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
if ( ndops <= msops ) {
   etree = etreeND ;
   ETree_free(etreeMS) ;
} else {
   etree = etreeMS ;
   ETree_free(etreeND) ;
}
/*
   ------------------------
   transform the front tree
   ------------------------
*/
MARKTIME(t1) ;
etree = ETree_transform(etree, graph->vwghts, maxzeros, maxsize, seed);
MARKTIME(t2) ;
transformCPU = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : transform tree", t2 - t1) ;
   fflush(msgFile) ;
}
nfront = ETree_nfront(etree) ;
nind   = ETree_nFactorIndices(etree) ;
nzf    = ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC) ;
ops    = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n real symmetric : final ops %.0f",
           ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC)) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   expand the front tree if the graph was compressed
   -------------------------------------------------
*/
if ( eqmapIV != NULL ) {
   etree2 = ETree_expand(etree, eqmapIV) ;
   ETree_free(etree) ;
   etree = etree2 ;
   Graph_free(graph) ;
   IV_free(eqmapIV) ;
} else {
   IVL_sortUp(graph->adjIVL) ;
}
DSTree_free(dstree) ;
MARKTIME(t3) ;
/*
   ------------------------
   print out the statistics
   ------------------------
*/
if ( msglvl >= 1 ) {
   fprintf(msgFile, "\n\n----------------------------------------"
           "\n\n Order the graph via best of ND and MS") ;
   fprintf(msgFile,
"\n\n                    # fronts  # indices    # entries         # ops"
"\n nested dissection   %7d %10d %12.0f  %12.0f"
"\n multisection        %7d %10d %12.0f  %12.0f"
"\n final ordering      %7d %10d %12.0f  %12.0f",
           ndnfront, ndnind, ndnzf, ndops,
           msnfront, msnind, msnzf, msops,
           nfront, nind, nzf, ops) ;
   totalCPU = t3 - t0 ;
   miscCPU  = totalCPU - (eqmapCPU + compressCPU + dstreeCPU 
                                   + ndCPU + msCPU + transformCPU) ;

   if ( totalCPU > 0.0 ) {
      fprintf(msgFile, 
              "\n\n CPU breakdown                            CPU %%"
              "\n    make equivalence map             %8.3f %6.2f"
              "\n    compress graph                   %8.3f %6.2f"
              "\n    construct domain/separator tree  %8.3f %6.2f"
              "\n    evaluate nested dissection       %8.3f %6.2f"
              "\n    evaluate multisection            %8.3f %6.2f"
              "\n    transform final tree             %8.3f %6.2f"
              "\n    miscellaneous time               %8.3f %6.2f"
              "\n    total time                       %8.3f",
       eqmapCPU, 100.*eqmapCPU/totalCPU,
       compressCPU, 100.*compressCPU/totalCPU,
       dstreeCPU, 100.*dstreeCPU/totalCPU,
       ndCPU, 100.*ndCPU/totalCPU,
       msCPU, 100.*msCPU/totalCPU,
       transformCPU, 100.*transformCPU/totalCPU,
       miscCPU, 100.*miscCPU/totalCPU,
       totalCPU) ;
   }
   fprintf(msgFile, "\n\n----------------------------------------") ;
}
return(etree) ; }

/*--------------------------------------------------------------------*/
