/*  init.c  */

#include "../MSMD.h"

#define MYDEBUG 0

#define BE_CAUTIOUS 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   initialization procedure


   created -- 96feb25, cca
   ---------------------------------------------------
*/
void
MSMD_init ( 
   MSMD       *msmd,
   Graph      *g, 
   int        stages[],
   MSMDinfo   *info 
) {
int             ierr, ii, iv, nstage, nvtx, stage ;
int             *vwghts ;
MSMDstageInfo   *stageInfo ;
MSMDvtx         *v ;
/*
   --------------------
   check the input data
   --------------------
*/
if ( msmd == NULL || g == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_init(%p,%p,%p,%p)"
           "\n bad input\n", msmd, g, stages, info) ;
   exit(-1) ;
}
/*
   ---------------------
   clear the data fields 
   ---------------------
*/
MSMD_clearData(msmd) ;
/*
   ----------------------------
   store the number of vertices
   ----------------------------
*/
msmd->nvtx = nvtx = g->nvtx ;
/*
   --------------------------
   allocate the IIheap object
   --------------------------
*/
msmd->heap = IIheap_new() ;
IIheap_init(msmd->heap, nvtx) ;
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n heap initialized") ;
   fflush(info->msgFile) ;
}
info->nbytes += IIheap_sizeOf(msmd->heap) ;
/*
   --------------------------
   allocate the IP structures
   --------------------------
*/
msmd->incrIP = nvtx ;
msmd->baseIP = IP_init(2*nvtx, IP_FORWARD) ;
msmd->freeIP = msmd->baseIP + 1 ;
msmd->baseIP->next = NULL ;
info->nbytes += nvtx*sizeof(struct _IP) ;
/*
   ---------------------
   allocate the vertices
   ---------------------
*/
ALLOCATE(msmd->vertices, struct _MSMDvtx, nvtx) ;
info->nbytes += nvtx*sizeof(struct _MSMDvtx) ;
for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
   v->id       =  iv  ;
   v->mark     =  'O' ;
   v->status   =  'R' ;
   v->bndwght  =   0  ;
   v->par      = NULL ;
   v->subtrees = NULL ;
   Graph_adjAndSize(g, iv, &v->nadj, &v->adj) ;
#if BE_CAUTIOUS
   for ( ii = 0 ; ii < v->nadj ; ii++ ) {
      if ( v->adj[ii] < 0 || v->adj[ii] > nvtx ) {
         fprintf(stderr, "\n bad adj, v = %d :", iv) ;
         IVfp80(stderr, v->nadj, v->adj, 20, &ierr) ;
         exit(-1) ;
      }
   }
#endif
}
if ( (vwghts = g->vwghts) == NULL ) {
   for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
      v->wght = 1 ;
   }
} else {
   for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
      v->wght = vwghts[iv] ;
   }
}
if ( stages == NULL ) {
   for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
      v->stage = 0 ;
   }
} else {
   for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
      v->stage = stages[iv] ;
   }
}
/*
   -------------------------
   allocate the work vectors
   -------------------------
*/
IV_init1(&msmd->ivtmpIV, nvtx) ;
IV_init1(&msmd->reachIV, nvtx) ;
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n vectors initialized") ;
   fprintf(info->msgFile, "\n ivtmpIV = %p", &msmd->ivtmpIV) ;
   IV_writeForHumanEye(&msmd->ivtmpIV, info->msgFile) ;
   fprintf(info->msgFile, "\n reachIV = %p", &msmd->reachIV) ;
   IV_writeForHumanEye(&msmd->reachIV, info->msgFile) ;
   fflush(info->msgFile) ;
}
info->nbytes += 2*nvtx*sizeof(int) ;
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n nvtx = %d, nvtx = %d", nvtx, nvtx) ;
   fflush(info->msgFile) ;
}
/*
   ---------------------------------------------------------------
   set the number of stages and allocate the MSMDinfoStages vector
   ---------------------------------------------------------------
*/
nstage = (stages == NULL) ? 0 : IVmax(nvtx, stages, &iv) ;
info->nstage = nstage ;
ALLOCATE(info->stageInfo, struct _MSMDstageInfo, 3+nstage) ;
for ( stage = 0, stageInfo = info->stageInfo ;
      stage <= 2 + nstage ;
      stage++, stageInfo++ ) {
   stageInfo->nstep    =  0  ;
   stageInfo->nfront   =  0  ;
   stageInfo->welim    =  0  ;
   stageInfo->nfind    =  0  ;
   stageInfo->nzf      =  0  ;
   stageInfo->ops      = 0.0 ;
   stageInfo->nexact2  =  0  ;
   stageInfo->nexact3  =  0  ;
   stageInfo->napprox  =  0  ;
   stageInfo->ncheck   =  0  ;
   stageInfo->nindst   =  0  ;
   stageInfo->noutmtch =  0  ;
}

return ; }

/*--------------------------------------------------------------------*/
