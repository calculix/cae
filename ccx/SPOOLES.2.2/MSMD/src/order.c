/*  order.c  */

#include "../MSMD.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- to order the graph using multi-stage minimum degree

   g      -- Graph object
   stages -- stage vector for vertices, 
      if NULL then
         all vertices on stage zero.
      otherwise 
         vertices with stage istage are eliminated 
         before any vertices with stage > istage

   working storage is free'd,
   statistics can be accessed through their variables or printed
   via the void MSMD_printStats(MSMD*,FILE*) method.

   created -- 96feb25, cca
   ---------------------------------------------------------------------
*/
void
MSMD_order ( 
   MSMD       *msmd,
   Graph      *g, 
   int        stages[],
   MSMDinfo   *info
) {
double          t0, t1, t2, t3 ;
int             istage, iv, nstage, nvtx ;
IP              *ip ;
MSMDstageInfo   *now, *total ;
MSMDvtx         *v ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if (  msmd == NULL || g == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_order(%p,%p,%p,%p)"
           "\n bad input\n", msmd, g, stages, info) ;
   exit(-1) ;
}
if ( info->msglvl > 2 ) {
   fprintf(info->msgFile, "\n\n inside MSMD_order()") ;
   if ( stages != NULL ) {
      int ierr ;
      fprintf(info->msgFile, "\n stages[%d]", g->nvtx) ;
      IVfp80(info->msgFile, g->nvtx, stages, 80, &ierr) ;
   }
   fflush(info->msgFile) ;
}
/*
   -------------------------------
   check the information structure
   -------------------------------
*/
if ( MSMDinfo_isValid(info) != 1 ) {
   fprintf(stderr, "\n fatal error in MSMD_order(%p,%p,%p,%p)"
           "\n bad MSMDinfo object\n", msmd, g, stages, info) ;
   exit(-1) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n\n trying to initialize MSMD object ") ;
   Graph_writeForHumanEye(g, info->msgFile) ;
   fflush(info->msgFile) ;
}
MSMD_init(msmd, g, stages, info) ;
nvtx   = g->nvtx ;
nstage = info->nstage ;
if ( info->msglvl > 2 ) {
   fprintf(info->msgFile, 
           "\n\n MSMD object initialized, %d stages", nstage) ;
   fflush(info->msgFile) ;
}
/*
   ------------------------------------
   load the reach set with all vertices
   ------------------------------------
*/
if ( info->compressFlag / 4 >= 1 ) {
/*
   ------------------
   compress the graph
   ------------------
*/
   if ( info->msglvl > 2 ) {
      fprintf(info->msgFile, "\n\n initial compression") ;
      fflush(info->msgFile) ;
   }
   IV_setSize(&msmd->reachIV, nvtx) ;
   IV_ramp(&msmd->reachIV, 0, 1) ;
   MSMD_findInodes(msmd, info) ;
   if ( info->msglvl > 2 ) {
      fprintf(info->msgFile, 
              "\n\n %d checked, %d found indistinguishable",
              info->stageInfo->ncheck, info->stageInfo->nindst) ;
      fflush(info->msgFile) ;
   }
   MSMD_cleanReachSet(msmd, info) ;
/*
   for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
      MSMD_cleanEdgeList(msmd, v, info) ;
   }
*/
}
IV_setSize(&msmd->reachIV, 0) ;
/*
   --------------------
   loop over the stages
   --------------------
*/
for ( info->istage = 0 ; info->istage <= nstage ; info->istage++ ) {
   if ( info->msglvl > 2 ) {
      fprintf(info->msgFile, 
              "\n\n ##### elimination stage %d", info->istage) ;
      fflush(info->msgFile) ;
   }
/*
   if ( info->istage == nstage ) {
      info->msglvl = 5 ;
   }
*/
   MARKTIME(t1) ;
   MSMD_eliminateStage(msmd, info) ;
   MARKTIME(t2) ;
   info->stageInfo->cpu = t2 - t1 ;
   info->stageInfo++ ;
}
/*
   -------------
   final cleanup 
   -------------
*/
{
MSMDvtx   *first, *last, *v ;

IV_setSize(&msmd->reachIV, 0) ;
first = msmd->vertices ;
last  = first + nvtx - 1 ;
for ( v = first ; v <= last ; v++ ) {
   switch ( v->status ) {
   case 'E' :
   case 'L' :
   case 'I' :
      break ;
   default :
      IV_push(&msmd->reachIV, v->id) ;
      break ;
   }
}
/*
fprintf(stdout, "\n reach set, %d entries", IV_size(&msmd->reachIV)) ;
IV_writeForHumanEye(&msmd->reachIV, stdout) ;
*/
MSMD_findInodes(msmd, info) ;
}
/*
   ---------------------------------------------------------
   make info->stagInfo point back to beginning of the stages
   ---------------------------------------------------------
*/
info->stageInfo -= nstage + 1 ;
/*
   ------------------------
   get the total statistics
   ------------------------
*/
for ( istage = 0, now = info->stageInfo,
      total = info->stageInfo + nstage + 1 ; 
      istage <= nstage ; 
      istage++, now++ ) {
   total->nstep    += now->nstep    ;
   total->nfront   += now->nfront   ;
   total->welim    += now->welim    ;
   total->nfind    += now->nfind    ;
   total->nzf      += now->nzf      ;
   total->ops      += now->ops      ;
   total->nexact2  += now->nexact2  ;
   total->nexact3  += now->nexact3  ;
   total->napprox  += now->napprox  ;
   total->ncheck   += now->ncheck   ;
   total->nindst   += now->nindst   ;
   total->noutmtch += now->noutmtch ;
}
/*
   -----------------------------------------------------
   free some working storage (leave the MSMDvtx objects
   so the user can extract the permutations and/or ETree
   -----------------------------------------------------
*/
IIheap_free(msmd->heap) ;
msmd->heap = NULL ;
IV_clearData(&msmd->ivtmpIV) ;
IV_clearData(&msmd->reachIV) ;
/*
while ( (ip = msmd->baseIP) != NULL ) {
   msmd->baseIP = ip->next ;
   IP_free(ip) ;
}
*/
MARKTIME(t3) ;
info->totalCPU = t3 - t0 ;

return ; }

/*--------------------------------------------------------------------*/
