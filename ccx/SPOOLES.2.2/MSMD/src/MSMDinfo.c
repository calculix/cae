/*  basics.C  */

#include "../MSMD.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor

   created -- 96feb25, cca
   -----------------------
*/
MSMDinfo *
MSMDinfo_new ( 
   void 
) {
MSMDinfo   *info ;

ALLOCATE(info, struct _MSMDinfo, 1) ;
MSMDinfo_setDefaultFields(info) ;

return(info) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   set the default data fields
   
   created -- 96feb25, cca
   ---------------------------
*/
void
MSMDinfo_setDefaultFields(
   MSMDinfo   *info
) {
info->compressFlag =    1    ;
info->prioType     =    1    ;
info->stepType     =  1.0    ;
info->seed         =    0    ;
info->msglvl       =    0    ;
info->msgFile      = stdout  ;
info->maxnbytes    =    0    ;
info->nbytes       =    0    ;
info->istage       =    0    ;
info->nstage       =    0    ;
info->stageInfo    =  NULL   ;
info->totalCPU     =   0.0   ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 96feb25, cca
   -----------------------
*/
void
MSMDinfo_clearData ( 
   MSMDinfo   *info
) {
if ( info->stageInfo  != NULL ) { 
   FREE(info->stageInfo ) ; 
}
MSMDinfo_setDefaultFields(info) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 96feb25, cca
   -----------------------
*/
void
MSMDinfo_free (
   MSMDinfo   *info
) {
MSMDinfo_clearData(info) ;
FREE(info) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- print the MSMDinfo object

   created -- 96feb25, cca
   ------------------------------------
*/
void
MSMDinfo_print ( 
   MSMDinfo    *info,
   FILE        *fp 
) {
int             istage ;
MSMDstageInfo   *stageinfo ;

if ( info == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in MSMDinfo_print(%p,%p)"
           "\n bad input\n", info, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n\n MSMDinfo :") ;
fprintf(fp, "\n    compressFlag = %d : ", info->compressFlag) ;
if ( info->compressFlag / 4  >= 1 ) {
   fprintf(fp, "compress graph, ") ;
}
switch ( info->compressFlag % 4 ) {
case 0 : 
   fprintf(fp, "during elimination do not compress") ; 
   break ;
case 1 : 
   fprintf(fp, "during elimination compress 2-adj nodes") ; 
   break ;
case 2 : 
   fprintf(fp, "during elimination compress all nodes") ; 
   break ;
default :
   fprintf(fp, "\n unknown type") ;
   break ;
}
fprintf(fp, "\n    prioType = %d : ", info->prioType) ;
switch ( info->prioType ) {
case 1 : 
   fprintf(fp, " true updates") ; 
   break ;
case 2 : 
   fprintf(fp, " approximate updates") ; 
   break ;
case 3 : 
   fprintf(fp, " true updates for 2-adj nodes, others approximate") ;
   break ;
default :
   fprintf(fp, " unknown type") ;
   break ;
}
fprintf(fp, "\n    stepType = %f : ", info->stepType) ;
if ( info->stepType < 1.0 ) {
   fprintf(fp, " single elimination") ;
} else if ( info->stepType == 1.0 ) {
   fprintf(fp, " multiple elimination of nodes of mininum degree") ;
} else {
   fprintf(fp, " multiple elimination in range [mindeg, %f*mindeg]",
          info->stepType) ;
}
fprintf(fp, "\n    msglvl       = %d ", info->msglvl) ;
fprintf(fp, "\n    maxnbytes    = %d ", info->maxnbytes) ;
fprintf(fp, "\n    ordering cpu = %8.3f ", info->totalCPU) ;
fprintf(fp, "\n    stage information") ;
fprintf(fp, 
"\n\n stage #steps #fronts #weight #frontind     nzf          ops    CPU") ;
for ( istage = 0, stageinfo = info->stageInfo  ;
      istage <= info->nstage ; istage++, stageinfo++ ) {
   fprintf(fp, "\n   %3d %5d %6d %7d %9d %10d %12.0f %8.3f",
           istage, stageinfo->nstep, stageinfo->nfront,
           stageinfo->welim, stageinfo->nfind, stageinfo->nzf,
           stageinfo->ops, stageinfo->cpu) ;
}
fprintf(fp, "\n total %5d %6d %7d %9d %10d %12.0f ",
        stageinfo->nstep, stageinfo->nfront,
        stageinfo->welim, stageinfo->nfind, stageinfo->nzf,
        stageinfo->ops) ;
fprintf(fp, 
    "\n\n stage #nexact2 #exact3 #approx #check #indst #outmatched") ;
for ( istage = 0, stageinfo = info->stageInfo  ;
      istage <= info->nstage ; istage++, stageinfo++ ) {
   fprintf(fp, "\n   %3d %6d %7d %6d %7d %8d %8d",
           istage, stageinfo->nexact2, stageinfo->nexact3,
           stageinfo->napprox, stageinfo->ncheck, stageinfo->nindst,
           stageinfo->noutmtch) ;
}
fprintf(fp, "\n total %6d %7d %6d %7d %8d %8d",
        stageinfo->nexact2, stageinfo->nexact3,
        stageinfo->napprox, stageinfo->ncheck, stageinfo->nindst,
        stageinfo->noutmtch) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   determine if the MSMDinfo object is valid

   created -- 96feb25, cca
   -----------------------------------------
*/
int
MSMDinfo_isValid (
   MSMDinfo   *info
) {
int   rc ;

if (  info == NULL 
   || info->compressFlag < 0 
   || info->compressFlag == 3 
   || info->compressFlag > 6
   || info->prioType < 1
   || info->prioType > 4 ) {
   rc = 0 ;
} else {
   rc = 1 ;
}

return(rc) ; }

/*--------------------------------------------------------------------*/
