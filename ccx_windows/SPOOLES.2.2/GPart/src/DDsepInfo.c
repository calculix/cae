/*  DDsepInfo.c  */

#include "../DDsepInfo.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   construct a new instance of the DDsepInfo object

   created -- 96feb24, cca
   ------------------------------------------------
*/
DDsepInfo *
DDsepInfo_new (
   void
) {
DDsepInfo   *info ;

ALLOCATE(info, struct _DDsepInfo, 1) ;

DDsepInfo_setDefaultFields(info) ;

return(info) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   set the default fields of the DDsepInfo object

   created  -- 96feb24, cca
   ---------------------------------------------
*/
void
DDsepInfo_setDefaultFields (
   DDsepInfo   *info
) {
if ( info == NULL ) {
   fprintf(stderr, "\n fatal error in DDsepInfo_setDefaultFields(%p)"
           "\n bad input\n", info) ;
   exit(-1) ;
}
info->seed          =      1 ;
info->minweight     =     40 ;
info->maxweight     =     80 ;
info->freeze        =    4.0 ;
info->alpha         =    1.0 ;
info->maxcompweight =    500 ;
info->ntreeobj      =      0 ;
info->DDoption      =      1 ;
info->nlayer        =      3 ;
info->cpuDD         =    0.0 ;
info->cpuMap        =    0.0 ;
info->cpuBPG        =    0.0 ;
info->cpuBKL        =    0.0 ;
info->cpuSmooth     =    0.0 ;
info->cpuSplit      =    0.0 ;
info->cpuTotal      =    0.0 ;
info->msglvl        =      0 ;
info->msgFile       = stdout ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   clear the data fields for a DDsepInfo object

   created  -- 96feb24, cca
   ---------------------------------------------
*/
void
DDsepInfo_clearData (
   DDsepInfo   *info
) {
if ( info == NULL ) {
   fprintf(stderr, "\n fatal error in DDsepInfo_clearData(%p)"
           "\n bad input\n", info) ;
   exit(-1) ;
}
DDsepInfo_setDefaultFields(info) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------
   free the DDsepInfo object

   created  -- 96feb24, cca
   ------------------------
*/
void
DDsepInfo_free (
   DDsepInfo   *info
) {
if ( info == NULL ) {
   fprintf(stderr, "\n fatal error in DDsepInfo_free(%p)"
           "\n bad input\n", info) ;
   exit(-1) ;
}
DDsepInfo_clearData(info) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n trying to free info") ;
   fflush(stdout) ;
#endif
FREE(info) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   write the CPU times
  
   created -- 97nov06, cca
   -----------------------
*/
void
DDsepInfo_writeCpuTimes (
   DDsepInfo   *info,
   FILE        *msgFile
) {
double   cpuMisc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( info == NULL || msgFile == NULL ) {
   fprintf(stderr, "\n fatal error in DDsepInfo_writeCpuTimes(%p,%p)"
           "\n bad input\n", info, msgFile) ;
   exit(-1) ;
}
cpuMisc = info->cpuTotal - info->cpuDD  - info->cpuSplit - info->cpuMap 
        - info->cpuBPG   - info->cpuBKL - info->cpuSmooth ;
if ( info->cpuTotal > 0 ) {
fprintf(msgFile, 
        "\n\n CPU breakdown for graph partition"
        "\n               raw CPU   per cent"
        "\n misc       : %9.2f %6.1f%%"
        "\n Split      : %9.2f %6.1f%%"
        "\n find DD    : %9.2f %6.1f%%"
        "\n DomSeg Map : %9.2f %6.1f%%"
        "\n DomSeg BPG : %9.2f %6.1f%%"
        "\n BKL        : %9.2f %6.1f%%"
        "\n Smooth     : %9.2f %6.1f%%"
        "\n Total      : %9.2f %6.1f%%",
        cpuMisc,         100.*cpuMisc/info->cpuTotal,
        info->cpuSplit,  100.*info->cpuSplit/info->cpuTotal,
        info->cpuDD,     100.*info->cpuDD/info->cpuTotal,
        info->cpuMap,    100.*info->cpuMap/info->cpuTotal,
        info->cpuBPG,    100.*info->cpuBPG/info->cpuTotal,
        info->cpuBKL,    100.*info->cpuBKL/info->cpuTotal,
        info->cpuSmooth, 100.*info->cpuSmooth/info->cpuTotal,
        info->cpuTotal,  100.) ;
}
return ; }

/*--------------------------------------------------------------------*/
