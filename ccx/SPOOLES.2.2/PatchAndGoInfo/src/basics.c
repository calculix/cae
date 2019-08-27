/*  basics.c  */

#include "../PatchAndGoInfo.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor method

   created -- 98aug26, cca
   -----------------------
*/
PatchAndGoInfo *
PatchAndGoInfo_new ( 
   void 
) {
PatchAndGoInfo   *info ;

ALLOCATE(info, struct _PatchAndGoInfo, 1) ;

PatchAndGoInfo_setDefaultFields(info) ;

return(info) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98aug26, cca
   -----------------------
*/
void
PatchAndGoInfo_setDefaultFields ( 
   PatchAndGoInfo   *info
) {
if ( info == NULL ) {
   fprintf(stderr, "\n fatal error in PatchAndGoInfo_setDefaultFields()"
           "\n bad input\n") ;
   exit(-1) ;
}
info->strategy =  -1  ;
info->toosmall =  0.0 ;
info->fudge    =  0.0 ;
info->fudgeIV  = NULL ;
info->fudgeDV  = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 98aug26, cca
   -----------------------
*/
void
PatchAndGoInfo_clearData ( 
   PatchAndGoInfo   *info
) {
if ( info == NULL ) {
   fprintf(stderr, "\n fatal error in PatchAndGoInfo_clearData()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( info->fudgeIV != NULL ) {
   IV_free(info->fudgeIV) ;
}
if ( info->fudgeDV != NULL ) {
   DV_free(info->fudgeDV) ;
}
PatchAndGoInfo_setDefaultFields(info) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 98aug26, cca
   -----------------------
*/
void
PatchAndGoInfo_free ( 
   PatchAndGoInfo   *info
) {
if ( info == NULL ) {
   fprintf(stderr, "\n fatal error in PatchAndGoInfo_free()"
           "\n bad input\n") ;
   exit(-1) ;
}
PatchAndGoInfo_clearData(info) ;
FREE(info) ;

return ; }

/*--------------------------------------------------------------------*/
