/*  basics.c  */

#include "../Chv.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98apr30, cca
   -----------------------
*/
Chv *
Chv_new ( 
   void 
) {
Chv   *chv ;

ALLOCATE(chv, struct _Chv, 1) ;
Chv_setDefaultFields(chv) ;

return(chv) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98apr30, cca
   -----------------------
*/
void
Chv_setDefaultFields (
   Chv   *chv
) {
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_setDefaultFields(%p)"
           "\n bad input", chv) ;
   exit(-1) ;
}
chv->id      =  -1  ;
chv->nD      =   0  ;
chv->nL      =   0  ;
chv->nU      =   0  ;
chv->type    = SPOOLES_REAL ;
chv->symflag = SPOOLES_SYMMETRIC ;
DV_setDefaultFields(&chv->wrkDV) ;
chv->rowind  = NULL ;
chv->colind  = NULL ;
chv->entries = NULL ;
chv->next    = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98apr30, cca
   --------------------------------------------------
*/
void
Chv_clearData ( 
   Chv   *chv 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_clearData(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_clearData(&chv->wrkDV) ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
Chv_setDefaultFields(chv) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98apr30, cca
   ------------------------------------------
*/
void
Chv_free ( 
   Chv   *chv 
) {
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_free(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
Chv_clearData(chv) ;
FREE(chv) ;

return ; }

/*--------------------------------------------------------------------*/
