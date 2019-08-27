/*  basics.c  */

#include "../SolveMap.h"

#define   MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simplest constructor

   created -- 98mar19, cca
   -----------------------
*/
SolveMap *
SolveMap_new ( 
   void 
) {
SolveMap   *solvemap ;

ALLOCATE(solvemap, struct _SolveMap, 1) ;
SolveMap_setDefaultFields(solvemap) ;

return(solvemap) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98mar19, cca
   -----------------------
*/
void
SolveMap_setDefaultFields (
   SolveMap   *solvemap
) {
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_setDefaultFields(%p)"
           "\n bad input", solvemap) ;
   exit(-1) ;
}
solvemap->symmetryflag = SPOOLES_SYMMETRIC ;
solvemap->nfront       = 0 ;
solvemap->nproc        = 0 ;
solvemap->owners       = NULL ;
solvemap->nblockUpper  = 0 ;
solvemap->rowidsUpper  = NULL ;
solvemap->colidsUpper  = NULL ;
solvemap->mapUpper     = NULL ;
solvemap->nblockLower  = 0 ;
solvemap->rowidsLower  = NULL ;
solvemap->colidsLower  = NULL ;
solvemap->mapLower     = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 98mar19, cca
   --------------------------------------------------
*/
void
SolveMap_clearData ( 
   SolveMap   *solvemap 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_clearData(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
/*
   -----------------------------------------------
   free any storage held in the int vector objects
   -----------------------------------------------
*/
if ( solvemap->owners != NULL ) {
   IVfree(solvemap->owners) ;
}
if ( solvemap->rowidsUpper != NULL ) {
   IVfree(solvemap->rowidsUpper) ;
}
if ( solvemap->colidsUpper != NULL ) {
   IVfree(solvemap->colidsUpper) ;
}
if ( solvemap->mapUpper != NULL ) {
   IVfree(solvemap->mapUpper) ;
}
if ( solvemap->rowidsLower != NULL ) {
   IVfree(solvemap->rowidsLower) ;
}
if ( solvemap->colidsLower != NULL ) {
   IVfree(solvemap->colidsLower) ;
}
if ( solvemap->mapLower != NULL ) {
   IVfree(solvemap->mapLower) ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
SolveMap_setDefaultFields(solvemap) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 98mar19, cca
   ------------------------------------------
*/
SolveMap *
SolveMap_free ( 
   SolveMap   *solvemap 
) {
if ( solvemap == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_free(%p)"
           "\n bad input\n", solvemap) ;
   exit(-1) ;
}
SolveMap_clearData(solvemap) ;
FREE(solvemap) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
