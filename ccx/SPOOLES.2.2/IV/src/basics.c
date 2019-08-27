/*  basics.C  */

#include "../IV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor method

   created -- 95oct06, cca
   -----------------------
*/
IV *
IV_new ( 
   void 
) {
IV   *iv ;

ALLOCATE(iv, struct _IV, 1) ;

IV_setDefaultFields(iv) ;

return(iv) ; }
/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 95oct06, cca
   -----------------------
*/
void
IV_setDefaultFields ( 
   IV   *iv
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_setDefaultFields(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
iv->size    =   0  ;
iv->maxsize =   0  ;
iv->owned   =   0  ;
iv->vec     = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 95oct06, cca
   -----------------------
*/
void
IV_clearData ( 
   IV   *iv
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_clearData(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
if ( iv->vec != NULL && iv->owned == 1 ) {
   IVfree(iv->vec) ;
}
IV_setDefaultFields(iv) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   destructor

   created -- 95oct06, cca
   -----------------------
*/
void
IV_free ( 
   IV   *iv
) {
if ( iv == NULL ) {
   fprintf(stderr, "\n fatal error in IV_free(%p)"
           "\n bad input\n", iv) ;
   exit(-1) ;
}
IV_clearData(iv) ;
FREE(iv) ;

return ; }

/*--------------------------------------------------------------------*/
