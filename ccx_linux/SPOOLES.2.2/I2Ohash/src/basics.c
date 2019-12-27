/*  basics.c  */

#include "../I2Ohash.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   create and return a new instance of the I2Ohash object

   created -- 98jan29, cca
   ------------------------------------------------------
*/
I2Ohash *
I2Ohash_new (
   void
) {
I2Ohash   *hashtable ;

ALLOCATE(hashtable, struct _I2Ohash, 1) ;

I2Ohash_setDefaultFields(hashtable) ;

return(hashtable) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   set the default fields for an I2Ohash object

   created -- 98jan29, cca
   --------------------------------------------
*/
void
I2Ohash_setDefaultFields (
   I2Ohash   *hashtable
) {
if ( hashtable == NULL ) {
   fprintf(stderr, "\n fatal error in I2Ohash_setDefaultFields(%p)"
           "\n hashtable is NULL\n", hashtable) ;
   exit(-1) ;
}
hashtable->nlist    =   0  ;
hashtable->grow     =   0  ;
hashtable->nitem    =   0  ;
hashtable->baseI2OP = NULL ;
hashtable->freeI2OP = NULL ;
hashtable->heads    = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 98jan29, cca
   -----------------------
*/
void
I2Ohash_clearData (
   I2Ohash   *hashtable
) {
I2OP   *i2op ;

if ( hashtable == NULL ) {
   fprintf(stderr, "\n fatal error in I2Ohash_clearData(%p)"
           "\n hashtable is NULL\n", hashtable) ;
   exit(-1) ;
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n I2Ohash_clearData") ;
   fflush(stdout) ;
#endif
while ( (i2op = hashtable->baseI2OP) != NULL ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n    baseI2OP = %p", i2op) ;
   fflush(stdout) ;
#endif
   hashtable->baseI2OP = i2op->next ;
   I2OP_free(i2op) ;
}
if ( hashtable->heads != NULL ) {
   FREE(hashtable->heads) ;
}
I2Ohash_setDefaultFields(hashtable) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   free the I2Ohash object

   created -- 98jan29, cca
   -----------------------
*/
void
I2Ohash_free (
   I2Ohash   *hashtable
) {
if ( hashtable == NULL ) {
   fprintf(stderr, "\n fatal error in I2Ohash_free(%p)"
           "\n hashtable is NULL\n", hashtable) ;
   exit(-1) ;
}
I2Ohash_clearData(hashtable) ;
FREE(hashtable) ;

return ; }

/*--------------------------------------------------------------------*/
