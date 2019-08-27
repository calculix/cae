/*  init.c  */

#include "../I2Ohash.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   initializer, 
   (1) set hashtable maximum size to nlist
       and allocate the heads[] vector
   (2) initialize nobj I2OP objects
   (2) set item growth factor to grow

   created -- 98jan29, cca
   ------------------------------------------
*/
void
I2Ohash_init ( 
   I2Ohash   *hashtable,
   int       nlist,
   int       nobj,
   int       grow 
) {
int   ii ;
if ( hashtable == NULL || nlist <= 0 ) {
   fprintf(stderr, "\n\n error in I2Ohash_init(%p,%d,%d,%d)"
           "\n hashtable is NULL or nlist = %d is nonpositive\n",
           hashtable, nlist, nobj, grow, nlist) ;
   exit(-1) ;
}
if ( nobj <= 0 && grow <= 0 ) {
   fprintf(stderr, "\n\n error in I2Ohash_init(%p,%d,%d,%d)"
           "\n nobj = %d, grow = %d\n",
           hashtable, nlist, nobj, grow, nobj, grow) ;
   exit(-1) ;
}
I2Ohash_clearData(hashtable) ;
hashtable->nlist = nlist ;
hashtable->grow  = grow ;
if ( nobj > 0 ) {
   hashtable->baseI2OP = I2OP_init(nobj+1, I2OP_FORWARD) ;
   hashtable->freeI2OP = hashtable->baseI2OP + 1 ;
   hashtable->baseI2OP->next = NULL ;
}
ALLOCATE(hashtable->heads, struct _I2OP *, nlist) ;
for ( ii = 0 ; ii < nlist ; ii++ ) {
   hashtable->heads[ii] = NULL ;
}

return ; }

/*--------------------------------------------------------------------*/
