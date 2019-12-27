/*  IO.c  */

#include "../I2Ohash.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to write the I2Ohash object to a file

   created -- 98jan29, cca
   ------------------------------------------------
*/
void
I2Ohash_writeForHumanEye ( 
   I2Ohash   *hashtable,
   FILE     *fp 
) {
double    measure ;
int       count, loc, nfull ;
I2OP      *i2op ;

if ( hashtable == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in I2Ohash_writeForHumanEye(%p,%p)"
           "\n hashtable is NULL or file pointer is NULL", 
           hashtable, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n\n I2Ohash : %d lists, %d items", 
        hashtable->nlist, hashtable->nitem) ;
nfull   = 0 ;
measure = 0.0 ;
for ( loc = 0 ; loc < hashtable->nlist ; loc++ ) {
   if ( (i2op = hashtable->heads[loc]) != NULL ) {
      fprintf(fp, "\n %4d : ", loc) ;
      count = 0 ;
      while ( i2op != NULL ) {
         if ( ++count % 4 == 0 ) {
            fprintf(fp, "\n") ;
         }
         fprintf(fp, " < %6d, %6d, %p >", 
                 i2op->value0, i2op->value1, i2op->value2) ;
         i2op = i2op->next ;
      }
      measure += count*count ;
      nfull++ ;
   }
}
measure = sqrt(measure) ;
fprintf(fp, "\n %d empty lists, %d items, %.3f ratio",
        nfull, hashtable->nitem, 
        measure*sqrt((double) hashtable->nlist)/hashtable->nitem) ;

return ; }

/*--------------------------------------------------------------------*/
