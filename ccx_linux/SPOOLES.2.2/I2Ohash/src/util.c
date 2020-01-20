/*  util.c  */

#include "../I2Ohash.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   insert the (key1, key2, value) triple into the hash table

   created -- 98jan29, cca
   ---------------------------------------------------------
*/
void
I2Ohash_insert ( 
   I2Ohash   *hashtable,
   int       key1, 
   int       key2, 
   void      *value 
) {
int   loc, loc1, loc2 ;
I2OP  *i2op, *j2op, *prev ;

if ( hashtable == NULL ) {
   fprintf(stderr, "\n error in I2Ohash_insert(%p,%d,%d,%p)"
           "\n hashtable is NULL \n", hashtable, key1, key2, value) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n inside I2Ohash_insert(%p,%d,%d,%p)",
        hashtable, key1, key2, value) ;
fflush(stdout) ;
#endif
/*
   ---------------------------------------------------
   get the list to hold the (key1, key2, value) triple
   ---------------------------------------------------
*/
loc1 = (key1 + 1) % hashtable->nlist ;
loc2 = (key2 + 1) % hashtable->nlist ;
/*loc  = (loc1*loc2) % hashtable->nlist ;*/
long int loc3  = (long int)loc1*(long int)loc2 % hashtable->nlist ;
loc=(int)loc3;
#if MYDEBUG > 0
fprintf(stdout, "\n loc1 = %d, loc2 = %d, loc3 = %d", loc1, loc2, loc) ;
fflush(stdout) ;
#endif
/*
   --------------------------------------------------------
   get the list item to hold the (key1, key2, value) triple
   --------------------------------------------------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n hashtable->freeI2OP = %p", hashtable->freeI2OP) ;
fflush(stdout) ;
#endif
if ( (i2op = hashtable->freeI2OP) == NULL ) {
   if ( hashtable->grow <= 0 ) { 
      fprintf(stderr, "\n fatal error in I2Ohash_insert(%p,%d,%d,%p)"
              "\n no free list items, grow = %d",
              hashtable, key1, key2, value, hashtable->grow) ;
      exit(-1) ;
   }
   i2op = I2OP_init(hashtable->grow, I2OP_FORWARD) ;
   hashtable->freeI2OP = i2op + 1 ;
   i2op->next = hashtable->baseI2OP ;
   hashtable->baseI2OP = i2op ;
   i2op = hashtable->freeI2OP ;
}
hashtable->freeI2OP = i2op->next ;
#if MYDEBUG > 0
fprintf(stdout, "\n i2op = %p, hashtable->freeI2OP = %p", 
        i2op, hashtable->freeI2OP) ;
fflush(stdout) ;
#endif
i2op->value0 = key1  ;
i2op->value1 = key2  ;
i2op->value2 = value ;
i2op->next   = NULL  ;
#if MYDEBUG > 0
fprintf(stdout, "\n i2op : value0 %d, value1 %d, value2 %p, next %p",
        i2op->value0, i2op->value1, i2op->value2, i2op->next) ;
fflush(stdout) ;
#endif
/*
   --------------------------------------------
   insert the (key1, key2, value) into its list
   --------------------------------------------
*/
for ( j2op = hashtable->heads[loc], prev = NULL ;
      j2op != NULL ; j2op = j2op->next ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n j2op : value0 %d, value1 %d, value2 %p, next %p",
           j2op->value0, j2op->value1, j2op->value2, j2op->next) ;
   fflush(stdout) ;
#endif
   if (  j2op->value0 > key1
      || (j2op->value0 == key1 && j2op->value1 >= key2 ) ) {
      break ;
   }
   prev = j2op ;
}
if ( prev == NULL ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n   heads[%d] = %p", loc, i2op) ;
#endif
   hashtable->heads[loc] = i2op ;
} else {
#if MYDEBUG > 0
   fprintf(stdout, "\n   %p->next = %p", prev, i2op) ;
#endif
   prev->next = i2op ;
}
#if MYDEBUG > 0
   fprintf(stdout, "\n   %p->next = %p", i2op, j2op) ;
#endif
i2op->next = j2op ;
hashtable->nitem++ ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   locate the first (key1,key2,*) in the hash table

   return value --
      0 -- no (key1,key2,*) value triple
      1 -- (key1,key2,value) value triple found,
           value put into *pvalue

   created -- 98jan29, cca
   ------------------------------------------------
*/
int
I2Ohash_locate ( 
   I2Ohash   *hashtable,
   int       key1,
   int       key2,
   void      **pvalue 
) {
int   loc, loc1, loc2, rc ;
I2OP   *i2op ;
/*
   ---------------
   check the input
   ---------------
*/
if ( hashtable == NULL || pvalue == NULL ) {
   fprintf(stderr, "\n error in I2Ohash_locate(%p,%d,%d,%p)"
           "\n hashtable or pvalue is NULL\n",
           hashtable, key1, key2, pvalue) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n inside I2Ohash_locate(%p,%d,%d,%p)",
        hashtable, key1, key2, pvalue) ;
fflush(stdout) ;
#endif
loc1 = (key1 + 1) % hashtable->nlist ;
loc2 = (key2 + 1) % hashtable->nlist ;
/*loc  = (loc1*loc2) % hashtable->nlist ;*/
long int loc3  = (long int)loc1*(long int)loc2 % hashtable->nlist ;
loc=(int)loc3;
#if MYDEBUG > 0
fprintf(stdout, "\n loc1 = %d, loc2 = %d, loc3 = %d", loc1, loc2, loc) ;
fflush(stdout) ;
#endif
/*
   ---------------------------------------------
   find the location of the first (key,*) triple
   ---------------------------------------------
*/
for ( i2op = hashtable->heads[loc] ; 
      i2op != NULL ; 
      i2op = i2op->next ) {
#if MYDEBUG > 0
   fprintf(stdout, 
      "\n  i2op = %p, value0 = %d, value1 = %d, value2 = %p, next = %p",
           i2op, i2op->value0, i2op->value1, i2op->value2, i2op->next) ;
   fflush(stdout) ;
#endif
   if (  i2op->value0 > key1 
     || (i2op->value0 == key1 && i2op->value1 >= key2) ) {
      break ;
   }
}
rc = 0 ;
if ( i2op != NULL && i2op->value0 == key1 && i2op->value1 == key2 ) {
/*
   --------------------------
   (key,*) found, set *pvalue
   --------------------------
*/
   *pvalue = i2op->value2 ;
   rc = 1 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   remove the first (key1,key2,*) from the hash table

   return value --
      0 -- no (key1,key2,*) value triple
      1 -- (key1,key2,value) value triple found,
           value put into *pvalue

   created -- 98jan29, cca
   --------------------------------------------------
*/
int
I2Ohash_remove ( 
   I2Ohash   *hashtable,
   int       key1,
   int       key2,
   void      **pvalue 
) {
int   loc, loc1, loc2, rc ;
I2OP   *i2op, *prev ;
/*
   ---------------
   check the input
   ---------------
*/
if ( hashtable == NULL || pvalue == NULL ) {
   fprintf(stderr, "\n error in I2Ohash_remove(%p,%d,%d,%p)"
           "\n hashtable or pvalue is NULL\n",
           hashtable, key1, key2, pvalue) ;
   exit(-1) ;
}
loc1 = (key1 + 1) % hashtable->nlist ;
loc2 = (key2 + 1) % hashtable->nlist ;
/*loc  = (loc1*loc2) % hashtable->nlist ;*/
long int loc3  = (long int)loc1*(long int)loc2 % hashtable->nlist ;
loc=(int)loc3;
/*
   ---------------------------------------------------
   find the location of the first (key1,key2,*) triple
   ---------------------------------------------------
*/
for ( i2op = hashtable->heads[loc], prev = NULL ;
      i2op != NULL ;
      i2op = i2op->next ) {
   if (  i2op->value0 > key1 
     || (i2op->value0 == key1 && i2op->value1 >= key2) ) {
      break ;
   }
   prev = i2op ;
}
rc = 0 ;
if ( i2op != NULL && i2op->value0 == key1 && i2op->value1 == key2 ) {
/*
   ----------------------------------
   (key,*) found, remove, set *pvalue
   ----------------------------------
*/
   if ( prev == NULL ) {
      hashtable->heads[loc] = i2op->next ;
   } else {
      prev->next = i2op->next ;
   }
   i2op->next = hashtable->freeI2OP ;
   hashtable->freeI2OP = i2op ;
   hashtable->nitem-- ;
   *pvalue = i2op->value2 ;
   rc = 1 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   return a measure of the nonuniform distribution of the entries.
   a value of 1.0 is best.

   created -- 98jan29, cca
   ---------------------------------------------------------------
*/
double
I2Ohash_measure (
   I2Ohash   *hashtable
) {
double   measure ;
int      count, loc ;
I2OP      *i2op ;

if ( hashtable == NULL ) {
   fprintf(stderr, "\n fatal error in I2Ohash_measure(%p)"
           "\n hashtable is NULL\n", hashtable) ;
   exit(-1) ;
}
measure = 0.0 ;
for ( loc = 0 ; loc < hashtable->nlist ; loc++ ) {
   if ( (i2op = hashtable->heads[loc]) != NULL ) {
      count = 0 ;
      while ( i2op != NULL ) {
         count++ ;
         i2op = i2op->next ;
      }
      measure += count*count ;
   }
}
measure = sqrt(measure) ;
measure *= sqrt((double) hashtable->nlist)/hashtable->nitem ;

return(measure) ; }

/*--------------------------------------------------------------------*/
