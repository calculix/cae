/*  I2OP.c  */

#include "../Utilities.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   initializer. 
   create and return an array of n I2OP structures.
   the structures are linked in one of three ways.
   flag = 0 (I2OP_NULL)     --> ip->next = NULL
   flag = 1 (I2OP_FORWARD)  --> ip->next = successor in list
   flag = 2 (I2OP_BACKWARD) --> ip->next = predecessor in list
   
   created -- 98feb06, cca
   ---------------------------------------------------------
*/
I2OP *
I2OP_init ( 
   int   n, 
   int   flag 
) {
I2OP    *base ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n <= 0 || flag < I2OP_NULL || flag > I2OP_BACKWARD ) {
   fprintf(stderr, "\n fatal error in I2OP_init(%d,%d)"
           "\n bad input\n", n, flag) ;
   exit(-1) ;
}
/*
   --------------------
   allocate the storage
   --------------------
*/
ALLOCATE(base, struct _I2OP, n) ;
/*
   -----------------------
   initialize the elements
   -----------------------
*/
I2OP_initStorage(n, flag, base) ;

return(base) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   initializer. 
   create and return an array of n I2OP structures.
   the structures are linked in one of three ways.
   flag = 0 (I2OP_NULL)     --> ip->next = NULL
   flag = 1 (I2OP_FORWARD)  --> ip->next = successor in list
   flag = 2 (I2OP_BACKWARD) --> ip->next = predecessor in list
   
   created -- 98feb06, cca
   ---------------------------------------------------------
*/
void
I2OP_initStorage ( 
   int   n, 
   int   flag,
   I2OP   *base
) {
I2OP   *elem, *firstelem, *lastelem ;
/*
   ---------------
   check the input
   ---------------
*/
if (  n <= 0 || flag < I2OP_NULL || flag > I2OP_BACKWARD 
   || base == NULL ) {
   fprintf(stderr, "\n fatal error in I2OP_initStorage(%d,%d,%p)"
   "\n bad input\n", n, flag, base) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n inside I2OP_initStorage(%d, %d, %p)",
        n, flag, base) ;
#endif
/*
   -----------------------
   initialize the elements
   -----------------------
*/
firstelem = base ;
lastelem  = firstelem + n - 1 ;
#if MYDEBUG > 0
fprintf(stdout, "\n firstelem %p, lastelem %p",
        firstelem, lastelem) ;
#endif
switch ( flag ) {
case I2OP_FORWARD :
   for ( elem = firstelem ; elem < lastelem ; elem++ ) {
      elem->value0 = elem->value1 = -1 ;
      elem->value2 = NULL ;
      elem->next   = elem + 1 ;
#if MYDEBUG > 0
      fprintf(stdout, 
           "\n elem %p, value0 %d, value1 %d, value2 %p, next %p",
           elem, elem->value0, elem->value1, elem->value2, elem->next) ;
#endif
   }
   lastelem->value0 = lastelem->value1 = -1 ;
   lastelem->value2 = NULL ;
   lastelem->next   = NULL ;
#if MYDEBUG > 0
   fprintf(stdout, 
        "\n elem %p, value0 %d, value1 %d, value2 %p, next %p",
        lastelem, lastelem->value0, lastelem->value1, 
        lastelem->value2, lastelem->next) ;
#endif
   break ;
case I2OP_BACKWARD :
   for ( elem = firstelem + 1 ; elem <= lastelem ; elem++ ) {
      elem->value0 = elem->value1 = -1 ;
      elem->value2 = NULL ;
      elem->next   = elem - 1 ;
   }
   firstelem->value0 = firstelem->value1 = -1 ;
   firstelem->value2 = NULL ;
   firstelem->next   = NULL ;
   break ;
case I2OP_NULL :
   for ( elem = firstelem ; elem <= lastelem ; elem++ ) {
      elem->value0 = elem->value1 = -1 ;
      elem->value2 = NULL ;
      elem->next   = NULL ;
   }
   break ;
default :
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   free the storage for an array of I2OP structures,
   must have been allocated by I2OP_init

   created -- 98feb06, cca
   -----------------------------------------------
*/
void
I2OP_free ( 
   I2OP   *ip
) {
if ( ip != NULL ) {
   FREE(ip) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   purpose -- to print out a I2OP list
   
   created -- 98feb06, cca
   ----------------------------------
*/
void
I2OP_fprintf ( 
   FILE   *fp, 
   I2OP    *elem 
) {
if ( fp != NULL && elem != NULL ) {
   int   i = 0 ;
   while ( elem != NULL ) {
      if ( i % 4 == 0 ) fprintf(fp, "\n ") ;
      fprintf(fp, " <%4d,%4d,%p>", 
              elem->value0, elem->value1, elem->value2) ;
      elem = elem->next ;
      i++ ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
