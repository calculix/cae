/*  init.c  */

#include "../A2.h"

/*====================================================================*/
/*
   ------------------------------------------------------------------
   initializer. sets the n1, n2, inc1 and inc2 fields.
   must have
      mtx != NULL
      type = SPOOLES_REAL or SPOOLES_COMPLEX
      n1, n2, inc1, inc2 > 0
      (inc1 = 1 and inc2 = nrow) or (inc1 = ncol and inc2 = 1)

   if entries is NULL then
      entries are allocated and zeroed.
   else
      mtx->nowned is set to 0
      mtx->entries is set to entries.
   endif
   n1, n2, inc1, inc2 set

   created -- 98apr15, cca
   ------------------------------------------------------------------
*/
void
A2_init ( 
   A2       *mtx, 
   int      type, 
   int      n1, 
   int      n2, 
   int      inc1, 
   int      inc2, 
   double   *entries 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || n1 <= 0 || n2 <= 0 || inc1 <= 0 || inc2 <= 0 ) {
   fprintf(stderr, "\n fatal error in A2_init(%p,%d,%d,%d,%d,%p)"
           "\n bad input\n", mtx, n1, n2, inc1, inc2, entries) ;
   exit(-1) ;
}
if ( type != SPOOLES_REAL && type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in A2_init(%p,%d,%d,%d,%d,%p)"
           "\n bad type %d\n", 
           mtx, n1, n2, inc1, inc2, entries, type) ;
   exit(-1) ;
}
if ( entries == NULL
   && !( (inc1 == 1 && inc2 == n1) || (inc1 == n2 && inc2 == 1) ) ) {
/*
   ---------------------------------------------------------------
   whoa, when we own the data we can only initialize a full matrix
   ---------------------------------------------------------------
*/
   fprintf(stderr, "\n fatal error in A2_init(%p,%d,%d,%d,%d,%p)"
           "\n entries is not NULL and we have bad increments"
           "\n inc1 = %d, inc2 = %d, nrow = %d, ncol = %d\n", 
           mtx, n1, n2, inc1, inc2, entries, 
           inc1, inc2, n1, n2) ;
   exit(-1) ;
}
if ( entries != NULL ) {
/*
   ----------------------------
   set pointer to these entries
   ----------------------------
*/
   if ( mtx->entries != NULL ) {
      DVfree(mtx->entries) ;
   }
   mtx->nowned  = 0 ;
   mtx->entries = entries ;
} else {
   int   nbytesNeeded, nbytesPresent ;
/*
   ----------------------------
   allocate and own the entries
   ----------------------------
*/
   if ( mtx->type == SPOOLES_REAL ) {
      nbytesPresent = mtx->nowned * sizeof(double) ;
   } else if ( mtx->type == SPOOLES_COMPLEX ) {
      nbytesPresent = 2 * mtx->nowned * sizeof(double) ;
   } else {
      nbytesPresent = 0 ;
   }
   if ( type == SPOOLES_REAL ) {
      nbytesNeeded = n1 * n2 * sizeof(double) ; 
   } else if ( type == SPOOLES_COMPLEX ) { 
      nbytesNeeded = 2 * n1 * n2 * sizeof(double) ; 
   }
   if ( nbytesNeeded > nbytesPresent ) {
      DVfree(mtx->entries) ;
      mtx->nowned  = n1 * n2 ;
      if ( type == SPOOLES_REAL ) {
         mtx->entries = DVinit(n1*n2, 0.0) ;
      } else if ( type == SPOOLES_COMPLEX ) { 
         mtx->entries = DVinit(2*n1*n2, 0.0) ;
      }
   }
}
mtx->type = type ;
mtx->n1   = n1   ;
mtx->n2   = n2   ;
mtx->inc1 = inc1 ;
mtx->inc2 = inc2 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   submatrix initializer 

   A(0:lastrow-firstrow,0:lastcol-firstcol) 
              = B(firstrow:lastrow, firstcol:lastcol)

   created -- 98apr15, cca
   --------------------------------------------------
*/
void
A2_subA2 (
   A2   *mtxA,
   A2   *mtxB,
   int   firstrow,
   int   lastrow,
   int   firstcol,
   int   lastcol
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtxA == NULL || mtxB == NULL 
   || firstrow < 0 || lastrow >= mtxB->n1 
   || firstcol < 0 || lastcol >= mtxB->n2 ) {
   fprintf(stderr, "\n fatal error in A2_subA2(%p,%p,%d,%d,%d,%d)"
           "\n bad input\n", 
           mtxA, mtxB, firstrow, lastrow, firstcol, lastcol) ;
   if ( mtxA != NULL ) {
      fprintf(stderr, "\n first A2") ;
      A2_writeForHumanEye(mtxA, stderr) ;
   }
   if ( mtxB != NULL ) {
      fprintf(stderr, "\n second A2") ;
      A2_writeForHumanEye(mtxB, stderr) ;
   }
   exit(-1) ;
}
if ( ! (A2_IS_REAL(mtxB) || A2_IS_COMPLEX(mtxB)) ) {
   fprintf(stderr, "\n fatal error in A2_subA2(%p,%p,%d,%d,%d,%d)"
           "\n bad type %d\n", mtxA, mtxB, firstrow, lastrow, 
           firstcol, lastcol, mtxB->type) ;
   exit(-1) ;
}
mtxA->type    = mtxB->type ;
mtxA->inc1    = mtxB->inc1 ;
mtxA->inc2    = mtxB->inc2 ;
mtxA->n1      = lastrow - firstrow + 1 ;
mtxA->n2      = lastcol - firstcol + 1 ;
if ( A2_IS_REAL(mtxB) ) {
   mtxA->entries = mtxB->entries 
                 + firstrow*mtxB->inc1 + firstcol*mtxB->inc2 ;
} else if ( A2_IS_COMPLEX(mtxB) ) {
   mtxA->entries = mtxB->entries 
                 + 2*(firstrow*mtxB->inc1 + firstcol*mtxB->inc2) ;
}
mtxA->nowned  = 0 ;

return ; }

/*--------------------------------------------------------------------*/
