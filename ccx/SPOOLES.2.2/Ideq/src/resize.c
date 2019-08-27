/*  resize.c  */

#include "../Ideq.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   resize the deque
   if the new size is large enough then
      copy the old data
      return 1
   else
      error, return -1
   endif

   created -- 96jun06, cca
   ------------------------------------
*/
int
Ideq_resize (
   Ideq   *deq,
   int    newsize
) {
int   count, head, j, size, tail ;
int   *ivec, *tmp ;
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL || newsize < 0 ) {
   fprintf(stderr, "\n fatal error in Ideq_resize(%p,%d)"
           "\n bad input\n", deq, newsize) ;
   exit(-1) ;
}
/*
   ---------------------------------------
   check that the new size is large enough
   ---------------------------------------
*/
if ( deq->tail >= deq->head ) {
   size = deq->tail - deq->head + 1 ;
} else {
   size = deq->tail + deq->iv.size - deq->head + 1 ;
}
if ( size > newsize ) {
   return(-1) ;
}
/*
   ---------------------
   create the new vector
   ---------------------
*/
tmp = IVinit(size, -1) ;
if ( (j = deq->head) != -1 ) {
   ivec = deq->iv.vec ;
   count = 0 ;
   while ( j != deq->tail ) {
      tmp[count++] = ivec[j] ;
      if ( j == size - 1 ) {
         j = 0 ;
      } else {
         j++ ;
      }
   }
   tmp[count++] = ivec[j] ;
   head = 0 ;
   tail = count - 1 ;
} else {
   head = tail = -1 ;
}
Ideq_clearData(deq) ;
IV_init(&deq->iv, newsize, NULL) ;
if ( size > 0 ) {
   IVcopy(size, deq->iv.vec, tmp) ;
}
IVfree(tmp) ;
deq->head    = head ;
deq->tail    = tail ;
deq->maxsize = newsize ;

return(1) ; }

/*--------------------------------------------------------------------*/
