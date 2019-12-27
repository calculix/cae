/*  util.c  */

#include "../Ideq.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the dequeue,
  
   created -- 96jun06, cca
   -----------------------
*/
void
Ideq_clear (
   Ideq   *deq
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_clear(%p)"
           "\n bad input\n", deq) ;
   exit(-1) ;
}
deq->head = deq->tail = -1 ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   return the head of the dequeue,
   return -1 if the dequeue is empty
  
   created -- 96jun06, cca
   ---------------------------------
*/
int
Ideq_head (
   Ideq   *deq
) {
int   val ;
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_head(%p)"
           "\n bad input\n", deq) ;
   exit(-1) ;
}
if ( deq->head == -1 ) {
   val = -1 ;
} else {
   val =  deq->iv.vec[deq->head] ;
} 
return(val) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return and remove the head of the dequeue,
   return -1 if the dequeue is empty
  
   created -- 96jun06, cca
   ------------------------------------------
*/
int
Ideq_removeFromHead (
   Ideq   *deq
) {
int   val ;
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_removeFromHead(%p)"
           "\n bad input\n", deq) ;
   exit(-1) ;
}
if ( deq->head == -1 ) {
   val = -1 ;
} else {
   val =  deq->iv.vec[deq->head] ;
   if ( deq->head == deq->tail ) {
      deq->head = deq->tail = -1 ;
   } else if ( deq->head == deq->iv.size - 1 ) {
      deq->head = 0 ;
   } else {
      deq->head++ ;
   }
} 
return(val) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   insert value at head of dequeue
   return value
     1 --> value inserted
    -1 --> no room in dequeue, must resize
  
   created -- 96jun06, cca
   ---------------------------------------
*/
int
Ideq_insertAtHead (
   Ideq   *deq,
   int    val 
) {
int   rc, size ;
int   *ivec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_insertAtHead(%p,%d)"
           "\n bad input\n", deq, val) ;
   exit(-1) ;
}
ivec = deq->iv.vec ;
size = deq->iv.size ;
if ( deq->head == -1 ) {
   ivec[0] = val ;
   deq->head = deq->tail = 0 ;
   rc = 1 ;
} else {
   if ( deq->head == 0 ) {
      if ( deq->tail == size - 1 ) {
         rc = -1 ;
      } else {
         ivec[(deq->head = size - 1)] = val ;
         rc = 1 ;
      }
   } else if ( deq->tail == deq->head - 1 ) {
      rc = -1 ;
   } else {
      ivec[--deq->head] = val ;
      rc = 1 ;
   }
} 
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   return the tail of the dequeue,
   return -1 if the dequeue is empty
  
   created -- 96jun06, cca
   ---------------------------------
*/
int
Ideq_tail (
   Ideq   *deq
) {
int   val ;
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_tail(%p)"
           "\n bad input\n", deq) ;
   exit(-1) ;
}
if ( deq->tail == -1 ) {
   val = -1 ;
} else {
   val =  deq->iv.vec[deq->tail] ;
} 
return(val) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return and remove the tail of the dequeue,
   return -1 if the dequeue is empty
  
   created -- 96jun06, cca
   ------------------------------------------
*/
int
Ideq_removeFromTail (
   Ideq   *deq
) {
int   val ;
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_removeFromTail(%p)"
           "\n bad input\n", deq) ;
   exit(-1) ;
}
if ( deq->tail == -1 ) {
   val = -1 ;
} else {
   val =  deq->iv.vec[deq->tail] ;
   if ( deq->head == deq->tail ) {
      deq->head = deq->tail = -1 ;
   } else if ( deq->tail == 0 ) {
      deq->tail = deq->iv.size - 1 ;
   } else {
      deq->tail-- ;
   }
} 
return(val) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   insert value at tail of dequeue
   return value
     1 --> value inserted
    -1 --> no room in dequeue, must resize
  
   created -- 96jun06, cca
   ---------------------------------------
*/
int
Ideq_insertAtTail (
   Ideq   *deq,
   int    val 
) {
int   rc, size ;
int   *ivec ;
/*
   ---------------
   check the input
   ---------------
*/
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_insertAtTail(%p,%d)"
           "\n bad input\n", deq, val) ;
   exit(-1) ;
}
ivec = deq->iv.vec ;
size = deq->iv.size ;
if ( deq->tail == -1 ) {
   ivec[0] = val ;
   deq->head = deq->tail = 0 ;
   rc = 1 ;
} else {
   if ( deq->tail == size - 1 ) {
      if ( deq->head == 0 ) {
         rc = -1 ;
      } else {
         ivec[(deq->tail = 0)] = val ;
         rc = 1 ;
      }
   } else if ( deq->tail + 1 == deq->head ) {
      rc = -1 ;
   } else {
      ivec[++deq->tail] = val ;
      rc = 1 ;
   }
} 
return(rc) ; }

/*--------------------------------------------------------------------*/
