/*  IO.c  */

#include "../Ideq.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- write the contents of the dequeue to a file
              in a human readable format

   created -- 98feb11, cca
   ------------------------------------------------------
*/
void
Ideq_writeForHumanEye (
   Ideq   *dequeue,
   FILE   *fp
) {
int   ii ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dequeue == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_writeForHumanEye(%p,%p)"
           "\n bad input\n", dequeue, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n\n Ideq : maxsize = %d, head = %d, tail = %d\n",
        dequeue->maxsize, dequeue->head, dequeue->tail) ;
if ( dequeue->head != -1 && dequeue->tail != -1) {
   if ( dequeue->head <= dequeue->tail ) {
      for ( ii = dequeue->head ; ii <= dequeue->tail ; ii++ ) {
         fprintf(fp, " %d", IV_entry(&dequeue->iv, ii)) ;
      }
   } else {
      for ( ii = dequeue->head ; ii < dequeue->maxsize ; ii++ ) {
         fprintf(fp, " %d", IV_entry(&dequeue->iv, ii)) ;
      }
      for ( ii = 0 ; ii <= dequeue->tail ; ii++ ) {
         fprintf(fp, " %d", IV_entry(&dequeue->iv, ii)) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
