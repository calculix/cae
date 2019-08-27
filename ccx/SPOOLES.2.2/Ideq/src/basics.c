/*  basics.c  */

#include "../Ideq.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   create and return a new instance of the Ideq object

   created -- 96jun06, cca
   -----------------------------------------------------
*/
Ideq *
Ideq_new (
   void
) {
Ideq   *deq ;

ALLOCATE(deq, struct _Ideq, 1) ;

Ideq_setDefaultFields(deq) ;

return(deq) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   set the default fields for an Ideq object

   created -- 96jun06, cca
   -------------------------------------------
*/
void
Ideq_setDefaultFields (
   Ideq   *deq
) {
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_setDefaultFields(%p)"
           "\n deq is NULL\n", deq) ;
   exit(-1) ;
}
deq->maxsize =  0 ;
deq->head    = -1 ;
deq->tail    = -1 ;
IV_setDefaultFields(&deq->iv) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data fields

   created -- 96jun06, cca
   -----------------------
*/
void
Ideq_clearData (
   Ideq   *deq
) {
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_clearData(%p)"
           "\n deq is NULL\n", deq) ;
   exit(-1) ;
}
IV_clearData(&deq->iv) ;
Ideq_setDefaultFields(deq) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   free the Ideq object

   created -- 96jun06, cca
   -----------------------
*/
void
Ideq_free (
   Ideq   *deq
) {
if ( deq == NULL ) {
   fprintf(stderr, "\n fatal error in Ideq_free(%p)"
           "\n deq is NULL\n", deq) ;
   exit(-1) ;
}
Ideq_clearData(deq) ;
FREE(deq) ;

return ; }

/*--------------------------------------------------------------------*/
