/*  instance.c  */

#include "../IVL.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   return the storage type of the object

   created -- 96dec06, cca
   -------------------------------------
*/
int
IVL_type (
   IVL   *ivl
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_type(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
return(ivl->type) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return the maximum number of lists

   created -- 96dec06, cca
   ----------------------------------
*/
int
IVL_maxnlist (
   IVL   *ivl
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_maxnlist(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
return(ivl->maxnlist) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------
   return the number of lists

   created -- 96dec06, cca
   --------------------------
*/
int
IVL_nlist (
   IVL   *ivl
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_nlist(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
return(ivl->nlist) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return the total size of the lists

   created -- 96dec06, cca
   ----------------------------------
*/
int
IVL_tsize (
   IVL   *ivl
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_tsize(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
return(ivl->tsize) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   return the storage increment

   created -- 96dec06, cca
   ----------------------------
*/
int
IVL_incr (
   IVL   *ivl
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_incr(%p)"
           "\n bad input\n", ivl) ;
   exit(-1) ;
}
return(ivl->incr) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------
   set the storage increment

   created -- 96dec06, cca
   -------------------------
*/
void
IVL_setincr (
   IVL   *ivl,
   int   incr
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || incr < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_setincr(%p,%d)"
           "\n bad input\n", ivl, incr) ;
   exit(-1) ;
}
ivl->incr = incr ;

return ; }

/*--------------------------------------------------------------------*/
