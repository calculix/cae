/*  findMaxFlow.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   find the maximum flow of a network

   created -- 96jun08, cca
   ----------------------------------
*/
void
Network_findMaxFlow (
   Network   *network
) {
Arc    *arc ;
Arc    **inheads, **outheads ;
FILE   *msgFile ;
Ideq   *deq ;
int    delta, msglvl, nnode, source, tag ;
int    *deltas, *pred, *tags ;
/*
   ---------------
   check the input
   ---------------
*/
if ( network == NULL || (nnode = network->nnode) <= 0 ) {
   fprintf(stderr, "\n fatal error in findMaxFlow(%p)"
           "\n bad input\n", network) ;
   exit(-1) ;
}
outheads = network->outheads ;
inheads  = network->inheads  ;
msglvl   = network->msglvl   ;
msgFile  = network->msgFile  ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n findMaxFlow :\n") ;
}
/*
   ----------------------------------
   set up the working data structures
   ----------------------------------
*/
deq = Ideq_new() ;
Ideq_resize(deq, nnode) ;
pred   = IVinit(nnode, -1) ;
tags   = IVinit(nnode, -1) ;
deltas = IVinit(nnode,  0) ;
/*
   -----------------
   find the max flow
   -----------------
*/
tag    = 0 ;
source = 0 ;
for ( arc = outheads[source] ; arc != NULL ; arc = arc->nextOut ) {
   network->ntrav++ ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n checking out node %d", arc->second) ;
   }
   while ( (delta = arc->capacity - arc->flow) > 0 ) {
      delta = Network_findAugmentingPath(network, arc->second, delta, 
                                         tag, deq, tags, deltas, pred) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n    delta = %d from findAugmentPath(%d)",
                 delta, arc->second) ;
      }
      if ( delta == 0 ) {
          break ;
       }
       Network_augmentPath(network, delta, pred) ;
       tag++ ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Ideq_free(deq) ;
IVfree(pred)   ;
IVfree(tags)   ;
IVfree(deltas) ;
 
return ; }

/*--------------------------------------------------------------------*/
