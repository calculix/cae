/*  addArc.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------
   add an arc to the network

   created -- 96jun08, cca
   -------------------------
*/
void
Network_addArc (
   Network   *network,
   int       firstNode,
   int       secondNode,
   int       capacity,
   int       flow
) {
Arc        *arc ;
Arc        **inheads, **outheads ;
ArcChunk   *chunk ;
int        nnode ;
/*
   ---------------
   check the input
   ---------------
*/
if (  network == NULL 
   || (nnode = network->nnode) <= 0
   || firstNode  < 0 || firstNode  >= nnode
   || secondNode < 0 || secondNode >= nnode
   || capacity <= 0
   || flow < 0 || flow > capacity ) {
   fprintf(stderr, "\n fatal error in Network_addArc(%p,%d,%d,%d,%d)"
           "\n bad input\n", network, firstNode, secondNode,
           capacity, flow) ;
   if ( network != NULL ) {
      fprintf(stderr, "\n nnode = %d", nnode) ;
   }
   exit(-1) ;
}
inheads  = network->inheads  ;
outheads = network->outheads ;
if ( (chunk = network->chunk) == NULL || chunk->inuse == chunk->size ) {
   ALLOCATE(chunk, struct _ArcChunk, 1) ;
   ALLOCATE(chunk->base, struct _Arc, nnode) ;
   chunk->size    = nnode ;
   chunk->inuse   =   0   ;
   chunk->next    = network->chunk ;
   network->chunk = chunk ;
}
arc = chunk->base + chunk->inuse++ ;
arc->first     = firstNode       ;
arc->second    = secondNode      ;
arc->capacity  = capacity        ;
arc->flow      = flow            ;
arc->nextOut   = outheads[firstNode] ;
outheads[firstNode] = arc ;
arc->nextIn    = inheads[secondNode] ;
inheads[secondNode] = arc ;
network->narc++ ;

return ; }

/*--------------------------------------------------------------------*/
