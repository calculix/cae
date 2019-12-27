/*  init.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   initialize the Network object

   nnode -- number of nodes in the network, must be > 2,
            source node is 0, sink node is nnode - 1
   narc  -- number of arcs. this value can be zero if the number 
            of arcs is not known at initialization time. 
            storage for the arcs will grow as needed

   created -- 96jun08, cca
   -------------------------------------------------------------
*/
void
Network_init (
   Network   *network,
   int       nnode,
   int       narc
) {
int   v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( network == NULL || nnode <= 2 || narc < 0 ) {
   fprintf(stderr, "\n fatal error in Network_init(%p,%d,%d)"
           "\n bad input\n", network, nnode, narc) ;
   exit(-1) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Network_clearData(network) ;
network->nnode = nnode ;
/*
   -------------------
   allocate the arrays
   -------------------
*/
ALLOCATE(network->inheads,  struct _Arc *, nnode) ;
ALLOCATE(network->outheads, struct _Arc *, nnode) ;
for ( v = 0 ; v < nnode ; v++ ) {
   network->inheads[v] = network->outheads[v] = NULL ;
}
if ( narc > 0 ) {
   ArcChunk   *chunk ;
/*
   -------------------------------
   allocate the ArcChunk structure
   -------------------------------
*/
   ALLOCATE(chunk, struct _ArcChunk, 1) ;
   ALLOCATE(chunk->base, struct _Arc, narc) ;
   chunk->size = narc ;
   chunk->inuse = 0 ;
   chunk->next = NULL ;
   network->chunk = chunk ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   purpose -- set the message fields

   created -- 96oct23, cca
   ---------------------------------
*/
void
Network_setMessageInfo (
   Network   *network,
   int       msglvl,
   FILE      *msgFile
) {
/*
   --------------
   check the data
   --------------
*/
if ( network == NULL ) {
   fprintf(stderr, "\n fatal error in Network_setMessageInfo(%p,%d,%p)"
           "\n bad input\n", network, msglvl, msgFile) ;
   exit(-1) ;
}
if ( msglvl >= 0 ) {
   network->msglvl = msglvl ;
}
if ( msgFile != NULL ) {
   network->msgFile = msgFile ;
} else {
   network->msgFile = stdout ;
}

return ; }

/*--------------------------------------------------------------------*/
