/*  basics.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   construct a new instance of the Network object

   created -- 96jun08, cca
   --------------------------------------------
*/
Network *
Network_new (
   void
) {
Network   *network ;

ALLOCATE(network, struct _Network, 1) ;

Network_setDefaultFields(network) ;

return(network) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   set the default fields of the Network object

   created  -- 96jun08, cca
   ---------------------------------------------
*/
void
Network_setDefaultFields (
   Network   *network
) {
if ( network == NULL ) {
   fprintf(stderr, "\n fatal error in Network_setDefaultFields(%p)"
           "\n bad input\n", network) ;
   exit(-1) ;
}
network->nnode    =   0  ;
network->narc     =   0  ;
network->ntrav    =   0  ;
network->inheads  = NULL ;
network->outheads = NULL ;
network->chunk    = NULL ;
network->msglvl   =   0  ;
network->msgFile  = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   clear the data fields for a Network object

   created  -- 96jun08, cca
   ---------------------------------------------
*/
void
Network_clearData (
   Network   *network
) {
ArcChunk   *chunk ;

if ( network == NULL ) {
   fprintf(stderr, "\n fatal error in Network_clearData(%p)"
           "\n bad input\n", network) ;
   exit(-1) ;
}
if ( network->inheads != NULL ) {
   FREE(network->inheads) ;
}
if ( network->outheads != NULL ) {
   FREE(network->outheads) ;
}
while ( (chunk = network->chunk) != NULL ) {
   network->chunk = chunk->next ;
   FREE(chunk->base) ;
   FREE(chunk) ;
}
Network_setDefaultFields(network) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------
   free the Network object

   created  -- 96jun08, cca
   ------------------------
*/
void
Network_free (
   Network   *network
) {
if ( network == NULL ) {
   fprintf(stderr, "\n fatal error in Network_free(%p)"
           "\n bad input\n", network) ;
   exit(-1) ;
}
Network_clearData(network) ;
FREE(network) ;

return ; }

/*--------------------------------------------------------------------*/
