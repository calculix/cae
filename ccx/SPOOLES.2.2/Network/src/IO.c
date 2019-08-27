/*  IO.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   print the network for debugging purposes

   created -- 96jun08, cca
   ----------------------------------------
*/
void
Network_writeForHumanEye (
   Network   *network,
   FILE      *fp
) {
Arc   *arc ;
Arc   **inheads, **outheads ;
int   jnode, nnode ;

Network_writeStats(network, fp) ;
nnode    = network->nnode    ;
inheads  = network->inheads  ;
outheads = network->outheads ;
for ( jnode = 0 ; jnode < nnode ; jnode++ ) {
   fprintf(fp, "\n in list for %d :", jnode) ;
   fflush(fp) ;
   for ( arc = inheads[jnode] ; arc != NULL ; arc = arc->nextIn ) {
      fprintf(fp, " <%d,%d,%d>", 
              arc->first, arc->flow, arc->capacity) ;
      fflush(fp) ;
   }
   fprintf(fp, "\n out list for %d :", jnode) ;
   fflush(fp) ;
   for ( arc = outheads[jnode] ; arc != NULL ; arc = arc->nextOut ) {
      fprintf(fp, " <%d,%d,%d>", 
              arc->second, arc->flow, arc->capacity) ;
      fflush(fp) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   print the network statistics for debugging purposes

   created -- 96jun08, cca
   ---------------------------------------------------
*/
void
Network_writeStats (
   Network   *network,
   FILE      *fp
) {
fprintf(fp, "\n\n Network : %d nodes, %d arcs, %d arc traversals",
        network->nnode, network->narc, network->ntrav) ;

return ; }

/*--------------------------------------------------------------------*/
