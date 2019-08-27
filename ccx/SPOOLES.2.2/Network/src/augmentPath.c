/*  augmentPath.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   given an augmenting path defined by the pred[] vector and delta,
   the change in flow, augment the path from the source to the sink

   created -- 86jun08, cca
   ----------------------------------------------------------------
*/
void
Network_augmentPath (
   Network   *network,
   int       delta,
   int       pred[]
) {
Arc    *arc ;
Arc    **inheads, **outheads ;
FILE   *msgFile ;
int    msglvl, nnode, sink, source, v, w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( network == NULL || (nnode = network->nnode) <= 0
   || delta <= 0 || pred == NULL ) {
   fprintf(stderr, "\n fatal error in Network_augmentPath(%p,%d,%p)"
           "\n bad input\n",
           network, delta, pred) ;
   exit(-1) ;
}
source   = 0 ;
sink     = nnode - 1 ;
inheads  = network->inheads  ;
outheads = network->outheads ;
msglvl   = network->msglvl   ;
msgFile  = network->msgFile  ;
/*
   -----------------------
   work back from the sink
   -----------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n augment path") ;
   fflush(msgFile) ;
}
w = sink ;
while ( 1 ) {
   v = pred[w] ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n    w = %d, v = %d", w, v) ;
   }
   for ( arc = inheads[v] ; arc != NULL ; arc = arc->nextIn ) {
      network->ntrav++ ;
      if ( arc->first == w ) {
         arc->flow -= delta ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n   backward arc(%d,%d), flow = %d",
                    w, v, arc->flow) ;
         }
         break ;
      }
   }
   if ( arc == NULL ) {
      for ( arc = outheads[v] ; arc != NULL ; arc = arc->nextOut ) {
         network->ntrav++ ;
         if ( arc->second == w ) {
            arc->flow += delta ;
            if ( msglvl > 2 ) {
               fprintf(msgFile, "\n   forward arc(%d,%d), flow = %d",
                    v, w, arc->flow) ;
            }
            break ;
         }
      }
   }
   if ( v == source ) {
      break ;
   }
   w = v ;
}
return ; }

/*--------------------------------------------------------------------*/
