/*  findAugmentingPath.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   find an augmenting path from the source to the sink through node.

   node   -- (source, node) is an arc below capacity
   delta  -- capacity(source, node) - flow(source, node)
   tag    -- tag used to build this traversal tree
   deq    -- dequeue object to handle traversal,
             out-nodes get put at head of list,
             in-nodes get put at tail of list
   tags   -- tags vector for the nodes
   deltas -- flow increment vector for the nodes
   pred   -- predecessor vector for the nodes

   created -- 96jun08, cca
   -----------------------------------------------------------------
*/
int
Network_findAugmentingPath (
   Network   *network,
   int       node,
   int       delta,
   int       tag,
   Ideq      *deq,
   int       tags[],
   int       deltas[],
   int       pred[]
) {
Arc    *arc ;
Arc    **inheads, **outheads ;
FILE   *msgFile ;
int    msglvl, nnode, resid, sink, source, v, w ;
/*
   ---------------
   check the input
   ---------------
*/ 
if (  network == NULL || (nnode = network->nnode) <= 0
   || node <= 0 || node >= nnode -1
   || deq == NULL || tags == NULL || deltas == NULL || pred == NULL ) {
   fprintf(stderr, 
"\n fatal error in Network_findAugmentingPath(%p,%d,%d,%d,%p,%p,%p,%p)"
"\n bad input\n",
           network, node, delta, tag, deq, tags, deltas, pred) ;
   exit(-1) ;
}
inheads  = network->inheads  ;
outheads = network->outheads ;
msglvl   = network->msglvl   ;
msgFile  = network->msgFile  ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n try to find augmenting path starting at node %d",
           node) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   load the dequeue with the seed node and
   set the tags, pred and delta fields
   ---------------------------------------
*/
Ideq_clear(deq) ;
Ideq_insertAtHead(deq, node) ;
source = 0 ;
sink   = nnode - 1 ;
tags[source] = tags[node] = tag ;
deltas[node] = delta ;
pred[node]   = source ;
/*
   ------------------------------------------
   try to find an augmenting path to the sink
   ------------------------------------------
*/
while ( tags[sink] != tag ) {
   v = Ideq_removeFromHead(deq) ;
   if ( v < 0 ) {
/*
      -----------------------------------
      dequeue is empty, break out of loop
      -----------------------------------
*/
      break ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n node %d removed from head of dequeue", v) ;
   }
/*
   ----------------------
   check out the out-arcs
   ----------------------
*/
   for ( arc = outheads[v] ; arc != NULL ; arc = arc->nextOut ) {
      network->ntrav++ ;
      w = arc->second ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n    out-arc (%d,%d), flow %d, capacity %d",
                 arc->first, arc->second, arc->flow, arc->capacity) ;
      }
      if ( tags[w] != tag && (resid = arc->capacity - arc->flow) > 0 ) {
         if ( resid > deltas[v] ) {
            resid = deltas[v] ;
         }
         deltas[w] = resid ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
                    ", now a tree arc, delta = %d", deltas[w]) ;
         }
         tags[w] = tag ;
         pred[w] =  v  ;
         if ( w == sink ) {
            return(resid) ;
         }
         Ideq_insertAtHead(deq, w) ;
      }
   }
/*
   ---------------------
   check out the in-arcs
   ---------------------
*/
   for ( arc = inheads[v] ; arc != NULL ; arc = arc->nextIn ) {
      network->ntrav++ ;
      w = arc->first ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n    in-arc (%d,%d), flow %d, capacity %d",
                arc->first, arc->second, arc->flow, arc->capacity) ;
      }
      if ( tags[w] != tag && (resid = arc->flow) > 0 ) {
         if ( resid > deltas[v] ) {
            resid = deltas[v] ;
         }
         deltas[w] = resid ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
                    ", now a tree arc, delta = %d", deltas[w]) ;
         }
         tags[w] = tag ;
         pred[w] =  v  ;
         Ideq_insertAtTail(deq, w) ;
      }
   }
}
/*
   ----------------------------------------------------
   if we got to here the sink was not reached, return 0
   ----------------------------------------------------
*/
return(0) ; }

/*--------------------------------------------------------------------*/
