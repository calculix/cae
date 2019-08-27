/*  findMincut.c  */

#include "../Network.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   mark the nodes 1 or 2 to define a min-cut,
   where source in X and sink in Y.
   start the search from the source node.
   for x in X and y in Y 
      arc (x,y) is in the min-cut when flow(x,y) == capacity(x,y)
      arc (y,x) is in the min-cut when flow(y,x) == 0
   on return, mark[*] is filled with 1 or 2, 
   where the mark[source] = 1 and mark[sink] = 2

   created -- 96jun08, cca
   --------------------------------------------------------------
*/
void
Network_findMincutFromSource (
   Network   *network,
   Ideq      *deq,
   int       mark[]
) {
Arc    *arc ;
Arc    **inheads, **outheads ;
FILE   *msgFile ;
int    msglvl, nnode, source, x, z ;
/*
   ---------------
   check the input
   ---------------
*/
if (  network == NULL || (nnode = network->nnode) <= 0 
   || deq == NULL || mark == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Network_findMincutFromSource(%p,%p,%p)"
           "\n bad input\n", network, deq, mark) ;
   exit(-1) ;
}
source   = 0 ;
inheads  = network->inheads  ;
outheads = network->outheads ;
msglvl   = network->msglvl   ;
msgFile  = network->msgFile  ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n Network_findMincutFromSource") ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------
   load all the nodes into Y except for the source
   -----------------------------------------------
*/
IVfill(nnode, mark, 2) ;
mark[source] = 1 ;
/*
   ---------------------------------------------------------
   do a breadth first traversal from the source
   visit x in X
      out edge (x,z), add z to X if flow(x,z) < capacity(x,z)
      in edge (z,x), add z to X if flow(z,x) > 0
   ---------------------------------------------------------
*/
Ideq_clear(deq) ;
Ideq_insertAtHead(deq, source) ;
while ( (x = Ideq_removeFromHead(deq)) != -1 ) {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n checking out node %d", x) ;
      fflush(msgFile) ;
   }
   for ( arc = outheads[x] ; arc != NULL ; arc = arc->nextOut ) {
      z = arc->second ;
      if ( mark[z] != 1 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile,
                    "\n    out-arc (%d,%d), flow %d, capacity %d",
                    x, z, arc->flow, arc->capacity) ;
            fflush(msgFile) ;
         }
         if ( arc->flow < arc->capacity ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, ", adding %d to X", z) ;
               fflush(msgFile) ;
            }
            Ideq_insertAtTail(deq, z) ;
            mark[z] = 1 ;
         }
      }
   }
   for ( arc = inheads[x] ; arc != NULL ; arc = arc->nextIn ) {
      z = arc->first ;
      if ( mark[z] != 1 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n    in-arc (%d,%d), flow %d",
                    z, x, arc->flow) ;
            fflush(msgFile) ;
         }
         if ( arc->flow > 0 ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, ", adding %d to X", z) ;
               fflush(msgFile) ;
            }
            Ideq_insertAtTail(deq, z) ;
            mark[z] = 1 ;
         }
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n leaving FindMincutFromSource") ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   mark the nodes 1 or 2 to define a min-cut,
   where source in X and sink in Y.
   start the search from the sink node.
   for x in X and y in Y 
      arc (x,y) is in the min-cut when flow(x,y) == capacity(x,y)
      arc (y,x) is in the min-cut when flow(y,x) == 0
   on return, mark[*] is filled with 1 or 2, 
   where the mark[source] = 1 and mark[sink] = 2

   created -- 96jun08, cca
   --------------------------------------------------------------
*/
void
Network_findMincutFromSink (
   Network   *network,
   Ideq      *deq,
   int       mark[]
) {
Arc    *arc ;
Arc    **inheads, **outheads ;
FILE   *msgFile ;
int    msglvl, nnode, sink, x, z ;
/*
   ---------------
   check the input
   ---------------
*/
if (  network == NULL || (nnode = network->nnode) <= 0 
   || deq == NULL || mark == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Network_findMincutFromSink(%p,%p,%p)"
           "\n bad input\n", network, deq, mark) ;
   exit(-1) ;
}
sink     = nnode - 1 ;
inheads  = network->inheads  ;
outheads = network->outheads ;
msglvl   = network->msglvl   ;
msgFile  = network->msgFile  ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n Network_findMincutFromSink") ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   load all the nodes into X except for the sink
   ---------------------------------------------
*/
IVfill(nnode, mark, 1) ;
mark[sink] = 2 ;
/*
   ---------------------------------------------------------
   do a breadth first traversal from the sink
   visit y in Y
      out edge (x,z), add z to X if flow(x,z) < capacity(x,z)
      in edge (z,x), add z to X if flow(z,x) > 0
   ---------------------------------------------------------
*/
Ideq_clear(deq) ;
Ideq_insertAtHead(deq, sink) ;
while ( (x = Ideq_removeFromHead(deq)) != -1 ) {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n checking out node %d", x) ;
      fflush(msgFile) ;
   }
   for ( arc = outheads[x] ; arc != NULL ; arc = arc->nextOut ) {
      z = arc->second ;
      if ( mark[z] != 2 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile,
                    "\n    out-arc (%d,%d), flow %d, capacity %d",
                    x, z, arc->flow, arc->capacity) ;
            fflush(msgFile) ;
         }
         if ( arc->flow > 0 ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, ", adding %d to X", z) ;
               fflush(msgFile) ;
            }
            Ideq_insertAtTail(deq, z) ;
            mark[z] = 2 ;
         }
      }
   }
   for ( arc = inheads[x] ; arc != NULL ; arc = arc->nextIn ) {
      z = arc->first ;
      if ( mark[z] != 2 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n    in-arc (%d,%d), flow %d",
                    z, x, arc->flow) ;
            fflush(msgFile) ;
         }
         if ( arc->flow < arc->capacity ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, ", adding %d to X", z) ;
               fflush(msgFile) ;
            }
            Ideq_insertAtTail(deq, z) ;
            mark[z] = 2 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
