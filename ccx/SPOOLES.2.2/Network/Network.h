/*  Network.h  */

#include "../cfiles.h"
#include "../Ideq.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   forward declarations of typedef's
   ---------------------------------
*/
typedef struct _Network   Network ;
typedef struct _Arc   Arc ;
typedef struct _ArcChunk   ArcChunk ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   structure to hold a network.

   nnode    -- number of nodes in the network.
               the source node is 0, the sink node is nnode - 1
   narc     -- number of arcs in the network
   ntrav    -- number of arcs that are traversed,
               used as a measure of complexity
   inheads  -- array of pointers used to point to the first arc 
               in the in list for the nodes
   outheads -- array of pointers used to point to the first arc 
               in the out list for the nodes
   chunk    -- pointer to the first ArcChunk structure
   msglvl   -- message level, default is 0
   msgFile  -- message file, default is stdout
   ------------------------------------------------------------
*/
struct _Network {
   int        nnode      ;
   int        narc       ;
   int        ntrav      ;
   Arc        **inheads  ;
   Arc        **outheads ;
   ArcChunk   *chunk     ;
   int        msglvl     ;
   FILE       *msgFile   ;
} ;
/*
   ----------------------------------------------------------------
   arc structure

   vertices in the arc are (first, second). there is one arc that
   is linked into the out list for vertex first and the in list for
   vertex second via the nextOut and nextIn fields, respectively.
   ----------------------------------------------------------------
*/
struct _Arc {
   int   first    ;
   int   second   ;
   int   capacity ;
   int   flow     ;
   Arc   *nextOut ;
   Arc   *nextIn  ;
} ;
/*
   ---------------------------------------------------------------
   structure to hold a vector of arcs, used in forming the network
   when the total number of arcs is not known ahead of time.

   size  -- dimension of base[] array
   inuse -- number of arc structures in use, 
            next free arc is located at &base[inuse]
   base  -- vector of Arc structures
   next  -- link to the next ArcChunk structure,
            used to free the storage when finished
   ---------------------------------------------------------------
*/
struct _ArcChunk {
   int        size  ;
   int        inuse ;
   Arc        *base ;
   ArcChunk   *next ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c   --------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------
   construct a new instance of the Network object

   created -- 96jun08, cca
   --------------------------------------------
*/
Network *
Network_new (
   void
) ;
/*
   ---------------------------------------------
   set the default fields of the Network object

   created  -- 96jun08, cca
   ---------------------------------------------
*/
void
Network_setDefaultFields (
   Network   *network
) ;
/*
   ---------------------------------------------
   clear the data fields for a Network object

   created  -- 96jun08, cca
   ---------------------------------------------
*/
void
Network_clearData (
   Network   *network
) ;
/*
   ------------------------
   free the Network object

   created  -- 96jun08, cca
   ------------------------
*/
void
Network_free (
   Network   *network
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c   ----------------------------------------
------------------------------------------------------------------------
*/
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in addArc.c   --------------------------------------
------------------------------------------------------------------------
*/
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in findAugmentingPath.c   --------------------------
------------------------------------------------------------------------
*/
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in augmentingPath.c   ------------------------------
------------------------------------------------------------------------
*/
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in findMaxFlow.c   ---------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------
   find the maximum flow of a network
 
   created -- 96jun08, cca
   ----------------------------------
*/
void
Network_findMaxFlow (
   Network   *network
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in findMincut.c   ----------------------------------
------------------------------------------------------------------------
*/
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c   ------------------------------------------
------------------------------------------------------------------------
*/
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
) ;
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
) ;
/*--------------------------------------------------------------------*/
