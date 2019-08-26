/*  Graph.h  */

#include "../IVL.h"
#include "../IV.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   The Graph object represents one of three types of graphs,
   defined by its "type" field.

   The Graph object uses IVL objects to store the vertex adjacency
   lists and edge weight lists (for types 2 and 3).

   The Graph object also supports graphs with "boundaries".
   This need arose from ordering graphs using constrained
   minimum degree where boundary vertices contribute to the
   degree but are not eliminated.

   created : 95-sep27, cca
   ---------------------------------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   data fields

   type -- type of graph
      type   vertices weighted?   edges weighted?
        0           no                  no
        1          yes                  no
        2           no                 yes
        3          yes                 yes
   nvtx     -- number of vertices, also equals adjIVL->nlist
               and ewghtIVL->nlist
   nvbnd    -- number of boundary vertices
   nedges   -- number of edges, also equals adjIVL->tsize
               and ewghtIVL->tsize
   totvwght -- total vertex weight
   totewght -- total edge weight
   adjIVL   -- IVL object that holds the adjacency lists
   vwghts   -- int vector that holds the vertex weights, 
               size nvtx + nvbnd, not NULL if type % 2 == 1
   ewghtIVL -- IVL object that holds the edge weight lists,
               not NULL if type >= 2
   ---------------------------------------------------------
*/
typedef struct _Graph   Graph ;
struct _Graph {
   int   type      ;
   int   nvtx      ;
   int   nvbnd     ;
   int   nedges    ;
   int   totvwght  ;
   int   totewght  ;
   IVL   *adjIVL   ;
   int   *vwghts   ;
   IVL   *ewghtIVL ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in basics.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- create and return a new Graph object

   created -- 95sep27, cca
   -----------------------------------------------
*/
Graph *
Graph_new ( 
   void
) ;
/*
   ------------------------------------------------------
   purpose -- set the default fields for the Graph object

   created -- 95sep27, cca
   ------------------------------------------------------
*/
void
Graph_setDefaultFields (
   Graph   *g
) ;
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 95sep27, cca
   --------------------------------
*/
void
Graph_clearData (
   Graph   *g
) ;
/*
   --------------------------------
   purpose -- free the Graph object

   created -- 95sep27, cca
   --------------------------------
*/
void
Graph_free (
   Graph   *g
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in compress.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------
   given a graph g and a fine-to-coarse map vector cmap[],
   return a compressed graph with type coarseType.
   note, the compressed graph will have no trivial boundary
         even if the original graph did have a boundary.

   created -- 95sep29, cca
   ---------------------------------------------------------------
*/
Graph *
Graph_compress (
   Graph   *g,
   int     cmap[],
   int     coarseType
) ;
/*
   ---------------------------------------------------------------
   given a graph g and a fine-to-coarse map vector *mapIV,
   return a compressed graph with type coarseType.
   note, the compressed graph will have no trivial boundary
         even if the original graph did have a boundary.
 
   created -- 96mar02, cca
   ---------------------------------------------------------------
*/
Graph *
Graph_compress2 (
   Graph   *g,
   IV      *mapIV,
   int     coarseType
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in equivMap.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   fill an IV object with an equivalence map
   if ( map[u] == map[v] ) then
      then u and v are adjacent to the same vertices
   endif
   NOTE : each empty list is mapped to a different value

   return value -- IV object that contains the equivalence map

   created -- 96feb23, cca
   -----------------------------------------------------------
*/
IV *
Graph_equivMap (
   Graph   *g
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in expand.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------------
   purpose -- take a Graph object and a map to expand it, create and
              return a bigger unit weight Graph object. this is useful 
              for expanding a compressed graph into a unit weight graph.
 
   created -- 96mar02, cca
   ---------------------------------------------------------------------
*/
Graph *
Graph_expand (
   Graph   *g,
   int     nvtxbig,
   int     map[]
) ;
/*
   ---------------------------------------------------------------------
   purpose -- take a Graph object and a map to expand it, create and
              return a bigger unit weight Graph object. this is useful 
              for expanding a compressed graph into a unit weight graph.
 
   created -- 96mar02, cca
   ---------------------------------------------------------------------
*/
Graph *
Graph_expand2 (
   Graph   *g,
   IV      *mapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in fillFromOffsets.c ----------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   purpose -- take an adjacency structure in the
              (offsets[neqns+1], adjncy[*]) form
              and load the Graph object

   g -- pointer to Graph object, must be initialized with nvtx = neqns
   neqns -- # of equations
   offsets -- offsets vector
   adjncy  -- big adjacency vector
      note, the adjacency for list v is found in
            adjncy[offsets[v]:offsets[v+1]-1]
      also note, offsets[] and adjncy[] must be zero based,
      if (offsets,adjncy) come from a harwell-boeing file, they use
      the fortran numbering, so each value must be decremented to
      conform with C's zero based numbering
   flag -- task flag
      flag = 0 --> just set the adjacency list for v to be that
                   found in adjncy[offsets[v]:offsets[v+1]-1]
      flag = 1 --> the input adjancency is just the upper triangle
                   (or strict upper triangle) as from a harwell-boeing
                   file. fill the Graph object with the full adjacency 
                   structure, including (v,v) edges

   created -- 96mar16, cca
   -------------------------------------------------------------------
*/
void
Graph_fillFromOffsets (
   Graph   *g,
   int     neqns,
   int     offsets[],
   int     adjncy[],
   int     flag
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in init.c ---------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------
   basic initializer for the Graph object

   type      -- graph type
   nvtx      -- # of vertices
   nvbnd     -- # of boundary vertices
   nedges    -- # of edges
   adjType   -- IVL type for adjacency object
   ewghtType -- IVL type for edge weights object

   created -- 95sep27, cca
   ---------------------------------------------
*/
void
Graph_init1 (
   Graph   *g,
   int     type,
   int     nvtx,
   int     nvbnd,
   int     nedges,
   int     adjType,
   int     ewghtType
) ;
/*
   --------------------------------------------------------
   second initializer for the Graph object.
   this function is used in the I/O routines
   Graph_readFromFormattedFile(Graph *g, FILE *fp) and
   Graph_readFromBinaryFile(Graph *g, FILE *fp) where
   the IVL object(s) and vwghts[] vector are created
   independently.

   type     -- graph type
   nvtx     -- # of vertices
   nvbnd    -- # of boundary vertices
   nedges   -- # of edges
   totvwght -- total vertex weight
   totewght -- total edge weight
   adjIVL   -- IVL object for adjacency structure
   vwghts   -- pointer to integer vector for vertex weights
   ewghtIVL -- IVL object for edge weights 

   created -- 95sep27, cca
   --------------------------------------------------------
*/
void
Graph_init2 (
   Graph   *g,
   int     type,
   int     nvtx,
   int     nvbnd,
   int     nedges,
   int     totvwght,
   int     totewght,
   IVL     *adjIVL,
   int     *vwghts,
   IVL     *ewghtIVL
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in util.c ---------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------------
   return the external degree (in terms of vertex weight) of vertex v

   created -- 95oct05, cca
   ---------------------------------------------------------------------
*/
int
Graph_externalDegree (
      Graph   *g,
      int     v
) ;
/*
   -----------------------------------
   method to access the adjacency list

   created -- 95oct05, cca
   -----------------------------------
*/
void
Graph_adjAndSize (
   Graph   *g,
   int     jvtx,
   int     *psize,
   int     **padj
) ;
/*
   -----------------------------------------------
   method to access the adjacency list 
   and possibly the edge weight list for a vertex.

   created -- 95sep29, cca
   -----------------------------------------------
*/
void
Graph_adjAndEweights (
   Graph   *g,
   int     jvtx,
   int     *psize,
   int     **padj,
   int     **pewghts
) ;
/*
   ----------------------------------------------------
   return the number of bytes taken by the Graph object

   created -- 95oct05, cca
   ----------------------------------------------------
*/
int
Graph_sizeOf (
   Graph   *g
) ;
/*
   --------------------------------------
   create and return an IV object filled 
   with a map from vertices to components

   created -- 96feb25, cca
   --------------------------------------
*/
IV *
Graph_componentMap (
   Graph   *g
) ;
/*
   -----------------------------------------------------------------
   given a Graph g and map from vertices to components,
   fill counts[icomp] with the number of vertices in component icomp
   and fill weight[icomp] with their weight

   created -- 96feb25, cca
   -----------------------------------------------------------------
*/
void
Graph_componentStats (
   Graph   *g,
   int     map[],
   int     counts[],
   int     weights[]
) ;
/*
   -------------------------------------------------------------------
   create and return a subgraph and a map from 
   its vertices to the vertices of the graph.

   g       -- graph from which to extract the subgraph
   icomp   -- component from which comes the vertices of the subgraph,
              icomp > 0
   compids -- component ids vector of graph
   pmap    -- pointer to hold address of map vector, the map from
              the subgraph's vertices to the graph's vertices

   return value -- pointer to subgraph Graph object

   created -- 95nov10, cca
   -------------------------------------------------------------------
*/
Graph *
Graph_subGraph (
   Graph   *g,
   int     icomp,
   int     compids[],
   int     **pmap
) ;
/*
   ----------------------------------
   return 1 if the graph is symmetric
   return 0 otherwise
 
   created -- 96oct31, cca
   ----------------------------------
*/
int
Graph_isSymmetric (
   Graph   *graph 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in polar.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   perform a breadth-first search
 
   v        -- vertex to check out
   levelsIV -- IV object that holds the levels vector
   listIV   -- IV object to hold order of search
   parIV    -- IV object to hold parent vector
 
   return value -- depth of the search tree
 
   created -- 96jul13, cca
   -----------------------------------------------------------
*/
int
Graph_BFS (
   Graph   *graph,
   int     v,
   IV      *levelsIV,
   IV      *listIV,
   IV      *parIV
) ;
/*
   ------------------------------------------------------------
   perform a breadth-first search from a list of seed vertices.
 
   seedsIV  -- IV object that holds the seed vertices
   levelsIV -- IV object that holds the levels vector
   listIV   -- IV object to hold order of search
   parIV    -- IV object to hold parent vector
 
   return value -- depth of the search tree
 
   created -- 96jul19, cca
   ------------------------------------------------------------
*/
int
Graph_BFS2 (
   Graph   *graph,
   IV      *seedsIV,
   IV      *levelsIV,
   IV      *listIV,
   IV      *parIV
) ;
/*
   ----------------------------------------------------------
   find a pseudodiameter of a graph.
   fill *pstart and *pend with the endpoint vertices.
   fill startLevelsIV and endLevelsIV with the level vectors
   for the BFS from the start and end vertices, respectively.
 
   created -- 96jul13, cca
   ----------------------------------------------------------
*/
void
Graph_pseudodiameter (
   Graph   *g,
   int     root,
   IV      *startLevelsIV,
   IV      *endLevelsIV,
   int     *pstart,
   int     *pend
) ;
/*
   ----------------------------------------------------------
   create and return an IV object that contains a list of
   vertices that are in the major axis of the graph.
   we have dropped two levels structures from end vertices
   to approximate the diameter of the graph.
   a vertex v is in the major axis if
      levels1[v] + levels2[v] - tolerance <= d(s,e)
   fill *pstart and *pend with the endpoint vertices.
 
   created -- 96jul13, cca
   ----------------------------------------------------------
*/
IV *
Graph_majorAxis (
   Graph   *graph,
   int     vstart,
   int     vend,
   IV      *startLevelsIV,
   IV      *endLevelsIV,
   int     tolerance
) ;
/*
   -----------------------------------------------------------
   compute and return the polar moment of inertia for a vertex
 
   v        -- vertex to check out
   levelsIV -- IV object that holds the levels vector
   listIV   -- IV object for workspace
 
   created -- 96jul13, cca
   -----------------------------------------------------------
*/
double
Graph_polarMoment (
   Graph   *graph,
   int     v,
   IV      *levelsIV,
   IV      *listIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in setListsFromOffsets.c ------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   purpose -- 

   take an adjacency structure in the (offsets[neqns+1], adjncy[*]) 
   form and load the Graph object. note, pointers to the lists are
   set, no new storage is allocated for the adjacency lists.

   however, during the ordering process each adjacency lists 
   may be shuffled.

   g       -- pointer to Graph object, 
              must be initialized with nvtx = neqns
   neqns   -- # of equations
   offsets -- offsets vector
   adjncy  -- big adjacency vector
      note, the adjacency for list v is found in
            adjncy[offsets[v]:offsets[v+1]-1]
      also note, offsets[] and adjncy[] must be zero based,
      if (offsets,adjncy) come from a harwell-boeing file, they use
      the fortran numbering, so each value must be decremented to
      conform with C's zero based numbering

   created -- 96oct24, cca
   -------------------------------------------------------------------
*/
void
Graph_setListsFromOffsets (
   Graph   *g,
   int     neqns,
   int     offsets[],
   int     adjncy[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in wirebasket.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   on input
      stages[v] = 0 --> vertex is in a domain
      stages[v] > 0 --> vertex is in the multisector
   this remains true on output, but the stage of a multisector
   vertex is equal to the number of different domains that
   are found within radius edges from itself
   
   created -- 97jul30, cca
   -----------------------------------------------------------
*/
void
Graph_wirebasketStages (
   Graph   *graph,
   IV      *stagesIV,
   int     radius
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in IO.c -----------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   purpose -- to read in a Graph object from a file

   input --

      fn -- filename, must be *.graphb or *.graphf

   return value -- 1 if success, 0 if failure

   created -- 95sep29, cca
   -------------------------------------------------
*/
int
Graph_readFromFile ( 
   Graph   *graph, 
   char    *fn 
) ;
/*
   -------------------------------------------------------
   purpose -- to read in a Graph object from a CHACO file
 
   input --
 
      fn -- filename
 
   return value -- 1 if success, 0 if failure
 
   created -- 98sep20, jjs
   --------------------------------------------------------
*/
int
Graph_readFromChacoFile (
   Graph   *graph,
   char    *fn
) ;
/*
   --------------------------------------------------------
   purpose -- to read a Graph object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95sep29, cca
   --------------------------------------------------------
*/
int
Graph_readFromFormattedFile ( 
   Graph   *graph, 
   FILE    *fp 
) ;
/*
   ----------------------------------------------------
   purpose -- to read a Graph object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95sep29, cca
   ----------------------------------------------------
*/
int
Graph_readFromBinaryFile ( 
   Graph   *graph, 
   FILE    *fp 
) ;
/*
   --------------------------------------------
   purpose -- to write a Graph object to a file

   input --

      fn -- filename
        *.graphb -- binary
        *.graphf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   --------------------------------------------
*/
int
Graph_writeToFile ( 
   Graph   *graph, 
   char    *fn 
) ;
/*
   ------------------------------------------------------
   purpose -- to write a Graph object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   ------------------------------------------------------
*/
int
Graph_writeToFormattedFile ( 
   Graph   *graph, 
   FILE    *fp 
) ;
/*
   ---------------------------------------------------
   purpose -- to write a Graph object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   ---------------------------------------------------
*/
int
Graph_writeToBinaryFile ( 
   Graph    *graph, 
   FILE   *fp 
) ;
/*
   -------------------------------------------------
   purpose -- to write a Graph object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   -------------------------------------------------
*/
int
Graph_writeForHumanEye ( 
   Graph    *graph, 
   FILE   *fp 
) ;
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the Graph object

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   -----------------------------------------------------------
*/
int
Graph_writeStats ( 
   Graph    *graph, 
   FILE   *fp 
) ;
/*
   ---------------------------------
   write out a graph to a METIS file

   created -- 95oct18, cca
   ---------------------------------
*/
int
Graph_writeToMetisFile (
   Graph   *g,
   FILE    *fp
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----   functions in getTree.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   purpose -- to extract a maximal tree from a graph
 
   on input,
      if mark[v] == 'N' then
         v is eligible to be placed in the tree
      endif
      i.e., mark[] serves as a mask vector
 
   on return,
      if ( mark[v] == 'Y' ) then
         v is in the tree, parent of v = par[v]
      else
         mark[v] = 'N', not in tree
      endif
 
   created -- 98apr10, cca
   -------------------------------------------------
*/
void
Graph_getTree (
   Graph   *graph,
   int     par[],
   char    mark[],
   int     seed,
   int     msglvl,
   FILE    *msgFile
) ;
/*
   --------------------------------------------------------
   purpose -- to extract a maximal nlevel-tree from a graph
      (a vertex v can have an edge in the graph with itself
      and any ancestor within nlevel levels in the tree.)
 
   on input, 
      if mark[v] == 'N' then 
         v is eligible to be placed in the tree
      endif
      i.e., mark[] serves as a mask vector
 
   on return,
      if ( mark[v] == 'Y' ) then
         v is in the tree, parent of v = par[v]
      else
         mark[v] = 'N', not in tree
      endif
 
   created -- 98apr10, cca
   --------------------------------------------------------
*/
void
Graph_getTree2 (
   Graph   *graph,
   int     nlevel,
   int     par[],
   char    mark[],
   int     seed,
   int     msglvl,
   FILE    *msgFile
) ;
/*--------------------------------------------------------------------*/
