/*  MSMD.h  */

#include "../cfiles.h"
#include "../Graph.h"
#include "../IIheap.h"
#include "../ETree.h"
#include "../IV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   typedef definitions as forward references
   -----------------------------------------
*/
typedef  struct _MSMD           MSMD          ;
typedef  struct _MSMDinfo       MSMDinfo      ;
typedef  struct _MSMDstageInfo  MSMDstageInfo ;
typedef  struct _MSMDvtx        MSMDvtx       ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   MSMD -- multistage constrained minimum degree object

   nvint    -- # of internal vertices
   nvbnd    -- # of boundary vertices
               nvbnd >  0 --> constrained md used
               nvbnd == 0 --> unconstrained md used 
                              for graph has no boundary
   heap     -- pointer to IIheap object, used as a priority queue
   incrIP   -- increment for new IP structures
   baseIP   -- pointer to base storage for working IP structures
   freeIP   -- pointer to free list of IP structures
   vertices -- array of Vtx objects
   ivtmpIV  -- IV object to hold a working vector
   reachIV  -- IV object to hold the reach set
   --------------------------------------------------------------------
*/
struct _MSMD {
   int      nvtx      ;
   IIheap   *heap     ;
   int      incrIP    ;
   IP       *baseIP   ;
   IP       *freeIP   ;
   MSMDvtx  *vertices ;
   IV       ivtmpIV   ;
   IV       reachIV   ;
} ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   the MSMDinfo object contains information for and about a run of
   the MSMD algorithm.

   user supplied or default information

      compressFlag -- flag for initial compression
         compressFlag / 4 >= 1 --> compress before elimination
         compressFlag % 4 == 2 --> compress at each elimination step,
                                   consider all nodes
         compressFlag % 4 == 1 --> compress at each elimination step,
                                   but only consider 2-adj nodes
         compressFlag % 4 == 0 --> do not perform any compression
         default value = 1
      prioType -- priority type
         prioType == 0 --> zero priority for all vertices,
                           i.e., random and independent set elimination
         prioType == 1 --> true external degree
         prioType == 2 --> approximate external updates
         prioType == 3 --> half and half
         default value = 1
      stepType -- definition of nodes in a step
         when stepType < 1 --> only one node is eliminated at a step,
            e.g., like QMD from SPARSPAK and YSMP
         stepType == 1 --> regular multiple elimination of nodes
            with minimum priority, e.g., GENMMD
         stepType >  1 --> extended multiple elimination
            an independent set of nodes is selected for elimination 
            whose degree satisfies minprio <= prio <= stepType*minprio
         default value = 1
      seed    -- random number seed
      msglvl  -- message level
         default value = 0, no statistics
      msgFile -- message file
         default value = stdout

   storage information

      maxnbytes -- maximum number of bytes
      nbytes    -- present number of bytes

   information available about the ordering 

      istage    -- present stage
      nstage    -- number of stages, supplied by the stages[] vector,
                   an input parameter to the MSMD_order method
      stageInfo -- pointer to vector of stageInfo objects, 
                   nstage + 1 in size
   ---------------------------------------------------------------------
*/
struct _MSMDinfo {
   int             compressFlag ;
   int             prioType     ;
   double          stepType     ;
   int             seed         ;
   int             msglvl       ;
   FILE            *msgFile     ;
   int             maxnbytes    ;
   int             nbytes       ;
   int             istage       ;
   int             nstage       ;
   MSMDstageInfo   *stageInfo   ;
   double          totalCPU     ;
} ;
/*--------------------------------------------------------------------*/
/*
   the MSMDstageInfo structure contains statistics 
   about the elimination at a certain stage

      nstep     -- # of elimination steps
      nfront    -- # of fronts
      welim     -- weight of vertices eliminated
      nfind     -- # of front indices
      nzf       -- number of factor entries
      ops       -- number of factor operations
      nexact2   -- # of 2-adjacent exact updates
      nexact3   -- # of 3-adjacent exact updates
      napprox   -- # of approximate degree updates
      ncheck    -- number of indistinguishable node checks
      nindst    -- number of indistinguishable nodes found
      noutmtch  -- number of outmatched nodes found
*/
struct _MSMDstageInfo {
   int      nstep    ;
   int      nfront   ;
   int      welim    ;
   int      nfind    ;
   int      nzf      ;
   double   ops      ;
   int      nexact2  ;
   int      nexact3  ;
   int      napprox  ;
   int      ncheck   ;
   int      nindst   ;
   int      noutmtch ;
   double   cpu      ;
} ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   this object holds information about a vertex

   id     -- id of the vertex
   mark   -- mark flag, 'O' or 'X'
   status -- status flag
      'L' -- eliminated leaf node
      'E' -- eliminated interior node
      'O' -- outmatched node
      'D' -- vertex is on degree heap
      'R' -- vertex is on reach list
      'I' -- interior segment
      'B' -- boundary segment
   stage -- stage of node, stage 0 nodes are eliminated first,
            stage 1 nodes next, etc
   wght  -- weight of the node
   nadj  -- size of adjacency
   adj   -- pointer to adjacency vector
     when the vertex has not yet been eliminated, adj[nadj]
        contains the uncovered edges
     when the vertex has been eliminated, adj[nadj] contains the
        list of boundary segments
   bndwght -- weight of boundary, used for a root segment
   par     -- pointer to a parent segment
     for a root segment, 
        pointer to the root segment of the parent front
     for a interior segment, 
        pointer to the parent segment in the same front
   subtrees -- pointer to head of IP list with adjacent subtrees
   --------------------------------------------------------------------
*/
struct _MSMDvtx {
   int       id        ;
   char      mark      ;
   char      status    ;
   int       stage     ;
   int       wght      ;
   int       nadj      ;
   int       *adj      ;
   int       bndwght   ;
   MSMDvtx   *par      ;
   IP        *subtrees ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in MSMDinfo.c  -------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   constructor

   created -- 96feb25, cca
   -----------------------
*/
MSMDinfo *
MSMDinfo_new ( 
   void 
) ;
/*
   ---------------------------
   set the default data fields
   
   created -- 96feb25, cca
   ---------------------------
*/
void
MSMDinfo_setDefaultFields(
   MSMDinfo   *info
) ;
/*
   -----------------------
   clear the data fields

   created -- 96feb25, cca
   -----------------------
*/
void
MSMDinfo_clearData ( 
   MSMDinfo   *info
) ;
/*
   -----------------------
   destructor

   created -- 96feb25, cca
   -----------------------
*/
void
MSMDinfo_free (
   MSMDinfo   *info
) ;
/*
   ------------------------------------
   purpose -- print the MSMDinfo object

   created -- 96feb25, cca
   ------------------------------------
*/
void
MSMDinfo_print ( 
   MSMDinfo    *info,
   FILE        *fp 
) ;
/*
   -----------------------------------------
   determine if the MSMDinfo object is valid
 
   created -- 96feb25, cca
   -----------------------------------------
*/
int
MSMDinfo_isValid (
   MSMDinfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c  ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   constructor

   created -- 96feb25, cca
   -----------------------
*/
MSMD *
MSMD_new ( 
   void 
) ;
/*
   ---------------------------
   set the default data fields
   
   created -- 96feb25, cca
   ---------------------------
*/
void
MSMD_setDefaultFields(
   MSMD   *msmd
) ;
/*
   -----------------------
   clear the data fields

   created -- 96feb25, cca
   -----------------------
*/
void
MSMD_clearData ( 
   MSMD   *msmd
) ;
/*
   -----------------------
   destructor

   created -- 96feb25, cca
   -----------------------
*/
void
MSMD_free (
   MSMD   *msmd
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in MSMDvtx.c  --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   print a vertex  object
 
   created -- 96feb25, cca
   -----------------------
*/
void
MSMDvtx_print (
   MSMDvtx    *v,
   FILE   *fp
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in order.c  ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------------
   purpose -- to order the graph using multi-stage minimum degree

   g      -- Graph object
   stages -- stage vector for vertices, 
      if NULL then
         all vertices on stage zero.
      otherwise 
         vertices with stage istage are eliminated 
         before any vertices with stage > istage

   working storage is free'd,
   statistics can be accessed through their variables or printed
   via the void MSMD_printStats(MSMD*,FILE*) method.

   created -- 96feb25, cca
   ---------------------------------------------------------------------
*/
void
MSMD_order ( 
   MSMD       *msmd,
   Graph      *g, 
   int        stages[],
   MSMDinfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in fillPerms.c  ------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------
   fill the two permutation vectors
 
   created -- 96feb24, cca
   --------------------------------
*/
void
MSMD_fillPerms (
   MSMD   *msmd,
   IV     *newToOldIV,
   IV     *oldToNewIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in frontETree.c  -----------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   create and return an ETree object that holds the front tree.
 
   created  -- 96jun23, cca
   ------------------------------------------------------------
*/
ETree *
MSMD_frontETree (
   MSMD   *msmd
) ;
/*--------------------------------------------------------------------*/
/*
========================================================================
===== start of private methods =========================================
========================================================================
*/
/*
------------------------------------------------------------------------
----- methods found in init.c  -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------
   initialization procedure


   created -- 96feb25, cca
   ---------------------------------------------------
*/
void
MSMD_init ( 
   MSMD       *msmd,
   Graph      *g, 
   int        stages[],
   MSMDinfo   *info 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in clearReachSet.c  --------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   clean the vertices in the reach set

   for each v in reach set
      clean subtree list
      clean edge list
   end for

   created -- 95nov08, cca
   -------------------------------------------------
*/
void
MSMD_cleanReachSet ( 
   MSMD       *msmd,
   MSMDinfo   *info
) ;
/*
   ----------------------------------
   clean v's subtree list of children

   created -- 95nov08, cca
   ----------------------------------
*/
void
MSMD_cleanSubtreeList ( 
   MSMD      *msmd,
   MSMDvtx   *v,
   MSMDinfo   *info
) ;
/*
   ----------------------------------------------
   for each uncovered (v,w)
      if v->subtrees \cap w->subtrees != emptyset
         remove (v,w) from uncovered edges
      end if
   end for

   created -- 95nov08, cca
   ----------------------------------------------
*/
void
MSMD_cleanEdgeList ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in eliminate.c  ------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------------
   eliminate all nodes in a stage

   created -- 96feb25, cca
   ---------------------------------------------------------------------
*/
void
MSMD_eliminateStage ( 
   MSMD       *msmd,
   MSMDinfo   *info
) ;
/*
   ------------------------------------------------------
   purpose -- to eliminate an independent set of vertices
 
   created -- 95feb25, cca
   ------------------------------------------------------
*/
int
MSMD_eliminateStep ( 
   MSMD       *msmd,
   MSMDinfo   *info
) ;
/*
   -----------------------------------------
   purpose -- eliminate vertex v
      1) create v's boundary list
      2) merge boundary list onto reach list
      3) for each vertex in the boundary
         3.1) add v to the subtree list

   created -- 96feb25, cca
   -----------------------------------------
*/
void 
MSMD_eliminateVtx ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in findInodes.c  -----------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- to find indistinguishable nodes in the reach set

   flag = 0 --> return
   flag = 1 --> check out nodes that are 2-adj
   flag = 2 --> check out nodes that are both 2-adj and not

   created -- 96feb15, cca
   -----------------------------------------------------------
*/
void
MSMD_findInodes ( 
   MSMD       *msmd,
   MSMDinfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in update.c  ---------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   purpose -- to update vertices in the reach set

   created -- 96feb25, cca
   ----------------------------------------------
*/
void
MSMD_update ( 
   MSMD       *msmd,
   MSMDinfo   *info
) ;
/*
   -------------------------------------------------------
   purpose -- to compute the exact boundary size of a node 
              adjacent to only two eliminated vertices

   created -- 96feb25, cca
   -------------------------------------------------------
*/
int
MSMD_exactDegree2 ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) ;
/*
   ------------------------------------------------------------
   purpose -- to compute the exact boundary size of a node that
              is not adjacent to only two eliminated vertices

   created -- 96feb25, cca
   ------------------------------------------------------------
*/
int
MSMD_exactDegree3 ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) ;
/*
   --------------------------------------------------------
   purpose -- to compute the approximate degree of a vertex

   created -- 96feb25, cca
   --------------------------------------------------------
*/
int
MSMD_approxDegree ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in makeSchurComplement.c  --------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------------
   purpose -- 

   if the elimination has halted before all the stages have been 
   eliminated, then create the schur complement graph and the map 
   from the original vertices those in the schur complement graph.

   schurGraph -- Graph object to contain the schur complement graph
   VtoPhi     -- IV object to contain the map from vertices in V
                 to schur complement vertices in Phi

   created -- 97feb01, cca
   ----------------------------------------------------------------
*/
void
MSMD_makeSchurComplement (
   MSMD    *msmd,
   Graph   *schurGraph,
   IV      *VtoPhiIV
) ;
/*--------------------------------------------------------------------*/
