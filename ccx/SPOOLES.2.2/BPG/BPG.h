/*  BPG.h  */

#include "../Graph.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   the BPG object holds a bipartite graph.
   the edge set E \subseteq X \times Y,
   where X and Y are two sets of vertices,
      X = { 0, 1, ..., nX - 1}
      Y = { nX, nX + 1, ..., nX + nY - 1 }

   the BPG object contains a Graph object graph,
   a hack because C does not possess inheritance.

   created -- 95oct06, cca
   ----------------------------------------------
*/
typedef struct _BPG   BPG ;
struct _BPG {
   int     nX       ;
   int     nY       ;
   Graph   *graph   ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- create and return a new BPG object

   created -- 95oct06, cca
   -----------------------------------------------
*/
BPG *
BPG_new ( 
   void
) ;
/*
   ------------------------------------------------------
   purpose -- set the default fields for the BPG object

   created -- 95oct06, cca
   ------------------------------------------------------
*/
void
BPG_setDefaultFields (
   BPG   *bpg
) ;
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 95oct06, cca
   --------------------------------
*/
void
BPG_clearData (
   BPG   *bpg
) ;
/*
   --------------------------------
   purpose -- free the BPG object

   created -- 95oct06, cca
   --------------------------------
*/
void
BPG_free (
   BPG   *bpg
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in DM.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------
   compute the generalized Dulmadge-Mendolsohn decomposition

   bpg     -- BPG bipartite graph object
   dmflags -- flags vector
                   / 0 if x in X_R
      dmflags[x] = | 1 if x in X_I
                   \ 2 if x in X_E
                   / 0 if y in Y_R
      dmflags[y] = | 1 if y in Y_I
                   \ 2 if y in Y_E
   stats -- statistics vector
      stats[0] -- weight of X_I
      stats[1] -- weight of X_E
      stats[2] -- weight of X_R
      stats[3] -- weight of Y_I
      stats[4] -- weight of Y_E
      stats[5] -- weight of Y_R

   created -- 95oct12, cca
   ---------------------------------------------------------
*/
void
BPG_DMdecomposition (
   BPG    *bpg,
   int    dmflags[],
   int    stats[],
   int    msglvl,
   FILE   *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in maxFlow.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------
   find the DM decomposition of a weighted bipartite graph
   using a simple max flow algorithm
 
   bpg     -- BPG bipartite graph object
   dmflags -- flags vector
                   / 0 if x in X_R
      dmflags[x] = | 1 if x in X_I
                   \ 2 if x in X_E
                   / 0 if y in Y_R
      dmflags[y] = | 1 if y in Y_I
                   \ 2 if y in Y_E
   stats -- statistics vector
      stats[0] -- weight of X_I
      stats[1] -- weight of X_E
      stats[2] -- weight of X_R
      stats[3] -- weight of Y_I
      stats[4] -- weight of Y_E
      stats[5] -- weight of Y_R
 
   created -- 96mar08, cca
   -------------------------------------------------------
*/
void
BPG_DMviaMaxFlow (
   BPG    *bpg,
   int    dmflags[],
   int    stats[],
   int    msglvl,
   FILE   *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   initialize a BPG object, set bpg->nX = nX, bpg->nY = nY,
   and bpg->graph = graph. convert graph into bipartite form, i.e.,
      adj(v) := adj(v) \cap {nX, ..., nX + nY - 1} for v in [0, nX)
      adj(v) := adj(v) \cap {0, ..., nX - 1}       for v in [nX, nX+nY)

   created -- 95oct06, cca
   --------------------------------------------------------------------
*/
void
BPG_init (
   BPG     *bpg,
   int     nX,
   int     nY,
   Graph   *graph
) ;
/*
   -----------------------------------------------------------
   initialize from a graph given a color vector and two target
   vectors. this method can be used to get the bipartite graph
   of the boundary of two components.

   graph  -- graph from which to extract the bipartite graph
   colors -- vector of colors for the vertices
   cX     -- color for the X vertices
   cY     -- color for the Y vertices
   cmap   -- map from vertices in g to vertices in bpg
   indX   -- vector to hold the global indices in X
   indY   -- vector to hold the global indices in Y

   created -- 95oct06, cca
   -----------------------------------------------------------
*/
void
BPG_initFromColoring (
   BPG     *bpg,
   Graph   *graph,
   int     colors[],
   int     cX,
   int     cY,
   int     cmap[],
   int     indX[],
   int     indY[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   purpose -- to read in a BPG object from a file

   input --

      fn -- filename, must be *.bpgb or *.bpgf

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
BPG_readFromFile ( 
   BPG    *bpg, 
   char   *fn 
) ;
/*
   --------------------------------------------------------
   purpose -- to read a BPG object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   --------------------------------------------------------
*/
int
BPG_readFromFormattedFile ( 
   BPG   *bpg, 
   FILE    *fp 
) ;
/*
   ----------------------------------------------------
   purpose -- to read a BPG object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95oct06, cca
   ----------------------------------------------------
*/
int
BPG_readFromBinaryFile ( 
   BPG   *bpg, 
   FILE    *fp 
) ;
/*
   --------------------------------------------
   purpose -- to write a BPG object to a file

   input --

      fn -- filename
        *.bpgb -- binary
        *.bpgf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   --------------------------------------------
*/
int
BPG_writeToFile ( 
   BPG   *bpg, 
   char    *fn 
) ;
/*
   ------------------------------------------------------
   purpose -- to write a BPG object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   ------------------------------------------------------
*/
int
BPG_writeToFormattedFile ( 
   BPG   *bpg, 
   FILE  *fp 
) ;
/*
   ---------------------------------------------------
   purpose -- to write a BPG object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   ---------------------------------------------------
*/
int
BPG_writeToBinaryFile ( 
   BPG    *bpg, 
   FILE   *fp 
) ;
/*
   -------------------------------------------------
   purpose -- to write a BPG object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -------------------------------------------------
*/
int
BPG_writeForHumanEye ( 
   BPG    *bpg, 
   FILE   *fp 
) ;
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the BPG object

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -----------------------------------------------------------
*/
int
BPG_writeStats ( 
   BPG    *bpg, 
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in makeGraphs.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- return the X by X graph where (x1,x2) is an edge
              if there exists a y in Y such that (x1,y) and
              (x2,y) are edges in the bipartite graph.

   created -- 95dec06, cca
              to be used in the BKL object.
   -----------------------------------------------------------
*/
Graph *
BPG_makeGraphXbyX (
   BPG   *bpg
) ;
/*
   -----------------------------------------------------------
   purpose -- return the Y by Y graph where (y1,y2) is an edge
              if there exists a x in X such that (x,y1) and
              (x,y2) are edges in the bipartite graph.

   created -- 95dec07, cca
   -----------------------------------------------------------
*/
Graph *
BPG_makeGraphYbyY (
   BPG   *bpg
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in pseudo.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------
   return a pseudoperipheral node

   created -- 95oct07, cca
   ------------------------------
*/
int
BPG_pseudoperipheralnode (
   BPG   *bpg,
   int   seed
) ;
/*
   ----------------------------------------------------
   return value -- # of vertices in the level structure

   created -- 95oct07, cca
   ----------------------------------------------------
*/
int
BPG_levelStructure (
   BPG   *bpg,
   int   root,
   int   list[],
   int   dist[],
   int   mark[],
   int   tag
) ;
/*--------------------------------------------------------------------*/
