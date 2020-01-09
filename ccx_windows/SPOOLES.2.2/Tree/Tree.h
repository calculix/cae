/*  Tree.h  */

#include "../cfiles.h"
#include "../IV.h"
#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------
   simple tree object

   created -- 95nov15, cca
   -----------------------
*/
typedef struct _Tree   Tree ;
struct _Tree {
   int   n    ;
   int   root ;
   int   *par ;
   int   *fch ;
   int   *sib ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- create and return a new Tree object

   created -- 95nov15, cca
   -----------------------------------------------
*/
Tree *
Tree_new ( 
   void
) ;
/*
   ------------------------------------------------------
   purpose -- set the default fields for the Tree object

   created -- 95nov15, cca
   ------------------------------------------------------
*/
void
Tree_setDefaultFields (
   Tree   *tree
) ;
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 95nov15, cca
   --------------------------------
*/
void
Tree_clearData (
   Tree   *tree
) ;
/*
   --------------------------------
   purpose -- free the Tree object

   created -- 95nov15, cca
   --------------------------------
*/
void
Tree_free (
   Tree   *tree
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   purpose -- return the number of nodes in the tree
 
   created -- 98jun12, cca
   -------------------------------------------------
*/
int
Tree_nnodes (
   Tree   *tree
) ;
/*
   --------------------------------------
   purpose -- return the root of the tree
 
   created -- 98jun12, cca
   --------------------------------------
*/
int
Tree_root (
   Tree   *tree
) ;
/*
   ------------------------------------------------
   purpose -- return a pointer to the parent vector
 
   created -- 98jun12, cca
   ------------------------------------------------
*/
int *
Tree_par (
   Tree   *tree
) ;
/*
   -----------------------------------------------------
   purpose -- return a pointer to the first child vector
 
   created -- 98jun12, cca
   -----------------------------------------------------
*/
int *
Tree_fch (
   Tree   *tree
) ;
/*
   -------------------------------------------------
   purpose -- return a pointer to the sibling vector
 
   created -- 98jun12, cca
   -------------------------------------------------
*/
int *
Tree_sib (
   Tree   *tree
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor

   created -- 95nov15, cca
   -----------------------
*/
void
Tree_init1 (
   Tree   *tree,
   int    size
) ;
/*
   --------------------------------
   initialize given a parent vector

   created -- 95nov15, cca
   --------------------------------
*/
void
Tree_init2 (
   Tree   *tree,
   int    size,
   int    par[]
) ;
/*
   ---------------------------------
   initialize given the tree vectors

   created -- 95nov15, cca
   ---------------------------------
*/
void
Tree_init3 (
   Tree   *tree,
   int    size,
   int    par[],
   int    fch[],
   int    sib[]
) ;
/*
   ------------------------------------
   set the fch[], sib[] and root fields

   created -- 95nov15, cca
   ------------------------------------
*/
void
Tree_setFchSibRoot (
   Tree   *tree
) ;
/*
   -----------------------
   set the root field

   created -- 95nov15, cca
   -----------------------
*/
void
Tree_setRoot (
   Tree   *tree
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in maximizeGain.c ----------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------
   purpose -- 
 
   given a gain value assigned to each node,
   find a set of nodes, no two in a child-ancestor
   relationship, that maximizes the total gain.
 
   this problem arises in finding the optimal domain/schur 
   complement partition for a semi-implicit factorization.
 
   created -- 98jun20, cca
   -------------------------------------------------------
*/
IV *
Tree_maximizeGainIV (
   Tree   *tree,
   IV     *gainIV,
   int    *ptotalgain,
   int    msglvl,
   FILE   *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in subtree.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   purpose -- to initialize subtree with the subtree 
              of tree using nodes in nodeidsIV
 
   return values ---
      1 -- normal return
     -1 -- subtree is NULL
     -2 -- nodeidsIV is NULL
     -3 -- tree is NULL
     -4 -- nodeidsIV is invalid
 
   created -- 98oct15, cca
   -------------------------------------------------
*/
int
Tree_initFromSubtree (
   Tree   *subtree,
   IV     *nodeidsIV,
   Tree   *tree
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   return the first vertex in a post-order traversal
   
   created -- 95nov15, cca
   -------------------------------------------------
*/
int
Tree_postOTfirst (
   Tree   *tree
) ;
/*
   ----------------------------------------------------------
   return the vertex that follows v in a post-order traversal
   ----------------------------------------------------------
*/
int
Tree_postOTnext (
   Tree   *tree,
   int    v
) ;
/*
   ------------------------------------------------
   return the first vertex in a pre-order traversal
   
   created -- 95nov15, cca
   ------------------------------------------------
*/
int
Tree_preOTfirst (
   Tree   *tree
) ;
/*
   ---------------------------------------------------------
   return the vertex that follows v in a pre-order traversal
   
   created -- 95nov15, cca
   ---------------------------------------------------------
*/
int
Tree_preOTnext (
   Tree   *tree,
   int    v
) ;
/*
   ---------------------------------------
   return the number of leaves in the tree

   created -- 95nov15, cca
   ---------------------------------------
*/
int
Tree_nleaves (
   Tree   *tree
) ;
/*
   -----------------------------------------------
   return the number of roots of the tree (forest)

   created -- 95nov15, cca
   -----------------------------------------------
*/
int
Tree_nroots (
   Tree   *tree
) ;
/*
   -----------------------------------------
   return the number of children of a vertex

   created -- 95nov15, cca
   -----------------------------------------
*/
int
Tree_nchild (
   Tree   *tree,
   int    v
) ;
/*
   -------------------------------------------
   this method returns an IV object that holds 
   the number of children for the tree nodes.
 
   created -- 96dec18, cca
   -------------------------------------------
*/
IV *
Tree_nchildIV (
   Tree   *tree
) ;
/*
   -----------------------------
   return the height of the tree
 
   created -- 96aug23, cca
   -----------------------------
*/
int
Tree_height (
   Tree   *tree
) ;
/*
   -------------------------------------------------------------
   return the maximum number of children of any node in the tree
 
   created -- 96sep05, cca
   -------------------------------------------------------------
*/
int
Tree_maxNchild (
   Tree   *tree
) ;
/*
   ---------------------------------------------
   return the number of bytes used by the object
   ---------------------------------------------
*/
int
Tree_sizeOf (
   Tree   *tree
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in metrics.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------
   create and return a subtree metric IV object
   input  : vmetricIV -- a metric defined on the vertices
   return : tmetricIV -- a metric defined on the subtrees
  
   created -- 96jun23, cca
   ------------------------------------------------------
*/
IV *
Tree_setSubtreeImetric (
   Tree   *tree,
   IV     *vmetricIV
) ;
/*
   ------------------------------------------------------
   create and return a subtree metric DV object
   input  : vmetricDV -- a metric defined on the vertices
   return : tmetricDV -- a metric defined on the subtrees
 
   created -- 96jun23, cca
   ------------------------------------------------------
*/
DV *
Tree_setSubtreeDmetric (
   Tree   *tree,
   DV     *vmetricDV
) ;
/*
   ------------------------------------------------------------
   create and return a depth metric IV object
   input  : vmetricIV -- a metric defined on the vertices
   output : dmetric[] -- a depth metric defined on the vertices
 
   dmetric[u] = vmetric[u] + dmetric[par[u]] if par[u] != -1
              = vmetric[u]                   if par[u] == -1
 
   created -- 96jun23, cca
   ------------------------------------------------------------
*/
IV *
Tree_setDepthImetric (
   Tree   *tree,
   IV     *vmetricIV
) ;
/*
   ------------------------------------------------------------
   create and return a depth metric DV object
   input  : vmetricDV -- a metric defined on the vertices
   output : dmetric[] -- a depth metric defined on the vertices
 
   dmetric[u] = vmetric[u] + dmetric[par[u]] if par[u] != -1
              = vmetric[u]                   if par[u] == -1
 
   created -- 96jun23, cca
   ------------------------------------------------------------
*/
DV *
Tree_setDepthDmetric (
   Tree   *tree,
   DV     *vmetricDV
) ;
/*
   ------------------------------------------------------------------
   create and return a height metric IV object
   input  : vmetricIV -- a metric defined on the vertices
   output : dmetricIV -- a depth metric defined on the vertices
 
   hmetric[v] = vmetric[v] + max{p(u) = v} hmetric[u] if fch[v] != -1
              = vmetric[v]                            if fch[v] == -1
 
   created -- 96jun23, cca
   ------------------------------------------------------------------
*/
IV *
Tree_setHeightImetric (
   Tree   *tree,
   IV     *vmetricIV
) ;
/*
   ------------------------------------------------------------------
   create and return a height metric DV object
   input  : vmetricDV -- a metric defined on the vertices
   output : dmetricDV -- a depth metric defined on the vertices
 
   hmetric[v] = vmetric[v] + max{p(u) = v} hmetric[u] if fch[v] != -1
              = vmetric[v]                            if fch[v] == -1
 
   created -- 96jun23, cca
   ------------------------------------------------------------------
*/
DV *
Tree_setHeightDmetric (
   Tree   *tree,
   DV     *vmetricDV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in justify.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   left-justify a tree by subtree size
   children are linked in ascending order of their subtree size
 
   created -- 95nov15, cca
   ------------------------------------------------------------
*/
void
Tree_leftJustify (
   Tree   *tree
) ;
/*
   ------------------------------------------------------
   left-justify a tree by a metric
   children are linked in ascending order of their metric
 
   created -- 96jun23, cca
   ------------------------------------------------------
*/
void
Tree_leftJustifyI (
   Tree   *tree,
   IV     *metricIV
) ;
/*
   ------------------------------------------------------
   left-justify a tree by a metric
   children are linked in ascending order of their metric
 
   created -- 96jun23, cca
   ------------------------------------------------------
*/
void
Tree_leftJustifyD (
   Tree   *tree,
   DV     *metricDV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in perms.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------
   fill the new-to-old permutation vector
 
   created -- 95nov15, cca
   --------------------------------------
*/
void
Tree_fillNewToOldPerm (
   Tree   *tree,
   int    newToOld[]
) ;
/*
   --------------------------------------
   fill the old-to-new permutation vector
 
   created -- 95nov15, cca
   --------------------------------------
*/
void
Tree_fillOldToNewPerm (
   Tree   *tree,
   int    oldToNew[]
) ;
/*
   ------------------------------------------------------
   fill the new-to-old and old-to-new permutation vectors
 
   created -- 95nov15, cca
   ------------------------------------------------------
*/
void
Tree_fillBothPerms (
   Tree   *tree,
   int    newToOld[],
   int    oldToNew[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   purpose -- to read in an Tree object from a file
 
   input --
 
      fn -- filename, must be *.treeb or *.treef
 
   return value -- 1 if success, 0 if failure
 
   created -- 95nov15, cca
   ------------------------------------------------
*/
int
Tree_readFromFile ( 
   Tree   *tree, 
   char   *fn 
) ;
/*
   -------------------------------------------------------
   purpose -- to read an Tree object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 95nov15, cca
   -------------------------------------------------------
*/
int
Tree_readFromFormattedFile (
   Tree   *tree,
   FILE   *fp
) ;
/*
   ---------------------------------------------------
   purpose -- to read an Tree object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 95nov15, cca
   ---------------------------------------------------
*/
int
Tree_readFromBinaryFile (
   Tree    *tree,
   FILE   *fp
) ;
/*
   --------------------------------------------
   purpose -- to write an Tree object to a file
 
   input --
 
      fn -- filename
        *.treeb -- binary
        *.treef -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov15, cca
   --------------------------------------------
*/
int
Tree_writeToFile (
   Tree   *tree,
   char   *fn
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write an Tree object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov15, cca
   ------------------------------------------------------
*/
int
Tree_writeToFormattedFile (
   Tree   *tree,
   FILE   *fp
) ;
/*
   ---------------------------------------------------
   purpose -- to write an Tree object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov15, cca
   ---------------------------------------------------
*/
int
Tree_writeToBinaryFile (
   Tree    *tree,
   FILE   *fp
) ;
/*
   --------------------------------------------------
   purpose -- to write an Tree object for a human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov15, cca
   --------------------------------------------------
*/
int
Tree_writeForHumanEye (
   Tree    *tree,
   FILE   *fp
) ;
/*
   ----------------------------------------------------------
   purpose -- to write out the statistics for the Tree object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov15, cca
   ----------------------------------------------------------
*/
int
Tree_writeStats (
   Tree    *tree,
   FILE   *fp
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in compress.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------
   create and return an IV object that contains
   the map from vertices to fundamental chains.
 
   return value -- # of fundamental chains
 
   created -- 96jun23, cca
   -------------------------------------------
*/
IV *
Tree_fundChainMap (
   Tree   *tree
) ;
/*
   -----------------------------------------------------------------
   compress a tree based on a map from old vertices to new vertices.
   the restriction on the map is that the set {u | map[u] = U} must
   be connected for all U.
 
   created  -- 95nov15, cca
   modified -- 96jan04, cca
      bug fixed, N computed incorrectly
   modified -- 96jun23, cca
      in calling sequence, int map[] converted to IV *mapIV 
   -----------------------------------------------------------------
*/
Tree *
Tree_compress (
   Tree   *tree,
   IV     *mapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in permute.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   return a permuted tree
 
   created -- 96jan04, cca
   -----------------------
*/
Tree *
Tree_permute (
   Tree   *tree,
   int    newToOld[],
   int    oldToNew[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in setBoxes.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   purpose -- fill boxes arrays for display of a tree
 
   vmetric[] -- vector of metric on the nodes
   tmetric[] -- vector of metric on the subtrees
   xmin, xmax, ymin, ymax -- bounds on box for root
   west[], east[], south[], north[] -- vector to hold bounds for
                                       the nodes in the tree
 
   return value --
     1 --> success
     2 --> no success, maxnchild > 3
 
   created -- 96dec20, cca
   -------------------------------------------------------------
*/
int
Tree_setBoxesII (
   Tree     *tree,
   int      vmetric[],
   int      tmetric[],
   double   xmin,
   double   xmax,
   double   ymin,
   double   ymax,
   double   west[],
   double   east[],
   double   south[],
   double   north[]
) ;
/*
   -------------------------------------------------------------
   purpose -- fill boxes arrays for display of a tree
 
   vmetric[] -- vector of metric on the nodes
   tmetric[] -- vector of metric on the subtrees
   xmin, xmax, ymin, ymax -- bounds on box for root
   west[], east[], south[], north[] -- vector to hold bounds for
                                       the nodes in the tree
 
   return value --
     1 --> success
     2 --> no success, maxnchild > 3
 
   created -- 96dec20, cca
   -------------------------------------------------------------
*/
int
Tree_setBoxesDD (
   Tree     *tree,
   double   vmetric[],
   double   tmetric[],
   double   xmin,
   double   xmax,
   double   ymin,
   double   ymax,
   double   west[],
   double   east[],
   double   south[],
   double   north[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in getCoords.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   purpose -- to get simple x[] and y[] coordinates
              for the tree vertices

   return values --
      1 -- normal return
     -1 -- tree is NULL
     -2 -- heightflag is invalid
     -3 -- coordflag is invalid
     -4 -- xDV is NULL
     -5 -- yDV is NULL

   created -- 99jan07, cca
   ------------------------------------------------
*/
int
Tree_getSimpleCoords (
   Tree   *tree,
   char   heightflag,
   char   coordflag,
   DV     *xDV,
   DV     *yDV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in draw.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   purpose -- to write an EPS file with a picture of a tree.
              each node can have its own radius and label

   filename  -- name of the file to be written
   xDV       -- x coordinates
   yDV       -- y coordinates
   rscale    -- scaling factor for radius of nodes
   radiusDV  -- radius of nodes, if NULL then radius = 1
   labelflag -- flag to specify whether labels are to be drawn
           1     -- draw labels
       otherwise -- do not draw labels
   fontscale -- scaling factor for font
   labelsIV  -- IV object that contains the labels of the nodes.
       if NULL then the node ids are used
   bbox[] -- bounding box for figure
      bbox[0] -- x_min
      bbox[1] -- y_min
      bbox[2] -- x_max
      bbox[3] -- y_max
   frame[] -- frame to hold tree
      frame[0] -- x_min
      frame[1] -- y_min
      frame[2] -- x_max
      frame[3] -- y_max
   bounds[] -- bounds for local coordinates
      if bounds is NULL then
         the tree fills the frame. note, this is a nonlinear process
         when the nodes have non-constant radii, and may not converge
         when the maximum radius is large when compared to the frame.
         if the process does not converge, a message is printed and
        the program exits.
      else
         bounds[0] -- xi_min
         bounds[1] -- eta_min
         bounds[2] -- xi_max
         bounds[3] -- eta_max
      endif

   recommendations,
      bbox[] = { 0, 0, 500, 200 } for tall skinny trees
               { 0, 0, 500, 500 } for wide trees
      frame[0] = bbox[0] + 10
      frame[1] = bbox[1] + 10
      frame[2] = bbox[2] - 10
      frame[3] = bbox[3] - 10

   return value
      1 -- normal return
     -1 -- tree is NULL
     -2 -- filename is NULL
     -3 -- xDV is NULL
     -4 -- yDV is NULL
     -5 -- rscale is negative
     -6 -- fontscale is negative
     -7 -- bbox is NULL
     -8 -- frame is NULL

   created -- 99jan07, cca
   ------------------------------------------------------------
*/
int
Tree_drawToEPS (
   Tree     *tree,
   char     *filename,
   DV       *xDV,
   DV       *yDV,
   double   rscale,
   DV       *radiusDV,
   int      labelflag,
   double   fontscale,
   IV       *labelsIV,
   double   bbox[],
   double   frame[],
   double   bounds[]
) ;
/*--------------------------------------------------------------------*/
