/*  ETree.h  */

#include "../SPOOLES.h"
#include "../cfiles.h"
#include "../Graph.h"
#include "../Tree.h"
#include "../IV.h"
#include "../DV.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   the ETree object is a tree that has a weight associated with
   each node and a weight associated with each node's boundary.
   it is useful to model:
   (1) a vertex elimination tree (for a unit weight graph), 
   (2) a compressed vertex elimination tree (for a compressed graph),
   (3) a front tree (for a factor graph)

   nfront       -- # of fronts
   nvtx         -- # of vertices
   tree         -- pointer to a Tree object, size nfront 
   nodwghtsIV   -- IV object of node weights, size nfront
   bnwwghtsIV   -- IV object of node boundary weights, size nfront
   vtxToFrontIV -- IV object that holds the map from vertices to fronts,
                 size nvtx

   created -- 96jun23, cca
   ---------------------------------------------------------------------
*/
typedef struct _ETree   ETree ;
struct _ETree {
   int    nfront        ;
   int    nvtx          ;
   Tree   *tree         ;
   IV     *nodwghtsIV   ;
   IV     *bndwghtsIV   ;
   IV     *vtxToFrontIV ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in basics.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- create and return a new ETree object

   created -- 95nov15, cca
   -----------------------------------------------
*/
ETree *
ETree_new ( 
   void
) ;
/*
   ------------------------------------------------------
   purpose -- set the default fields for the ETree object

   created -- 95nov15, cca
   ------------------------------------------------------
*/
void
ETree_setDefaultFields (
   ETree   *etree
) ;
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 95nov15, cca
   --------------------------------
*/
void
ETree_clearData (
   ETree   *etree
) ;
/*
   --------------------------------
   purpose -- free the ETree object

   created -- 95nov15, cca
   --------------------------------
*/
void
ETree_free (
   ETree   *etree
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in instance.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------
   return the number of fronts
 
   created -- 97feb28, cca
   ---------------------------
*/
int
ETree_nfront (
   ETree   *etree
) ;
/*
   -----------------------------
   return the number of vertices
 
   created -- 97feb28, cca
   -----------------------------
*/
int
ETree_nvtx (
   ETree   *etree
) ;
/*
   -----------------------------------
   return a pointer to the Tree object
 
   created -- 97feb28, cca
   -----------------------------------
*/
Tree *
ETree_tree (
   ETree   *etree
) ;
/*
   ---------------------------
   return the root of the tree
 
   created -- 97feb28, cca
   ---------------------------
*/
int
ETree_root (
   ETree   *etree
) ;
/*
   -------------------------------------
   return a pointer to the parent vector
 
   created -- 97feb28, cca
   -------------------------------------
*/
int *
ETree_par (
   ETree   *etree
) ;
/*
   ------------------------------------------
   return a pointer to the first child vector
 
   created -- 97feb28, cca
   ------------------------------------------
*/
int *
ETree_fch (
   ETree   *etree
) ;
/*
   --------------------------------------
   return a pointer to the sibling vector
 
   created -- 97feb28, cca
   --------------------------------------
*/
int *
ETree_sib (
   ETree   *etree
) ;
/*
   ------------------------------------------
   return a pointer to the nodwghts IV object
 
   created -- 97feb28, cca
   ------------------------------------------
*/
IV *
ETree_nodwghtsIV (
   ETree   *etree
) ;
/*
   -------------------------------------------
   return a pointer to the nodwghts int vector
 
   created -- 97feb28, cca
   -------------------------------------------
*/
int *
ETree_nodwghts (
   ETree   *etree
) ;
/*
   ------------------------------------------
   return a pointer to the bndwghts IV object
 
   created -- 97feb28, cca
   ------------------------------------------
*/
IV *
ETree_bndwghtsIV (
   ETree   *etree
) ;
/*
   -------------------------------------------
   return a pointer to the bndwghts int vector
 
   created -- 97feb28, cca
   -------------------------------------------
*/
int *
ETree_bndwghts (
   ETree   *etree
) ;
/*
   --------------------------------------------
   return a pointer to the vtxToFront IV object
 
   created -- 97feb28, cca
   --------------------------------------------
*/
IV *
ETree_vtxToFrontIV (
   ETree   *etree
) ;
/*
   ---------------------------------------------
   return a pointer to the vtxToFront int vector
 
   created -- 97feb28, cca
   ---------------------------------------------
*/
int *
ETree_vtxToFront (
   ETree   *etree
) ;
/*
   ------------------------------------------------
   purpose -- return the number of internal degrees
              of freedom in front J
 
   created -- 97may23, cca
   ------------------------------------------------
*/
int
ETree_frontSize (
   ETree   *etree,
   int     J
) ;
/*
   ------------------------------------------------
   purpose -- return the number of external degrees 
              of freedom in front J
 
   created -- 97may23, cca
   ------------------------------------------------
*/
int
ETree_frontBoundarySize (
   ETree   *etree,
   int     J
) ;
/*
   ------------------------------------------------------------
   purpose -- compute the maximum number of indices and entries
              in a front
 
   symflag = 1 -->
      count only column indices
      count upper entries in (1,1) block and (1,2) block
   symflag = 2 --> 
      count row and column indices
      count entries in (1,1), (1,2) and (2,1) blocks
 
   created -- 97may23, cca
   ------------------------------------------------------------
*/
void
ETree_maxNindAndNent (
   ETree   *etree,
   int     symflag,
   int     *pmaxnind,
   int     *pmaxnent
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in util.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object
 
   created -- 95nov15, cca
   ----------------------------------------------
*/
int
ETree_sizeOf (
   ETree   *etree
) ;
/*
   ----------------------------------------
   return the number of factor indices
 
   created  -- 95nov15, cca
   modified -- 96jan11, cca
   ----------------------------------------
*/
int
ETree_nFactorIndices (
   ETree   *etree
) ;
/*
   ------------------------------------------
   return the number of factor entries
 
   symflag -- symmetry flag
     0 (SPOOLES_SYMMETRIC)    -- symmetric
     1 (SPOOLES_HERMITIAN)    -- hermitian
     2 (SPOOLES_NONSYMMETRIC) -- nonsymmetric
 
   created -- 98jun05, cca
   ------------------------------------------
*/
int
ETree_nFactorEntries (
   ETree   *etree,
   int     symflag
) ;
/*
   ------------------------------------------
   return the number of factor operations
 
   type -- type of matrix entries
     1 (SPOOLES_REAL)    -- real entries
     2 (SPOOLES_COMPLEX) -- complex entries
   symflag -- symmetry flag
     0 (SPOOLES_SYMMETRIC)    -- symmetric
     1 (SPOOLES_HERMITIAN)    -- hermitian
     2 (SPOOLES_NONSYMMETRIC) -- nonsymmetric
 
   created -- 98jun05, cca
   ------------------------------------------
*/
double
ETree_nFactorOps (
   ETree   *etree,
   int     type,
   int     symflag
) ;
/*
   ----------------------------------------
   return the number of entries an LU front
 
   created -- 96dec04, cca
   ----------------------------------------
*/
double
ETree_nFactorEntriesInFront (
   ETree   *etree,
   int     symflag,
   int     J
) ;
/*
   -------------------------------------------------------
   return the number of internal LU operations for a front
 
   created -- 96dec04, cca
   -------------------------------------------------------
*/
double
ETree_nInternalOpsInFront (
   ETree   *etree,
   int     type,
   int     symflag,
   int     J
) ;
/*
   -------------------------------------------------------
   return the number of external LU operations for a front
 
   created -- 96dec04, cca
   -------------------------------------------------------
*/
double
ETree_nExternalOpsInFront (
   ETree   *etree,
   int     type,
   int     symflag,
   int     J
) ;
/*
   ------------------------------------
   return an IV object that contains
   the number of entries for each front
 
   created -- 98jan30, cca
   ------------------------------------
*/
IV *
ETree_factorEntriesIV (
   ETree   *etree,
   int     symflag
) ;
/*
   ---------------------------------------------------------
   return a DV object that contains the number of operations
   for each front using a backward looking algorithm
 
   created -- 96dec04, cca
   ---------------------------------------------------------
*/
DV *
ETree_backwardOps (
   ETree   *etree,
   int     type,
   int     symflag,
   int     *vwghts,
   IVL     *symbfacIVL
) ;
/*
   ---------------------------------------------------------
   return a DV object that contains the number of operations
   for each front using a forward-looking algorithm
 
   created -- 96dec04, cca
   ---------------------------------------------------------
*/
DV *
ETree_forwardOps (
   ETree   *etree,
   int     type,
   int     symflag
) ;
/*
   ---------------------------------------------------------------
   given an IV object that maps uncompressed vertices to vertices,
   create and return an ETree object that is relative to the 
   uncompressed graph.

   created -- 97feb13, cca
   ---------------------------------------------------------------
*/
ETree *
ETree_expand (
   ETree   *etree,
   IV      *eqmapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in init.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   initialize the object given the number of nodes
 
   created -- 95nov15, cca
   -----------------------------------------------
*/
void
ETree_init1 (
   ETree   *etree,
   int     nfront,
   int     nvtx
) ;
/*
   ----------------------------------------
   initialize the ETree object from a graph
 
   created -- 95nov15, cca
   ----------------------------------------
*/
void
ETree_initFromGraph (
   ETree   *etree,
   Graph   *g
) ;
/*
   --------------------------------------------------------
   initialize the ETree object from a graph and permutation
 
   created -- 95nov15, cca
   --------------------------------------------------------
*/
void
ETree_initFromGraphWithPerms (
   ETree   *etree,
   Graph   *g,
   int     newToOld[],
   int     oldToNew[]
) ;
/*
   --------------------------------------------------------------
   purpose -- initialize the front tree for a dense matrix
   
   n -- size of the matrix
   option -- mapping option
      1 --> have all fronts (save the last) contain the same
            number of vertices
      2 --> have all fronts have roughly equal numbers of entries
 
   created -- 96aug19, cca
   --------------------------------------------------------------
*/
void
ETree_initFromDenseMatrix (
   ETree   *etree,
   int     n,
   int     option,
   int     param
) ;
/*
   -------------------------------------
   initialize the ETree object
   (1) read ETree object from file
   (2) get the old-to-new permutation
   (3) permute the ETree
   (4) return the old-to-new permutation
 
   created -- 97jul13, cca
   -------------------------------------
*/
IV *
ETree_initFromFile (
   ETree   *frontETree,
   char    *inETreeFileName,
   int     msglvl,
   FILE    *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in ms.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   returns a compidsIV IV object that maps the
   vertices to a domain (compids[v] > 1)
   or to the multisector (compids[v] = 0).
   the vertices in the multisector is specified
   by their depth of their front in the tree.
 
   created -- 96jan04, cca
   ------------------------------------------------
*/
IV *
ETree_msByDepth (
   ETree   *etree,
   int     depth
) ;
/*
   ----------------------------------------------------------------
  construct a multisector based on vertices found in a subtree.
 
   created -- 96jan04, cca
   ----------------------------------------------------------------
*/
IV *
ETree_msByNvtxCutoff (
   ETree    *etree,
   double   cutoff
) ;
/*
   --------------------------------------------------
   construct a multisector based on the number 
   of factor entries found in a subtree.
 
   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC
     SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   created -- 96jan04, cca
   --------------------------------------------------
*/
IV *
ETree_msByNentCutoff (
   ETree    *etree,
   double   cutoff,
   int      symflag
) ;
/*
   --------------------------------------------------
   construct a multisector based on the number
   of factor operations found in a subtree.
 
   type -- type of entries,
     SPOOLES_REAL or SPOOLES_COMPLEX
 
   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC
     SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   created -- 96jan04, cca
   --------------------------------------------------
*/
IV *
ETree_msByNopsCutoff (
   ETree    *etree,
   double   cutoff,
   int      type,
   int      symflag
) ;
/*
   --------------------------------------------------------------
   purpose -- given a front tree and a multisector map vector,
     fill the map vector with domain ids and the three statistics
     arrays with domain and schur complement statistics.
 
   frontETree -- front tree object, unchanged on output
   msIV -- map from fronts to domains or schur complement
     on input, ms[J] = 0 --> J is in the schur complement
               ms[J] = 1 --> J is not in the schur complement
     on output, ms[J] =  0 --> J is in the schur complement
                ms[J] != 0 --> J is in domain ms[J]
   on output
      nvtxIV -- nvtx[ireg] = # of dof in region ireg
      nzfIV  -- nzf[ireg] = # of factor entries in region ireg
      opsIV  -- ops[ireg] = # of factor ops in region ireg
 
   type -- type of entries, SPOOLES_REAL or SPOOLES_COMPLEX
 
   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC
     SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   created -- 98jan30, cca
   --------------------------------------------------------------
*/
void
ETree_msStats (
   ETree   *frontETree,
   IV      *msIV,
   IV      *nvtxIV,
   IV      *nzfIV,
   DV      *opsDV,
   int     type,
   int     symlag
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in permute.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   fill the new-to-old permutation vector for the fronts
 
   created -- 96jun23, cca
   -----------------------------------------------------
*/
IV *
ETree_newToOldFrontPerm (
   ETree  *etree
) ;
/*
   -----------------------------------------------------
   fill the old-to-new permutation vector for the fronts
 
   created -- 96jun23, cca
   -----------------------------------------------------
*/
IV *
ETree_oldToNewFrontPerm (
   ETree  *etree
) ;
/*
   -------------------------------------------------------
   fill the new-to-old permutation vector for the vertices
 
   created -- 96jun23, cca
   -------------------------------------------------------
*/
IV *
ETree_newToOldVtxPerm (
   ETree  *etree
) ;
/*
   -------------------------------------------------------
   fill the old-to-new permutation vector for the vertices
 
   created -- 96jun23, cca
   -------------------------------------------------------
*/
IV *
ETree_oldToNewVtxPerm (
   ETree  *etree
) ;
/*
   -------------------------------------------------------
   purpose -- permute the vertices, 
              overwrite entries in the vertex-to-front map
 
   created -- 96oct03, cca
   -------------------------------------------------------
*/
void
ETree_permuteVertices (
   ETree   *etree,
   IV      *vtxOldToNewIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in compress.c -------------------------------------
------------------------------------------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- 
   to create and return an IV object that contains the map 
   from old to new fronts that are fundamental chains.

   created  -- 96jun23, cca
   -------------------------------------------------------
*/
IV *
ETree_fundChainMap (
   ETree   *etree
) ;
/*
   -------------------------------------------------------
   purpose -- 
   to create and return an IV object that contains the map 
   from old to new fronts that are fundamental supernodes.

   created  -- 96jun23, cca
   -------------------------------------------------------
*/
IV *
ETree_fundSupernodeMap (
   ETree   *etree
) ;
/*
   -----------------------------------------------------------
   compress an ETree object given a map from old to new nodes.
   note, a new node must be a connected set of the old nodes.

   return value -- pointer to new ETree object

   created -- 96jun23, cca.
   -----------------------------------------------------------
*/
ETree *
ETree_compress (
   ETree   *etree,
   IV      *frontmapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in justify.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   left-justify a tree by subtree size
   children are linked in ascending order of their subtree size

   created -- 96jan11, cca
   ------------------------------------------------------------
*/
void
ETree_leftJustify (
   ETree   *etree
) ;
/*
   ------------------------------------------------------
   left-justify a etree by a metric
   children are linked in ascending order of their metric

   created -- 96jan11, cca
   ------------------------------------------------------
*/
void
ETree_leftJustifyI (
   ETree   *etree,
   IV      *metricIV
) ;
/*
   ------------------------------------------------------
   left-justify a etree by a metric
   children are linked in ascending order of their metric

   created -- 96jan11, cca
   ------------------------------------------------------
*/
void
ETree_leftJustifyD (
   ETree   *etree,
   DV      *metricDV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in metrics.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------
   return an IV object with the weights 
   of the vertices in each front.

   created -- 96jun23, cca
   ------------------------------------
*/
IV *
ETree_nvtxMetric (
   ETree   *etree
) ;
/*
   ---------------------------------------------------------------
   return an IV object with the number 
   of factor entries in each front.
 
   symflag -- symmetryflag 
      SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   created -- 96jun23, cca
   ---------------------------------------------------------------
*/
IV *
ETree_nentMetric (
   ETree   *etree,
   int     flag
) ;
/*
   ---------------------------------------------------------------
   return a DV object with the number
   of factor operations in each front.
 
   type -- type of entries,
      SPOOLES_REAL or SPOOLES_COMPLEX
 
   symflag -- symmetryflag,
      SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   created -- 96jun23, cca
   ---------------------------------------------------------------
*/
DV *
ETree_nopsMetric (
   ETree   *etree,
   int     type,
   int     symflag
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in stages.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   generate a stages vector to be used by the CMD object.
   (1) if v = par[u] then
          stages[v] >= stages[u]
       endif
   (2) if v is a leaf then
          stages[v] = 0
       endif
   (3) if u and v belong to the same fundamental supernode then
          stages[v] = stages[u]
       endif
 
   basically, all nodes in a domain (a subtree) have stage zero,
   and the stages of all fundamental supernodes ancestor to that
   subtree are distinct.
 
   input --
 
      msIV -- IV object that contains the vertices in the
              multisector (non-domain vertices)
 
   return value --
 
      stagesIV -- an IV object that contains the stage for each vertex
 
   created -- 96feb19, cca
   -------------------------------------------------------------------
*/
IV *
ETree_stagesViaMS (
   ETree    *etree,
   IV       *msIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in transform.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------
   transform an ETree object by
   (1) merging small fronts into larger fronts
       using the ETree_mergeFrontsOne() method
   (2) merging small fronts into larger fronts
       using the ETree_mergeFrontsAll() method
   (3) merging small fronts into larger fronts
       using the ETree_mergeFrontsAny() method
   (4) split a large front into a chain of smaller fronts
       using the ETree_splitFronts() method
 
   created  -- 96jun27, cca
   ------------------------------------------------------
*/
ETree *
ETree_transform (
   ETree   *etree,
   int     vwghts[],
   int     maxzeros,
   int     maxfrontsize,
   int     seed
) ;
/*
   ------------------------------------------------------
   transform an ETree object by 
   (1) merging small fronts into larger fronts
       using the ETree_mergeFrontsOne() method
   (2) merging small fronts into larger fronts
       using the ETree_mergeFrontsAll() method
   (3) split a large front into a chain of smaller fronts
       using the ETree_splitFronts() method

   created  -- 96jun27, cca
   ------------------------------------------------------
*/
ETree *
ETree_transform2 (
   ETree   *etree,
   int     vwghts[],
   int     maxzeros,
   int     maxfrontsize,
   int     seed
) ;
/*
   --------------------------------------------------------------------
   purpose -- merge the front tree allowing only chains of nodes to
     merge that create at most maxzeros zero entries inside a front

   return --
      IV object that has the old front to new front map
 
   created -- 98jan29, cca
   --------------------------------------------------------------------
*/
ETree *
ETree_mergeFrontsOne (
   ETree   *etree,
   int     maxzeros,
   IV      *nzerosIV
) ;
/*
   -------------------------------------------------------
   purpose -- merge the front tree allowing a parent
              to absorb all children when that creates
              at most maxzeros zero entries inside a front
 
   return --
      IV object that has the old front to new front map
 
   created -- 98jan29, cca
   -------------------------------------------------------
*/
ETree *
ETree_mergeFrontsAll (
   ETree   *etree,
   int     maxzeros,
   IV      *nzerosIV
) ;
/*
   --------------------------------------------------------------------
   purpose -- merge the front tree allowing at most
              maxzeros zero entries inside a front
 
   return -- 
      IV object that has the old front to new front map
 
   created -- 96jun23, cca
   modified -- 97dec18, cca
      bug fixed that incorrectly counted the number of zeros in a front
   --------------------------------------------------------------------
*/
ETree *
ETree_mergeFrontsAny (
   ETree   *etree,
   int     maxzeros,
   IV      *nzerosIV
) ;
/*
   -------------------------------------------------
   expand an ETree object by splitting a large front
   into a chain of smaller fronts.
 
   created -- 96jun27, cca
   -------------------------------------------------
*/
ETree *
ETree_splitFronts (
   ETree   *etree,
   int     vwghts[],
   int     maxfrontsize,
   int     seed
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in maps.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   this method constructs and returns an IV object that holds the
   map from fronts to threads for a wrap map of the front tree.
 
   created -- 96dec12, cca
   --------------------------------------------------------------
*/
IV *
ETree_wrapMap (
   ETree   *frontTree,
   int     type,
   int     symflag,
   DV      *cumopsDV
) ;
/*
   ----------------------------------------------------------------
  this method constructs and returns an IV object that holds the
   map from fronts to threads for a balanced map of the front tree.
  the fronts are visited in the post-order traversal.
 
   created -- 96dec12, cca
   ----------------------------------------------------------------
*/
IV *
ETree_balancedMap (
   ETree   *frontTree,
   int     type,
   int     symflag,
   DV      *cumopsDV
) ;
/*
   -----------------------------------------------
   this method constructs and returns an IV object
   that holds the map from fronts to threads for a
   "subtree-subset" map of the front tree.
 
   created -- 97jan15, cca
   -----------------------------------------------
*/
IV *
ETree_subtreeSubsetMap (
   ETree   *frontTree,
   int     type,
   int     symflag,
   DV      *cumopsDV
) ;
/*
   ----------------------------------------------------------------
   this method constructs and returns an IV object that holds the
   map from fronts to threads for a domain decomposition 
   balanced map of the front tree.
   the domains are mapped to threads using a balanced map,
   and the schur complement fronts are mapped to threads 
   using a balanced map, but the two balanced maps are independent.

   created -- 97jan17, cca
   ----------------------------------------------------------------
*/
IV *
ETree_ddMap (
   ETree    *frontTree,
   int      type,
   int      symflag,
   DV       *cumopsDV,
   double   cutoff
) ;
/*
   ----------------------------------------------------------------
   this method constructs and returns an IV object that holds the
   map from fronts to threads for a domain decomposition 
   balanced map of the front tree.
   the domains are mapped to threads using a balanced map,
   and the schur complement fronts are mapped to threads 
   using a balanced map, but the two balanced maps are independent.

   created -- 97jan17, cca
   ----------------------------------------------------------------
*/
IV *
ETree_ddMapNew (
   ETree   *frontTree,
   int     type,
   int     symflag,
   IV      *msIV,
   DV      *cumopsDV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in splice.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   this method is used to splice together two front trees
   when the domain vertices and schur complement vertices
   have been ordered separately.

   etree0 -- the lower front tree is for vertices in the domain.
   graph0 -- graph for all the vertices
   mapIV  -- IV object that maps vertices to schur complement
             vertices, if IV_entry(mapIV, v) < 0 then v is 
             a domain vertex.
   etree1 -- the upper front tree is for vertices in the schur 
             complement.

   created -- 97feb01, cca
   -------------------------------------------------------------
*/
ETree *
ETree_spliceTwoETrees (
   ETree   *etree0,
   Graph   *graph0,
   IV      *mapIV,
   ETree   *etree1
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in initFromSubtree.c ------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- to initialize subtree with the subtree 
              of the front tree using nodes in nodeidsIV.
              vtxIV is filled with the vertices in the subtree
 
   return values ---
      1 -- normal return
     -1 -- subtree is NULL
     -2 -- nodeidsIV is NULL
     -3 -- etree is NULL
     -4 -- nodeidsIV is invalid
     -5 -- vtxIV is NULL
 
   created -- 98oct15, cca
   -----------------------------------------------------------
*/
int
ETree_initFromSubtree (
   ETree   *subtree,
   IV      *nodeidsIV,
   ETree   *etree,
   IV      *vtxIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in semi.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   purpose -- 
 
   to find the optimal domain/schur complement partition
   for a semi-implict factorization.
 
   the gain of a subtree sbt(J) is equal to
 
   |L_{bnd{J},sbt{J}}| - |A_{bnd{J},sbt{J}}|
      - alpha *|L_{sbt{J},sbt{J}}| 
 
   when alpha = 0 we minimize active storage
   when alpha = 1 we minimize solve operations
 
   *ptotalgain is filled with the total gain
 
   the return value is compidsIV,
      compids[J] = 0 --> J is in the schur complement
      compids[J] != 0 --> J is in domain compids[J]
 
   created -- 98jun20, cca
   -----------------------------------------------------
*/
IV *
ETree_optPart (
   ETree    *etree,
   Graph    *graph,
   IVL      *symbfacIVL,
   double   alpha, 
   int      *ptotalgain,
   int      msglvl,
   FILE     *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in storage.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the active storage to eliminate J
               using the multifrontal method
 
   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC,
              SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   created -- 97may21, cca
   ---------------------------------------------------------------
*/
void
ETree_MFstackProfile (
   ETree    *etree,
   int      symflag,
   double   dvec[]
) ;
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the active storage to eliminate J
               using the left-looking general sparse method
 
   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC,
              SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   created -- 97may21, cca
   ---------------------------------------------------------------
*/
void
ETree_GSstorageProfile (
   ETree    *etree,
   int      symflag,
   IVL      *symbfacIVL,
   int      *vwghts,
   double   dvec[]
) ;
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the active storage to eliminate J
               using the right-looking general sparse method

   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC,
              SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 98dec19, cca
   ---------------------------------------------------------------
*/
void
ETree_FSstorageProfile (
   ETree    *etree,
   int      symflag,
   IVL      *symbfacIVL,
   double   dvec[]
) ;
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the stack storage to solve for J
               in a forward solve
 
   created -- 97nov30, cca
   ---------------------------------------------------------------
*/
void
ETree_forwSolveProfile (
   ETree    *etree,
   double   dvec[]
) ;
/*
   ---------------------------------------------------------------
   purpose --  fill dvec[J] with the stack storage to solve for J
               in a backward solve
 
   created -- 97nov30, cca
   ---------------------------------------------------------------
*/
void
ETree_backSolveProfile (
   ETree    *etree,
   double   dvec[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in IO.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   purpose -- to read in an ETree object from a file

   input --

      fn -- filename, must be *.etreeb or *.etreef

   return value -- 1 if success, 0 if failure

   created -- 95nov15, cca
   -------------------------------------------------
*/
int
ETree_readFromFile ( 
   ETree   *etree, 
   char    *fn 
) ;
/*
   --------------------------------------------------------
   purpose -- to read an ETree object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95nov15, cca
   --------------------------------------------------------
*/
int
ETree_readFromFormattedFile ( 
   ETree   *etree, 
   FILE    *fp 
) ;
/*
   ----------------------------------------------------
   purpose -- to read an ETree object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95nov15, cca
   ----------------------------------------------------
*/
int
ETree_readFromBinaryFile ( 
   ETree    *etree, 
   FILE   *fp 
) ;
/*
   --------------------------------------------
   purpose -- to write an ETree object to a file

   input --

      fn -- filename
        *.etreeb -- binary
        *.etreef -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   --------------------------------------------
*/
int
ETree_writeToFile ( 
   ETree   *etree, 
   char   *fn 
) ;
/*
   ------------------------------------------------------
   purpose -- to write an ETree object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ------------------------------------------------------
*/
int
ETree_writeToFormattedFile ( 
   ETree   *etree, 
   FILE    *fp 
) ;
/*
   ---------------------------------------------------
   purpose -- to write an ETree object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ---------------------------------------------------
*/
int
ETree_writeToBinaryFile ( 
   ETree    *etree, 
   FILE   *fp 
) ;
/*
   ---------------------------------------------------
   purpose -- to write an ETree object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ---------------------------------------------------
*/
int
ETree_writeForHumanEye ( 
   ETree    *etree, 
   FILE   *fp 
) ;
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the ETree object

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   -----------------------------------------------------------
*/
int
ETree_writeStats ( 
   ETree    *etree, 
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
