/*  GPart.h  */

#include "../Graph.h"
#include "../BPG.h"
#include "../DSTree.h"
#include "../Tree.h"
#include "../IV.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   The GPart object is used to generate and store information
   about a partition of the graph.

   id      -- id for the object
   g       -- pointer to Graph object, not free'd by the destructor
   nvtx    -- # of vertices in the graph, internal vertices only
   nvbnd   -- # of vertices in the boundary of the graph, 
   ncomp   -- # of components in the graph
   compidsIV -- IV object that contains a map from
                vertices to components, size nvtx
      compids[v] == 0 --> v is on the interface
      compids[v] != 0 --> v belongs to domain compids[v]
   cweightsIV -- IV object that contains the component weights
   par        -- pointer to parent GPart object
   fch        -- pointer to first child GPart object
   sib        -- pointer to sibling GPart object
   vtxMapIV   -- IV object that contains a map map from local vertices 
                 to parent's vertices or global vertices
   msglvl     -- message level, default is zero
   msgFile    -- message file pointer, default is stdout

   created  -- 95oct21, cca
   ----------------------------------------------------------------
*/
typedef struct _GPart  GPart ;
struct _GPart {
   int     id         ;
   Graph   *g         ;
   int     nvtx       ;
   int     nvbnd      ;
   int     ncomp      ;
   IV      compidsIV  ;
   IV      cweightsIV ;
   GPart   *par       ;
   GPart   *fch       ;
   GPart   *sib       ;
   IV      vtxMapIV   ;
   int     msglvl     ;
   FILE    *msgFile   ;
} ;
/*
   ------------------------------------------
   include the DDsepInfo object's header file
   ------------------------------------------
*/
#include "DDsepInfo.h"
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c   --------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------
   construct a new instance of the GPart object

   created -- 95oct05, cca
   --------------------------------------------
*/
GPart *
GPart_new (
   void
) ;
/*
   ---------------------------------------------
   set the default fields of the GPart object

   created  -- 95oct05, cca
   modified -- 95nov29, cca
      par, fch, sib and vtxMap fields included
   ---------------------------------------------
*/
void
GPart_setDefaultFields (
   GPart   *gpart
) ;
/*
   ---------------------------------------------
   clear the data fields for a GPart object

   created  -- 95oct05, cca
   modified -- 95nov29, cca
      par, fch, sib and vtxMap fields included
   ---------------------------------------------
*/
void
GPart_clearData (
   GPart   *gpart
) ;
/*
   ------------------------
   free the GPart object

   created  -- 95oct05, cca
   modified -- 95nov29, cca
      gpart now free'd
   ------------------------
*/
void
GPart_free (
   GPart   *gpart
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in DDviaFishnet.c   --------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------
   purpose -- to construct and return a multisector
              using the fishnet algorithm

   freeze    -- parameter used to keep nodes of high degree
                in the interface
   minweight -- minimum weight for a domain
   maxweight -- maximum weight for a domain
   seed      -- random number seed

   created -- 96feb16, cca
   ---------------------------------------------------------------
*/
void
GPart_DDviaFishnet (
   GPart    *gpart,
   double   freeze,
   int      minweight,
   int      maxweight,
   int      seed
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in DDviaProjection.c   -----------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------
   set the compids[] vector using a global map from vertices
   to domains and interface nodes.

   DDmapIV -- IV object that contains the map from vertices
              to domains and interface nodes

   created -- 96mar17, cca
   ---------------------------------------------------------
*/
void
GPart_DDviaProjection (
   GPart   *gpart,
   IV      *DDmapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in DDsepInfo.c   -----------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   construct a new instance of the DDsepInfo object

   created -- 96feb24, cca
   ------------------------------------------------
*/
DDsepInfo *
DDsepInfo_new (
   void
) ;
/*
   ---------------------------------------------
   set the default fields of the DDsepInfo object

   created  -- 96feb24, cca
   ---------------------------------------------
*/
void
DDsepInfo_setDefaultFields (
   DDsepInfo   *info
) ;
/*
   ---------------------------------------------
   clear the data fields for a DDsepInfo object

   created  -- 96feb24, cca
   ---------------------------------------------
*/
void
DDsepInfo_clearData (
   DDsepInfo   *info
) ;
/*
   ------------------------
   free the DDsepInfo object

   created  -- 96feb24, cca
   ------------------------
*/
void
DDsepInfo_free (
   DDsepInfo   *info
) ;
/*
   -----------------------
   write the CPU times
  
   created -- 97nov06, cca
   -----------------------
*/
void
DDsepInfo_writeCpuTimes (
   DDsepInfo   *info,
   FILE        *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in domSegMap.c   -----------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   fill *pndom with ndom, the number of domains.
   fill *pnseg with nseg, the number of segments.
   domains are numbered in [0, ndom), segments in [ndom,ndom+nseg).

   return -- an IV object that contains the map 
             from vertices to segments

   created -- 99feb29, cca
   --------------------------------------------------------------------
*/
IV *
GPart_domSegMap (
   GPart   *gpart,
   int     *pndom,
   int     *pnseg
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in identifyWideSep.c   -----------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   identify the wide separator
 
   return -- IV object that holds the nodes in the wide separator

   created -- 96oct21, cca
   --------------------------------------------------------------
*/
IV *
GPart_identifyWideSep (
   GPart   *gpart,
   int     nlevel1,
   int     nlevel2
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c   ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------
   initialize the GPart object

   created -- 95oct05, cca
   ---------------------------
*/
void
GPart_init (
   GPart   *gpart,
   Graph   *g
) ;
/*
   -----------------------
   set the message fields

   created -- 96oct21, cca
   -----------------------
*/
void 
GPart_setMessageInfo (
   GPart   *gpart,
   int     msglvl,
   FILE    *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in makeYCmap.c   -----------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------
   make the map from wide separator vertices Y 
   to components {0, 1, 2, 3}.

   YCmap[y] == 0 --> y is not adjacent to either component
   YCmap[y] == 1 --> y is adjacent to only component 1
   YCmap[y] == 2 --> y is adjacent to only component 2
   YCmap[y] == 3 --> y is adjacent to components 1 and 2

   created -- 96jun09, cca
   -------------------------------------------------------
*/
IV *
GPart_makeYCmap (
   GPart   *gpart,
   IV      *YVmapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in RBviaDDsep.c   ----------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------
   this method constructs recursive partition of a graph.
   it returns a DSTree object to represent the splitting of
   the tree into subgraphs.
   the info object contains the information needed by the
   DDsep algorithm.

   created -- 96feb24, cca
   --------------------------------------------------------
*/
DSTree *
GPart_RBviaDDsep (
   GPart       *gpart,
   DDsepInfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in smoothBisector.c   ------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------
   smooth a wide separator
 
   nlevel -- number of levels one each side of the separator
             to include into the separator
   alpha  -- partition cost function parameter
   
   created -- 96jun02, cca
   ---------------------------------------------------------
*/
float
GPart_smoothBisector (
   GPart   *gpart,
   int     nlevel,
   float   alpha
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in smoothBy2layers.c   -----------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   smooth a bisector using the alternating two-layer algorithm

   option -- network flag
      1 --> make the network bipartite as for the
            Dulmage-Mendelsohn decomposition
      otherwise -- use network induced by the wide separator
   alpha -- cost function parameter

   created -- 96jun08, cca
   -----------------------------------------------------------
*/
void
GPart_smoothBy2layers (
   GPart   *gpart,
   int     option,
   float   alpha
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in smoothYSep.c   ----------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   smooth the wide separator given by YVmapIV by forming 
   the associated network, solving the max flow problem, 
   and examining two min-cuts.

   YVmapIV -- map from wide vertices Y to vertices V
   YCmapIV -- map from wide vertices Y to {0,1,2,3}
     YCmap[y] = 0 --> treat y as an internal node
                      adjacent to neither component
     YCmap[y] = 1 --> treat y as adjacent only to component 1
     YCmap[y] = 2 --> treat y as adjacent only to component 2
     YCmap[y] = 3 --> treat y as adjacent to components 1 and 2
   alpha -- cost function parameter

   created -- 96jun08, cca
   ------------------------------------------------------------
*/
float
GPart_smoothYSep (
   GPart   *gpart,
   IV      *YVmapIV,
   IV      *YCmapIV,
   float   alpha
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in split.c   ---------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------
   split the graph partition object into pieces

   created -- 95nov29, cca
   --------------------------------------------
*/
void
GPart_split (
   GPart   *gpart
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in TwoSetViaBKL.c   --------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------------
   given a domain decomposition, find a bisector
   1. construct the domain/segment graph
   2. use block kernihan-lin to get an initial bisector

   alpha   -- cost function parameter for BKL
   seed    -- random number seed
   cpus    -- array to store CPU times
              cpus[0] -- time to find domain/segment map
              cpus[1] -- time to find domain/segment bipartite graph
              cpus[2] -- time to find two-set partition

   return value -- cost of the partition

   created  -- 96mar09, cca
   -----------------------------------------------------------------
*/
double
GPart_TwoSetViaBKL (
   GPart       *gpart,
   double      alpha,
   int         seed,
   double      cpus[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c   ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------
   set the component weights from the compids[] vector

   created  -- 95oct05, cca
   modified -- 95nov29, cca
   ----------------------------------------------------
*/
void
GPart_setCweights (
   GPart   *gpart
) ;
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95oct05, cca
   ----------------------------------------------
*/
int 
GPart_sizeOf (
   GPart   *gpart
) ;
/*
   ----------------------------------------------------
   return 1 if vertex is adjacent to only one domain
            and fill *pdomid with the domain's id
   return 0 otherwise
 
   created -- 95oct19, cca
   -------------------------------------------------
*/
int
GPart_vtxIsAdjToOneDomain (
   GPart   *gpart,
   int     v,
   int     *pdomid
) ;
/*
   ------------------------------------------------------
   return 1 if the partition has a valid vertex separator
          0 otherwise

   created -- 95oct18, cca
   ------------------------------------------------------
*/
int
GPart_validVtxSep (
   GPart   *gpart
) ;
/*
   -------------------------------------
   return an IV object filled with the
   weights of the component's boundaries

   created -- 96oct21, cca
   -------------------------------------
*/
IV *
GPart_bndWeightsIV (
   GPart   *gpart 
) ;
/*--------------------------------------------------------------------*/
