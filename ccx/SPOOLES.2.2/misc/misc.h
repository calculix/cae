/*  misc.h  */

#include "../cfiles.h"
#include "../DenseMtx.h"
#include "../InpMtx.h"
#include "../ETree.h"
#include "../Coords.h"
#include "../Graph.h"
#include "../MSMD.h"
#include "../GPart.h"

/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in ND.c --------------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   purpose -- main procedure. if the input region is 1 x 1 x 1,
              the node is put into the permutation vector.
              otherwise the region is split into three pieces,
              two subregions and a separator, and recursive 
              calls are made to order the subregions

   input --

      n1         -- number of points in the first direction
      n2         -- number of points in the second direction
      n3         -- number of points in the third direction
      new_to_old -- pointer to the permutation vector
      west       -- west coordinate
      east       -- east coordinate
      south      -- south coordinate
      north      -- north coordinate
      bottom     -- bottom coordinate
      top        -- top coordinate

   created -- 95nov15, cca
   ------------------------------------------------------------
*/
void
mkNDperm ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   new_to_old[], 
   int   west, 
   int   east, 
   int   south, 
   int   north, 
   int   bottom, 
   int   top 
) ;
/*
   -----------------------------------------------------
   purpose -- to print a vector on a 2-d grid

   input --

      n1   -- number of points in the first direction
      n2   -- number of points in the second direction
      ivec -- integer vector to be printed in %4d format
      fp   -- file pointer

   created -- 95nov16, cca
   -----------------------------------------------------
*/
void
fp2DGrid ( 
   int    n1, 
   int    n2, 
   int    ivec[], 
   FILE   *fp 
) ;
/*
   -----------------------------------------------------
   purpose -- to print a vector on a 3-d grid

   input --

      n1   -- number of points in the first direction
      n2   -- number of points in the second direction
      n3   -- number of points in the third direction
      ivec -- integer vector to be printed in %4d format
      fp   -- file pointer

   created -- 95nov16, cca
   -----------------------------------------------------
*/
void
fp3DGrid ( 
   int    n1, 
   int    n2, 
   int    n3, 
   int    ivec[], 
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in localND.c ---------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   compute an old-to-new ordering for 
   local nested dissection in two dimensions

   n1       -- number of grid points in first direction
   n2       -- number of grid points in second direction
   p1       -- number of domains in first direction
   p2       -- number of domains in second direction
   dsizes1  -- domain sizes in first direction, size p1
               if NULL, then we construct our own
   dsizes2  -- domain sizes in second direction, size p2
               if NULL, then we construct our own
   oldToNew -- old-to-new permutation vector

   note : the following must hold
      n1 > 0, n2 >0, n1 >= 2*p1 - 1, n2 >= 2*p2 - 1, p2 > 1
      sum(dsizes1) = n1 - p1 + 1 and sum(dsizes2) = n2 - p2 + 1

   created -- 95nov16, cca
   ------------------------------------------------------------
*/
void
localND2D ( 
   int   n1, 
   int   n2, 
   int   p1, 
   int   p2, 
   int   dsizes1[], 
   int   dsizes2[], 
   int   oldToNew[] 
) ;
/*
   ------------------------------------------------------------
   compute an old-to-new ordering for 
   local nested dissection in three dimensions

   n1       -- number of grid points in first direction
   n2       -- number of grid points in second direction
   n3       -- number of grid points in third direction
   p1       -- number of domains in first direction
   p2       -- number of domains in second direction
   p3       -- number of domains in third direction
   dsizes1  -- domain sizes in first direction, size p1
               if NULL, then we construct our own
   dsizes2  -- domain sizes in second direction, size p2
               if NULL, then we construct our own
   dsizes3  -- domain sizes in third direction, size p3
               if NULL, then we construct our own
   oldToNew -- old-to-new permutation vector

   note : the following must hold
      n1 > 0, n2 >0, n3 > 0,
      n1 >= 2*p1 - 1, n2 >= 2*p2 - 1, n3 >= 2*p3 - 1, p3 > 1
      sum(dsizes1) = n1 - p1 + 1, sum(dsizes2) = n2 - p2 + 1
      sum(dsizes3) = n3 - p3 + 1

   created -- 95nov16, cca
   ------------------------------------------------------------
*/
void
localND3D ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   p1, 
   int   p2, 
   int   p3,
   int   dsizes1[], 
   int   dsizes2[], 
   int   dsizes3[], 
   int   oldToNew[] 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in orderViaBestOfNDandMS.c -------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- return an ETree object for the better of a 
              nested dissection and multisection orderings
 
   graph -- graph to order
   maxdomainsize -- used to control the incomplete nested dissection 
     process. any subgraph whose weight is less than maxdomainsize
     is not split further.
   maxzeros -- maximum number of zero entries allowed in a front
   maxsize  -- maximum number of internal columns in a front
   seed     -- random number seed
   msglvl  -- message level
      1 -- timings and statistics
      2 -- more timings and statistics
      3 -- lots of output
   msgFile -- message file
 
   created -- 98aug26, cca
   ------------------------------------------------------------------
*/
ETree *
orderViaBestOfNDandMS (
   Graph    *graph,
   int      maxdomainsize,
   int      maxzeros,
   int      maxsize,
   int      seed,
   int      msglvl,
   FILE     *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in orderViaND.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- return an ETree object for a nested dissection ordering
 
   graph -- graph to order
   maxdomainsize -- used to control the incomplete nested dissection 
     process. any subgraph whose weight is less than maxdomainsize
    is not split further.
   seed    -- random number seed
   msglvl  -- message level, 0 --> no output, 1 --> timings
   msgFile -- message file
 
   created -- 97nov06, cca
   ------------------------------------------------------------------
*/
ETree *
orderViaND (
   Graph   *graph,
   int     maxdomainsize,
   int     seed,
   int     msglvl,
   FILE    *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in orderViaMS.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- return an ETree object for a multisection ordering
 
   graph -- graph to order
   maxdomainsize -- used to control the incomplete nested dissection 
     process. any subgraph whose weight is less than maxdomainsize
    is not split further.
   seed    -- random number seed
   msglvl  -- message level, 0 --> no output, 1 --> timings
   msgFile -- message file
 
   created -- 97nov06, cca
   ------------------------------------------------------------------
*/
ETree *
orderViaMS (
   Graph   *graph,
   int     maxdomainsize,
   int     seed,
   int     msglvl,
   FILE    *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in orderViaMMD.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------
   purpose -- return an ETree object for 
              a multiple minimum degree ordering
 
   graph -- graph to order
   seed    -- random number seed
   msglvl  -- message level, 0 --> no output, 1 --> timings
   msgFile -- message file
 
   created -- 97nov08, cca
   --------------------------------------------------------
*/
ETree *
orderViaMMD (
   Graph   *graph,
   int     seed,
   int     msglvl,
   FILE    *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in drawGraphEPS.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   draw a graph to an EPS file
 
   (1) read a Graph object
   (2) read a Coords object
   (3) read an IV object that contains a tag vector.
       if (v,w) is an edge in the graph
       and tags[v] == tags[w] then
          draw edge (v,w)
   (4) bbox[4] is the bounding box for the plot.
       bbox = { xsw, ysw, xne, yne }
       try bbox = { 0, 0, 500, 500 }
       because coordinates are measured in points,
       72 points per inch.
   (5) rect[4] contains the frame for the plot.
       to put a 20 point margin around the plot,
       rect[0] = bbox[0] + 20
       rect[1] = bbox[1] + 20
       rect[2] = bbox[2] - 20
       rect[3] = bbox[3] - 20
 
   created -- 98apr11, cca
   -------------------------------------------------
*/
void
drawGraphEPS (
   Graph    *graph,
   Coords   *coords,
   IV       *tagsIV,
   double   bbox[],
   double   rect[],
   double   linewidth1,
   double   linewidth2,
   double   radius,
   char     *epsFileName,
   int      msglvl,
   FILE     *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in mkNDlinsys.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------------
   purpose -- to create a linear system A*X = B
      for a nested dissection ordering of a 2-d or 3-d regular grid
   input --
 
      n1 -- # of nodes in first direction
      n2 -- # of nodes in second direction
      n3 -- # of nodes in third direction
      maxzeros -- relaxation factor for fronts,
         maximum number of zero entries in a front
      maxsize  -- split parameter for large fronts,
         maximum number of internal vertices in a front
      type -- type of entries
         SPOOLES_REAL or SPOOLES_COMPLEX
      symmetryflag -- symmetry of the matrix
         SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
      nrhs    -- number of right hand sides
      seed    -- seed for random number generator
      msglvl  -- message level
      msgFile -- message file
 
   output --
 
      pfrontETree -- to be filled with address of front tree 
      psymbfacIVL -- to be filled with address of symbolic factorization
      pmtxA       -- to be filled with address of matrix object A
      pmtxX       -- to be filled with address of matrix object X
      pmtxB       -- to be filled with address of matrix object B
 
   created -- 98may16, cca
   ---------------------------------------------------------------------
*/
void
mkNDlinsys (
   int        n1,
   int        n2,
   int        n3,
   int        maxzeros,
   int        maxsize,
   int        type,
   int        symmetryflag,
   int        nrhs,
   int        seed,
   int        msglvl,
   FILE       *msgFile,
   ETree      **pfrontETree,
   IVL        **psymbfacIVL,
   InpMtx     **pmtxA,
   DenseMtx   **pmtxX,
   DenseMtx   **pmtxB
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods in mkNDlinsysQR.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
  ---------------------------------------------------------------------
   purpose -- to create an overdetermined linear system A*X = B
      for a nested dissection ordering of a 2-d or 3-d regular grid
 
   input --
 
      n1 -- # of nodes in first direction
      n2 -- # of nodes in second direction
      n3 -- # of nodes in third direction
      type -- type of entries
         SPOOLES_REAL or SPOOLES_COMPLEX
      nrhs    -- number of right hand sides
      seed    -- seed for random number generator
      msglvl  -- message level
      msgFile -- message file
 
   output --
 
      pfrontETree -- to be filled with address of front tree 
      psymbfacIVL -- to be filled with address of symbolic factorization
      pmtxA       -- to be filled with address of matrix object A
      pmtxX       -- to be filled with address of matrix object X
      pmtxB       -- to be filled with address of matrix object B
 
   created -- 98may29, cca
   ---------------------------------------------------------------------
*/
void
mkNDlinsysQR (
   int        n1,
   int        n2,
   int        n3,
   int        type,
   int        nrhs,
   int        seed,
   int        msglvl,
   FILE       *msgFile,
   ETree      **pfrontETree,
   IVL        **psymbfacIVL,
   InpMtx     **pmtxA,
   DenseMtx   **pmtxX,
   DenseMtx   **pmtxB
) ;
/*--------------------------------------------------------------------*/
