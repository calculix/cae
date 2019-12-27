/*  DSTree.h  */

#include "../cfiles.h"
#include "../Tree.h"
#include "../IV.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   the DSTree object models a domain/separator tree.
   it contains a Tree object to represent the connectivity of
   the domains and separators and an IV object to contain the
   map from vertices to domains and separators.

   tree  -- pointer to a Tree object that contains the tree adjacency
   mapIV -- pointer to IV object that contains the map 
            from vertices to domains and separators

   created -- 96mar10, cca
   ---------------------------------------------------------------------
*/
typedef struct _DSTree   DSTree ;
struct _DSTree {
   Tree   *tree  ;
   IV     *mapIV ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in basics.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- create and return a new DSTree object

   created -- 96mar10, cca
   -----------------------------------------------
*/
DSTree *
DSTree_new ( 
   void
) ;
/*
   ------------------------------------------------------
   purpose -- set the default fields for the DSTree object

   created -- 96mar10, cca
   ------------------------------------------------------
*/
void
DSTree_setDefaultFields (
   DSTree   *dstree
) ;
/*
   --------------------------------
   purpose -- clear the data fields

   created -- 96mar10, cca
   --------------------------------
*/
void
DSTree_clearData (
   DSTree   *dstree
) ;
/*
   --------------------------------
   purpose -- free the DSTree object

   created -- 96mar10, cca
   --------------------------------
*/
void
DSTree_free (
   DSTree   *dstree
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in instance.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------
   return a pointer to the internal Tree object
 
   created -- 97jun21, cca
   --------------------------------------------
*/
Tree *
DSTree_tree (
   DSTree   *dstree
) ;
/*
   -------------------------------------
   return a pointer to the map IV object
 
   created -- 97jun21, cca
   -------------------------------------
*/
IV *
DSTree_mapIV (
   DSTree   *dstree
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
 
   created -- 96mar10, cca
   -----------------------------------------------
*/
void
DSTree_init1 (
   DSTree   *dstree,
   int      ndomsep,
   int      nvtx
) ;
/*
   -----------------------------------------
   initialize the object given a Tree object
 
   created -- 96mar10, cca
   -----------------------------------------
*/
void
DSTree_init2 (
   DSTree   *dstree,
   Tree     *tree,
   IV       *mapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in stages.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------
   create a stages vector for a nested dissection ordering

   created -- 96mar10, cca
   -------------------------------------------------------
*/
IV *
DSTree_NDstages (
   DSTree   *dstree
) ;
/*
   ----------------------------------------------------------------
   create a stages vector for a ``half'' nested dissection ordering

   created -- 96mar10, cca
   ----------------------------------------------------------------
*/
IV *
DSTree_ND2stages (
   DSTree   *dstree
) ;
/*
   ------------------------------------------------------------
   create a stages vector for a two-level multisection ordering

   created -- 96mar10, cca
   ------------------------------------------------------------
*/
IV *
DSTree_MS2stages (
   DSTree   *dstree
) ;
/*
   --------------------------------------------------------------
   create a stages vector for a three-level multisection ordering

   created -- 96mar10, cca
   --------------------------------------------------------------
*/
IV *
DSTree_MS3stages (
   DSTree   *dstree
) ;
/*
   ---------------------------------------------------
   create a stages vector via cutoff on domain weights
 
   created -- 96mar10, cca
   ---------------------------------------------------
*/
IV *
DSTree_stagesViaDomainWeight (
   DSTree   *dstree,
   int      *vwghts,
   DV       *cutoffDV
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
 
   created -- 96mar10, cca
   ----------------------------------------------
*/
int
DSTree_sizeOf (
   DSTree   *dstree
) ;
/*
   ---------------------------------------------
   renumber the fronts by a post-order traversal
   
   created -- 96apr13, cca
   ---------------------------------------------
*/
void
DSTree_renumberViaPostOT (
   DSTree * dstree 
) ;
/*
   -----------------------------------------------------------
   purpose -- return the weight of the vertices in the domains
 
   created -- 97jun21, cca
   -----------------------------------------------------------
*/
int
DSTree_domainWeight (
   DSTree   *dstree,
   int      vwghts[]
) ;
/*
   -----------------------------------------------------------
   purpose -- return the weight of the vertices in the domains
 
   created -- 97jun21, cca
   -----------------------------------------------------------
*/
int
DSTree_domainWeight (
   DSTree   *dstree,
   int      vwghts[]
) ;
/*
   --------------------------------------------------------------
   purpose -- return the weight of the vertices in the separators
 
   created -- 97jun21, cca
   --------------------------------------------------------------
*/
int
DSTree_separatorWeight (
   DSTree   *dstree,
   int      vwghts[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods founds in IO.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------
   purpose -- to read in an DSTree object from a file

   input --

      fn -- filename, must be *.dstreeb or *.dstreef

   return value -- 1 if success, 0 if failure

   created -- 96mar10, cca
   --------------------------------------------------
*/
int
DSTree_readFromFile ( 
   DSTree   *dstree, 
   char    *fn 
) ;
/*
   ---------------------------------------------------------
   purpose -- to read an DSTree object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 96mar10, cca
   ---------------------------------------------------------
*/
int
DSTree_readFromFormattedFile ( 
   DSTree   *dstree, 
   FILE    *fp 
) ;
/*
   ------------------------------------------------------
   purpose -- to read an DSTree object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 96mar10, cca
   ------------------------------------------------------
*/
int
DSTree_readFromBinaryFile ( 
   DSTree    *dstree, 
   FILE   *fp 
) ;
/*
   ----------------------------------------------
   purpose -- to write an DSTree object to a file

   input --

      fn -- filename
        *.dstreeb -- binary
        *.dstreef -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   ----------------------------------------------
*/
int
DSTree_writeToFile ( 
   DSTree   *dstree, 
   char   *fn 
) ;
/*
   --------------------------------------------------------
   purpose -- to write an DSTree object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   --------------------------------------------------------
*/
int
DSTree_writeToFormattedFile ( 
   DSTree   *dstree, 
   FILE    *fp 
) ;
/*
   -----------------------------------------------------
   purpose -- to write an DSTree object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   -----------------------------------------------------
*/
int
DSTree_writeToBinaryFile ( 
   DSTree    *dstree, 
   FILE   *fp 
) ;
/*
   ----------------------------------------------------
   purpose -- to write an DSTree object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   ----------------------------------------------------
*/
int
DSTree_writeForHumanEye ( 
   DSTree    *dstree, 
   FILE   *fp 
) ;
/*
   ------------------------------------------------------------
   purpose -- to write out the statistics for the DSTree object

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   ------------------------------------------------------------
*/
int
DSTree_writeStats ( 
   DSTree    *dstree, 
   FILE   *fp 
) ;
/*--------------------------------------------------------------------*/
