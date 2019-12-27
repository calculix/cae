/*  EGraph.h  */

#include "../cfiles.h"
#include "../Graph.h"
#include "../IVL.h"

/*====================================================================*/
/*
   ---------------------------------------------------------------
   the EGraph object stores an element graph

   type   -- type of element graph
      0 --> vertices have unit weight
      1 --> vertices are weighted
   nelem  -- # of elements
   nvtx   -- # of vertices
   adj    -- IVL object that holds the lists 
             of elements and their vertices
   vwghts -- pointer to vector of vertex weights, NULL if type = 0

   created -- 95nov03, cca
   ---------------------------------------------------------------
*/
typedef struct _EGraph {
   int   type    ;
   int   nelem   ;
   int   nvtx    ;
   IVL   *adjIVL ;
   int   *vwghts ;
} EGraph ;

/*====================================================================*/
/*
------------------------------------------------------------------------
------ methods found in basics.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   constructor
 
   created -- 95nov03, cca
   -----------------------
*/
EGraph *
EGraph_new (
   void
) ;
/*
   -----------------------
   set the default fields
 
   created -- 95nov03, cca
   -----------------------
*/
void
EGraph_setDefaultFields (
   EGraph   *eg
) ;/*
   -----------------------
   clear the data fields
 
   created -- 95nov03, cca
   -----------------------
*/
void
EGraph_clearData (
   EGraph   *eg
) ;
/*
   -----------------------
   destructor
 
   created -- 95nov03, cca
   -----------------------
*/
void
EGraph_free (
   EGraph   *eg
) ;
/*====================================================================*/
/*
------------------------------------------------------------------------
------ methods found in init.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------
   purpose -- initialize the EGraph object
 
   created -- 96oct24, cca
   ---------------------------------------
*/
void
EGraph_init (
   EGraph   *egraph,
   int      type,
   int      nelem,
   int      nvtx,
   int      IVL_type
) ;
/*====================================================================*/
/*
------------------------------------------------------------------------
------ methods found in mkAdjGraph.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------
   create a Graph object that holds the adjacency graph 
   of the assembled elements. 
 
   created -- 95nov03, cca
   ----------------------------------------------------
*/
Graph *
EGraph_mkAdjGraph ( 
   EGraph   *eg 
) ;
/*====================================================================*/
/*
------------------------------------------------------------------------
------ methods found in misc.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   make an element graph for a n1 x n2 grid with ncomp components
 
   created -- 95nov03, cca
   --------------------------------------------------------------
*/
EGraph *
EGraph_make9P ( 
   int   n1, 
   int   n2, 
   int   ncomp 
) ;
/*

-------------------------------------------------------------------
  make an element graph for a n1 x n2 x n3 grid with ncomp
components
 
   created -- 95nov03, cca

-------------------------------------------------------------------
*/
EGraph *
EGraph_make27P (
   int   n1,
   int   n2,
   int   n3,
   int   ncomp
) ;
/*====================================================================*/
/*
   ------------------------
   IO methods found in IO.c
   -----------------------
*/
/*
   -------------------------------------------------
   purpose -- to read in a EGraph object from a file
 
   input --
 
      fn -- filename, must be *.egraphb or *.egraphf
 
   return value -- 1 if success, 0 if failure
 
   created -- 95nov03, cca
   -------------------------------------------------
*/
int
EGraph_readFromFile ( 
   EGraph   *egraph, 
   char     *fn 
) ;
/*
   --------------------------------------------------------
   purpose -- to read a EGraph object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 95nov03, cca
   --------------------------------------------------------
*/
int
EGraph_readFromFormattedFile (
   EGraph   *egraph,
   FILE     *fp
) ;
/*
   ----------------------------------------------------
   purpose -- to read a EGraph object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 95nov03, cca
   ----------------------------------------------------
*/
int
EGraph_readFromBinaryFile (
   EGraph   *egraph,
   FILE    *fp
) ;
/*
   --------------------------------------------
   purpose -- to write a EGraph object to a file
 
   input --
 
      fn -- filename
        *.egraphb -- binary
        *.egraphf -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov03, cca
   --------------------------------------------
*/
int
EGraph_writeToFile (
   EGraph   *egraph,
   char     *fn
) ;
/*
   ------------------------------------------------------
   purpose -- to write a EGraph object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov03, cca
   ------------------------------------------------------
*/
int
EGraph_writeToFormattedFile (
   EGraph   *egraph,
   FILE     *fp
) ;
/*
   ---------------------------------------------------
   purpose -- to write a EGraph object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov03, cca
   ---------------------------------------------------
*/
int
EGraph_writeToBinaryFile (
   EGraph   *egraph,
   FILE     *fp
) ;
/*
   -------------------------------------------------
   purpose -- to write a EGraph object for a human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov03, cca
   -------------------------------------------------
*/
int
EGraph_writeForHumanEye (
   EGraph    *egraph,
   FILE   *fp
) ;
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the EGraph object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95nov03, cca
   -----------------------------------------------------------
*/
int
EGraph_writeStats (
   EGraph    *egraph,
   FILE   *fp
) ;
/*====================================================================*/
/*
   -------------------------------------
   miscellaneous methods found in misc.c
   -------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   return the element graph for a 9 point operator on a 
   n1 x n2 grid with ncomp components at each grid point
   -----------------------------------------------------
*/
EGraph *
EGraph_make9P ( int n1, int n2, int ncomp ) ;
/*
   ----------------------------------------------------------
   return the element graph for a 27 point operator on a 
   n1 x n2 x n3 grid with ncomp components at each grid point
   ----------------------------------------------------------
*/
EGraph *
EGraph_make27P ( int n1, int n2, int n3, int ncomp ) ;

/*--------------------------------------------------------------------*/
