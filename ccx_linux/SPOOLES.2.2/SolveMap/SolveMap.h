/*  SolveMap.h  */

#include "../IVL.h"
#include "../Tree.h"
#include "../SPOOLES.h"

/*--------------------------------------------------------------------*/
typedef struct _SolveMap   SolveMap ;
struct _SolveMap {
   int   symmetryflag ;
   int   nfront       ;
   int   nproc        ;
   int   *owners      ;
   int   nblockUpper  ;
   int   *rowidsUpper ;
   int   *colidsUpper ;
   int   *mapUpper    ;
   int   nblockLower  ;
   int   *rowidsLower ;
   int   *colidsLower ;
   int   *mapLower    ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98mar19, cca
   -----------------------
*/
SolveMap *
SolveMap_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98mar19, cca
   -----------------------
*/
void
SolveMap_setDefaultFields (
   SolveMap   *solvemap
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98mar19, cca
   --------------------------------------------------
*/
void
SolveMap_clearData (
   SolveMap   *solvemap
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98mar19, cca
   ------------------------------------------
*/
SolveMap *
SolveMap_free (
   SolveMap   *solvemap
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------
   purpose -- to return the symmetry flag of the object
 
   created -- 98mar19, cca
   ----------------------------------------------------
*/
int
SolveMap_symmetryflag (
   SolveMap   *solvemap
) ;
/*
   -----------------------------------------
   purpose -- to return the number of fronts
 
   created -- 98mar19, cca
   -----------------------------------------
*/
int
SolveMap_nfront (
   SolveMap   *solvemap
) ;
/*
   ---------------------------------------------
   purpose -- to return the number of processors
 
   created -- 98mar19, cca
   ---------------------------------------------
*/
int
SolveMap_nproc (
   SolveMap   *solvemap
) ;
/*
   ----------------------------------------------------------------
  purpose -- to return the number of blocks in the upper adjacency
 
   created -- 98mar19, cca
   ----------------------------------------------------------------
*/
int
SolveMap_nblockUpper (
   SolveMap   *solvemap
) ;
/*
   ----------------------------------------------------------------
  purpose -- to return the number of blocks in the lower adjacency
 
   created -- 98mar19, cca
   ----------------------------------------------------------------
*/
int
SolveMap_nblockLower (
   SolveMap   *solvemap
) ;
/*
   ---------------------------------------------------
   purpose -- to return a pointer to the owners vector
 
   created -- 98mar19, cca
   ----------------------------------------------------
*/
int *
SolveMap_owners (
   SolveMap   *solvemap
) ;
/*
   ----------------------------------------------------
   purpose -- to return a pointer to the row ids vector
              for the upper adjacency structure
 
   created -- 98mar19, cca
   -----------------------------------------------------
*/
int *
SolveMap_rowidsUpper (
   SolveMap   *solvemap
) ;
/*
   -------------------------------------------------------
   purpose -- to return a pointer to the column ids vector
              for the upper adjacency structure
 
   created -- 98mar19, cca
   --------------------------------------------------------
*/
int *
SolveMap_colidsUpper (
   SolveMap   *solvemap
) ;
/*
   ------------------------------------------------
   purpose -- to return a pointer to the map vector
              for the upper adjacency structure
 
   created -- 98mar19, cca
   -------------------------------------------------
*/
int *
SolveMap_mapUpper (
   SolveMap   *solvemap
) ;
/*
   ----------------------------------------------------
   purpose -- to return a pointer to the row ids vector
              for the upper adjacency structure
 
   created -- 98mar19, cca
   -----------------------------------------------------
*/
int *
SolveMap_rowidsLower (
   SolveMap   *solvemap
) ;
/*
   -------------------------------------------------------
   purpose -- to return a pointer to the column ids vector
              for the upper adjacency structure
 
   created -- 98mar19, cca
   --------------------------------------------------------
*/
int *
SolveMap_colidsLower (
   SolveMap   *solvemap
) ;
/*
   ------------------------------------------------
   purpose -- to return a pointer to the map vector
              for the upper adjacency structure
 
   created -- 98mar19, cca
   -------------------------------------------------
*/
int *
SolveMap_mapLower (
   SolveMap   *solvemap
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- set the scalars and allocate the vectors for the object
 
   created -- 98mar19, cca
   ------------------------------------------------------------------ */
void
SolveMap_init (
   SolveMap   *solvemap,
   int        symmetryflag,
   int        nfront,
   int        nproc,
   int        nblockUpper,
   int        nblockLower
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in maps.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------
   purpose -- map the off diagonal blocks 
      to processes in a random fashion
 
   created -- 98mar19, cca
   --------------------------------------
*/
void
SolveMap_randomMap (
   SolveMap   *solvemap,
   int        symmetryflag,
   IVL        *upperBlockIVL,
   IVL        *lowerBlockIVL,
   int        nproc,
   IV         *ownersIV,
   int        seed,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------------
   purpose -- map the off diagonal blocks to
      processes in a domain decomposition fashion
 
   created -- 98mar28, cca
   ----------------------------------------------
*/
void
SolveMap_ddMap (
   SolveMap   *solvemap,
   int        symmetryflag,
   IVL        *upperBlockIVL,
   IVL        *lowerBlockIVL,
   int        nproc,
   IV         *ownersIV,
   Tree       *tree,
   int        seed,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in setup.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- to set up the linked lists for the forward solve
              local to process myid
 
   created -- 98mar19, cca
   -----------------------------------------------------------
*/
IP **
SolveMap_forwardSetup (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ------------------------------------------------------------
   purpose -- to set up the linked lists for the backward solve
              local to process myid
 
   created -- 98mar19, cca
   ------------------------------------------------------------
*/
IP **
SolveMap_backwardSetup (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------
   purpose -- return the owner of block (rowid, colid).
 
   created -- 98mar19, cca
   ----------------------------------------------------
*/
int
SolveMap_owner (
   SolveMap   *solvemap,
   int        rowid,
   int        colid
) ;
/*
   ---------------------------------------------------------------
   purpose -- return an IVL object whose list K contains all
      processes who do not own U(K,K) but own a U(J,K) for some J.
 
   if myid == -1 then
      the entire IVL object is created and returned
   else
      only the portion of the IVL object pertinent
      to myid is created and returned
   endif
 
   created -- 98may24, cca
   ---------------------------------------------------------------
*/
IVL *
SolveMap_upperSolveIVL (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   -----------------------------------------------------
   purpose -- return an IV object whose entry J contains
     the number of all processes who do not own U(J,J)
     but own a U(J,K) for some J.
 
   if myid == -1 then
      all entries in the vector are filled
   else
      all those entries pertinent to myid are filled
   endif
 
   created -- 98mar19, cca
   -----------------------------------------------------
*/
IV *
SolveMap_upperAggregateIV (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   --------------------------------------------------------
   purpose -- return an IV object whose J'th entry contains
      the number of processes who do not own L(J,J) but own
      a L(J,I) for some I < J.
 
   if myid == -1 then
      all entries in the vector are filled
   else
      all those entries pertinent to myid are filled
   endif
 
   created -- 98mar20, cca
   --------------------------------------------------------
*/
IV *
SolveMap_lowerAggregateIV (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   -------------------------------------------------------
   purpose -- return an IVL object whose list J contains
      the processes who do not own L(J,J) but own a L(K,J)
      for some K > J.
 
   if myid == -1 then
      the entire IVL object is created and returned
   else
      only the portion of the IVL object pertinent
      to myid is created and returned
   endif
 
   created -- 98may24, cca
   -------------------------------------------------------
*/
IVL *
SolveMap_lowerSolveIVL (
   SolveMap   *solvemap,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   purpose -- to read in an SolveMap object from a file
 
   input --
 
      fn -- filename, must be *.solvemapb or *.solvemapf
 
   return value -- 1 if success, 0 if failure
 
   created -- 98apr09, cca
   -----------------------------------------------------
*/
int
SolveMap_readFromFile ( 
   SolveMap   *solvemap, 
   char       *fn 
) ;
/*
   -----------------------------------------------------------
   purpose -- to read an SolveMap object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 98apr09, cca
   -----------------------------------------------------------
*/
int
SolveMap_readFromFormattedFile (
   SolveMap   *solvemap,
   FILE       *fp
) ;
/*
   ---------------------------------------------------
   purpose -- to read an SolveMap object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 98apr09, cca
   ---------------------------------------------------
*/
int
SolveMap_readFromBinaryFile (
   SolveMap    *solvemap,
   FILE   *fp
) ;
/*
   ------------------------------------------------
   purpose -- to write an SolveMap object to a file
 
   input --
 
      fn -- filename
        *.solvemapb -- binary
        *.solvemapf -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98apr09, cca
   ------------------------------------------------
*/
int
SolveMap_writeToFile (
   SolveMap    *solvemap,
   char   *fn
) ;
/*
   -----------------------------------------------------
   purpose -- to write an SolveMap object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98apr09, cca
   -----------------------------------------------------
*/
int
SolveMap_writeToFormattedFile (
   SolveMap    *solvemap,
   FILE   *fp
) ;
/*
   --------------------------------------------------
   purpose -- to write an SolveMap object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98apr09, cca
   --------------------------------------------------
*/
int
SolveMap_writeToBinaryFile (
   SolveMap    *solvemap,
   FILE   *fp
) ;
/*
   ------------------------------------------------------
   purpose -- to write an SolveMap object for a human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98apr09, cca
   ------------------------------------------------------
*/
int
SolveMap_writeForHumanEye (
   SolveMap   *solvemap,
   FILE       *fp
) ;
/*
   --------------------------------------------------------------
   purpose -- to write out the statistics for the SolveMap object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98apr09, cca
   --------------------------------------------------------------
*/
int
SolveMap_writeStats (
   SolveMap   *solvemap,
   FILE       *fp
) ;
/*--------------------------------------------------------------------*/
