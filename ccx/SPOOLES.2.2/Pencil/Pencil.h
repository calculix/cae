/*  Pencil.h  */

#include "../InpMtx.h"
#include "../DenseMtx.h"
#include "../Drand.h"
#include "../SPOOLES.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   this object stores a matrix pencil
      A + sigma * B

   created -- 98may02, cca
   ------------------------------------------------------------------
*/
typedef struct _Pencil   Pencil ;
struct _Pencil {
   int      type     ;
   int      symflag  ;
   InpMtx   *inpmtxA ;
   InpMtx   *inpmtxB ;
   double   sigma[2] ;
} ;

#define PENCIL_IS_REAL(pencil)    ((pencil)->type == SPOOLES_REAL)
#define PENCIL_IS_COMPLEX(pencil) ((pencil)->type == SPOOLES_COMPLEX)

#define PENCIL_IS_SYMMETRIC(pencil) \
        ((pencil)->symflag == SPOOLES_SYMMETRIC)
#define PENCIL_IS_HERMITIAN(pencil) \
        ((pencil)->symflag == SPOOLES_HERMITIAN)
#define PENCIL_IS_NONSYMMETRIC(pencil) \
        ((pencil)->symflag == SPOOLES_NONSYMMETRIC)
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98may02, cca
   -----------------------
*/
Pencil *
Pencil_new (
   void
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may02, cca
   -----------------------
*/
void
Pencil_setDefaultFields (
   Pencil   *pencil
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may02, cca
   --------------------------------------------------
*/
void
Pencil_clearData (
   Pencil   *pencil
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may02, cca
   ------------------------------------------
*/
Pencil *
Pencil_free (
   Pencil   *pencil
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   initialize the object
 
   created -- 98may02, cca
   -----------------------
*/
void
Pencil_init (
  Pencil   *pencil,
  int      type,
  int      symflag,
  InpMtx   *inpmtxA,
  double   sigma[],
  InpMtx   *inpmtxB
) ;
/*
   --------------------------
   change the coordinate type
 
   created -- 98may02, cca
   --------------------------
*/
void
Pencil_changeCoordType (
   Pencil   *pencil,
   int       newType
) ;
/*
   -----------------------
   change the storage mode
 
   created -- 98may02, cca
   -----------------------
*/
void
Pencil_changeStorageMode (
   Pencil   *pencil,
   int       newMode
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in permute.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------
   permute the matrix pencil
 
   created -- 98may02, cca
   -------------------------
*/
void
Pencil_permute (
   Pencil   *pencil,
   IV        *rowOldToNewIV,
   IV        *colOldToNewIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in mmm.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------
   compute Y := Y + (A + sigma*B)*X
 
   created -- 98may02, cca
   --------------------------------
*/
void
Pencil_mmm (
   Pencil     *pencil,
   DenseMtx   *Y,
   DenseMtx   *X
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in setup.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------------
  initialize the matrix pencil A + sigma*B
 
   myid -- id of process, used in MPI implementation
      if myid = 0 then 
         the pencil is loaded with the matrices read from the files
     else
         the pencil is loaded with the empty matrices 
      endif
   symflag -- symmetry flag, 
      PENCIL_SYMMETRIC    -- symmetric
      PENCIL_HERMITIAN    -- hermitian
      PENCIL_NONSYMMETRIC -- nonsymmetric
      if symmetric or hermitian, drop entries in lower triangle 
   inpmtxAfile  -- filename for A
   sigma        -- scaling factor
   inpmtxBfile  -- filename for B
   randomflag   -- random flag, 
      if 1 then fill with random numbers
   drand        -- random number generator
   msglvl       -- message level
   msgFile      -- message file
 
   return value -- pointer to a Pencil object
 
   created -- 98may02, cca
   ----------------------------------------------------------------
*/
Pencil *
Pencil_setup (
   int        myid,
   int        symflag,
   char       *inpmtxAfile,
   double     sigma[],
   char       *inpmtxBfile,
   int        randomflag,
   Drand      *drand,
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
   --------------------------------------
   sort and compress the pencil's entries
 
   created -- 98may02, cca
   --------------------------------------
*/
void
Pencil_sortAndCompress (
   Pencil   *pencil 
) ;
/*
   ------------------------------
   convert the storage to vectors
 
   created -- 98may02, cca
   ------------------------------
*/
void
Pencil_convertToVectors (
   Pencil   *pencil 
) ;
/*
   ----------------------------------------------
   map entries to the lower triangle,
   used after a permutation of a symmetric matrix
 
   created -- 98may02, cca
   ----------------------------------------------
*/
void
Pencil_mapToLowerTriangle (
   Pencil   *pencil 
) ;
/*
   ----------------------------------------------
   map entries to the upper triangle,
   used after a permutation of a symmetric matrix
 
   created -- 98may02, cca
   ----------------------------------------------
*/
void
Pencil_mapToUpperTriangle (
   Pencil   *pencil
) ;
/*
   -------------------------------------------------------------
   purpose -- to return the full, symmetric adjacency IVL object
              for the graph of (A + B) + sigma * (A + B)^T
 
   created -- 98may02, cca
   -------------------------------------------------------------
*/
IVL *
Pencil_fullAdjacency (
   Pencil  *pencil
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------
   purpose -- to read in a Pencil object from a file
 
   input --
 
      fn -- filename, must be *.inpmtxb or *.inpmtxf
 
   return value -- 1 if success, 0 if failure
 
   created -- 98may02, cca
   --------------------------------------------------
*/
int
Pencil_readFromFiles (
   Pencil   *pencil,
   char      *inpmtxAfileName,
   char      *inpmtxBfileName
) ;
/*
   ----------------------------------------------------
   purpose -- to write a Pencil object for a human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98may02, cca
   ----------------------------------------------------
*/
int
Pencil_writeForHumanEye (
   Pencil   *pencil,
   FILE      *fp
) ;
/*
   -------------------------------------------------------------
   purpose -- to write out the statistics for the Pencil object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98may02, cca
   -------------------------------------------------------------
*/
int
Pencil_writeStats ( 
   Pencil   *pencil, 
   FILE      *fp 
) ;
/*--------------------------------------------------------------------*/
