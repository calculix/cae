#include "../MT.h"
#include "../InpMtx.h"
#include "../ETree.h"
#include "../FrontMtx.h"
#include "../SymbFac.h"
#include "../misc.h"
#include "../timings.h"
/*
   ----------------------------------------------------------------
   this object is the bridge between the lanczos and spooles codes.

   NOTE: multithreaded version

   graph statistics
      neqns          -- number of vertices in the uncompressed graph
      nedges         -- number of edges in the uncompressed graph
      Neqns          -- number of vertices in the compressed graph
      Nedges         -- number of edges in the compressed graph

   ordering parameters
      maxdomainsize -- maximum domain size for recursive dissection
      maxnzeros     -- maximum number of zeros in a front
      maxsize       -- maximum size of a front
      seed          -- random number seed 
      compressCutoff -- cutoff for graph compression, 
         if Neqns <= compressCutoff * neqns then the compressed
         graph is created, ordered and used to create the symbolic
         factorization

   matrix parameters
      neqns        -- number of equations
      type         -- SPOOLES_REAL or SPOOLES_COMPLEX 
      symmetryflag -- SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN 
         or SPOOLES_NONSYMMETRIC

   factorization parameters
      sparsityflag -- FRONTMTX_DENSE_FRONTS or FRONTMTX_SPARSE_FRONTS
      pivotingflag -- SPOOLES_PIVOTING or SPOOLES_NO_PIVOTING
      tau          -- upper bound on magnitude of factor entries
         when pivoting is used
      droptol      -- lower bound on magnitude of factor entries
         when sparse fronts are stored
      lookahead    -- lookahead parameter for factorization
                      recommended value = nthread

   message info
      msglvl -- message level for SPOOLES programs
         set msglvl =  0 for no output
         set msglvl =  1 for scalar and timing output
         set msglvl >= 2 for lots of output
      msgFile -- message file for debug and diagnostic output

   internal objects, free'd in cleanup
      frontETree -- object that contains the front tree information
      symbfacIVL -- object that contains the symbolic factorization
      mtxmanager -- SubMtx manager object that handles storage 
                    for the submatrices of the factors.
      frontmtx   -- object that contains the factorization
      oldToNewIV -- object that contains the old-to-new permutation
      newToOldIV -- object that contains the new-to-old permutation

   multithreaded information
      nthread   -- # of threads
      ownersIV  -- map from fronts to owning threads
      solvemap  -- map from submatrices to owning threads
      cumopsDV  -- vector that holds operations for each thread

   statistics
      stats[ 0] -- # of pivots
      stats[ 1] -- # of pivot tests
      stats[ 2] -- # of delayed rows and columns
      stats[ 3] -- # of entries in D
      stats[ 4] -- # of entries in L
      stats[ 5] -- # of entries in U

   timings
      cpus[ 0] -- time to construct Graph
      cpus[ 1] -- time to compress Graph
      cpus[ 2] -- time to order Graph
      cpus[ 3] -- time for symbolic factorization
      cpus[ 4] -- total setup time
      cpus[ 5] -- parallel factor setup time
      cpus[ 6] -- time to permute original matrix
      cpus[ 7] -- time to initialize the front matrix
      cpus[ 8] -- time to factor the front matrix
      cpus[ 9] -- time to post-process the front matrix
      cpus[10] -- total factor time
      cpus[11] -- parallel solve setup time
      cpus[12] -- permute rhs
      cpus[13] -- solve time
      cpus[14] -- permute solution
      cpus[15] -- total solve time

   created -- 98sep17, cca
   ----------------------------------------------------------------
*/
typedef struct _BridgeMT BridgeMT ;
struct _BridgeMT {
/*
   ----------------
   graph parameters
   ----------------
*/
  int           neqns          ;
  int           nedges         ;
  int           Neqns          ;
  int           Nedges         ;
/*
   -------------------
   ordering parameters
   -------------------
*/
  int           maxdomainsize  ;
  int           maxnzeros      ;
  int           maxsize        ;
  int           seed           ;
  double        compressCutoff ;
/*
   -----------------
   matrix parameters
   -----------------
*/
  int           type           ;
  int           symmetryflag   ;
/*
   ------------------------
   factorization parameters
   ------------------------
*/
  int              sparsityflag   ;
  int              pivotingflag   ;
  double           tau            ;
  double           droptol        ;
  int              lookahead      ;
  PatchAndGoInfo   *patchinfo    ;
/*
   -------------------
   pointers to objects
   -------------------
*/
  ETree         *frontETree    ;
  IVL           *symbfacIVL    ;
  SubMtxManager *mtxmanager    ;
  FrontMtx      *frontmtx      ;
  IV            *oldToNewIV    ;
  IV            *newToOldIV    ;
/*
   ------------------------
   multithreaded parameters
   ------------------------
*/
  int           nthread        ;
  IV            *ownersIV      ;
  SolveMap      *solvemap      ;
  DV            *cumopsDV      ;
/*
   ---------------------------------
   message info, statistics and cpus
   ---------------------------------
*/
  int           msglvl         ;
  FILE          *msgFile       ;
  int           stats[6]       ;
  double        cpus[16]       ;
} ;

/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   constructor method
 
   created -- 98sep18, cca
   -----------------------
*/
BridgeMT *
BridgeMT_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
 
   created -- 98sep18, cca
   -----------------------
*/
int
BridgeMT_setDefaultFields ( 
   BridgeMT   *bridge
) ;
/*
   -----------------------
   clear the data fields
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
 
   created -- 98sep18, cca
   -----------------------
*/
int
BridgeMT_clearData (
   BridgeMT   *bridge
) ;
/*
   -----------------------
   destructor
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
 
   created -- 98sep18, cca
   -----------------------
*/
int
BridgeMT_free (
   BridgeMT   *bridge
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in setparams.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------
   purpose -- to set the matrix parameters
 
   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- neqns <= 0
    -3 -- type is invalid
    -4 -- symmetryflag is invalid
    -5 -- matrix is hermitian but type is real
 
   created -- 98sep25, cca
   -------------------------------------------
*/
int
BridgeMT_setMatrixParams (
   BridgeMT   *bridge,
   int        neqns,
   int        type,
   int        symmetryflag
) ;
/*
   -------------------------------------------
   purpose -- to set the ordering parameters
 
   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- maxdomainsize <= 0 
    -3 -- maxsize <= 0 
    -4 -- compressCutoff > 1.0
 
   created -- 98sep25, cca
   -------------------------------------------
*/
int
BridgeMT_setOrderingParams (
   BridgeMT   *bridge,
   int         maxdomainsize, 
   int         maxnzeros, 
   int         maxsize,
   int         seed,
   int         compressCutoff
) ;
/*
   -------------------------------------
   purpose -- to set the message info
 
   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- msglvl > 0 and msgFile is NULL
 
   created -- 98sep18, cca
   -------------------------------------
*/
int
BridgeMT_setMessageInfo (
   BridgeMT   *bridge,
   int        msglvl,
   FILE       *msgFile 
) ;
/*
   ----------------------------------------------
   purpose -- to set the factorization parameters
 
   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- sparsityflag is invalid
    -3 -- pivotingflag is invalid
    -4 -- tau < 2.0
    -5 -- droptol < 0.0
    -6 -- lookahead < 0
 
   created -- 98sep25, cca
   ----------------------------------------------
*/
int
BridgeMT_setFactorParams (
   BridgeMT         *bridge,
   int              sparsityflag,
   int              pivotingflag,
   double           tau,
   double           droptol,
   int              lookahead,
   PatchAndGoInfo   *patchinfo
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------
   purpose -- load *pobj with the address of the
              old-to-new permutation IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep18, cca
   ---------------------------------------------
*/
int
BridgeMT_oldToNewIV (
   BridgeMT   *bridge,
   IV         **pobj
) ;
/*
   ---------------------------------------------
   purpose -- load *pobj with the address of the
              new-to-old permutation IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep18, cca
   ---------------------------------------------
*/
int
BridgeMT_newToOldIV (
   BridgeMT   *bridge,
   IV         **pobj
) ;
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the front ETree object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep18, cca
   --------------------------------------
*/
int
BridgeMT_frontETree (
   BridgeMT   *bridge,
   ETree      **pobj
) ;
/*
   ---------------------------------------------
   purpose -- load *pobj with the address of the
              symbolic factorization IVL object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep18, cca
   ---------------------------------------------
*/
int
BridgeMT_symbfacIVL (
   BridgeMT   *bridge,
   IVL        **pobj
) ;
/*
   -----------------------------------------
   purpose -- load *pobj with the address of
              the submatrix manager object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep18, cca
   -----------------------------------------
*/
int
BridgeMT_mtxmanager (
   BridgeMT        *bridge,
   SubMtxManager   **pobj
) ;
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the front matrix object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep18, cca
   --------------------------------------
*/
int
BridgeMT_frontmtx (
   BridgeMT   *bridge,
   FrontMtx   **pobj
) ;
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the owners IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep24, cca
   --------------------------------------
*/
int
BridgeMT_ownersIV (
   BridgeMT   *bridge,
   IV         **pobj
) ;
/*
   -----------------------------------------
   purpose -- load *pobj with the address of
              the solve map SolveMap object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep24, cca
   -----------------------------------------
*/
int
BridgeMT_solvemap (
   BridgeMT   *bridge,
   SolveMap   **pobj
) ;
/*
   ----------------------------------------------------
   purpose -- load *pnthread with the number of threads
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pnthread is NULL
 
   created -- 98sep24, cca
   ----------------------------------------------------
*/
int
BridgeMT_nthread (
   BridgeMT   *bridge,
   int        *pnthread
) ;
/*
   --------------------------------------------------------
   purpose -- load *plookahead with the lookahead parameter
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- plookahead is NULL
 
   created -- 98sep24, cca
   --------------------------------------------------------
*/
int
BridgeMT_lookahead (
   BridgeMT   *bridge,
   int        *plookahead
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in info.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   purpose --  generate and return some statistics 
               about the factor and solve
 
   type -- type of entries
     SPOOLES_REAL or SPOOLES_COMPLEX
   symmetryflag -- symmetry type
     SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
 
   on return ---
      *pnfront     -- # of fronts
      *pnfactorind -- # of factor indices
      *pnfactorent -- # of factor entries
      *pnsolveops  -- # of solve operations 
      *pnfactorops -- # of factor operations 
 
   return values --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- type is bad, must be SPOOLES_REAL or SPOOLES_COMPLEX
     -3 -- symmetryflag is bad, must be SPOOLES_SYMMETRIC,
           SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
     -4 -- type and symmetryflag mismatch
     -5 -- front tree is not present
     -6 -- pnfront is NULL
     -7 -- pnfactorind is NULL
     -8 -- pnfactorent is NULL
     -9 -- pnsolveops is NULL
    -10 -- pnfactorops is NULL
 
   created -- 98oct01, cca
   --------------------------------------------------------------
*/
int
BridgeMT_factorStats (
   BridgeMT   *bridge,
   int        type,
   int        symmetryflag,
   int        *pnfront,
   int        *pnfactorind,
   int        *pnfactorent,
   int        *pnsolveops,
   double     *pnfactorops
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in setup.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   purpose --
 
   given an InpMtx object that contains the structure of A, initialize
     the bridge data structure for the serial factor's and solve's.

   note: all parameters are pointers to be compatible with
         fortran's call by reference.
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- mtxA is NULL
 
   created -- 98sep17, cca
   -------------------------------------------------------------------
*/
int
BridgeMT_setup (
   BridgeMT   *bridge,
   InpMtx     *mtxA
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in factorSetup.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------
   purpose -- to construct the map from fronts to threads,
      and compute operations for each thread.

   nthread -- number of threads
   maptype -- type of map for parallel factorization
     maptype = 1 --> wrap map
     maptype = 2 --> balanced map
     maptype = 3 --> subtree-subset map
     maptype = 4 --> domain decomposition map
   cutoff -- used when maptype = 4 as upper bound on
     relative domain size
  default is maptype = 4 and cutoff = 1/(2*nthread)

   return value --
      1 -- success
     -1 -- bridge is NULL
     -2 -- nthread is invalid, must be > 0
     -3 -- front tree is NULL

   created -- 98sep24, cca
   -------------------------------------------------------
*/
int
BridgeMT_factorSetup (
   BridgeMT   *bridge,
   int        nthread,
   int        maptype,
   double     cutoff
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in factor.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   purpose -- to permute (if necessary) the original matrix,
      and to initialize, factor and postprocess the factor matrix
 
   return value ---
      1 -- normal return, factorization complete
      0 -- factorization did not complete, see error flag
     -1 -- bridge is NULL
     -2 -- mtxA is NULL
     -3 -- perror is NULL 
 
   created -- 98sep18, cca
   --------------------------------------------------------------
*/
int
BridgeMT_factor (
   BridgeMT   *bridge,
   InpMtx     *mtxA,
   int        permuteflag,
   int        *perror
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solveSetup.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- to setup for the parallel solve
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- frontmtx is NULL
     -3 -- frontmtx has not yet been postprocessed
 
   created -- 98sep24, cca
   -----------------------------------------------
*/
int
BridgeMT_solveSetup (
   BridgeMT   *bridge
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solve.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------
   purpose -- to solve the linear system
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- X is NULL
     -3 -- Y is NULL
     -4 -- frontmtx is NULL
     -5 -- mtxmanager is NULL
     -6 -- oldToNewIV not available
     -7 -- newToOldIV not available
 
   created -- 98sep18, cca
   -------------------------------------
*/
int
BridgeMT_solve (
   BridgeMT   *bridge,
   int        permuteflag,
   DenseMtx   *X,
   DenseMtx   *Y
) ;
/*--------------------------------------------------------------------*/
