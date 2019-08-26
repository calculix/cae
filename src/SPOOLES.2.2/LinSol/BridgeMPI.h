/*  BridgeMPI.h  */

#include "../MPI.h"
#include "../InpMtx.h"
#include "../ETree.h"
#include "../FrontMtx.h"
#include "../SymbFac.h"
#include "../misc.h"
#include "../timings.h"
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   this object is the bridge between the nastran and spooles codes.

   NOTE: MPI version

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
      patchinfo    -- pointer to PatchAndGoInfo object, can be NULL
      lookahead    -- lookahead parameter for factorization
                      recommended value = nthread

   message info
      msglvl -- message level for SPOOLES programs
         set msglvl =  0 for no output
         set msglvl =  1 for scalar and timing output
         set msglvl >= 2 for lots of output
      msgFile -- message file for debug and diagnostic output

   internal objects, free'd in cleanup
      frontETree     -- object that contains the front tree information
      symbfacIVL     -- object that contains the symbolic factorization
      mtxmanager     -- SubMtx manager object that handles storage 
                        for the submatrices of the factors.
      frontmtx       -- object that contains the factorization
      oldToNewIV     -- object that contains the old-to-new permutation
      newToOldIV     -- object that contains the new-to-old permutation
      vtxmapIV       -- object that contains the map from vertices 
                        to processors
      rowmapIV       -- object that contains the map from rows 
                        to processors in the solve
      ownedColumnsIV -- object that contains the owned columns
                        during the solve
      Aloc           -- InpMtx object to hold local owned matrix
      Xloc           -- DenseMtx object to hold local solution
      Yloc           -- DenseMtx object to hold local right hand side

   MPI information
      nproc     -- # of processors
      myid      -- rank of this processors
      comm      -- MPI communicator
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
      cpus[ 4] -- time to broadcast the front tree
      cpus[ 5] -- time to broadcast the symbolic factorization
      cpus[ 6] -- total setup time
      cpus[ 7] -- parallel factor setup time
      cpus[ 8] -- time to permute original matrix
      cpus[ 9] -- time to distribute the original matrix
      cpus[10] -- time to initialize the front matrix
      cpus[11] -- time to factor the front matrix
      cpus[12] -- time to post-process the front matrix
      cpus[13] -- total factor time
      cpus[14] -- parallel solve setup time
      cpus[15] -- time to permute rhs
      cpus[16] -- time to distribute rhs
      cpus[17] -- time to create the solution matrix
      cpus[18] -- time to solve the linear system
      cpus[19] -- time to gather solution on processor zero
      cpus[20] -- time to permute solution
      cpus[21] -- total solve time

   created -- 98sep17, cca
   ----------------------------------------------------------------
*/
typedef struct _BridgeMPI BridgeMPI ;
struct _BridgeMPI {
/*
   ----------------
   graph parameters
   ----------------
*/
   int   neqns  ;
   int   nedges ;
   int   Neqns  ;
   int   Nedges ;
/*
   -------------------
   ordering parameters
   -------------------
*/
   double   compressCutoff ;
   int      maxdomainsize  ;
   int      maxnzeros      ;
   int      maxsize        ;
   int      seed           ;
/*
   -----------------
   matrix parameters
   -----------------
*/
   int   type         ;
   int   symmetryflag ;
/*
   ------------------------
   factorization parameters
   ------------------------
*/
   int              sparsityflag ;
   int              pivotingflag ;
   double           tau          ;
   double           droptol      ;
   int              lookahead    ;
   PatchAndGoInfo   *patchinfo   ;
/*
   -------------------
   pointers to objects
   -------------------
*/
   ETree           *frontETree     ;
   IVL             *symbfacIVL     ;
   SubMtxManager   *mtxmanager     ;
   FrontMtx        *frontmtx       ;
   IV              *oldToNewIV     ;
   IV              *newToOldIV     ;
/*
   -----------------------
   pointers to MPI objects
   -----------------------
*/
   IV              *ownersIV       ;
   SolveMap        *solvemap       ;
   DV              *cumopsDV       ;
   IV              *vtxmapIV       ;
   IV              *rowmapIV       ;
   IV              *ownedColumnsIV ;
   InpMtx          *Aloc           ;
   DenseMtx        *Xloc           ;
   DenseMtx        *Yloc           ;
/*
   --------------
   MPI parameters
   --------------
*/
   int        nproc ;
   int        myid  ;
   MPI_Comm   comm  ;
/*
   ---------------------------------
   message info, statistics and cpus
   ---------------------------------
*/
   int      msglvl   ;
   FILE     *msgFile ;
   int      stats[6] ;
   double   cpus[22] ;
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
 
   created -- 98sep25, cca
   -----------------------
*/
BridgeMPI *
BridgeMPI_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
 
   created -- 98sep25, cca
   -----------------------
*/
int
BridgeMPI_setDefaultFields ( 
   BridgeMPI   *bridge
) ;
/*
   -----------------------
   clear the data fields
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
 
   created -- 98sep25, cca
   -----------------------
*/
int
BridgeMPI_clearData (
   BridgeMPI   *bridge
) ;
/*
   -----------------------
   destructor
 
   return value ---
      1 -- normal return
     -1 -- bridge is NULL
 
   created -- 98sep25, cca
   -----------------------
*/
int
BridgeMPI_free (
   BridgeMPI   *bridge
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
BridgeMPI_setMatrixParams (
   BridgeMPI   *bridge,
   int         neqns,
   int         type,
   int         symmetryflag
) ;
/*
   -------------------------------------------
   purpose -- to set the MPI parameters
 
   return value --
     1 -- normal return
    -1 -- bridge object is NULL
    -2 -- nproc <= 0
    -3 -- myid < 0 or myid >= nproc
 
   created -- 98sep25, cca
   -------------------------------------------
*/
int
BridgeMPI_setMPIparams (
   BridgeMPI   *bridge,
   int         nproc,
   int         myid,
   MPI_Comm    comm
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
BridgeMPI_setOrderingParams (
   BridgeMPI   *bridge,
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
 
   created -- 98sep25, cca
   -------------------------------------
*/
int
BridgeMPI_setMessageInfo (
   BridgeMPI   *bridge,
   int         msglvl,
   FILE        *msgFile
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
BridgeMPI_setFactorParams (
   BridgeMPI        *bridge,
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
BridgeMPI_oldToNewIV (
   BridgeMPI   *bridge,
   IV          **pobj
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
BridgeMPI_newToOldIV (
   BridgeMPI   *bridge,
   IV          **pobj
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
BridgeMPI_frontETree (
   BridgeMPI   *bridge,
   ETree       **pobj
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
BridgeMPI_symbfacIVL (
   BridgeMPI   *bridge,
   IVL         **pobj
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
BridgeMPI_mtxmanager (
   BridgeMPI       *bridge,
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
BridgeMPI_frontmtx (
   BridgeMPI   *bridge,
   FrontMtx    **pobj
) ;
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the owners IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep25, cca
   --------------------------------------
*/
int
BridgeMPI_ownersIV (
   BridgeMPI   *bridge,
   IV          **pobj
) ;
/*
   -----------------------------------------
   purpose -- load *pobj with the address of
              the solve map SolveMap object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98sep25, cca
   -----------------------------------------
*/
int
BridgeMPI_solvemap (
   BridgeMPI   *bridge,
   SolveMap    **pobj
) ;
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the vtxmap IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   --------------------------------------
*/
int
BridgeMPI_vtxmapIV (
   BridgeMPI   *bridge,
   IV          **pobj
) ;
/*
   --------------------------------------
   purpose -- load *pobj with the address
              of the rowmap IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   --------------------------------------
*/
int
BridgeMPI_rowmapIV (
   BridgeMPI   *bridge,
   IV          **pobj
) ;
/*
   ----------------------------------------
   purpose -- load *pobj with the address
              of the ownedColumns IV object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   ----------------------------------------
*/
int
BridgeMPI_ownedColumns (
   BridgeMPI   *bridge,
   IV          **pobj
) ;
/*
   ----------------------------------------
   purpose -- load *pobj with the address
              of the Xloc DenseMtx object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   ----------------------------------------
*/
int
BridgeMPI_Xloc (
   BridgeMPI   *bridge,
   DenseMtx    **pobj
) ;
/*
   ----------------------------------------
   purpose -- load *pobj with the address
              of the Yloc DenseMtx object
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pobj is NULL
 
   created -- 98oct01, cca
   ----------------------------------------
*/
int
BridgeMPI_Yloc (
   BridgeMPI   *bridge,
   DenseMtx    **pobj
) ;
/*
   -----------------------------------------------------
   purpose -- load *pnproc with the number of processors
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pnproc is NULL
 
   created -- 98sep18, cca
   -----------------------------------------------------
*/
int
BridgeMPI_nproc (
   BridgeMPI   *bridge,
   int         *pnproc
) ;
/*
   ----------------------------------------------------
   purpose -- load *pmyid with the id of this processor
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- pmyid is NULL
 
   created -- 98sep18, cca
   ----------------------------------------------------
*/
int
BridgeMPI_myid (
   BridgeMPI   *bridge,
   int         *pmyid
) ;
/*
   ----------------------------------------------------
   purpose -- load *plookahead with the lookahead value 
              for the factorization
 
   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- plookahead is NULL
 
   created -- 98sep18, cca
   ----------------------------------------------------
*/
int
BridgeMPI_lookahead (
   BridgeMPI   *bridge,
   int         *plookahead
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
BridgeMPI_factorStats (
   BridgeMPI   *bridge,
   int         type,
   int         symmetryflag,
   int         *pnfront,
   int         *pnfactorind,
   int         *pnfactorent,
   int         *pnsolveops,
   double      *pnfactorops
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
     -2 -- myid is zero and mtxA is NULL
 
   created -- 98sep25, cca
   -------------------------------------------------------------------
*/
int
BridgeMPI_setup (
   BridgeMPI   *bridge,
   InpMtx      *mtxA
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in factorSetup.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------
   purpose -- to construct the map from fronts to processors,
      and compute operations for each processor.
 
   maptype -- type of map for parallel factorization
      maptype = 1 --> wrap map
      maptype = 2 --> balanced map
      maptype = 3 --> subtree-subset map
      maptype = 4 --> domain decomposition map
   cutoff -- used when maptype = 4 as upper bound on
      relative domain size
 
   return value --
      1 -- success
     -1 -- bridge is NULL
     -2 -- front tree is NULL
 
   created -- 98sep25, cca
   ----------------------------------------------------------
*/
int
BridgeMPI_factorSetup (
   BridgeMPI   *bridge,
   int         maptype,
   double      cutoff
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
      if permuteflag == 1 then
         matrix is permuted into new ordering
      endif
 
   return value ---
      1 -- normal return, factorization complete
      0 -- factorization did not complete, see error flag
     -1 -- bridge is NULL
     -2 -- perror is NULL
 
   created -- 98sep18, cca
   --------------------------------------------------------------
*/
int
BridgeMPI_factor (
   BridgeMPI        *bridge,
   InpMtx           *mtxA,
   int              permuteflag,
   int              *perror
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
BridgeMPI_solveSetup (
   BridgeMPI   *bridge
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solve.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------
   purpose -- to solve the linear system
      MPI version
      if permuteflag is 1 then
         rhs is permuted into new ordering
         solution is permuted into old ordering
 
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
   --------------------------------------------
*/
int
BridgeMPI_solve (
   BridgeMPI  *bridge,
   int        permuteflag,
   DenseMtx   *X,
   DenseMtx   *Y
) ;
/*--------------------------------------------------------------------*/
