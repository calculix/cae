/*  FrontMtx.h  */

#include "../Pencil.h"
#include "../ETree.h"
#include "../IVL.h"
#include "../PatchAndGoInfo.h"
#include "../Chv.h"
#include "../ChvManager.h"
#include "../ChvList.h"
#include "../SubMtx.h"
#include "../SubMtxList.h"
#include "../SubMtxManager.h"
#include "../DenseMtx.h"
#include "../Ideq.h"
#include "../SolveMap.h"
#include "../Lock.h"
#include "../I2Ohash.h"

#include "../SPOOLES.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   the FrontMtx object is used to compute and store a matrix
   factorization in serial, multithreaded and MPI modes.

   there are two data modes:
     1 --> data is stored by fronts, a 1-D data decomposition.
           we use five pointer vectors (p_mtxDJJ[], p_mtxUJJ[],
           p_mtxUJN[], p_mtxLJJ[] and p_mtxLNJ[]) to store pointers
           to the factor submatrices
     2 --> data is stored by submatrices, a 2-D data decomposition.
           we use two hash objects (lowerhash and upperhash), a
           pointer vector p_mtxDJJ[] and two IVL objects (lowerblockIVL
           and upperblockIVL) to store the structure of the matrix.

   nfront -- number of fronts in the matrix
   neqns  -- number of rows and columns in the matrix
   type   -- type of entries
      1 -- real
      2 -- complex
   symmetryflag -- flag to specify symmetry of the matrix
      0 --> symmetric structure, symmetric entries
      1 --> symmetric structure, nonsymmetric entries
      2 --> nonsymmetric structure, nonsymmetric entries
   sparsityflag -- flag to specify dense or sparse fronts
      0 --> dense fronts
      1 --> sparse fronts
   pivotingflag -- flag to specify pivoting enabled
      0 --> pivoting not used
      1 --> pivoting used
   dataMode -- flag to specify storage mode
      1 --> 1-dimensional data decomposition, used for factorization
      2 --> 2-dimensional data decomposition, used for solves
   nentD -- number of entries in the diagonal matrix
   nentL -- number of entries in the lower triangular matrix
   nentU -- number of entries in the upper triangular matrix

   tree       -- pointer to an Tree object that holds the tree of fronts
   frontETree -- pointer to an ETree object that holds the front tree
   symbfacIVL -- pointer to an IVL object that holds 
      the symbolic factorization

   frontsizesIV -- pointer to an IV object that holds the number of
      internal rows and columns in each front
   rowadjIVL  -- pointer to an IVL object that holds the row ids
      of the fronts, used only for nonsymmetric matrices with pivoting
   coladjIVL  -- pointer to an IVL object that holds the column ids
      of the fronts, used only with pivoting

   p_mtxDJJ -- vector of pointers to the diagonal D_{J,J} objects
   p_mtxUJJ -- vector of pointers to the upper U_{J,J} objects
   p_mtxUJN -- vector of pointers to the upper U_{J,N} objects
   p_mtxLJJ -- vector of pointers to the lower L_{J,J} objects
   p_mtxLNJ -- vector of pointers to the lower L_{N,J} objects

   lowerblockIVL -- pointer to an IVL object that holds the sparsity
      structure of the block L matrix, i.e., front-front edges,
      used only for nonsymmetric matrices with pivoting 
   upperblockIVL -- pointer to an IVL object that holds the sparsity
      structure of the block U matrix, i.e., front-front edges
   lowerhash -- pointer to a hash table object that holds the
      submatrices in L, used only for nonsymmetric factorizations
   upperhash -- pointer to a hash table object that holds the
      submatrices in U

   manager -- object to manage instances of Mtx objects
   lock -- mutex object that controls access to allocating
      storage in IVL and DVL objects, can be NULL
   nlocks -- number of times the lock was locked
   patchinfo -- information object needed for special factorizations

   created -- 98feb27, cca
   -------------------------------------------------------------------
*/
typedef struct _FrontMtx   FrontMtx ;
struct _FrontMtx {
   int             nfront         ;
   int             neqns          ;
   int             type           ;
   int             symmetryflag   ;
   int             sparsityflag   ;
   int             pivotingflag   ;
   int             dataMode       ;
   int             nentD          ;
   int             nentL          ;
   int             nentU          ;
   Tree            *tree          ;
   ETree           *frontETree    ;
   IV              *frontsizesIV  ;
   IVL             *symbfacIVL    ;
   IVL             *rowadjIVL     ;
   IVL             *coladjIVL     ;
   IVL             *lowerblockIVL ;
   IVL             *upperblockIVL ;
   SubMtx          **p_mtxDJJ     ;
   SubMtx          **p_mtxUJJ     ;
   SubMtx          **p_mtxUJN     ;
   SubMtx          **p_mtxLJJ     ;
   SubMtx          **p_mtxLNJ     ;
   I2Ohash         *lowerhash     ;
   I2Ohash         *upperhash     ;
   SubMtxManager   *manager       ;
   Lock            *lock          ;
   PatchAndGoInfo  *patchinfo     ;
   int             nlocks         ;
} ;

#define FRONTMTX_IS_REAL(mtx)     ((mtx)->type == SPOOLES_REAL)
#define FRONTMTX_IS_COMPLEX(mtx)  ((mtx)->type == SPOOLES_COMPLEX)

#define FRONTMTX_IS_SYMMETRIC(mtx)    \
   ((mtx)->symmetryflag == SPOOLES_SYMMETRIC)
#define FRONTMTX_IS_HERMITIAN(mtx)    \
   ((mtx)->symmetryflag == SPOOLES_HERMITIAN)
#define FRONTMTX_IS_NONSYMMETRIC(mtx) \
   ((mtx)->symmetryflag == SPOOLES_NONSYMMETRIC)

#define FRONTMTX_DENSE_FRONTS  0
#define FRONTMTX_SPARSE_FRONTS 1
#define FRONTMTX_IS_DENSE_FRONTS(mtx)  \
   ((mtx)->sparsityflag == FRONTMTX_DENSE_FRONTS)
#define FRONTMTX_IS_SPARSE_FRONTS(mtx) \
   ((mtx)->sparsityflag == FRONTMTX_SPARSE_FRONTS)

#define FRONTMTX_IS_PIVOTING(mtx)  \
   ((mtx)->pivotingflag == SPOOLES_PIVOTING)

#define FRONTMTX_1D_MODE      1
#define FRONTMTX_2D_MODE      2
#define FRONTMTX_IS_1D_MODE(mtx)  \
   ((mtx)->dataMode == FRONTMTX_1D_MODE)
#define FRONTMTX_IS_2D_MODE(mtx)  \
   ((mtx)->dataMode == FRONTMTX_2D_MODE)

#define NO_LOCK               0
#define LOCK_IN_PROCESS       1
#define LOCK_IN_ALL_PROCESSES 2
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98may04, cca
   -----------------------
*/
FrontMtx *
FrontMtx_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may04, cca
   -----------------------
*/
void
FrontMtx_setDefaultFields (
   FrontMtx   *frontmtx
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may04, cca
   --------------------------------------------------
*/
void
FrontMtx_clearData (
   FrontMtx   *frontmtx
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may04, cca
   ------------------------------------------
*/
void
FrontMtx_free (
   FrontMtx   *frontmtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------
   purpose -- return the number of fronts
 
   created -- 98may04, cca
   --------------------------------------
*/
int
FrontMtx_nfront (
   FrontMtx   *frontmtx
) ;
/*
   -----------------------------------------
   purpose -- return the number of equations
 
   created -- 98may04, cca
   -----------------------------------------
*/
int
FrontMtx_neqns (
   FrontMtx   *frontmtx
) ;
/*
   ----------------------------------------------------
   purpose -- return a pointer to the front Tree object
 
   created -- 98may04, cca
   ----------------------------------------------------
*/
Tree *
FrontMtx_frontTree (
   FrontMtx   *frontmtx
) ;
/*
   ----------------------------------------------------------------
  simple method to return the dimensions of front J and the number
   of bytes necessary for the Chv object to hold the front.
 
   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_initialFrontDimensions (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnD,
   int         *pnL,
   int         *pnU,
   int         *pnbytes
) ;
/*
   ---------------------------------------------------------
   return the number of internal rows and columns in front J
 
   created -- 98may04, cca
   ---------------------------------------------------------
*/
int
FrontMtx_frontSize (
   FrontMtx   *frontmtx,
   int         J
) ;
/*
   ------------------------------------------------------
   set the number of internal rows and columns in front J
 
   created -- 98may04, cca
   ------------------------------------------------------
*/
void
FrontMtx_setFrontSize (
   FrontMtx   *frontmtx,
   int         J,
   int         size
) ;
/*
   ---------------------------------------------
   fill *pncol with the number of columns and
   *pcolind with a pointer to the column indices
 
   created -- 98may04, cca
   ---------------------------------------------
*/
void
FrontMtx_columnIndices (
   FrontMtx   *frontmtx,
   int         J,
   int         *pncol,
   int         **pcolind
) ;
/*
   -------------------------------------------
   fill *pnrow with the number of rows and
   *prowind with a pointer to the rows indices
 
   created -- 98may04, cca
   -------------------------------------------
*/
void
FrontMtx_rowIndices (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnrow,
   int         **prowind
) ;
/*
   -----------------------------------------------------------
   purpose -- return a pointer to the (J,J) diagonal submatrix
 
   created -- 98may04, cca
   -----------------------------------------------------------
*/
SubMtx *
FrontMtx_diagMtx (
   FrontMtx   *frontmtx,
   int         J
) ;
/*
   --------------------------------------------------------
   purpose -- return a pointer to the (J,K) upper submatrix
 
   created -- 98may04, cca
   --------------------------------------------------------
*/
SubMtx *
FrontMtx_upperMtx (
   FrontMtx   *frontmtx,
   int         J,
   int         K
) ;
/*
   --------------------------------------------------------
   purpose -- return a pointer to the (K,J) lower submatrix
 
   created -- 98may04, cca
   --------------------------------------------------------
*/
SubMtx *
FrontMtx_lowerMtx (
   FrontMtx   *frontmtx,
   int         K,
   int         J
) ;
/*
   --------------------------------------------------
   purpose -- fill *pnadj with the number of fronts K
              such that L_{K,J} != 0 and *padj with a
              pointer to a list of those fronts
 
   created -- 98may04, cca
   --------------------------------------------------
*/
void
FrontMtx_lowerAdjFronts (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnadj,
   int         **padj
) ;
/*
   --------------------------------------------------
   purpose -- fill *pnadj with the number of fronts K
              such that U_{J,K} != 0 and *padj with a
              pointer to a list of those fronts
 
   created -- 98may04, cca
   --------------------------------------------------
*/
void
FrontMtx_upperAdjFronts (
   FrontMtx   *frontmtx,
   int         J,
   int         *pnadj,
   int         **padj
) ;
/*
   ------------------------------------------------------
   purpose -- return the number of nonzero L_{K,J} blocks
 
   created -- 98may04, cca
   ------------------------------------------------------
*/
int
FrontMtx_nLowerBlocks (
   FrontMtx   *frontmtx
) ;
/*
   ------------------------------------------------------
   purpose -- return the number of nonzero U_{K,J} blocks
 
   created -- 98may04, cca
   ------------------------------------------------------
*/
int
FrontMtx_nUpperBlocks (
   FrontMtx   *frontmtx
) ;
/*
   ---------------------------------------------------------
   purpose -- return a pointer to the upper block IVL object
 
   created -- 98jun13, cca
   ---------------------------------------------------------
*/
IVL *
FrontMtx_upperBlockIVL (
   FrontMtx   *frontmtx
) ;
/*
   ---------------------------------------------------------
   purpose -- return a pointer to the lower block IVL object
 
   created -- 98jun13, cca
   ---------------------------------------------------------
*/
IVL *
FrontMtx_lowerBlockIVL (
   FrontMtx   *frontmtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- basic initializer
 
   frontETree -- ETree object that stores the front tree
   symbfacIVL -- IVL object that stores the symbolic factorization
   manager    -- SubMtxManager object to manage SubMtx objects
   type       -- type of entries
      SPOOLES_REAL --> real
      SPOOLES_COMPLEX --> complex
   symmetryflag -- symmetry flag,
      SPOOLES_SYMMETRIC --> symmetric structure and entries
      SPOOLES_HERMITIAN --> hermitian (complex only)
      SPOOLES_NONSYMMETRIC --> nonsymmetric entries
   sparsityflag -- flag to specify dense or sparse fronts
      FRONTMTX_DENSE_FRONTS --> dense fronts
      FRONTMTX_SPARSE_FRONTS --> sparse fronts
   pivotingflag -- flag to specify pivoting enabled
      SPOOLES_NO_PIVOTING --> pivoting not used
      SPOOLES_PIVOTING    --> pivoting used
 
   in a multithreaded environment, we need to protect the critical
   section where data is allocated. we use a lockflag to do this.
   in a serial or distributed environment, use lockflag = 0.
 
   lockflag -- flag to specify lock status
      NO_LOCK --> mutex lock is not allocated or initialized
      LOCK_IN_PROCESS --> mutex lock is allocated and
         it can synchronize only threads in this process.
      LOCK_OVER_ALL_PROCESSES --> mutex lock is allocated and
          it can synchronize only threads in this and other processes.
 
   in a distributed environment we need to specify which process
   owns each front. when we can preallocate data structures
   (when there is no pivoting and dense fronts are stored) we
   need each process to determine what parts of the data it
   can allocate and set up. in a serial or multithreaded 
   environment, use ownersIV = NULL.
 
      ownersIV -- map from fronts to owning processes
      myid     -- id of this process.
 
   submatrices (be they lower, diagonal, block diagonal, upper)
   are stored in SubMtx objects. the management of these objects,
   (allocation and deallocation) is managed by the SubMtxManager
   manager object.
 
      manager -- SubMtxManager object to handle the submatrices
 
   created  -- 98may04, cca
   ------------------------------------------------------------------
*/
void
FrontMtx_init (
   FrontMtx        *frontmtx,
   ETree           *frontETree,
   IVL             *symbfacIVL,
   int             type,
   int             symmetryflag,
   int             sparsityflag,
   int             pivotingflag,
   int             lockflag,
   int             myid,
   IV              *ownersIV,
   SubMtxManager   *manager,
   int             msglvl,
   FILE            *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in factor.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   compute an (U^T + I)D(I + U) or (L + I)D(I + L) factorization of A.
   this is a wrapper method around FrontMtx_factorPencil().
 
   input --
 
      frontmtx -- pointer to the FrontMtx object that will hold
                  the factorization
      pencil   -- pointer to the Pencil object that holds A + sigma*B
      tau      -- upper bound on entries in L and U,
                  used only when pivoting enabled
      droptol  -- lower bound on entries in L and U,
                  used only when sparsity enabled
      perror   -- error flag, on return
         *perror >= 0 --> factorization failed at front *perror
         *perror <  0 --> factorization succeeded
      cpus[]   -- timing array
         cpus[0] -- initialize fronts
         cpus[1] -- load original entries
         cpus[2] -- get updates from descendents
         cpus[3] -- assembled postponed data
         cpus[4] -- factor the front
         cpus[5] -- extract postponed data
         cpus[6] -- store factor entries
         cpus[7] -- miscellaneous time
         cpus[8] -- total time
      stats[] -- statistics array
         stats[0] -- # of pivots
         stats[1] -- # of pivot tests
         stats[2] -- # of delayed rows and columns
         stats[3] -- # of entries in D
         stats[4] -- # of entries in L
         stats[5] -- # of entries in U
      msglvl   -- message level
      msgFile  -- message file
 
   created  -- 98mar27, cca
   modified -- 98mar27, cca
      perror added to argument list
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_factorInpMtx (
   FrontMtx     *frontmtx,
   InpMtx       *inpmtx,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) ;
/*
   -------------------------------------------------------------------
   compute an (U^T + I)D(I + U) or (L + I)D(I + L)
   factorization of A + sigma*B.
 
   input --
 
      frontmtx -- pointer to the FrontMtx object that will hold
                  the factorization
      pencil   -- pointer to the Pencil object that holds A + sigma*B
      tau      -- upper bound on entries in L and U,
                  used only when pivoting enabled
      droptol  -- lower bound on entries in L and U,
                  used only when sparsity enabled
      perror   -- error flag, on return
         *perror >= 0 --> factorization failed at front *perror
         *perror <  0 --> factorization succeeded
      cpus[]   -- timing array
         cpus[0] -- initialize fronts
         cpus[1] -- load original entries
         cpus[2] -- get updates from descendents
         cpus[3] -- assembled postponed data
         cpus[4] -- factor the front
         cpus[5] -- extract postponed data
         cpus[6] -- store factor entries
         cpus[7] -- miscellaneous time
         cpus[8] -- total time
      stats[] -- statistics array
         stats[0] -- # of pivots
         stats[1] -- # of pivot tests
         stats[2] -- # of delayed rows and columns
         stats[3] -- # of entries in D
         stats[4] -- # of entries in L
         stats[5] -- # of entries in U
      msglvl   -- message level
      msgFile  -- message file
 
   created  -- 98mar27, cca
   modified -- 98mar27, cca
      perror added to argument list
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_factorPencil (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in factorUtil.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------
   purpose -- initialize a front
 
   created -- 98may04, cca
   -----------------------------
*/
void
FrontMtx_initializeFront (
   FrontMtx   *frontmtx,
   Chv        *frontJ,
   int        J
) ;
/*
   ------------------------------------------------------------------
   purpose -- to visit a front during a factorization.
              note: this method is called by the serial,
              multithreaded and MPI factorization codes.
 
   frontmtx     -- front matrix object
   pencil       -- matrix pencil for A + sigma*B
   J            -- id of front we are working on
   myid         -- id of thread or process 
   owners[]     -- map from fronts to owning threads,
                   in a serial environment, owners = NULL
   fronts[]     -- vector of pointers to fronts
   lookahead    -- parameter controls the partial upward visiting
                   of ancestors before returning
   tau          -- used when pivoting enabled,
                   |L_{j,i}| and |U_{i,j}| <= tau
   droptol      -- used for an approximate factorization
                   stored |L_{j,i}| and |U_{i,j}| > droptol
   status[]     -- status vector for the fronts,
                   'W' -- needs to be woke up
                   'A' -- active front
                   'F' -- front is finished
   heads[]      -- vector of pointers to IP objects that store the
                   front-to-front update lists
   pivotsizesIV -- IV object used during the factorization of a front
                   when pivoting is enabled
   workDV       -- DV object used for working storage
   parent       -- front parent vector
   aggList      -- list object used in parallel environment, used
                   to store aggregate fronts
   postList     -- list object used in pivoting and/or parallel
                   environments, used to stored delayed data and/or
                   to synchronize the factorization
   chvmanager   -- used to manage working storage for Chv objects
   stats[]      -- statistics vector
   cpus[]       -- vector to hold breakdown of cpu times
   msglvl       -- message level
   msgFil       -- message file
 
   created -- 98may04, cca
   ------------------------------------------------------------------
*/
char
FrontMtx_factorVisit (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   int          J,
   int          myid,
   int          owners[],
   Chv          *fronts[],
   int          lookahead,
   double       tau,
   double       droptol,
   char         status[],
   IP           *heads[],
   IV           *pivotsizesIV,
   DV           *workDV,
   int          parent[],
   ChvList      *aggList,
   ChvList      *postList,
   ChvManager   *chvmanager,
   int          stats[],
   double       cpus[],
   int          msglvl,
   FILE         *msgFile
) ;
/*
   --------------------------------------------
   purpose -- set up the front's data structure
 
   created -- 98mar27, cca
   --------------------------------------------
*/
Chv *
FrontMtx_setupFront (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   int          J,
   int          myid,
   int          owners[],
   ChvManager   *chvmanager,
   double       cpus[],
   int          msglvl,
   FILE         *msgFile
) ;
/*
   --------------------------------------------------------------------
   purpose -- to set up the link data structures
      for a parallel factorization for process myid
 
   return value -- pointer to IP *heads[nfront+2], which contains
      the beginning of a list of IP objects that store the remaining
      updates to the fronts.
      note, heads[nfront] is the first IP object in the free list.
      heads[nfront+1] is the base address of the allocated IP objects.
 
   created -- 98mar07, cca
   --------------------------------------------------------------------
*/
IP **
FrontMtx_factorSetup (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   -----------------------------------------------
   create and return the nactiveChild vector.
   nactiveChild[J] contains the number of children
   of J that belong to an active path
 
   created -- 97jul03, cca
   -----------------------------------------------
*/
int *
FrontMtx_nactiveChild (
   FrontMtx   *frontmtx,
   char       *status,
   int        myid
) ;
/*
   ---------------------------------------------------------
   create, initialize and return a Ideq object
   that will be used for a parallel factorization,
   forward solve, or backward solve.
 
   the status[] vector will be set as follows:
      status[J] = activeflag   if J is on an active path
      status[J] = inactiveflag if J is not on an active path
 
   created -- 98mar27, cca
   ---------------------------------------------------------
*/
Ideq *
FrontMtx_setUpDequeue (
   FrontMtx   *frontmtx,
   int        owners[],
   int        myid,
   char       status[],
   IP         *heads[],
   char       activeFlag,
   char       inactiveFlag,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   -----------------------------------------------------------------
   purpose -- load the dequeue with the leaves of the active subtree
              used for a factorization and forward solve
 
   created -- 98mar27, cca
   -----------------------------------------------------------------
*/
void
FrontMtx_loadActiveLeaves (
   FrontMtx   *frontmtx,
   char       status[],
   char       activeFlag,
   Ideq       *dequeue
) ;
/*
   -----------------------------------------------
   create, initialize and return a ChvList object
   to deal with postponed chevrons
 
   created -- 97jul03, cca
   -----------------------------------------------
*/
ChvList *
FrontMtx_postList (
   FrontMtx   *frontmtx,
   IV          *frontOwnersIV,
   int         lockflag
) ;
/*
   -----------------------------------------------
   create, initialize and return a ChvList object
   to deal with aggregate chevrons
 
   created -- 97jul03, cca
   -----------------------------------------------
*/
ChvList *
FrontMtx_aggregateList (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        lockflag
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in loadEntries.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   load entries from sigma*A
 
   chv     -- pointer to the Chv object that holds the front
   pencil  -- pointer to a Pencil that holds the matrix entries 
   msglvl  -- message level
   msgFile -- message file
 
   created  -- 97jul18, cca
   ------------------------------------------------------------
*/
void
FrontMtx_loadEntries (
   Chv      *chv,
   Pencil   *pencil,
   int      msglvl,
   FILE     *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in update.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   accumulate updates to front J, store them in the Chv object
 
   created -- 98may04, cca
   ------------------------------------------------------------
*/
void
FrontMtx_update (
   FrontMtx   *frontmtx,
   Chv        *frontJ,
   IP         *heads[],
   char       status[],
   DV         *tempDV,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in postponed.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- to assemble any postponed data into frontJ
 
   frontJ  -- pointer to Chv objec that contains current front
   chvlist -- pointer to a ChvList object that handles the
              lists of postponed Chv objects
   chvmanager -- pointer to a ChvManager object for the list
                 of free Chv objects
   pndelay -- pointer to address to contain the # of delayed rows
              and columns that were assembled into the front
 
   return value -- pointer to Chv object that contains the new front
 
   created -- 98may04, cca
   ------------------------------------------------------------------
*/
Chv *
FrontMtx_assemblePostponedData (
   FrontMtx     *frontmtx,
   Chv          *frontJ,
   ChvList      *chvlist,
   ChvManager   *chvmanager,
   int          *pndelay
) ;
/*
   ---------------------------------------------------------
   purpose -- extract and store the postponed data
 
   frontJ  -- pointer to present front object
   npost   -- # of postponed rows and columns in frontJ
   K       -- parent of J
   chvlist -- pointer to a ChvList object that handles the
              lists of postponed Chv objects
              a singly linked list to assemble
   chvmanager -- pointer to a ChvManager object for the list
                 of free Chv objects
 
   created -- 98may04, cca
   ---------------------------------------------------------
*/
void
FrontMtx_storePostponedData (
   FrontMtx     *frontmtx,
   Chv          *frontJ,
   int          npost,
   int          K,
   ChvList      *chvlist,
   ChvManager   *chvmanager
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in storeFront.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------------
   purpose -- to store the factor indicies and entries from front J
 
   pivotsizesIV -- used for symmetric or hermitian and pivoting
   droptol      -- used for drop tolerance factorization,
                   an entry is stored if its magnitude > droptol
 
   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_storeFront (
   FrontMtx   *frontmtx,
   Chv        *frontJ,
   IV         *pivotsizesIV,
   double     droptol,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in postProcess.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   purpose -- post-process the factorization
      (1) permute row and column adjacency objects if necessary
      (2) permute lower and upper matrices if necessary
      (3) update the block adjacency objects if necessary
      (4) split the chevron submatrices into submatrices
          and make the submatrix indices local w.r.t their fronts
 
   created -- 98mar05, cca
   --------------------------------------------------------------
*/
void
FrontMtx_postProcess (
   FrontMtx   *frontmtx,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in permute.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- permute the upper adjacency structure so that the
      indices in bnd{J} are in ascending order w.r.t. their ancestors
 
   created -- 98mar05, cca
   ------------------------------------------------------------------
*/
void
FrontMtx_permuteUpperAdj (
   FrontMtx   *frontmtx,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ------------------------------------------------------------------
   purpose -- permute the lower adjacency structure so that the
      indices in bnd{J} are in ascending order w.r.t. their ancestors
 
   created -- 98mar05, cca
   ------------------------------------------------------------------
*/
void
FrontMtx_permuteLowerAdj (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) ;
/*
   ------------------------------------------------------------
   purpose -- if the columns indices of the front matrix are in
      different order than the column indices of U_{J,bnd{J}},
      sort the columns of U_{J,bnd{J}} into ascending order
      w.r.t the column indices of the front matrix.
 
   created -- 98mar05, cca
   ------------------------------------------------------------
*/
void
FrontMtx_permuteUpperMatrices (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) ;
/*
   --------------------------------------------------------
   purpose -- if the row indices of the front matrix are in
      different order than the row indices of L_{bnd{J},J},
      sort the rows of L_{bnd{J},J} into ascending order
      w.r.t the row indices of the front matrix.
 
   created -- 98mar05, cca
   --------------------------------------------------------
*/
void
FrontMtx_permuteLowerMatrices (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in split.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------------
   purpose -- for each U_{J,bnd{J}} matrix, remove from hash table,
              split into their U_{J,K} submatrices and insert
              into the hash table.
 
   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_splitUpperMatrices (
   FrontMtx   *frontmtx,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------------------------------
   purpose -- for each L_{bnd{J},J} matrix, remove from hash table,
              split into their L_{K,J} submatrices and insert
              into the hash table.
 
   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_splitLowerMatrices (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solve.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   purpose -- to solve a linear system 
 
   frontmtx -- FrontMtx object that holds the factor matrices
   solmtx   -- DenseMtx that holds the solution
   rhsmtx   -- DenseMtx that holds the right hand side matrix
      note: right hand side entries are copied from rhsmtx,
            and solution entries are copied into solmtx.
            when the row indices of rhsmtx and solmtx are
            identical, rhsmtx and solmtx can be the same object.
   cpus -- vector to hold cpu breakdown time
      cpus[0] -- set up solves
      cpus[1] -- fetch rhs and store solution
      cpus[2] -- forward solve
      cpus[3] -- diagonal solve
      cpus[4] -- backward solve
      cpus[5] -- total time
   mtxmanager -- object that manages working storage
   msglvl  -- message level
   msgFile -- message file
 
   created -- 98may04, cca
   ------------------------------------------------------------
*/
void
FrontMtx_solve (
   FrontMtx        *frontmtx,
   DenseMtx        *solmtx,
   DenseMtx        *rhsmtx,
   SubMtxManager   *mtxmanager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solveUtil.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------
   load the right hand side for the owned fronts
 
   created -- 98mar19, cca
   ---------------------------------------------
*/
SubMtx **
FrontMtx_loadRightHandSide ( 
   FrontMtx        *frontmtx,
   DenseMtx        *rhsmtx,
   int             owners[],
   int             myid,
   SubMtxManager   *mtxmanager,
   int             msglvl,
   FILE            *msgFile
) ;
/*
   --------------------------------------
   visit front J during the forward solve
 
   created -- 98mar27, cca
   --------------------------------------
*/
void
FrontMtx_forwardVisit (
   FrontMtx        *frontmtx,
   int             J,
   int             nrhs,
   int             *owners,
   int             myid,
   SubMtxManager   *mtxmanager,
   SubMtxList      *aggList,
   SubMtx          *p_mtx[],
   char            frontIsDone[],
   IP              *heads[],
   SubMtx          *p_agg[],
   char            status[],
   int             msglvl,
   FILE            *msgFile
) ;
/*
   --------------------------------------------------
   purpose -- visit front J during the diagonal solve
 
   created -- 98feb20, cca
   --------------------------------------------------
*/
void
FrontMtx_diagonalVisit (
   FrontMtx   *frontmtx,
   int        J,
   int        owners[],
   int        myid,
   SubMtx     *p_mtx[],
   char       frontIsDone[],
   SubMtx     *p_agg[],
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ---------------------------------------
   visit front J during the backward solve
 
   created -- 98mar27, cca
   ---------------------------------------
*/
void
FrontMtx_backwardVisit (
   FrontMtx        *frontmtx,
   int             J,
   int             nrhs,
   int             *owners,
   int             myid,
   SubMtxManager   *mtxmanager,
   SubMtxList      *aggList,
   SubMtx          *p_mtx[],
   char            frontIsDone[],
   IP              *heads[],
   SubMtx          *p_agg[],
   char            status[],
   int             msglvl,
   FILE            *msgFile
) ;
/*
   ---------------------------------------------------
   purpose -- move the solution from the individual
     SubMtx objects into the global solution SubMtx object
 
   created -- 98feb20
   ---------------------------------------------------
*/
void
FrontMtx_storeSolution (
   FrontMtx        *frontmtx,
   int             owners[],
   int             myid,
   SubMtxManager   *manager,
   SubMtx          *p_mtx[],
   DenseMtx        *solmtx,
   int             msglvl,
   FILE            *msgFile
) ;
/*
   ---------------------------------------------------
   purpose -- load the dequeue with the roots of the
              active subtree used for a backward solve
 
   created -- 98mar27, cca
   ---------------------------------------------------
*/
void
FrontMtx_loadActiveRoots (
   FrontMtx   *frontmtx,
   char        status[],
   char        activeFlag,
   Ideq        *dequeue
) ;
/*
   --------------------------------------------------------------------
   purpose -- to set up the link data structures for a forward solve
 
   return value -- pointer to IP *heads[nfront+2], which contains
      the beginning of a list of IP objects that store the remaining
      updates to the fronts.
      note, heads[nfront] is the first IP object in the free list.
      heads[nfront+1] is the base address of the allocated IP objects.
 
   created -- 98feb20, cca
   --------------------------------------------------------------------
*/
IP **
FrontMtx_forwardSetup (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) ;
/*
   ---------------------------------------------------------
   purpose -- set up the linked lists for the backward solve
 
   created -- 98feb20, cca
   ---------------------------------------------------------
*/
IP **
FrontMtx_backwardSetup (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in QRfactor.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   purpose -- to compute the (U+I)D(I+U) factorization of A^TA,
              where A = QR, R = (D^{1/2} + D^{1/2}U)

   cpus[0] -- setup time
   cpus[1] -- initialize and load staircase matrix
   cpus[2] -- permute the rows into staircase form
   cpus[3] -- factor the matrix
   cpus[4] -- scale and store factor entries
   cpus[5] -- store update entries
   cpus[6] -- miscellaneous time
 
   created -- 98may28, cca
   ------------------------------------------------------------
*/
void
FrontMtx_QR_factor (
   FrontMtx     *frontmtx,
   InpMtx       *mtxA,
   ChvManager   *chvmanager,
   double       cpus[],
   double       *pfacops,
   int          msglvl,
   FILE         *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in QRsolve.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------
   minimize ||b - Ax||_2 by
   solving (U^T + I) * D * (I + U) X = A^T * B
      where A = QR = QD(I + U)
   by calling FrontMtx_solve()
 
   mtxmanager -- object that manages working storage
   cpus -- vector of cpu time breakdowns
      cpus[0] -- set up solves
      cpus[1] -- fetch rhs and store solution
      cpus[2] -- forward solve
      cpus[3] -- diagonal solve
      cpus[4] -- backward solve
      cpus[5] -- total solve time
      cpus[6] -- time to compute A^T * B
      cpus[7] -- total time
 
   created  -- 97may27, dkw
   modified -- 98may28, cca
   -------------------------------------------------
*/
void
FrontMtx_QR_solve (
   FrontMtx        *frontmtx,
   InpMtx          *mtxA,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   SubMtxManager   *manager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in QRutil.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   purpose -- to setup two data structures for a QR serial
              or multithreaded factorization
 
   rowsIVL[J]    -- list of rows of A to be assembled into front J
   firstnz[irow] -- column with location of leading nonzero of row in A
 
   created -- 98may29, cca
   --------------------------------------------------------------------
*/
void
FrontMtx_QR_setup (
   FrontMtx   *frontmtx,
   InpMtx     *mtxA,
   IVL        **prowsIVL,
   int        **pfirstnz,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   -----------------------------------------------
   purpose --  visit front J during a serial 
               or multithreaded QR factorization
 
   cpus[1] -- initialize and load staircase matrix
   cpus[2] -- factor the matrix
   cpus[3] -- scale and store factor entries
   cpus[4] -- store update entries
 
   created -- 98may28, cca
   -----------------------------------------------
*/
void
FrontMtx_QR_factorVisit (
   FrontMtx     *frontmtx,
   int          J,
   InpMtx       *mtxA,
   IVL          *rowsIVL,
   int          firstnz[],
   ChvList      *updlist,
   ChvManager   *chvmanager,
   char         status[],
   int          colmap[],
   DV           *workDV,
   double       cpus[],
   double       *pfacops,
   int          msglvl,
   FILE         *msgFile
) ;
/*
   --------------------------------------------------------------
   purpose -- create and return an A2 object that contains rows
              of A and rows from update matrices of the children.
              the matrix may not be in staircase form
 
   created -- 98may25, cca
   --------------------------------------------------------------
*/
A2 *
FrontMtx_QR_assembleFront (
   FrontMtx   *frontmtx,
   int        J,
   InpMtx     *mtxA,
   IVL        *rowsIVL,
   int        firstnz[],
   int        colmap[],
   Chv        *firstchild,
   DV         *workDV,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------------------
   store the factor entries of the reduced front matrix
 
   created -- 98may25, cca
   ----------------------------------------------------
*/
void
FrontMtx_QR_storeFront (
   FrontMtx   *frontmtx,
   int        J,
   A2         *frontJ,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   -------------------------------------------------
   purpose -- to create and return a Chv object that
              holds the update matrix for front J
 
   created -- 98may25, cca
   -------------------------------------------------
*/
Chv *
FrontMtx_QR_storeUpdate (
   FrontMtx     *frontmtx,
   int          J,
   A2           *frontJ,
   ChvManager   *chvmanager,
   int          msglvl,
   FILE         *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------
   purpose -- produce a map from each column
              to the front that contains it
 
   created -- 98may04, cca
   -----------------------------------------
*/
IV *
FrontMtx_colmapIV (
   FrontMtx   *frontmtx
) ;
/*
   --------------------------------------------------------------------
   purpose -- produce a map from each row to the front that contains it
 
   created -- 98may04, cca
   --------------------------------------------------------------------
*/
IV *
FrontMtx_rowmapIV (
   FrontMtx   *frontmtx
) ;
/*
   -------------------------------------------------------------
   compute the inertia of a symmetric matrix
 
   fill *pnnegative with the number of negative eigenvalues of A
   fill *pnzero     with the number of zero eigenvalues of A
   fill *pnpositive with the number of positive eigenvalues of A
 
   created -- 98may04, cca
   -------------------------------------------------------------
*/
void
FrontMtx_inertia (
   FrontMtx   *frontmtx,
   int         *pnnegative,
   int         *pnzero,
   int         *pnpositive
) ;
/*
   -------------------------------------------------------
   purpose -- create and return an IV object that contains 
              all the row ids owned by process myid.
 
   created -- 98jun13, cca
   -------------------------------------------------------
*/
IV *
FrontMtx_ownedRowsIV (
   FrontMtx   *frontmtx,
   int        myid,
   IV         *ownersIV,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   -------------------------------------------------------
   purpose -- create and return an IV object that contains 
              all the column ids owned by process myid.
 
   created -- 98jun13, cca
   -------------------------------------------------------
*/
IV *
FrontMtx_ownedColumnsIV (
   FrontMtx   *frontmtx,
   int        myid,
   IV         *ownersIV,
   int        msglvl,
   FILE       *msgFile
) ;
/*
   ----------------------------------------------------------------
  purpose -- to create and return an IVL object that holds the
      submatrix nonzero pattern for the upper triangular factor.
 
   NOTE: this method supercedes calling IVL_mapEntries() on
         the column adjacency structure. that gave problems when
         pivoting forced some fronts to have no eliminated columns.
        in some cases, solve aggregates were expected when none
         were forthcoming.
 
   created -- 98aug20, cca
   ----------------------------------------------------------------
*/
IVL *
FrontMtx_makeUpperBlockIVL (
   FrontMtx   *frontmtx,
   IV         *colmapIV
) ;
/*
   ----------------------------------------------------------------
  purpose -- to create and return an IVL object that holds the
      submatrix nonzero pattern for the lower triangular factor.
 
   NOTE: this method supercedes calling IVL_mapEntries() on
         the row adjacency structure. that gave problems when
         pivoting forced some fronts to have no eliminated columns.
        in some cases, solve aggregates were expected when none
         were forthcoming.
 
   created -- 98aug20, cca
   ----------------------------------------------------------------
*/
IVL *
FrontMtx_makeLowerBlockIVL (
   FrontMtx   *frontmtx,
   IV         *rowmapIV
) ;
/*
   ---------------------------------------------------------------
   purpose -- to compute and return the number of floating point
      operations to perform a solve with a single right hand side.
 
   created -- 98sep26, cca
   ---------------------------------------------------------------
*/
int
FrontMtx_nSolveOps (
   FrontMtx   *frontmtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   purpose -- to read in an FrontMtx object from a file
 
   input --
 
      fn -- filename, must be *.frontmtxb or *.frontmtxf
 
   return value -- 1 if success, 0 if failure
 
   created -- 98may04, cca
   -----------------------------------------------------
*/
int
FrontMtx_readFromFile (
   FrontMtx   *frontmtx,
   char       *fn
) ;
/*
   ------------------------------------------------------------
   purpose -- to read an FrontMtx object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 98may04, cca
   ------------------------------------------------------------
*/
int
FrontMtx_readFromFormattedFile (
   FrontMtx   *frontmtx,
   FILE        *fp
) ;
/*
   --------------------------------------------------------
   purpose -- to read an FrontMtx object from a binary file
 
   return value -- 1 if success, 0 if failure
 
   created -- 98may04, cca
   --------------------------------------------------------
*/
int
FrontMtx_readFromBinaryFile (
   FrontMtx   *frontmtx,
   FILE        *fp
) ;
/*
   -------------------------------------------------
   purpose -- to write an FrontMtx object to a file
 
   input --
 
      fn -- filename
        *.frontmtxb -- binary
        *.frontmtxf -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98may04, cca
   -------------------------------------------------
*/
int
FrontMtx_writeToFile (
   FrontMtx   *frontmtx,
   char        *fn
) ;
/*
   -----------------------------------------------------------
   purpose -- to write an FrontMtx object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98may04, cca
   -----------------------------------------------------------
*/
int
FrontMtx_writeToFormattedFile (
   FrontMtx   *frontmtx,
   FILE       *fp
) ;
/*
   -------------------------------------------------------
   purpose -- to write an FrontMtx object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98may04, cca
   -------------------------------------------------------
*/
int
FrontMtx_writeToBinaryFile (
   FrontMtx   *frontmtx,
   FILE        *fp
) ;
/*
   ---------------------------------------------------------------
   purpose -- to write out the statistics for the FrontMtx object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98may04, cca
   ---------------------------------------------------------------
*/
int
FrontMtx_writeStats (
   FrontMtx   *frontmtx,
   FILE        *fp
) ;
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form
 
   created -- 98may04, cca
   ----------------------------------------
*/
int
FrontMtx_writeForHumanEye (
   FrontMtx   *frontmtx,
   FILE       *fp
) ;
/*
   -------------------------------------------------------------
   purpose -- to write the factor matrices out for a matlab file
 
      Lname -- name for lower triangular matrix
      Dname -- name for diagonal matrix
      Uname -- name for upper triangular matrix
 
   presently works only with 1-d data decomposition
 
   created -- 98sep23, cca
   -------------------------------------------------------------
*/
int
FrontMtx_writeForMatlab (
   FrontMtx   *frontmtx,
   char       *Lname,
   char       *Dname,
   char       *Uname,
   FILE       *fp
) ;
/*--------------------------------------------------------------------*/
