/*  spoolesMPI.h  */

#ifndef _SPOOLES_MPI_

#define _SPOOLES_MPI_
#include <mpi.h>
#include "../FrontMtx.h"
#include "../misc.h"
#include "../timings.h"
#include "../SymbFac.h"
#include "../ETree.h"
#include "../SolveMap.h"
#include "../IVL.h"
#include "../IV.h"
#include "../Ideq.h"
#include "../Pencil.h"
#include "../Drand.h"

#endif

/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in IVallgather.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------------
   purpose -- 
 
   the IV objects objIV and ownersIV are found on each process.
   the ownersIV object is identical over all the processes, and
   owners[ii] tells which processes owns location ii of the obj[]
   vector. on return from this entry, the obj[] vector is replicated
   over all the processes. each process sends the (ii,obj[ii]) pairs
   that it owns to all the other processes.
 
   created -- 98apr02, cca
   -----------------------------------------------------------------
*/
void
IV_MPI_allgather (
   IV         *objIV,
   IV         *ownersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        tag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in IVLallgather.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   purpose -- 
 
   the IVL object ivl and IV object ownersIV are both found on 
   each process.  the ownersIV object is identical over all the 
   processes, and owners[ii] tells which processes owns list ii 
   of the ivl object. on return from this method, the ivl object 
   is replicated over all the processes. each process sends 
   the lists that it owns to all the other processes.
 
   created -- 98apr03, cca
   -------------------------------------------------------------
*/
void
IVL_MPI_allgather (
   IVL        *ivl,
   IV         *ownersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        tag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in IVL_alltoall.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   this method is used during the setup for matrix-vector multiplies.
   each processor has computed the vertices it needs from other
   processors, these lists are contained in sendIVL. on return,
   recvIVL contains the lists of vertices this processor must send
   to all others.
 
   sendIVL -- on input, list[q] contains the vertices needed by
              this processor that are owned by q
   recvIVL -- on output, list[q] contains the vertices owned by
              this processor that are needed by q. note, if NULL
              on input, a new IVL object is allocated
   stats[] -- statistics vector
     stats[0] -- contains # of sends
     stats[1] -- contains # of receives
     stats[2] -- contains # of bytes sent
     stats[3] -- contains # of bytes received
   firsttag -- first tag for messages,
     tags in range [firsttag, firsttag+nproc-1] are used
 
   return value -- recvIVL
 
   created -- 98jul26, cca
   ------------------------------------------------------------------
*/
IVL *
IVL_MPI_alltoall (
   IVL        *sendIVL,
   IVL        *recvIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in aggListMPI.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   create, initialize and return a ChvList object
   to deal with aggregate chevrons
      
   created  -- 98may21, cca
   -----------------------------------------------
*/
ChvList *
FrontMtx_MPI_aggregateList (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        tag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in splitDenseMtx.c ----------------------------------
------------------------------------------------------------------------
*/
/*

-----------------------------------------------------------------
   purpose -- to split a DenseMtx object by rows
 
   mtx         -- DenseMtx object
   rowmapIV    -- map from rows to owning processes
   firsttag    -- first tag to be used in these messages
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl      -- message level
   msgFile     -- message file
   comm        -- MPI communicator
 
   return value -- a new DenseMtx object filled with the owned rows

 
   created  -- 98may16, cca

-----------------------------------------------------------------
*/
DenseMtx *
DenseMtx_MPI_splitByRows (
   DenseMtx   *mtx,
   IV         *rowmapIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*
   -----------------------------------------------------------------
   purpose -- to scatter a DenseMtx object by rows
 
   Xglobal     -- global DenseMtx object, significant only for root
   Xlocal      -- local DenseMtx object, if NULL on input, will
                  be created if necessary
   rowmapIV    -- map from rows to owning processes
   firsttag    -- first tag to be used in these messages
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl      -- message level
   msgFile     -- message file
   comm        -- MPI communicator
 
   return value -- Xlocal, a local DenseMtx object, may be NULL
 
   created  -- 98may16, cca
   modified -- 98sep26, cca
      mtx is not modified
   -----------------------------------------------------------------
*/
DenseMtx *
DenseMtx_MPI_splitFromGlobalByRows (
   DenseMtx   *Xglobal,
   DenseMtx   *Xlocal,
   IV         *rowmapIV,
   int        root,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*
   -----------------------------------------------------------------
   purpose -- to merge a DenseMtx object by rows
 
   Xlocal      -- DenseMtx object, can be NULL
   Xglobal     -- DenseMtx object, can be NULL
                  significant only for root
   firsttag    -- first tag to be used in these messages
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl      -- message level
   msgFile     -- message file
   comm        -- MPI communicator
 
   return value -- 
      if processor is root
          Xglobal is returned, if was NULL on input, it is created
      else
          NULL
      endif
 
   Xlocal is not modified
 
   created  -- 98sep27, cca
   -----------------------------------------------------------------
*/
DenseMtx *
DenseMtx_MPI_mergeToGlobalByRows (
   DenseMtx   *Xglobal,
   DenseMtx   *Xlocal,
   int        root,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in splitInpMtx.c -----------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   purpose -- to split a distributed InpMtx object
 
   inpmtx     -- pointer to the local InpMtx object
   mapIV      -- pointer to the map from vertices to processes
   firsttag   -- first tag value, one will be used
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl     -- local message level
   msgFile    -- local message file
   comm       -- MPI communication structure
 
   return value -- pointer to the new InpMtx object
      that contains the owned entries.
 
   created  -- 97jun20, cca
   modified -- 97oct17, cca
      stats added
   ------------------------------------------------------------
*/
InpMtx *
InpMtx_MPI_split (
   InpMtx     *inpmtx,
   IV         *mapIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*
   ------------------------------------------------------------
   purpose -- to split a distributed InpMtx object
 
   Aglobal  -- pointer to the global InpMtx object,
               significant only for root processor
   Alocal   -- pointer to the local InpMtx object,
               if NULL and will be nonempty, it will be created
   mapIV    -- pointer to the map from vertices to processes
   root     -- processor that presently owns the global matrix
   firsttag -- first tag value, one will be used
   stats[4] -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl  -- local message level
   msgFile -- local message file
   comm    -- MPI communication structure
 
   return value -- pointer to the local InpMtx object 
      that contains the owned entries.
 
   created  -- 98may16, cca
   modified -- 98sep26, cca
      inpmtx is not modified
   ------------------------------------------------------------
*/
InpMtx *
InpMtx_MPI_splitFromGlobal (
   InpMtx     *Aglobal,
   InpMtx     *Alocal,
   IV         *mapIV,
   int        root,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in splitPencil.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   purpose -- to split a distributed Pencil object
 
   pencil     -- pointer to the local Pencil object
   mapIV      -- pointer to the map from vertices to processes
   firsttag   -- first tag value, two will be used, tag and tag+1
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl     -- local message level
   msgFile    -- local message file
   comm       -- MPI communication structure
 
   created  -- 98may20, cca
   --------------------------------------------------------------
*/
void
Pencil_MPI_split (
   Pencil     *pencil,
   IV         *mapIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in symbfacMPI.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   perform the symbolic factorization in parallel
 
   input --
      frontETree    -- front tree for the factorization
      frontOwnersIV -- map from fronts to owning process
      pencil        -- matrix pencil object
      firsttag      -- first tag to be used for messages,
                       will use tag, ..., tag + nfront - 1
      msglvl        -- message level
      msgFile       -- message file
      comm          -- MPI communicator
 
   return value --
      symbfacIVL -- contains front indices for supported fronts
 
   created -- 98may20, cca
   ------------------------------------------------------------
*/
IVL *
SymbFac_MPI_initFromPencil (
   ETree      *frontETree,
   IV         *frontOwnersIV,
   Pencil     *pencil,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*
   ------------------------------------------------------------
   perform the symbolic factorization in parallel
 
   input --
      frontETree    -- front tree for the factorization
      frontOwnersIV -- map from fronts to owning process
      inpmtx        -- matrix object
      firsttag      -- first tag to be used for messages,
                       will use tag, ..., tag + nfront - 1
      msglvl        -- message level
      msgFile       -- message file
      comm          -- MPI communicator
 
   return value --
      symbfacIVL -- contains front indices for supported fronts
 
   created -- 98may20, cca
   ------------------------------------------------------------
*/
IVL *
SymbFac_MPI_initFromInpMtx (
   ETree      *frontETree,
   IV         *frontOwnersIV,
   InpMtx     *inpmtx,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in postProcess.c ------------------------------------
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
 
   created -- 98may20, cca
   --------------------------------------------------------------
*/
void
FrontMtx_MPI_postProcess (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*
   -------------------------------------------------------------
   purpose -- to permute the indices of the upper adjacency
      structure so that for each front J, bnd{J} is in 
      ascending order w.r.t. K cup bnd{K} for par[J] = K.
      process q sends to process r one message that contains
      J cup bnd{J} for all J owned by q and needed by r.
      once all the indices for the supported fronts are present,
      the indices in the upper adjacency structure are reordered
      as necessary.
 
   created -- 98may20, cca
   -------------------------------------------------------------
*/
void
FrontMtx_MPI_permuteUpperAdj (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*
   -------------------------------------------------------------
   purpose -- to permute the indices of the lower adjacency
      structure so that for each front J, bnd{J} is in
      ascending order w.r.t. K cup bnd{K} for par[J] = K.
      process q sends to process r one message that contains
      J cup bnd{J} for all J owned by q and needed by r.
      once all the indices for the supported fronts are present,
      the indices in the lower adjacency structure are reordered
      as necessary.
 
   created -- 98may20, cca
   -------------------------------------------------------------
*/
void
FrontMtx_MPI_permuteLowerAdj (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in fullAdjMPI.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- given a distributed InpMtx object,
      gather all the indices onto each process and create
      an IVL object that contains the full adjacency structure
 
   created -- 97dec17, cca
   -----------------------------------------------------------
*/
IVL *
InpMtx_MPI_fullAdjacency (
   InpMtx    *inpmtx,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) ;
/*
   -----------------------------------------------------------
   purpose -- given a distributed Pencil object,
      gather all the indices onto each process and create
      an IVL object that contains the full adjacency structure
 
   created -- 97dec18, cca
   -----------------------------------------------------------
*/
IVL *
Pencil_MPI_fullAdjacency (
   Pencil    *pencil,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in factorMPI.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   this is the method that computes the
   parallel factorization of a sparse matrix A using MPI.
 
   input ---
 
      frontmtx      -- front matrix object that holds the factors
      inpmtx        -- matrix object for A
      tau           -- tolerance used in pivoting
      droptol       -- tolerance used in the approximate factorization
      chvmanager    -- ChvManager object to handle working storage
      frontOwnersIV -- map from fronts to owning processes
      lookahead     -- parameter that governs the ``lookahead''
                       behavior, use lookahead = 0 as a default
      cpus          -- vector to store breakdown of cpu times
         cpus[ 0] -- initialize fronts
         cpus[ 1] -- load original entries
         cpus[ 2] -- update fronts
         cpus[ 3] -- insert aggregate data
         cpus[ 4] -- assemble aggregate data
         cpus[ 5] -- assemble postponed data
         cpus[ 6] -- factor fronts
         cpus[ 7] -- extract postponed data
         cpus[ 8] -- store factor entries
         cpus[ 9] -- post initial receives
         cpus[10] -- check for received messages
         cpus[11] -- post initial sends
         cpus[12] -- check for sent messages
      stats         -- vector to store statistics
         stats[ 0] -- # of aggregate sends
         stats[ 1] -- # of bytes in the aggregate sends
         stats[ 2] -- # of aggregate received
         stats[ 3] -- # of bytes in the aggregate received
         stats[ 4] -- # of postponed data sends
         stats[ 5] -- # of bytes in the postponed data sends
         stats[ 6] -- # of postponed data received
         stats[ 7] -- # of bytes in the postponed data received
         stats[ 8] -- # of active Chv objects (working storage)
         stats[ 9] -- # of active bytes in working storage
         stats[10] -- # of requested bytes in working storage
      msglvl        -- message level
      msgFile       -- message file
      firsttag      -- first tag to use during the factorization,
                       reserved tags are tag, ..., tag + 3*nfront + 2
      comm          -- MPI communicator
 
   return value -- if process id is zero, a pointer to the first
   Chv object in a list that contains postponed rows and columns
   that could not be eliminated.
 
   created  -- 98may21, cca
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_MPI_factorInpMtx (
   FrontMtx     *frontmtx,
   InpMtx       *inpmtx,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   IV           *frontOwnersIV,
   int          lookahead,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile,
   int          firsttag,
   MPI_Comm     comm
) ;
/*
   -------------------------------------------------------------------
   this is the method that computes the 
   parallel factorization of A + sigma*B using MPI.
 
   input ---
 
      frontmtx      -- front matrix object that holds the factors
      pencil        -- matrix pencil object, contains A + sigma*B
      tau           -- tolerance used in pivoting
      droptol       -- tolerance used in the approximate factorization
      chvmanager    -- ChvManager object to handle working storage
      frontOwnersIV -- map from fronts to owning processes
      lookahead     -- parameter that governs the ``lookahead''
                       behavior, use lookahead = 0 as a default
      cpus          -- vector to store breakdown of cpu times
         cpus[ 0] -- initialize fronts
         cpus[ 1] -- load original entries
         cpus[ 2] -- update fronts
         cpus[ 3] -- insert aggregate data
         cpus[ 4] -- assemble aggregate data
         cpus[ 5] -- assemble postponed data
         cpus[ 6] -- factor fronts
         cpus[ 7] -- extract postponed data
         cpus[ 8] -- store factor entries
         cpus[ 9] -- post initial receives
         cpus[10] -- check for received messages
         cpus[11] -- post initial sends
         cpus[12] -- check for sent messages
      stats         -- vector to store statistics
         stats[ 0] -- # of aggregate sends
         stats[ 1] -- # of bytes in the aggregate sends
         stats[ 2] -- # of aggregate received
         stats[ 3] -- # of bytes in the aggregate received
         stats[ 4] -- # of postponed data sends
         stats[ 5] -- # of bytes in the postponed data sends
         stats[ 6] -- # of postponed data received
         stats[ 7] -- # of bytes in the postponed data received
         stats[ 8] -- # of active Chv objects (working storage)
         stats[ 9] -- # of active bytes in working storage
         stats[10] -- # of requested bytes in working storage
      msglvl        -- message level
      msgFile       -- message file
      firsttag      -- first tag to use during the factorization,
                       reserved tags are tag, ..., tag + 3*nfront + 2
      comm          -- MPI communicator
 
   return value -- if process id is zero, a pointer to the first
   Chv object in a list that contains postponed rows and columns
   that could not be eliminated.
 
   created  -- 98may21, cca
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_MPI_factorPencil (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   IV           *frontOwnersIV,
   int          lookahead,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile,
   int          firsttag,
   MPI_Comm     comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in splitFrontMtx.c ----------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   purpose -- after the factorization has been computed and the
      front matrices have been split into submatrices, and after
      a solve map object has been computed, the submatrices
      are sent to the process that owns them.
 
   frontmtx -- stores the factor matrix
   solvemap -- stores the map from submatrices to processes
 
   created -- 98may21, cca
   -------------------------------------------------------------
*/
void
FrontMtx_MPI_split (
   FrontMtx   *frontmtx,
   SolveMap   *solvemap,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in Graph_Bcast.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- to broadcast a Graph IVL object from 
              one processor to all the others
 
   created -- 98sep10, cca
   -----------------------------------------------
*/
Graph *
Graph_MPI_Bcast (
   Graph      *graph,
   int        root,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in IVL_Bcast.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------
   purpose -- to broadcast an IVL object from
              one processor to all the others
 
   created -- 98sep10, cca
   ------------------------------------------
*/
IVL *
IVL_MPI_Bcast (
   IVL        *ivl,
   int        root,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in ETreeBcast.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------
   purpose -- to broadcast a front tree object 
              from one process to all the others
 
   created -- 98may21, cca
   ---------------------------------------------
*/
ETree *
ETree_MPI_Bcast (
   ETree      *etree,
   int        root,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in rowmapMPI.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- after pivoting for a nonsymmetric factorization,
              some delayed rows may belong to a process other
              than its original owner. this method returns an
              IV object that maps rows to owning processes.
 
   created -- 98may22, cca
   -----------------------------------------------------------
*/
IV *
FrontMtx_MPI_rowmapIV (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in colmapMPI.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   purpose -- after pivoting for a nonsymmetric factorization,
              some delayed columns may belong to a process other
              than its original owner. this method returns an
              IV object that maps columns to owning processes.
 
   created -- 98may22, cca
   -------------------------------------------------------------
*/
IV *
FrontMtx_MPI_colmapIV (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in utilities.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------
   return the maximum tag value
 
   created -- 98jan08, cca
   ----------------------------
*/
int
maxTagMPI (
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in solveMPI.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   MPI solve method for (L + I)D(I + U)X = B or (U^T + I)D(I + U)X = B
 
   created -- 98may21, ca
   -------------------------------------------------------------------
*/
void
FrontMtx_MPI_solve (
   FrontMtx        *frontmtx,
   DenseMtx        *solmtx,
   DenseMtx        *rhsmtx,
   SubMtxManager   *mtxmanager,
   SolveMap        *solvemap,
   double          cpus[],
   int             stats[],
   int             msglvl,
   FILE            *msgFile,
   int             firstTag,
   MPI_Comm        comm
) ;
/*--------------------------------------------------------------------*/
#define MMM_WITH_A   0
#define MMM_WITH_AT  1
#define MMM_WITH_AH  2
#define MMM_LOCAL    1
#define MMM_GLOBAL   2
/*
   --------------------------------------------------------------------
   symflag -- symmetry flag
      SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC
   opflag -- operations flag
      0 (MMM_WITH_A)  -- multiply with A
      1 (MMM_WITH_AT) -- multiply with A^T
      2 (MMM_WITH_AH) -- multiply with A^H
   XownedIV -- list of rows of X owned by this processor
   XsupIV   -- vector of rows of X that are supported by this processor
   XmapIV   -- vector that maps global rows of X to local rows of X
               that are supported by this processor,
       if jj = Xsup[ii] then
          Xmap[jj] = ii 
       else
          Xmap[jj] = -1
       endif
   XsendIVL -- list jproc contains the local entries of Xloc that get
               sent to processor jproc.
   XrecvIVL -- list jproc contains the local entries of Xsupp that will
               be received from processor jproc.
   YownedIV -- list of rows of Y owned by this processor
   Xsupp    -- DenseMtx object used in the matrix-matrix multiply
   YsupIV -- vector of rows of Y that are supported by this processor
   YmapIV -- vector that maps global rows of Y to local rows of Y
             that are supported by this processor,
       if jj = Ysup[ii] then
          Ymap[jj] = ii 
       else
          Ymap[jj] = -1
       endif
   YsendIVL -- list jproc contains the local entries of Yloc that get
               sent to processor jproc.
   YrecvIVL -- list jproc contains the local entries of Ysupp that will
               be received from processor jproc.
   Ysupp    -- DenseMtx object used in the matrix-matrix multiply

   created -- 98aug21, cca
   --------------------------------------------------------------------
*/
typedef struct _MatMulInfo   MatMulInfo ;
struct _MatMulInfo {
   int        symflag   ;
   int        opflag    ;
   IV         *XownedIV ;
   IV         *XsupIV   ;
   IV         *XmapIV   ;
   IVL        *XsendIVL ;
   IVL        *XrecvIVL ;
   DenseMtx   *Xsupp    ;
   IV         *YownedIV ;
   IV         *YsupIV   ;
   IV         *YmapIV   ;
   IVL        *YsendIVL ;
   IVL        *YrecvIVL ;
   DenseMtx   *Ysupp    ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in MMM.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------
   purpose -- setup the distributed matrix-matrix multiply
              information object
 
   created -- 98aug21, cca
   -------------------------------------------------------
*/
MatMulInfo *
MatMul_MPI_setup (
   InpMtx     *A,
   int        symflag,
   int        opflag,
   IV         *XownersIV,
   IV         *YownersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        tag,
   MPI_Comm   comm
) ;
/*
   ------------------------------------------------------------
   set the indices of A to be local with respect to its support
 
   created -- 98aug21, cca
   ------------------------------------------------------------
*/
void
MatMul_setLocalIndices (
   MatMulInfo   *info,
   InpMtx       *A
) ;
/*
   -------------------------------------------------------------
   set the indices of A to be global with respect to its support
 
   created -- 98aug21, cca
   -------------------------------------------------------------
*/
void
MatMul_setGlobalIndices (
   MatMulInfo   *info,
   InpMtx       *A
) ;
/*
   ---------------------------------------------------------
   purpose -- compute the distributed matrix-matrix multiply
              Y := Y - alpha * A * X
      where A, Xloc and Yloc are distributed
 
   created -- 98aug21, cca
   ---------------------------------------------------------
*/
void
MatMul_MPI_mmm (
   MatMulInfo   *info,
   DenseMtx     *Yloc,
   double       alpha[],
   InpMtx       *A,
   DenseMtx     *Xloc,
   int          stats[],
   int          msglvl,
   FILE         *msgFile,
   int          tag,
   MPI_Comm     comm
) ;
/*
   --------------------------------------------------
   free the MatMulInfo object and its data structures 
 
   created -- 98aug21, cca
   --------------------------------------------------
*/
void
MatMul_cleanup (
   MatMulInfo   *info
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in DenseMtx_gather.c --------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   purpose -- to gather entries from a distributed X into Y .
      what rows of X to be sent to other processors are 
      found in sendIVL. what rows of Y to be received 
      from other processors are found in recvIVL.
 
   Y  -- on return, contains the rows specified by recvIVL.
      row indices of Y will be in ascending order.
   X  -- this processor's part of the distributed partitioned
      DenseMtx object. row indices of X are assumed to be 
      in ascending order.
   sendIVL -- list jproc contains the global ids of rows in X
      that need to be sent to processor jproc. note, lists are
      assumed to be in ascending order and are local with respect to X.
   recvIVL -- list jproc contains the global ids of rows in jproc's
     part of X that will to be sent to this processor. note,
      lists are assumed to be in ascending order and are local
      with respect to Y.
 
   created -- 98jul31, cca
   --------------------------------------------------------------------
*/
void
DenseMtx_MPI_gatherRows (
   DenseMtx   *Y,
   DenseMtx   *X,
   IVL        *sendIVL,
   IVL        *recvIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in DenseMtx_scatterAdd.c ----------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   purpose -- to scatter/add entries from a distributed X into Y.
      what rows of X to be sent to other processors are 
      found in sendIVL. what rows of Y to be received 
      from other processors are found in recvIVL.
 
   Y  -- on return, contains the rows specified by recvIVL.
      row indices of Y will be in ascending order.
   X  -- this processor's part of the distributed partitioned
      DenseMtx object. row indices of X are assumed to be 
      in ascending order.
   sendIVL -- list jproc contains the global ids of rows in X
      that need to be sent to processor jproc. note, lists are
      assumed to be in ascending order and are local with respect to X.
   recvIVL -- list jproc contains the global ids of rows in jproc's
     part of X that will to be sent to this processor. note,
      lists are assumed to be in ascending order and are local
      with respect to Y.
 
   created -- 98jul31, cca
   --------------------------------------------------------------------
*/
void
DenseMtx_MPI_scatterAddRows (
   DenseMtx   *Y,
   DenseMtx   *X,
   IVL        *sendIVL,
   IVL        *recvIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---- methods found in makeSendRecvIVLs.c -------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------------
   purpose -- to analyze and organize communication. it was written
     in support of a distributed matrix-vector multiply but can be
     used for other applications.
 
   each processor has a list of items it "supports" or needs found
   in the supportedIV object. the globalmapIV object contains the
   map from items to owning processors. we need to figure out what
   items this processor will send to and receive from each other
   processor. this information is found in the sendIVL and recvIVL
   objects. 
 
   on return, list jproc of sendIVL contains the items owned by
   this processor and needed by jproc.
   on return, list jproc of recvIVL contains the items needed by
   this processor and owned by jproc.
 
   as a concrete example, consider a distributed Y = A * X.
   the matrix A, the right hand side X and the vector Y are
   distributed among processors. 
 
   consider the case where the supportedIV object contains the rows
   of X that are needed by this processor to perform its part of the
   matrix-vector multiply. globalmapIV contains the map from rows
   of X to the owning processors. on return, list jproc of sendIVL
   contains the row indices of X owned by this processor that are
   needed by processor jproc. on return, list jproc of recvIVL 
   contains the row indices of X needed by this processor that are
   owned by processor jproc.
 
   consider the case where the supportedIV object contains the rows
   of Y that will be updated by this processor when it performs it
   part of the matrix-vector multiply. globalmapIV contains the map
   from rows of Y to their owning processors. on return, list jproc
   of recvIVL contains the row indices of Y on this processor that
   need to be sent to their owner, processor jproc. on return, list
   jproc of sendIVL contains the row indices of Y owned by this 
   processor that will be sent by processor jproc to this
   processor.
 
   created -- 98aug01, cca
   -----------------------------------------------------------------
*/
void
makeSendRecvIVLs (
   IV         *supportedIV,
   IV         *globalmapIV,
   IVL        *sendIVL,
   IVL        *recvIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) ;
/*--------------------------------------------------------------------*/
