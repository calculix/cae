/*  SetupMPI.c  */

#include "../BridgeMPI.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_Setup = 0 ;
static double time_Setup = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose --

   given InpMtx objects that contain A and B, initialize the bridge
   data structure for the MPI factor's, solve's and mvm's.
  
   NOTE: all the input arguments are pointers 
         to allow calls from Fortran

   data -- pointer to a Bridge object
   pprbtype -- pointer to value containing problem type
     *prbtype = 1 --> A X = B X Lambda, vibration problem
     *prbtype = 2 --> A X = B X Lambda, buckling problem
     *prbtype = 3 --> A X = X Lambda, simple eigenvalue problem
   pneqns  -- pointer to value containing number of equations
   pmxbsz  -- pointer to value containing blocksize
   A       -- pointer to InpMtx object containing A
   B       -- pointer to InpMtx object containing B
   pseed   -- pointer to value containing a random number seed
   pmsglvl -- pointer to value containing a message level
   msgFile -- message file pointer

   return value --
      1 -- normal return
     -1 -- data is NULL
     -2 -- pprbtype is NULL
     -3 -- *pprbtype is invalid
     -4 -- pneqns is NULL
     -5 -- *pneqns is invalid
     -6 -- pmxbsz is NULL
     -7 -- *pmxbsz is invalid
     -8 -- A and B are NULL
     -9 -- pseed is NULL
    -10 -- pmsglvl is NULL
    -11 -- *pmsglvl > 0 and msgFile is NULL
    -12 -- comm is NULL

   created -- 98aug10, cca
   ----------------------------------------------------------------
*/
int
SetupMPI (
   void       *data,
   int        *pprbtype,
   int        *pneqns,
   int        *pmxbsz,
   InpMtx     *A,
   InpMtx     *B,
   int        *pseed,
   int        *pmsglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
BridgeMPI   *bridge = (BridgeMPI *) data ;
double      cutoff, minops ;
double      *opcounts ;
double      sigma[2] ;
DV          *cumopsDV ;
Graph       *graph ;
int         maxdomainsize, maxsize, maxzeros, msglvl, myid, mxbsz,
            nedges, neqns, nmyowned, nproc, prbtype, root, seed, tag ;
int         stats[4] ;
IVL         *adjIVL ;
#if MYDEBUG > 0
double      t1, t2 ;
#endif
/*
   ---------------------------------
   get number of processors and rank
   ---------------------------------
*/
MPI_Comm_rank(comm, &myid)  ;
MPI_Comm_size(comm, &nproc) ;
#if MYDEBUG > 0
count_Setup++ ;
MARKTIME(t1) ;
if ( myid == 0 ) {
   fprintf(stdout, "\n (%d) SetupMPI()", count_Setup) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, "\n (%d) SetupMPI()", count_Setup) ;
fflush(bridge->msgFile) ;
#endif
/*
   --------------------
   check the input data
   --------------------
*/
if ( data == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n data is NULL\n") ;
   return(-1) ;
}
if ( pprbtype == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n prbtype is NULL\n") ;
   return(-2) ;
}
prbtype = *pprbtype ;
if ( prbtype < 1 || prbtype > 3 ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n prbtype = %d, is invalid\n", prbtype) ;
   return(-3) ;
}
if ( pneqns == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n pneqns is NULL\n") ;
   return(-4) ;
}
neqns = *pneqns ;
if ( neqns <= 0 ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n neqns = %d, is invalid\n", neqns) ;
   return(-5) ;
}
if ( pmxbsz == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n pmxbsz is NULL\n") ;
   return(-6) ;
}
mxbsz = *pmxbsz ;
if ( mxbsz <= 0 ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n *pmxbsz = %d, is invalid\n", mxbsz) ;
   return(-7) ;
}
if ( A == NULL && B == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n A and B are NULL\n") ;
   return(-8) ;
}
if ( pseed == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n pseed is NULL\n") ;
   return(-9) ;
}
seed = *pseed ;
if ( pmsglvl == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n pmsglvl is NULL\n") ;
   return(-10) ;
}
msglvl = *pmsglvl ;
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n msglvl = %d, msgFile = NULL\n", msglvl) ;
   return(-11) ;
}
if ( comm == NULL ) {
   fprintf(stderr, "\n fatal error in SetupMPI()"
           "\n comm = NULL\n") ;
   return(-12) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside SetupMPI()"
           "\n neqns = %d, prbtype = %d, mxbsz = %d, seed = %d",
           neqns, prbtype, mxbsz, seed) ;
   if ( A != NULL ) {
      fprintf(msgFile, "\n\n matrix A") ;
      InpMtx_writeForHumanEye(A, msgFile) ;
   }
   if ( B != NULL ) {
      fprintf(msgFile, "\n\n matrix B") ;
      InpMtx_writeForHumanEye(B, msgFile) ;
   }
   fflush(msgFile) ;
}
bridge->myid    = myid    ;
bridge->nproc   = nproc   ;
bridge->prbtype = prbtype ;
bridge->neqns   = neqns   ;
bridge->mxbsz   = mxbsz   ;
bridge->A       = A       ;
bridge->B       = B       ;
bridge->seed    = seed    ;
bridge->msglvl  = msglvl  ;
bridge->msgFile = msgFile ;
bridge->comm    = comm    ;
/*
   ------------------------------------
   STEP 1: create and initialize pencil
   ------------------------------------
*/
sigma[0] = 1.0;
sigma[1] = 0.0;
bridge->pencil = Pencil_new() ;
Pencil_init(bridge->pencil, SPOOLES_REAL, SPOOLES_SYMMETRIC,
            A, sigma, B) ;
/*
   ----------------------------------------
   STEP 2: convert to row or column vectors
   ----------------------------------------
*/
if ( A != NULL ) {
   if ( ! INPMTX_IS_BY_ROWS(A) && ! INPMTX_IS_BY_COLUMNS(A) ) {
      InpMtx_changeCoordType(A, INPMTX_BY_ROWS) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(A) ) {
      InpMtx_changeStorageMode(A, INPMTX_BY_VECTORS) ;
   }
}
if ( B != NULL ) {
   if ( ! INPMTX_IS_BY_ROWS(B) && ! INPMTX_IS_BY_COLUMNS(B) ) {
      InpMtx_changeCoordType(B, INPMTX_BY_ROWS) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(B) ) {
      InpMtx_changeStorageMode(B, INPMTX_BY_VECTORS) ;
   }
}
/*
   ---------------------------------------
   STEP 3: create a Graph object for A + B
   ---------------------------------------
*/
graph  = Graph_new() ;
IVfill(4, stats, 0) ;
adjIVL = Pencil_MPI_fullAdjacency(bridge->pencil, 
                                  stats, msglvl, msgFile, comm);
nedges = IVL_tsize(adjIVL),
Graph_init2(graph, 0, bridge->neqns, 0, nedges,
            bridge->neqns, nedges, adjIVL, NULL, NULL) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n graph of the input pencil") ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------------
   STEP 4: order the graph. each processor orders the graph
           independently. then the processors determine who 
           has the best ordering, and that ordering information 
           (contained in the FrontETree object) is broadcast 
           to all processors.
   ------------------------------------------------------------
*/
maxdomainsize = neqns / 64 ;
if ( maxdomainsize == 0 ) {
   maxdomainsize = 1 ;
}
maxzeros  = (int) (0.01*neqns) ;
maxsize   = 64 ;
bridge->frontETree = orderViaBestOfNDandMS(graph, maxdomainsize,
                                 maxzeros, maxsize, bridge->seed + myid,
                                 msglvl, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree from ordering") ;
   fflush(msgFile) ;
}
Graph_free(graph) ;
opcounts = DVinit(nproc, 0.0) ;
opcounts[myid] = ETree_nFactorOps(bridge->frontETree, 
                                  SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n before gather \n" );
   fflush(msgFile) ;
}

MPI_Allgather((void *) &opcounts[myid], 1, MPI_DOUBLE,
              (void *) opcounts, 1, MPI_DOUBLE, comm) ;

if ( msglvl > 2 ) {
   fprintf(msgFile, "\n after gather \n" );
   fflush(msgFile) ;
}

minops = DVmin(nproc, opcounts, &root) ;
root = 0 ;
DVfree(opcounts) ;
bridge->frontETree = ETree_MPI_Bcast(bridge->frontETree, root, 
				     msglvl, msgFile, comm) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n best front tree") ;
   ETree_writeForHumanEye(bridge->frontETree, msgFile) ;
   fflush(msgFile) ;
}
/*
ETree_writeToFile(bridge->frontETree, "stk35.etreef") ;
*/
/*
   ------------------------------------------------------
   STEP 5: get the old-to-new and new-to-old permutations
   ------------------------------------------------------
*/
bridge->oldToNewIV = ETree_oldToNewVtxPerm(bridge->frontETree) ;
bridge->newToOldIV = ETree_newToOldVtxPerm(bridge->frontETree) ;
IV_writeToFile(bridge->oldToNewIV, "oldToNew.ivf") ;
IV_writeToFile(bridge->newToOldIV, "newToOld.ivf") ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n old-to-new permutation") ;
   fprintf(msgFile, "\n\n new-to-old permutation") ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   STEP 6: permute the vertices in the front tree
   ----------------------------------------------
*/
ETree_permuteVertices(bridge->frontETree, bridge->oldToNewIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted front etree") ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------------
   STEP 7: generate 
           (1) ownersIV  -- the map from fronts to processors
           (2) vtxmapIV  -- the map from vertices to processors
           (3) myownedIV -- vertices owned by this processor
   ------------------------------------------------------------
*/
cutoff   = 1./(2*nproc) ;
cumopsDV = DV_new() ;
DV_init(cumopsDV, nproc, NULL) ;
bridge->ownersIV = ETree_ddMap(bridge->frontETree, SPOOLES_REAL, 
                               SPOOLES_SYMMETRIC, cumopsDV, cutoff) ;
DV_free(cumopsDV) ;
bridge->vtxmapIV = IV_new() ;
IV_init(bridge->vtxmapIV, bridge->neqns, NULL) ;
IVgather(bridge->neqns, IV_entries(bridge->vtxmapIV), 
         IV_entries(bridge->ownersIV), 
         ETree_vtxToFront(bridge->frontETree)) ;
bridge->myownedIV = IV_targetEntries(bridge->vtxmapIV, myid) ;
nmyowned = IV_size(bridge->myownedIV) ;
bridge->rowmapIV     = NULL ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n map from fronts to owning processes") ;
   IV_writeForHumanEye(bridge->ownersIV, msgFile) ;
   fprintf(msgFile, "\n\n map from vertices to owning processes") ;
   IV_writeForHumanEye(bridge->vtxmapIV, msgFile) ;
   fprintf(msgFile, "\n\n vertices owned by this process") ;
   IV_writeForHumanEye(bridge->myownedIV, msgFile) ;
   fflush(msgFile) ;
}
/*
IV_writeToFile(bridge->ownersIV, "owners.ivf") ;
Tree_writeToFile(bridge->frontETree->tree, "stk35.treef") ;
*/
/*
   ---------------------------------------------------
   STEP 8: permute the entries in the pencil. 
           note, after the permutation the 
           entries are mapped into the upper triangle.
   ---------------------------------------------------
*/
Pencil_permute(bridge->pencil, bridge->oldToNewIV, bridge->oldToNewIV) ;
Pencil_mapToUpperTriangle(bridge->pencil) ;
Pencil_changeCoordType(bridge->pencil, INPMTX_BY_CHEVRONS) ;
Pencil_changeStorageMode(bridge->pencil, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted pencil") ;
   fflush(msgFile) ;
}
tag = 0 ;
Pencil_MPI_split(bridge->pencil, bridge->vtxmapIV, 
                 stats, msglvl, msgFile, tag, comm) ;
Pencil_changeCoordType(bridge->pencil, INPMTX_BY_CHEVRONS) ;
Pencil_changeStorageMode(bridge->pencil, INPMTX_BY_VECTORS) ;
bridge->A = A = bridge->pencil->inpmtxA ;
bridge->B = B = bridge->pencil->inpmtxB ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted and split pencil") ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------
   STEP 9: compute the symbolic factorization
   ------------------------------------------
*/
bridge->symbfacIVL = SymbFac_MPI_initFromPencil(bridge->frontETree, 
                                    bridge->ownersIV, bridge->pencil, 
                                    stats, msglvl, msgFile, tag, comm) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n symbolic factorization") ;
   IVL_writeForHumanEye(bridge->symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   STEP 10: create a FrontMtx object to hold the factorization
   -----------------------------------------------------------
*/
bridge->frontmtx = FrontMtx_new() ;
/*
   ---------------------------------------------------------------------
   STEP 11: create a SubMtxManager object to hold the factor submatrices
   ---------------------------------------------------------------------
*/
bridge->mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(bridge->mtxmanager, NO_LOCK, 0) ;
/*
   ---------------------------------------------
   STEP 12: create a SolveMap object to hold the
            map from submatrices to processes
   ---------------------------------------------
*/
bridge->solvemap = SolveMap_new() ;
/*
   -----------------------------------
   STEP 13: set up the distributed mvm
   -----------------------------------
*/
bridge->coordFlag = GLOBAL ;
if ( *pprbtype == 3 ) {
/*
   ----------------------------------
   matrix is the identity, do nothing
   ----------------------------------
*/
   bridge->info = NULL ;
   bridge->Xloc = NULL ;
   bridge->Yloc = NULL ;
} else {
   InpMtx   *mtx ;
   int      n ;
   int      *list ;
/*
   ----------------------------------------------
   generalized eigenvalue problem, find out which
   matrix to use with the matrix-vector multiply
   ----------------------------------------------
*/
   if ( *pprbtype == 1 ) {
      mtx = B ;
   } else if ( *pprbtype == 2 ) {
      mtx = A ;
   } 

   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n before MatMul_MPI_setup \n" );
      fflush(msgFile) ;
   }

   bridge->info = MatMul_MPI_setup(mtx, SPOOLES_SYMMETRIC, MMM_WITH_A, 
                                   bridge->vtxmapIV, bridge->vtxmapIV,
                                   stats, msglvl, msgFile, tag, comm) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n after MatMul_MPI_setup \n" );
      fflush(msgFile) ;
   }

/*
   ---------------------------------
   set up the Xloc and Yloc matrices
   ---------------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n SPOOLES_REAL = %d", SPOOLES_REAL) ;
      fprintf(msgFile, "\n bridge->mxbsz = %d", bridge->mxbsz) ;
      fprintf(msgFile, "\n initializing Xloc") ;
      fflush(msgFile) ;
   }
   bridge->Xloc = DenseMtx_new() ;
   DenseMtx_init(bridge->Xloc, SPOOLES_REAL, 0, 0, 
                 nmyowned, bridge->mxbsz, 1, nmyowned) ;
   DenseMtx_rowIndices(bridge->Xloc, &n, &list) ;
   IVcopy(n, list, IV_entries(bridge->myownedIV)) ;
   bridge->Yloc = DenseMtx_new() ;
   DenseMtx_zero(bridge->Xloc) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n initializing Yloc") ;
      fflush(msgFile) ;
   }
   DenseMtx_init(bridge->Yloc, SPOOLES_REAL, 0, 0, 
                 nmyowned, bridge->mxbsz, 1, nmyowned) ;
   DenseMtx_rowIndices(bridge->Yloc, &n, &list) ;
   IVcopy(n, list, IV_entries(bridge->myownedIV)) ;
   DenseMtx_zero(bridge->Yloc) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n initial Xloc matrix") ;
      fprintf(msgFile, "\n\n initial Yloc matrix") ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n leaving SetupMPI()\n") ;
   fflush(msgFile) ;
}

#if MYDEBUG > 0
MARKTIME(t2) ;
time_Setup += t2 - t1 ;
if ( bridge->myid == 0 ) {
   fprintf(stdout, ", %8.3f seconds, %8.3f total time",
           t2 - t1, time_Setup) ;
   fflush(stdout) ;
}
#endif
#if MYDEBUG > 1
fprintf(bridge->msgFile, ", %8.3f seconds, %8.3f total time",
        t2 - t1, time_Setup) ;
fflush(bridge->msgFile) ;
#endif

return(1) ; }

/*--------------------------------------------------------------------*/
