/*  init.c  */

#include "../FrontMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
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
) {
SubMtx   *mtx ;
int      J, nbytes, nentD, nentL, nentU, neqns, nD, nent, nfront, nU ;
int      *bndwghts, *nodwghts, *owners, *vtxToFront ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || frontETree == NULL || symbfacIVL == NULL
   || (ownersIV != NULL && myid < 0) 
   || manager == NULL ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_init()"
           "\n frontmtx %p, frontETree %p, symbfacIVL %p"
           "\n myid %d, ownersIV %p, manager %p"
           "\n bad input\n", 
           frontmtx, frontETree, symbfacIVL, myid, ownersIV, manager) ;
   exit(-1) ;
}
if ( type != SPOOLES_REAL && type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in FrontMtx_init()"
           "\n type %d must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           type) ;
   exit(-1) ;
}
if ( type == SPOOLES_REAL && 
   ! (symmetryflag == SPOOLES_SYMMETRIC
      || symmetryflag == SPOOLES_NONSYMMETRIC) ) {
   fprintf(stderr, 
  "\n fatal error in FrontMtx_init()"
  "\n type is real"
  "\n symmetryflag is not SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC") ;
   exit(-1) ;
}
if ( type == SPOOLES_COMPLEX && 
   ! (symmetryflag == SPOOLES_SYMMETRIC
      || symmetryflag == SPOOLES_HERMITIAN
      || symmetryflag == SPOOLES_NONSYMMETRIC) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_init()"
           "\n type is real, symmetryflag is not SPOOLES_SYMMETRIC,"
           "\n SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC") ;
   exit(-1) ;
}
if ( ! ( pivotingflag == SPOOLES_PIVOTING 
      || pivotingflag == SPOOLES_NO_PIVOTING) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_init()"
  "\n pivotingflag must be SPOOLES_PIVOTING or SPOOLES_NO_PIVOTING\n") ;
   exit(-1) ;
}
if ( !  (lockflag == NO_LOCK
      || lockflag == LOCK_IN_PROCESS
      || lockflag == LOCK_OVER_ALL_PROCESSES) ) {
   fprintf(stderr, 
       "\n fatal error in FrontMtx_init()"
       "\n invalid lockflag, must be NO_LOCK"
       "\n LOCK_IN_PROCESS or LOCK_OVER_ALL_PROCESSES") ;
   exit(-1) ;
}
nfront     = frontETree->nfront         ;
neqns      = frontETree->nvtx           ;
nodwghts   = ETree_nodwghts(frontETree) ;
bndwghts   = ETree_bndwghts(frontETree) ;
vtxToFront = ETree_vtxToFront(frontETree) ;
if ( ownersIV != NULL ) {
   owners = IV_entries(ownersIV) ;
} else {
   owners = NULL ;
}
/*
   ----------------------
   set the default fields
   ----------------------
*/
FrontMtx_setDefaultFields(frontmtx) ;
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
frontmtx->nfront       = nfront       ;
frontmtx->neqns        = neqns        ;
frontmtx->type         = type         ;
frontmtx->symmetryflag = symmetryflag ;
frontmtx->sparsityflag = sparsityflag ;
frontmtx->pivotingflag = pivotingflag ;
frontmtx->dataMode     = FRONTMTX_1D_MODE ;
/*
   ---------------------------------------------------------------
   set the front tree ETree and symbolic factorization IVL objects
   ---------------------------------------------------------------
*/
frontmtx->tree       = frontETree->tree ;
frontmtx->frontETree = frontETree ;
frontmtx->symbfacIVL = symbfacIVL ;
/*
   ----------------------------
   set the frontsizes IV object
   ----------------------------
*/
frontmtx->frontsizesIV = IV_new() ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   IV_init(frontmtx->frontsizesIV, nfront, NULL) ;
   IV_fill(frontmtx->frontsizesIV, 0) ;
} else {
   IV_init(frontmtx->frontsizesIV, nfront, nodwghts) ;
}
/*
   ------------------------------------------------------
   set the row and column adjacency objects, if necessary
   ------------------------------------------------------
*/
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   frontmtx->coladjIVL = IVL_new() ;
   IVL_init1(frontmtx->coladjIVL, IVL_CHUNKED, nfront) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      frontmtx->rowadjIVL = IVL_new() ;
      IVL_init1(frontmtx->rowadjIVL, IVL_CHUNKED, nfront) ;
   }
}
/*
   ---------------------------------
   allocate the five pointer vectors
   ---------------------------------
*/
ALLOCATE(frontmtx->p_mtxDJJ, struct _SubMtx *, nfront) ;
ALLOCATE(frontmtx->p_mtxUJJ, struct _SubMtx *, nfront) ;
ALLOCATE(frontmtx->p_mtxUJN, struct _SubMtx *, nfront) ;
for ( J = 0 ; J < nfront ; J++ ) {
   frontmtx->p_mtxDJJ[J] = NULL ;
   frontmtx->p_mtxUJJ[J] = NULL ;
   frontmtx->p_mtxUJN[J] = NULL ;
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   ALLOCATE(frontmtx->p_mtxLJJ, struct _SubMtx *, nfront) ;
   ALLOCATE(frontmtx->p_mtxLNJ, struct _SubMtx *, nfront) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      frontmtx->p_mtxLJJ[J] = NULL ;
      frontmtx->p_mtxLNJ[J] = NULL ;
   }
}
/*
   -------------------------------
   set the submatrix manager field
   -------------------------------
*/
frontmtx->manager = manager ;
/*
   -------------------------------------
   initialize submatrices where possible
   -------------------------------------
*/
if (  ! FRONTMTX_IS_PIVOTING(frontmtx) 
   && FRONTMTX_IS_DENSE_FRONTS(frontmtx) ) {
   double   *entries ;
   int      ii, jj, ncol, nrow ;
   int      *firstlocs, *sizes ;

   nentD = nentL = nentU = 0 ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( owners == NULL || owners[J] == myid ) {
         nD = nodwghts[J] ;
         nU = bndwghts[J] ;
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n\n J %d, nD %d, nU %d", J, nD, nU) ;
            fflush(msgFile) ;
         }
/*
         ---------------
         diagonal matrix
         ---------------
*/
         nbytes = SubMtx_nbytesNeeded(type, SUBMTX_DIAGONAL, 
                                      nD, nD, nD) ;
         mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
         SubMtx_init(mtx, type, SUBMTX_DIAGONAL, J, J, nD, nD, nD);
         SubMtx_diagonalInfo(mtx, &ncol, &entries) ;
         SubMtx_zero(mtx) ;
         nentD += nD ;
         frontmtx->p_mtxDJJ[J] = mtx ;
         if ( msglvl > 3 ) {
            fprintf(msgFile, 
                  "\n diagonal (%d,%d) matrix %p, %d entries, %d bytes",
                  J, J, mtx, nD, nbytes) ;
            fflush(msgFile) ;
         }
         if ( (nent = (nD*(nD-1))/2) > 0 ) {
/*
            -------------------------------------------
            U_{J,J} and possibly lower L_{J,J} matrices
            -------------------------------------------
*/
            nbytes = SubMtx_nbytesNeeded(type, SUBMTX_DENSE_SUBCOLUMNS,
                                         nD, nD, nent) ;
            mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
            SubMtx_init(mtx, type, SUBMTX_DENSE_SUBCOLUMNS, 
                        J, J, nD, nD, nent);
            SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, &firstlocs,
                                       &sizes, &entries) ;
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               firstlocs[jj] = 0 ;
               sizes[jj]     = jj ;
            }
            SubMtx_zero(mtx) ;
            nentU += nent ;
            frontmtx->p_mtxUJJ[J] = mtx ;
            if ( msglvl > 3 ) {
               fprintf(msgFile, 
                "\n upper (%d,%d) matrix %p, %d entries, %d bytes",
                J, J, mtx, nent, nbytes) ;
               fflush(msgFile) ;
            }
            if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
               nbytes = SubMtx_nbytesNeeded(type, SUBMTX_DENSE_SUBROWS,
                                            nD, nD, nent) ;
               mtx = SubMtxManager_newObjectOfSizeNbytes(manager, 
                                                         nbytes);
               SubMtx_init(mtx, type, SUBMTX_DENSE_SUBROWS, 
                           J, J, nD, nD, nent);
               SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, &firstlocs,
                                     &sizes, &entries) ;
               for ( ii = 0 ; ii < nrow ; ii++ ) {
                  firstlocs[ii] = 0 ;
                  sizes[ii]     = ii ;
               }
               SubMtx_zero(mtx) ;
               nentL += nent ;
               frontmtx->p_mtxLJJ[J] = mtx ;
               if ( msglvl > 3 ) {
                  fprintf(msgFile, 
                     "\n lower (%d,%d) matrix %p, %d entries, %d bytes",
                     J, J, mtx, nent, nbytes) ;
                  fflush(msgFile) ;
               }
            }
         }
         if ( (nent = nD*nU) > 0 ) {
/*
            -----------------------------------------------
            U_{J,bnd{J}} and possibly L_{bnd{J},J} matrices
            -----------------------------------------------
*/
            nbytes = SubMtx_nbytesNeeded(type, SUBMTX_DENSE_COLUMNS,
                                         nD, nU, nent) ;
            mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
            SubMtx_init(mtx, type, SUBMTX_DENSE_COLUMNS, 
                        J, nfront, nD, nU, nent);
            SubMtx_zero(mtx) ;
            nentU += nent ;
            frontmtx->p_mtxUJN[J] = mtx ;
            if ( msglvl > 3 ) {
               fprintf(msgFile, 
                     "\n upper (%d,%d) matrix %p, %d entries, %d bytes",
                     J, nfront, mtx, nent, nbytes) ;
               fflush(msgFile) ;
            }
            if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
               nbytes = SubMtx_nbytesNeeded(type, SUBMTX_DENSE_ROWS, 
                                            nU, nD, nent);
               mtx = SubMtxManager_newObjectOfSizeNbytes(manager, 
                                                         nbytes);
               SubMtx_init(mtx, type, SUBMTX_DENSE_ROWS, 
                           nfront, J, nU, nD, nent);
               SubMtx_zero(mtx) ;
               nentL += nent ;
               frontmtx->p_mtxLNJ[J] = mtx ;
               if ( msglvl > 3 ) {
                  fprintf(msgFile, 
                     "\n lower (%d,%d) matrix %p, %d entries, %d bytes",
                     nfront, J, mtx, nent, nbytes) ;
                  fflush(msgFile) ;
               }
            }
         }
      }
   }
   frontmtx->nentD = nentD ;
   frontmtx->nentL = nentL ;
   frontmtx->nentU = nentU ;
}
if (  lockflag == LOCK_OVER_ALL_PROCESSES 
   || lockflag == LOCK_IN_PROCESS ) {
/*
   -----------------
   allocate the lock
   -----------------
*/
   frontmtx->lock = Lock_new() ;
   Lock_init(frontmtx->lock, lockflag) ;
}
/*
   --------------------------------------
   set the patch-and-go information field
   --------------------------------------
*/
frontmtx->patchinfo = NULL ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n frontmtx->lock = %p", frontmtx->lock) ;
   fflush(msgFile) ;
}
   
return ; }

/*--------------------------------------------------------------------*/
