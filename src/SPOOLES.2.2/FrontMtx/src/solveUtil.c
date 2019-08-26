/*  solveUtil.c  */

#include "../FrontMtx.h"
#include "../../SubMtxList.h"

/*--------------------------------------------------------------------*/
static SubMtx * initBJ ( int type, int J, int nJ, int nrhs, 
   SubMtxManager *mtxmanager, int msglvl, FILE *msgFile ) ;
static void computeForwardUpdates ( FrontMtx *frontmtx, SubMtx *BJ,
   int J, IP *heads[], char frontIsDone[], SubMtx *p_mtx[],
   int msglvl, FILE *msgFile ) ;
static void assembleAggregates ( int J, SubMtx *BJ, SubMtxList *aggList,
   SubMtxManager *mtxmanager, int msglvl, FILE *msgFile ) ;
static void computeBackwardUpdates ( FrontMtx *frontmtx, SubMtx *ZJ,
   int J, IP *heads[], char frontIsDone[], SubMtx *p_mtx[],
   int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
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
) {
char     localrhs ;
SubMtx   *BJ ;
SubMtx   **p_agg ;
double   *bJ, *rhs ;
int      inc1, inc2, irow, jrhs, J, kk, nbytes, ncolJ, 
         neqns, nfront, nJ, nrhs, nrowInRhs, nrowJ ;
int      *rowind, *rowindJ, *rowmap ;

nrowInRhs = rhsmtx->nrow ;
neqns     = frontmtx->neqns ;
if ( nrowInRhs != neqns ) {
/*
   ----------------------------------------------------
   the rhs matrix is only part of the total rhs matrix.
   (this happens in an MPI environment where the rhs
   is partitioned among the processors.)
   create a map from the global row indices to the
   indices local to this rhs matrix.
   ----------------------------------------------------
*/
   rowmap = IVinit(neqns, -1) ;
   rowind = rhsmtx->rowind ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n rhsmtx->rowind") ;
      IVfprintf(msgFile, rhsmtx->nrow, rowind) ;
      fflush(msgFile) ;
   }
   for ( irow = 0 ; irow < nrowInRhs ; irow++ ) {
      rowmap[rowind[irow]] = irow ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n rowmap") ;
      IVfprintf(msgFile, neqns, rowmap) ;
      fflush(msgFile) ;
   }
   localrhs = 'T' ;
} else {
   localrhs = 'F' ;
}
/*
   --------------------------------------
   allocate the vector of SubMtx pointers
   --------------------------------------
*/ 
nfront = FrontMtx_nfront(frontmtx) ;
ALLOCATE(p_agg, struct _SubMtx *, nfront) ;
for ( J = 0 ; J < nfront ; J++ ) {
   p_agg[J] = NULL ;
}
DenseMtx_dimensions(rhsmtx, &neqns, &nrhs) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if (  (owners == NULL || owners[J] == myid)
      && (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
      if ( localrhs == 'T' ) {
/*
        ------------------------------------------------------
         map the global row indices into the local row indices
        ------------------------------------------------------
*/
         for ( irow = 0 ; irow < nJ ; irow++ ) {
            rowindJ[irow] = rowmap[rowindJ[irow]] ;
         }
      }
/*
      ------------------------------------------------
      initialize bmtx, a SubMtx object to hold b_{J,*}
      ------------------------------------------------
*/
      nbytes = SubMtx_nbytesNeeded(frontmtx->type, SUBMTX_DENSE_COLUMNS,
                                   nJ, nrhs, nJ*nrhs);
      BJ = SubMtxManager_newObjectOfSizeNbytes(mtxmanager, nbytes) ;
      SubMtx_init(BJ, frontmtx->type, SUBMTX_DENSE_COLUMNS, 
                  J, 0, nJ, nrhs, nJ*nrhs) ;
      p_agg[J] = BJ ;
      if ( BJ == NULL ) {
         fprintf(stderr,
            "\n fatal error in load rhs(%d), BJ = NULL", J) ;
         exit(-1) ;
      }
/*
      -----------------------
      load BJ with b_{J,*}
      -----------------------
*/
      rhs = DenseMtx_entries(rhsmtx) ;
      SubMtx_denseInfo(BJ, &nrowJ, &ncolJ, &inc1, &inc2, &bJ) ;
      if ( FRONTMTX_IS_REAL(frontmtx) ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( irow = 0 ; irow < nJ ; irow++ ) {
               kk = rowindJ[irow] ;
               bJ[irow] = rhs[kk] ;
            }
            bJ  += nJ ;
            rhs += nrowInRhs ;
         }
      } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( irow = 0 ; irow < nJ ; irow++ ) {
               kk = rowindJ[irow] ;
               bJ[2*irow]   = rhs[2*kk]   ;
               bJ[2*irow+1] = rhs[2*kk+1] ;
            }
            bJ  += 2*nJ ;
            rhs += 2*nrowInRhs ;
         }
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n rhs for J = %d, BJ = %p", J, BJ) ;
         fflush(msgFile) ;
         SubMtx_writeForHumanEye(BJ, msgFile) ;
         fflush(msgFile) ;
      }
      if ( localrhs == 'T' ) {
/*
        -----------------------------------------------------------
         map the local row indices back into the global row indices
        -----------------------------------------------------------
*/
         for ( irow = 0 ; irow < nJ ; irow++ ) {
            rowindJ[irow] = rowind[rowindJ[irow]] ;
         }
      }
   }
}
if ( localrhs == 'T' ) {
   IVfree(rowmap) ;
}
return(p_agg) ; }

/*--------------------------------------------------------------------*/
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
) {
char     aggDone, updDone ;
SubMtx   *BJ, *LJJ, *UJJ ;
int      nJ ;
 
if ( (nJ = FrontMtx_frontSize(frontmtx, J)) == 0 ) {
/*
   -----------------------------------------------------
   front has no eliminated rows or columns, quick return
   -----------------------------------------------------
*/
   if ( owners == NULL || owners[J] == myid ) {
      frontIsDone[J] = 'Y'  ;
   }
   status[J] = 'F' ;
   return ;
}
if ( heads[J] != NULL ) {
/*
   -------------------------------------
   there are internal updates to perform
   -------------------------------------
*/
   if ( (BJ = p_agg[J]) == NULL ) {
/*
      ---------------------------
      create the aggregate object
      ---------------------------
*/
      BJ = p_agg[J] = initBJ(frontmtx->type, J, nJ, nrhs, 
                             mtxmanager, msglvl, msgFile) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n BJ = %p", BJ) ;
      fflush(msgFile) ;
      SubMtx_writeForHumanEye(BJ, msgFile) ;
      fflush(msgFile) ;
   }
/*
   ---------------------------
   compute any waiting updates
   ---------------------------
*/
   computeForwardUpdates(frontmtx, BJ, J, heads, frontIsDone, p_mtx,
                         msglvl, msgFile) ;
}
if ( heads[J] == NULL ) {
   updDone = 'Y' ;
} else {
   updDone = 'N' ;
}
if ( aggList != NULL && owners[J] == myid ) {
/*
   -----------------------
   assemble any aggregates
   -----------------------
*/
   aggDone = 'N' ;
   if ( (BJ = p_agg[J]) == NULL ) {
      fprintf(stderr,
             "\n 2. fatal error in forwardVisit(%d), BJ = NULL", J) ;
      exit(-1) ;
   }
   assembleAggregates(J, BJ, aggList, mtxmanager, 
                             msglvl, msgFile) ;
   if ( SubMtxList_isCountZero(aggList, J) == 1 ) {
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n aggregate count is zero") ;
         fflush(msgFile) ;
      }
      assembleAggregates(J, BJ, aggList, mtxmanager, 
                                msglvl, msgFile) ;
      aggDone = 'Y' ;
   }
} else {
   aggDone = 'Y' ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n updDone = %c, aggDone = %c", 
           updDone, aggDone) ;
   fflush(msgFile) ;
}
if ( updDone == 'Y' && aggDone == 'Y' ) {
   BJ = p_agg[J] ;
   if ( owners == NULL || owners[J] == myid ) {
/*
      -------------------------------------
      owned front, ready for interior solve
      -------------------------------------
*/
      if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
         LJJ = FrontMtx_lowerMtx(frontmtx, J, J) ;
         if ( LJJ != NULL ) {
            SubMtx_solve(LJJ, BJ) ;
         }
      } else {
         UJJ = FrontMtx_upperMtx(frontmtx, J, J) ;
         if ( UJJ != NULL ) {
            if ( FRONTMTX_IS_SYMMETRIC(frontmtx) ) {
               SubMtx_solveT(UJJ, BJ) ;
            } else if ( FRONTMTX_IS_HERMITIAN(frontmtx) ) {
               SubMtx_solveH(UJJ, BJ) ;
            }
         }
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n after forward solve") ;
         SubMtx_writeForHumanEye(BJ, msgFile) ;
         fflush(msgFile) ;
      }
/*
      ------------------------------------------------
      move YJ (stored in BJ) into p_mtx[],
      signal front as done, and set status to finished
      ------------------------------------------------
*/
      p_agg[J]       = NULL ;
      p_mtx[J]       = BJ   ;
      frontIsDone[J] = 'Y'  ;
   } else if ( BJ != NULL ) {
/*
      --------------------------------------
      unowned front, put into aggregate list
      --------------------------------------
*/
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n putting BJ into aggregate list") ;
         fflush(msgFile) ;
      }
      SubMtxList_addObjectToList(aggList, BJ, J) ;
      p_agg[J]  = NULL ;
   }
   status[J] = 'F' ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
if ( owners == NULL || owners[J] == myid ) {
   SubMtx   *BJ, *DJJ ;
 
   if ( (BJ = p_mtx[J]) != NULL ) {
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n BJ = %p", BJ) ;
         SubMtx_writeForHumanEye(BJ, msgFile) ;
         fflush(msgFile) ;
      }
      DJJ = FrontMtx_diagMtx(frontmtx, J) ;
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n DJJ = %p", DJJ) ;
         SubMtx_writeForHumanEye(DJJ, msgFile) ;
         fflush(msgFile) ;
      }
      SubMtx_solve(DJJ, BJ) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n b_{%d,*} after diagonal solve", J) ;
         SubMtx_writeForHumanEye(BJ, msgFile) ;
         fflush(msgFile) ;
      }
      p_mtx[J] = NULL ;
      p_agg[J] = BJ   ;
   }
   frontIsDone[J] = 'Y' ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
char     aggDone, updDone ;
SubMtx     *UJJ, *ZJ ;
int      nJ ;
 
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n inside FrontMtx_backwardVisit(%d), nJ = %d",
           J, FrontMtx_frontSize(frontmtx, J)) ;
   fflush(msgFile) ;
}
if ( (nJ = FrontMtx_frontSize(frontmtx, J)) == 0 ) {
/*
   -----------------------------------------------------
   front has no eliminated rows or columns, quick return
   -----------------------------------------------------
*/
   if ( owners == NULL || owners[J] == myid ) {
      frontIsDone[J] = 'Y'  ;
   }
   status[J] = 'F' ;
   return ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n heads[%d] = %p", J, heads[J]) ;
   fflush(msgFile) ;
}
if ( heads[J] != NULL ) {
/*
   -------------------------------------
   there are internal updates to perform
   -------------------------------------
*/
   if ( (ZJ = p_agg[J]) == NULL ) {
/*
      ---------------------------
      create the aggregate object
      ---------------------------
*/
      ZJ = p_agg[J] = initBJ(frontmtx->type, J, nJ, nrhs, 
                             mtxmanager, msglvl, msgFile) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n ZJ = %p", ZJ) ;
      SubMtx_writeForHumanEye(ZJ, msgFile) ;
      fflush(msgFile) ;
   }
/*
   ---------------------------
   compute any waiting updates
   ---------------------------
*/
   computeBackwardUpdates(frontmtx, ZJ, J, heads, frontIsDone, p_mtx,
                          msglvl, msgFile) ;
}
if ( heads[J] == NULL ) {
   updDone = 'Y' ;
} else {
   updDone = 'N' ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n updDone = %c", updDone) ;
   fflush(msgFile) ;
}
if ( aggList != NULL && owners[J] == myid ) {
/*
   -----------------------
   assemble any aggregates
   -----------------------
*/
   aggDone = 'N' ;
   if ( (ZJ = p_agg[J]) == NULL ) {
      fprintf(stderr,
             "\n 2. fatal error in backwardVisit(%d), ZJ = NULL", J) ;
      exit(-1) ;
   }
   assembleAggregates(J, ZJ, aggList, mtxmanager, 
                             msglvl, msgFile) ;
   if ( SubMtxList_isCountZero(aggList, J) == 1 ) {
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n aggregate count is zero") ;
         fflush(msgFile) ;
      }
      assembleAggregates(J, ZJ, aggList, mtxmanager, 
                                msglvl, msgFile) ;
      aggDone = 'Y' ;
   }
} else {
   aggDone = 'Y' ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n aggDone = %c", aggDone) ;
   fflush(msgFile) ;
}
if ( updDone == 'Y' && aggDone == 'Y' ) {
   ZJ = p_agg[J] ;
   if ( owners == NULL || owners[J] == myid ) {
/*
      -------------------------------------
      owned front, ready for interior solve
      -------------------------------------
*/
      UJJ = FrontMtx_upperMtx(frontmtx, J, J) ;
      if ( UJJ != NULL ) {
         SubMtx_solve(UJJ, ZJ) ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n after backward solve") ;
         SubMtx_writeForHumanEye(ZJ, msgFile) ;
         fflush(msgFile) ;
      }
/*
      ------------------------------------------------
      move YJ (stored in BJ) into p_mtx[],
      signal front as done, and set status to finished
      ------------------------------------------------
*/
      p_agg[J]       = NULL ;
      p_mtx[J]       = ZJ   ;
      frontIsDone[J] = 'Y'  ;
   } else if ( ZJ != NULL ) {
/*
      --------------------------------------
      unowned front, put into aggregate list
      --------------------------------------
*/
      SubMtxList_addObjectToList(aggList, ZJ, J) ;
      p_agg[J]  = NULL ;
   }
   status[J] = 'F'  ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n status[%d] = %c", J, status[J]) ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
char     localsol ;
SubMtx   *xmtxJ ;
double   *sol, *xJ ;
int      inc1, inc2, irow, jrhs, J, kk,
         ncolJ, neqns, nfront, nJ, nrhs, nrowInSol, nrowJ ;
int      *colindJ, *colmap, *rowind ;

if ( (nrowInSol = solmtx->nrow) != (neqns = frontmtx->neqns) ) {
/*
   --------------------------------------------------------------
   the solution matrix is only part of the total solution matrix.
   (this happens in an MPI environment where the rhs
   is partitioned among the processors.)
   create a map from the global row indices to the
   indices local to this solution matrix.
   --------------------------------------------------------------
*/
   colmap = IVinit(neqns, -1) ;
   rowind = solmtx->rowind ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n solmtx->rowind") ;
      IVfprintf(msgFile, solmtx->nrow, rowind) ;
      fflush(msgFile) ;
   }
   for ( irow = 0 ; irow < nrowInSol ; irow++ ) {
      colmap[rowind[irow]] = irow ;
   }
   localsol = 'T' ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n colmap") ;
      IVfprintf(msgFile, neqns, colmap) ;
      fflush(msgFile) ;
   }
} else {
   localsol = 'F' ;
}
DenseMtx_dimensions(solmtx, &neqns, &nrhs) ;
nfront = FrontMtx_nfront(frontmtx) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if (  (owners == NULL || owners[J] == myid)
      && (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
      xmtxJ = p_mtx[J] ;
      if ( xmtxJ == NULL ) {
         fprintf(stderr,
            "\n fatal error in storeSolution(%d)"
            "\n thread %d, xmtxJ = NULL", J, myid) ;
         exit(-1) ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n storing solution for front %d", J) ;
         SubMtx_writeForHumanEye(xmtxJ, msgFile) ;
         fflush(msgFile) ;
      }
      if ( localsol == 'T' ) {
/*
        ------------------------------------------------------
         map the global row indices into the local row indices
        ------------------------------------------------------
*/
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n global row indices") ;
            IVfprintf(msgFile, nJ, colindJ) ;
            fflush(msgFile) ;
         }
         for ( irow = 0 ; irow < nJ ; irow++ ) {
            colindJ[irow] = colmap[colindJ[irow]] ;
         }
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n local row indices") ;
            IVfprintf(msgFile, nJ, colindJ) ;
            fflush(msgFile) ;
         }
      }
/*
      ----------------------------------
      store x_{J,*} into solution matrix
      ----------------------------------
*/
      sol = DenseMtx_entries(solmtx) ;
      SubMtx_denseInfo(xmtxJ, &nrowJ, &ncolJ, &inc1, &inc2, &xJ) ;
      if ( FRONTMTX_IS_REAL(frontmtx) ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( irow = 0 ; irow < nJ ; irow++ ) {
               kk = colindJ[irow] ;
               sol[kk] = xJ[irow] ;
            }
            sol += neqns ;
            xJ  += nJ ;
         }
      } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( irow = 0 ; irow < nJ ; irow++ ) {
               kk = colindJ[irow] ;
               sol[2*kk]   = xJ[2*irow]   ;
               sol[2*kk+1] = xJ[2*irow+1] ;
            }
            sol += 2*neqns ;
            xJ  += 2*nJ ;
         }
      }
/*
fprintf(msgFile, "\n solution for front %d stored", J) ;
*/
      SubMtxManager_releaseObject(manager, xmtxJ) ;
      if ( localsol == 'T' ) {
/*
        -----------------------------------------------------------
         map the local row indices back into the global row indices
        -----------------------------------------------------------
*/
         for ( irow = 0 ; irow < nJ ; irow++ ) {
            colindJ[irow] = rowind[colindJ[irow]] ;
         }
      }
   }
}
if ( localsol == 'T' ) {
   IVfree(colmap) ;
}
/*
fprintf(msgFile, "\n\n SOLUTION") ;
DenseMtx_writeForHumanEye(solmtx, msgFile) ;
*/

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- initialize an aggregate SubMtx object

   created -- 98mar27, cca
   ----------------------------------------------
*/
static SubMtx *
initBJ (
   int             type,
   int             J,
   int             nJ,
   int             nrhs,
   SubMtxManager   *mtxmanager,
   int             msglvl,
   FILE            *msgFile
) {
SubMtx     *BJ ;
double   *entries ;
int      inc1, inc2, nbytes ;
/*
   ------------------------------------------
   B_J not yet allocated (must not be owned),
   create and zero the entries
   ------------------------------------------
*/
nbytes = SubMtx_nbytesNeeded(type, SUBMTX_DENSE_COLUMNS, 
                             nJ, nrhs, nJ*nrhs);
BJ = SubMtxManager_newObjectOfSizeNbytes(mtxmanager, nbytes) ;
if ( BJ == NULL ) {
   fprintf(stderr,
          "\n 1. fatal error in forwardVisit(%d), BJ = NULL", J) ;
   exit(-1) ;
}
SubMtx_init(BJ, type, SUBMTX_DENSE_COLUMNS, J, 0, nJ, nrhs, nJ*nrhs) ;
SubMtx_denseInfo(BJ, &nJ, &nrhs, &inc1, &inc2, &entries) ;
if ( type == SPOOLES_REAL ) {
   DVzero(nJ*nrhs, entries) ;
} else if ( type == SPOOLES_COMPLEX ) {
   DVzero(2*nJ*nrhs, entries) ;
}
return(BJ) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- compute any updates to BJ

   created -- 98mar26, cca
   ------------------------------------
*/
static void
computeForwardUpdates (
   FrontMtx   *frontmtx,
   SubMtx     *BJ,
   int        J,
   IP         *heads[],
   char       frontIsDone[],
   SubMtx     *p_mtx[],
   int        msglvl,
   FILE       *msgFile
) {
SubMtx   *LJI, *UIJ, *YI ;
int    I ;
IP     *ip, *nextip ;
/*
   -------------------------------
   loop over the remaining updates
   -------------------------------
*/
for ( ip = heads[J], heads[J] = NULL ;
      ip != NULL ;
      ip = nextip ) {
   I = ip->val ; nextip = ip->next ;
   if ( msglvl > 3 ) {
      fprintf(msgFile,
              "\n\n frontIsDone[%d] = %c, p_mtx[%d] = %p", 
              I, frontIsDone[I], I, p_mtx[I]) ;
      fflush(msgFile) ;
   }
   if ( frontIsDone[I] == 'Y' ) {
      if ( (YI = p_mtx[I]) != NULL ) {
/*
         --------------------------------
         Y_I exists and has been computed
         --------------------------------
*/
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n\n before solve: YI = %p", YI) ;
            SubMtx_writeForHumanEye(YI, msgFile) ;
            fflush(msgFile) ;
         }
         if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
            if ( (LJI = FrontMtx_lowerMtx(frontmtx, J, I)) != NULL ) {
               if ( msglvl > 3 ) {
                  fprintf(msgFile, "\n\n LJI = %p", LJI) ;
                  SubMtx_writeForHumanEye(LJI, msgFile) ;
                  fflush(msgFile) ;
               }
               SubMtx_solveupd(BJ, LJI, YI) ;
            }
         } else {
            if ( (UIJ = FrontMtx_upperMtx(frontmtx, I, J)) != NULL ) {
               if ( msglvl > 3 ) {
                  fprintf(msgFile, "\n\n UIJ = %p", UIJ) ;
                  SubMtx_writeForHumanEye(UIJ, msgFile) ;
                  fflush(msgFile) ;
               }
               if ( FRONTMTX_IS_SYMMETRIC(frontmtx) ) {
                  SubMtx_solveupdT(BJ, UIJ, YI) ;
               } else if ( FRONTMTX_IS_HERMITIAN(frontmtx) ) {
                  SubMtx_solveupdH(BJ, UIJ, YI) ;
               }
            }
         }
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n\n after update: BJ = %p", BJ) ;
            SubMtx_writeForHumanEye(BJ, msgFile) ;
            fflush(msgFile) ;
         }
      }
   } else {
/*
      ------------------------
      Y_I is not yet available
      ------------------------
*/
      ip->next = heads[J] ;
      heads[J] = ip ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- assemble any aggregates in the aggregate list

   created -- 98mar26, cca
   --------------------------------------------------------
*/
static void
assembleAggregates (
   int             J,
   SubMtx          *BJ,
   SubMtxList      *aggList,
   SubMtxManager   *mtxmanager,
   int             msglvl, 
   FILE            *msgFile
) {
SubMtx     *BJhat, *BJhead ;
double   *entBJ, *entBJhat ;
int      inc1, inc1hat, inc2, inc2hat, ncol, ncolhat, nrow, nrowhat ;
 
if ( BJ == NULL || aggList == NULL ) {
   fprintf(stderr,
          "\n fatal error in assembleAggregates()"
          "\n BJ = %p, aggList = %p", BJ, aggList) ;
   exit(-1) ;
}
if ( SubMtxList_isListNonempty(aggList, BJ->rowid) ) {
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n aggregate list is not-empty") ;
      fflush(msgFile) ;
   }
   SubMtx_denseInfo(BJ, &nrow, &ncol, &inc1, &inc2, &entBJ) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile,
          "\n\n BJ(%d,%d) : nrow %d, ncol %d, inc1 %d, inc2 %d, ent %p",
          BJ->rowid, BJ->colid, nrow, ncol, inc1, inc2, entBJ) ;
      fflush(msgFile) ;
   }
   BJhead = SubMtxList_getList(aggList, J) ;
   for ( BJhat = BJhead ; BJhat != NULL ; BJhat = BJhat->next ) {
      if ( BJhat == NULL ) {
         fprintf(stderr,
                 "\n 3. fatal error in forwardVisit(%d)"
                 "\n BJhat = NULL", J) ;
         exit(-1) ;
      }
      SubMtx_denseInfo(BJhat, &nrowhat, &ncolhat, &inc1hat, &inc2hat,
                       &entBJhat) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile,
         "\n BJhat(%d,%d) : nrow %d, ncol %d, inc1 %d, inc2 %d, ent %p",
           BJhat->rowid, BJhat->colid,
           nrowhat, ncolhat, inc1hat, inc2hat, entBJhat) ;
         fflush(msgFile) ;
      }
      if ( nrow != nrowhat || ncol != ncolhat
         || inc1 != inc1hat || inc2 != inc2hat || entBJhat == NULL ) {
         fprintf(stderr, "\n fatal error") ;
         exit(-1) ;
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n BJ") ;
         SubMtx_writeForHumanEye(BJ, msgFile) ;
         fprintf(msgFile, "\n\n BJhat") ;
         SubMtx_writeForHumanEye(BJhat, msgFile) ;
         fflush(msgFile) ;
      }
      if ( SUBMTX_IS_REAL(BJhat) ) {
         DVadd(nrow*ncol, entBJ, entBJhat) ;
      } else if ( SUBMTX_IS_COMPLEX(BJhat) ) {
         DVadd(2*nrow*ncol, entBJ, entBJhat) ;
      }
   }
   SubMtxManager_releaseListOfObjects(mtxmanager, BJhead) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n BJ after assembly") ;
      SubMtx_writeForHumanEye(BJ, msgFile) ;
      fflush(msgFile) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- compute any updates to ZJ

   created -- 98mar26, cca
   ------------------------------------
*/
static void
computeBackwardUpdates (
   FrontMtx   *frontmtx,
   SubMtx     *ZJ,
   int        J,
   IP         *heads[],
   char       frontIsDone[],
   SubMtx     *p_mtx[],
   int        msglvl,
   FILE       *msgFile
) {
SubMtx   *UJK, *XK ;
int    K ;
IP     *ip, *nextip ;
/*
   -------------------------------
   loop over the remaining updates
   -------------------------------
*/
for ( ip = heads[J], heads[J] = NULL ;
      ip != NULL ;
      ip = nextip ) {
   K = ip->val ; nextip = ip->next ;
   if ( msglvl > 3 ) {
      fprintf(msgFile,
              "\n\n frontIsDone[%d] = %c", K, frontIsDone[K]) ;
      fflush(msgFile) ;
   }
   if ( frontIsDone[K] == 'Y' ) {
      if ( (XK = p_mtx[K]) != NULL ) {
/*
         --------------------------------
         X_K exists and has been computed
         --------------------------------
*/
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n\n before solve: XK = %p", XK) ;
            SubMtx_writeForHumanEye(XK, msgFile) ;
            fflush(msgFile) ;
         }
         if ( (UJK = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            if ( msglvl > 3 ) {
               fprintf(msgFile, "\n\n UJK = %p", UJK) ;
               SubMtx_writeForHumanEye(UJK, msgFile) ;
               fflush(msgFile) ;
            }
            SubMtx_solveupd(ZJ, UJK, XK) ;
         }
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n\n after update: ZJ = %p", ZJ) ;
            SubMtx_writeForHumanEye(ZJ, msgFile) ;
            fflush(msgFile) ;
         }
      }
   } else {
/*
      ------------------------
      X_K is not yet available
      ------------------------
*/
      ip->next = heads[J] ;
      heads[J] = ip ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
) {
int    J ;
int    *sib ;
 
sib    = frontmtx->tree->sib ;
for ( J = frontmtx->tree->root ; J != -1 ; J = sib[J] ) {
   if ( status[J] == activeFlag ) {
      Ideq_insertAtTail(dequeue, J) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
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
   int        msglvl,
   FILE       *msgFile
) {
int   ii, J, K, nadj, nblock, nfront ;
int   *adj ;
IP    *ip ;
IP    **heads ;
/*
   --------------------------------------------------
   set up the update head/links for the forward solve
   --------------------------------------------------
*/
nblock = FrontMtx_nLowerBlocks(frontmtx) ;
nfront = FrontMtx_nfront(frontmtx) ;
ALLOCATE(heads, struct _IP *, nfront + 2) ;
for ( J = 0 ; J <= nfront + 1 ; J++ ) {
   heads[J] = NULL ;
}
heads[nfront] = heads[nfront+1] = IP_init(nblock, 1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   FrontMtx_lowerAdjFronts(frontmtx, J, &nadj, &adj) ;
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      if ( (K = adj[ii]) > J ) {
         ip = heads[nfront] ; heads[nfront] = ip->next ;
         ip->val = J ; ip->next = heads[K] ; heads[K] = ip ;
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n linking L(%d,%d) to L(%d,%d)",
                    K, J, K, (ip->next == NULL) ? -1 : ip->next->val) ;
            fflush(msgFile) ;
         }
      }
   }
}
return(heads) ; }
 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- set up the linked lists for the backward solve
 
   created -- 98feb20, cca
   ---------------------------------------------------------
*/
IP **
FrontMtx_backwardSetup (
   FrontMtx   *frontmtx,
   int        msglvl,
   FILE       *msgFile
) {
IP    *ip ;
IP    **heads ;
int   ii, J, K, nadj, nblock, nfront ;
int   *adj ;
 
nfront = FrontMtx_nfront(frontmtx) ;
nblock = FrontMtx_nLowerBlocks(frontmtx) ;
ALLOCATE(heads, struct _IP *, nfront + 2) ;
for ( J = 0 ; J <= nfront + 1 ; J++ ) {
   heads[J] = NULL ;
}
heads[nfront] = heads[nfront+1] = IP_init(nblock, 1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   FrontMtx_upperAdjFronts(frontmtx, J, &nadj, &adj) ;
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      if ( (K = adj[ii]) > J ) {
         if ( heads[nfront] == NULL ) {
            fprintf(stderr,
                    "\n WHOA, heads[nfront] = %p", heads[nfront]) ;
            exit(-1) ;
         }
         ip = heads[nfront] ; heads[nfront] = ip->next ;
         ip->val = K ; ip->next = heads[J] ; heads[J] = ip ;
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n linking U(%d,%d) to U(%d,%d)",
                    J, K, J, (ip->next == NULL) ? -1 : ip->next->val) ;
            fflush(msgFile) ;
         }
      }
   }
}
return(heads) ; }

/*--------------------------------------------------------------------*/
