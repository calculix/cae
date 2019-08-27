/*  factor.c  */

#include "../BridgeMPI.h"

/*--------------------------------------------------------------------*/
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
   BridgeMPI   *bridge,
   InpMtx      *mtxA,
   int         permuteflag,
   int         *perror
) {
Chv             *rootchv ;
ChvManager      *chvmanager ;
double          cputotal, nfops, t0, t1, t2 ;
double          cpus[20] ;
int             firsttag, msglvl, myid, neqns, nproc, nzf ;
int             stats[20] ;
InpMtx          *Aloc ;
FILE            *msgFile ;
FrontMtx        *frontmtx ;
SubMtxManager   *mtxmanager ;

/*--------------------------------------------------------------------*/

MARKTIME(t0) ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_factor()"
           "\n bridge is NULL\n") ;
   return(-1) ;
}
if ( perror == NULL ) {
   fprintf(stderr, "\n error in BridgeMPI_factor()"
           "\n perror is NULL\n") ;
   return(-2) ;
}
msglvl  = bridge->msglvl  ;
msgFile = bridge->msgFile ;
neqns   = bridge->neqns   ;
myid    = bridge->myid    ;
nproc   = bridge->nproc   ;

if ( msglvl > 1 ) {
   fprintf(msgFile, "\n all set in BridgeMPI_factor()") ;
   fprintf(msgFile, "\n myid = %d", myid) ;
   fprintf(msgFile, "\n permuteflag = %d", permuteflag) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/

MARKTIME(t1) ;
if ( myid == 0 ) {
   if ( permuteflag == 1 ) {
      int   *oldToNew = IV_entries(bridge->oldToNewIV) ;
/*
      ------------------------------------------------
      permute the input matrix and convert to chevrons
      ------------------------------------------------
*/
      InpMtx_permute(mtxA, oldToNew, oldToNew) ;
      if (  bridge->symmetryflag == SPOOLES_SYMMETRIC
         || bridge->symmetryflag == SPOOLES_HERMITIAN ) {
         InpMtx_mapToUpperTriangle(mtxA) ;
      }
   }
   if ( ! INPMTX_IS_BY_CHEVRONS(mtxA) ) {
      InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(mtxA) ) {
      InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
   }
   MARKTIME(t2) ;
   bridge->cpus[8] += t2 - t1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : permute and format A", t2 - t1) ;
      fflush(msgFile) ;
   }
}
/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   distribute the front matrix.
   processor 0 scatters the entries 
   to all the other processors.
   --------------------------------
*/
MARKTIME(t1) ;
if ( myid != 0 ) {
   mtxA = NULL ;
}
firsttag = 0 ;
IVzero(4, stats) ;
Aloc = bridge->Aloc = InpMtx_MPI_splitFromGlobal(mtxA, NULL, 
                                    bridge->vtxmapIV, 0, stats, msglvl, 
                                    msgFile, firsttag, bridge->comm) ;
if ( Aloc != NULL ) {
   InpMtx_changeStorageMode(Aloc, INPMTX_BY_VECTORS) ;
}
MARKTIME(t2) ;
bridge->cpus[9] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : distribute A", t2 - t1) ;
   fflush(msgFile) ;
}
if ( Aloc != NULL && msglvl > 1 ) {
   fprintf(msgFile, "\n\n permuted and split InpMtx") ;
   InpMtx_writeForHumanEye(Aloc, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ---------------------------
   initialize the front matrix
   ---------------------------
*/
MARKTIME(t1) ;
if ( (mtxmanager = bridge->mtxmanager) == NULL ) {
   mtxmanager = bridge->mtxmanager = SubMtxManager_new() ;
   SubMtxManager_init(mtxmanager, NO_LOCK, 0) ;
}
if ( (frontmtx = bridge->frontmtx) == NULL ) {
   frontmtx = bridge->frontmtx = FrontMtx_new() ;
} else {
   FrontMtx_clearData(frontmtx) ;
}
FrontMtx_init(frontmtx, bridge->frontETree, bridge->symbfacIVL,
              bridge->type, bridge->symmetryflag, bridge->sparsityflag,
              bridge->pivotingflag, NO_LOCK, myid, bridge->ownersIV, 
              mtxmanager, bridge->msglvl, bridge->msgFile) ;
frontmtx->patchinfo = bridge->patchinfo ;
MARKTIME(t2) ;
bridge->cpus[10] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : initialize front matrix", t2 - t1) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n bridge ownersIV") ;
   IV_writeForHumanEye(bridge->ownersIV, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -----------------
   factor the matrix
   -----------------
*/
nzf   = ETree_nFactorEntries(bridge->frontETree, bridge->symmetryflag) ;
nfops = ETree_nFactorOps(bridge->frontETree, 
                         bridge->type, bridge->symmetryflag) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n %d factor entries, %.0f factor ops, %8.3f ratio",
           nzf, nfops, nfops/nzf) ;
   fflush(msgFile) ;
}
IVzero(20, stats) ;
DVzero(20, cpus) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, NO_LOCK, 1) ;
MARKTIME(t1) ;
rootchv = FrontMtx_MPI_factorInpMtx(frontmtx, Aloc, bridge->tau, 
             bridge->droptol, chvmanager, bridge->ownersIV, 
             bridge->lookahead, perror, cpus, stats, bridge->msglvl, 
             bridge->msgFile, firsttag, bridge->comm) ;
MARKTIME(t2) ;
firsttag += 4*bridge->frontmtx->nfront ;
IVcopy(6, bridge->stats, stats) ;
bridge->cpus[11] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n CPU %8.3f : factor matrix, %8.3f mflops",
           t2 - t1, 1.e-6*nfops/(t2-t1)) ;
   fprintf(msgFile,
           "\n %8d pivots, %8d pivot tests, %8d delayed vertices"
           "\n %d entries in D, %d entries in L, %d entries in U",
           stats[0], stats[1], stats[2], stats[3], stats[4], stats[5]) ;
   cputotal = cpus[8] ;
   if ( cputotal > 0.0 ) {
      fprintf(msgFile,
      "\n    initialize fronts       %8.3f %6.2f"
      "\n    load original entries   %8.3f %6.2f"
      "\n    update fronts           %8.3f %6.2f"
      "\n    assemble postponed data %8.3f %6.2f"
      "\n    factor fronts           %8.3f %6.2f"
      "\n    extract postponed data  %8.3f %6.2f"
      "\n    store factor entries    %8.3f %6.2f"
      "\n    miscellaneous           %8.3f %6.2f"
      "\n    total time              %8.3f",
      cpus[0], 100.*cpus[0]/cputotal,
      cpus[1], 100.*cpus[1]/cputotal,
      cpus[2], 100.*cpus[2]/cputotal,
      cpus[3], 100.*cpus[3]/cputotal,
      cpus[4], 100.*cpus[4]/cputotal,
      cpus[5], 100.*cpus[5]/cputotal,
      cpus[6], 100.*cpus[6]/cputotal,
      cpus[7], 100.*cpus[7]/cputotal, cputotal) ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n submatrix mananger after factorization") ;
   SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
   fprintf(msgFile, "\n\n chevron mananger after factorization") ;
   ChvManager_writeForHumanEye(chvmanager, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n front factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
ChvManager_free(chvmanager) ;
if ( *perror >= 0 ) {
   return(0) ;
}
/*
   -----------------------------
   post-process the front matrix
   -----------------------------
*/
MARKTIME(t1) ;
IVzero(4, stats) ;
FrontMtx_MPI_postProcess(frontmtx, bridge->ownersIV, stats,
                         msglvl, msgFile, firsttag, bridge->comm) ;
MARKTIME(t2) ;
bridge->cpus[12] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n\n CPU %8.3f : post-process the matrix", t2 - t1) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n submatrix mananger after post-processing") ;
   SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n front factor matrix after post-processing") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/

MARKTIME(t2) ;
bridge->cpus[13] += t2 - t0 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n CPU %8.3f : total factor time", t2 - t0) ;
   fflush(msgFile) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
