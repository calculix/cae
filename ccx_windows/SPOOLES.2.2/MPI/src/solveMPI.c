/*  solveMPI.c  */

#include "../spoolesMPI.h"
#include "../../timings.h"

#define AGGREGATE 1
#define SOLUTION  2

/*--------------------------------------------------------------------*/
typedef struct _SolveMsg  SolveMsg ;
struct _SolveMsg {
   int           info[3] ; /* type, frontid, nbytes */
   void          *base   ;
   SubMtx        *mtx    ;
   SolveMsg      *next   ;
   MPI_Request   req     ;
} ;
/*--------------------------------------------------------------------*/
static SolveMsg * SolveMsg_new ( void ) ;
static void checkForAggMessages ( FrontMtx *frontmtx, int J, int nrhs, 
   int owners[], int myid, char recvsPosted[], int inAggCounts[], 
   SolveMsg *p_msg[], SubMtxManager *mtxmanager, SubMtxList *aggList, 
   int stats[], int firsttag, MPI_Comm comm, int msglvl, FILE *msgFile);
static SubMtx * checkForSolMessages ( FrontMtx *frontmtx, int nrhs, 
   int J, int owners[], int myid, IP *heads[], char frontIsDone[],
   char recvsPosted[], SolveMsg *p_msg[], SubMtx *p_mtx[], int nUsed[],
   SubMtxManager *mtxmanager, int stats[], int firsttag, 
   MPI_Comm comm, int msglvl, FILE *msgFile ) ;
static SolveMsg * checkSentMessages ( SolveMsg *firstmsg,
   SubMtxManager *mtxmanager, int msglvl, FILE *msgFile ) ;
static SolveMsg * sendMessages ( int J, int owners[], int myid,
   SubMtx *p_mtx[], IVL *solveIVL, SubMtxList *aggList, 
   SubMtxManager *mtxmanager, SolveMsg *firstmsg, int stats[],
   int firsttag, MPI_Comm comm, int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   MPI solve method for (L + I)D(I + U)X = B or (U^T + I)D(I + U)X = B

   created -- 98may21, ca
   -------------------------------------------------------------------
*/
void
FrontMtx_MPI_solve (
   FrontMtx        *frontmtx,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   SubMtxManager   *mtxmanager,
   SolveMap        *solvemap,
   double          cpus[],
   int             stats[],
   int             msglvl,
   FILE            *msgFile,
   int             firsttag,
   MPI_Comm        comm
) {
char         *frontIsDone, *recvsPosted, *status ;
SubMtx       *mtx, *releaseHead ;
SubMtx       **p_mtx ;
SubMtxList   *aggList ;
double       t0, t1, t2, t3 ;
Ideq         *dequeue ;
int          I, J, K, lasttag, myid, nfront, nJ, nproc, nrhs, tagbound ;
int          *fch, *inAggCounts, *nactiveChild, *nUsed, 
             *owners, *par, *sib ;
IP           *ip ;
IP           **heads ;
IV           *aggCountsIV ;
IVL          *solveIVL ;
SolveMsg     *firstmsg ;
SolveMsg     **p_msg ;
Tree         *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || mtxX == NULL || mtxB == NULL 
   || mtxmanager == NULL || solvemap == NULL 
   || cpus == NULL || stats == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_solve()"
           "\n bad input\n") ;
   exit(-1) ;
}
MARKTIME(t0) ;
nfront = frontmtx->nfront ;
tree   = frontmtx->tree ;
par    = tree->par ;
fch    = tree->fch ;
sib    = tree->sib ;
nrhs   = mtxB->ncol ;
owners = SolveMap_owners(solvemap) ;
/*
   ------------------------------
   get id and number of processes
   ------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
/*
   -------------------
   check the tag range
   -------------------
*/
lasttag  = firsttag + 2*nproc ;
tagbound = maxTagMPI(comm) ;
if ( firsttag < 0 || lasttag > tagbound ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_solve()"
           "\n firsttag = %d, lasttag = %d, tag_bound = %d", 
           firsttag, lasttag, tagbound) ;
   exit(-1) ;
}
/*
   ---------------------------------------------------------------
   allocate the working data that lasts for the entire solve

   recvsPosted[J] == 'Y' --> any receives for J have been posted,
                             receives can only be posted once
   frontIsDone[J] == 'Y' --> front J is done w.r.t. this processor
                             note, if |J| = 0, front is done
   status[J] = 'W' --> J is asleep or active
               'F' --> J is finished
   p_msg[J] = pointer to first message for J
      if ( J owned ) {
         p_msg[J] = first in list of incoming aggregates
      } if ( J used ) {
         p_msg[J] = message to hold X_J
      }
   nUsed[*] -- used to manage scope of unowned X_J
      if ( J not owned but X_J used ) {
         nUsed[J] = # of remaining times to use X_J
      }
   ---------------------------------------------------------------
*/
recvsPosted = CVinit(nfront, 'N') ;
frontIsDone = CVinit(nfront, 'W') ;
status      = CVinit(nfront, 'W') ;
ALLOCATE(p_msg, SolveMsg *, nfront) ;
for ( J = 0 ; J < nfront ; J++ ) {
   p_msg[J] = NULL ;
}
nUsed = IVinit(nfront, 0) ;
/*
   ----------------------------------------------------------------
   set up the working data specific to the forward solve

   solveIVL     -- holds lists of processors that need XJ objects
   aggCountsIV  -- vector of aggregate counts
   aggList      -- object to manage aggregate instances
   heads[*]     -- linked lists for updates local to this processor
   dequeue      -- object to manager bottom-up traversal
   nactiveChild -- manages placement of parents on dequeue
   ----------------------------------------------------------------
*/
MARKTIME(t1) ;
solveIVL = SolveMap_lowerSolveIVL(solvemap, myid, msglvl, msgFile) ;
aggCountsIV = SolveMap_lowerAggregateIV(solvemap, myid, 
                                        msglvl, msgFile) ;
inAggCounts = IV_entries(aggCountsIV) ;
aggList = SubMtxList_new() ;
SubMtxList_init(aggList, nfront, inAggCounts, NO_LOCK, NULL) ;
heads   = SolveMap_forwardSetup(solvemap, myid, msglvl, msgFile) ;
dequeue = FrontMtx_setUpDequeue(frontmtx, owners, myid, status,
                                heads, 'W', 'F', msglvl, msgFile) ;
FrontMtx_loadActiveLeaves(frontmtx, status, 'W', dequeue) ;
nactiveChild = FrontMtx_nactiveChild(frontmtx, status, myid) ;
for ( J = 0 ; J < nfront ; J++ ) {
   for ( ip = heads[J] ; ip != NULL ; ip = ip->next ) {
      I = ip->val ;
      if ( owners[I] != myid ) {
         nUsed[I]++ ;
      }
   }
}
if ( msglvl > 2 ) {
   IP   *ip ;
   fprintf(msgFile, "\n\n ### setup for forward solve") ;
   fprintf(msgFile, "\n\n solveIVL ") ;
   IVL_writeForHumanEye(solveIVL, msgFile) ;
   fprintf(msgFile, "\n\n aggList ") ;
   SubMtxList_writeForHumanEye(aggList, msgFile) ;
   fprintf(msgFile, "\n\n initial dequeue") ;
   Ideq_writeForHumanEye(dequeue, msgFile) ;
   fprintf(msgFile, "\n\n forward update lists") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( heads[J] != NULL ) {
         fprintf(msgFile, "\n %d : ", J) ;
         for ( ip = heads[J] ; ip != NULL ; ip = ip->next ) {
            fprintf(msgFile, " %d", ip->val) ;
         }
      }
   }
   fprintf(msgFile, "\n\n nactiveChild") ;
   IVfprintf(msgFile, nfront, nactiveChild) ;
   fflush(msgFile) ;
}
MARKTIME(t2) ;
cpus[0] += t2 - t1 ;
/*
   --------------------------------------------------------
   load the right hand side, p_mtx[*] is returned
   if ( J owned ) {
      p_mtx[J] is initially B_J and overwritten with X_J
   } else {
      if ( J is supported by q ) {
         p_mtx[J] contains Y_J^q until sent to owner of Y_J
      }
      if ( XJ is used by q ) {
         p_mtx[J] contains X_J when received
      }
   }
   --------------------------------------------------------
*/
MARKTIME(t1) ;
p_mtx = FrontMtx_loadRightHandSide(frontmtx, mtxB, owners, myid,
                                   mtxmanager, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[1] += t2 - t1 ;
/*
   -------------
   forward solve
   -------------
*/
MARKTIME(t1) ;
firstmsg = NULL ;
while ( (J = Ideq_removeFromHead(dequeue)) != -1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile,
              "\n\n ### forward solve, checking out J = %d, status %c",
              J, status[J]) ;
      fflush(msgFile) ;
   }
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
/*
      -------------------------------------
      check for necessary incoming solution 
      matrices needed to compute X_J
      -------------------------------------
*/
      releaseHead = checkForSolMessages(frontmtx, nrhs, J, owners, myid,
                              heads, frontIsDone, recvsPosted, p_msg, 
                              p_mtx, nUsed, mtxmanager, stats, firsttag,
                              comm, msglvl, msgFile) ;
      if ( owners[J] == myid ) {
/*
         -------------------------------------
         check for incoming aggregate messages
         -------------------------------------
*/
         checkForAggMessages(frontmtx, J, nrhs, owners, myid, 
                      recvsPosted, inAggCounts, p_msg, mtxmanager, 
                      aggList, stats, firsttag, comm, msglvl, msgFile) ;
      }
/*
      ---------------
      visit the front
      ---------------
*/
      FrontMtx_forwardVisit(frontmtx, J, nrhs, owners, myid, mtxmanager,
                            aggList, p_mtx, frontIsDone, heads, 
                            p_mtx, status, msglvl, msgFile) ;
   } else {
/*
      ----------------------------------------
      front has no eliminated rows and columns
      ----------------------------------------
*/
      status[J] = 'F' ;
      releaseHead = NULL ;
   }
   if ( status[J] == 'F' ) {
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n front J is finished") ;
         fflush(msgFile) ;
      }
      if ( nJ > 0 ) {
/*
         ------------------------------------------------------------
         front J is finished, send any solution or aggregate messages
         ------------------------------------------------------------
*/
         firstmsg = sendMessages(J, owners, myid, p_mtx, solveIVL, 
                                 aggList, mtxmanager, firstmsg, stats,
                                 firsttag, comm, msglvl, msgFile) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n after sendMessages, p_mtx[%d] = %p",
                    J, p_mtx[J]) ;
            fflush(msgFile) ;
         }
      }
      if ( (K = par[J]) != -1 && --nactiveChild[K] == 0 ) {
/*
         -----------------------------------
         all children of parent are finished
         add parent to head of dequeue
         -----------------------------------
*/
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n placing %d on head of dequeue", K) ;
            fflush(msgFile) ;
         }
         Ideq_insertAtHead(dequeue, K) ;
      }
   } else {
/*
      -------------------------------------------------
      front J is not finished, place on tail of dequeue
      -------------------------------------------------
*/
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n front J is not finished") ;
         fflush(msgFile) ;
      }
      Ideq_insertAtTail(dequeue, J) ;
   }
/*
   --------------------------------------------
   check all posted sends, any message that 
   has been received, recycle its SubMtx object
   --------------------------------------------
*/
   firstmsg = checkSentMessages(firstmsg, mtxmanager, msglvl, msgFile) ;
/*
   ----------------------------------------------------------
   release all external X_J objects that are no longer needed
   ----------------------------------------------------------
*/
   while ( (mtx = releaseHead) != NULL ) {
      p_mtx[mtx->rowid] = NULL ;
      releaseHead = mtx->next ;
      SubMtxManager_releaseObject(mtxmanager, mtx) ;
   }
}
while ( firstmsg != NULL ) {
/*
   ------------------------------------------
   check all posted sends, any message that 
   has been received, recycle its SubMtx object
   ------------------------------------------
*/
   firstmsg = checkSentMessages(firstmsg, mtxmanager, msglvl, msgFile) ;
}
IV_free(aggCountsIV) ;
IP_free(heads[nfront+1]) ;
FREE(heads) ;
IVfree(nactiveChild) ;
Ideq_free(dequeue) ;
SubMtxList_free(aggList) ;
IVL_free(solveIVL) ;
MARKTIME(t2) ;
cpus[2] += t2 - t1 ;
/*
   ---------------------
   do the diagonal solve
   ---------------------
*/
MARKTIME(t1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( owners[J] == myid ) {
      FrontMtx_diagonalVisit(frontmtx, J, owners, myid, p_mtx,
                              frontIsDone, p_mtx, msglvl, msgFile) ;
   }
}
MARKTIME(t2) ;
cpus[3] += t2 - t1 ;
/*
   ------------------------------------
   set up for the backward solve
   (1) IP lists
   (2) set up the dequeue
   (3) set up the aggregate list object
   (4) nUsed[] vector
   ------------------------------------
*/
/*
   ----------------------------------------------------------------
   set up the working data specific to the forward solve

   solveIVL     -- holds lists of processors that need XJ objects
   aggCountsIV  -- vector of aggregate counts
   aggList      -- object to manage aggregate instances
   heads[*]     -- linked lists for updates local to this processor
   dequeue      -- object to manager bottom-up traversal
   ----------------------------------------------------------------
*/
MARKTIME(t1) ;
solveIVL = SolveMap_upperSolveIVL(solvemap, myid, msglvl, msgFile) ;
aggCountsIV = SolveMap_upperAggregateIV(solvemap, myid,
                                        msglvl, msgFile) ;
inAggCounts = IV_entries(aggCountsIV) ;
aggList = SubMtxList_new() ;
SubMtxList_init(aggList, nfront, inAggCounts, 0, NULL) ;
heads = SolveMap_backwardSetup(solvemap, myid, msglvl, msgFile) ;
dequeue = FrontMtx_setUpDequeue(frontmtx, owners, myid, status, 
                                heads, 'W', 'F', msglvl, msgFile) ;
FrontMtx_loadActiveRoots(frontmtx, status, 'W', dequeue) ;
CVfill(nfront, recvsPosted, 'N') ;
CVfill(nfront, frontIsDone, 'N') ;
for ( J = 0 ; J < nfront ; J++ ) {
   for ( ip = heads[J] ; ip != NULL ; ip = ip->next ) {
      I = ip->val ;
      if ( owners[I] != myid ) {
         nUsed[I]++ ;
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n ### setup for backward solve") ;
   fprintf(msgFile, "\n\n solveIVL ") ;
   IVL_writeForHumanEye(solveIVL, msgFile) ;
   fprintf(msgFile, "\n\n aggList ") ;
   SubMtxList_writeForHumanEye(aggList, msgFile) ;
   fprintf(msgFile, "\n\n initial dequeue") ;
   Ideq_writeForHumanEye(dequeue, msgFile) ;
   fprintf(msgFile, "\n\n forward update lists") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      IP   *ip ;
      fprintf(msgFile, "\n %d :", J) ;
      fflush(msgFile) ;
      for ( ip = heads[J] ; ip != NULL ; ip = ip->next ) {
         fprintf(msgFile, " %d", ip->val) ;
         fflush(msgFile) ;
      }
   }
}
MARKTIME(t2) ;
cpus[0] += t2 - t1 ;
/*
   --------------
   backward solve
   --------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n thread %d : starting backward solve", myid) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
firsttag += nfront ;
while ( (J = Ideq_removeFromHead(dequeue)) != -1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile,
              "\n\n ### backward solve, checking out J = %d, status %c",
              J, status[J]) ;
      fflush(msgFile) ;
   }
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
/*
      ----------------------------------------------
      check for necessary incoming solution matrices
      ----------------------------------------------
*/
      releaseHead = checkForSolMessages(frontmtx, nrhs, J, owners, myid,
                              heads, frontIsDone, recvsPosted, p_msg, 
                              p_mtx, nUsed, mtxmanager, stats, firsttag,
                              comm, msglvl, msgFile) ;
      if ( owners[J] == myid ) {
/*
         -------------------------------------
         check for incoming aggregate messages
         -------------------------------------
*/
         checkForAggMessages(frontmtx, J, nrhs, owners, myid, 
                      recvsPosted, inAggCounts, p_msg, mtxmanager, 
                      aggList, stats, firsttag, comm, msglvl, msgFile) ;
      }
/*
      ---------------
      visit the front
      ---------------
*/
      FrontMtx_backwardVisit(frontmtx, J, nrhs, owners, myid, 
                             mtxmanager, aggList, p_mtx, frontIsDone, 
                             heads, p_mtx, status, msglvl, msgFile) ;
   } else {
/*
      ----------------------------------------
      front has no eliminated rows and columns
      ----------------------------------------
*/
      status[J] = 'F' ;
      releaseHead = NULL ;
   }
   if ( status[J] == 'F' ) {
      if ( nJ > 0 ) {
/*
         ------------------------------------------------------------
         front J is finished, send any solution or aggregate messages
         ------------------------------------------------------------
*/
         firstmsg = sendMessages(J, owners, myid, p_mtx, solveIVL, 
                                 aggList, mtxmanager, firstmsg, stats,
                                 firsttag, comm, msglvl, msgFile) ;
      }
/*
      ---------------------------------
      load any children on active paths
      ---------------------------------
*/
      for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
         if ( status[I] == 'W' ) {
            Ideq_insertAtHead(dequeue, I) ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, "\n placing %d at head of dequeue", I) ;
               fflush(msgFile) ;
            }
         }
      }
   } else {
/*
      -------------------------------------------------
      front J is not finished, place on tail of dequeue
      -------------------------------------------------
*/
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n front J is not finished") ;
         fflush(msgFile) ;
      }
      Ideq_insertAtTail(dequeue, J) ;
   }
/*
   --------------------------------------------
   check all posted sends, any message that 
   has been received, recycle its SubMtx object
   --------------------------------------------
*/
   firstmsg = checkSentMessages(firstmsg, mtxmanager, msglvl, msgFile) ;
/*
   ----------------------------------------------------------
   release all external X_J objects that are no longer needed
   ----------------------------------------------------------
*/
   while ( (mtx = releaseHead) != NULL ) {
      p_mtx[mtx->rowid] = NULL ;
      releaseHead = mtx->next ;
      SubMtxManager_releaseObject(mtxmanager, mtx) ;
   }
}
while ( firstmsg != NULL ) {
/*
   ------------------------------------------
   check all posted sends, any message that 
   has been received, recycle its SubMtx object
   ------------------------------------------
*/
   firstmsg = checkSentMessages(firstmsg, mtxmanager, msglvl, msgFile) ;
}
IV_free(aggCountsIV) ;
IP_free(heads[nfront+1]) ;
FREE(heads) ;
Ideq_free(dequeue) ;
SubMtxList_free(aggList) ;
IVL_free(solveIVL) ;
MARKTIME(t2) ;
cpus[4] += t2 - t1 ;
/*
   --------------------------
   store the solution entries
   --------------------------
*/
MARKTIME(t1) ;
FrontMtx_storeSolution(frontmtx, owners, myid, mtxmanager,
                        p_mtx, mtxX, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[1] += t2 - t1 ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
FREE(p_mtx) ;
FREE(p_msg) ;
IVfree(nUsed) ;
CVfree(status) ;
CVfree(recvsPosted) ;
CVfree(frontIsDone) ;

MARKTIME(t3) ;
cpus[5] = t3 - t0 - cpus[0] - cpus[1] - cpus[2] - cpus[3] - cpus[4] ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   this method is called inside the event loop for the forward 
   and backward solves. 
   (1) it posts and checks for external X_I messages
       where X_I is needed to update Y_J
   (2) it returns a pointer to the first SubMtx object in a list
       of external X_I objects that can be free'd.

   created -- 98apr02, cca
   -----------------------------------------------------------
*/
static SubMtx *
checkForSolMessages (
   FrontMtx        *frontmtx,
   int             nrhs,
   int             J,
   int             owners[],
   int             myid,
   IP              *heads[],
   char            frontIsDone[],
   char            recvsPosted[],
   SolveMsg        *p_msg[],
   SubMtx          *p_mtx[],
   int             nUsed[],
   SubMtxManager   *mtxmanager,
   int             stats[],
   int             firsttag,
   MPI_Comm        comm,
   int             msglvl,
   FILE            *msgFile
) {
SubMtx       *firstToRelease, *XI ;
int          flag, I, nbytes, nI, source, tag ;
IP           *ip ;
SolveMsg     *msg ;
MPI_Status   status ;

if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside checkForSolMessages(%d)", J) ;
   fflush(msgFile) ;
}
firstToRelease = NULL ;
for ( ip = heads[J] ; ip != NULL ; ip = ip->next ) {
   I = ip->val ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n I = %d, source = %d", I, owners[I]) ;
      fflush(msgFile) ;
   }
   if ( (source = owners[I]) != myid ) {
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n nI = %d", 
                 FrontMtx_frontSize(frontmtx, I)) ;
         fflush(msgFile) ;
      }
      if ( (nI = FrontMtx_frontSize(frontmtx, I)) == 0 ) {
/*
         -------------------------------------------------
         there is no X_I to receive, set local status flag
         -------------------------------------------------
*/
         frontIsDone[I] = 'Y' ;
      } else {
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n recvsPosted[%d] = %c", 
                    I, recvsPosted[I]) ;
            fflush(msgFile) ;
         }
         if ( recvsPosted[I] == 'N' ) {
/*
            --------------------
            post receive for X_I
            --------------------
*/
            msg = SolveMsg_new() ;
            msg->info[0] = SOLUTION ;
            msg->info[1] = I ;
            msg->info[2] = 0 ;
            nbytes = SubMtx_nbytesNeeded(frontmtx->type,
                         SUBMTX_DENSE_COLUMNS, nI, nrhs, nI*nrhs) ;
            XI = SubMtxManager_newObjectOfSizeNbytes(mtxmanager, 
                                                     nbytes);
            SubMtx_init(XI, frontmtx->type, SUBMTX_DENSE_COLUMNS, 
                        I, 0, nI, nrhs, nI*nrhs) ;
            msg->mtx  = XI ;
            msg->base = SubMtx_workspace(XI) ;
            tag = firsttag + I ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, 
                       "\n posting solution Irecv,"
                       "msg %p, I %d, tag %d, nbytes %d",
                       msg, I, tag, nbytes) ;
               fflush(msgFile) ;
            }
            MPI_Irecv(msg->base, nbytes, MPI_BYTE, source, tag,
                      comm, &msg->req) ;
            p_msg[I] = msg ;
            recvsPosted[I] = 'Y' ;
            stats[4]++ ; stats[6] += nbytes ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, 
        "\n post Irecv XI : I %d, source %d, stats[4] %d, stats[6] %d",
        I, source, stats[4], stats[6]) ;
               fflush(msgFile) ;
            }
         }
         if ( (msg = p_msg[I]) != NULL ) {
/*
            --------------------------------
            check for message containing X_I
            --------------------------------
*/
            if ( msglvl > 1 ) {
               fprintf(msgFile, "\n testing for X_{%d}", I) ;
               fflush(msgFile) ;
            }
            MPI_Test(&msg->req, &flag, &status) ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, ", flag = %d", flag) ;
               fflush(msgFile) ;
            }
            if ( flag != 0 ) {
               source = status.MPI_SOURCE ;
               tag    = status.MPI_TAG    ;
               if ( msglvl > 1 ) {
                  fprintf(msgFile, 
                          "\n X_{%d} received from source %d", 
                          I, source) ;
                  fflush(msgFile) ;
               }
               if ( msg->mtx == NULL ) {
                  fprintf(stderr, 
            "\n proc %d: fatal error in checkForSolMsg, J = %d, I = %d",
                  myid, J, I) ;
                  exit(-1) ;
               }
               SubMtx_initFromBuffer(msg->mtx) ;
               p_msg[I] = NULL ;
               p_mtx[I] = msg->mtx ;
               frontIsDone[I] = 'Y' ;
               FREE(msg) ;
            }
         }
         if ( (XI = p_mtx[I]) != NULL && --nUsed[I] == 0 ) {
/*
            ----------------------------------
            we can release X_I after this step
            ----------------------------------
*/
            XI->next = firstToRelease ;
            firstToRelease = XI ;
         }
      }
   }
}
return(firstToRelease) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   this method is called inside the event loop for the forward and 
   backward solves. it posts and checks for aggregate Y_J^q messages

   created -- 98apr02, cca
   -----------------------------------------------------------------
*/
static void
checkForAggMessages (
   FrontMtx        *frontmtx,
   int             J,
   int             nrhs,
   int             owners[],
   int             myid,
   char            recvsPosted[],
   int             inAggCounts[],
   SolveMsg        *p_msg[],
   SubMtxManager   *mtxmanager,
   SubMtxList      *aggList,
   int             stats[],
   int             firsttag,
   MPI_Comm        comm,
   int             msglvl,
   FILE            *msgFile
) {
SubMtx         *YJ ;
int          count, flag, ii, nbytes, nincoming, nJ, tag ;
SolveMsg     *msg, *nextmsg ;
MPI_Status   status ;

if ( owners[J] == myid && recvsPosted[J] == 'N' ) {
/*
   ---------------------------------------------
   post the receives for the incoming aggregates
   ---------------------------------------------
*/
   if ( (nincoming = inAggCounts[J]) > 0 ) {
/*
      -----------------------
      post receives for Y_J^q
      -----------------------
*/
      nJ = FrontMtx_frontSize(frontmtx, J) ;
      for ( ii = 0 ; ii < nincoming ; ii++ ) {
         nbytes = SubMtx_nbytesNeeded(frontmtx->type,
                              SUBMTX_DENSE_COLUMNS, nJ, nrhs, nJ*nrhs) ;
         YJ = SubMtxManager_newObjectOfSizeNbytes(mtxmanager, nbytes);
         SubMtx_init(YJ, frontmtx->type, SUBMTX_DENSE_COLUMNS, 
                     J, 0, nJ, nrhs, nJ*nrhs) ;
         msg = SolveMsg_new() ;
         msg->info[0] = AGGREGATE ;
         msg->info[1] = J ;
         msg->info[2] = nbytes ;
         msg->base = SubMtx_workspace(YJ) ;
         msg->mtx  = YJ ;
         tag = firsttag + J ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                    "\n posting aggregate Irecv,"
                    "msg %p, J %d, tag %d, nbytes %d",
                    msg, J, tag, nbytes) ;
            fflush(msgFile) ;
         }
         MPI_Irecv(msg->base, nbytes, MPI_BYTE, MPI_ANY_SOURCE, tag,
                   comm, &msg->req) ;
         msg->next = p_msg[J] ;
         p_msg[J] = msg ;
         stats[5]++ ; stats[7] += nbytes ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                    "\n post Irecv YJ : J %d, stats[5] %d, stats[7] %d",
                    J, stats[5], stats[7]) ;
            fflush(msgFile) ;
         }
      }
   }
   recvsPosted[J] = 'Y' ;
}
if ( p_msg[J] != NULL ) {
/*
   ---------------------------------------
   check for unreceived aggregate messages
   ---------------------------------------
*/
   for ( msg = p_msg[J], p_msg[J] = NULL ; 
         msg != NULL ;
         msg = nextmsg ) {
      nextmsg = msg->next ;
      if ( msglvl > 1 ) {
         fprintf(msgFile,
                 "\n checking aggregate message %p, J %d",
                 msg, J) ;
         fflush(msgFile) ;
      }
      MPI_Test(&msg->req, &flag, &status) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, ", flag = %d", flag) ;
         fflush(msgFile) ;
      }
      if ( flag == 0 ) {
/*
         --------------------------------------
         message not received, put back on list
         --------------------------------------
*/
         msg->next = p_msg[J] ;
         p_msg[J] = msg ;
      } else {
/*
         -------------------------
         message received, 
         initialize SubMtx object 
         and put on aggregate list
         -------------------------
*/
         if ( msglvl > 1 ) {
            MPI_Get_count(&status, MPI_BYTE, &count) ;
            fprintf(msgFile,
                    "\n message received, source %d, tag %d, nbytes %d",
                    status.MPI_SOURCE, status.MPI_TAG, count) ;
            fflush(msgFile) ;
         }
         if ( msg->mtx == NULL ) {
            fprintf(stderr, 
                   "\n proc %d: fatal error in checkForAggMsg, J = %d", 
                   myid, J) ;
            exit(-1) ;
         }
         SubMtx_initFromBuffer(msg->mtx) ;
         SubMtxList_addObjectToList(aggList, msg->mtx, J) ;
         FREE(msg) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   this method is called inside the event loop for the forward and 
   backward solves when front J is now complete. if J is owned, then
   X_J is sent to all processors that need it. (this information is
   found in the solveIVL object.) otherwise, Y_J^q exists in the
   aggregate list object. Y_J^q is removed and sent to the owner of J.

   return value -- any new messages are prepended onto the list of
                   sent messages and the head of the list is returned.

   created -- 98apr02, cca
   -------------------------------------------------------------------
*/
static SolveMsg *
sendMessages (
   int             J,
   int             owners[],
   int             myid,
   SubMtx          *p_mtx[] ,
   IVL             *solveIVL,
   SubMtxList      *aggList,
   SubMtxManager   *mtxmanager,
   SolveMsg        *firstmsg,
   int             stats[],
   int             firsttag,
   MPI_Comm        comm,
   int             msglvl,
   FILE            *msgFile
) {
SubMtx       *mtx, *XJ, *YJ ;
int        dest, ii, mproc, nbytes, tag ;
int        *procids ;
SolveMsg   *msg ;
void       *work, *XJbuff ;

if ( (dest = owners[J]) == myid ) {
/*
   ---------------------------------------
   send X_J to all processors that need it
   ---------------------------------------
*/
   IVL_listAndSize(solveIVL, J, &mproc, &procids) ;
   if ( mproc > 0 ) {
      XJ = p_mtx[J] ;
      if ( XJ == NULL ) {
         fprintf(stderr, 
                "\n proc %d: error in sendMessages, J = %d, XJ is NULL",
                myid, J) ;
         exit(-1) ;
      }
      nbytes = SubMtx_nbytesInUse(XJ) ;
      XJbuff = SubMtx_workspace(XJ) ;
      tag = firsttag + J ;
      for ( ii = 0 ; ii < mproc ; ii++ ) {
         dest = procids[ii] ;
         mtx = SubMtxManager_newObjectOfSizeNbytes(mtxmanager, nbytes) ;
         work = SubMtx_workspace(mtx) ;
         memcpy(work, XJbuff, nbytes) ;
         msg = SolveMsg_new() ;
         msg->info[0] = SOLUTION ;
         msg->info[1] = J ;
         msg->info[2] = nbytes ;
         msg->base    = work ;
         msg->mtx     = mtx  ;
         msg->next    = firstmsg ;
         firstmsg  = msg ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                    "\n posting X_J Isend, "
                    "msg %p, J %d, dest %d, tag %d, nbytes %d",
                    msg, J, dest, tag, nbytes) ;
            fflush(msgFile) ;
         }
         MPI_Isend(msg->base, nbytes, MPI_BYTE, dest, tag, 
                   comm, &msg->req) ;
         stats[0]++, stats[2] += nbytes ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
           "\n post Isend XJ : J %d, dest %d, stats[0] %d, stats[2] %d",
           J, dest, stats[0], stats[2]) ;
            fflush(msgFile) ;
         }
      }
   }
} else {
/*
   --------------------------
   send Y_J^myid to the owner
   --------------------------
*/
   YJ = SubMtxList_getList(aggList, J) ;
/*
   if ( YJ == NULL ) {
      fprintf(stderr, 
             "\n proc %d: error in sendMessages, J = %d, YJ is NULL",
             myid, J) ;
      exit(-1) ;
   }
*/
   if ( YJ != NULL ) {
      msg = SolveMsg_new() ;
      msg->info[0] = AGGREGATE ;
      msg->info[1] = J ;
      msg->info[2] = nbytes = SubMtx_nbytesInUse(YJ) ;
      msg->mtx  = YJ ;
      msg->base = SubMtx_workspace(YJ) ;
      msg->next = firstmsg ;
      tag = firsttag + J ;
      firstmsg  = msg ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
                 "\n posting Y_J^%d Isend, "
                 "msg %p, J %d, dest %d, tag %d, nbytes %d",
                 myid, msg, J, dest, tag, nbytes) ;
         fflush(msgFile) ;
      }
      MPI_Isend(msg->base, nbytes, MPI_BYTE, dest, tag, 
                comm, &msg->req) ;
      p_mtx[J] = NULL ;
      stats[1]++, stats[3] += nbytes ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
           "\n post Isend YJ : J %d, dest %d, stats[1] %d, stats[3] %d",
           J, dest, stats[1], stats[3]) ;
         fflush(msgFile) ;
      }
   }
}
return(firstmsg) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   this method is called inside the event loop for the forward and 
   backward solves when front J is now complete. it checks the list
   of sent messages. any messages that have been received have their
   SubMtx object given back to the manager and their SolveMsg object
   free'd.

   return value -- head of list of unreceived messages

   created -- 98apr02, cca
   -----------------------------------------------------------------
*/
static SolveMsg *
checkSentMessages (
   SolveMsg      *firstmsg,
   SubMtxManager   *mtxmanager,
   int           msglvl,
   FILE          *msgFile
) {
int          flag ;
MPI_Status   status ;
SolveMsg     *msg, *nextmsg ;

for ( msg = firstmsg, firstmsg = NULL ; msg != NULL ; msg = nextmsg ) {
   nextmsg = msg->next ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, 
           "\n checking sent message %p : type %d, front %d, nbytes %d",
           msg, msg->info[0], msg->info[1], msg->info[2]) ;
      fflush(msgFile) ;
   }
   MPI_Test(&msg->req, &flag, &status) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, ", flag = %d", flag) ;
      fflush(msgFile) ;
   }
   if ( flag == 1 ) {
      if ( msg->mtx == NULL ) {
         fprintf(msgFile, "\n WHOA!, msg = %p, msg->mtx = NULL", msg) ;
         fflush(msgFile) ;
      } else {
         SubMtxManager_releaseObject(mtxmanager, msg->mtx) ;
      }
      FREE(msg) ;
   } else {
      msg->next = firstmsg ;
      firstmsg  = msg ;
   }
}
return(firstmsg) ; }

/*--------------------------------------------------------------------*/
static SolveMsg *
SolveMsg_new (
   void 
) {
SolveMsg   *msg ;

ALLOCATE(msg, struct _SolveMsg, 1) ;
msg->info[0] = 0 ;
msg->info[1] = 0 ;
msg->info[2] = 0 ;
msg->base    = NULL ;
msg->next    = NULL ;
msg->mtx     = NULL ;

return(msg) ; }

/*--------------------------------------------------------------------*/
