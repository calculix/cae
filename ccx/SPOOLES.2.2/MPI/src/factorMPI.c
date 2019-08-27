/*  factorMPI.c  */

#include "../spoolesMPI.h"
#include "../../timings.h"

#define MYDIAGNOSTIC 0

#define AGGREGATE_CHV          1
#define POSTPONED_NOTIFICATION 2
#define POSTPONED_CHV          3
#define AGGREGATE_CHV_TAG(K)          (firsttag + (K))
#define POSTPONED_NOTIFICATION_TAG(K) (firsttag + nfront + (K))
#define POSTPONED_CHV_TAG(K)          (firsttag + 2*nfront + 1 + (K))
#define ERROR_TAG                     (firsttag + 3*nfront + 2)

/*--------------------------------------------------------------------*/
typedef struct _FactorMsg FactorMsg ;
struct _FactorMsg {
   int           info[3] ;
   void          *base   ;
   Chv           *chv    ;
   MPI_Request   req     ;
   FactorMsg     *next   ;
} ;
/*--------------------------------------------------------------------*/
static FactorMsg * FactorMsg_new ( void ) ;
static void FactorMsg_setDefaultFields ( FactorMsg *msg ) ;
static void FactorMsg_clearData ( FactorMsg *msg ) ;
static void FactorMsg_free ( FactorMsg *msg ) ;
static FactorMsg * wakeupFront ( FrontMtx *frontmtx, int J, int myid, 
   int frontOwners[], int firsttag, ChvManager *chvmanager, 
   ChvList *aggList, int stats[], double cpus[], MPI_Comm comm, 
   int msglvl, FILE *msgFile ) ;
static FactorMsg * checkMessages ( FrontMtx *frontmtx, int J,
   FactorMsg *firstmsg, ChvList *aggList, ChvList *postList,
   ChvManager *chvmanager, int firsttag, int stats[], MPI_Comm comm,
   int msglvl, FILE *msgFile ) ;
static FactorMsg * sendMessages ( FrontMtx *frontmtx, int J, int myid, 
   int nfront, int par[], int owners[], int firsttag, 
   ChvManager*chvmanager, ChvList *aggList, ChvList *postList, 
   FactorMsg *firstmsgsent, int stats[], MPI_Comm comm, 
   int msglvl, FILE *msgFile ) ;
static FactorMsg * checkSentMessages ( FactorMsg *firstmsgsent,
   ChvManager *chvmanager, int msglvl, FILE *msgFile ) ;
static FactorMsg * sendErrorMessages ( int nfront, int J, int nproc,
   int myid, int firsttag, FactorMsg *firstmsgsent, int stats[],
   MPI_Comm comm, int msglvl, FILE *msgFile ) ;
static void cancelAndFreeMessages ( FactorMsg *head,
   ChvManager *chvmanager, int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
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
         stats[ 0] -- # of pivots
         stats[ 1] -- # of pivot tests
         stats[ 2] -- # of delayed rows and columns
         stats[ 3] -- # of entries in D
         stats[ 4] -- # of entries in L
         stats[ 5] -- # of entries in U
         stats[ 6] -- # of aggregate sends
         stats[ 7] -- # of bytes in the aggregate sends
         stats[ 8] -- # of aggregate received
         stats[ 9] -- # of bytes in the aggregate received
         stats[10] -- # of postponed data sends
         stats[11] -- # of bytes in the postponed data sends
         stats[12] -- # of postponed data received
         stats[13] -- # of bytes in the postponed data received
         stats[14] -- # of active Chv objects (working storage)
         stats[15] -- # of active bytes in working storage
         stats[16] -- # of requested bytes in working storage
      msglvl        -- message level
      msgFile       -- message file
      firsttag      -- first tag to use during the factorization,
                       reserved tags are tag, ..., tag + 4*nfront + 2
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
) {
Chv      *rootchv ;
double   zero[2] = {0.0, 0.0} ;
int      lasttag, tagbound ;
Pencil   pencil ;

lasttag = 3*frontmtx->nfront + 2 ;
tagbound = maxTagMPI(comm) ;
if ( firsttag < 0 || lasttag > tagbound ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_factorInpMtx()"
           "\n tag range is [%d,%d], tag_bound = %d",
           firsttag, lasttag, tagbound) ;
   exit(-1) ;
}
Pencil_setDefaultFields(&pencil) ;
Pencil_init(&pencil, frontmtx->type, frontmtx->symmetryflag,
            inpmtx, zero, NULL) ;
rootchv = FrontMtx_MPI_factorPencil(frontmtx, &pencil, tau, droptol,
             chvmanager, frontOwnersIV, lookahead, perror, cpus, stats, 
             msglvl, msgFile, firsttag, comm) ;

return(rootchv) ; }

/*--------------------------------------------------------------------*/
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
         stats[ 0] -- # of pivots
         stats[ 1] -- # of pivot tests
         stats[ 2] -- # of delayed rows and columns
         stats[ 3] -- # of entries in D
         stats[ 4] -- # of entries in L
         stats[ 5] -- # of entries in U
         stats[ 6] -- # of aggregate sends
         stats[ 7] -- # of bytes in the aggregate sends
         stats[ 8] -- # of aggregate received
         stats[ 9] -- # of bytes in the aggregate received
         stats[10] -- # of postponed data sends
         stats[11] -- # of bytes in the postponed data sends
         stats[12] -- # of postponed data received
         stats[13] -- # of bytes in the postponed data received
         stats[14] -- # of active Chv objects (working storage)
         stats[15] -- # of active bytes in working storage
         stats[16] -- # of requested bytes in working storage
      msglvl   -- message level
      msgFile  -- message file
      firsttag -- first tag to use during the factorization,
                  reserved tags are tag, ..., tag + 4*nfront + 2
      comm     -- MPI communicator

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
) {
char         *status, *recvsPosted ;
Chv          *rootchv ;
Chv          **fronts ;
ChvList      *aggList, *postList ;
double       t1, t2 ;
DV           workDV ;
FactorMsg    *firstmsgsent, *head, *msg ;
FactorMsg    **p_msg ;
Ideq         *dequeue ;
int          flag, iproc, J, K, lasttag, myid, nfront, nproc, tagbound ;
int          *nactiveChild, *owners, *par, *term_flags ;
int          commstats[4] ;
IP           **heads ;
IV           pivotsizesIV ;
MPI_Status   mpistatus ;
/*
   --------------
   check the data
   --------------
*/
if ( frontmtx == NULL || pencil == NULL || tau < 1.0 || droptol < 0.0 
   || chvmanager == NULL || frontOwnersIV == NULL || cpus == NULL 
   || stats == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_factorPencil()"
           "\n frontmtx = %p, pencil = %p, tau = %f, droptol = %f"
           "\n chvmanager = %p, frontOwnersIV = %p, firsttag = %d"
           "\n cpus = %p, stats = %p, msglvl = %d, msgFile = %p"
           "\n bad input\n",
           frontmtx, pencil, tau, droptol, chvmanager, frontOwnersIV, 
           firsttag, cpus, stats, msglvl, msgFile) ;
   exit(-1) ;
}
nfront  = frontmtx->nfront ;
lasttag = 3*nfront + 2 ;
tagbound  = maxTagMPI(comm) ;
if ( firsttag < 0 || lasttag > tagbound ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_factorPencil()"
           "\n tag range is [%d,%d], tag_bound = %d",
           firsttag, lasttag, tagbound) ;
   exit(-1) ;
}
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
par    = frontmtx->tree->par ;
owners = IV_entries(frontOwnersIV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n chvmanager = %p, firsttag = %d", 
           chvmanager, firsttag) ;
   fflush(msgFile) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   --------------------------------------------
   create a ChvList object to handle the 
   postponed data.  it does not require a lock.
   --------------------------------------------
*/
   postList = FrontMtx_postList(frontmtx, frontOwnersIV, NO_LOCK) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n postList = %p", postList) ;
      ChvList_writeForHumanEye(postList, msgFile) ;
      fflush(msgFile) ;
   }
} else {
   postList = NULL ;
}
/*
   ------------------------------------------------------------
   create a ChvList object to handle aggregate fronts. no lock
   is necessary. this object must be created cooperatively, 
   since the aggregate counts are a function of the owners map 
   and the symbolic factorization, and the latter is not global 
   ------------------------------------------------------------
*/
IVzero(4, commstats) ;
aggList = FrontMtx_MPI_aggregateList(frontmtx, frontOwnersIV, 
             commstats, msglvl, msgFile, firsttag, comm) ;
firsttag += 2 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n aggList = %p", aggList) ;
   ChvList_writeForHumanEye(aggList, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------------
   initialize the heads[] vector for the owned updates
   ---------------------------------------------------
*/
heads = FrontMtx_factorSetup(frontmtx, frontOwnersIV, myid,
                             msglvl, msgFile) ;
/*
   ----------------------------------------------------------------
  initialize the Ideq object that holds the initial fronts
   of the active paths, owned fronts with no children that
   are owned or updates by this thread. 
   status[J] == 'W' --> J belongs to an active path for this thread
  ----------------------------------------------------------------
*/
status = CVinit(nfront, 'F') ;
dequeue = FrontMtx_setUpDequeue(frontmtx, owners, myid, status, heads,
                                 'W', 'F', msglvl, msgFile) ;
FrontMtx_loadActiveLeaves(frontmtx, status, 'W', dequeue) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n status") ;
   CVfprintf(msgFile, nfront, status) ;
   fprintf(msgFile, "\n fronts in initial dequeue") ;
   Ideq_writeForHumanEye(dequeue, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------
   initialize the nactiveChild[] vector,
   nactiveChild[J] measures the number of children 
   that belong to active paths of this thread
   -----------------------------------------------
*/
nactiveChild = FrontMtx_nactiveChild(frontmtx, status, myid) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n nactiveChild") ;
   IVfprintf(msgFile, nfront, nactiveChild) ;
   fflush(msgFile) ;
}
/*
   ------------------------------
   initialize the working storage
   ------------------------------
*/
ALLOCATE(fronts, Chv *, nfront) ;
for ( J = 0 ; J < nfront ; J++ ) {
   fronts[J] = NULL ;
}
DV_setDefaultFields(&workDV) ;
IV_setDefaultFields(&pivotsizesIV) ;
ALLOCATE(p_msg, struct _FactorMsg *, nfront+1) ;
for ( J = 0 ; J <= nfront ; J++ ) {
   p_msg[J] = NULL ;
}
recvsPosted = CVinit(nfront, 'N') ;
/*
   ---------------------------
   loop while a path is active
   ---------------------------
*/
*perror = -1 ;
firstmsgsent = NULL ;
while ( (J = Ideq_removeFromHead(dequeue)) != -1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n ### checking out front %d", J) ;
      fprintf(msgFile, "\n remaining fronts in dequeue") ;
      Ideq_writeForHumanEye(dequeue, msgFile) ;
      fflush(msgFile) ;
   }
   if ( recvsPosted[J] == 'N' && owners[J] == myid ) {
      MARKTIME(t1) ;
/*
      -----------------------------------
      post any receives for owned front J
      -----------------------------------
*/
      p_msg[J] = wakeupFront(frontmtx, J, myid, owners, firsttag,
                             chvmanager, aggList, stats, cpus,
                             comm, msglvl, msgFile) ;
      recvsPosted[J] = 'Y' ;
      MARKTIME(t2) ;
      cpus[10] += t2 - t1 ;
   }
   if ( p_msg[J] != NULL ) {
/*
      ---------------------------------------
      check for received messages for front J
      ---------------------------------------
*/
      MARKTIME(t1) ;
      p_msg[J] = checkMessages(frontmtx, J, p_msg[J], aggList, 
                               postList, chvmanager, firsttag, 
                               stats, comm, msglvl, msgFile) ;
      MARKTIME(t2) ;
      cpus[11] += t2 - t1 ;
   }
/*
   ----------------------------------
   check for numeric work for front J
   ----------------------------------
*/
   FrontMtx_factorVisit(frontmtx, pencil, J, myid, owners, fronts,
                        lookahead, tau, droptol, status, heads, 
                        &pivotsizesIV, &workDV, par, aggList, postList, 
                        chvmanager, stats, cpus, msglvl, msgFile) ;
   if ( status[J] == 'E' ) {
/*
      ---------------------------------------
      error detected while processing front J
      ---------------------------------------
*/
      *perror = J ;
      firstmsgsent = sendErrorMessages(nfront, J, nproc, myid, firsttag,
                          firstmsgsent, stats, comm, msglvl, msgFile) ;
      break ;
   } else if ( status[J] == 'F' ) {
#if MYDIAGNOSTIC > 0
{
      double   t1 ;
      MARKTIME(t1) ;
      fprintf(stdout, 
              "\n proc %d, time %8.3f : front %d is finished", 
              myid, t1, J) ;
      fflush(stdout) ;
      fprintf(msgFile, 
              "\n proc %d, time %8.3f : front %d is finished", 
              myid, t1, J) ;
      fflush(msgFile) ;
}
#endif
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n front %d is finished", J) ;
         fflush(msgFile) ;
      }
      if ( (K = par[J]) != -1 && --nactiveChild[K] == 0 ) {
/*
         ------------------------
         load parent into dequeue
         ------------------------
*/
         Ideq_insertAtHead(dequeue, K) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n\n 2. adding front %d to dequeue", K) ;
            fprintf(msgFile, "\n fronts in dequeue") ;
            Ideq_writeForHumanEye(dequeue, msgFile) ;
            fflush(msgFile) ;
         }
      }
/*
      ------------------------------------
      send any messages to other processes
      ------------------------------------
*/
      MARKTIME(t1) ;
      firstmsgsent = sendMessages(frontmtx, J, myid, nfront, par, 
                                  owners, firsttag, chvmanager, aggList,
                                  postList, firstmsgsent, stats,
                                  comm, msglvl, msgFile) ;
      MARKTIME(t2) ;
      cpus[12] += t2 - t1 ;
   } else {
/*
      ---------------------------------------------
      front is still active, place on tail of queue
      ---------------------------------------------
*/
      Ideq_insertAtTail(dequeue, J) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n front %d is not finished", J) ;
         fprintf(msgFile, "\n fronts in dequeue") ;
         Ideq_writeForHumanEye(dequeue, msgFile) ;
         fflush(msgFile) ;
      }
   }
   if ( firstmsgsent != NULL ) {
/*
      ------------------------------------------------------
      check the list of sent messages to recycle Chv objects
      ------------------------------------------------------
*/
      MARKTIME(t1) ;
      firstmsgsent = checkSentMessages(firstmsgsent, chvmanager,
                                       msglvl, msgFile) ;
      MARKTIME(t2) ;
      cpus[13] += t2 - t1 ;
   }
/* 
   --------------------------------------
   probe for an error termination message
   --------------------------------------
*/
   MPI_Iprobe(MPI_ANY_SOURCE, ERROR_TAG, comm, &flag, &mpistatus) ;
   if ( flag == 1 ) {
/*
      ---------------------------------------------------------
      an error message has been received, break out of the loop
      ---------------------------------------------------------
*/
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n ERROR MESSAGE RECEIVED FROM %d",
                 mpistatus.MPI_SOURCE) ;
         fflush(msgFile) ;
      }
      break ;
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n done with main loop") ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------
   now do an allgather on the termination flags
   --------------------------------------------
*/
term_flags = IVinit(nproc, -1) ;
MPI_Allgather((void *) perror, 1, MPI_INT, 
              (void *) term_flags, 1, MPI_INT, comm) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n allgather on error code complete") ;
   IVfprintf(msgFile, nproc, term_flags) ;
   fflush(msgFile) ;
}
if ( *perror < 0 ) {
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      if ( term_flags[iproc] >= 0 ) {
         *perror = term_flags[iproc] ;
         break ;
      }
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n error code = %d", *perror) ;
   fflush(msgFile) ;
}
rootchv = NULL ;
if ( *perror >= 0 ) {
/*
   ----------------------------------------
   error termination has occured,
   link all messages together into one list
   cancel all pending messages
   ----------------------------------------
*/
   head = firstmsgsent ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (msg = p_msg[J]) != NULL ) {
         while ( msg->next != NULL ) {
            msg = msg->next ;
         }
         msg->next = head ;
         head = p_msg[J] ;
      }
   }
   cancelAndFreeMessages(head, chvmanager, msglvl, msgFile) ;
} else {
/*
   -------------------------------------------
   normal termination, wait until all messages 
   have been received and free their storage
   -------------------------------------------
*/
   while ( firstmsgsent != NULL ) {
/*
      -------------------------------
      check the list of sent messages
      -------------------------------
*/
      MARKTIME(t1) ;
      firstmsgsent = checkSentMessages(firstmsgsent, chvmanager,
                                       msglvl, msgFile) ;
      MARKTIME(t2) ;
      cpus[13] += t2 - t1 ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n all sent messages have been free'd") ;
      fflush(msgFile) ;
   }
   if ( postList != NULL ) {
/*
      -----------------------------------------
      set the root chevron of the factorization
      -----------------------------------------
*/
      rootchv = ChvList_getList(postList, nfront) ;
   }
}
/*
   --------------------
   load some statistics
   --------------------
*/
stats[3]  = frontmtx->nentD ;
stats[4]  = frontmtx->nentL ;
stats[5]  = frontmtx->nentU ;
stats[14] = chvmanager->nactive ;
stats[15] = chvmanager->nbytesactive ;
stats[16] = chvmanager->nbytesrequested ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( postList != NULL ) {
   ChvList_free(postList) ;
}
ChvList_free(aggList) ;
CVfree(status) ;
IVfree(nactiveChild) ;
IP_free(heads[nfront+1]) ;
FREE(heads) ;
FREE(fronts) ;
Ideq_free(dequeue) ;
DV_clearData(&workDV) ;
IV_clearData(&pivotsizesIV) ;
CVfree(recvsPosted) ;
FREE(p_msg) ;
IVfree(term_flags) ;

return(rootchv) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------
   initialize a new instance
 
   created -- 97nov13, cca
   -------------------------
*/
static FactorMsg *
FactorMsg_new (
   void
) {
FactorMsg   *msg ;

ALLOCATE(msg, struct _FactorMsg, 1) ;
FactorMsg_setDefaultFields(msg) ;

return(msg) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields
 
   created -- 97nov13, cca
   -----------------------
*/
static void
FactorMsg_setDefaultFields (
   FactorMsg   *msg
) {
msg->info[0] =   0  ;
msg->info[1] =   0  ;
msg->info[2] =   0  ;
msg->base    = NULL ;
msg->chv     = NULL ;
msg->req     = NULL ;
msg->next    = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data
 
   created -- 97nov13, cca
   -----------------------
*/
static void
FactorMsg_clearData (
   FactorMsg   *msg
) {
FactorMsg_setDefaultFields(msg) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   free the object
 
   created -- 97nov13, cca
   -----------------------
*/
static void
FactorMsg_free (
   FactorMsg   *msg
) {
FactorMsg_clearData(msg) ;
FREE(msg) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   when a front is initially active, 
   post its necessary receive messages

   created -- 97nov13, cca
   -----------------------------------
*/
static FactorMsg *
wakeupFront ( 
   FrontMtx     *frontmtx,
   int          J,
   int          myid,
   int          frontOwners[],
   int          firsttag,
   ChvManager   *chvmanager,
   ChvList      *aggList,
   int          stats[],
   double       cpus[],
   MPI_Comm     comm,
   int          msglvl,
   FILE         *msgFile
) {
Chv         *chv ;
FactorMsg   *first, *msg ;
int         aggcount, ii, I, nbytes, nD, nfront, nL, nU, tag ;
int         *fch, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || J < 0 || J > (nfront = frontmtx->nfront)
   || chvmanager == NULL || aggList == NULL
   || (msglvl > 1 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in wakeupFront()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside wakeupFront, J = %d", J) ;
   fflush(msgFile) ;
}
first    = NULL ;
aggcount = aggList->counts[J] ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n aggcount for %d = %d", J, aggcount) ;
   fflush(msgFile) ;
}
if ( aggcount > 0 ) {
/*
   -----------------------------------------------
   post receives for the incoming aggregate fronts
   ----------------------------------------------
*/
   tag = AGGREGATE_CHV_TAG(J) ;
   for ( ii = 0 ; ii < aggcount ; ii++ ) {
      chv = FrontMtx_setupFront(frontmtx, NULL, J, myid, frontOwners,
                                chvmanager, cpus, msglvl, msgFile) ;
      msg          = FactorMsg_new() ;
      msg->info[0] = AGGREGATE_CHV ;
      msg->info[1] =   J    ;
      msg->info[2] = nbytes = Chv_nbytesInWorkspace(chv) ;
      msg->chv     = chv    ;
      msg->base    = (void *) Chv_workspace(chv) ;
      msg->next    = first  ;
      first        = msg    ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
          "\n 1. posting Irecv, msg %p, type = %d, J = %d, tag = %d",
           msg, AGGREGATE_CHV, J, tag) ;
         fflush(msgFile) ;
      }
      MPI_Irecv(msg->base, nbytes, MPI_BYTE, MPI_ANY_SOURCE,
                tag, comm, &msg->req) ;
      FrontMtx_initialFrontDimensions(frontmtx, J, 
                                      &nD, &nL, &nU, &nbytes) ;
      stats[8]++ ;
      stats[9] += nbytes ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, ", return") ;
         fflush(msgFile) ;
      }
   }
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) && frontOwners[J] == myid ) {
/*
   ---------------------------------------------
   post receives for child notification messages
   ---------------------------------------------
*/
   fch = frontmtx->tree->fch ;
   if ( fch[J] != -1 ) {
      sib  = frontmtx->tree->sib ;
      for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
         if ( frontOwners[I] != myid ) {
            msg          = FactorMsg_new() ;
            msg->info[0] = POSTPONED_NOTIFICATION ;
            msg->info[1] =   I   ;
            msg->info[2] =   0   ;
            msg->next    = first ;
            first        = msg   ;
            stats[12]++ ;
            stats[13] += 3*sizeof(int) ;
            tag = POSTPONED_NOTIFICATION_TAG(I) ;
            if ( msglvl > 1 ) {
            fprintf(msgFile, 
             "\n 2. posting Irecv, msg %p, type = %d, I = %d, tag = %d",
              msg, POSTPONED_NOTIFICATION, I, tag) ;
               fflush(msgFile) ;
            }
            MPI_Irecv((void *) msg->info, 3, MPI_INT, frontOwners[I], 
                      tag, comm, &msg->req) ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, ", return") ;
               fflush(msgFile) ;
            }
         }
      }
   }
}
return(first) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   check for completion of messages for front J

   created -- 97nov13, cca
   --------------------------------------------
*/
static FactorMsg *
checkMessages (
   FrontMtx     *frontmtx,
   int          J,
   FactorMsg    *firstmsg,
   ChvList      *aggList,
   ChvList      *postList,
   ChvManager   *chvmanager,
   int          firsttag,
   int          stats[],
   MPI_Comm     comm,
   int          msglvl,
   FILE         *msgFile
) {
FactorMsg    *msg, *nextmsg ;
int          flag, nbytes, nfront, source, tag, type ;
MPI_Status   status ;
/*
   ----------------------------------
   loop over the messages in the list
   ----------------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside checkMessages, J = %d, firstmsg = %p", 
           J, firstmsg) ;
   fflush(msgFile) ;
}
nfront = frontmtx->nfront ;
for ( msg = firstmsg, firstmsg = NULL ; 
      msg != NULL ; 
      msg = nextmsg ) {
/*
   ---------------------------
   set next message to look at
   ---------------------------
*/
   nextmsg   = msg->next ;
   msg->next = NULL ;
/*
   --------------------------------------
   extract the type, front id and # bytes
   --------------------------------------
*/
   type   = msg->info[0] ;
   nbytes = msg->info[2] ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, 
              "\n checking msg %p : type %d, front %d, nbytes %d",
              msg, type, msg->info[1], nbytes) ;
      fflush(msgFile) ;
   }
   MPI_Test(&msg->req, &flag, &status) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, ", flag = %d", flag) ;
      fflush(msgFile) ;
   }
   if ( flag != 1 ) {
/*
      -----------------------------------------------
      message has not yet been received, keep in list
      -----------------------------------------------
*/
      msg->next = firstmsg ;
      firstmsg  = msg ;
   } else {
/*
      -------------------------
      message has been received
      -------------------------
*/
      source = status.MPI_SOURCE ;
      tag    = status.MPI_TAG ;
      nbytes = msg->info[2] ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
                 "\n    message received: source %d, tag %d, nbytes %d",
                 source, tag, nbytes) ;
         fflush(msgFile) ;
      }
      switch ( type ) {
      case AGGREGATE_CHV :
/*
         -----------------
         aggregate chevron
         -----------------
*/
         Chv_initFromBuffer(msg->chv) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n    received aggregate front %d", J) ;
            Chv_writeForHumanEye(msg->chv, msgFile) ;
            fflush(msgFile) ;
         }
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n    adding chv to aggregate list"
                    "\n    before count is %d", aggList->counts[J]) ;
            fflush(msgFile) ;
         }
         ChvList_addObjectToList(aggList, msg->chv, J) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile,
                    "\n    count is now %d", aggList->counts[J]) ;
            fflush(msgFile) ;
         }
         FactorMsg_free(msg) ;
         break ;
      case POSTPONED_NOTIFICATION :
         if ( nbytes > 0 ) {
/*
            --------------------------------------------------
            post receive for message to contain postponed data
            --------------------------------------------------
*/
            msg->info[0] = POSTPONED_CHV ;
            msg->chv     = ChvManager_newObjectOfSizeNbytes(chvmanager,
                                                             nbytes) ;
            msg->base    = (void *) Chv_workspace(msg->chv) ;
            msg->next    = nextmsg ;
            nextmsg      = msg ;
            tag          += nfront + 1 ;
            stats[12]++ ;
            stats[13] += nbytes ;
            if ( msglvl > 1 ) {
               fprintf(msgFile,
          "\n    2. posting Irecv, msg %p, type = %d, I = %d, tag = %d",
          msg, POSTPONED_NOTIFICATION, msg->info[1], tag) ;
               fflush(msgFile) ;
            }
            MPI_Irecv(msg->base, nbytes, MPI_BYTE,
                      source, tag, comm, &msg->req) ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, ", return") ;
               fflush(msgFile) ;
            }
         } else {
/*
            ------------------------
            postponed data is empty
            ------------------------
*/
            if ( msglvl > 1 ) {
               fprintf(msgFile,
                     "\n    adding NULL object to postponed list"
                     "\n    before, count is %d", postList->counts[J]) ;
               fflush(msgFile) ;
            }
            ChvList_addObjectToList(postList, NULL, J) ;
            if ( msglvl > 1 ) {
               fprintf(msgFile,
                       "\n   count is now %d", postList->counts[J]) ;
               fflush(msgFile) ;
            }
            FactorMsg_free(msg) ;
         }
         break ;
      case POSTPONED_CHV :
/*
         ------------------------------
         put postponed data on its list
         ------------------------------
*/
         Chv_initFromBuffer(msg->chv) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n    received postponed data %d", J) ;
            Chv_writeForHumanEye(msg->chv, msgFile) ;
            fflush(msgFile) ;
         }
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                     "\n    before, count is %d", postList->counts[J]) ;
            fflush(msgFile) ;
         }
         ChvList_addObjectToList(postList, msg->chv, J) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile,
                    "\n    count is now %d", postList->counts[J]) ;
            fflush(msgFile) ;
         }
         FactorMsg_free(msg) ;
         break ;
      default :
         break ;
      }
   }
}
return(firstmsg) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   after front J has completed, remove relevant data from the 
   aggregage and/or postponed lists and send to the owning processor.

   created -- 97jul10,  cca
   ------------------------------------------------------------------
*/
static FactorMsg *
sendMessages (
   FrontMtx     *frontmtx,
   int          J,
   int          myid,
   int          nfront,
   int          par[],
   int          owners[],
   int          firsttag,
   ChvManager   *chvmanager,
   ChvList      *aggList,
   ChvList      *postList,
   FactorMsg    *firstmsgsent,
   int          stats[],
   MPI_Comm     comm,
   int          msglvl,
   FILE         *msgFile
) {
Chv        *chv ;
FactorMsg   *msg ;
int         destination, K, nbytes, tag ;

destination = owners[J] ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
      "\n\n inside sendMessages, J = %d, firsttag = %d, destination %d",
      J, firsttag, destination) ;
   fflush(msgFile) ;
}
if ( destination != myid ) {
/*
   ------------------------------------------
   this process does not own the front, so
   an aggregate front will have been created.
   ------------------------------------------
*/
   chv = ChvList_getList(aggList, J) ;
   if ( chv != NULL ) {
/*
      ------------------------------------
      an aggregate front has been created,
      send to the owning processor
      ------------------------------------
*/
      nbytes = Chv_nbytesNeeded(chv->nD, chv->nL, chv->nU, 
                                frontmtx->type, chv->symflag) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
                 "\n proc %d does not own front %d, nbytes = %d", 
                 myid, J, nbytes) ;
         fflush(msgFile) ;
      }
      msg = FactorMsg_new() ;
      msg->info[0] = AGGREGATE_CHV ;
      msg->info[1] = J ;
      msg->info[2] = nbytes ;
      msg->chv     = chv ;
      msg->base    = (void *) Chv_workspace(chv) ;
      tag          = AGGREGATE_CHV_TAG(J) ;
      msg->next    = firstmsgsent ;
      firstmsgsent = msg  ;
      stats[6]++ ;
      stats[7] += nbytes ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
  "\n posting Isend, msg %p, type %d, J %d, nbytes %d, dest %d, tag %d",
  msg, msg->info[0], msg->info[1], msg->info[2], destination, tag) ;
         fprintf(msgFile, 
                 "\n    wrk %p, nD %d, nL %d, nU %d, symflag %d",
                 Chv_workspace(chv), chv->nD, chv->nL, 
                 chv->nU, chv->symflag) ;
         fflush(msgFile) ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n aggregate front %d, chv %p", J, chv) ;
         if ( chv != NULL ) {
            Chv_writeForHumanEye(msg->chv, msgFile) ;
         }
         fflush(msgFile) ;
      }
      MPI_Isend(msg->base, nbytes, MPI_BYTE, destination, tag, 
                comm, &msg->req) ;
   }
} else if ( postList != NULL ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n    checking for postponed messages") ;
      fflush(msgFile) ;
   }
/*
   ----------------------------------------------
   pivoting enabled, look for a postponed chevron
   ----------------------------------------------
*/
   K = par[J] ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n    J = %d, K = %d", J, K) ;
      fflush(msgFile) ;
   }
   if ( K != -1 && (destination = owners[K]) != myid ) {
      if ( msglvl > 1 ) {
         fprintf(msgFile, ", destination = %d", destination) ;
         fflush(msgFile) ;
      }
/*
      ---------------------------------------
      owner of the parent is not this process
      ---------------------------------------
*/
      if ( ChvList_isListNonempty(postList, K) == 1 ) {
/*
         ----------------------------------------
         there is a postponed chevron on the list
         ----------------------------------------
*/
         chv    = ChvList_getList(postList, K) ;
         nbytes = Chv_nbytesNeeded(chv->nD, chv->nL, chv->nU, 
                                   frontmtx->type, chv->symflag) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                    "\n    nD = %d, nL = %d, nU = %d, nbytes = %d",
                    chv->nD, chv->nL, chv->nU, nbytes) ;
            fflush(msgFile) ; 
         }
      } else {
         chv    = NULL ;
         nbytes = 0 ;
      }
/*
      ----------------------------------------
      owner of the parent is not this process,
      send a notification message
      ----------------------------------------
*/
      msg = FactorMsg_new() ;
      msg->info[0] = POSTPONED_NOTIFICATION ;
      msg->info[1] = J ;
      msg->info[2] = nbytes ;
      tag          = POSTPONED_NOTIFICATION_TAG(J) ;
      msg->next    = firstmsgsent ;
      firstmsgsent = msg  ;
      stats[10]++ ;
      stats[11] += 3*sizeof(int) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
  "\n posting Isend, msg %p, type %d, J %d, nbytes %d, dest %d, tag %d",
  msg, msg->info[0], msg->info[1], msg->info[2], destination, tag) ;
         fflush(msgFile) ;
      }
      MPI_Isend((void *) msg->info, 3, MPI_INT, destination, tag, 
                comm, &msg->req) ;
      if ( nbytes > 0 ) {
/*
         ----------------------------------------------------------
         there is postponed data, send the data in a second message
         ----------------------------------------------------------
*/
         msg = FactorMsg_new() ;
         msg->info[0] = POSTPONED_CHV ;
         msg->info[1] = J ;
         msg->info[2] = nbytes ;
         msg->chv     = chv ;
         msg->base    = (void *) Chv_workspace(chv) ;
         tag          = POSTPONED_CHV_TAG(J) ;
         msg->next    = firstmsgsent ;
         firstmsgsent = msg  ;
         stats[10]++ ;
         stats[11] += nbytes ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
  "\n posting Isend, msg %p, type %d, J %d, nbytes %d, dest %d, tag %d",
  msg, msg->info[0], msg->info[1], msg->info[2], destination, tag) ;
         fprintf(msgFile, 
                 "\n    wrk %p, nD %d, nL %d, nU %d, symflag %d",
                 Chv_workspace(chv), chv->nD, chv->nL, 
                 chv->nU, chv->symflag) ;
            fflush(msgFile) ;
         }
         MPI_Isend(msg->base, nbytes, MPI_BYTE, destination, 
                   tag, comm, &msg->req) ;
      }
   }
}
return(firstmsgsent) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   check the list of posted sends.
   if a send has completed, free the message structure 
   and return the Chv object (if present) to its manager.

   created -- 97nov13, cca
   -------------------------------------------------------
*/
static FactorMsg *
checkSentMessages (
   FactorMsg     *firstmsgsent,
   ChvManager   *chvmanager,
   int           msglvl,
   FILE          *msgFile
) {
FactorMsg    *msg, *nextmsg ;
int          flag, J, nbytes, type ;
MPI_Status   status ;

for ( msg = firstmsgsent, firstmsgsent = NULL ;
      msg != NULL ;
      msg = nextmsg ) {
   nextmsg = msg->next ;
   type    = msg->info[0] ;
   J       = msg->info[1] ;
   nbytes  = msg->info[2] ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, 
              "\n checking sent msg %p : type %d, front %d, nbytes %d",
              msg, type, J, nbytes) ;
      fflush(msgFile) ;
   }
   MPI_Test(&msg->req, &flag, &status) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, ", flag = %d", flag) ;
      fflush(msgFile) ;
   }
   if ( flag == 1 ) {
      if ( msg->chv != NULL ) {
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n    chv %p released", msg->chv) ;
            fflush(msgFile) ;
         }
         ChvManager_releaseObject(chvmanager, msg->chv) ;
      }
      FactorMsg_free(msg) ;
   } else {
      msg->next    = firstmsgsent ;
      firstmsgsent = msg ;
   }
}
return(firstmsgsent) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- an error has been detected at front J, 
              send a message to all other processors

   created -- 98aug29, cca
   -------------------------------------------------
*/
static FactorMsg *
sendErrorMessages (
   int         nfront,
   int         J,
   int         nproc,
   int         myid,
   int         firsttag,
   FactorMsg   *firstmsgsent,
   int         stats[],
   MPI_Comm    comm,
   int         msglvl,
   FILE        *msgFile
) {
int         iproc, tag ;
FactorMsg   *msg ;

for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( iproc != myid ) {
      msg = FactorMsg_new() ;
      msg->info[0] = tag = ERROR_TAG ;
      msg->info[1] = J ;
      msg->info[2] = 0 ;
      msg->chv     = NULL ;
      msg->base    = (void *) msg->info ;
      msg->next    = firstmsgsent ;
      firstmsgsent = msg ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
         "\n posting error Isend, info = {%d, %d, %d}, dest %d, tag %d",
        msg->info[0], msg->info[1], msg->info[2], iproc, tag) ;
         fflush(msgFile) ;
      }
      MPI_Isend(msg->base, 2, MPI_INT, iproc, tag, comm, &msg->req) ;
   }
}
return(firstmsgsent) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- when an error has been detected,  
              cancel all outstanding sends and 
              receives and free their storage
  
   created -- 98aug29, cca
   -------------------------------------------
*/
static void
cancelAndFreeMessages (
   FactorMsg    *head,
   ChvManager   *chvmanager,
   int          msglvl,
   FILE         *msgFile
) {
FactorMsg    *msg, *nextmsg ;
int          flag ;
MPI_Status   status ;
/*
   ------------------------------------------------
   cancel all the messages and test for termination
   ------------------------------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside cancelAndFreeMessages()") ;
   fflush(msgFile) ;
}
for ( msg = head, head = NULL ; msg != NULL ; msg = nextmsg ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n msg %p, info = {%d,%d,%d}",
              msg, msg->info[0], msg->info[1], msg->info[2]) ;
      fflush(msgFile) ;
   }
   nextmsg = msg->next ;
   MPI_Cancel(&msg->req) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, ", canceled") ;
      fflush(msgFile) ;
   }
   MPI_Test(&msg->req, &flag, &status) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, ", tested, flag = %d", flag) ;
      fflush(msgFile) ;
   }
   if ( flag == 1 ) {
/*
      ----------------------
      message has terminated
      ----------------------
*/
      if ( msg->chv != NULL ) {
         if ( msglvl > 1 ) {
            fprintf(msgFile, ", releasing chv %p", msg->chv) ;
            fflush(msgFile) ;
         }
         ChvManager_releaseObject(chvmanager, msg->chv) ;
      }
      FactorMsg_free(msg) ;
   } else {
      msg->next = head ;
      head = msg ;
   }
}
/*
   --------------------------------------------
   loop until all messages have been terminated
   --------------------------------------------
*/
while ( head != NULL ) {
   for ( msg = head, head = NULL ; msg != NULL ; msg = nextmsg ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n msg %p, info = {%d,%d,%d}",
              msg, msg->info[0], msg->info[1], msg->info[2]) ;
      fflush(msgFile) ;
   }
      nextmsg = msg->next ;
      MPI_Test(&msg->req, &flag, &status) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, ", tested, flag = %d", flag) ;
         fflush(msgFile) ;
      }
      if ( flag == 1 ) {
/*
         ----------------------
         message has terminated
         ----------------------
*/
         if ( msg->chv != NULL ) {
            if ( msglvl > 1 ) {
               fprintf(msgFile, ", releasing chv %p", msg->chv) ;
               fflush(msgFile) ;
            }
            ChvManager_releaseObject(chvmanager, msg->chv) ;
         }
         FactorMsg_free(msg) ;
      } else {
         msg->next = head ;
         head = msg ;
      }
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n leaving cancelAndFreeMessages()") ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
