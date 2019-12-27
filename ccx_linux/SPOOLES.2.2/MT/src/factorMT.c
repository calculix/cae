/*  factorMT.c  */

#include "../spoolesMT.h"
#include "../../Ideq.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   worker method for each thread
   -----------------------------
*/
static void * FrontMtx_workerFactor ( void *arg ) ;
/*
   ----------------------------------------------------------------
  this object is used during a threaded factorization
 
   pencil    -- matrix pencil A + sigma * B
   tau       -- upper bound on factor entries when pivoting enabled
   droptol   -- drop tolerance used for sparse fronts
   ownersIV  -- map from fronts to threads
   lookahead -- parameter used to control computation lookahead
      0   --> no lookahead
      > 0 --> look up a number of levels up the tree
 
   frontmtx  -- object used to store factorization
   manager   -- object used to manage working Chv objects
   aggList   -- object used to store aggregate data
   postList  -- object used to store postponed data
   perror    -- pointer to external error flag
 
   heads   -- IP vector to maintain update lists
   myid    -- thread id
   fronts  -- vector of pointers to active fronts
   stats   -- int statistics vector
   cpus    -- double vector to store breakdown of cpu times
   msglvl  -- message level
   msgFile -- message file
 
   created -- 97may30, cca
   ----------------------------------------------------------------
*/
typedef struct _FactorData   FactorData ;
struct _FactorData {
/*
   -------------------------
   global data, not modified
   -------------------------
*/
   Pencil      *pencil   ;
   double      tau       ;
   double      droptol   ;
   IV          *ownersIV ;
   int         lookahead ;
/*
   ---------------------
   shared data, modified
   ---------------------
*/
   FrontMtx    *frontmtx   ;
   ChvManager  *chvmanager ;
   ChvList     *aggList    ;
   ChvList     *postList   ;
   int         *perror     ;
/*
   ----------
   local data
   ----------
*/
   int          myid      ;
   int          stats[10] ;
   double       cpus[20]  ;
   int          msglvl    ;
   FILE         *msgFile  ;
} ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   parallel factorization method for A.
   all but two input parameters are the same as the serial method.
   this is a wrapper method around FrontMtx_MT_factorInpMtx().

   ownersIV  -- pointer to IV object that holds map from fronts 
                to threads
   lookahead -- lookahead parameter that allows computation at
                higher levels of the front tree to proceed when
                lower fronts are not yet finish. use lookahead = 0
                to turn off this option. otherwise lookahead ancestors
                of an active unfinished front can be active.

   return value -- pointer to the first Chv object in a list
                   that contains postponed data
 
   created -- 98may16, cca
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_MT_factorInpMtx (
   FrontMtx     *frontmtx,
   InpMtx       *inpmtx,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   IV           *ownersIV,
   int          lookahead,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) {
Chv      *rootchv ;
double   sigma[2] = {0.0, 0.0} ;
Pencil   pencil ;

Pencil_setDefaultFields(&pencil) ;
Pencil_init(&pencil, frontmtx->type, frontmtx->symmetryflag,
            inpmtx, sigma, NULL) ;
rootchv = FrontMtx_MT_factorPencil(frontmtx, &pencil, tau, droptol,
                                 chvmanager, ownersIV, lookahead, 
                                 perror, cpus, stats, msglvl, msgFile) ;

return(rootchv) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   parallel factorization method for A + sigma*B.
   all but two input parameters are the same 
   as the FrontMtx_factorPencil method.

   ownersIV  -- pointer to IV object that holds map from fronts 
                to threads
   lookahead -- lookahead parameter that allows computation at
                higher levels of the front tree to proceed when
                lower fronts are not yet finish. use lookahead = 0
                to turn off this option. otherwise lookahead ancestors
                of an active unfinished front can be active.

   return value -- pointer to the first Chv object in a list
                   that contains postponed data
 
   created -- 98may16, cca
   -------------------------------------------------------------------
*/
Chv *
FrontMtx_MT_factorPencil (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   double       tau,
   double       droptol,
   ChvManager   *chvmanager,
   IV           *ownersIV,
   int          lookahead,
   int          *perror,
   double       cpus[],
   int          stats[],
   int          msglvl,
   FILE         *msgFile
) {
char         buffer[20] ;
Chv          *rootchv ;
ChvList      *aggList ;
ChvList      *postList ;
double       t0, t1, t2 ;
FactorData   *data, *dataObjects ;
FILE         *fp ;
int          ierr, ii, myid, nfront, nthread, rc ;
int          *owners ;
/*
   --------------
   check the data
   --------------
*/
MARKTIME(t0) ;
if (  frontmtx == NULL || pencil == NULL || tau < 1.0 || droptol < 0.0
   || ownersIV == NULL || lookahead < 0 || cpus == NULL || stats == NULL
   || msglvl < 0 || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MT_factorPencil()"
           "\n frontmtx = %p, pencil = %p"
           "\n tau = %f, droptol = %f, ownersIV = %p, lookahead = %d"
           "\n cpus = %p, stats = %p, msglvl = %d, msgFile = %p"
           "\n bad input\n", frontmtx, pencil, tau, droptol, 
           ownersIV, lookahead, cpus, stats, msglvl, msgFile) ;
   exit(-1) ;
}
IV_sizeAndEntries(ownersIV, &nfront, &owners) ;
nthread = 1 + IV_max(ownersIV) ;
/*
   --------------------------------------------------------------------
   create :
   (1) a ChvList object to handle the lists of aggregate Chv objects,
   (2) if pivoting is enabled, a ChvList object to handle 
       the lists of postponed Chv objects.
   --------------------------------------------------------------------
*/
MARKTIME(t1) ;
aggList = FrontMtx_aggregateList(frontmtx, ownersIV, 1) ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   postList = FrontMtx_postList(frontmtx, ownersIV, 1) ;
} else {
   postList = NULL ;
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : initialize lists and manager", 
           t2 - t1) ;
}
/*
   ------------------
   set the error flag
   ------------------
*/
*perror = -1 ;
/*
   ----------------------------------------------------------
   create nthread FactorData objects and load with their data
   ----------------------------------------------------------
*/
MARKTIME(t1) ;
ALLOCATE(dataObjects, struct _FactorData, nthread) ;
for ( myid = 0, data = dataObjects ; myid < nthread ; myid++, data++ ) {
   data->pencil     = pencil     ;
   data->tau        = tau        ;
   data->droptol    = droptol    ;
   data->ownersIV   = ownersIV   ;
   data->lookahead  = lookahead  ;
   data->frontmtx   = frontmtx   ;
   data->chvmanager = chvmanager ;
   data->aggList    = aggList    ;
   data->postList   = postList   ;
   data->perror     = perror     ;
   data->myid       = myid       ;
   IVzero(10, data->stats) ;
   DVzero(20, data->cpus) ;
   data->msglvl = msglvl ;
   if ( msglvl > 0 ) {
      sprintf(buffer, "res.%d", myid) ;
      if ( (fp = fopen(buffer, "w")) == NULL ) {
         fprintf(stderr, "\n fatal error, unable to open file %s",
                 buffer) ;
         exit(-1) ;
      }
      data->msgFile = fp ;
   } else {
      data->msgFile = NULL ;
   }
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : initialize data objects", t2 - t1) ;
}
/*
   ---------------------------
   switch over the thread type
   ---------------------------
*/
#if THREAD_TYPE == TT_SOLARIS
/*
   ---------------
   solaris threads
   ---------------
*/
MARKTIME(t1) ;
thr_setconcurrency(nthread) ;
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : set concurrency time", t2 - t1) ;
}
MARKTIME(t1) ;
for ( myid = 0, data = dataObjects ; 
      myid < nthread - 1 ;
      myid++, data++ ) {
   rc = thr_create(NULL, 0, FrontMtx_workerFactor, data, 0, NULL) ;
   if ( rc != 0 ) {
      fprintf(stderr, 
              "\n fatal error, myid = %d, rc = %d from thr_create",
              myid, rc) ;
      exit(-1) ;
   }
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : thread creation time", t2 - t1) ;
}
MARKTIME(t1) ;
FrontMtx_workerFactor(data) ;
for ( myid = 0 ; myid < nthread - 1 ; myid++ ) {
   thr_join(0, 0, 0) ;
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : thread join time", t2 - t1) ;
}
#endif
#if THREAD_TYPE == TT_POSIX
/*
   -------------
   POSIX threads
   -------------
*/
{
pthread_t        *tids ;
pthread_attr_t   attr  ;
void             *status ;
/*
#####   NOTE: for SGI machines, this command must be present
#####         for the thread scheduling to be efficient.
#####         this is NOT a POSIX call, but SGI needs it anyway
pthread_setconcurrency(nthread) ;
*/
pthread_attr_init(&attr) ;
/*
pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM) ;
*/
pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS) ;
ALLOCATE(tids, pthread_t, nthread) ;
MARKTIME(t1) ;
for ( myid = 0, data = dataObjects ; myid < nthread ; myid++, data++ ) {
/*
   rc = pthread_create(&tids[myid], &attr, 
                       FrontMtx_workerFactor, data) ;
*/
   rc = pthread_create(&tids[myid], NULL, FrontMtx_workerFactor, data) ;
   if ( rc != 0 ) {
      fprintf(stderr, 
              "\n fatal error, myid = %d, rc = %d from pthread_create",
              myid, rc) ;
      exit(-1) ;
   } else if ( msglvl > 1 ) {
      fprintf(stderr, "\n thread %d created", tids[myid]) ;
   }
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : thread creation time", t2 - t1) ;
}
MARKTIME(t1) ;
for ( myid = 0 ; myid < nthread ; myid++ ) {
   pthread_join(tids[myid], &status) ;
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : thread join time", t2 - t1) ;
}
FREE(tids) ;
pthread_attr_destroy(&attr) ;
}
#endif
/*
   --------------------------------------------
   get a pointer to any postponed root chevrons
   --------------------------------------------
*/
if ( postList != NULL ) {
   rootchv = ChvList_getList(postList, nfront) ;
} else {
   rootchv = NULL ;
}
/*
   -------------------
   fill the statistics
   -------------------
*/
for ( myid = 0, data = dataObjects ; 
      myid < nthread ;
      myid++, data++ ) {
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n thread %d stats", myid) ;
      IVfp80(msgFile, 10, data->stats, 20, &ierr) ;
      fprintf(msgFile, "\n thread %d cpus", myid) ;
      DVfprintf(msgFile, 10, data->cpus) ;
   }
   for ( ii = 0 ; ii < 10 ; ii++ ) {
      stats[ii] += data->stats[ii] ;
   }
   for ( ii = 0 ; ii <= 10 ; ii++ ) {
      cpus[ii] += data->cpus[ii] ;
   }
}
stats[3] = frontmtx->nentD ;
stats[4] = frontmtx->nentL ;
stats[5] = frontmtx->nentU ;
stats[6] = frontmtx->nlocks ;
stats[7] = aggList->nlocks ;
if ( postList != NULL ) {
   stats[8] = postList->nlocks ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, 
           "\n\n factorization has finished" 
           "\n %d locks of the front matrix"
           "\n %d locks of the aggregate list",
           frontmtx->nlocks, aggList->nlocks) ;
   if ( postList != NULL ) {
      fprintf(msgFile, "\n %d locks of the aggregate list",
          aggList->nlocks) ;
   }
}
/*
   -------------
   free the data
   -------------
*/
MARKTIME(t1) ;
ChvList_free(aggList) ;
if ( postList != NULL ) {
   ChvList_free(postList) ;
}
FREE(dataObjects) ;
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : total time", t2 - t1) ;
}

return(rootchv) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- worker method to factor the matrix


   created -- 97may30, cca
   ---------------------------------------------
*/
static void *
FrontMtx_workerFactor (
   void   *arg
) {
char          *status ;
Chv           **fronts ;
ChvList       *aggList, *postList ;
ChvManager    *chvmanager ;
FactorData    *data ;
FrontMtx      *frontmtx ;
Pencil        *pencil ;
double        droptol, tau, t1, t2 ;
double        *cpus ;
DV            *workDV ;
ETree         *frontETree ;
FILE          *msgFile ;
Ideq          *deq ;
int           J, K, lookahead, msglvl, myid, nfront ;
int           *nactiveChild, *owners, *par, *stats ;
IP            **heads ;
IV            *ownersIV, pivotsizesIV ;
Tree          *tree ;
#if THREAD_TYPE == TT_POSIX
/*
fprintf(stdout, "\n ### inside workerFactor, pthread_self() = %d", 
        pthread_self()) ;
fflush(stdout) ;
*/
#endif
/*
   -------------------------------
   extract pointers and dimensions
   -------------------------------
*/
MARKTIME(t1) ;
data = (FactorData *) arg ;
pencil     = data->pencil   ;
tau        = data->tau      ;
droptol    = data->droptol  ;
msglvl     = data->msglvl   ;
msgFile    = data->msgFile  ;
myid       = data->myid     ;
frontmtx   = data->frontmtx ;
chvmanager = data->chvmanager ;
aggList    = data->aggList    ;
postList   = data->postList   ;
frontETree = frontmtx->frontETree ;
tree       = ETree_tree(frontETree) ;
nfront     = ETree_nfront(frontETree) ;
par        = ETree_par(frontETree) ;
lookahead  = data->lookahead ;
ownersIV   = data->ownersIV ;
owners     = IV_entries(ownersIV) ;
stats      = data->stats;
cpus       = data->cpus ;
#if THREAD_TYPE == TT_SOLARIS
if ( msglvl > 2 ) {
   fprintf(stdout, 
           "\n ### inside workerFactor, myid = %d, thr_self() = %d", 
           myid, thr_self()) ;
   fflush(stdout) ;
}
#endif
#if THREAD_TYPE == TT_POSIX
if ( msglvl > 2 ) {
   fprintf(stdout, "\n ### inside workerFactor, myid = %d" 
                   ", pthread_self() = %d", myid, pthread_self()) ;
   fflush(stdout) ;
}
#endif
/*
   ---------------------------------------------------
   initialize the heads[] vector for the owned updates
   ---------------------------------------------------
*/
heads = FrontMtx_factorSetup(frontmtx, ownersIV, myid,
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
deq = FrontMtx_setUpDequeue(frontmtx, IV_entries(ownersIV), myid, 
                             status, heads, 'W', 'F', msglvl, msgFile) ;
FrontMtx_loadActiveLeaves(frontmtx, status, 'W', deq) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n status") ;
   CVfprintf(msgFile, nfront, status) ;
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
IV_setDefaultFields(&pivotsizesIV) ;
/*
if ( FRONTMTX_IS_PIVOTING(frontmtx) 
   && (  FRONTMTX_IS_SYMMETRIC(frontmtx) 
      || FRONTMTX_IS_HERMITIAN(frontmtx) ) ) {
   pivotsizesIV = IV_new() ;
} else {
   fprintf(msgFile, "\n no pivotsizesIV needed") ;
   fflush(msgFile) ;
   pivotsizesIV = NULL ;
}
*/
workDV = DV_new() ;
/*
   ---------------------------
   loop while a path is active
   ---------------------------
*/
IVzero(10, stats) ;
while ( (J = Ideq_removeFromHead(deq)) != -1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n ### checking out front %d", J) ;
      fflush(msgFile) ;
   }
   FrontMtx_factorVisit(frontmtx, pencil, J, myid, owners, fronts,
      lookahead, data->tau, data->droptol, status, heads, &pivotsizesIV,
      workDV, par, data->aggList, data->postList, 
      data->chvmanager, stats, cpus, msglvl, msgFile) ;
   if ( status[J] == 'E' ) {
/*
      ----------------------------
      error encountered at front J
      ----------------------------
*/
      *(data->perror) = J ;
      break ;
   } else if ( status[J] == 'F' ) {
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n front %d is finished", J) ;
         fflush(msgFile) ;
      }
      if ( (K = par[J]) != -1 && --nactiveChild[K] == 0 ) {
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n\n adding front %d to dequeue", K) ;
            fflush(msgFile) ;
         }
         Ideq_insertAtHead(deq, K) ;
      }
   } else {
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n front %d not yet done", J) ;
         fflush(msgFile) ;
      }
      Ideq_insertAtTail(deq, J) ;
   }
   if ( *(data->perror) >= 0 ) {
/*
      ------------------------------------
      error encountered, break out of loop
      ------------------------------------
*/
      break ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_clearData(&pivotsizesIV) ;
if ( workDV != NULL ) {
   DV_free(workDV) ;
}
CVfree(status) ;
IVfree(nactiveChild) ;
IP_free(heads[nfront+1]) ;
FREE(heads) ;
FREE(fronts) ;
Ideq_free(deq) ;
MARKTIME(t2) ;
cpus[10] = t2 - t1 ;
cpus[9] = t2 - t1 - cpus[0] - cpus[1] - cpus[2] - cpus[3] 
        - cpus[4] - cpus[5] - cpus[6] - cpus[7] - cpus[8] ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
