/*  QRfactorMT.c  */

#include "../spoolesMT.h"
#include "../../Ideq.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   worker method for each thread
   -----------------------------
*/
static void * FrontMtx_QR_workerFactor ( void *arg ) ;
/*
   ----------------------------------------------------------------
   this object is used during a threaded factorization
 
   mtxA      -- InpMtx matrix A 
   rowsIVL   -- list[J] contains rows of A to be assembled into front J
   firstnz   -- firstnz[irow] contains first column 
                of first nonzero entry in row
   ownersIV  -- map from fronts to threads
 
   frontmtx   -- object used to store factorization
   chvmanager -- object used to manage working Chv objects
   updlist    -- object used to store update Chv objects
 
   myid      -- thread id
   facops    -- number of factor operations
   cpus      -- double vector to store breakdown of cpu times
   msglvl    -- message level
   msgFile   -- message file
 
   created -- 97may30, cca
   ----------------------------------------------------------------
*/
typedef struct _QR_factorData   QR_factorData ;
struct _QR_factorData {
/*
   -------------------------
   global data, not modified
   -------------------------
*/
   InpMtx      *mtxA     ;
   IVL         *rowsIVL  ;
   int         *firstnz  ;
   IV          *ownersIV ;
/*
   ---------------------
   shared data, modified
   ---------------------
*/
   FrontMtx    *frontmtx   ;
   ChvManager  *chvmanager ;
   ChvList     *updlist    ;
/*
   ----------
   local data
   ----------
*/
   int          myid      ;
   double       facops    ;
   double       cpus[7]   ;
   int          msglvl    ;
   FILE         *msgFile  ;
} ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- compute a QR factorization using multiple threads

   created -- 98may29, cca
   ------------------------------------------------------------
*/
void
FrontMtx_MT_QR_factor (
   FrontMtx     *frontmtx,
   InpMtx       *mtxA,
   ChvManager   *chvmanager,
   IV           *ownersIV,
   double       cpus[],
   double       *pfacops,
   int          msglvl,
   FILE         *msgFile
) {
ChvList         *updlist ;
double          t0, t1 ;
IVL             *rowsIVL ;
int             ithread, myid, nthread, rc ;
int             *firstnz ;
QR_factorData   *data, *dataObjects ;

/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || mtxA == NULL || chvmanager == NULL 
   || ownersIV == NULL || cpus == NULL || pfacops == NULL 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MT_QR_factor()"
           "\n bad input\n") ;
   exit(-1) ;
}
nthread = 1 + IV_max(ownersIV) ;
/*
   ----------------------------------------------------------------
   create the update Chv list object
   create the rowsIVL object, where
      list(J) = list of rows that are assembled in front J
   firstnz[irowA] = first column with nonzero element in A(irowA,*)
   ----------------------------------------------------------------
*/
MARKTIME(t0) ;
updlist = FrontMtx_postList(frontmtx, ownersIV, LOCK_IN_PROCESS) ;
FrontMtx_QR_setup(frontmtx, mtxA, &rowsIVL, &firstnz, msglvl, msgFile) ;
MARKTIME(t1) ;
cpus[0] = t1 - t0 ;
/*
   ------------------------------------
   create and load nthread data objects
   ------------------------------------
*/
ALLOCATE(dataObjects, struct _QR_factorData, nthread) ;
for ( myid = 0, data = dataObjects ; myid < nthread ; myid++, data++ ) {
   data->mtxA       = mtxA       ;
   data->rowsIVL    = rowsIVL    ;
   data->firstnz    = firstnz    ;
   data->ownersIV   = ownersIV   ;
   data->frontmtx   = frontmtx   ;
   data->chvmanager = chvmanager ;
   data->updlist    = updlist    ;
   data->myid       = myid       ;
   DVzero(7, data->cpus) ;
   data->facops = 0.0 ;
   data->msglvl  = msglvl ;
   if ( msglvl > 0 ) {
      char   buffer[20] ;
      sprintf(buffer, "res.%d", myid) ;
      if ( (data->msgFile = fopen(buffer, "w")) == NULL ) {
         fprintf(stderr, "\n fatal error in FrontMtx_MT_QR_factor()"
                 "\n unable to open file %s", buffer) ;
         exit(-1) ;
      }
   } else {
      data->msgFile = NULL ;
   }
}
#if THREAD_TYPE == TT_SOLARIS
/*
   ----------------------------------
   Solaris threads.
   (1) set the concurrency
   (2) create nthread - 1 new threads
   (3) execute own thread
   (4) join the threads
   ----------------------------------
*/
thr_setconcurrency(nthread) ;
for ( myid = 0, data = dataObjects ; 
      myid < nthread - 1 ; 
      myid++, data++ ) {
   rc = thr_create(NULL, 0, FrontMtx_QR_workerFactor, data, 0, NULL) ;
   if ( rc != 0 ) {
      fprintf(stderr, 
              "\n fatal error, myid = %d, rc = %d from thr_create()", 
              myid, rc) ;
      exit(-1) ;
   }
}
FrontMtx_QR_workerFactor(data) ;
for ( myid = 0 ; myid < nthread - 1 ; myid++ ) {
   thr_join(0, 0, 0) ;
}
#endif
#if THREAD_TYPE == TT_POSIX
/*
   ----------------------------------
   POSIX threads.
   (1) if SGI, set the concurrency
   (2) create nthread new threads
   (3) join the threads
   ----------------------------------
*/
{
pthread_t        *tids ;
pthread_attr_t   attr  ;
void             *status ;
/*
   ---------------------------------------------------------
   #### NOTE: for SGI machines, this command must be present
   ####       for the thread scheduling to be efficient.
   ####       this is NOT a POSIX call, but necessary
   ---------------------------------------------------------
   pthread_setconcurrency(nthread) ;
*/
   pthread_attr_init(&attr) ;
/*
   pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM) ;
*/
   pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS) ;
   ALLOCATE(tids, pthread_t, nthread) ;
   for ( myid = 0 ; myid < nthread ; myid++ ) {
      tids[myid] = 0 ;
   }
   for ( myid = 0, data = dataObjects ; 
         myid < nthread ; 
         myid++, data++ ) {
      rc = pthread_create(&tids[myid], &attr,
                          FrontMtx_QR_workerFactor, data) ;
      if ( rc != 0 ) {
         fprintf(stderr,
                 "\n fatal error in FrontMtx_MT_QR_factor()"
                 "\n myid = %d, rc = %d from pthread_create()",
                 myid, rc) ;
         exit(-1) ;
      } else if ( msglvl > 2 ) {
         fprintf(stderr, "\n thread %d created", tids[myid]) ;
      }
   }
   for ( myid = 0 ; myid < nthread ; myid++ ) {
      pthread_join(tids[myid], &status) ;
   }
   FREE(tids) ;
   pthread_attr_destroy(&attr) ;
}
#endif
/*
   ----------------------------------------------
   fill the cpu vector and factor operation count
   ----------------------------------------------
*/
*pfacops = 0 ;
for ( myid = 0, data = dataObjects ; myid < nthread ; myid++, data++ ) {
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n thread %d cpus", myid) ;
      DVfprintf(msgFile, 7, data->cpus) ;
   }
   for ( ithread = 0 ; ithread < 7 ; ithread++ ) {
      cpus[ithread] += data->cpus[ithread] ;
   }
   *pfacops += data->facops ;
}
/*
   -------------
   free the data
   -------------
*/
ChvList_free(updlist) ;
IVL_free(rowsIVL) ;
IVfree(firstnz) ;
FREE(dataObjects) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- worker method to factor the matrix
 
 
   created -- 98may29, cca
   ----------------------------------------------------
*/
static void *
FrontMtx_QR_workerFactor (
   void   *arg
) {
char            *status ;
ChvList         *updlist ;
ChvManager      *chvmanager ;
double          facops, t0, t1 ;
double          *cpus ;
DV              workDV ;
FILE            *msgFile ;
FrontMtx        *frontmtx ;
Ideq            *dequeue ;
InpMtx          *mtxA ;
int             J, K, myid, neqns, nfront, msglvl ;
int             *colmap, *firstnz, *nactiveChild, *owners, *par ;
IVL             *rowsIVL ;
QR_factorData   *data ;

MARKTIME(t0) ;
data = (QR_factorData *) arg ;
mtxA       = data->mtxA     ;
rowsIVL    = data->rowsIVL  ;
firstnz    = data->firstnz  ;
IV_sizeAndEntries(data->ownersIV, &nfront, &owners) ;
frontmtx   = data->frontmtx   ;
chvmanager = data->chvmanager ;
updlist    = data->updlist    ;
myid       = data->myid       ;
cpus       = data->cpus       ;
msglvl     = data->msglvl     ;
msgFile    = data->msgFile    ;
par        = frontmtx->tree->par ;
neqns      = FrontMtx_neqns(frontmtx) ;
/*
   --------------------------------------------------------
   status[J] = 'F' --> J finished
             = 'W' --> J waiting to be finished
   create the Ideq object to handle the bottom-up traversal
   nactiveChild[K] = # of unfinished children of K,
      when zero, K can be placed on the dequeue
   --------------------------------------------------------
*/
status = CVinit(nfront, 'F') ;
dequeue = FrontMtx_setUpDequeue(frontmtx, owners, myid, status,
                                NULL, 'W', 'F', msglvl, msgFile) ;
FrontMtx_loadActiveLeaves(frontmtx, status, 'W', dequeue) ;
nactiveChild = FrontMtx_nactiveChild(frontmtx, status, myid) ;
colmap = IVinit(neqns, -1) ;
DV_setDefaultFields(&workDV) ;
facops = 0.0 ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n owners") ;
   IVfprintf(msgFile, nfront, owners) ;
   fprintf(msgFile, "\n Ideq") ;
   Ideq_writeForHumanEye(dequeue, msgFile) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
cpus[0] += t1 - t0 ;
/*
   ---------------------------
   loop while a path is active
   ---------------------------
*/
while ( (J = Ideq_removeFromHead(dequeue)) != -1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n ### checking out front %d, owner %d",
              J, owners[J]) ;
   }
   if ( owners[J] == myid ) {
/*
      --------------------------------
      front J is ready to be processed
      --------------------------------
*/
      FrontMtx_QR_factorVisit(frontmtx, J, mtxA, rowsIVL, firstnz,
                              updlist, chvmanager, status, colmap,
                              &workDV, cpus, &facops, msglvl, msgFile) ;
      if ( status[J] == 'F' ) {
/*
         ------------------------------------------
         front J is finished, put parent on dequeue 
         if it exists or all children are finished
         ------------------------------------------
*/
         if ( (K = par[J]) != -1 && --nactiveChild[K] == 0 ) {
            Ideq_insertAtHead(dequeue, K) ;
         }
      } else {
/*
         -----------------------------------------------
         front J is not complete, put on tail of dequeue
         -----------------------------------------------
*/
         Ideq_insertAtTail(dequeue, J) ;
      }
   } else {
/*
      -------------------------------------------
      front J is not owned, put parent on dequeue 
      if it exists and all children are finished
      -------------------------------------------
*/
      if ( (K = par[J]) != -1 && --nactiveChild[K] == 0 ) {
         Ideq_insertAtHead(dequeue, K) ;
      }
   }
}
data->facops = facops ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
CVfree(status) ;
Ideq_free(dequeue) ;
IVfree(nactiveChild) ;
IVfree(colmap) ;
DV_clearData(&workDV) ;
MARKTIME(t1) ;
cpus[6] = t1 - t0 ;
cpus[5] = t1 - t0 - cpus[0] - cpus[1]
         - cpus[2] - cpus[3] - cpus[4] ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
