/*  solveMT.c  */

#include "../spoolesMT.h"
#include "../../Ideq.h"
#include "../../timings.h"
#include "../../SolveMap.h"
#include "../../SubMtxList.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   declarations for static functions
   ---------------------------------
*/
static void * FrontMtx_workerSolve ( void *arg ) ;
static void setupListObjects ( FrontMtx *frontmtx, SolveMap *solvemap,
   SubMtxList *forwAggList, SubMtxList *backAggList,
   int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   data structure for each thread

   global data

      frontmtx    -- FrontMtx object that contains the matrix 
                     factorization
      solmtx      -- DenseMtx object, on output it holds the solution
      rhsmtx      -- DenseMtx object, on input it holds the right 
                     hand side
      solvemap    -- SolveMap object that holds the map of computations
                     to processes for the forward and backward solve
      mtxmanager  -- SubMtxManager object to manage working storage
      forwAggList -- SubMtxList to hold aggregate objects 
                     during the forward solve
      backAggList -- SubMtxList to hold aggregate objects 
                     during the backward solve
      p_mtx       -- array of pointers to SubMtx objects that point
                     to global solution SubMtx objects for the fronts
      frontIsDone -- char vector used for synchronization
        there are three steps to the solve
           (1) (L + I) Y = B
           (2) D Z = Y
           (3) (I + U) X = B
        frontIsDone[J] = 'N' --> J at step is not yet complete
        frontIsDone[J] = 'Y' --> J at step is complete
      lowerIsDone -- vector used for synchronization
         'Y' --> forward solve for iproc complete
         'N' --> forward solve for iproc not complete
      diagIsDone -- vector used for synchronization
         'Y' --> diagonal solve for iproc complete
         'N' --> diagonal solve for iproc not complete
      upperIsDone -- vector used for synchronization
         'Y' --> backward solve for iproc complete
         'N' --> backward solve for iproc not complete

   local data

      myid    -- thread id
      cpus    -- vector to hold cpu breakdown
      msglvl  -- local message level
      msgFile -- local message file
   --------------------------------------------------------------------
*/
typedef struct _SolveData   SolveData ;
struct _SolveData {
   FrontMtx        *frontmtx    ;
   DenseMtx        *solmtx      ;
   DenseMtx        *rhsmtx      ;
   SolveMap        *solvemap    ;
   SubMtxManager   *mtxmanager  ;
   int             myid         ;
   SubMtxList      *forwAggList ;
   SubMtxList      *backAggList ;
   SubMtx          **p_mtx      ;
   char            *frontIsDone ;
   char            *lowerIsDone ;
   char            *diagIsDone  ;
   char            *upperIsDone ;
   double          *cpus        ;
   int             msglvl       ;
   FILE            *msgFile     ;
} ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   multithreaded solve method for (L + I)D(I + U) X = B
   or (U^T + I)D(I + U) X = B.

   created -- 98mar19, cca
   ----------------------------------------------------
*/
void
FrontMtx_MT_solve (
   FrontMtx        *frontmtx,
   DenseMtx        *solmtx,
   DenseMtx        *rhsmtx,
   SubMtxManager   *mtxmanager,
   SolveMap        *solvemap,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) {
char         buffer[20] ;
char         *diagIsDone, *frontIsDone, *lowerIsDone, *upperIsDone ;
SubMtx       **p_mtx ;
SubMtxList   *backAggList, *forwAggList ;
SolveData    *data, *dataObjects ;
double       t1, t2 ;
FILE         *fp ;
int          ii, J, myid, nfront, nthread, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || solmtx == NULL || rhsmtx == NULL 
   || mtxmanager == NULL || solvemap == NULL || cpus == NULL 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MT_solve()"
        "\n frontmtx = %p, solmtx = %p, rhsmtx = %p, mtxmanager = %p"
        "\n solvemap = %p, cpus = %p"
        "\n  msglvl = %d, msgFile = %p"
        "\n\n bad input\n", frontmtx, solmtx, rhsmtx, mtxmanager,
        solvemap, cpus, msglvl, msgFile) ; 
   exit(-1) ;
}
MARKTIME(t1) ;
nfront  = FrontMtx_nfront(frontmtx) ;
nthread = SolveMap_nproc(solvemap) ;
/*
   --------------------------------------------------
   set up the forward and backward solve list objects
   --------------------------------------------------
*/
backAggList = SubMtxList_new() ;
forwAggList = SubMtxList_new() ;
setupListObjects(frontmtx, solvemap, forwAggList, backAggList,
                 msglvl, msgFile) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n forward aggregate list") ;
   SubMtxList_writeForHumanEye(forwAggList, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   allocate the global vectors
   ---------------------------
*/
ALLOCATE(p_mtx,  struct _SubMtx *, nfront) ;
for ( J = 0 ; J < nfront ; J++ ) {
   p_mtx[J] = NULL ;
}
frontIsDone = CVinit(nfront, 'N') ;
lowerIsDone = CVinit(solvemap->nproc, 'N') ;
diagIsDone  = CVinit(solvemap->nproc, 'N') ;
upperIsDone = CVinit(solvemap->nproc, 'N') ;
/*
   ----------------------------------
   allocate the local data structures
   ----------------------------------
*/
ALLOCATE(dataObjects, struct _SolveData, nthread) ;
for ( myid = 0, data = dataObjects ; myid < nthread ; myid++, data++ ) {
   data->frontmtx     = frontmtx    ;
   data->solmtx       = solmtx      ;
   data->rhsmtx       = rhsmtx      ;
   data->solvemap     = solvemap    ;
   data->mtxmanager   = mtxmanager  ;
   data->myid         = myid        ;
   data->forwAggList  = forwAggList ;
   data->backAggList  = backAggList ;
   data->p_mtx        = p_mtx       ;
   data->frontIsDone  = frontIsDone ;
   data->lowerIsDone  = lowerIsDone ;
   data->diagIsDone   = diagIsDone  ;
   data->upperIsDone  = upperIsDone ;
   data->cpus         = cpus        ;
   data->msglvl       = msglvl      ;
   if ( msglvl > 0 ) {
      sprintf(buffer, "solve.res.%d", myid) ;
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
   fprintf(msgFile, "\n CPU %8.3f : initialization time", t2 - t1) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------
   switch over the type of thread library
   --------------------------------------
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
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n CPU %8.3f : set concurrency time", t2 - t1) ;
}
MARKTIME(t1) ;
for ( myid = 0, data = dataObjects ;
      myid < nthread - 1 ;
      myid++, data++ ) {
   rc = thr_create(NULL, 0, FrontMtx_workerSolve, data, 0, NULL) ;
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
   fflush(msgFile) ;
}
MARKTIME(t1) ;
FrontMtx_workerSolve(data) ;
for ( myid = 0 ; myid < nthread - 1 ; myid++ ) {
   thr_join(0, 0, 0) ;
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : thread join time", t2 - t1) ;
   fflush(msgFile) ;
}
#endif
#if THREAD_TYPE == TT_POSIX
/*
   -------------
   POSIX threads
   -------------
*/
{
pthread_attr_t   attr ;
pthread_t        *tids ;
void             *status ;
/*
#####   NOTE: for SGI machines, this command must be present
#####         for the thread scheduling to be efficient.
#####         this is NOT a POSIX call, but SGI needs it anyway
pthread_setconcurrency(nthread) ;
*/
MARKTIME(t1) ;
pthread_attr_init(&attr) ;
/*
pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM) ;
*/
pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS) ;
ALLOCATE(tids, pthread_t, nthread) ;
for ( myid = 0, data = dataObjects ; myid < nthread ; myid++, data++ ) {
   rc = pthread_create(&tids[myid], &attr,
                       FrontMtx_workerSolve, data) ;
   if ( rc != 0 ) {
      fprintf(stderr,
              "\n fatal error, myid = %d, rc = %d from pthread_create",
              myid, rc) ;
      exit(-1) ;
   }
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : thread creation time", t2 - t1) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
for ( myid = 0 ; myid < nthread ; myid++ ) {
   pthread_join(tids[myid], &status) ;
}
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : thread join time", t2 - t1) ;
   fflush(msgFile) ;
}
FREE(tids) ;
pthread_attr_destroy(&attr) ;
}
#endif
/*
   ---------------------
   free the working data
   ---------------------
*/
MARKTIME(t1) ;
DVzero(8, cpus) ;
for ( myid = 0, data = dataObjects ; myid < nthread ; myid++, data++ ) {
   for ( ii = 0 ; ii < 8 ; ii++ ) {
      cpus[ii] += data->cpus[ii] ;
   }
}
FREE(dataObjects) ;
FREE(p_mtx) ;
CVfree(frontIsDone) ;
CVfree(lowerIsDone) ;
CVfree(diagIsDone) ;
CVfree(upperIsDone) ;
SubMtxList_free(forwAggList) ;
SubMtxList_free(backAggList) ;
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %8.3f : free working data", t2 - t1) ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   worker method to solve the linear system
 
   created -- 97jun27, cca
   ----------------------------------------
*/
static void *
FrontMtx_workerSolve (
   void   *arg
) {
char            *diagIsDone, *frontIsDone, *lowerIsDone, 
                *status, *upperIsDone ;
DenseMtx        *rhsmtx, *solmtx ;
SubMtx          **p_agg, **p_mtx ;
SubMtxList      *backAggList, *forwAggList ;
SubMtxManager   *mtxmanager ;
FrontMtx        *frontmtx ;
double          t0, t1, t2, t3 ;
double          *cpus ;
SolveData       *data ;
FILE            *msgFile ;
Ideq            *dequeue ;
int             everybodyDone, iproc, I, J, K, msglvl, myid, 
                nfront, nproc, nrhs ;
int             *fch, *nactiveChild, *owners, *par, *sib ;
IP              **heads ;
SolveMap        *solvemap ;
Tree            *tree ;
/*
   -------------------------------
   extract pointers and dimensions
   -------------------------------
*/
MARKTIME(t0) ;
data        = (SolveData *) arg ;
frontmtx    = data->frontmtx ;
rhsmtx      = data->rhsmtx ;
solmtx      = data->solmtx ;
solvemap    = data->solvemap ;
mtxmanager  = data->mtxmanager  ;
forwAggList = data->forwAggList ;
backAggList = data->backAggList ;
p_mtx       = data->p_mtx     ;
frontIsDone = data->frontIsDone ;
lowerIsDone = data->lowerIsDone ;
diagIsDone  = data->diagIsDone  ;
upperIsDone = data->upperIsDone ;
myid        = data->myid ;
cpus        = data->cpus ;
msglvl      = data->msglvl ;
msgFile     = data->msgFile ;
nfront      = frontmtx->nfront ;
nproc       = solvemap->nproc ;
tree        = frontmtx->frontETree->tree ;
par         = tree->par ;
fch         = tree->fch ;
sib         = tree->sib ;
nrhs        = rhsmtx->ncol ;
owners      = SolveMap_owners(solvemap) ;
/*
   -----------------------------------------------
   set up for the forward solve
   (1) IP lists for the forward solve
   (2) local status[] vector 
   (3) create the dequeue for the forward solve
   (4) create the number of active children vector
   -----------------------------------------------
*/
MARKTIME(t1) ;
heads   = SolveMap_forwardSetup(solvemap, myid, msglvl, msgFile) ;
status  = CVinit(nfront, 'F') ;
dequeue = FrontMtx_setUpDequeue(frontmtx, owners, myid, status, heads,
                                'W', 'F', msglvl, msgFile) ;
FrontMtx_loadActiveLeaves(frontmtx, status, 'W', dequeue) ;
nactiveChild = FrontMtx_nactiveChild(frontmtx, status, myid) ;
MARKTIME(t2) ;
cpus[0] += t2 - t1 ;
if ( msglvl > 2 ) {
   IP   *ip ;
   fprintf(msgFile, "\n\n forward update lists") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( heads[J] != NULL ) {
         fprintf(msgFile, "\n %d : ", J) ;
         for ( ip = heads[J] ; ip != NULL ; ip = ip->next ) {
            fprintf(msgFile, " %d", ip->val) ;
         }
      }
   }
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n dequeue set up for forward solve") ;
   Ideq_writeForHumanEye(dequeue, msgFile) ;
   fprintf(msgFile, "\n\n nactiveChild") ;
   IVfprintf(msgFile, nfront, nactiveChild) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   load the right hand side
   ------------------------
*/
MARKTIME(t1) ;
p_agg = FrontMtx_loadRightHandSide(frontmtx, rhsmtx, owners, myid, 
                                   mtxmanager, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[1] += t2 - t1 ;
/*
   -------------
   forward solve
   -------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n\n thread %d : starting forward solve", myid) ;
fflush(stdout) ;
#endif
MARKTIME(t1) ;
while ( (J = Ideq_removeFromHead(dequeue)) != -1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile,
              "\n\n ### forward solve, checking out J = %d, status %c",
              J, status[J]) ;
      fflush(msgFile) ;
   }
   FrontMtx_forwardVisit(frontmtx, J, nrhs, owners, myid, mtxmanager, 
                         forwAggList, p_mtx, frontIsDone, heads, 
                         p_agg, status, msglvl, msgFile) ;
   if ( status[J] == 'F' ) {
/*
      -------------------
      front J is finished
      -------------------
*/
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n front J is finished") ;
         fflush(msgFile) ;
      }
      if ( (K = par[J]) != -1 && --nactiveChild[K] == 0 ) {
/*
         -----------------------------------
         all children of parent are finished
         add parent to head of dequeue
         -----------------------------------
*/
         if ( msglvl > 1 ) {
            fprintf(msgFile, ", placing %d on head of dequeue", K) ;
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
}
IP_free(heads[nfront+1]) ;
FREE(heads) ;
IVfree(nactiveChild) ;
Ideq_free(dequeue) ;
MARKTIME(t2) ;
cpus[2] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n forward solve finished") ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------
   all processes wait here until all are done
   ------------------------------------------
*/
lowerIsDone[myid] = 'Y' ;
do {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n 1. checking out lowerIsDone[] vector: ") ;
      fflush(msgFile) ;
   }
   everybodyDone = 1 ;
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, " %c", lowerIsDone[iproc]) ;
         fflush(msgFile) ;
      }
      if ( lowerIsDone[iproc] != 'Y' ) {
         everybodyDone = 0 ;
         break ;
      }
   }
} while ( everybodyDone == 0 ) ;
/*
   ---------------------
   do the diagonal solve
   ---------------------
*/
MARKTIME(t1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( owners[J] == myid ) {
      FrontMtx_diagonalVisit(frontmtx, J, owners, myid, p_mtx,
                              frontIsDone, p_agg, msglvl, msgFile) ;
/*
      ------------------------------------------------------------
      set frontIsDone[J] to 'N' to be ready for the backward solve
      ------------------------------------------------------------
*/
      frontIsDone[J] = 'N' ;
   }
}
MARKTIME(t2) ;
cpus[3] += t2 - t1 ;
/*
   ------------------------------------------
   all processes wait here until all are done
   ------------------------------------------
*/
diagIsDone[myid] = 'Y' ;
do {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n 2. checking out diagIsDone[] vector: ") ;
      fflush(msgFile) ;
   }
   everybodyDone = 1 ;
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, " %c", diagIsDone[iproc]) ;
         fflush(msgFile) ;
      }
      if ( diagIsDone[iproc] != 'Y' ) {
         everybodyDone = 0 ;
         break ;
      }
   }
} while ( everybodyDone == 0 ) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n everybody is done") ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   set up for the backward solve
   (1) IP lists for the backward solve
   (2) set up the dequeue for the backward solve
   ---------------------------------------------
*/
MARKTIME(t1) ;
heads = SolveMap_backwardSetup(solvemap, myid, msglvl, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n ") ;
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
dequeue = FrontMtx_setUpDequeue(frontmtx, owners, myid, status, heads,
                                 'W', 'F', msglvl, msgFile) ;
FrontMtx_loadActiveRoots(frontmtx, status, 'W', dequeue) ;
MARKTIME(t2) ;
cpus[1] += t2 - t1 ;
/*
   --------------
   backward solve
   --------------
*/
#if MYDEBUG > 0
fflush(stdout) ;
#endif
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n thread %d : starting backward solve", myid) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
while ( (J = Ideq_removeFromHead(dequeue)) != -1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile,
              "\n\n ### backward solve, checking out J = %d, status %c",
              J, status[J]) ;
      fflush(msgFile) ;
   }
   FrontMtx_backwardVisit(frontmtx, J, nrhs, owners, myid, mtxmanager, 
                           backAggList, p_mtx, frontIsDone, heads, 
                           p_agg, status, msglvl, msgFile) ;
   if ( status[J] == 'F' ) {
/*
      -------------------
      front J is finished
      -------------------
*/
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n front J is finished") ;
         fflush(msgFile) ;
      }
      for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
         if ( status[I] == 'W' ) {
            Ideq_insertAtHead(dequeue, I) ;
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
}
CVfree(status) ;
IP_free(heads[nfront+1]) ;
FREE(heads) ;
Ideq_free(dequeue) ;
MARKTIME(t2) ;
cpus[4] += t2 - t1 ;
/*
   ------------------------------------------
   all processes wait here until all are done
   ------------------------------------------
*/
upperIsDone[myid] = 'Y' ;
do {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n 2. checking out upperIsDone[] vector: ") ;
      fflush(msgFile) ;
   }
   everybodyDone = 1 ;
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, " %c", upperIsDone[iproc]) ;
         fflush(msgFile) ;
      }
      if ( upperIsDone[iproc] != 'Y' ) {
         everybodyDone = 0 ;
         break ;
      }
   }
} while ( everybodyDone == 0 ) ;
/*
   --------------------------
   store the solution entries
   --------------------------
*/
MARKTIME(t1) ;
FrontMtx_storeSolution(frontmtx, owners, myid, mtxmanager, 
                        p_mtx, solmtx, msglvl, msgFile) ;
MARKTIME(t2) ;
cpus[5] += t2 - t1 ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
FREE(p_agg) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   set up the aggregate list objects for the forward and back solves

   created -- 98mar20, cca
   -----------------------------------------------------------------
*/
static void
setupListObjects (
   FrontMtx      *frontmtx,
   SolveMap      *solvemap,
   SubMtxList    *forwAggList,
   SubMtxList    *backAggList,
   int           msglvl, 
   FILE          *msgFile
) {
char   *flags ;
int    J, nfront ;
int    *aggcounts ;
IV     *aggcountsIV ;

nfront = SolveMap_nfront(solvemap) ;
flags  = CVinit(nfront, 'N') ;
/*
   ------------------------------------------------------
   set up the aggregate list object for the forward solve
   ------------------------------------------------------
*/
aggcountsIV = SolveMap_lowerAggregateIV(solvemap, -1, msglvl, msgFile) ;
aggcounts = IV_entries(aggcountsIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( aggcounts[J] > 1 ) {
      flags[J] = 'Y' ;
   }
}
SubMtxList_init(forwAggList, nfront, aggcounts, 1, flags) ;
IV_free(aggcountsIV) ;
/*
   -------------------------------------------------------
   set up the aggregate list object for the backward solve
   -------------------------------------------------------
*/
aggcountsIV = SolveMap_upperAggregateIV(solvemap, -1, msglvl, msgFile) ;
aggcounts = IV_entries(aggcountsIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( aggcounts[J] > 1 ) {
      flags[J] = 'Y' ;
   }
}
SubMtxList_init(backAggList, nfront, aggcounts, 1, flags) ;
IV_free(aggcountsIV) ;
CVfree(flags) ;

return ; }

/*--------------------------------------------------------------------*/
