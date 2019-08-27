/*  factorUtil.c  */

#include "../FrontMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   purpose -- initialize a front
 
   created -- 98may04, cca
   -----------------------------
*/
void
FrontMtx_initializeFront (
   FrontMtx   *frontmtx,
   Chv        *frontJ,
   int        J
) {
int   ncolJ, nD, nrowJ ;
int   *colindJ, *ivec ;
/*
   -----------------------------------
   get the number of internal vertices
   -----------------------------------
*/
nD = ETree_frontSize(frontmtx->frontETree, J) ;
/*
   --------------------------------------
   get the internal and boundary vertices
   --------------------------------------
*/
IVL_listAndSize(frontmtx->symbfacIVL, J, &ncolJ, &colindJ) ;
#if MYDEBUG > 0
fprintf(stdout, "\n front %d, column indices", J) ;
IVfprintf(stdout, ncolJ, colindJ) ;
fflush(stdout) ;
#endif
/*
   ------------------------------------------------------
   initialize the front's dimensions, indices and entries
   ------------------------------------------------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n front %d, nD %d, nU %d", J, nD, ncolJ - nD) ;
fflush(stdout) ;
#endif
Chv_init(frontJ, J, nD, ncolJ - nD, ncolJ - nD,
         frontmtx->type, frontmtx->symmetryflag) ;
/*
   -----------------------
   fill the column indices
   -----------------------
*/
Chv_columnIndices(frontJ, &ncolJ, &ivec) ;
IVcopy(ncolJ, ivec, colindJ) ;
#if MYDEBUG > 0
fprintf(stdout, "\n front %d, colind = %p", J, ivec) ;
IVfprintf(stdout, ncolJ, ivec) ;
fflush(stdout) ;
#endif
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
   --------------------
   fill the row indices
   --------------------
*/
   Chv_rowIndices(frontJ, &nrowJ, &ivec) ;
   IVcopy(nrowJ, ivec, colindJ) ;
}
/*
   ------------------------
   zero the front's entries
   ------------------------
*/
Chv_zero(frontJ) ;
#if MYDEBUG > 0
fprintf(stdout, "\n leaving Front_initializeFront()") ;
Chv_writeForHumanEye(frontJ, stdout) ;
fflush(stdout) ;
#endif
 
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------
   static functions
   ----------------
*/
static void assembleAggregates ( Chv *frontJ, ChvList *aggList,
   ChvManager *chvmanager, double cpus[], int msglvl, FILE *msgFile ) ;
static char assemblePostponedData ( FrontMtx *frontmtx, Chv *frontJ,
   int *pndelay, Chv *fronts[], ChvList *postList,
   ChvManager *chvmanager, double cpus[], int msglvl, FILE *msgFile ) ;
static int factorFront ( FrontMtx *frontmtx, Chv *frontJ,
   int ndelay, double tau, IV *pivotsizesIV, DV *workDV, int stats[],
   double cpus[], int msglvl, FILE *msgFile ) ;
static void storeEntries ( FrontMtx *frontmtx, Chv *frontJ, int nelim,
   double droptol, IV *pivotsizesIV, ChvList *postList,
   ChvManager *chvmanager, int parent[], int stats[], double cpus[],
   int msglvl, FILE *msgFile ) ;
/*
   ------------------------------------------------------------------
   purpose -- to visit a front during a factorization.
              note: this method is called by the serial,
              multithreaded and MPI factorization codes.

   frontmtx     -- front matrix object
   pencil       -- matrix pencil for A + sigma*B
   J            -- id of front we are working on
   myid         -- id of thread or process 
   owners[]     -- map from fronts to owning threads,
                   in a serial environment, owners = NULL
   fronts[]     -- vector of pointers to fronts
   lookahead    -- parameter controls the partial upward visiting
                   of ancestors before returning
   tau          -- used when pivoting enabled,
                   |L_{j,i}| and |U_{i,j}| <= tau
   droptol      -- used for an approximate factorization
                   stored |L_{j,i}| and |U_{i,j}| > droptol
   status[]     -- status vector for the fronts,
                   'W' -- needs to be woke up
                   'A' -- active front
                   'F' -- front is finished
   heads[]      -- vector of pointers to IP objects that store the
                   front-to-front update lists
   pivotsizesIV -- IV object used during the factorization of a front
                   when pivoting is enabled
   workDV       -- DV object used for working storage
   parent       -- front parent vector
   aggList      -- list object used in parallel environment, used
                   to store aggregate fronts
   postList     -- list object used in pivoting and/or parallel
                   environments, used to stored delayed data and/or
                   to synchronize the factorization
   chvmanager   -- used to manage working storage for Chv objects
   stats[]      -- statistics vector
   cpus[]       -- vector to hold breakdown of cpu times
   msglvl       -- message level
   msgFil       -- message file

   created -- 98may04, cca
   ------------------------------------------------------------------
*/
char
FrontMtx_factorVisit (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   int          J,
   int          myid,
   int          owners[],
   Chv          *fronts[],
   int          lookahead,
   double       tau,
   double       droptol,
   char         status[],
   IP           *heads[],
   IV           *pivotsizesIV,
   DV           *workDV,
   int          parent[],
   ChvList      *aggList,
   ChvList      *postList,
   ChvManager   *chvmanager,
   int          stats[],
   double       cpus[],
   int          msglvl,
   FILE         *msgFile
) {
/*
   ---------------
   local variables
   ---------------
*/
char     allAggregatesHere, allPostponedAssmb, allUpdatesDone ;
Chv     *frontJ ;
double   t1, t2 ;
int      K, ndelay, nelim ;

if ( status[J] == 'F' ) {
   return('F') ;
}
allUpdatesDone    = 'N' ;
allAggregatesHere = 'N' ;
allPostponedAssmb = 'N' ;
frontJ = NULL ;
if ( heads[J] != NULL ) {
/*
   --------------------------------
   internal updates to be performed
   --------------------------------
*/
   if ( (frontJ = fronts[J]) == NULL ) {
/*
      --------------------
      initialize the front
      --------------------
*/
      fronts[J] = FrontMtx_setupFront(frontmtx, pencil, J, myid, 
                      owners, chvmanager, cpus, msglvl, msgFile) ;
      frontJ = fronts[J] ;
      status[J] = 'A' ;
   }
/*
   ---------------------------------
   compute updates from owned fronts
   ---------------------------------
*/
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n performing updates") ;
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   FrontMtx_update(frontmtx, frontJ, heads,
                    status, workDV, msglvl, msgFile) ;
   MARKTIME(t2) ;
   cpus[2] += t2 - t1 ;
}
if ( heads[J] == NULL ) {
   allUpdatesDone = 'Y' ;
}
if ( owners == NULL || owners[J] == myid ) {
/*
   --------------
   front is owned
   --------------
*/
   if ( (frontJ = fronts[J]) == NULL ) {
/*
      --------------------
      initialize the front
      --------------------
*/
      fronts[J] = FrontMtx_setupFront(frontmtx, pencil, J, myid, 
                            owners, chvmanager, cpus, msglvl, msgFile) ;
      frontJ = fronts[J] ;
      status[J] = 'A' ;
   }
   if ( aggList != NULL ) {
/*
      ------------------------------------------
      we are operating in a parallel environment
      ------------------------------------------
*/
      if ( ChvList_isListNonempty(aggList, J) == 1 ) {
/*
         -------------------------------------------------
         assemble any aggregate updates from other threads
         -------------------------------------------------
*/
         assembleAggregates(frontJ, aggList, chvmanager, 
                            cpus, msglvl, msgFile) ;
      }
      if ( ChvList_isCountZero(aggList, J) == 1 ) {
/*
         -------------------------------------------------------
         all aggregates are assembled or in the list to assemble
         -------------------------------------------------------
*/
         if ( ChvList_isListNonempty(aggList, J) == 1 ) {
/*
            -------------------------------------------------
            assemble any aggregate updates from other threads
            -------------------------------------------------
*/
            assembleAggregates(frontJ, aggList, chvmanager, 
                               cpus, msglvl, msgFile) ;
         }
         allAggregatesHere = 'Y' ;
      }
   } else {
      allAggregatesHere = 'Y' ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n allUpdatesDone %c, allAggregatesHere %c",
              allUpdatesDone, allAggregatesHere) ;
      fflush(msgFile) ;
   }
   if ( allUpdatesDone == 'Y' && allAggregatesHere == 'Y' ) {
/*
      -------------------------------------
      all internal updates and all external 
      aggregates have been incorporated
      -------------------------------------
*/
      if ( postList != NULL ) {
/*
         ---------------------------------------------
         front is ready to assemble any postponed data
         ---------------------------------------------
*/
         allPostponedAssmb = assemblePostponedData(frontmtx, frontJ, 
                                    &ndelay, fronts, postList, 
                                    chvmanager, cpus, msglvl, msgFile) ;
         frontJ = fronts[J] ;
      } else {
         allPostponedAssmb = 'Y' ;
         ndelay = 0 ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
                 "\n\n allPostponedAssmb %c", allPostponedAssmb) ;
         fflush(msgFile) ;
      }
      if ( allPostponedAssmb == 'Y' ) {
/*
         -----------------------------
         front is ready to be factored
         -----------------------------
*/
         nelim = factorFront(frontmtx, frontJ, ndelay, tau, 
                             pivotsizesIV, workDV, stats, 
                             cpus, msglvl, msgFile) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                    "\n\n J = %d, nelim = %d", frontJ->id, nelim) ;
            fflush(msgFile) ;
         }
         if ( (! FRONTMTX_IS_PIVOTING(frontmtx)) && nelim < frontJ->nD){
/*
            ------------------------------------------------
            a singular front matrix has been found.
            release the front and set the status to error
            ------------------------------------------------
*/
            ChvManager_releaseObject(chvmanager, frontJ) ;
            fronts[J] = NULL ;
            status[J] = 'E' ;
         } else {
/*
            -------------------------------------------
            store any factor entries and postponed data
            -------------------------------------------
*/
            storeEntries(frontmtx, frontJ, nelim, droptol, pivotsizesIV,
                         postList, chvmanager, parent, stats, 
                         cpus, msglvl, msgFile) ;
/*
            ------------------------------------------------
            release the front and set the status to finished
            ------------------------------------------------
*/
            ChvManager_releaseObject(chvmanager, frontJ) ;
            fronts[J] = NULL ;
            status[J] = 'F' ;
         }
      }
   }
} else if ( allUpdatesDone == 'Y' ) {
   if ( frontJ != NULL ) {
/*
      --------------------------------------------------------
      front is not owned and all internal updates are complete
      put aggregate update front on the aggregate list
      --------------------------------------------------------
*/
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n done with unowned front %3d", J) ;
         fflush(msgFile) ;
      }
      if ( msglvl > 3 ) {
         Chv_writeForHumanEye(frontJ, msgFile) ;
         fflush(msgFile) ;
      }
      MARKTIME(t1) ;
      ChvList_addObjectToList(aggList, frontJ, J) ;
      MARKTIME(t2) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n thread %2d, done with unowned front %3d",
              myid, J) ;
      fflush(stdout) ;
#endif
      cpus[7] += t2 - t1 ;
   }
   status[J] = 'F' ;
}
if ( --lookahead >= 0 && (K = parent[J]) != -1 ) {
/*
   -----------------------------
   visit parent before returning
   -----------------------------
*/
   FrontMtx_factorVisit(frontmtx, pencil, J, myid, owners, fronts,
                   lookahead, tau, droptol, status, heads, pivotsizesIV,
                   workDV, parent, aggList, postList, chvmanager,
                   stats, cpus, msglvl, msgFile) ;
}
return(status[J]) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- set up the front's data structure

   created -- 98mar27, cca
   --------------------------------------------
*/
Chv *
FrontMtx_setupFront (
   FrontMtx     *frontmtx,
   Pencil       *pencil,
   int          J,
   int          myid,
   int          owners[],
   ChvManager   *chvmanager,
   double       cpus[],
   int          msglvl,
   FILE         *msgFile
) {
Chv     *frontJ ;
double   t1, t2 ;
int      nbytes, nD, nL, nU ;

if ( msglvl > 4 ) {
   fprintf(msgFile, 
           "\n\n inside FrontMtx_setupFront()"
           "\n frontmtx %p, pencil %p, J %d, myid %d"
           "\n owners %p, chvmanager %p, cpus %p"
           "\n msglvl %d, msgFile %p",
           frontmtx, pencil, J, myid, owners, chvmanager,
           cpus, msglvl, msgFile) ; 
   fflush(msgFile) ;
}
MARKTIME(t1) ;
FrontMtx_initialFrontDimensions(frontmtx, J,
                                 &nD, &nL, &nU, &nbytes) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nD %d, nL %d, nU %d, nbytes %d",
           nD, nL, nU, nbytes) ;
   fflush(msgFile) ;
}
frontJ = ChvManager_newObjectOfSizeNbytes(chvmanager, nbytes) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n frontJ = %p", frontJ) ;
   fflush(msgFile) ;
}
Chv_init(frontJ, J, nD, nL, nU, frontmtx->type, frontmtx->symmetryflag);
FrontMtx_initializeFront(frontmtx, frontJ, J) ;
MARKTIME(t2) ;
cpus[0] += t2 - t1 ;
if ( pencil != NULL && (owners == NULL || owners[J] == myid) ) {
/*
   -------------------------------------------------
   thread owns this front, load the original entries
   -------------------------------------------------
*/
   MARKTIME(t1) ;
   FrontMtx_loadEntries(frontJ, pencil, msglvl, msgFile) ;
   MARKTIME(t2) ;
   cpus[1] += t2 - t1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n original entries loaded") ;
      fflush(msgFile) ;
   }
/*
   if ( FRONTMTX_IS_HERMITIAN(frontmtx) && J == frontmtx->nfront - 1 ) {
      fprintf(msgFile, "\n last front after entries loaded") ;
      Chv_writeForHumanEye(frontJ, msgFile) ;
      fflush(msgFile) ;
   }
*/
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n front initialized") ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   Chv_writeForHumanEye(frontJ, msgFile) ;
   fflush(msgFile) ;
}
return(frontJ) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- assemble any aggregate updates from other threads
   ------------------------------------------------------------
*/
static void
assembleAggregates (
   Chv          *frontJ,
   ChvList      *aggList,
   ChvManager   *chvmanager,
   double       cpus[],
   int          msglvl,
   FILE         *msgFile
) {
Chv     *chv, *headchv ;
double   t1, t2 ;

MARKTIME(t1) ;
headchv = ChvList_getList(aggList, frontJ->id) ;
for ( chv = headchv ; chv != NULL ; chv = chv->next ) {
   Chv_assembleChv(frontJ, chv) ;
}
MARKTIME(t2) ;
cpus[8] += t2 - t1 ;
ChvManager_releaseListOfObjects(chvmanager, headchv) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n after assembly") ;
   Chv_writeForHumanEye(frontJ, msgFile) ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- assemble any postponed data

   created -- 98mar27, cca
   --------------------------------------
*/
static char
assemblePostponedData (
   FrontMtx     *frontmtx,
   Chv          *frontJ,
   int          *pndelay,
   Chv          *fronts[],
   ChvList      *postList,
   ChvManager   *chvmanager,
   double       cpus[],
   int          msglvl,
   FILE         *msgFile
) {
char     rc ;
double   t1, t2 ;
int      J ;

if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n frontmtx %p, frontJ %p, pndelay %p"
           "\n fronts %p, postList %p, chvmanager %p, cpus %p",
           frontmtx, frontJ, pndelay, 
           fronts, postList, chvmanager, cpus) ;
   fflush(msgFile) ;
}


J = frontJ->id ;
if ( msglvl > 1 ) {
   if ( ChvList_isCountZero(postList, J) == 1) {
      fprintf(msgFile, "\n all postponed data is here") ;
      fflush(msgFile) ;
   } else {
      fprintf(msgFile, "\n still waiting on postponed data") ;
      fflush(msgFile) ;
   }
}
if ( ChvList_isCountZero(postList, J) == 1 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n assembling postponed data") ;
      fflush(msgFile) ;
   }
/*
   ---------------------------
   assemble any postponed data
   ---------------------------
*/
   MARKTIME(t1) ;
   fronts[J] = FrontMtx_assemblePostponedData(frontmtx, frontJ, 
                                   postList, chvmanager, pndelay) ;
   if ( frontJ != fronts[J] ) {
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n releasing old front") ;
         fflush(msgFile) ;
      }
      ChvManager_releaseObject(chvmanager, frontJ) ;
   }
   MARKTIME(t2) ;
   cpus[3] += t2 - t1 ;
   rc = 'Y' ;
} else {
   rc = 'N' ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   purpose -- factor the front

   created -- 98mar27, cca
   ---------------------------
*/
static int
factorFront (
   FrontMtx   *frontmtx,
   Chv        *frontJ,
   int        ndelay,
   double     tau,
   IV         *pivotsizesIV,
   DV         *workDV,
   int        stats[],
   double     cpus[],
   int        msglvl,
   FILE       *msgFile
) {
double   t1, t2 ;
int      nelim, npost ;

if ( msglvl > 2 ) {
   fprintf(msgFile, "\n frontJ = %p, ndelay = %d", 
           frontJ, ndelay) ;
   fprintf(msgFile, "\n tau = %12.4e", tau) ;
   fprintf(msgFile, "\n stats %p, cpus %p", stats, cpus) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   Chv_writeForHumanEye(frontJ, msgFile) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   nelim = Chv_factorWithPivoting(frontJ, ndelay, 
                                  frontmtx->pivotingflag, pivotsizesIV, 
                                  workDV, tau, &stats[1]) ;
} else {
   nelim = Chv_factorWithNoPivoting(frontJ, frontmtx->patchinfo) ;
}
npost = frontJ->nD - nelim ;
stats[2] += npost ;
if (  !(FRONTMTX_IS_PIVOTING(frontmtx))
   || FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   stats[0] += nelim ;
} else {
   stats[0] += IV_size(pivotsizesIV) ;
}
MARKTIME(t2) ;
cpus[4] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n front %d, nelim = %d, npost = %d",
           frontJ->id, nelim, npost) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   Chv_writeForHumanEye(frontJ, msgFile) ;
   fflush(msgFile) ;
}
return(nelim) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- store the factor entries in the front matrix object,
              and extract any postponed data and put it in the
              list object

   created -- 98mar27, cca
   ---------------------------------------------------------------
*/
static void
storeEntries (
   FrontMtx     *frontmtx,
   Chv          *frontJ,
   int          nelim,
   double       droptol,
   IV           *pivotsizesIV,
   ChvList      *postList,
   ChvManager   *chvmanager,
   int          parent[],
   int          stats[],
   double       cpus[],
   int          msglvl,
   FILE         *msgFile
) {
double   t1, t2 ;
int      npost ;

npost = frontJ->nD - nelim ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n storing factor data, nelim = %d", nelim) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
frontJ->nD -= npost ;
frontJ->nL += npost ;
frontJ->nU += npost ;
FrontMtx_storeFront(frontmtx, frontJ, pivotsizesIV, droptol,
                     msglvl, msgFile) ;
frontJ->nD += npost ;
frontJ->nL -= npost ;
frontJ->nU -= npost ;
MARKTIME(t2) ;
cpus[6] += t2 - t1 ;
if ( postList != NULL ) {
   Chv   *postponedchv ;
   if ( npost > 0 ) {
/*
      ---------------------------------
      there is postponed data,
      extract and put on postponed list
      ---------------------------------
*/
      postponedchv = frontJ ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n postponed data for front %d", 
                 frontJ->id) ;
         Chv_writeForHumanEye(postponedchv, msgFile) ;
         fflush(msgFile) ;
      }
   } else {
      postponedchv = NULL ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n storing postponed data %p", postponedchv);
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   FrontMtx_storePostponedData(frontmtx, postponedchv, npost, 
                       parent[frontJ->id], postList, chvmanager) ;
   MARKTIME(t2) ;
   cpus[5] += t2 - t1 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- to set up the link data structures 
      for a parallel factorization for process myid
 
   return value -- pointer to IP *heads[nfront+2], which contains
      the beginning of a list of IP objects that store the remaining
      updates to the fronts. 
      note, heads[nfront] is the first IP object in the free list.
      heads[nfront+1] is the base address of the allocated IP objects.
 
   created -- 98mar07, cca
   --------------------------------------------------------------------
*/
IP **
FrontMtx_factorSetup (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        myid,
   int        msglvl,
   FILE       *msgFile
) {
int   count, ii, J, K, nadj, nfront ;
int   *adj, *mark, *owners, *vtxToFront ;
IP    *ip ;
IP    **heads ;
IVL   *symbfacIVL ;
/*
   ---------------------------------------------------
   count the number of updates this front will perform
   ---------------------------------------------------
*/
nfront = FrontMtx_nfront(frontmtx) ;
if ( frontOwnersIV != NULL ) {
   owners = IV_entries(frontOwnersIV) ;
} else {
   owners = NULL ;
}
symbfacIVL = frontmtx->symbfacIVL ;
vtxToFront = ETree_vtxToFront(frontmtx->frontETree) ;
mark = IVinit(nfront, -1) ;
for ( J = count = 0 ; J < nfront ; J++ ) {
   if ( owners == NULL || owners[J] == myid ) {
      IVL_listAndSize(symbfacIVL, J, &nadj, &adj) ;
      mark[J] = J ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         K = vtxToFront[adj[ii]] ;
         if ( mark[K] != J ) {
            mark[K] = J ;
            count++ ;
         }
      }
   }
}
/*
   --------------------------------------------------
   set up the update head/links for the factorization
   --------------------------------------------------
*/
ALLOCATE(heads, struct _IP *, nfront + 2) ;
for ( J = 0 ; J <= nfront + 1 ; J++ ) {
   heads[J] = NULL ;
}
heads[nfront] = heads[nfront+1] = IP_init(count, 1) ;
IVfill(nfront, mark, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( owners == NULL || owners[J] == myid ) {
      IVL_listAndSize(symbfacIVL, J, &nadj, &adj) ;
      mark[J] = J ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         K = vtxToFront[adj[ii]] ;
         if ( mark[K] != J ) {
            mark[K] = J ;
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
}
IVfree(mark) ;

return(heads) ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   create and return the nactiveChild vector.
   nactiveChild[J] contains the number of children 
   of J that belong to an active path
                  
   created -- 97jul03, cca
   -----------------------------------------------
*/
int *
FrontMtx_nactiveChild (
   FrontMtx   *frontmtx,
   char       *status,
   int        myid
) {
int    J, K, nfront ;
int    *nactiveChild, *par ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || status == NULL || myid < 0 ) {
   fprintf(stderr, "\n fatal error in FrontMtx_nativeChild(%p,%p,%d)"
           "\n bad input\n", frontmtx, status, myid) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
par    = ETree_par(frontmtx->frontETree) ;
/*
   ---------------------------------------
   create and fill the nactiveChild vector
   ---------------------------------------
*/
nactiveChild = IVinit(nfront, 0) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( status[J] == 'W' && (K = par[J]) != -1 ) {
      nactiveChild[K]++ ;
   }
}
return(nactiveChild) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   create, initialize and return a Ideq object
   that will be used for a parallel factorization,
   forward solve, or backward solve.

   the status[] vector will be set as follows:
      status[J] = activeflag   if J is on an active path
      status[J] = inactiveflag if J is not on an active path

   created -- 98mar27, cca
   ---------------------------------------------------------
*/
Ideq *
FrontMtx_setUpDequeue (
   FrontMtx   *frontmtx,
   int        owners[],
   int        myid,
   char       status[],
   IP         *heads[],
   char       activeFlag,
   char       inactiveFlag,
   int        msglvl,
   FILE       *msgFile
) {
Ideq   *dequeue ;
int    J, K, nfront, npath ;
int    *par ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || owners == NULL 
   || status == NULL || myid < 0 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_setUpDequeue()"
           "\n frontmtx %p, owners %p, myid %d, status %p"
           "\n heads %p, activeFlag %c, inactiveFlag %c"
           "\n bad input\n", 
           frontmtx, owners, myid, status, heads,
           activeFlag, inactiveFlag) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
tree   = frontmtx->tree ;
par    = tree->par ;
/*
   ------------------------------------
   count the number of active paths in
   the tree and set the status[] vector
   ------------------------------------
*/
CVfill(nfront, status, inactiveFlag) ;
for ( J = npath = 0 ; J < nfront ; J++ ) {
   if ( status[J] == inactiveFlag ) {
      if ( owners[J] == myid || (heads != NULL && heads[J] != NULL) ) {
         npath++ ;
         for ( K = J ; 
               K != -1 && status[K] != activeFlag ; 
               K = par[K] ) {
            status[K] = activeFlag ;
         }
      }
   }
}
dequeue = Ideq_new() ;
Ideq_resize(dequeue, npath) ;

return(dequeue) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- load the dequeue with the leaves of the active subtree
              used for a factorization and forward solve

   created -- 98mar27, cca
   -----------------------------------------------------------------
*/
void
FrontMtx_loadActiveLeaves (
   FrontMtx   *frontmtx,
   char       status[],
   char       activeFlag,
   Ideq       *dequeue
) {
int   I, J, nactivechild, nfront ;
int   *fch, *sib ;

nfront = frontmtx->nfront ;
fch    = frontmtx->tree->fch ;
sib    = frontmtx->tree->sib ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( status[J] == activeFlag ) {
      if ( fch[J] == -1 ) {
         Ideq_insertAtTail(dequeue, J) ;
      } else {
         nactivechild = 0 ;
         for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
            if ( status[I] == activeFlag ) {
               nactivechild++ ;
               break ;
            }
         }
         if ( nactivechild == 0 ) {
            Ideq_insertAtTail(dequeue, J) ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   create, initialize and return a ChvList object
   to deal with postponed chevrons
                  
   created -- 97jul03, cca
   -----------------------------------------------
*/
ChvList *
FrontMtx_postList (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        lockflag
) {
char      *flags ;
ChvList   *postList ;
int       count, I, J, jthread, nchild, nfront, nthread ;
int       *counts, *fch, *frontOwners, *mark, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || frontOwnersIV == NULL 
   || lockflag < 0 || lockflag > 2 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_postList(%p,%p,%d)"
           "\n bad input\n", frontmtx, frontOwnersIV, lockflag) ;
   exit(-1) ;
}
fch = ETree_fch(frontmtx->frontETree) ;
sib = ETree_sib(frontmtx->frontETree) ;
IV_sizeAndEntries(frontOwnersIV, &nfront, &frontOwners) ;
counts = IVinit(nfront+1, 0) ;
if ( lockflag > 0 ) {
   flags = CVinit(nfront+1, 'N') ;
} else {
   flags = NULL ;
}
nthread = 1 + IV_max(frontOwnersIV) ;
mark    = IVinit(nthread, -1) ;
/*
   --------------------
   loop over the fronts
   --------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   count = nchild = 0 ;
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      nchild++ ;
      jthread = frontOwners[I] ;
      if ( mark[jthread] != J ) {
         mark[jthread] = J ;
         count++ ;
      }
   }
   counts[J] = nchild ;
   if ( flags != NULL ) {
      if ( count > 1 ) {
         flags[J] = 'Y' ;
      } else {
         flags[J] = 'N' ;
      }
   }
}
count = nchild = 0 ;
for ( J = ETree_root(frontmtx->frontETree) ; J != -1 ; J = sib[J] ) {
   nchild++ ;
   jthread = frontOwners[J] ;
   if ( mark[jthread] != J ) {
      mark[jthread] = J ;
      count++ ;
   }
}
counts[nfront] = nchild ;
if ( flags != NULL ) {
   if ( count > 1 ) {
      flags[nfront] = 'Y' ;
   } else {
      flags[nfront] = 'N' ;
   }
}
/*
   -----------------------------------------
   create and initialize the ChvList object
   -----------------------------------------
*/
postList = ChvList_new() ;
ChvList_init(postList, nfront+1, counts, lockflag, flags) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(mark) ;
IVfree(counts) ;
if ( flags != NULL ) {
   CVfree(flags) ;
}

return(postList) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   create, initialize and return a ChvList object
   to deal with aggregate chevrons
                  
   created -- 97jul03, cca
   -----------------------------------------------
*/
ChvList *
FrontMtx_aggregateList (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        lockflag
) {
char       *flags ;
ChvList   *aggList ;
int        count, ii, I, J, jthread, K, myid, nfront, nthread, size ;
int        *counts, *frontOwners, *head, *indices, *link, 
           *mark, *offsets, *vtxToFront ;
IVL        *symbfacIVL ;

/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || frontOwnersIV == NULL 
   || lockflag < 0 || lockflag > 2 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_aggregateList(%p,%p,%d)"
           "\n bad input\n", frontmtx, frontOwnersIV, lockflag) ;
   exit(-1) ;
}
symbfacIVL = frontmtx->symbfacIVL ;
vtxToFront = ETree_vtxToFront(frontmtx->frontETree) ;
IV_sizeAndEntries(frontOwnersIV, &nfront, &frontOwners) ;
nthread = 1 + IV_max(frontOwnersIV) ;
mark    = IVinit(nthread, -1) ;
head    = IVinit(nfront,  -1) ;
link    = IVinit(nfront,  -1) ;
offsets = IVinit(nfront,   0) ;
counts  = IVinit(nfront,   0) ;
if ( lockflag > 0 ) {
   flags = CVinit(nfront, 'N') ;
} else {
   flags = NULL ;
}
/*
   --------------------
   loop over the fronts
   --------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   myid = frontOwners[J] ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n front %d, owner %d", J, myid) ;
   fflush(stdout) ;
#endif
   mark[myid] = J ;
   count = 0 ;
/*
   ---------------------------------------------------
   loop over all descendent fronts that might update J
   ---------------------------------------------------
*/
   while ( (I = head[J]) != -1 ) {
      head[J] = link[I] ;
      jthread = frontOwners[I] ;
#if MYDEBUG > 0
      fprintf(stdout, "\n descendent front %d, owner %d", I,
jthread) ;
      fflush(stdout) ;
#endif
      if ( mark[jthread] != J ) {
/*
         --------------------------------
         expect an aggregate from jthread
         --------------------------------
*/
#if MYDEBUG > 0
         fprintf(stdout, ", incoming aggregate") ;
         fflush(stdout) ;
#endif
         mark[jthread] = J ;
         count++ ;
      }
/*
      --------------------------------------------------
      link front I to next ancestor that it does not own
      --------------------------------------------------
*/
      IVL_listAndSize(symbfacIVL, I, &size, &indices) ;
      for ( ii = offsets[I] ; ii < size ; ii++ ) {
         if (  (K = vtxToFront[indices[ii]]) > J
            && frontOwners[K] != jthread ) {
#if MYDEBUG > 0
            fprintf(stdout, ", link to %d", K) ;
            fflush(stdout) ;
#endif
            offsets[I] = ii ;
            link[I] = head[K] ;
            head[K] = I ;
            break ;
         }
      }
   }
/*
   -------------------------------------
   set the number of incoming aggregates
   -------------------------------------
*/
   counts[J] = count ;
#if MYDEBUG > 0
   fprintf(stdout, "\n counts[%d] = %d", J, counts[J]) ;
   fflush(stdout) ;
#endif
/*
   ---------------------------------------------------
   set the flags to see if the list needs to be locked
   ---------------------------------------------------
*/
   if ( flags != NULL ) {
      if ( count > 1 ) {
         flags[J] = 'Y' ;
      } else {
         flags[J] = 'N' ;
      }
#if MYDEBUG > 0
      fprintf(stdout, ", flags[%d] = %c", J, flags[J]) ;
      fflush(stdout) ;
#endif
   }
/*
   --------------------------------------------------
   link front J to next ancestor that it does not own
   --------------------------------------------------
*/
   IVL_listAndSize(symbfacIVL, J, &size, &indices) ;
   for ( ii = 0 ; ii < size ; ii++ ) {
      if (  (K = vtxToFront[indices[ii]]) > J && frontOwners[K] != myid
) {
#if MYDEBUG > 0
         fprintf(stdout, ", linking to %d", K) ;
         fflush(stdout) ;
#endif
         offsets[J] = ii ;
         link[J] = head[K] ;
         head[K] = J ;
         break ;
      }
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n counts") ;
IVfprintf(stdout, nfront, counts) ;
fflush(stdout) ;
#endif
/*
   -----------------------------------------
   create and initialize the ChvList object
   -----------------------------------------
*/
aggList = ChvList_new() ;
ChvList_init(aggList, nfront, counts, lockflag, flags) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(counts) ;
IVfree(head) ;
IVfree(link) ;
IVfree(offsets) ;
IVfree(mark) ;
if ( flags != NULL ) {
   CVfree(flags) ;
}
return(aggList) ; }

/*--------------------------------------------------------------------*/
