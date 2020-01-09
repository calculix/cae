/*  symbfac.c  */

#include "../spoolesMPI.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
typedef struct _Msg Msg ;
struct _Msg {
   int           id    ;
   int           size  ;
   void          *base ;
   MPI_Request   req   ;
   Msg           *next ;
} ;
static IVL * initSupportedIVL ( ETree *frontETree, IV *frontOwnersIV,
   int myid, int msglvl, FILE *msgFile ) ;
static void loadInternalIndices ( ETree *frontETree, InpMtx *inpmtxA,
   InpMtx *inpmtxB, IV *frontOwnersIV, int myid, IVL *symbfacIVL,
   int msglvl, FILE *msgFile ) ;
static void doCooperativeWork ( ETree *frontETree, IV *frontOwnersIV, 
   IVL *symbfacIVL, int stats[], int msglvl, 
   FILE *msgFile, int firsttag, MPI_Comm comm ) ;
static int mergeIndices ( int sizeJ, int indJ[], 
                          int sizeK, int indK[], int nleft ) ;
static Msg * Msg_new ( void ) ;
static void Msg_setDefaultFields ( Msg *msg ) ;
static void Msg_clearData ( Msg *msg ) ;
static void Msg_free ( Msg *msg ) ;
static Msg * wakeupFront ( int J, int myid, int owners[], 
   char supportTable[], ETree *frontETree, IVL *symbfacIVL,
   int firsttag, int stats[], MPI_Comm comm, 
   int msglvl, FILE *msgFile ) ;
static Msg * checkRecvMessages ( int J, Msg *firstmsg, int nLeft[],
   int nodwghts[], IVL *symbfacIVL, int msglvl, FILE *msgFile ) ;
static Msg * sendMessages ( int J, int myid, int owners[], int nfront, 
   int nproc, char supportTable[], int par[], IVL *symbfacIVL,
   Msg *firstMsgSent, int firsttag, int stats[], MPI_Comm comm, 
   int msglvl, FILE *msgFile ) ;
static Msg * checkSendMessages ( Msg *firstMsgSent,
   int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   perform the symbolic factorization in parallel

   input --
      frontETree    -- front tree for the factorization
      frontOwnersIV -- map from fronts to owning process
      inpmtx        -- matrix object
      firsttag      -- first tag to be used for messages,
                       will use tag, ..., tag + nfront - 1
      msglvl        -- message level
      msgFile       -- message file
      comm          -- MPI communicator

   return value --
      symbfacIVL -- contains front indices for supported fronts

   created -- 98may20, cca
   ------------------------------------------------------------
*/
IVL *
SymbFac_MPI_initFromInpMtx (
   ETree      *frontETree,
   IV         *frontOwnersIV,
   InpMtx     *inpmtx,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int   lasttag, myid, nproc, tagbound; 
IVL   *symbfacIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || frontOwnersIV == NULL
   || inpmtx == NULL || stats == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in SymbFac_MPI_initFromInpMtx()"
           "\n comm = %p, frontETree = %p, frontOwnersIV = %p"
           "\n inpmtx = %p, firsttag = %d, msglvl = %d, msgFile = %p"
           "\n bad input\n", comm, frontETree, frontOwnersIV, inpmtx,
           firsttag, msglvl, msgFile) ;
   exit(-1) ;
}
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n myid = %d, nproc = %d", myid, nproc) ;
   fflush(msgFile) ;
}
if ( firsttag < 0 ) {
   fprintf(stderr, "\n fatal error in SymbFac_MPI_initFromInpMtx()"
           "\n firsttag = %d\n", firsttag) ;
   exit(-1) ;
}
lasttag = firsttag + frontETree->nfront ;
if ( lasttag > (tagbound = maxTagMPI(comm)) ) {
   fprintf(stderr, "\n fatal error in SymbFac_MPI_initFromInpMtx()"
           "\n lasttag = %d, tag_bound = %d", lasttag, tagbound) ;
   exit(-1) ;
}
/*
   --------------------------------------
   first, set up the supported IVL object
   --------------------------------------
*/
symbfacIVL = initSupportedIVL(frontETree, frontOwnersIV, 
                              myid, msglvl, msgFile) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n local supported IVL after initialization") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------------------
   second, fill the index lists of the owned fronts with the internal
   vertices and boundary vertices available from the matrix pencil
   ------------------------------------------------------------------
*/
loadInternalIndices(frontETree, inpmtx, NULL, frontOwnersIV, 
                    myid, symbfacIVL, msglvl, msgFile ) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n after loading internal indices") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   third, cooperate among the processors to generate 
   the rest of the supported symbolic factorization
   -------------------------------------------------
*/
doCooperativeWork(frontETree, frontOwnersIV, symbfacIVL, 
                  stats, msglvl, msgFile, firsttag, comm) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n final local supported IVL ") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
return(symbfacIVL) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   perform the symbolic factorization in parallel

   input --
      frontETree    -- front tree for the factorization
      frontOwnersIV -- map from fronts to owning process
      pencil        -- matrix pencil object
      firsttag      -- first tag to be used for messages,
                       will use tag, ..., tag + nfront - 1
      msglvl        -- message level
      msgFile       -- message file
      comm          -- MPI communicator

   return value --
      symbfacIVL -- contains front indices for supported fronts

   created -- 98may20, cca
   ------------------------------------------------------------
*/
IVL *
SymbFac_MPI_initFromPencil (
   ETree      *frontETree,
   IV         *frontOwnersIV,
   Pencil     *pencil,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int   lasttag, myid, nproc, tagbound; 
IVL   *symbfacIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || frontOwnersIV == NULL
   || pencil == NULL || stats == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in SymbFac_MPI_initFromPencil()"
           "\n comm = %p, frontETree = %p, frontOwnersIV = %p"
           "\n pencil = %p, firsttag = %d, msglvl = %d, msgFile = %p"
           "\n bad input\n", comm, frontETree, frontOwnersIV, pencil,
           firsttag, msglvl, msgFile) ;
   exit(-1) ;
}
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n myid = %d, nproc = %d", myid, nproc) ;
   fflush(msgFile) ;
}
if ( firsttag < 0 ) {
   fprintf(stderr, "\n fatal error in SymbFac_MPI_initFromPencil()"
           "\n firsttag = %d\n", firsttag) ;
   exit(-1) ;
}
lasttag = firsttag + frontETree->nfront ;
if ( lasttag > (tagbound = maxTagMPI(comm)) ) {
   fprintf(stderr, "\n fatal error in SymbFac_MPI_initFromPencil()"
           "\n lasttag = %d, tag_bound = %d", lasttag, tagbound) ;
   exit(-1) ;
}
/*
   --------------------------------------
   first, set up the supported IVL object
   --------------------------------------
*/
symbfacIVL = initSupportedIVL(frontETree, frontOwnersIV, 
                              myid, msglvl, msgFile) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n local supported IVL after initialization") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------------------
   second, fill the index lists of the owned fronts with the internal
   vertices and boundary vertices available from the matrix pencil
   ------------------------------------------------------------------
*/
loadInternalIndices(frontETree, pencil->inpmtxA, pencil->inpmtxB, 
                    frontOwnersIV, myid, symbfacIVL, msglvl, msgFile ) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n after loading internal indices") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   third, cooperate among the processors to generate 
   the rest of the supported symbolic factorization
   -------------------------------------------------
*/
doCooperativeWork(frontETree, frontOwnersIV, symbfacIVL, 
                  stats, msglvl, msgFile, firsttag, comm) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n final local supported IVL ") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
return(symbfacIVL) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose -- initialize the symbolic factorization IVL for process myid

   created -- 98may20
   ---------------------------------------------------------------------
*/
static IVL *
initSupportedIVL (
   ETree   *frontETree,
   IV      *frontOwnersIV,
   int     myid,
   int     msglvl,
   FILE    *msgFile
) {
char   *frontSupported ;
int    J, K, nfront ;
int    *bndwghts, *nodwghts, *owners, *par, *sizes ;
IVL    *symbfacIVL ;

nfront   = ETree_nfront(frontETree) ;
nodwghts = ETree_nodwghts(frontETree),
bndwghts = ETree_bndwghts(frontETree),
par      = ETree_par(frontETree) ;
owners   = IV_entries(frontOwnersIV) ;
/*
   -----------------------------------------------------------------
   create the frontSupported[] vector,
   frontSupported[J] = 'Y' --> this process supports front J
   frontSupported[J] = 'N' --> this process does not support front J
   -----------------------------------------------------------------
*/
frontSupported = CVinit(nfront, 'N') ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( owners[J] == myid ) {
      K = J ;
      while ( K != -1 && frontSupported[K] == 'N' ) {
         frontSupported[K] = 'Y' ;
         K = par[K] ;
      }
   }
}
/*
   -----------------------------------------------
   initialize the local symbolic factorization IVL
   -----------------------------------------------
*/
symbfacIVL = IVL_new() ;
sizes = IVinit(nfront, 0) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( frontSupported[J] == 'Y' ) {
      sizes[J] = nodwghts[J] + bndwghts[J] ;
   }
}
IVL_init3(symbfacIVL, IVL_CHUNKED, nfront, sizes) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(sizes)  ;
CVfree(frontSupported) ;

return(symbfacIVL) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   perform the symbolic factorization in parallel

   input --
      frontETree    -- front tree for the factorization
      frontOwnersIV -- map from fronts to owning process
      pencil        -- matrix pencil
      firsttag      -- first tag to be used for messages,
                       will use tag, ..., tag + nfront - 1
      msglvl        -- message level
      msgFile       -- message file
      comm          -- MPI communicator

   return value --
      symbfacIVL -- contains front indices for supported fronts

   created -- 98may20, cca
   --------------------------------------------------------------------
*/
static void
doCooperativeWork (
   ETree      *frontETree,
   IV         *frontOwnersIV,
   IVL        *symbfacIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
char   *frontIsSupported, *status, *supportTable ;
Ideq   *dequeue ;
int    count, ii, iproc, J, K, myid, nDJ, nfront, npath, nproc, 
       nvtx, sizeJ, sizeK ;
int    *bndwghts, *frontOwners, *indicesJ, *indicesK, 
       *nactiveChild, *nLeft, *nodwghts, *par, *vtxToFront ;
Msg    *firstMsgSent ;
Msg    **p_msg ;
/*
   -------------------------------
   extract pointers and dimensions
   -------------------------------
*/
nfront     = ETree_nfront(frontETree) ;
nvtx       = ETree_nvtx(frontETree) ;
nodwghts   = ETree_nodwghts(frontETree),
bndwghts   = ETree_bndwghts(frontETree),
vtxToFront = ETree_vtxToFront(frontETree),
par        = ETree_par(frontETree) ;
IV_sizeAndEntries(frontOwnersIV, &nfront, &frontOwners) ;
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n myid = %d, nproc = %d", myid, nproc) ;
   fflush(msgFile) ;
}
/*
   --------------------------
   generate the support table
   --------------------------
*/
supportTable = CVinit(nfront*nproc, 'N') ;
for ( J = 0 ; J < nfront ; J++ ) {
   iproc = frontOwners[J] ;
   frontIsSupported = supportTable + iproc*nfront ;
   for ( K = J ; K != -1 && frontIsSupported[K] == 'N' ; K = par[K] ) {
      frontIsSupported[K] = 'Y' ;
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n\n supportTable") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      fprintf(msgFile, "\n") ;
      for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
         fprintf(msgFile, " %c", supportTable[J + iproc*nfront]) ;
      }
   }
   fflush(msgFile) ;
}
/*
   ---------------------------------------------------------
   initialize the status[] and nactiveChild[] vectors and 
   load the dequeue with the leaves of the supported subtree
   ---------------------------------------------------------
*/
frontIsSupported = supportTable + myid*nfront ;
status = CVinit(nfront, 'F') ;
for ( J = npath = 0 ; J < nfront ; J++ ) {
   if ( status[J] == 'F' && frontIsSupported[J] == 'Y' ) {
      npath++ ;
      for ( K = J ; K != -1 && status[K] == 'F' ; K = par[K] ) {
         status[K] = 'W' ;
      }
   }
}
dequeue = Ideq_new() ;
Ideq_resize(dequeue, npath) ;
nactiveChild = IVinit(nfront, 0) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( status[J] == 'W' ) {
      if ( (K = par[J]) != -1 ) {
         nactiveChild[K]++ ;
      }
   }
}
for ( J = 0 ; J < nfront ; J++ ) {
   if ( frontOwners[J] == myid && nactiveChild[J] == 0 ) {
      Ideq_insertAtTail(dequeue, J) ;
   }
}
/*
   ---------------------------
   generate the nLeft[] vector
   ---------------------------
*/
nLeft = IVinit(nfront, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( frontOwners[J] == myid ) {
      IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
      for ( ii = sizeJ - 1, count = 0 ; ii >= 0 ; ii-- ) {
         if ( indicesJ[ii] != -1 ) {
            break ;
         }
         count++ ;
      }
      nLeft[J] = count ;
   }
}
/*
   -------------------------------------------
   generate the vector of pointers to messages
   -------------------------------------------
*/
firstMsgSent = NULL ;
ALLOCATE(p_msg, struct _Msg *, nfront) ;
for ( J = 0 ; J < nfront ; J++ ) {
   p_msg[J] = NULL ;
}
/*
   ----------------------------------
   loop while the dequeue is nonempty
   ----------------------------------
*/
while ( (J = Ideq_removeFromHead(dequeue)) != -1 ) {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n ### popping %d from dequeue", J) ;
      fflush(msgFile) ;
   }
   if ( status[J] == 'W' ) {
/*
      --------------------------------------------
      wake up the front, post its receive messages
      --------------------------------------------
*/
      p_msg[J] = wakeupFront(J, myid, frontOwners, supportTable,
                             frontETree, symbfacIVL, firsttag, 
                             stats, comm, msglvl, msgFile) ;
      status[J] = 'A' ;
   }
   if ( p_msg[J] != NULL ) {
/*
      -----------------------------
      try to receive messages for J
      -----------------------------
*/
      p_msg[J] = checkRecvMessages(J, p_msg[J], nLeft, nodwghts,  
                                   symbfacIVL, msglvl, msgFile) ;
   }
   if ( p_msg[J] == NULL ) {
      if ( msglvl > 3 ) {
         IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
         fprintf(msgFile, 
                 "\n adjacency list for J = %d, nJ = %d, bndJ = %d",
                 J, nodwghts[J], bndwghts[J]) ;
         IVfprintf(msgFile, sizeJ, indicesJ) ;
         fflush(msgFile) ;
      }
/*
      -------------------------------------------
      all messages received:
      check for the length of the indices vector,
      set the boundary size of the front
      ------------------------------------------
*/
      status[J] = 'D' ;
      IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
      for ( ii = 0 ; ii < sizeJ ; ii++ ) {
         if ( indicesJ[ii] == -1 ) {
            break ;
         }
      }
      if ( ii != sizeJ ) {
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n WHOA, front J %d, sizeJ %d", J, sizeJ);
            fprintf(msgFile, "\n old bndwght = %d, new bndwght = %d",
                    bndwghts[J], ii - nodwghts[J]) ;
            IVfprintf(msgFile, sizeJ, indicesJ) ;
            fflush(msgFile) ;
         }
         bndwghts[J] = ii - nodwghts[J] ;
         IVL_setList(symbfacIVL, J, ii, indicesJ) ;
      }
      if ( frontOwners[J] == myid ) {
/*       
         ------------------------------------------
         front is owned, send messages if necessary
         ------------------------------------------
*/
         firstMsgSent = sendMessages(J, myid, frontOwners, nfront, 
                                     nproc, supportTable, par, 
                                     symbfacIVL, firstMsgSent, firsttag,
                                     stats, comm, msglvl, msgFile) ;
      }
      if ( (K = par[J]) != -1 ) {
         if ( frontOwners[K] == myid && nLeft[K] > 0 ) {
/*
            ------------------------------
            merge bnd{J} into K cup bnd{K}
            ------------------------------
*/
            nDJ = nodwghts[J] ;
            IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
            IVL_listAndSize(symbfacIVL, K, &sizeK, &indicesK) ;
            if ( msglvl > 3 ) {
               fprintf(msgFile, "\n merging J = %d into K = %d", J, K) ;
               fprintf(msgFile, ", nleft = %d", nLeft[K]) ;
               fprintf(msgFile, "\n before, boundary indices of J") ;
               IVfprintf(msgFile, sizeJ - nDJ, indicesJ + nDJ) ;
               fprintf(msgFile, "\n before, indices of K") ;
               IVfprintf(msgFile, sizeK, indicesK) ;
            }
            nLeft[K] = mergeIndices(sizeJ - nDJ, indicesJ + nDJ, 
                                    sizeK, indicesK, nLeft[K]) ;
            if ( msglvl > 3 ) {
               fprintf(msgFile, 
                       "\n after, indices of K, nleft = %d", nLeft[K]) ;
               IVfprintf(msgFile, sizeK, indicesK) ;
            }
         }
         if ( --nactiveChild[K] == 0 ) {
/*
            ----------------------
            put K on head of queue
            ----------------------
*/
            Ideq_insertAtHead(dequeue, K) ;
         }
      }
   } else {
/*
      ----------------------
      put J on tail of queue
      ----------------------
*/
      Ideq_insertAtTail(dequeue, J) ;
   }
   firstMsgSent = checkSendMessages(firstMsgSent, msglvl, msgFile) ;
}
/*
   ------------------------------------------
   check for all sent messages to be received
   ------------------------------------------
*/
do {
   firstMsgSent = checkSendMessages(firstMsgSent, msglvl, msgFile) ;
} while ( firstMsgSent != NULL ) ;
/*
   --------------------------------
   post-process the adjacency lists
   --------------------------------
*/
/*
if ( msglvl > 2 ) {
      fprintf(msgFile, "\n bndwghts") ;
   IVfprintf(msgFile, nfront, bndwghts) ;
   fprintf(msgFile, "\n nodwghts") ;
   IVfprintf(msgFile, nfront, nodwghts) ;
}
frontIsSupported = supportTable + myid*nfront ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( frontIsSupported[J] == 'Y' ) {
      IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
      for ( ii = 0 ; ii < sizeJ ; ii++ ) {
         if ( indicesJ[ii] == -1 ) {
            break ;
         }
      }
      if ( ii != sizeJ ) {
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n WHOA, front J %d, sizeJ %d", J, sizeJ);
            fprintf(msgFile, "\n old bndwght = %d, new bndwght = %d",
                    bndwghts[J], ii - nodwghts[J]) ;
            IVfprintf(msgFile, sizeJ, indicesJ) ;
            fflush(msgFile) ;
         }
         bndwghts[J] = ii - nodwghts[J] ;
         IVL_setList(symbfacIVL, J, ii, indicesJ) ;
      }
   }
}
*/
/*
   ------------------------
   free the working storage
   ------------------------
*/
CVfree(supportTable) ;
CVfree(status) ;
IVfree(nactiveChild) ;
IVfree(nLeft) ;
Ideq_free(dequeue) ;
FREE(p_msg) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   constructor

   created -- 98may20, cca
   -----------------------
*/
static Msg *
Msg_new (
   void
) {
Msg   *msg ;

ALLOCATE(msg, struct _Msg, 1) ;
Msg_setDefaultFields(msg) ;

return(msg) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   set the default fields

   created -- 98may20, cca
   -----------------------
*/
static void
Msg_setDefaultFields (
   Msg   *msg
) {
msg->id   =  -1  ;
msg->size =   0  ;
msg->base = NULL ;
msg->next = NULL ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   clear the data

   created -- 98may20, cca
   -----------------------
*/
static void
Msg_clearData (
   Msg   *msg
) {

Msg_setDefaultFields(msg) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   free the object

   created -- 98may20, cca
   -----------------------
*/
static void
Msg_free (
   Msg   *msg
) {

Msg_clearData(msg) ;
FREE(msg) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   when front J is first examined, post Irecv's for its messages
  
   created -- 98may20, cca
   -------------------------------------------------------------
*/
static Msg *
wakeupFront (
   int        J,
   int        myid,
   int        owners[],
   char       supportTable[],
   ETree      *frontETree,
   IVL        *symbfacIVL,
   int        firsttag,
   int        stats[],
   MPI_Comm   comm,
   int        msglvl,
   FILE       *msgFile
) {
Msg   *first, *msg ;
int   I, nfront ;
int   *bndwghts, *fch, *indices, *nodwghts, *sib ;

nfront = ETree_nfront(frontETree) ;
first = NULL ;
if ( owners[J] == myid ) {
/*
   ---------------------------------------------------
   owned front, post messages for unsupported children
   ---------------------------------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n waking up owned front %d", J) ;
      fflush(msgFile) ;
   }
   fch      = ETree_fch(frontETree) ;
   sib      = ETree_sib(frontETree) ;
   nodwghts = ETree_nodwghts(frontETree) ;
   bndwghts = ETree_bndwghts(frontETree) ;
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      if ( supportTable[I + myid*nfront] != 'Y' ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n    unsupported child %d", I) ;
            fflush(msgFile) ;
         }
/*
         --------------------------------------------------------
         create a message, allocate temporary storage for indices
         --------------------------------------------------------
*/
         msg = Msg_new() ;
         msg->id = I ;
         msg->size = nodwghts[I] + bndwghts[I] ;
         msg->base = (void *) IVinit(msg->size, -1) ;
         msg->next = first ;
         first     = msg ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
              "\n    1. posting Irecv, msg %p, source = %d, tag = %d",
              msg, owners[I], firsttag + I) ;
            fflush(msgFile) ;
         }
         stats[1]++ ;
         stats[3] += msg->size * sizeof(int) ;
         MPI_Irecv(msg->base, msg->size, MPI_INT, owners[I], 
                   firsttag + I, comm, &msg->req) ;
      }
   }
} else if ( supportTable[J + myid*nfront] == 'Y' ) {
/*
   --------------------------------------------------------
   supported but not owned. create a message but point 
   into storage for the indices in the symbolic IVL object.
   --------------------------------------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n waking up supported front %d", J) ;
      fflush(msgFile) ;
   }
   msg = Msg_new() ;
   msg->id = J ;
   IVL_listAndSize(symbfacIVL, J, &msg->size, &indices) ;
   msg->base = (void *) indices ;
   msg->next = first ;
   first     = msg ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
        "\n    1. posting Irecv, msg %p, source = %d, tag = %d",
        msg, owners[J], firsttag + J) ;
      fflush(msgFile) ;
   }
   stats[1]++ ;
   stats[3] += msg->size * sizeof(int) ;
   MPI_Irecv(msg->base, msg->size, MPI_INT, owners[J], firsttag + J,
             comm, &msg->req) ;
}
return(first) ; }
   
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   check the messages waiting to be received for front J

   return value -- pointer to first message still to be received

   created -- 98may20, cca
   -------------------------------------------------------------
*/
static Msg *
checkRecvMessages (
   int    J,
   Msg    *firstmsg,
   int    nLeft[],
   int    nodwghts[],
   IVL    *symbfacIVL,
   int    msglvl,
   FILE   *msgFile
) {
int          flag, I, nDI, sizeI, sizeJ ;
int          *indicesI, *indicesJ ;
MPI_Status   status ;
Msg          *msg, *nextmsg ;

if ( msglvl > 2 ) {
   fprintf(msgFile, "\n checking for received messages") ;
   fflush(msgFile) ;
}
for ( msg = firstmsg, firstmsg = NULL ; msg != NULL ; msg = nextmsg ) {
   nextmsg   = msg->next ;
   msg->next = NULL ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n    checking msg %p : id %d, size %d",
              msg, msg->id, msg->size) ;
      fflush(msgFile) ;
   }
   MPI_Test(&msg->req, &flag, &status) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, ", flag = %d", flag) ;
      fflush(msgFile) ;
   }
   if ( flag != 0 ) {
      if ( msg->id == J ) {
/*
         --------------------------------
         J is supported but not owned, 
         indices are in the symbolic IVL, 
         free the message
         --------------------------------
*/
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n    supported front %d", J) ;
            fflush(msgFile) ;
         }
         Msg_free(msg) ;
      } else {
/*
         ---------------------------------------------------
         message is from an unsupported child
         merge indices, then free indices, then free message
         ---------------------------------------------------
*/
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
                    "\n    unsupported child front %d", msg->id) ;
            fflush(msgFile) ;
         }
         if ( nLeft[J] > 0 ) {
            I        = msg->id ;
            nDI      = nodwghts[I] ;
            sizeI    = msg->size - nDI ;
            indicesI = (int *) msg->base + nDI ;
/*
            sizeI    = msg->size ;
            indicesI = (int *) msg->base ;
*/
            IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
            if ( msglvl > 3 ) {
               fprintf(msgFile, "\n merging I = %d into J = %d", I, J) ;
               fprintf(msgFile, ", nleft = %d", nLeft[J]) ;
               fprintf(msgFile, "\n before, boundary indices of I") ;
               IVfprintf(msgFile, sizeI, indicesI) ;
               fprintf(msgFile, "\n before, indices of J") ;
               IVfprintf(msgFile, sizeJ, indicesJ) ;
            }
            nLeft[J] = mergeIndices(sizeI, indicesI, 
                                    sizeJ, indicesJ, nLeft[J]) ;
         }
         if ( msglvl > 3 ) {
            fprintf(msgFile, 
                    "\n after, indices of J, nleft = %d", nLeft[J]) ;
            IVfprintf(msgFile, sizeJ, indicesJ) ;
         }
         IVfree(msg->base) ;
         Msg_free(msg) ;
      }
   } else {
/*
      ----------------------------
      keep the message on the list
      ----------------------------
*/
      msg->next = firstmsg ;
      firstmsg  = msg ;
   }
}
return(firstmsg) ; }

/*--------------------------------------------------------------------*/
/*
*/
static Msg *
sendMessages (
   int        J, 
   int        myid, 
   int        owners[], 
   int        nfront, 
   int        nproc, 
   char       supportTable[], 
   int        par[], 
   IVL        *symbfacIVL,
   Msg        *firstMsgSent,
   int        firsttag, 
   int        stats[], 
   MPI_Comm   comm,
   int        msglvl, 
   FILE       *msgFile
) {
int   destination, K, sizeJ ;
int   *indicesJ ;
Msg   *msg ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n ready to send messages") ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   get size and pointer to indices
   -------------------------------
*/
IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
/*
   -----------------------------------
   send indices to supported processes
   -----------------------------------
*/
for ( destination = 0 ; destination < nproc ; destination++ ) {
   if (  destination != myid 
      && supportTable[J + destination*nfront] == 'Y' ) {
      msg          = Msg_new() ;
      msg->id      = J ;
      msg->size    = sizeJ ;
      msg->base    = (void *) indicesJ ;
      msg->next    = firstMsgSent ;
      firstMsgSent = msg ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
           "\n    1. posting Isend, msg %p, destination = %d, tag = %d",
           msg, destination, firsttag + J) ;
         fflush(msgFile) ;
      }
      stats[0]++ ;
      stats[2] += msg->size * sizeof(int) ;
      MPI_Isend(msg->base, msg->size, MPI_INT, destination, 
                firsttag + J, comm, &msg->req) ;
   }
}
if ( (K = par[J]) != -1 && supportTable[J + owners[K]*nfront] != 'Y' ) {
/*
   ----------------------------------------------------
   send indices to unsupported process that owns parent
   ----------------------------------------------------
*/
   msg          = Msg_new() ;
   msg->id      = J ;
   msg->size    = sizeJ ;
   msg->base    = (void *) indicesJ ;
   msg->next    = firstMsgSent ;
   firstMsgSent = msg ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
        "\n    2. posting Isend, msg %p, destination = %d, tag = %d",
        msg, owners[K], firsttag + J) ;
      fflush(msgFile) ;
   }
   stats[0]++ ;
   stats[2] += msg->size * sizeof(int) ;
   MPI_Isend(msg->base, msg->size, MPI_INT, owners[K], 
             firsttag + J, comm, &msg->req) ;
}
return(firstMsgSent) ; }
  
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   merge indices in indJ[sizeJ] into indK[sizeK],
   there are nleft free locations at the end of indK[].
   both indJ[sizeJ] and indK[sizeK] are in ascending order.
   return value is the number of present free locations in indK[].
 
   created -- 98may20, cca
   ---------------------------------------------------------------
*/
static int
mergeIndices (
   int   sizeJ,
   int   indJ[],
   int   sizeK,
   int   indK[],
   int   nleft
) {
int   firstK, ii, jj, jlast, kk, v ;
 
firstK = indK[0] ;
kk = jlast = sizeK - nleft ;
ii = 0, jj = 0 ; 
while ( ii < sizeJ && jj < jlast ) {
   v = indJ[ii] ;
   if ( v < firstK ) {
      ii++ ;
   } else if ( v < indK[jj] ) {
      indK[kk++] = v ;
      ii++ ;
   } else if ( v == indK[jj] ) {
      ii++, jj++ ;
   } else {
      jj++ ;
   }
}
for ( ; ii < sizeJ ; ii++ ) {
   v = indJ[ii] ;
   if ( v >= 0 ) {
      indK[kk++] = indJ[ii] ;
   }
}
if ( kk > jlast ) {
   IVqsortUp(kk, indK) ;
}
return(sizeK - kk) ; }
 
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   load the owned index lists with information from the InpMtx object
 
   created -- 98may20, cca
   ------------------------------------------------------------------
*/
static void
loadInternalIndices (
   ETree    *frontETree,
   InpMtx   *inpmtxA,
   InpMtx   *inpmtxB,
   IV       *frontOwnersIV,
   int      myid,
   IVL      *symbfacIVL,
   int      msglvl,
   FILE     *msgFile
) {
int   count, ii, J, nfront, nvtx, off, sizeJ, v, vsize, w ;
int   *frontOwners, *head, *indicesJ, *link, *mark, *vind, *vtxToFront ;

nfront      = ETree_nfront(frontETree) ;
nvtx        = ETree_nvtx(frontETree) ;
vtxToFront  = ETree_vtxToFront(frontETree) ;
frontOwners = IV_entries(frontOwnersIV) ;
/*
   -------------------------------------------------
   set up the lists of vertices for each owned front
   -------------------------------------------------
*/
head = IVinit(nfront, -1) ;
link = IVinit(nvtx,   -1) ;
mark = IVinit(nvtx,  -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   J = vtxToFront[v] ;
   if ( frontOwners[J] == myid ) {
      link[v] = head[J] ;
      head[J] = v ;
   }
}
/*
   --------------------------------------------------------
   load the internal indices into the owned adjacency lists
   --------------------------------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   if ( frontOwners[J] == myid ) {
      IVL_listAndSize(symbfacIVL, J, &sizeJ, &indicesJ) ;
      count = 0 ;
/*
      ---------------------
      load internal indices
      ---------------------
*/
      for ( v = head[J] ; v != -1 ; v = link[v] ) {
         mark[v] = J ;
         indicesJ[count++] = v ;
      }
/*
      -----------------------------------------
      load indices from the first InpMtx object
      -----------------------------------------
*/
      for ( v = head[J] ; v != -1 && count < sizeJ ; v = link[v] ) {
         InpMtx_vector(inpmtxA, v, &vsize, &vind) ;
         for ( ii = 0 ; ii < vsize && count < sizeJ ; ii++ ) {
            off = vind[ii] ;
            if ( off >= 0 ) {
               w = v + off ;
            } else {
               w = v - off ;
            }
            if ( vtxToFront[w] < J ) {
               fprintf(msgFile,
                       "\n error, J = %d, w = %d, map[%d] = %d",
                       J, w, w, vtxToFront[w]) ;
               exit(-1) ;
            }
            if ( mark[w] != J ) {
               mark[w] = J ;
               indicesJ[count++] = w ;
            }
         }
      }
      if ( inpmtxB != NULL ) {
/*
         -------------------
         load indices from B
         -------------------
*/
         for ( v = head[J] ; v != -1 && count < sizeJ ; v = link[v] ) {
            InpMtx_vector(inpmtxB, v, &vsize, &vind) ;
            for ( ii = 0 ; ii < vsize && count < sizeJ ; ii++ ) {
               off = vind[ii] ;
               if ( off >= 0 ) {
                  w = v + off ;
               } else {
                  w = v - off ;
               }
               if ( vtxToFront[w] < J ) {
                  fprintf(msgFile,
                          "\n error, J = %d, w = %d, map[%d] = %d",
                          J, w, w, vtxToFront[w]) ;
                  exit(-1) ;
               }
               if ( mark[w] != J ) {
                  mark[w] = J ;
                  indicesJ[count++] = w ;
               }
            }
         }
      }
/*
      -----------------------------------
      sort the indices in ascending order 
      -----------------------------------
*/
      IVqsortUp(count, indicesJ) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(head) ;
IVfree(link) ;
IVfree(mark) ;
 
return ; }
 
/*--------------------------------------------------------------------*/
/*
*/
static Msg * 
checkSendMessages (
   Msg    *firstMsgSent,
   int    msglvl,
   FILE   *msgFile
) {
int          flag ;
MPI_Status   status ;
Msg          *msg, *nextmsg ;

if ( firstMsgSent != NULL ) {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n checking for sent messages") ;
      fflush(msgFile) ;
   }
   for ( msg = firstMsgSent, firstMsgSent = NULL ; 
         msg != NULL ; 
         msg = nextmsg ) {
      nextmsg = msg->next ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n    checking msg %p : id %d, size %d",
                 msg, msg->id, msg->size) ;
         fflush(msgFile) ;
      }
      MPI_Test(&msg->req, &flag, &status) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, ", flag = %d", flag) ;
         fflush(msgFile) ;
      }
      if ( flag != 0 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, ", free'ing object") ;
            fflush(msgFile) ;
         }
         Msg_free(msg) ;
      } else {
         msg->next = firstMsgSent ;
         firstMsgSent = msg ;
      }
   }
}

return(firstMsgSent) ; }

/*--------------------------------------------------------------------*/
