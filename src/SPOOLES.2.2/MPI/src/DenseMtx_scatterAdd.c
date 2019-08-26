/*  DenseMtx_scatterAdd.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
typedef struct _Msg Msg ;
struct _Msg {
   int           id    ;
   int           size  ;
   double        *base ;
   MPI_Request   req   ;
   Msg           *next ;
} ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- to scatter/add entries from a distributed X into Y.
      what rows of X to be sent to other processors are 
      found in sendIVL. what rows of Y to be received 
      from other processors are found in recvIVL.

   Y  -- on return, contains the rows specified by recvIVL.
      row indices of Y will be in ascending order.
   X  -- this processor's part of the distributed partitioned
      DenseMtx object. row indices of X are assumed to be 
      in ascending order.
   sendIVL -- list jproc contains the global ids of rows in X
      that need to be sent to processor jproc. note, lists are
      assumed to be in ascending order and are local with respect to X.
   recvIVL -- list jproc contains the global ids of rows in jproc's
      part of X that will to be sent to this processor. note,
      lists are assumed to be in ascending order and are local
      with respect to Y.

   created -- 98jul31, cca
   --------------------------------------------------------------------
*/
void
DenseMtx_MPI_scatterAddRows (
   DenseMtx   *Y,
   DenseMtx   *X,
   IVL        *sendIVL,
   IVL        *recvIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
double       *recvvec, *sendvec ;
int          flag, iirow, irow, jproc, jrow, myid, ncolX, ncolY, nproc,
             nrecv, nrowX, nrowY, nsend, nword, tag ;
int          *colindX, *colindY, *recvrowids, 
             *rowindX, *rowindY, *sendrowids ;
Msg          *msg, *nextmsg, *recvhead, *sendhead ;
MPI_Status   status ;
/*
   ---------------
   check the input
   ---------------
*/
if ( Y == NULL || X == NULL || sendIVL == NULL || recvIVL == NULL
     || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_scatterAddRows()"
           "\n bad input\n") ;
   exit(-1) ;
}
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( DENSEMTX_IS_REAL(X) ) {
   nword = 1 ;
} else if ( DENSEMTX_IS_COMPLEX(X) ) {
   nword = 2 ;
} else {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_scatterAddRows()"
           "\n X->type = %d\n", X->type) ;
   exit(-1) ;
}
DenseMtx_columnIndices(Y, &ncolY,  &colindY) ;
DenseMtx_rowIndices(Y, &nrowY, &rowindY) ;
DenseMtx_columnIndices(X, &ncolX,  &colindX) ;
DenseMtx_rowIndices(X, &nrowX, &rowindX) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n sendIVL ") ;
   IVL_writeForHumanEye(sendIVL, msgFile) ;
   fprintf(msgFile, "\n\n recvIVL ") ;
   IVL_writeForHumanEye(recvIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------
   scatter/add the internal rows
   -----------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n loading internal rows") ;
   fflush(msgFile) ;
}
IVL_listAndSize(sendIVL, myid, &nsend, &sendrowids) ;
IVL_listAndSize(recvIVL, myid, &nrecv, &recvrowids) ;
for ( iirow = 0 ; iirow < nsend ; iirow++ ) {
   irow = sendrowids[iirow] ;
   jrow = recvrowids[iirow] ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n irow %d, jrow %d", irow, jrow) ;
      fflush(msgFile) ;
   }
   DenseMtx_addRow(Y, jrow, X, irow) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after adding internal rows") ;
   DenseMtx_writeForHumanEye(Y, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   post the sends and the receives
   -------------------------------
*/
recvhead = sendhead = NULL ;
for ( jproc = 0 ; jproc < nproc ; jproc++ ) {
   if ( jproc != myid ) {
      IVL_listAndSize(sendIVL, jproc, &nsend, &sendrowids) ;
      IVL_listAndSize(recvIVL, jproc, &nrecv, &recvrowids) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n jproc %d, nsend %d, nrecv %d",
                 jproc, nsend, nrecv) ;
         fflush(msgFile) ;
      }
      if ( nsend > 0 ) {
/*
         -------------------------
         create the message object
         -------------------------
*/
         ALLOCATE(msg, struct _Msg, 1) ;
         msg->id   = jproc ;
         msg->size = nword * nsend * ncolY ;
         msg->base = sendvec = DVinit(msg->size, 0.0) ;
         msg->next = sendhead, sendhead = msg ;
         tag       = firsttag + myid*nproc + jproc ;
/*
         -----------------------------------
         fill the buffer with matrix entries
         -----------------------------------
*/
         for ( iirow = 0 ; iirow < nsend ; iirow++ ) {
            irow = sendrowids[iirow] ;
            DenseMtx_copyRowIntoVector(X, irow, sendvec) ;
            sendvec += nword*ncolY ;
         }
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n sendvec") ;
            DVfprintf(msgFile, msg->size, msg->base) ;
            fflush(msgFile) ;
         }
/*
         -------------
         post the send
         -------------
*/
         stats[0]++ ;
         stats[2] += msg->size * sizeof(double) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
                    "\n posting Isend to %d, size %d, tag %d",
                    jproc, msg->size, tag) ;
            fflush(msgFile) ;
         }
         MPI_Isend(msg->base, msg->size, MPI_DOUBLE, 
                   jproc, tag, comm, &msg->req) ;

      }
      if ( nrecv > 0 ) {
/*
         -------------------------
         create the message object
         -------------------------
*/
         ALLOCATE(msg, struct _Msg, 1) ;
         msg->id   = jproc ;
         msg->size = nword * nrecv * ncolY ;
         msg->base = (void *) DVinit(msg->size, 0.0) ;
         msg->next = recvhead, recvhead = msg ;
         tag       = firsttag + jproc*nproc + myid ;
/*
         ----------------
         post the receive
         ----------------
*/
         stats[1]++ ;
         stats[3] += msg->size * sizeof(double) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
                    "\n posting Irecv from %d, size %d, tag %d",
                    jproc, msg->size, tag) ;
            fflush(msgFile) ;
         }
         MPI_Irecv(msg->base, msg->size, MPI_DOUBLE, 
                   jproc, tag, comm, &msg->req) ;

      }
   }
}
/*
   -------------------------------------------
   loop while there are messages to receive or
   sent messages that have not been received
   -------------------------------------------
*/
while ( sendhead != NULL || recvhead != NULL ) {
   for ( msg = sendhead, sendhead = NULL ; 
         msg != NULL ; 
         msg = nextmsg ) {
      nextmsg = msg->next ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n msg %p to %d", msg, msg->id) ; 
         fflush(msgFile) ;
      }
      MPI_Test(&msg->req, &flag, &status) ;
      if ( flag == 1 ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, ", received") ;
            fflush(msgFile) ;
         }
         DVfree((double *) msg->base) ;
         FREE(msg) ;
      } else {
         msg->next = sendhead, sendhead = msg ;
      }
   }
   for ( msg = recvhead, recvhead = NULL ; 
         msg != NULL ; 
         msg = nextmsg ) {
      nextmsg = msg->next ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n msg %p from %d", msg, msg->id) ; 
         fflush(msgFile) ;
      }
      MPI_Test(&msg->req, &flag, &status) ;
      if ( flag == 1 ) {
         jproc = msg->id ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, ", received") ;
            fflush(msgFile) ;
         }
         IVL_listAndSize(recvIVL, jproc, &nrecv, &recvrowids) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n recvrowids") ;
            IVfprintf(msgFile, nrecv, recvrowids) ;
            fflush(msgFile) ;
         }
         recvvec = msg->base ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n recvvec") ;
            DVfprintf(msgFile, nword*nrecv*ncolY, recvvec) ;
            fflush(msgFile) ;
         }
         for ( iirow = 0 ; iirow < nrecv ; iirow++ ) {
            irow = recvrowids[iirow] ;
            DenseMtx_addVectorIntoRow(Y, irow, recvvec) ;
            recvvec += nword*ncolY ;
         }
         DVfree((double *) msg->base) ;
         FREE(msg) ;
      } else {
         msg->next = recvhead, recvhead = msg ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
