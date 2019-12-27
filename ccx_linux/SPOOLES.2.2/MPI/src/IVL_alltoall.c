/*  IVL_alltoall.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   this method is used during the setup for matrix-vector multiplies.
   each processor has computed the vertices it needs from other
   processors, these lists are contained in sendIVL. on return,
   recvIVL contains the lists of vertices this processor must send
   to all others.

   sendIVL -- on input, list[q] contains the vertices needed by
              this processor that are owned by q
   recvIVL -- on output, list[q] contains the vertices owned by
              this processor that are needed by q. note, if NULL
              on input, a new IVL object is allocated
   stats[] -- statistics vector
     stats[0] -- contains # of sends
     stats[1] -- contains # of receives
     stats[2] -- contains # of bytes sent
     stats[3] -- contains # of bytes received
   firsttag -- first tag for messages,
     tags in range [firsttag, firsttag+nproc-1] are used

   return value -- recvIVL

   created -- 98jul26, cca
   ------------------------------------------------------------------
*/
IVL *
IVL_MPI_alltoall (
   IVL        *sendIVL,
   IVL        *recvIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int          count, destination, lasttag, left, myid, nproc, offset, q, 
             recvcount, right, sendcount, source, tag, tagbound ;
int          *incounts, *outcounts, *recvvec, *sendvec ;
MPI_Status   status ;
/*
   ---------------
   check the input
   ---------------
*/
if ( sendIVL == NULL || stats == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(msgFile, "\n fatal error in IVL_MPI_alltoall()"
           "\n bad input\n") ;
   exit(-1) ;
}
/*
   ---------------------------------------
   get id of self and number of processors
   ---------------------------------------
*/
MPI_Comm_rank(comm, &myid)  ;
MPI_Comm_size(comm, &nproc) ;
if ( sendIVL->nlist != nproc ) {
   fprintf(msgFile, "\n fatal error in IVL_MPI_alltoall()"
           "\n sendIVL: nproc = %d, nlist = %d\n", 
           nproc, sendIVL->nlist) ;
   exit(-1) ;
}
lasttag = firsttag + nproc ;
if ( lasttag > (tagbound = maxTagMPI(comm)) ) {
   fprintf(stderr, "\n fatal error in IVL_MPI_alltoall()"
           "\n lasttag = %d, tag_bound = %d", lasttag, tagbound) ;
   exit(-1) ;
}

if ( recvIVL == NULL ) {
   recvIVL = IVL_new() ;
} else {
   IVL_clearData(recvIVL) ;
}
IVL_init1(recvIVL, IVL_CHUNKED, nproc) ;
/*
   ------------------------------------------
   outcounts[] is sendIVL->sizes[]
   incounts[] will be recvIVL->sizes[]
   fill incounts via a call to MPI_Alltoall()
   and then initialize the recvIVL lists.
   ------------------------------------------
*/
outcounts = sendIVL->sizes ;
incounts  = IVinit(nproc, 0) ;
MPI_Alltoall((void *) outcounts, 1, MPI_INT,
             (void *) incounts,  1, MPI_INT, comm) ;
for ( q = 0 ; q < nproc ; q++ ) {
   IVL_setList(recvIVL, q, incounts[q], NULL) ;
}
IVfree(incounts) ;
/*
   ---------------------------------------------------
   load list myid of sendIVL into list myid of recvIVL
   ---------------------------------------------------
*/
IVL_listAndSize(sendIVL, myid, &sendcount, &sendvec) ;
IVL_setList(recvIVL, myid, sendcount, sendvec) ;
/*
   ---------------------------------------------------------
   now loop over the processes, send and receive information
   ---------------------------------------------------------
*/
for ( offset = 1, tag = firsttag ; offset < nproc ; offset++, tag++ ) {
   right = (myid + offset) % nproc ;
   if ( offset <= myid ) {
      left = myid - offset ;
   } else {
      left = nproc + myid - offset ;
   }
   IVL_listAndSize(sendIVL, right, &sendcount, &sendvec) ;
   IVL_listAndSize(recvIVL, left,  &recvcount, &recvvec) ;
   if ( sendcount > 0 ) {
      destination = right ;
      stats[0]++ ;
      stats[2] += sendcount*sizeof(int) ;
   } else {
      destination = MPI_PROC_NULL ;
   }
   if ( recvcount > 0 ) {
      source = left ;
      stats[1]++ ;
      stats[3] += recvcount*sizeof(int) ;
   } else {
      source = MPI_PROC_NULL ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
  "\n offset %d, recvcount %d, source %d, sendcount %d, destination %d",
              offset, recvcount, source, sendcount, destination) ;
      fflush(msgFile) ;
   }
/*
   -----------------
   do a send/receive
   -----------------
*/
   MPI_Sendrecv((void *) sendvec, sendcount, MPI_INT, destination, tag,
                (void *) recvvec, recvcount, MPI_INT, source,      tag,
                comm, &status) ;
   if ( source != MPI_PROC_NULL ) {
      MPI_Get_count(&status, MPI_INT, &count) ;
      if ( count != recvcount ) {
         fprintf(stderr, "\n fatal error in IVL_MPI_alltoall()"
                 "\n proc %d : source %d, count %d, recvcount %d\n",
                 myid, source, count, recvcount) ;
         exit(-1) ;
      }
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n send/recv completed") ;
      fflush(msgFile) ;
   }
}
return(recvIVL) ; }

/*--------------------------------------------------------------------*/
