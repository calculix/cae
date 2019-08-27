/*  IVallgather.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- 

   the IV objects objIV and ownersIV are found on each process.
   the ownersIV object is identical over all the processes, and
   owners[ii] tells which processes owns location ii of the obj[]
   vector. on return from this entry, the obj[] vector is replicated
   over all the processes. each process sends the (ii,obj[ii]) pairs
   that it owns to all the other processes.

   created -- 98apr02, cca
   -----------------------------------------------------------------
*/
void
IV_MPI_allgather (
   IV         *objIV,
   IV         *ownersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int          count, destination, ii, incount, iproc, jj, lasttag, left, 
             maxcount, myid, nowners, nproc, nvec, offset, 
             outcount, right, source, tag, tagbound, value ;
int          *counts, *inbuffer, *outbuffer, *owners, *vec ;
MPI_Status   status ;
/*
   ---------------
   check the input
   ---------------
*/
if ( objIV == NULL || ownersIV == NULL ) {
   fprintf(stderr, "\n fatal error in IV_MPI_allgather()"
           "\n objIV = %p, ownersIV = %p\n",
           objIV, ownersIV) ;
   exit(-1) ;
}
/*
   ----------------------------------------------
   get id of self, # of processes and # of fronts
   ----------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
IV_sizeAndEntries(objIV, &nvec, &vec) ;
IV_sizeAndEntries(ownersIV, &nowners, &owners) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside IV_MPI_allgather"
           "\n nproc = %d, myid = %d, nvec = %d, nowners = %d",
           nproc, myid, nvec, nowners) ;
   fflush(msgFile) ;
}
if ( nvec != nowners || vec == NULL || owners == NULL ) {
   fprintf(stderr, "\n fatal error in IV_MPI_allgather()"
           "\n nvec = %d, nowners = %d, vec = %p, owners = %p\n",
           nvec, nowners, vec, owners) ;
   exit(-1) ;
}
/*
   -------------------
   check the tag range
   -------------------
*/
lasttag = firsttag + nproc ;
tagbound = maxTagMPI(comm) ;
if ( firsttag < 0 || lasttag > tagbound ) {
   fprintf(stderr, "\n fatal error in IV_MPI_allgather()"
           "\n firsttag = %d, lasttag = %d, tagbound = %d\n",
           firsttag, lasttag, tagbound) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n objIV") ;
   IV_writeForHumanEye(objIV, msgFile) ;
   fprintf(msgFile, "\n\n ownersIV") ;
   IV_writeForHumanEye(ownersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------------------
   step 1 : determine the number of entries owned by each vector
   -------------------------------------------------------------
*/
counts = IVinit(nproc, 0) ;
for ( ii = 0 ; ii < nvec ; ii++ ) {
   if ( owners[ii] < 0 || owners[ii] >= nproc ) {
      fprintf(stderr, "\n owners[%d] = %d", ii, owners[ii]) ;
      exit(-1) ;
   }
   counts[owners[ii]]++ ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n counts") ;
   IVfprintf(msgFile, nproc, counts) ;
   fflush(msgFile) ;
}
/*
   -----------------------------
   set up the in and out buffers
   -----------------------------
*/
if ( counts[myid] > 0 ) {
   outbuffer = IVinit(2*counts[myid], -1) ;
   for ( ii = jj = 0 ; ii < nvec ; ii++ ) {
      if ( owners[ii] == myid ) {
         outbuffer[jj++] = ii ;
         outbuffer[jj++] = vec[ii] ;
      }
   }
   if ( jj != 2*counts[myid] ) {
      fprintf(msgFile, "\n jj = %d, 2*counts[%d] = %d",
              jj, myid, 2*counts[myid]) ;
      fprintf(stderr, "\n jj = %d, 2*counts[%d] = %d",
              jj, myid, 2*counts[myid]) ;
      exit(-1) ;
   }
} else {
   outbuffer = NULL ;
}
maxcount = IVmax(nproc, counts, &iproc) ;
if ( maxcount > 0 ) {
   inbuffer = IVinit(2*maxcount, -1) ;
} else {
   inbuffer = NULL ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n outbuffer %p, maxcount %d, inbuffer %p",
           outbuffer, maxcount, inbuffer) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------
   step 2: loop over the other processes
      send and receive information
   -------------------------------------
*/
outcount = 2*counts[myid] ;
for ( offset = 1, tag = firsttag ; offset < nproc ; offset++, tag++ ) {
   right = (myid + offset) % nproc ;
   if ( offset <= myid ) {
      left = myid - offset ;
   } else {
      left = nproc + myid - offset ;
   }
   if ( outcount > 0 ) {
      destination = right ;
      stats[0]++ ;
      stats[2] += outcount*sizeof(int) ;
   } else {
      destination = MPI_PROC_NULL ;
   }
   incount = 2*counts[left] ;
   if ( incount > 0 ) {
      source = left ;
      stats[1]++ ;
      stats[3] += incount*sizeof(int) ;
   } else {
      source = MPI_PROC_NULL ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n offset %d, source %d, destination %d",
              offset, source, destination) ;
      fflush(msgFile) ;
   }
/*
   -----------------
   do a send/receive
   -----------------
*/
   MPI_Sendrecv(outbuffer, outcount, MPI_INT, destination, tag,
                inbuffer,  incount,  MPI_INT, source,      tag,
                comm, &status) ;
   if ( source != MPI_PROC_NULL ) {
      MPI_Get_count(&status, MPI_INT, &count) ;
      if ( count != incount ) {
         fprintf(stderr,
                 "\n 1. fatal error in IV_MPI_allgather()"
                 "\n proc %d : source = %d, count = %d, incount = %d\n",
                 myid, source, count, incount) ;
         exit(-1) ;
      }
   }
/*
   ----------------------------
   set the values in the vector
   ----------------------------
*/
   for ( jj = 0 ; jj < incount ; jj += 2 ) {
      ii    = inbuffer[jj] ;
      value = inbuffer[jj+1] ;
      vec[ii] = value ;
   }
   if ( jj != incount ) {
      fprintf(msgFile, "\n jj = %d, incount = %d", jj, incount) ;
      fprintf(stderr, "\n jj = %d, incount = %d", jj, incount) ;
      exit(-1) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n after setting values") ;
      IVfprintf(msgFile, nvec, vec) ;
      fflush(msgFile) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(counts) ;
if ( outbuffer != NULL ) {
   IVfree(outbuffer) ;
}
if ( inbuffer != NULL ) {
   IVfree(inbuffer) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n leaving IV_MPI_gatherall()") ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
