/*  IVLallgather.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- 

   the IVL object ivl and IV object ownersIV are both found on 
   each process.  the ownersIV object is identical over all the 
   processes, and owners[ii] tells which processes owns list ii 
   of the ivl object. on return from this method, the ivl object 
   is replicated over all the processes. each process sends 
   the lists that it owns to all the other processes.

   created -- 98apr03, cca
   -------------------------------------------------------------
*/
void
IVL_MPI_allgather (
   IVL        *ivl,
   IV         *ownersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int          count, destination, ii, ilist, incount, jlist, 
             jproc, left, maxcount, myid, nlist, nmylists, 
             notherlists, nowners, nproc, offset, outcount, 
             right, size, source, tag ;
int          *counts, *inbuffer, *list, *outbuffer, *owners ;
MPI_Status   status ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || ownersIV == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_MPI_allgather()"
           "\n ivl = %p, ownersIV = %p\n",
           ivl, ownersIV) ;
   exit(-1) ;
}
/*
   ----------------------------------------------
   get id of self, # of processes and # of fronts
   ----------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
nlist = ivl->nlist ;
IV_sizeAndEntries(ownersIV, &nowners, &owners) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside IVL_MPI_allgather()"
           "\n nproc = %d, myid = %d, nlist = %d, nowners = %d",
           nproc, myid, nlist, nowners) ;
   fflush(msgFile) ;
}
if ( nlist != nowners || owners == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_MPI_allgather()"
           "\n nlist = %d, nowners = %d, owners = %p\n",
           nlist, nowners, owners) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n ivl") ;
   IVL_writeForHumanEye(ivl, msgFile) ;
   fprintf(msgFile, "\n\n ownersIV") ;
   IV_writeForHumanEye(ownersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------
   step 1 : determine the size of the message that
            this process will send to the others
   -----------------------------------------------
*/
for ( ilist = 0, outcount = 1 ; ilist < nlist ; ilist++ ) {
   if ( owners[ilist] < 0 || owners[ilist] >= nproc ) {
      fprintf(stderr, "\n owners[%d] = %d", ilist, owners[ilist]) ;
      exit(-1) ;
   }
   if ( owners[ilist] == myid ) {
      outcount += 2 ;
      IVL_listAndSize(ivl, ilist, &size, &list) ;
      outcount += size ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n outcount = %d", outcount) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------
   do an all-to-all gather/scatter
   counts[jproc] = # of int's in the message from jproc
   ----------------------------------------------------
*/
counts = IVinit(nproc, 0) ;
counts[myid] = outcount ;
MPI_Allgather((void *) &counts[myid], 1, MPI_INT,
              (void *) counts,  1, MPI_INT, comm) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n counts") ;
   IVfprintf(msgFile, nproc, counts) ;
   fflush(msgFile) ;
}
/*
   -----------------------------
   set up the in and out buffers
   -----------------------------
*/
if ( outcount > 0 ) {
   outbuffer = IVinit(outcount, -1) ;
   for ( ilist = nmylists = 0, ii = 1 ; ilist < nlist ; ilist++ ) {
      if ( owners[ilist] == myid ) {
         nmylists++ ;
         IVL_listAndSize(ivl, ilist, &size, &list) ;
         outbuffer[ii++] = ilist ;
         outbuffer[ii++] = size  ;
         if ( size > 0 ) {
            IVcopy(size, &outbuffer[ii], list) ;
            ii += size ;
         }
      }
   }
   outbuffer[0] = nmylists ;
   if ( ii != outcount ) {
      fprintf(stderr, "\n myid = %d, ii = %d, outcount = %d",
              myid, ii, outcount) ;
      fprintf(msgFile, "\n myid = %d, ii = %d, outcount = %d",
              myid, ii, outcount) ;
      exit(-1) ;
   }
} else {
   outbuffer = NULL ;
}
maxcount = IVmax(nproc, counts, &jproc) ;
if ( maxcount > 0 ) {
   inbuffer = IVinit(maxcount, -1) ;
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
   incount = counts[left] ;
   if ( incount > 0 ) {
      source = left ;
      stats[1]++ ;
      stats[3] += incount*sizeof(int) ;
   } else {
      source = MPI_PROC_NULL ;
   }
   if ( msglvl > 2 ) {
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
                 "\n 1. fatal error in IVL_MPI_allgather()"
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
   notherlists = inbuffer[0] ;
   for ( ilist = 0, ii = 1 ; ilist < notherlists ; ilist++ ) {
      jlist = inbuffer[ii++] ;
      size  = inbuffer[ii++] ;
      if ( size > 0 ) {
         IVL_setList(ivl, jlist, size, &inbuffer[ii]) ;
         ii += size ;
      }
   }
   if ( ii != incount ) {
      fprintf(msgFile, "\n ii = %d, incount = %d", ii, incount) ;
      fprintf(stderr, "\n ii = %d, incount = %d", ii, incount) ;
      exit(-1) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n after setting values") ;
      IVL_writeForHumanEye(ivl, msgFile) ;
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
   fprintf(msgFile, "\n\n leaving IVL_MPI_gatherall()") ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
