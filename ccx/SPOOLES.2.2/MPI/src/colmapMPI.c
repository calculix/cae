/*  colmapMPI.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- after pivoting for a nonsymmetric factorization,
              some delayed columns may belong to a process other
              than its original owner. this method returns an
              IV object that maps columns to owning processes.

   created -- 98may22, cca
   -------------------------------------------------------------
*/
IV *
FrontMtx_MPI_colmapIV (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
int   buffersize, ii, iproc, J, myid, nDJ, neqns, nfront, nproc, 
      ncolJ, nToSend, v ;
int   *buffer, *counts, *frontOwners, *inbuffer, *outbuffer, 
      *colindJ, *colmap, *vtxToFront ;
IV    *colmapIV ;
/*
   -------------------------------------------
   get the process id and number of processors
   -------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
neqns      = frontmtx->neqns ;
vtxToFront = ETree_vtxToFront(frontmtx->frontETree) ;
IV_sizeAndEntries(frontOwnersIV, &nfront, &frontOwners) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside FrontMtx_MPI_colmapIV()"
           "\n myid = %d, nproc = %d, nfront = %d, neqns = %d",
           myid, nproc, nfront, neqns) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------------
   loop through the owned fronts and store each column in an 
   owned front that was originally owned by another processor
   ----------------------------------------------------------
*/
outbuffer = IVinit(neqns, -1) ;
for ( J = nToSend = 0 ; J < nfront ; J++ ) {
   if (  frontOwners[J] == myid 
      && (nDJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n front %d owned, nDJ = %d, ncolJ = %d",
                 J, nDJ, ncolJ) ;
         fflush(msgFile) ;
      }
      for ( ii = 0 ; ii < nDJ ; ii++ ) {
         v = colindJ[ii] ;
         if ( frontOwners[vtxToFront[v]] != myid ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, "\n column %d originally owned by %d",
                       v, frontOwners[vtxToFront[v]]) ;
               fflush(msgFile) ;
            }
            outbuffer[nToSend++] = v ;
         }
      }
   }
}
IVqsortUp(nToSend, outbuffer) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n shifted vertices") ;
   IVfprintf(msgFile, nToSend, outbuffer) ;
   fflush(msgFile) ;
}
counts = IVinit(nproc, 0) ;
/*
   --------------------------------------------
   use an all-gather call to get the number of
   moved columns that are owned by each process
   --------------------------------------------
*/
MPI_Allgather((void *) &nToSend, 1, MPI_INT, counts, 1, MPI_INT, comm) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n after the all-gather operation, counts") ;
   IVfprintf(msgFile, nproc, counts) ;
   fflush(msgFile) ;
}
buffersize = IVmax(nproc, counts, &iproc) ;
inbuffer   = IVinit(buffersize, -1) ;
/*
   -----------------------------------
   initialize the column map IV object
   -----------------------------------
*/
colmapIV = IV_new() ;
IV_init(colmapIV, neqns, NULL) ;
colmap = IV_entries(colmapIV) ;
IVgather(neqns, colmap, frontOwners, vtxToFront) ;
/*
   --------------------------------------------------------------
   loop over the other processes, receive vector of moved columns
   --------------------------------------------------------------
*/
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( counts[iproc] > 0 ) {
      if ( iproc == myid ) {
/*
        -------------------------------------
        send buffer vector to other processes
        -------------------------------------
*/
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n sending outbuffer to all processes") ;
            IVfprintf(msgFile, nToSend, outbuffer) ;
            fflush(msgFile) ;
         }
         MPI_Bcast(outbuffer, nToSend, MPI_INT, iproc, comm) ;
         buffer = outbuffer ;
      } else {
/*
        -----------------------------------------
        receive the vector from the other process
        -----------------------------------------
*/
         MPI_Bcast(inbuffer, counts[iproc], MPI_INT, iproc, comm) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n received inbuffer from process %d",
                    iproc) ;
            IVfprintf(msgFile, counts[iproc], inbuffer) ;
            fflush(msgFile) ;
         }
         buffer = inbuffer ;
      }
/*
      -------------------------
      set the column map values
      -------------------------
*/
      for ( ii = 0 ; ii < counts[iproc] ; ii++ ) {
         v = buffer[ii] ;
         colmap[v] = iproc ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(inbuffer)  ;
IVfree(outbuffer) ;
IVfree(counts)    ;

return(colmapIV) ; }

/*--------------------------------------------------------------------*/
