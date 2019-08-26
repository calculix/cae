/*  split.c  */

#include "../spoolesMPI.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to split a DenseMtx object by rows

   mtx         -- DenseMtx object
   rowmapIV    -- map from rows to owning processes
   firsttag    -- first tag to be used in these messages
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl      -- message level
   msgFile     -- message file
   comm        -- MPI communicator

   return value -- a new DenseMtx object filled with the owned rows 

   created  -- 98may16, cca
   modified -- 98sep26, cca
      mtx is not modified
   -----------------------------------------------------------------
*/
DenseMtx *
DenseMtx_MPI_splitByRows (
   DenseMtx   *mtx,
   IV         *rowmapIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
DenseMtx     *inmtx, *keepmtx, *outmtx ;
double       *inbuffer, *outbuffer ;
int          destination, ii, inbuffersize, incount, iproc, irow,
             lasttag, left, myid, ncol, ndouble, neqns, nkeep, 
             nmoved, nowned, nproc, nrecv, nrow, nsend, tagbound, 
             offset, outbuffersize, outcount, right, source, tag, type ;
int          *head, *link, *rowind, *rowmap, *rowsToRecv, *rowsToSend ;
MPI_Status   status ;
/*
   -------------------------------------------------
   get id of self, # of processes and # of equations
   -------------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
/*--------------------------------------------------------------------*/
{
int   rc = 1 ;
int   *rcs = IVinit(nproc, -1) ;
/*
   --------------
   check the data
   --------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_splitByRows()"
          "\n mtx is NULL\n") ;
   rc = -1 ;
}
if ( rowmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_splitByRows()"
          "\n rowmapIV is NULL\n") ;
   rc = -2 ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_splitByRows()"
          "\n msglvl > 0 and msgFile is NULL\n") ;
   rc = -3 ;
}
if ( firsttag < 0 ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_splitByRows()"
           "\n firsttag = %d\n", firsttag) ;
   rc = -4 ;
}
lasttag = firsttag + nproc ;
if ( lasttag > (tagbound = maxTagMPI(comm)) ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_splitByRows()"
           "\n lasttag = %d, tag_bound = %d", lasttag, tagbound) ;
   rc = -5 ;
}
MPI_Allgather((void *) &rc, 1, MPI_INT,
              (void *) rcs, 1, MPI_INT, comm) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( rcs[iproc] != 1 ) {
      if ( msgFile != NULL ) {
         fprintf(msgFile,
                 "\n fatal error in DenseMtx_MPI_splitByRows()"
                 "\n trouble with return code") ;
         IVfprintf(msgFile, nproc, rcs) ;
         MPI_Finalize() ;
         exit(rc) ;
      }
   }
}
IVfree(rcs) ;
}
/*--------------------------------------------------------------------*/
/*
   -----------------------
   get type and dimensions
   -----------------------
*/
type = mtx->type ;
IV_sizeAndEntries(rowmapIV, &neqns, &rowmap) ;
DenseMtx_dimensions(mtx, &nrow, &ncol) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside DenseMtx_MPI_splitByRows"
           "\n nproc = %d, myid = %d, neqns = %d, nrow = %d, ncol = %d",
           nproc, myid, neqns, nrow, ncol) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------------
   communicate the type's and ncol's to all the processors
   -------------------------------------------------------
*/
{
int   *ivec = IVinit(nproc, -1) ;
MPI_Allgather((void *) &type, 1, MPI_INT,
              (void *) ivec, 1, MPI_INT, comm) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( ivec[iproc] != type ) {
      if ( msgFile != NULL ) {
         fprintf(msgFile,
                 "\n fatal error in DenseMtx_MPI_splitByRows()"
                 "\n trouble with types") ;
         IVfprintf(msgFile, nproc, ivec) ;
         MPI_Finalize() ;
         exit(-1) ;
      }
   }
}
MPI_Allgather((void *) &ncol, 1, MPI_INT,
              (void *) ivec, 1, MPI_INT, comm) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( ivec[iproc] != ncol ) {
      if ( msgFile != NULL ) {
         fprintf(msgFile,
                 "\n fatal error in DenseMtx_MPI_splitByRows()"
                 "\n trouble with ncols") ;
         IVfprintf(msgFile, nproc, ivec) ;
         MPI_Finalize() ;
         exit(-1) ;
      }
   }
}
IVfree(ivec) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   get the counts of the entries to send to the other processors
   -------------------------------------------------------------
*/
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
rowsToSend = IVinit(2*nproc, 0) ;
rowsToRecv = rowsToSend + nproc ;
head = IVinit(nproc, -1) ;
link = IVinit(nrow, -1) ;
for ( ii = 0, nkeep = 0 ; ii < nrow ; ii++ ) {
   irow  = rowind[ii] ;
   iproc = rowmap[irow] ;
   link[ii] = head[iproc] ;
   head[iproc] = ii ;
   if ( iproc != myid ) {
      rowsToSend[iproc]++ ;
   } else {
      nkeep++ ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nkeep = %d, row send counts ", nkeep) ;
   IVfprintf(msgFile, nproc, rowsToSend) ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   do an all-to-all gather/scatter
   -------------------------------
*/
MPI_Alltoall((void *) rowsToSend, 1, MPI_INT,
             (void *) rowsToRecv, 1, MPI_INT, comm) ;
nowned = nkeep + IVsum(nproc, rowsToRecv) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nkeep = %d, row receive counts ", nkeep) ;
   IVfprintf(msgFile, nproc, rowsToRecv) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   determine the buffer size
   -------------------------
*/
nsend = IVmax(nproc, rowsToSend, &iproc) ;
nrecv = IVmax(nproc, rowsToRecv, &iproc) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nsend = %d, nrecv = %d", nsend, nrecv) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   allocate the send/receive DenseMtx objects
   -------------------------------------------
*/
if ( nsend > 0 ) {
   outmtx = DenseMtx_new() ;
   if ( mtx->inc1 == 1 ) {
      DenseMtx_init(outmtx, type, myid, -1, nsend, ncol, 1, nsend) ;
   } else {
      DenseMtx_init(outmtx, type, myid, -1, nsend, ncol, ncol, 1) ;
   }
} else {
   outmtx = NULL ;
}
if ( nrecv > 0 ) {
   inmtx = DenseMtx_new() ;
   if ( mtx->inc1 == 1 ) {
      DenseMtx_init(inmtx, type, myid, -1, nrecv, ncol, 1, nrecv) ;
   } else {
      DenseMtx_init(inmtx, type, myid, -1, nrecv, ncol, ncol, 1) ;
   }
} else {
   inmtx = NULL ;
}
/*
   -------------------------------------
   allocate the DenseMtx object to keep
   -------------------------------------
*/
keepmtx = DenseMtx_new() ;
if ( mtx->inc1 == 1 ) {
   DenseMtx_init(keepmtx, type, myid, -1, nowned, ncol, 1, nowned) ;
} else {
   DenseMtx_init(keepmtx, type, myid, -1, nowned, ncol, ncol, 1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n keepmtx object allocated") ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------------------
   copy all rows to keep from the input matrix into the keep matrix
   ----------------------------------------------------------------
*/
for ( ii = head[myid], nmoved = 0 ; ii != -1 ; ii = link[ii] ) {
   DenseMtx_copyRowAndIndex(keepmtx, nmoved, mtx, ii) ;
   nmoved++ ;
}
if ( nmoved > 0 ) {
/*
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n keepmtx") ;
      DenseMtx_writeForHumanEye(keepmtx, msgFile) ;
      fflush(msgFile) ;
   }
*/
}
/*
   --------------------------------------------------------------
   loop over the processes, gather their values and send them off
   --------------------------------------------------------------
*/
for ( offset = 1, tag = firsttag ; offset < nproc ; offset++, tag++ ) {
   right    = (myid + offset) % nproc ;
   left     = (nproc + myid - offset) % nproc ;
   outcount = rowsToSend[right] ;
   incount  = rowsToRecv[left] ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
       "\n\n ### process %d, send %d to right %d, recv %d from left %d",
       myid, outcount, right, incount, left) ;
      fflush(msgFile) ;
   }
   if ( outcount > 0 ) {
/*
      --------------------------
      load the out matrix object
      --------------------------
*/
      if ( mtx->inc1 == 1 ) {
         DenseMtx_init(outmtx, type, 
                       myid, -1, outcount, ncol, 1, outcount) ;
      } else {
         DenseMtx_init(outmtx, type, 
                       myid, -1, outcount, ncol, ncol, 1) ;
      }
      for ( ii = head[right], nmoved = 0 ; ii != -1 ; ii = link[ii] ) {
         DenseMtx_copyRowAndIndex(outmtx, nmoved, mtx, ii) ;
         nmoved++ ;
      }
      destination = right ;
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n\n outmtx for process %d", destination) ;
         DenseMtx_writeForHumanEye(outmtx, msgFile) ;
         fflush(msgFile) ;
      }
      stats[0]++ ;
      stats[2] += sizeof(double)*DV_size(&outmtx->wrkDV) ;
   } else {
/*
      ------------------------------------------
      set the destination to be the NULL process
      ------------------------------------------
*/
      destination = MPI_PROC_NULL ;
   }
   if ( incount > 0 ) {
/*
      ----------------------------------
      initialize the input matrix object
      ----------------------------------
*/
      if ( mtx->inc1 == 1 ) {
         DenseMtx_init(inmtx, type,
                       myid, -1, incount, ncol, 1, incount) ;
      } else {
         DenseMtx_init(inmtx, type,
                       myid, -1, incount, ncol, ncol, 1) ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n inmtx initialized to have %d rows",
                 incount) ;
         fflush(msgFile) ;
      }
      source = left ;
      stats[1]++ ;
      stats[3] += sizeof(double)*DV_size(&inmtx->wrkDV) ;
   } else {
      source = MPI_PROC_NULL ;
   }
/*
   -----------------
   do a send/receive
   -----------------
*/
   inbuffersize = outbuffersize = 0 ;
   inbuffer     = outbuffer     = NULL ;
   if ( outmtx != NULL ) {
      outbuffersize = DV_size(&outmtx->wrkDV) ;
      outbuffer     = DV_entries(&outmtx->wrkDV) ;
   }
   if ( inmtx != NULL ) {
      inbuffersize = DV_size(&inmtx->wrkDV) ;
      inbuffer     = DV_entries(&inmtx->wrkDV) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
              "\n inbuffersize  = %d, inbuffer  = %p"
              "\n outbuffersize = %d, outbuffer = %p",
              inbuffersize, inbuffer, outbuffersize, outbuffer) ;
      fflush(msgFile) ;
   }
   MPI_Sendrecv((void*) outbuffer, outbuffersize, MPI_DOUBLE, 
                destination, tag, (void*) inbuffer, inbuffersize, 
                MPI_DOUBLE, source, tag, comm, &status) ;
   if ( msglvl > 2 ) {
      MPI_Get_count(&status, MPI_DOUBLE, &ndouble) ;
      fprintf(msgFile, 
            "\n\n message received, source %d, tag %d, double count %d",
            status.MPI_SOURCE, status.MPI_TAG, ndouble) ;
      fprintf(msgFile, "\n MPI_ERROR = %d", status.MPI_ERROR) ;
      fflush(msgFile) ;
   }
   if ( source != MPI_PROC_NULL ) {
/*
      -------------------------------------
      initialize the object from its buffer
      -------------------------------------
*/
      DenseMtx_initFromBuffer(inmtx) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile,
                 "\n DenseMtx object intialized from its buffer") ;
         fflush(msgFile) ;
      }
      if ( msglvl > 4 ) {
         DenseMtx_writeForHumanEye(inmtx, msgFile) ;
         fflush(msgFile) ;
      }
   }
   if ( incount > 0 ) {
      if ( nkeep + incount > nowned ) {
         fprintf(msgFile, "\n fatal error in DenseMtx_splitByRows()"
              "\n nkeep = %d, nrecv = %d, nowned = %d",
              nkeep, nrecv, nowned) ;
         exit(-1) ;
      }
      for ( irow = 0 ; irow < incount ; irow++, nkeep++ ) {
         DenseMtx_copyRowAndIndex(keepmtx, nkeep, inmtx, irow) ;
      }
   }
}
/*
   -------------------------
   sort the rows and columns
   -------------------------
*/
DenseMtx_sort(keepmtx) ;
/*
   ------------------------------------------------------
   check that the matrix contains only the rows it should
   ------------------------------------------------------
*/
nrow   = keepmtx->nrow ;
rowind = keepmtx->rowind ;
for ( ii = 0 ; ii < nrow ; ii++ ) {
   irow = rowind[ii] ;
   if ( irow < 0 || irow >= neqns ) {
      fprintf(stderr, 
         "\n process %d : local row %d, global row %d, neqns = %d\n",
         myid, ii, irow, neqns) ;
      exit(-1) ;
   }
   if ( rowmap[irow] != myid ) {
      fprintf(stderr, 
           "\n process %d : local row %d, global row %d, map = %d\n",
           myid, ii, irow, rowmap[irow]) ;
      exit(-1) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( outmtx != NULL ) {
   DenseMtx_free(outmtx) ;
}
if ( inmtx != NULL ) {
   DenseMtx_free(inmtx) ;
}
IVfree(rowsToSend) ;
IVfree(head) ;
IVfree(link) ;

return(keepmtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to scatter a DenseMtx object by rows

   Xglobal     -- global DenseMtx object, significant only for root
   Xlocal      -- local DenseMtx object, if NULL on input, will
                  be created if necessary
   rowmapIV    -- map from rows to owning processes
   firsttag    -- first tag to be used in these messages
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl      -- message level
   msgFile     -- message file
   comm        -- MPI communicator

   return value -- Xlocal, a local DenseMtx object

   created  -- 98may16, cca
   modified -- 98sep26, cca
      mtx is not modified
   -----------------------------------------------------------------
*/
DenseMtx *
DenseMtx_MPI_splitFromGlobalByRows (
   DenseMtx   *Xglobal,
   DenseMtx   *Xlocal,
   IV         *rowmapIV,
   int        root,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
DenseMtx     *tempmtx ;
double       *buffer ;
int          ii, iproc, irow, maxnrow, myid, ncolX, nkeep, 
             nproc, nrowloc, nrowmap, nrowX, nsend, size, type ;
int          *counts, *head, *link, *rowind, *rowmap ;
MPI_Status   status ;
/*
   -------------------------------------------------
   get id of self, # of processes and # of equations
   -------------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( root < 0 || root >= nproc ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_splitByRows()"
           "\n root = %d, nproc = %d\n", root, nproc) ;
   MPI_Finalize() ;
   exit(-1) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------
   check the data
   --------------
*/
{
int   rc = 1 ;
int   *rcs = IVinit(nproc, -1) ;

if ( myid == root ) {
   if ( Xglobal == NULL ) {
      fprintf(stderr, 
             "\n fatal error in DenseMtx_MPI_splitFromGlobalByRows()"
             "\n Xglobal is NULL\n") ;
      rc = -1 ;
   }
   if ( rowmapIV == NULL ) {
      fprintf(stderr, 
             "\n fatal error in DenseMtx_MPI_splitFromGlobalByRows()"
             "\n rowmapIV is NULL\n") ;
      rc = -2 ;
   }
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(msgFile, 
          "\n fatal error in DenseMtx_MPI_splitFromGlobalByRows()"
          "\n msglvl > 0 and msgFile = NULL\n") ;
   rc = -3 ;
}
if ( firsttag < 0 ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_MPI_splitFromGlobalByRows()"
           "\n firsttag = %d\n", firsttag) ;
   rc = -4 ;
}
MPI_Allgather((void *) &rc, 1, MPI_INT,
              (void *) rcs, 1, MPI_INT, comm) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( rcs[iproc] != 1 ) {
      if ( msgFile != NULL ) {
         fprintf(msgFile,
                "\n fatal error in DenseMtx_MPI_splitFromGlobalByRows()"
                "\n trouble with return code") ;
         IVfprintf(msgFile, nproc, rcs) ;
         MPI_Finalize() ;
         exit(rc) ;
      }
   }
}
IVfree(rcs) ;
}
/*--------------------------------------------------------------------*/

if ( myid == root ) {
   type = Xglobal->type ;
   IV_sizeAndEntries(rowmapIV, &nrowmap, &rowmap) ;
   DenseMtx_dimensions(Xglobal, &nrowX, &ncolX) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n inside DenseMtx_MPI_splitFromGlobalByRows"
       "\n nproc = %d, myid = %d, nrowmap = %d, nrowX = %d, ncolX = %d",
       nproc, myid, nrowmap, nrowX, ncolX) ;
      fflush(msgFile) ;
   }
}
/*
   ----------------------------------------
   broadcast the type of entries and number
   of right hand sides to all processors
   ----------------------------------------
*/
MPI_Bcast((void *) &type, 1, MPI_INT, root, comm) ;
MPI_Bcast((void *) &ncolX, 1, MPI_INT, root, comm) ;
stats[0] += 2 ;
stats[1] += 2 ;
stats[2] += 2*sizeof(int) ;
stats[3] += 2*sizeof(int) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside DenseMtx_MPI_splitFromGlobalByRows"
      "\n type %d, ncolX %d", type, ncolX) ;
   fflush(msgFile) ;
}
if ( myid == root ) {
/*
   ------------------------------------------------
   create a head/link structure for the matrix rows
   ------------------------------------------------
*/
   DenseMtx_rowIndices(Xglobal, &nrowX, &rowind) ;
   counts = IVinit(nproc,  0) ;
   head   = IVinit(nproc, -1) ;
   link   = IVinit(nrowX, -1) ;
   for ( ii = nrowX - 1 ; ii >= 0 ; ii-- ) {
      irow = rowind[ii] ;
      iproc = rowmap[irow] ;
      link[ii] = head[iproc] ;
      head[iproc] = ii ;
      counts[iproc]++ ;
   }
} else {
   counts = NULL ;
}
/*
   -------------------------------------------------
   communicate the number of rows for each processor
   -------------------------------------------------
*/
MPI_Scatter((void *) counts,   1, MPI_INT, 
            (void *) &nrowloc, 1, MPI_INT, root, comm) ;
stats[0] += 1 ;
stats[1] += 1 ;
stats[2] += (nproc-1)*sizeof(int) ;
stats[3] += (nproc-1)*sizeof(int) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n nrowloc = %d", nrowloc) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------------------
   if necessary, create the local Xloc matrix, then initialize
   -----------------------------------------------------------
*/
if ( Xlocal == NULL ) {
   Xlocal = DenseMtx_new() ;
}
DenseMtx_init(Xlocal, type, -1, -1, nrowloc, ncolX, 1, nrowloc) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Xlocal after initialization") ;
   DenseMtx_writeForHumanEye(Xlocal, msgFile) ;
   fflush(msgFile) ;
}
if ( myid == root ) {
/*
   ---------------------------------
   load local matrix with owned rows
   ---------------------------------
*/
   if ( nrowloc > 0 ) {
      int   iglob, iloc = 0 ;
      for ( iglob = head[root] ; iglob != -1 ; iglob = link[iglob] ) {
         DenseMtx_copyRowAndIndex(Xlocal, iloc, Xglobal, iglob) ;
         iloc++ ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n Xlocal on root") ;
         DenseMtx_writeForHumanEye(Xlocal, msgFile) ;
         fflush(msgFile) ;
      }
   }
/*
   -----------------------------------------------------
   create a temporary matrix to send to other processors
   -----------------------------------------------------
*/
   counts[myid] = 0 ;
   maxnrow = IVmax(nproc, counts, &iproc) ;
   if ( maxnrow > 0 ) {
      DenseMtx   *tempmtx = DenseMtx_new() ;
      DenseMtx_init(tempmtx, type, -1, -1, maxnrow, ncolX, 1, maxnrow) ;
/*
      ----------------------------------
      loop over the processors
         collect owned rows into tempmtx
         send tempmtx to processor
      ----------------------------------
*/
      for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
         if ( iproc != root && (nrowloc = counts[iproc]) > 0 ) {
            DenseMtx_init(tempmtx, type, -1, -1,
                          nrowloc, ncolX, 1, nrowloc) ;
            nsend = 0 ;
            for ( ii = head[iproc] ; ii != -1 ; ii = link[ii] ) {
               DenseMtx_copyRowAndIndex(tempmtx, nsend, Xglobal, ii) ;
               nsend++ ;
            }
            if ( msglvl > 2 ) {
               fprintf(msgFile, "\n\n tempmtx for proc %d", iproc) ;
               DenseMtx_writeForHumanEye(tempmtx, msgFile) ;
               fflush(msgFile) ;
            }
            size   = DV_size(&tempmtx->wrkDV) ;
            buffer = DV_entries(&tempmtx->wrkDV) ;
            MPI_Send((void *) buffer, size, MPI_DOUBLE, 
                     iproc, firsttag, comm) ;
            stats[0] += 1 ;
            stats[2] += size*sizeof(double) ;
         }
      }
      DenseMtx_free(tempmtx) ;
   }
/*
   ------------------------
   free the working storage
   ------------------------
*/
   IVfree(head) ;
   IVfree(link) ;
   IVfree(counts) ;
} else {
/*
   ------------------
   non-root processor
   ------------------
*/
   if ( nrowloc > 0 ) {
      size   = DV_size(&Xlocal->wrkDV) ;
      buffer = DV_entries(&Xlocal->wrkDV) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n size = %d, buffer = %p", size, buffer) ;
         fflush(msgFile) ;
      }
      MPI_Recv((void *) buffer, size, MPI_DOUBLE, 
               root, firsttag, comm, &status);
      stats[1] += 1 ;
      stats[3] += size*sizeof(double) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n Xlocal rec'd from root %d", root) ;
         fflush(msgFile) ;
      }
      DenseMtx_initFromBuffer(Xlocal) ;
      if ( msglvl > 2 ) {
         DenseMtx_writeForHumanEye(Xlocal, msgFile) ;
         fflush(msgFile) ;
      }
   } else {
      Xlocal = NULL ;
   }
}
if ( msglvl > 3 ) {
   if ( Xlocal != NULL ) {
      fprintf(msgFile, "\n\n Xlocal") ;
      DenseMtx_writeForHumanEye(Xlocal, msgFile) ;
   } else {
      fprintf(msgFile, "\n\n Xlocal is NULL") ;
   }
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n leaving DenseMtx_splitFromGlobalByRows()") ;
   fflush(msgFile) ;
}
return(Xlocal) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to merge a DenseMtx object by rows

   Xlocal      -- DenseMtx object, can be NULL
   Xglobal     -- DenseMtx object, can be NULL
                  significant only for root
   firsttag    -- first tag to be used in these messages
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl      -- message level
   msgFile     -- message file
   comm        -- MPI communicator

   return value -- 
      if processor is root
          Xglobal is returned, if was NULL on input, it is created
      else
          NULL
      endif

   Xlocal is not modified

   created  -- 98sep27, cca
   -----------------------------------------------------------------
*/
DenseMtx *
DenseMtx_MPI_mergeToGlobalByRows (
   DenseMtx   *Xglobal,
   DenseMtx   *Xlocal,
   int        root,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
double   *buffer ;
int      iproc, maxnrow, myid, ncolX, nproc, nrowXloc, size, type ;
int      *incounts ;
/*
   -------------------------------------------------
   get id of self, # of processes and # of equations
   -------------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( root < 0 || root >= nproc ) {
   fprintf(stderr, "\n fatal error in DenseMtx_MPI_splitByRows()"
           "\n root = %d, nproc = %d\n", root, nproc) ;
   MPI_Finalize() ;
   exit(-1) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------
   check the data
   --------------
*/
{
int   rc = 1 ;
int   *ivec = IVinit(nproc, -1) ;

if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(msgFile, 
          "\n fatal error in DenseMtx_MPI_mergeToGlobalByRows()"
          "\n msglvl > 0 and msgFile = NULL\n") ;
   rc = -2 ;
}
if ( firsttag < 0 ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_MPI_mergeToGlobalByRows()"
           "\n firsttag = %d\n", firsttag) ;
   rc = -3 ;
}
MPI_Allgather((void *) &rc, 1, MPI_INT,
              (void *) ivec, 1, MPI_INT, comm) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( ivec[iproc] != 1 ) {
      if ( msgFile != NULL ) {
         fprintf(msgFile,
                "\n fatal error in DenseMtx_MPI_mergeToGlobalByRows()"
                "\n trouble with return code") ;
         IVfprintf(msgFile, nproc, ivec) ;
         MPI_Finalize() ;
         exit(rc) ;
      }
   }
}
/*
   ------------------------------
   check the type of the matrices
   ------------------------------
*/
type = (Xlocal != NULL) ? Xlocal->type : -1 ;
MPI_Allgather((void *) &type, 1, MPI_INT,
              (void *) ivec, 1, MPI_INT, comm) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( ivec[iproc] != -1 ) {
      if ( type == -1 ) {
         type = ivec[iproc] ;
      } else if ( type != ivec[iproc] ) {
         if ( msgFile != NULL ) {
            fprintf(msgFile,
                  "\n fatal error in DenseMtx_MPI_mergeToGlobalByRows()"
                  "\n trouble with types\n") ;
            IVfprintf(msgFile, nproc, ivec) ;
            MPI_Finalize() ;
            exit(-1) ;
         }
      }
   }
}
/*
   --------------------------------------
   check the # of columns of the matrices
   --------------------------------------
*/
ncolX = (Xlocal != NULL) ? Xlocal->ncol : 0 ;
MPI_Allgather((void *) &ncolX, 1, MPI_INT,
              (void *) ivec, 1, MPI_INT, comm) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( ivec[iproc] != 0 ) {
      if ( ncolX == 0 ) {
         ncolX = ivec[iproc] ;
      } else if ( ncolX != ivec[iproc] ) {
         if ( msgFile != NULL ) {
            fprintf(msgFile,
                  "\n fatal error in DenseMtx_MPI_mergeToGlobalByRows()"
                  "\n trouble with ncolX\n") ;
            IVfprintf(msgFile, nproc, ivec) ;
            MPI_Finalize() ;
            exit(-1) ;
         }
      }
   }
}
IVfree(ivec) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   gather the number of incoming rows onto processor root
   ------------------------------------------------------
*/
nrowXloc = (Xlocal != NULL) ? Xlocal->nrow : 0 ;
if ( myid == root ) {
   incounts = IVinit(nproc, 0) ;
} else {
   incounts = NULL ;
}
MPI_Gather(&nrowXloc, 1, MPI_INT, 
           (void *) incounts, 1, MPI_INT, root, comm) ;
if ( myid == root ) {
   DenseMtx     *tempmtx ;
   int          count, iglob, iloc, incount, nrowXglobal ;
   MPI_Status   status ;
/*
   ----------------------------------------------------------
   if necessary, create the global matrix, then initialize it
   ----------------------------------------------------------
*/
   nrowXglobal = IVsum(nproc, incounts) ;
   if ( Xglobal == NULL ) {
      Xglobal = DenseMtx_new() ;
   }
   DenseMtx_init(Xglobal, type, -1, -1, 
                 nrowXglobal, ncolX, 1, nrowXglobal) ;
/*
   ---------------------------------
   load local matrix with owned rows
   ---------------------------------
*/
   for ( iloc = iglob = 0 ; iloc < nrowXloc ; iloc++, iglob++ ) {
      DenseMtx_copyRowAndIndex(Xglobal, iglob, Xlocal, iloc) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n after loading Xlocal on root") ;
      DenseMtx_writeForHumanEye(Xglobal, msgFile) ;
      fflush(msgFile) ;
   }
/*
   ----------------------------------------------------------
   create a temporary matrix to receive from other processors
   ----------------------------------------------------------
*/
   incounts[root] = 0 ;
   maxnrow = IVmax(nproc, incounts, &iproc) ;
   tempmtx = DenseMtx_new() ;
   DenseMtx_init(tempmtx, type, -1, -1, maxnrow, ncolX, 1, maxnrow) ;
   size   = DV_size(&tempmtx->wrkDV) ;
   buffer = DV_entries(&tempmtx->wrkDV) ;
/*
   ----------------------------
   loop over the processors
      receive rows into tempmtx
   ----------------------------
*/
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      if ( iproc != root && (incount = incounts[iproc]) > 0 ) {
         MPI_Recv((void *) buffer, size, MPI_DOUBLE, 
                  iproc, firsttag, comm, &status) ;
         MPI_Get_count(&status, MPI_DOUBLE, &count) ;
         stats[1] += 1 ;
         stats[3] += count*sizeof(double) ;
         DenseMtx_initFromBuffer(tempmtx) ;
         for ( iloc = 0 ; iloc < incount ; iloc++, iglob++ ) {
            DenseMtx_copyRowAndIndex(Xglobal, iglob, tempmtx, iloc) ;
         }
      }
   }
/*
   ------------------------
   free the working storage
   ------------------------
*/
   IVfree(incounts) ;
   DenseMtx_free(tempmtx) ;
} else {
/*
   ------------------
   non-root processor
   ------------------
*/
   if ( nrowXloc > 0 ) {
      size   = DV_size(&Xlocal->wrkDV) ;
      buffer = DV_entries(&Xlocal->wrkDV) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n size = %d, buffer = %p", size, buffer) ;
         fflush(msgFile) ;
      }
      MPI_Send((void *) buffer, size, MPI_DOUBLE, 
               root, firsttag, comm) ;
      stats[0] += 1 ;
      stats[2] += size*sizeof(double) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n Xlocal sent to root %d", root) ;
         fflush(msgFile) ;
      }
   }
   Xglobal = NULL ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n leaving DenseMtx_mergeToGlobalByRows()") ;
   fflush(msgFile) ;
}
return(Xglobal) ; }

/*--------------------------------------------------------------------*/
