/*  split.c  */

#include "../spoolesMPI.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
static int numberOfDoubles ( int count, int type ) ;
static void writeBuffers ( int count, int ibuf1[], int ibuf2[],
   double dbuf[], int type, FILE *fp ) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to split a distributed InpMtx object

   inpmtx     -- pointer to the local InpMtx object
   mapIV      -- pointer to the map from vertices to processes
   firsttag   -- first tag value, one will be used
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl     -- local message level
   msgFile    -- local message file
   comm       -- MPI communication structure

   return value -- pointer to the new InpMtx object 
      that contains the owned entries.

   created  -- 98may16, cca
   modified -- 98sep26, cca
      inpmtx is not modified
   ------------------------------------------------------------
*/
InpMtx *
InpMtx_MPI_split (
   InpMtx     *inpmtx,
   IV         *mapIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
double       *dvec, *inbuffer, *outbuffer ;
InpMtx       *keepmtx ;
int          destination, ii, ient, incount, iproc, left, maxincount, 
             maxoutcount, myid, nDblIn, nDblOut, nent, nkeep, nowned, 
             nproc, nvector, nvtx, offset, outcount, right, source, 
             tagbound, type, val, vecid ;
int          *head, *incounts, *ivec1, *ivec2, *link,
             *map, *outcounts ;
MPI_Status   status ;
/*
   --------------------------------------
   get id of self and number of processes
   --------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
/*
   --------------------------
   check the data and the map
   --------------------------
*/
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_MPI_split()"
           "\n inpmtx is NULL") ;
   exit(-1) ;
}
if ( inpmtx->storageMode != INPMTX_BY_VECTORS ) {
   fprintf(stderr, "\n fatal error in InpMtx_MPI_split()"
           "\n inpmtx must be stored by vectors") ;
   exit(-1) ;
}
if ( firsttag < 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_MPI_split()"
           "\n firsttag = %d\n", firsttag) ;
   exit(-1) ;
}
if ( firsttag > (tagbound = maxTagMPI(comm)) ) {
   fprintf(stderr, "\n fatal error in InpMtx_MPI_split()"
           "\n firsttag = %d, tagbound = %d\n", firsttag, tagbound) ;
   exit(-1) ;
}
type  = InpMtx_inputMode(inpmtx) ;
switch ( type ) {
case INPMTX_INDICES_ONLY :
   dvec = NULL ;
   break ;
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   dvec = InpMtx_dvec(inpmtx) ;
   break ;
default :
   fprintf(stderr, "\n fatal error in InpMtx_MPI_split()"
           "\n bad inputMode\n") ;
   exit(-1) ;
   break ;
}
nent    = InpMtx_nent(inpmtx)    ;
ivec1   = InpMtx_ivec1(inpmtx)   ;
ivec2   = InpMtx_ivec2(inpmtx)   ;
nvector = InpMtx_nvector(inpmtx) ;
IV_sizeAndEntries(mapIV, &nvtx, &map) ;
if ( nvtx <= 0 || map == NULL ) {
   fprintf(stderr, "\n process %d : nvtx = %d, map = %p",
           myid, nvtx, map) ;
   exit(-1) ;
}
if ( (val = IVmin(nent, ivec1, &ient)) < 0 ) {
   fprintf(stderr, "\n process %d : IV_min(ivec1) = %d", 
           myid, val) ;
   exit(-1) ;
}
if ( (val = IVmax(nent, ivec1, &ient)) >= nvtx ) {
   fprintf(stderr, "\n process %d : IV_max(ivec1) = %d", 
           myid, val) ;
   exit(-1) ;
}
if ( (val = IV_min(mapIV)) < 0 ) {
   fprintf(stderr, "\n process %d : IVmin(mapIV) = %d", 
           myid, val) ;
   exit(-1) ;
}
if ( (val = IV_max(mapIV)) >= nproc ) {
   fprintf(stderr, "\n process %d : IVmax(mapIV) = %d", 
           myid, val) ;
   exit(-1) ;
}
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   compute the number of entries that 
   must be sent to each other process
   ----------------------------------
*/
outcounts = IVinit(nproc, 0) ;
incounts  = IVinit(nproc, 0) ;
if ( nvector > 0 ) {
   int   *vecids = InpMtx_vecids(inpmtx)  ;
   int   *sizes  = InpMtx_sizes(inpmtx)   ;

   head = IVinit(nproc, -1) ;
   link = IVinit(nvector, -1) ;
   for ( ii = nvector - 1 ; ii >= 0 ; ii-- ) {
      vecid = vecids[ii] ;
      iproc = map[vecid] ;
      link[ii] = head[iproc] ;
      head[iproc] = ii ;
      outcounts[iproc] += sizes[ii] ;
   }
   nkeep = outcounts[myid] ;
   outcounts[myid] = 0 ;
} else {
   head = NULL ;
   link = NULL ;
   nkeep = 0 ;
}
/*
   -------------------------------
   do an all-to-all gather/scatter
   -------------------------------
*/
MPI_Alltoall((void *) outcounts, 1, MPI_INT,
             (void *) incounts,  1, MPI_INT, comm) ;
incounts[myid] = 0 ;
nowned = nkeep + IVsum(nproc, incounts) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n nkeep %d, nowned %d", nkeep, nowned) ;
   fprintf(msgFile, "\n\n incounts") ;
   IVfprintf(msgFile, nproc, incounts) ;
   fprintf(msgFile, "\n\n outcounts") ;
   IVfprintf(msgFile, nproc, outcounts) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------
   allocate the send/receive buffers
   ---------------------------------
*/
if ( (maxincount = IVmax(nproc, incounts, &iproc)) > 0 ) {
   nDblIn = numberOfDoubles(maxincount, type) ;
   inbuffer = DVinit2(nDblIn) ;
} else {
   nDblIn   =   0  ;
   inbuffer = NULL ;
}
if ( (maxoutcount = IVmax(nproc, outcounts, &iproc)) > 0 ) {
   nDblOut = numberOfDoubles(maxoutcount, type) ;
   outbuffer = DVinit2(nDblOut) ;
} else {
   nDblOut   =   0  ;
   outbuffer = NULL ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n\n nDblIn %d, inbuffer %p, nDblOut %d, outbuffer %p",
           nDblIn, inbuffer, nDblOut, outbuffer) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   set up the new InpMtx object to hold the owned entries
   -------------------------------------------------------
*/
keepmtx = InpMtx_new() ;
InpMtx_init(keepmtx, inpmtx->coordType, type, nowned, 0) ;
/*
   ---------------------------------------------
   extract the triples from the original matrix,
   and load the triples into the keep matrix
   ---------------------------------------------
*/
if ( nkeep > 0 ) {
   int   *offsets = InpMtx_offsets(inpmtx) ;
   int   *sizes   = InpMtx_sizes(inpmtx)   ;
   int   ii, ient, size ;

   for ( ii = head[myid] ; ii != -1 ; ii = link[ii] ) {
      ient = offsets[ii] ;
      size = sizes[ii]   ;
      switch ( type ) {
      case INPMTX_INDICES_ONLY :
         InpMtx_inputTriples(keepmtx, size, ivec1 + ient, ivec2 + ient);
         break ;
      case SPOOLES_REAL :
         InpMtx_inputRealTriples(keepmtx, size, ivec1 + ient, 
                                 ivec2 + ient, dvec + ient) ;
         break ;
      case SPOOLES_COMPLEX :
         InpMtx_inputComplexTriples(keepmtx, size, ivec1 + ient, 
                                    ivec2 + ient, dvec + 2*ient) ;
         break ;
      }
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, 
              "\n\n keepmtx after storing owned original entries") ;
      InpMtx_writeForHumanEye(keepmtx, msgFile) ;
      fflush(msgFile) ;
   }
}
/*
   ----------------------------------
   loop over the other processes, 
      gather values and send them off
      receive values
   ----------------------------------
*/
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n\n ready to split entries and send to other processes") ;
   fflush(msgFile) ;
}
for ( offset = 1 ; offset < nproc ; offset++ ) {
   right = (myid + offset) % nproc ;
   left  = (offset <= myid) ? (myid - offset) : (nproc + myid - offset);
   outcount = outcounts[right] ;
   incount  = incounts[left] ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
         "\n ### process %d, send %d to right %d, recv %d from left %d",
              myid, outcount, right, incount, left) ;
      fflush(msgFile) ;
   }
   if ( outcount > 0 ) {
      double   *dbuf ;
      int      *offsets = InpMtx_offsets(inpmtx) ;
      int      *sizes   = InpMtx_sizes(inpmtx)   ;
      int      *ibuf1, *ibuf2 ;
      int      ii, ient, jent, size ;
/*
      -----------------------------------------
      load the entries to send to process right
      -----------------------------------------
*/
      ibuf1 = (int *) outbuffer ;
      ibuf2 = ibuf1 + outcount ;
      dbuf  = (double *) (ibuf2 + outcount) ;
      for ( ii = head[right], jent = 0 ; ii != -1 ; ii = link[ii] ) {
         ient = offsets[ii] ;
         size = sizes[ii]   ;
         IVcopy(size, ibuf1 + jent, ivec1 + ient) ;
         IVcopy(size, ibuf2 + jent, ivec2 + ient) ;
         switch ( type ) {
         case SPOOLES_REAL :
            DVcopy(size, dbuf + jent, dvec + ient) ;
            break ;
         case SPOOLES_COMPLEX :
            DVcopy(2*size, dbuf + 2*jent, dvec + 2*ient) ;
            break ;
         }
         jent += size ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n outcount %d, ibuf1 %p, ibuf2 %p, dbuf %p",
                 outcount, ibuf1, ibuf2, dbuf) ;
         writeBuffers(outcount, ibuf1, ibuf2, dbuf, type, msgFile) ;
         fflush(msgFile) ;
      }
      destination = right ;
      nDblOut = numberOfDoubles(outcount, type) ;
      stats[0]++ ; stats[2] += nDblOut*sizeof(double) ;
   } else {
      destination = MPI_PROC_NULL ;
      nDblOut  = 0 ;
   }
   if ( incount > 0 ) {
      source = left ;
   } else {
      source = MPI_PROC_NULL ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n nDblOut = %d, nDblIn = %d",
              nDblOut, nDblIn) ;
      fflush(msgFile) ;
   }
/*
   -----------------
   do a send receive
   -----------------
*/
   MPI_Sendrecv((void *) outbuffer, nDblOut, MPI_DOUBLE, 
                destination, firsttag, (void *) inbuffer, nDblIn, 
                MPI_DOUBLE, source, firsttag, comm, &status) ;
   if ( incount > 0 ) {
/*
      -------------------------------------
      load the triples into the keep matrix
      -------------------------------------
*/
      int      count ;
      int      *ibuf1 = (int *) inbuffer ;
      int      *ibuf2 = ibuf1 + incount ;
      double   *dbuf  = (double *) (ibuf2 + incount) ;

      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n incoming buffer") ;
         writeBuffers(incount, ibuf1, ibuf2, dbuf, type, msgFile) ;
         fflush(msgFile) ;
      }
      switch ( type ) {
      case INPMTX_INDICES_ONLY :
         InpMtx_inputTriples(keepmtx, incount, ibuf1, ibuf2) ;
         break ;
      case SPOOLES_REAL :
         InpMtx_inputRealTriples(keepmtx, incount, ibuf1, ibuf2, dbuf) ;
         break ;
      case SPOOLES_COMPLEX :
         InpMtx_inputComplexTriples(keepmtx, incount, 
                                    ibuf1, ibuf2, dbuf) ;
         break ;
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, 
                 "\n\n keepmtx after storing received entries") ;
         InpMtx_writeForHumanEye(keepmtx, msgFile) ;
         fflush(msgFile) ;
      }
      MPI_Get_count(&status, MPI_DOUBLE, &count) ;
      stats[1]++ ; stats[3] += count*sizeof(double) ;
   }
}
if ( keepmtx != NULL ) {
/*
   ----------------------------------------------------
   sort and compress the entries and convert to vectors
   ----------------------------------------------------
*/
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n before changing storage mode to %d", 2) ;
      InpMtx_writeForHumanEye(keepmtx, msgFile) ;
      fflush(msgFile) ;
   }
   InpMtx_changeStorageMode(keepmtx, INPMTX_BY_VECTORS) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n after changing storage mode to %d", 2) ;
      InpMtx_writeForHumanEye(keepmtx, msgFile) ;
      fflush(msgFile) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
if ( inbuffer != NULL ) {
   DVfree(inbuffer) ;
}
if ( outbuffer != NULL ) {
   DVfree(outbuffer) ;
}
IVfree(outcounts) ;
IVfree(incounts) ;
if ( link != NULL ) {
   IVfree(head) ;
   IVfree(link) ;
}
return(keepmtx) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to split a distributed InpMtx object

   Aglobal  -- pointer to the global InpMtx object,
               significant only for root processor
   Alocal   -- pointer to the local InpMtx object,
               if NULL and will be nonempty, it will be created
   mapIV    -- pointer to the map from vertices to processes
   root     -- processor that presently owns the global matrix
   firsttag -- first tag value, one will be used
   stats[4] -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl  -- local message level
   msgFile -- local message file
   comm    -- MPI communication structure

   return value -- pointer to the local InpMtx object 
      that contains the owned entries.

   created  -- 98may16, cca
   modified -- 98sep26, cca
      inpmtx is not modified
   ------------------------------------------------------------
*/
InpMtx *
InpMtx_MPI_splitFromGlobal (
   InpMtx     *Aglobal,
   InpMtx     *Alocal,
   IV         *mapIV,
   int        root,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
double   *dvec, *inbuffer, *outbuffer ;
int      coordType, ii, ient, incount, iproc, maxoutcount, myid, nDblIn,
         nDblOut, nent, nowned, nproc, nvector, nvtx, 
         outcount, tagbound, type, val, vecid ;
int      *head, *ivec1, *ivec2, *link, *map, *outcounts ;
/*
   --------------------------------------
   get id of self and number of processes
   --------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
/*
   --------------------------
   check the data and the map
   --------------------------
*/
if ( myid == root ) {
   int   rc = 1 ;
   if ( Aglobal == NULL ) {
      fprintf(stderr, "\n fatal error in InpMtx_MPI_splitFromGlobal()"
              "\n Aglobal is NULL") ;
      rc = -1 ;
   }
   if ( Aglobal->storageMode != INPMTX_BY_VECTORS ) {
      fprintf(stderr, "\n fatal error in InpMtx_MPI_splitFromGlobal()"
              "\n Aglobal must be stored by vectors") ;
      rc = -2 ;
   }
   if ( firsttag < 0 ) {
      fprintf(stderr, "\n fatal error in InpMtx_MPI_splitFromGlobal()"
              "\n firsttag = %d\n", firsttag) ;
      rc = -3 ;
   }
   if ( firsttag > (tagbound = maxTagMPI(comm)) ) {
      fprintf(stderr, "\n fatal error in InpMtx_MPI_splitFromGlobal()"
              "\n firsttag = %d, tagbound = %d\n", firsttag, tagbound) ;
      rc = -4 ;
   }
   type  = InpMtx_inputMode(Aglobal) ;
   switch ( type ) {
   case INPMTX_INDICES_ONLY :
      dvec = NULL ;
      break ;
   case SPOOLES_REAL :
   case SPOOLES_COMPLEX :
      dvec = InpMtx_dvec(Aglobal) ;
      break ;
   default :
      fprintf(stderr, "\n fatal error in InpMtx_MPI_splitFromGlobal()"
              "\n bad inputMode\n") ;
      rc = -5 ;
      break ;
   }
   nent    = InpMtx_nent(Aglobal)    ;
   ivec1   = InpMtx_ivec1(Aglobal)   ;
   ivec2   = InpMtx_ivec2(Aglobal)   ;
   nvector = InpMtx_nvector(Aglobal) ;
   IV_sizeAndEntries(mapIV, &nvtx, &map) ;
   if ( nvtx <= 0 || map == NULL ) {
      fprintf(stderr, "\n process %d : nvtx = %d, map = %p",
              myid, nvtx, map) ;
      rc = -6 ;
   }
   if ( (val = IVmin(nent, ivec1, &ient)) < 0 ) {
      fprintf(stderr, "\n process %d : IV_min(ivec1) = %d", 
              myid, val) ;
      rc = -7 ;
   }
   if ( (val = IVmax(nent, ivec1, &ient)) >= nvtx ) {
      fprintf(stderr, "\n process %d : IV_max(ivec1) = %d", 
              myid, val) ;
      rc = -8 ;
   }
   if ( (val = IV_min(mapIV)) < 0 ) {
      fprintf(stderr, "\n process %d : IVmin(mapIV) = %d", 
              myid, val) ;
      rc = -9 ;
   }
   if ( (val = IV_max(mapIV)) >= nproc ) {
      fprintf(stderr, "\n process %d : IVmax(mapIV) = %d", 
              myid, val) ;
      rc = -10 ;
   }
/*
   --------------------------------------------------
   if the error code is not 1, broadcast it and exit,
   otherwise broadcast the entries type
   --------------------------------------------------
*/
   if ( rc != 1 ) {
      MPI_Bcast((void *) &rc, 1, MPI_INT, root, comm) ;
      MPI_Finalize() ;
      exit(rc) ;
   } else {
      MPI_Bcast((void *) &type,      1, MPI_INT, root, comm) ;
      coordType = Aglobal->coordType ;
      MPI_Bcast((void *) &coordType, 1, MPI_INT, root, comm) ;
   }
} else {
/*
   --------------------------------------
   receive the type of entries, 
   if < 0 then an error has been found
   by the root process, finalize and exit
   --------------------------------------
*/
   MPI_Bcast((void *) &type, 1, MPI_INT, root, comm) ;
   if ( type < 0 ) {
      MPI_Finalize() ;
      exit(type) ;
   }
   MPI_Bcast((void *) &coordType, 1, MPI_INT, root, comm) ;
}
/*--------------------------------------------------------------------*/
if ( myid == root ) {
/*
   ----------------------------------
   compute the number of entries that 
   must be sent to each other process
   ----------------------------------
*/
   outcounts = IVinit(nproc, 0) ;
   if ( nvector > 0 ) {
      int   *vecids = InpMtx_vecids(Aglobal)  ;
      int   *sizes  = InpMtx_sizes(Aglobal)   ;
   
      head = IVinit(nproc, -1) ;
      link = IVinit(nvector, -1) ;
      for ( ii = nvector - 1 ; ii >= 0 ; ii-- ) {
         vecid = vecids[ii] ;
         iproc = map[vecid] ;
         link[ii] = head[iproc] ;
         head[iproc] = ii ;
         outcounts[iproc] += sizes[ii] ;
      }
      nowned = outcounts[myid] ;
      outcounts[myid] = 0 ;
   } else {
      head = NULL ;
      link = NULL ;
      nowned = 0 ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n nowned %d", nowned) ;
      fprintf(msgFile, "\n\n outcounts") ;
      IVfprintf(msgFile, nproc, outcounts) ;
      fflush(msgFile) ;
   }
/*
   ---------------------------------------------
   scatter the outcounts to the other processors
   ---------------------------------------------
*/
   MPI_Scatter((void *) outcounts, 1, MPI_INT,
               (void *) &incount, 1, MPI_INT, root, comm) ;
/*
   ------------------------
   allocate the send buffer
   ------------------------
*/
   if ( (maxoutcount = IVmax(nproc, outcounts, &iproc)) > 0 ) {
      nDblOut = numberOfDoubles(maxoutcount, type) ;
      outbuffer = DVinit2(nDblOut) ;
   } else {
      nDblOut   =   0  ;
      outbuffer = NULL ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n nDblOut %d, outbuffer %p",
              nDblOut, outbuffer) ;
      fflush(msgFile) ;
   }
} else {
/*
   ----------------------------------
   receive the in-count from the root
   ----------------------------------
*/
   MPI_Scatter(NULL, 0, MPI_INT,
               (void *) &incount, 1, MPI_INT, root, comm) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n incount %d", incount) ;
      fflush(msgFile) ;
   }
/*
   ---------------------------
   allocate the receive buffer
   ---------------------------
*/
   if ( incount > 0 ) {
      nDblIn   = numberOfDoubles(incount, type) ;
      inbuffer = DVinit2(nDblIn) ;
   } else {
      nDblIn   =   0  ;
      inbuffer = NULL ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n nDblIn %d, inbuffer %p",
              nDblIn, inbuffer) ;
      fflush(msgFile) ;
   }
   nowned = incount ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   set up the new InpMtx object to hold the owned entries
   -------------------------------------------------------
*/
if ( Alocal == NULL ) {
   Alocal = InpMtx_new() ;
}
InpMtx_init(Alocal, coordType, type, nowned, 0) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n Alocal initialized, nowned %d", nowned) ;
   fflush(msgFile) ;
}
if ( myid == root ) {
/*
   ---------------------------------------------
   extract the triples from the original matrix,
   and load the triples into the keep matrix
   ---------------------------------------------
*/
   if ( nowned > 0 ) {
      int   *offsets = InpMtx_offsets(Aglobal) ;
      int   *sizes   = InpMtx_sizes(Aglobal)   ;
      int   ii, ient, size ;
   
      for ( ii = head[myid] ; ii != -1 ; ii = link[ii] ) {
         ient = offsets[ii] ;
         size = sizes[ii]   ;
         switch ( type ) {
         case INPMTX_INDICES_ONLY :
            InpMtx_inputTriples(Alocal, size, 
                                ivec1 + ient, ivec2 + ient);
            break ;
         case SPOOLES_REAL :
            InpMtx_inputRealTriples(Alocal, size, ivec1 + ient, 
                                    ivec2 + ient, dvec + ient) ;
            break ;
         case SPOOLES_COMPLEX :
            InpMtx_inputComplexTriples(Alocal, size, ivec1 + ient, 
                                       ivec2 + ient, dvec + 2*ient) ;
            break ;
         }
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, 
                 "\n\n Alocal after storing owned original entries") ;
         InpMtx_writeForHumanEye(Alocal, msgFile) ;
         fflush(msgFile) ;
      }
   }
/*
   ----------------------------------
   loop over the other processes, 
      gather values and send them off
      receive values
   ----------------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, 
              "\n\n ready to send entries to other processes") ;
      fflush(msgFile) ;
   }
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      if ( iproc != myid ) {
         outcount = outcounts[iproc] ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n ### process %d, send %d to %d",
                    myid, outcount, iproc) ;
            fflush(msgFile) ;
         }
         if ( outcount > 0 ) {
            double   *dbuf ;
            int      *offsets = InpMtx_offsets(Aglobal) ;
            int      *sizes   = InpMtx_sizes(Aglobal)   ;
            int      *ibuf1, *ibuf2 ;
            int      ii, ient, jent, size ;
/*
            -----------------------------------------
            load the entries to send to process iproc
            -----------------------------------------
*/
            ibuf1 = (int *) outbuffer ;
            ibuf2 = ibuf1 + outcount ;
            dbuf  = (double *) (ibuf2 + outcount) ;
            jent  = 0 ;
            for ( ii = head[iproc] ; ii != -1 ; ii = link[ii] ) {
               ient = offsets[ii] ;
               size = sizes[ii]   ;
               IVcopy(size, ibuf1 + jent, ivec1 + ient) ;
               IVcopy(size, ibuf2 + jent, ivec2 + ient) ;
               switch ( type ) {
               case SPOOLES_REAL :
                  DVcopy(size, dbuf + jent, dvec + ient) ;
                  break ;
               case SPOOLES_COMPLEX :
                  DVcopy(2*size, dbuf + 2*jent, dvec + 2*ient) ;
                  break ;
               }
               jent += size ;
            }
            if ( msglvl > 2 ) {
               fprintf(msgFile, 
                       "\n outcount %d, ibuf1 %p, ibuf2 %p, dbuf %p",
                       outcount, ibuf1, ibuf2, dbuf) ;
               writeBuffers(outcount, ibuf1, ibuf2, dbuf, type,msgFile);
               fflush(msgFile) ;
            }
            nDblOut = numberOfDoubles(outcount, type) ;
            stats[0]++ ; stats[2] += nDblOut*sizeof(double) ;
/*
            -----------------
            do a send receive
            -----------------
*/
            MPI_Send((void *) outbuffer, nDblOut, MPI_DOUBLE, 
                      iproc, firsttag, comm) ;
         }
      }
   }
/*
   ------------------------
   free the working storage
   ------------------------
*/
   if ( outbuffer != NULL ) {
      DVfree(outbuffer) ;
   }
   IVfree(outcounts) ;
   if ( link != NULL ) {
      IVfree(head) ;
      IVfree(link) ;
   }
} else {
   if ( incount > 0 ) {
      int          count ;
      int          *ibuf1 = (int *) inbuffer ;
      int          *ibuf2 = ibuf1 + incount ;
      double       *dbuf  = (double *) (ibuf2 + incount) ;
      MPI_Status   status ;
/*
      --------------------------------
      receive the buffer from the root
      --------------------------------
*/
      MPI_Recv((void *) inbuffer, nDblIn, MPI_DOUBLE, 
                root, firsttag, comm, &status) ;
      MPI_Get_count(&status, MPI_DOUBLE, &count) ;
      if ( count != nDblIn ) {
         fprintf(stderr, 
                 "\n fatal error in InpMtx_MPI_splitFromGlobal()"
                 "\n proc %d, nDblIn %d, count %d\n", 
                 myid, nDblIn, count) ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n incoming buffer") ;
         writeBuffers(incount, ibuf1, ibuf2, dbuf, type,msgFile);
         fflush(msgFile) ;
      }
/*
      -------------------------------------
      load the triples into the keep matrix
      -------------------------------------
*/
      switch ( type ) {
      case INPMTX_INDICES_ONLY :
         InpMtx_inputTriples(Alocal, incount, ibuf1, ibuf2) ;
         break ;
      case SPOOLES_REAL :
         InpMtx_inputRealTriples(Alocal, incount, ibuf1, ibuf2, dbuf) ;
         break ;
      case SPOOLES_COMPLEX :
         InpMtx_inputComplexTriples(Alocal, incount, 
                                    ibuf1, ibuf2, dbuf) ;
         break ;
      }
      if ( msglvl > 3 ) {
         fprintf(msgFile, 
                 "\n\n Alocal after storing received entries") ;
         InpMtx_writeForHumanEye(Alocal, msgFile) ;
         fflush(msgFile) ;
      }
      stats[1]++ ; stats[3] += count*sizeof(double) ;
/*
      ------------------------
      free the working storage
      ------------------------
*/
      if ( inbuffer != NULL ) {
         DVfree(inbuffer) ;
      }
   }
}
if ( Alocal != NULL ) {
/*
   ----------------------------------------------------
   sort and compress the entries and convert to vectors
   ----------------------------------------------------
*/
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n before changing storage mode to %d", 2) ;
      InpMtx_writeForHumanEye(Alocal, msgFile) ;
      fflush(msgFile) ;
   }
   InpMtx_changeStorageMode(Alocal, INPMTX_BY_VECTORS) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n after changing storage mode to %d", 2) ;
      InpMtx_writeForHumanEye(Alocal, msgFile) ;
      fflush(msgFile) ;
   }
}
return(Alocal) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   return the # of double words that 
   form a message of count triples 

   created -- 98sep27, cca
   ---------------------------------
*/
static int
numberOfDoubles (
   int   count,
   int   type
) {
int   nDbl ;

if ( sizeof(int) == sizeof(double) ) {
   nDbl = 2*count ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   nDbl = count ;
}
if ( type == SPOOLES_REAL ) {
   nDbl += count ;
} else if ( type == SPOOLES_COMPLEX ) {
   nDbl += 2*count ;
}
return(nDbl) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   write out the buffers

   created -- 98sep27, cca
   -----------------------
*/
static void
writeBuffers (
   int      count,
   int      ibuf1[],
   int      ibuf2[],
   double   dbuf[],
   int      type,
   FILE     *fp
) {
fprintf(fp, "\n ibuf1") ;
IVfprintf(fp, count, ibuf1) ;
fprintf(fp, "\n ibuf2") ;
IVfprintf(fp, count, ibuf2) ;
if ( type == SPOOLES_REAL ) {
   fprintf(fp, "\n dbuf") ;
   DVfprintf(fp, count, dbuf) ;
} else if ( type == SPOOLES_COMPLEX ) {
   fprintf(fp, "\n dbuf") ;
   DVfprintf(fp, 2*count, dbuf) ;
}

return ; }

/*--------------------------------------------------------------------*/
