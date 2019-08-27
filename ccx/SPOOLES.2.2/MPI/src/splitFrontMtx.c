/*  splitFrontMtx.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
typedef struct _msg   Msg ;
struct _msg {
   int           rowid  ;
   int           colid  ;
   int           nbytes ;
   SubMtx        *mtx   ;
   MPI_Request   req    ;
   Msg           *next  ;
} ;
static void split ( FrontMtx *frontmtx, SolveMap *solvemap, char cflag,
 int stats[], int msglvl, FILE *msgFile, int firsttag, MPI_Comm comm ) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- after the factorization has been computed and the
      front matrices have been split into submatrices, and after
      a solve map object has been computed, the submatrices
      are sent to the process that owns them.

   frontmtx -- stores the factor matrix
   solvemap -- stores the map from submatrices to processes
   created -- 98may21, cca
   -------------------------------------------------------------
*/
void
FrontMtx_MPI_split (
   FrontMtx   *frontmtx,
   SolveMap   *solvemap,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
/*
   --------------
   check the data
   --------------
*/
if ( frontmtx == NULL || solvemap == NULL 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(msgFile, "\n fatal error in FrontMtx_MPI_split()"
         "\n frontmtx %p, solvemap %p, firsttag %d"
         "\n stats %p, msglvl %d, msgFile %p"
         "\n bad input\n", 
         frontmtx, solvemap, firsttag, stats, msglvl, msgFile) ;
   exit(-1) ;
}
split(frontmtx, solvemap, 'U', stats, 
      msglvl, msgFile, firsttag, comm) ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   split(frontmtx, solvemap, 'L', stats, 
         msglvl, msgFile, firsttag, comm) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- after the factorization has been computed and the
      front matrices have been split into submatrices, and after
      a solve map object has been computed, the submatrices
      are sent to the process that owns them.

   frontmtx -- stores the factor matrix
   solvemap -- stores the map from submatrices to processes
   cflag    -- 'U' --> split the upper matrix
               'L' --> split the lower matrix

   created -- 98may21, cca
   -------------------------------------------------------------
*/
static void
split (
   FrontMtx   *frontmtx,
   SolveMap   *solvemap,
   char       cflag,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
SubMtx         *mtx ;
SubMtxManager  *mtxmanager ;
int            colid, count, destination, flag, incount, iproc, 
               J, JK, kk, K, left, maxcount, myid, nblock, nbytes, 
               nproc, offset, outcount, right, rowid, source, tag ;
int            *colids, *inbuff, *incounts, *map, *offsets, *owners, 
               *outbuff, *outcounts, *rowids ;
int            **p_inbuff, **p_outbuff ;
I2Ohash        *hash ;
Msg            *messages, *msg, *nextmsg, *recvhead, *sendhead ;
MPI_Status     status ;

if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n ##### inside FrontMtx_MPI_split()") ;
   fflush(msgFile) ;
}
/*
   --------------
   check the data
   --------------
*/
if ( frontmtx == NULL || solvemap == NULL 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(msgFile, "\n fatal error in FrontMtx_MPI_split()"
         "\n frontmtx %p, solvemap %p, firsttag %d"
         "\n stats %p, msglvl %d, msgFile %p"
         "\n bad input\n", 
         frontmtx, solvemap, firsttag, stats, msglvl, msgFile) ;
   exit(-1) ;
}
/*
   ---------------------------------
   get id of self and # of processes
   ---------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
mtxmanager = frontmtx->manager ;
if ( cflag == 'U' ) {
   hash  = frontmtx->upperhash  ;
} else {
   hash  = frontmtx->lowerhash  ;
}
if ( cflag == 'U' ) {
   nblock     = SolveMap_nblockUpper(solvemap) ;
   owners     = SolveMap_owners(solvemap) ;
   rowids     = SolveMap_rowidsUpper(solvemap) ;
   colids     = SolveMap_colidsUpper(solvemap) ;
   map        = SolveMap_mapUpper(solvemap) ;
} else if ( frontmtx->pivotingflag != 0 ) {
   nblock     = SolveMap_nblockLower(solvemap) ;
   owners     = SolveMap_owners(solvemap) ;
   rowids     = SolveMap_rowidsLower(solvemap) ;
   colids     = SolveMap_colidsLower(solvemap) ;
   map        = SolveMap_mapLower(solvemap) ;
} else {
   nblock     = SolveMap_nblockUpper(solvemap) ;
   owners     = SolveMap_owners(solvemap) ;
   colids     = SolveMap_rowidsUpper(solvemap) ;
   rowids     = SolveMap_colidsUpper(solvemap) ;
   map        = SolveMap_mapUpper(solvemap) ;
}
/*
   ------------------------------------------------
   step 1: compute the number of submatrices this 
           process will send to every other process
   ------------------------------------------------
*/
incounts  = IVinit(nproc, 0) ;
outcounts = IVinit(nproc, 0) ;
if ( cflag == 'U' ) {
   for ( JK = 0 ; JK < nblock ; JK++ ) {
      J = rowids[JK] ;
      K = colids[JK] ;
      if (  owners[J] == myid 
         && (iproc = map[JK]) != myid
         && (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
         outcounts[iproc]++ ;
      }
   }
} else {
   for ( JK = 0 ; JK < nblock ; JK++ ) {
      K = rowids[JK] ;
      J = colids[JK] ;
      if (  owners[J] == myid 
         && (iproc = map[JK]) != myid
         && (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
         outcounts[iproc]++ ;
      }
   }
}
/*
   --------------------------------------------------------
   step 2: do an all-to-all gather/scatter to get the 
           number of incoming submatrices from each process
   --------------------------------------------------------
*/
MPI_Alltoall((void *) outcounts, 1, MPI_INT,
             (void *) incounts,  1, MPI_INT, comm) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n incounts") ;
   IVfprintf(msgFile, nproc, incounts) ;
   fprintf(msgFile, "\n\n outcounts") ;
   IVfprintf(msgFile, nproc, outcounts) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------
   step 3: allocate the buffer vectors
   -----------------------------------
*/
ALLOCATE(p_outbuff, int *, nproc) ;
ALLOCATE(p_inbuff,  int *, nproc) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( outcounts[iproc] > 0 ) {
      p_outbuff[iproc] = IVinit(3*outcounts[iproc], 0) ;
   } else {
      p_outbuff[iproc] = NULL ;
   }
   if ( incounts[iproc] > 0 ) {
      p_inbuff[iproc] = IVinit(3*incounts[iproc], 0) ;
   } else {
      p_inbuff[iproc] = NULL ;
   }
}
/*
   -----------------------------------
   step 4: fill the outbuffer vectors
   ----------------------------------
*/
offsets = IVinit(nproc, 0) ;
if ( cflag == 'U' ) {
   for ( JK = kk = 0 ; JK < nblock ; JK++ ) {
      J = rowids[JK] ;
      K = colids[JK] ;
      if (  owners[J] == myid 
         && (iproc = map[JK]) != myid 
         && (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
         outbuff = p_outbuff[iproc] ;
         kk = offsets[iproc] ;
         outbuff[kk++] = J ;
         outbuff[kk++] = K ;
         outbuff[kk++] = SubMtx_nbytesInUse(mtx) ;
         offsets[iproc] = kk ;
      }
   }
} else {
   for ( JK = kk = 0 ; JK < nblock ; JK++ ) {
      K = rowids[JK] ;
      J = colids[JK] ;
      if (  owners[J] == myid 
         && (iproc = map[JK]) != myid 
         && (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
         outbuff = p_outbuff[iproc] ;
         kk = offsets[iproc] ;
         outbuff[kk++] = K ;
         outbuff[kk++] = J ;
         outbuff[kk++] = SubMtx_nbytesInUse(mtx) ;
         offsets[iproc] = kk ;
      }
   }
}
IVfree(offsets) ;
/*
   -------------------------------------------------
   step 5: send/receive the buffer vectors in stages
   -------------------------------------------------
*/
for ( offset = 1, tag = firsttag ; offset < nproc ; offset++, tag++ ) {
   right    = (myid + offset) % nproc ;
   left     = (nproc + myid - offset) % nproc ;
   outcount = outcounts[right] ;
   incount  = incounts[left] ;
   if ( msglvl > 1 ) {
      fprintf(msgFile,
         "\n ### process %d, send %d to right %d, recv %d from left %d",
              myid, outcount, right, incount, left) ;
      fflush(msgFile) ;
   }
   if ( outcount > 0 ) {
      destination = right ;
      stats[0]++ ;
      stats[2] += 3*outcount*sizeof(int) ;
      outbuff = p_outbuff[right] ;
   } else {
      destination = MPI_PROC_NULL ;
      outbuff = NULL ;
   }
   if ( incount > 0 ) {
      source = left ;
      stats[1]++ ;
      stats[3] += 3*incount*sizeof(int) ;
      inbuff = p_inbuff[left] ;
   } else {
      source = MPI_PROC_NULL ;
      inbuff = NULL ;
   }
/*
   -----------------
   do a send/receive
   -----------------
*/
   if ( msglvl > 1 ) {
      fprintf(msgFile,
         "\n ### Sendrecv: dest %d, tag %d, source %d, tag %d",
              destination, tag, source, tag) ;
      fflush(msgFile) ;
   }
   MPI_Sendrecv(outbuff, 3*outcount, MPI_INT, destination, tag,
                inbuff,  3*incount,  MPI_INT, source,      tag,
                comm, &status) ;
   if ( source != MPI_PROC_NULL ) {
      MPI_Get_count(&status, MPI_INT, &count) ;
      if ( count != 3*incount ) {
         fprintf(stderr,
                 "\n 1. fatal error in FrontMtx_MPI_split()"
                 "\n proc %d : source = %d, count = %d, incount = %d\n",
                 myid, source, count, 3*incount) ;
         exit(-1) ;
      }
   }
   if ( incount > 0 && msglvl > 1 ) {
      fprintf(msgFile, "\n inbuffer from proc %d", left) ;
      for ( kk = 0 ; kk < incount ; kk++ ) {
         fprintf(msgFile, "\n %8d %8d %8d", 
                 inbuff[3*kk], inbuff[3*kk+1], inbuff[3*kk+2]) ;
      }
      fflush(msgFile) ;
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile,
      "\n ### outbuff[] and inbuff[] sent and received\n") ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   step 6: compute the maximum number of messages 
           sent and received in any stage
   ----------------------------------------------
*/
maxcount = 0 ;
for ( offset = 1 ; offset < nproc ; offset++ ) {
   right    = (myid + offset) % nproc ;
   left     = (nproc + myid - offset) % nproc ;
   outcount = outcounts[right] ;
   incount  = incounts[left] ;
   if ( maxcount < incount + outcount ) {
      maxcount = incount + outcount ;
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n maxcount = %d", maxcount) ;
   fflush(msgFile) ;
}
ALLOCATE(messages, struct _msg, maxcount) ;
/*
   ----------------------------------------------
   step 8: send/receive the submatrices in stages
   ----------------------------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n starting to send/receive the submatrices") ;
   fflush(msgFile) ;
}
tag = firsttag ;
for ( offset = 1 ; offset < nproc ; offset++ ) {
   right    = (myid + offset) % nproc ;
   left     = (nproc + myid - offset) % nproc ;
   outcount = outcounts[right] ;
   incount  = incounts[left] ;
   if ( msglvl > 1 ) {
      fprintf(msgFile,
         "\n ### process %d, send %d to right %d, recv %d from left %d",
              myid, outcount, right, incount, left) ;
      fflush(msgFile) ;
   }
   msg = messages ;
/*
   ---------------------------------------------------
   post the receives for submatrices from process left
   ---------------------------------------------------
*/
   inbuff = p_inbuff[left] ;
   for ( JK = 0, recvhead = NULL ; JK < incount ; JK++, msg++ ) {
      rowid  = inbuff[3*JK]   ;
      colid  = inbuff[3*JK+1] ;
      nbytes = inbuff[3*JK+2] ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
                 "\n ready to receive mtx(%d,%d) with %d bytes from %d",
                 rowid, colid, nbytes, left) ;
         fflush(msgFile) ;
      }
      mtx = SubMtxManager_newObjectOfSizeNbytes(mtxmanager, nbytes);
      I2Ohash_insert(hash, rowid, colid, mtx) ;
      msg->rowid  = rowid    ;
      msg->colid  = colid    ;
      msg->mtx    = mtx      ;
      msg->nbytes = nbytes   ;
      msg->next   = recvhead ;
      recvhead    = msg ;
      if ( msglvl > 1 ) {
         fprintf(msgFile,
                 "\n ### posting Irecv, left %d, tag %d, nbytes %d",
                 left, JK + firsttag, nbytes) ;
         fflush(msgFile) ;
      }
      MPI_Irecv(SubMtx_workspace(mtx), nbytes, MPI_BYTE,
                left, JK + firsttag, comm, &msg->req) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, ", Irecv posted") ;
         fflush(msgFile) ;
      }
      stats[1]++ ;
      stats[3] += nbytes ;
   }
/*
   -----------------------------------------------------
   post the sends for submatrices going to process right
   -----------------------------------------------------
*/
   outbuff = p_outbuff[right] ;
   for ( JK = kk = 0, sendhead = NULL ; JK < outcount ; JK++, msg++ ) {
      rowid  = outbuff[3*JK]   ;
      colid  = outbuff[3*JK+1] ;
      nbytes = outbuff[3*JK+2] ;
      I2Ohash_remove(hash, rowid, colid, (void *) &mtx) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, 
                 "\n ready to send mtx(%d,%d) with %d bytes to %d",
                 rowid, colid, nbytes, right) ;
         fflush(msgFile) ;
      }
      msg->rowid  = rowid    ;
      msg->colid  = colid    ;
      msg->mtx    = mtx      ;
      msg->nbytes = nbytes   ;
      msg->next   = sendhead ;
      sendhead    = msg ;
      if ( msglvl > 1 ) {
         fprintf(msgFile,
                 "\n ### posting Isend, right %d, tag %d, nbytes %d",
                 right, JK + firsttag, nbytes) ;
         SubMtx_writeForHumanEye(msg->mtx, msgFile) ;
         fprintf(msgFile, "\n first seven entries of buffer") ;
         IVfprintf(msgFile, 7, (int *) SubMtx_workspace(msg->mtx)) ;
         fflush(msgFile) ;
      }
      MPI_Isend(SubMtx_workspace(mtx), nbytes, MPI_BYTE,
                right, JK + firsttag, comm, &msg->req) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n Isend posted") ;
         fflush(msgFile) ;
      }
      stats[0]++ ;
      stats[2] += nbytes ;
   }
/*
   -----------------------------------------
   loop while send/receives are not complete
   -----------------------------------------
*/
   while ( sendhead != NULL || recvhead != NULL ) {
      for ( msg = sendhead, sendhead = NULL ; 
            msg != NULL ; 
            msg = nextmsg ) {
         nextmsg = msg->next ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                    "\n testing sent msg %p, (%d,%d) with %d bytes",
                    msg, msg->rowid, msg->colid, msg->nbytes) ;
            fflush(msgFile) ;
         }
         MPI_Test(&msg->req, &flag, &status) ;
         if ( flag == 1 ) {
/*
            ---------------------------------------------
            message has been received, release the matrix
            ---------------------------------------------
*/
            if ( msglvl > 1 ) {
               fprintf(msgFile, ", send complete, releasing mtx") ;
               fflush(msgFile) ;
            }
            SubMtxManager_releaseObject(mtxmanager, msg->mtx) ;
         } else {
/*
            ---------------------------------------------------
            message has not been received, keep message in list
            ---------------------------------------------------
*/
            msg->next = sendhead ;
            sendhead  = msg ;
         }
      }
      for ( msg = recvhead, recvhead = NULL ; 
            msg != NULL ; 
            msg = nextmsg ) {
         nextmsg = msg->next ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, 
                    "\n testing recv msg %p, (%d,%d) with %d bytes",
                    msg, msg->rowid, msg->colid, msg->nbytes) ;
            fflush(msgFile) ;
         }
         MPI_Test(&msg->req, &flag, &status) ;
         if ( flag == 1 ) {
/*
            ------------------------------------------------
            message has been received, initialize the matrix
            ------------------------------------------------
*/
            if ( msglvl > 1 ) {
               fprintf(msgFile, ", recv complete, initializing mtx") ;
               fprintf(msgFile, "\n first seven entries of buffer") ;
               IVfprintf(msgFile, 7, (int *) SubMtx_workspace(msg->mtx)) ;
               fflush(msgFile) ;
            }
            SubMtx_initFromBuffer(msg->mtx) ;
            if ( msglvl > 1 ) {
               SubMtx_writeForHumanEye(msg->mtx, msgFile) ;
               fflush(msgFile) ;
            }
         } else {
/*
            ---------------------------------------------------
            message has not been received, keep message in list
            ---------------------------------------------------
*/
            msg->next = recvhead ;
            recvhead  = msg ;
         }
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(incounts) ;
IVfree(outcounts) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( (inbuff = p_inbuff[iproc]) != NULL ) {
      IVfree(inbuff) ;
   }
   if ( (outbuff = p_outbuff[iproc]) != NULL ) {
      IVfree(outbuff) ;
   }
}
FREE(p_inbuff) ;
FREE(p_outbuff) ;
FREE(messages) ;

return ; }

/*--------------------------------------------------------------------*/
