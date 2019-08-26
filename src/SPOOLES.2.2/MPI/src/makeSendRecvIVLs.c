/*  makeSendRecvIVLs.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to analyze and organize communication. it was written
      in support of a distributed matrix-vector multiply but can be 
      used for other applications.

   each processor has a list of items it "supports" or needs found
   in the supportedIV object. the globalmapIV object contains the
   map from items to owning processors. we need to figure out what
   items this processor will send to and receive from each other
   processor. this information is found in the sendIVL and recvIVL
   objects. 

   on return, list jproc of sendIVL contains the items owned by
   this processor and needed by jproc.
   on return, list jproc of recvIVL contains the items needed by
   this processor and owned by jproc.

   as a concrete example, consider a distributed Y = A * X.
   the matrix A, the right hand side X and the vector Y are
   distributed among processors. 

   consider the case where the supportedIV object contains the rows
   of X that are needed by this processor to perform its part of the
   matrix-vector multiply. globalmapIV contains the map from rows
   of X to the owning processors. on return, list jproc of sendIVL 
   contains the row indices of X owned by this processor that are
   needed by processor jproc. on return, list jproc of recvIVL 
   contains the row indices of X needed by this processor that are
   owned by processor jproc.

   consider the case where the supportedIV object contains the rows
   of Y that will be updated by this processor when it performs it
   part of the matrix-vector multiply. globalmapIV contains the map
   from rows of Y to their owning processors. on return, list jproc
   of recvIVL contains the row indices of Y on this processor that 
   need to be sent to their owner, processor jproc. on return, list 
   jproc of sendIVL contains the row indices of Y owned by this 
   processor that will be sent by processor jproc to this processor.

   created -- 98aug01, cca
   -----------------------------------------------------------------
*/
void
makeSendRecvIVLs (
   IV         *supportedIV,
   IV         *globalmapIV,
   IVL        *sendIVL,
   IVL        *recvIVL,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int   count, ii, item, jproc, maxitem, myid, nitem, nproc ;
int   *head, *items, *link, *list, *map ;
/*
   ---------------
   check the input
   ---------------
*/
if ( supportedIV == NULL || globalmapIV == NULL 
   || sendIVL == NULL || recvIVL == NULL 
   || stats == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, "\n fatal error in makeSendRecvIVLs()"
           "\n bad input\n") ;
   exit(-1) ;
}
IV_sizeAndEntries(supportedIV, &nitem,  &items) ;
if ( nitem == 0 ) {
   maxitem = 0 ;
} else { 
   maxitem = items[nitem-1] ;
}
map = IV_entries(globalmapIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside makeSendRecvIVLs()"
           "\n supportedIV") ;
   IV_writeForHumanEye(supportedIV, msgFile) ;
   fprintf(msgFile, "\n globalmapIV") ;
   IV_writeForHumanEye(globalmapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   get id of self and number of processors
   ---------------------------------------
*/
MPI_Comm_rank(comm, &myid)  ;
MPI_Comm_size(comm, &nproc) ;
/*
   ----------------------------------------------------
   link the items into lists via their owning processor
   ----------------------------------------------------
*/
head = IVinit(nproc, -1) ;
link = IVinit(1 + maxitem, -1) ;
for ( ii = 0 ; ii < nitem ; ii++ ) {
   item = items[ii] ;
   jproc = map[item] ;
   link[item] = head[jproc] ;
   head[jproc] = item ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n items linked by owning processor") ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   initialize and fill the recvIVL object.
   list jproc contains the items that 
   this processor needs from jproc.
   ---------------------------------------
*/
IVL_init1(recvIVL, IVL_CHUNKED, nproc) ;
if ( nitem > 0 ) {
   list = IVinit(nitem, -1) ;
   for ( jproc = 0 ; jproc < nproc ; jproc++ ) {
      count = 0 ;
      for ( item = head[jproc] ; item != -1 ; item = link[item] ) {
         list[count++] = item ;
      }
      IVqsortUp(count, list) ;
      IVL_setList(recvIVL, jproc, count, list) ;
   }
   IVfree(list) ;
   IVfree(head) ;
   IVfree(link) ;
}
if ( msglvl > 5 ) {
   fprintf(msgFile, "\n\n recvIVL") ;
   IVL_writeForHumanEye(recvIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------------
   compute the sendIVL object via an all-to-all communication
   ----------------------------------------------------------
*/
IVL_MPI_alltoall(recvIVL, sendIVL, stats, 
                 msglvl, msgFile, firsttag, comm) ;

return ; }

/*--------------------------------------------------------------------*/
