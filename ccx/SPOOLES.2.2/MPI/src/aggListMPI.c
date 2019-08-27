/*  aggList.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   create, initialize and return a ChvList object
   to deal with aggregate chevrons
      
   created  -- 98may21, cca
   modified -- 98jul31, cca
      now uses IVL_MPI_alltoall()
   -----------------------------------------------
*/
ChvList *
FrontMtx_MPI_aggregateList (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        tag,
   MPI_Comm   comm
) {
char      *mark ;
ChvList   *aggList ;
int       count, ierr, ii, jproc, J, K, myid, nfront, nproc, size ;
int       *aggcounts, *frontOwners, *head, *indices, *link, *list,
          *updated, *vtxToFront ;
IVL       *recvIVL, *symbfacIVL, *sendIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || frontOwnersIV == NULL ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_MPI_aggregateList(%p,%p,%p)"
           "\n bad input\n", frontmtx, frontOwnersIV, comm) ;
   exit(-1) ;
}
if ( tag < 0 || tag > maxTagMPI(comm) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_aggregateList()"
           "\n tag = %d, tag_bound = %d", tag, maxTagMPI(comm)) ;
   exit(-1) ;
}
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
symbfacIVL = frontmtx->symbfacIVL ;
vtxToFront = ETree_vtxToFront(frontmtx->frontETree) ;
IV_sizeAndEntries(frontOwnersIV, &nfront, &frontOwners) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n\n inside FrontMtx_aggListMPI, myid = %d, nproc = %d",
           myid, nproc) ; 
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------
   mark all fronts that are supported by an owned front
   collect lists of updated fronts by owning processor
   ----------------------------------------------------
*/
mark = CVinit(nfront, 'N') ;
head = IVinit(nproc,  -1) ;
link = IVinit(nfront, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   jproc = frontOwners[J] ;
   if ( jproc == myid ) {
      IVL_listAndSize(symbfacIVL, J, &size, &indices) ;
      for ( ii = 0 ; ii < size ; ii++ ) {
         K = vtxToFront[indices[ii]] ;
         if ( mark[K] == 'N' ) {
            mark[K] = 'Y' ;
            jproc   = frontOwners[K] ;
            link[K] = head[jproc] ;
            head[jproc] = K ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, "\n front %d supported", K) ;
               fflush(msgFile) ;
            }
         }
      }
   }
}
/*
   -------------------------
   set up the sendIVL object
   -------------------------
*/
list = IVinit(nfront, -1) ;
sendIVL = IVL_new() ;
IVL_init1(sendIVL, IVL_CHUNKED, nproc) ;
for ( jproc = 0 ; jproc < nproc ; jproc++ ) {
   for ( K = head[jproc], count = 0 ; K != -1 ; K = link[K] ) {
      list[count++] = K ;
   }
   IVL_setList(sendIVL, jproc, count, list) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n send IVL for aggregate lists") ;
   IVL_writeForHumanEye(sendIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------
   get the recvIVL object
   ----------------------
*/
recvIVL = IVL_MPI_alltoall(sendIVL, NULL, 
                           stats, msglvl, msgFile, tag, comm) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n receive IVL for aggregate lists") ;
   IVL_writeForHumanEye(recvIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   fill the aggcounts vector
   -------------------------
*/
aggcounts = IVinit(nfront, 0) ;
for ( jproc = 0 ; jproc < nproc ; jproc++ ) {
   if ( jproc != myid ) {
      IVL_listAndSize(recvIVL, jproc, &count, &updated) ;
      for ( ii = 0 ; ii < count ; ii++ ) {
         aggcounts[updated[ii]]++ ;
      }
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n aggcounts") ;
   IVfp80(msgFile, nfront, aggcounts, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------
   create and initialize the ChvList object
   -----------------------------------------
*/
aggList = ChvList_new() ;
ChvList_init(aggList, nfront, aggcounts, 0, NULL) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(aggcounts) ;
IVfree(head) ;
IVfree(link) ;
IVfree(list) ;
CVfree(mark) ;
IVL_free(sendIVL) ;
IVL_free(recvIVL) ;

return(aggList) ; }

/*--------------------------------------------------------------------*/
