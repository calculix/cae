/*  fullAdjMPI.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- given a distributed InpMtx object,
      gather all the indices onto each process and create
      an IVL object that contains the full adjacency structure

   created -- 97dec17, cca
   -----------------------------------------------------------
*/
IVL *
InpMtx_MPI_fullAdjacency (
   InpMtx    *inpmtx,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
InpMtx   *adjmtx ;
int       ierr, iproc, maxnent, myid, nent, nproc, oldtype, totalnent ;
int       *buffer, *counts, *ivec1, *ivec2 ;
IVL       *adjIVL ;
/*
   --------------------------------------
   get id of self and number of processes
   --------------------------------------
*/
MPI_Comm_rank(comm, &myid)  ;
MPI_Comm_size(comm, &nproc) ;
/*
   ----------------------
   convert to row storage
   ----------------------
*/
oldtype = InpMtx_coordType(inpmtx) ;
InpMtx_changeCoordType(inpmtx, INPMTX_BY_ROWS) ;
nent  = InpMtx_nent(inpmtx) ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %d internal entries", nent) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   find out how many entries each process owns
   -------------------------------------------
*/
counts = IVinit(nproc, 0) ;
counts[myid] = nent ;
MPI_Allgather((void *) &counts[myid], 1, MPI_INT,
              (void *) counts, 1, MPI_INT, comm) ;
totalnent = IVsum(nproc, counts) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %d total entries", totalnent) ;
   fprintf(msgFile, "\n\n counts vector") ;
   IVfp80(msgFile, nproc, counts, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   allocate a new InpMtx object to hold the indices
   -------------------------------------------------
*/
adjmtx = InpMtx_new() ;
InpMtx_init(adjmtx, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, totalnent, 0) ;
/*
   -----------------------------------------
   allocate a buffer to send/receive entries
   -----------------------------------------
*/
maxnent = IVmax(nproc, counts, &iproc) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %d maximum entries", maxnent) ;
   fflush(msgFile) ;
}
buffer = IVinit(2*maxnent, -1) ;
/*
   ----------------------------
   now send and receive entries
   ----------------------------
*/
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   nent = counts[iproc] ;
   if ( iproc == myid ) {
/*
      -----------------------
      load the entries buffer
      -----------------------
*/
      IVcopy(nent, buffer, ivec1) ;
      IVcopy(nent, buffer + nent, ivec2) ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n owned entries in buffer") ;
         fflush(msgFile) ;
      }
      if ( msglvl > 2 ) {
         IVfprintf(msgFile, 2*nent, buffer) ;
         fflush(msgFile) ;
      }
      stats[0]++ ;
      stats[2] += 2*nent*sizeof(int) ;
   } else {
      stats[1]++ ;
      stats[3] += 2*nent*sizeof(int) ;
   }
/*
   ----------------------------------------
   broadcast the entries from process iproc
   ----------------------------------------
*/
   MPI_Bcast((void *) buffer, 2*nent, MPI_INT, iproc, comm) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n after bcast, buffer") ;
      IVfprintf(msgFile, 2*nent, buffer) ;
      fflush(msgFile) ;
   }
/*
   ------------------------------------------
   load the entries into the adjacency matrix
   ------------------------------------------
*/
   InpMtx_inputTriples(adjmtx, nent, buffer, buffer + nent) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n raw InpMtxobject") ;
   InpMtx_writeForHumanEye(adjmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   sort and compress the adjacency entries
   and convert to vector form
   ---------------------------------------
*/
InpMtx_sortAndCompress(adjmtx) ;
InpMtx_changeStorageMode(adjmtx, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n sorted InpMtxobject") ;
   InpMtx_writeForHumanEye(adjmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------
   create the IVL object that contains the full adjacency
   ------------------------------------------------------
*/
adjIVL = InpMtx_fullAdjacency(adjmtx) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n full adjacency object") ;
   IVL_writeForHumanEye(adjIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------
   convert back to original storage
   --------------------------------
*/
InpMtx_changeCoordType(inpmtx, oldtype) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(counts) ;
IVfree(buffer) ;
InpMtx_free(adjmtx) ;

return(adjIVL) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- given a distributed Pencil object,
      gather all the indices onto each process and create
      an IVL object that contains the full adjacency structure

   created -- 97dec18, cca
   -----------------------------------------------------------
*/
IVL *
Pencil_MPI_fullAdjacency (
   Pencil    *pencil,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
InpMtx   *adjmtx, *inpmtxA, *inpmtxB ;
int       ierr, iproc, maxnent, myid, nent, nentA, nentB,
          nproc, oldtypeA, oldtypeB, totalnent ;
int       *buffer, *counts, *ivec1, *ivec2, *tempbuffer ;
IVL       *adjIVL ;
/*
   ------------------------------------------------------
   check for a simple pencil (only one nontrivial matrix)
   ------------------------------------------------------
*/
inpmtxA = pencil->inpmtxA ;
inpmtxB = pencil->inpmtxB ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n inside Pencil_MPI_fullAdjacency(), A = %p, B = %p",
           inpmtxA, inpmtxB) ;
   fprintf(msgFile, "\n sigma = [%12.4e,%12.4e]", 
           pencil->sigma[0], pencil->sigma[1]) ;
   fflush(msgFile) ;
}
if ( inpmtxA == NULL ) {
   if ( pencil->sigma[0] != 0.0 && inpmtxB != NULL ) {
      adjIVL = InpMtx_MPI_fullAdjacency(inpmtxB, stats, 
                                         msglvl, msgFile, comm) ;
   } else {
      adjIVL = NULL ;
   }
   return(adjIVL) ;
} else if ( (pencil->sigma[0] == 0.0 && pencil->sigma[1] == 0.0)
          || inpmtxB == NULL ) {
   adjIVL = InpMtx_MPI_fullAdjacency(inpmtxA, 
                                      stats, msglvl, msgFile, comm) ;
   return(adjIVL) ;
}
/*
   --------------------------------------
   get id of self and number of processes
   --------------------------------------
*/
MPI_Comm_rank(comm, &myid)  ;
MPI_Comm_size(comm, &nproc) ;
/*
   ---------------------------------------------
   allocate a send buffer, fill with the indices 
   from the two matrices, sort and compress
   ---------------------------------------------
*/
oldtypeA = InpMtx_coordType(inpmtxA) ;
InpMtx_changeCoordType(inpmtxA, INPMTX_BY_ROWS) ;
oldtypeB = InpMtx_coordType(inpmtxB) ;
InpMtx_changeCoordType(inpmtxB, INPMTX_BY_ROWS) ;
nentA  = InpMtx_nent(inpmtxA)  ;
nentB  = InpMtx_nent(inpmtxB)  ;
if ( nentA > 0 || nentB > 0 ) {
   tempbuffer = IVinit(2*(nentA + nentB), -1) ;
   ivec1 = tempbuffer ;
   ivec2 = ivec1 + nentA + nentB ;
   if ( nentA > 0 ) {
      IVcopy(nentA, ivec1, InpMtx_ivec1(inpmtxA)) ;
      IVcopy(nentA, ivec2, InpMtx_ivec2(inpmtxA)) ;
   }
   if ( nentB > 0 ) {
      IVcopy(nentB, ivec1 + nentA, InpMtx_ivec1(inpmtxB)) ;
      IVcopy(nentB, ivec2 + nentA, InpMtx_ivec2(inpmtxB)) ;
   }
   if ( msglvl > 5 ) {
      fprintf(msgFile, "\n\n before sort and compress") ;
      fprintf(msgFile, "\n ivec1") ;
      IVfprintf(msgFile, nentA + nentB, ivec1) ;
      fprintf(msgFile, "\n ivec2") ;
      IVfprintf(msgFile, nentA + nentB, ivec2) ;
      fflush(msgFile) ;
   }
   nent = IV2sortUpAndCompress(nentA + nentB, ivec1, ivec2) ;
   if ( msglvl > 5 ) {
      fprintf(msgFile, "\n\n after sort and compress, nent = %d", nent);
      fprintf(msgFile, "\n ivec1") ;
      IVfprintf(msgFile, nent, ivec1) ;
      fprintf(msgFile, "\n ivec2") ;
      IVfprintf(msgFile, nent, ivec2) ;
      fflush(msgFile) ;
   }
} else {
   nent       = 0 ;
   tempbuffer = NULL ;
   ivec1      = NULL ;
   ivec2      = NULL ;
}
if ( msglvl > 5 ) {
   fprintf(msgFile, "\n\n %d internal entries", nent) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   find out how many entries each process owns
   -------------------------------------------
*/
counts = IVinit(nproc, 0) ;
counts[myid] = nent ;
MPI_Allgather((void *) &counts[myid], 1, MPI_INT,
              (void *) counts, 1, MPI_INT, comm) ;
totalnent = IVsum(nproc, counts) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %d total entries", totalnent) ;
   fprintf(msgFile, "\n\n counts vector") ;
   IVfp80(msgFile, nproc, counts, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------
   allocate a new InpMtx object to hold the indices
   -------------------------------------------------
*/
adjmtx = InpMtx_new() ;
InpMtx_init(adjmtx, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, totalnent, 0) ;
/*
   -----------------------------------------
   allocate a buffer to send/receive entries
   and load with the owned entries
   -----------------------------------------
*/
maxnent = IVmax(nproc, counts, &iproc) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n %d maximum entries", maxnent) ;
   fflush(msgFile) ;
}
buffer = IVinit(2*maxnent, -1) ;
/*
   ----------------------------
   now send and receive entries
   ----------------------------
*/
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   nent = counts[iproc] ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n processor %d's turn with %d entries",
              iproc, nent) ;
      fflush(msgFile) ;
   }
   if ( nent > 0 ) {
      if ( iproc == myid ) {
/*
         -----------------------
         load the entries buffer
         -----------------------
*/
         IVcopy(nent, buffer, ivec1) ;
         IVcopy(nent, buffer + nent, ivec2) ;
         if ( msglvl > 1 ) {
            fprintf(msgFile, "\n\n owned entries in buffer") ;
            fflush(msgFile) ;
         }
         if ( msglvl > 2 ) {
            IVfprintf(msgFile, 2*nent, buffer) ;
            fflush(msgFile) ;
         }
         stats[0]++ ;
         stats[2] += 2*nent*sizeof(int) ;
      } else {
         stats[1]++ ;
         stats[3] += 2*nent*sizeof(int) ;
      }
/*
      ----------------------------------------
      broadcast the entries from process iproc
      ----------------------------------------
*/
      MPI_Bcast((void *) buffer, 2*nent, MPI_INT, iproc, comm) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n after bcast, buffer") ;
         IVfprintf(msgFile, 2*nent, buffer) ;
         fflush(msgFile) ;
      }
/*
      ------------------------------------------
      load the entries into the adjacency matrix
      ------------------------------------------
*/
      InpMtx_inputTriples(adjmtx, nent, buffer, buffer + nent) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n entries from buffer loaded") ;
         fflush(msgFile) ;
      }
   }
}
/*
   ---------------------------------------
   sort and compress the adjacency entries
   ---------------------------------------
*/
InpMtx_sortAndCompress(adjmtx) ;
InpMtx_changeStorageMode(adjmtx, INPMTX_BY_VECTORS) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n adjmtx") ;
   InpMtx_writeForHumanEye(adjmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------
   create the IVL object that contains the full adjacency
   ------------------------------------------------------
*/
adjIVL = InpMtx_fullAdjacency(adjmtx) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n full adjacency object") ;
   IVL_writeForHumanEye(adjIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------
   convert back to original storage
   --------------------------------
*/
InpMtx_changeCoordType(inpmtxA, oldtypeA) ;
InpMtx_changeCoordType(inpmtxB, oldtypeB) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(counts) ;
if ( tempbuffer != NULL ) {
   IVfree(tempbuffer) ;
}
IVfree(buffer) ;
InpMtx_free(adjmtx) ;

return(adjIVL) ; }

/*--------------------------------------------------------------------*/
