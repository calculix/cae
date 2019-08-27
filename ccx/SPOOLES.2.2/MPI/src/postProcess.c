/*  postProcess.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- post-process the factorization
      (1) permute row and column adjacency objects if necessary
      (2) permute lower and upper matrices if necessary
      (3) update the block adjacency objects if necessary
      (4) split the chevron submatrices into submatrices
          and make the submatrix indices local w.r.t their fronts

   created -- 98may20, cca
   --------------------------------------------------------------
*/
void
FrontMtx_MPI_postProcess (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int   lasttag, nfront, nproc, tagbound ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || frontOwnersIV == NULL || stats == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_MPI_postProcess()"
           "\n frontmtx %p, frontOwnersIV %p, firsttag %d"
           "\n stats %p, msglvl %d, msgFile %p, comm %p"
           "\n bad input\n", frontmtx, frontOwnersIV, firsttag, 
           stats, msglvl, msgFile, comm) ;
   exit(-1) ;
}
MPI_Comm_size(comm, &nproc) ;
lasttag  = firsttag + 5*nproc ;
tagbound = maxTagMPI(comm) ;
if ( firsttag < 0 || lasttag > tagbound ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_postProcess()"
           "\n firsttag = %d, tagbound = %d", firsttag, tagbound) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n nfront = %d, pivotingflag = %d",
           nfront, frontmtx->pivotingflag) ;
   fflush(msgFile) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   IV   *colmapIV, *rowmapIV ;
/*
   --------------------------------------
   gather the global frontsizes IV object
   --------------------------------------
*/
   IV_MPI_allgather(frontmtx->frontsizesIV, frontOwnersIV, 
                    stats, msglvl, msgFile, firsttag, comm) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n global frontsizes IV object") ;
      IV_writeForHumanEye(frontmtx->frontsizesIV, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -------------------------------
   permute the adjacency object(s)
   -------------------------------
*/
   FrontMtx_MPI_permuteUpperAdj(frontmtx, frontOwnersIV, stats, 
                                msglvl, msgFile, firsttag, comm) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n new column adjacency object") ;
      IVL_writeForHumanEye(frontmtx->coladjIVL, msgFile) ;
      fflush(msgFile) ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      FrontMtx_MPI_permuteLowerAdj(frontmtx, frontOwnersIV, stats,
                                   msglvl, msgFile, firsttag, comm) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n new row adjacency object") ;
         IVL_writeForHumanEye(frontmtx->rowadjIVL, msgFile) ;
         fflush(msgFile) ;
      }
   }
/*
   -------------------------------------------------------------
   permute the U_{J,bnd{J}} and L_{bnd{J},J} triangular matrices
   -------------------------------------------------------------
*/
   FrontMtx_permuteUpperMatrices(frontmtx, msglvl, msgFile) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      FrontMtx_permuteLowerMatrices(frontmtx, msglvl, msgFile) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n front factor matrix after pivoting") ;
      FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   }
/*
   -----------------------------------------------
   get the map from columns to owning fronts
   and create the new upper block adjacency object
   -----------------------------------------------
*/
   colmapIV = FrontMtx_colmapIV(frontmtx) ;
   frontmtx->upperblockIVL 
                      = FrontMtx_makeUpperBlockIVL(frontmtx, colmapIV) ;
   IV_free(colmapIV) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      -------------------------------------------
      get the map from rows to owning fronts and
      create the new lower block adjacency object
      -------------------------------------------
*/
      rowmapIV = FrontMtx_rowmapIV(frontmtx) ;
      frontmtx->lowerblockIVL 
                      = FrontMtx_makeLowerBlockIVL(frontmtx, rowmapIV) ;
      IV_free(rowmapIV) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n local upper block adjacency object") ;
      IVL_writeForHumanEye(frontmtx->upperblockIVL, msgFile) ;
      if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
         fprintf(msgFile, "\n\n local lower block adjacency object") ;
         IVL_writeForHumanEye(frontmtx->lowerblockIVL, msgFile) ;
      }
      fflush(msgFile) ;
   }
} else {
/*
   ---------------------------------------
   get the upper block adjacency structure
   ---------------------------------------
*/
   IV *vtxToFrontIV = ETree_vtxToFrontIV(frontmtx->frontETree) ;
   frontmtx->upperblockIVL 
                  = FrontMtx_makeUpperBlockIVL(frontmtx, vtxToFrontIV) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      frontmtx->lowerblockIVL 
                  = FrontMtx_makeUpperBlockIVL(frontmtx, vtxToFrontIV) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n local upper block adjacency object") ;
      IVL_writeForHumanEye(frontmtx->upperblockIVL, msgFile) ;
      if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
         fprintf(msgFile, "\n\n local lower block adjacency object") ;
         IVL_writeForHumanEye(frontmtx->lowerblockIVL, msgFile) ;
      }
      fflush(msgFile) ;
   }
}
/*
   --------------------------------------
   all-gather the block adjacency objects
   --------------------------------------
*/
IVL_MPI_allgather(frontmtx->upperblockIVL, frontOwnersIV, 
                  stats, msglvl, msgFile, firsttag, comm) ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   IVL_MPI_allgather(frontmtx->lowerblockIVL, frontOwnersIV, 
                     stats, msglvl, msgFile, firsttag, comm) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n global upper block adjacency object") ;
   IVL_writeForHumanEye(frontmtx->upperblockIVL, msgFile) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      fprintf(msgFile, "\n\n global lower block adjacency object") ;
      IVL_writeForHumanEye(frontmtx->lowerblockIVL, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   ------------------------
   allocate the hash tables
   ------------------------
*/
frontmtx->upperhash = I2Ohash_new() ;
I2Ohash_init(frontmtx->upperhash, nfront, nfront, nfront) ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   frontmtx->lowerhash = I2Ohash_new() ;
   I2Ohash_init(frontmtx->lowerhash, nfront, nfront, nfront) ;
} else {
   frontmtx->lowerhash = NULL ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n hash tables allocated") ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------------
   split the U_{J,bnd{J}} and L_{bnd{J},J} into submatrices
   put the U_{J,K} and L_{K,J} matrices into hash tables,
   free the p_mtx*[] vectors.
   --------------------------------------------------------
*/
FrontMtx_splitUpperMatrices(frontmtx, msglvl, msgFile) ;
FREE(frontmtx->p_mtxUJJ) ; frontmtx->p_mtxUJJ = NULL ;
FREE(frontmtx->p_mtxUJN) ; frontmtx->p_mtxUJN = NULL ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   FrontMtx_splitLowerMatrices(frontmtx, msglvl, msgFile) ;
   FREE(frontmtx->p_mtxLJJ) ; frontmtx->p_mtxLJJ = NULL ;
   FREE(frontmtx->p_mtxLNJ) ; frontmtx->p_mtxLNJ = NULL ;
}
frontmtx->dataMode = FRONTMTX_2D_MODE ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n matrices split into submatrices") ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to permute the indices of the upper adjacency
      structure so that for each front J, bnd{J} is in 
      ascending order w.r.t. K cup bnd{K} for par[J] = K.
      process q sends to process r one message that contains
      J cup bnd{J} for all J owned by q and needed by r.
      once all the indices for the supported fronts are present,
      the indices in the upper adjacency structure are reordered
      as necessary.

   created -- 98may20, cca
   -------------------------------------------------------------
*/
void
FrontMtx_MPI_permuteUpperAdj (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int          count, destination, incount, iproc, J, K, left, lasttag, 
             myid, ncolJ, nfront, nproc, offset, outcount, right,
             source, tag, tagbound ;
int          *colindJ, *inbuffer, *incounts, *mark, *owners, 
             *outbuffer, *outcounts, *par ;
IVL          *coladjIVL ;
MPI_Status   status ;
/*
   --------------
   check the data
   --------------
*/
if ( frontmtx == NULL || frontOwnersIV == NULL || stats == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(msgFile, "\n fatal error in FrontMtx_MPI_permuteUpperAdj()"
         "\n frontmtx %p, frontOwnersIV %p, firsttag %d"
         "\n stats %p, msglvl %d, msgFile %p"
         "\n bad input\n", 
         frontmtx, frontOwnersIV, firsttag, stats, msglvl, msgFile) ;
   exit(-1) ;
}
/*
   ----------------------------------------------
   get id of self, # of processes and # of fronts
   ----------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside FrontMtx_MPI_permuteUpperAdj"
           "\n nproc = %d, myid = %d", nproc, myid) ;
   fflush(msgFile) ;
}
/*
   -------------------
   check the tag value
   -------------------
*/
lasttag  = firsttag + nproc ;
tagbound = maxTagMPI(comm) ;
if ( firsttag < 0 || lasttag > tagbound ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_permuteUpperAdj()"
           "\n firsttag = %d, tagbound = %d", firsttag, tagbound) ;
   exit(-1) ;
}

nfront    = FrontMtx_nfront(frontmtx) ;
coladjIVL = frontmtx->coladjIVL ;
par       = frontmtx->frontETree->tree->par ;
owners    = IV_entries(frontOwnersIV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n coladjIVL") ;
   IVL_writeForHumanEye(coladjIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   step 1 : determine the message size that this
      process will send to each other process
   ---------------------------------------------
*/
incounts  = IVinit(2*nproc, 0) ;
outcounts = incounts + nproc ;
mark      = IVinit(nfront, -1) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( iproc != myid ) {
/*
      -----------------------------------
      set mark[J] = 1 if iproc supports J
      -----------------------------------
*/
      IVfill(nfront, mark, -1) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == iproc ) {
            for ( K = J ; K != -1 && mark[K] == -1 ; K = par[K] ) {
               mark[K] = 1 ;
            }
         }
      }
/*
      ------------------------------------------------
      compute the size of the message to send to iproc
      ------------------------------------------------
*/
      for ( J = count = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == myid && mark[J] == 1 ) {
            FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
            count += 2 + ncolJ ;
         }
      }
      outcounts[iproc] = count ;
   }
}
/*
   -------------------------------
   do an all-to-all gather/scatter
   -------------------------------
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
   -----------------------------
   set up the in and out buffers
   -----------------------------
*/
count = IVmax(nproc, incounts, &iproc) ;
inbuffer = IVinit(count, -1) ;
count = IVmax(nproc, outcounts, &iproc) ;
outbuffer = IVinit(count, -1) ;
/*
   ----------------------------------------
   step 2: loop over the other processes,
      gather information and send them off,
      receive information
   ----------------------------------------
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
/*
      -----------------------------------
      set mark[J] = 1 if right supports J
      -----------------------------------
*/
      IVfill(nfront, mark, -1) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == right ) {
            for ( K = J ; K != -1 && mark[K] == -1 ; K = par[K] ) {
               mark[K] = 1 ;
            }
         }
      }
/*
      ----------------------------------------------
      load the message with the owned column indices
      that are needed by process right
      ----------------------------------------------
*/
      for ( J = count = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == myid && mark[J] == 1 ) {
            FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
            if ( msglvl > 1 ) {
               fprintf(msgFile, "\n loading adj(%d) :", J) ;
               IVfprintf(msgFile, ncolJ, colindJ) ;
            }
            outbuffer[count++] = J ;
            outbuffer[count++] = ncolJ ;
            IVcopy(ncolJ, outbuffer + count, colindJ) ;
            count += ncolJ ;
         }
      }
      if ( count != outcount ) {
         fprintf(stderr, 
                 "\n 0. fatal error in FrontMtx_MPI_permuteUpperAdj()"
                 "\n proc %d : count = %d, outcount = %d\n",
                 myid, count, incount) ;
         exit(-1) ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n message to %d", right) ;
         IVfprintf(msgFile, count, outbuffer) ;
         fflush(msgFile) ;
      }
      destination = right ;
      stats[0]++ ;
      stats[2] += outcount ;
   } else {
      destination = MPI_PROC_NULL ;
   }
   if ( incount > 0 ) {
      source = left ;
      stats[1]++ ;
      stats[3] += incount ;
   } else {
      source = MPI_PROC_NULL ;
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
                 "\n 1. fatal error in FrontMtx_MPI_permuteUpperAdj()"
                 "\n proc %d : source = %d, count = %d, incount = %d\n",
                 myid, source, count, incount) ;
         exit(-1) ;
      }
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n\n message from %d", source) ;
         IVfprintf(msgFile, count, inbuffer) ;
         fflush(msgFile) ;
      }
   }
/*
   --------------------------------------------
   set up the new lists in the coladjIVL object
   --------------------------------------------
*/
   count = 0 ;
   while ( count < incount ) {
      J     = inbuffer[count++] ;
      ncolJ = inbuffer[count++] ;
      if ( msglvl > 1 ) {
         fprintf(msgFile, "\n setting list (%d) :", J) ;
         IVfprintf(msgFile, ncolJ, inbuffer + count) ;
      }
      IVL_setList(coladjIVL, J, ncolJ, inbuffer + count) ;
      count += ncolJ ;
   }
   if ( count != incount ) {
      fprintf(stderr, 
              "\n 2. fatal error in FrontMtx_MPI_permuteUpperAdj()"
              "\n proc %d : source = %d, count = %d, incount = %d\n",
              myid, source, count, incount) ;
      exit(-1) ;
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n necessary upper adj") ;
   IVL_writeForHumanEye(coladjIVL, msgFile) ;
}
/*
   ----------------------------------------------------------------
   now reorder the supported portion of the column adjacency object
   ----------------------------------------------------------------
*/
FrontMtx_permuteUpperAdj(frontmtx, msglvl, msgFile) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(mark) ;
IVfree(incounts) ;
IVfree(inbuffer) ;
IVfree(outbuffer) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to permute the indices of the lower adjacency
      structure so that for each front J, bnd{J} is in 
      ascending order w.r.t. K cup bnd{K} for par[J] = K.
      process q sends to process r one message that contains
      J cup bnd{J} for all J owned by q and needed by r.
      once all the indices for the supported fronts are present,
      the indices in the lower adjacency structure are reordered
      as necessary.

   created -- 98may20, cca
   -------------------------------------------------------------
*/
void
FrontMtx_MPI_permuteLowerAdj (
   FrontMtx   *frontmtx,
   IV         *frontOwnersIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
int          count, destination, incount, iproc, J, K, left, lasttag, 
             myid, nfront, nproc, nrowJ, offset, outcount, right,
             source, tag, tagbound ;
int          *inbuffer, *incounts, *mark, *owners, 
             *outbuffer, *outcounts, *par, *rowindJ ;
IVL          *rowadjIVL ;
MPI_Status   status ;
/*
   --------------
   check the data
   --------------
*/
if ( frontmtx == NULL || frontOwnersIV == NULL || stats == NULL
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(msgFile, "\n fatal error in FrontMtx_MPI_permuteLowerAdj()"
         "\n frontmtx %p, frontOwnersIV %p, firsttag %d"
         "\n stats %p, msglvl %d, msgFile %p"
         "\n bad input\n", 
         frontmtx, frontOwnersIV, firsttag, stats, msglvl, msgFile) ;
   exit(-1) ;
}
/*
   ----------------------------------------------
   get id of self, # of processes and # of fronts
   ----------------------------------------------
*/
MPI_Comm_rank(comm, &myid) ;
MPI_Comm_size(comm, &nproc) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n inside FrontMtx_MPI_permuteLowerAdj"
           "\n nproc = %d, myid = %d", nproc, myid) ;
   fflush(msgFile) ;
}
/*
   -------------------
   check the tag value
   -------------------
*/
lasttag  = firsttag + nproc ;
tagbound = maxTagMPI(comm) ;
if ( firsttag < 0 || lasttag > tagbound ) {
   fprintf(stderr, "\n fatal error in FrontMtx_MPI_permuteUpperAdj()"
           "\n firsttag = %d, tagbound = %d", firsttag, tagbound) ;
   exit(-1) ;
}

nfront    = FrontMtx_nfront(frontmtx) ;
rowadjIVL = frontmtx->rowadjIVL ;
par       = frontmtx->frontETree->tree->par ;
owners    = IV_entries(frontOwnersIV) ;
/*
   ---------------------------------------------
   step 1 : determine the message size that this
      process will send to each other process
   ---------------------------------------------
*/
incounts  = IVinit(2*nproc, 0) ;
outcounts = incounts + nproc ;
mark      = IVinit(nfront, -1) ;
for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
   if ( iproc != myid ) {
/*
      -----------------------------------
      set mark[J] = 1 if iproc supports J
      -----------------------------------
*/
      IVfill(nfront, mark, -1) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == iproc ) {
            for ( K = J ; K != -1 && mark[K] == -1 ; K = par[K] ) {
               mark[K] = 1 ;
            }
         }
      }
/*
      ------------------------------------------------
      compute the size of the message to send to iproc
      ------------------------------------------------
*/
      for ( J = count = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == myid && mark[J] == 1 ) {
            FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
            count += 2 + nrowJ ;
         }
      }
      outcounts[iproc] = count ;
   }
}
/*
   -------------------------------
   do an all-to-all gather/scatter
   -------------------------------
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
   -----------------------------
   set up the in and out buffers
   -----------------------------
*/
count = IVmax(nproc, incounts, &iproc) ;
inbuffer = IVinit(count, -1) ;
count = IVmax(nproc, outcounts, &iproc) ;
outbuffer = IVinit(count, -1) ;
/*
   ----------------------------------------
   step 2: loop over the other processes,
      gather information and send them off,
      receive information
   ----------------------------------------
*/
tag = firsttag ;
for ( offset = 1 ; offset < nproc ; offset++ ) {
   right = (myid + offset) % nproc ;
   if ( offset <= myid ) {
      left = myid - offset ;
   } else {
      left = nproc + myid - offset ;
   }
   outcount = outcounts[right] ;
   incount  = incounts[left] ;
   if ( msglvl > 1 ) {
      fprintf(msgFile,
         "\n ### process %d, send %d to right %d, recv %d from left %d",
              myid, outcount, right, incount, left) ;
      fflush(msgFile) ;
   }
   if ( outcount > 0 ) {
/*
      -----------------------------------
      set mark[J] = 1 if right supports J
      -----------------------------------
*/
      IVfill(nfront, mark, -1) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == right ) {
            for ( K = J ; K != -1 && mark[K] == -1 ; K = par[K] ) {
               mark[K] = 1 ;
            }
         }
      }
/*
      -------------------------------------------
      load the message with the owned row indices
      that are needed by process right
      -------------------------------------------
*/
      for ( J = count = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == myid && mark[J] == 1 ) {
            FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
            outbuffer[count++] = J ;
            outbuffer[count++] = nrowJ ;
            IVcopy(nrowJ, outbuffer + count, rowindJ) ;
            count += nrowJ ;
         }
      }
      destination = right ;
      stats[0]++ ;
      stats[2] += outcount ;
   } else {
      destination = MPI_PROC_NULL ;
   }
   if ( incount > 0 ) {
      source = left ;
      stats[1]++ ;
      stats[3] += incount ;
   } else {
      source = MPI_PROC_NULL ;
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
                 "\n 1. fatal error in FrontMtx_MPI_permuteLowerAdj()"
                 "\n proc %d : source = %d, count = %d, incount = %d\n",
                 myid, source, count, incount) ;
         exit(-1) ;
      }
   }
/*
   --------------------------------------------
   set up the new lists in the rowadjIVL object
   --------------------------------------------
*/
   count = 0 ;
   while ( count < incount ) {
      J     = inbuffer[count++] ;
      nrowJ = inbuffer[count++] ;
      IVL_setList(rowadjIVL, J, nrowJ, inbuffer + count) ;
      count += nrowJ ;
   }
   if ( count != incount ) {
      fprintf(stderr, 
              "\n 2. fatal error in FrontMtx_MPI_permuteLowerAdj()"
              "\n proc %d : source = %d, count = %d, incount = %d\n",
              myid, source, count, incount) ;
      exit(-1) ;
   }
}
/*
   -------------------------------------------------------------
   now reorder the supported portion of the row adjacency object
   -------------------------------------------------------------
*/
FrontMtx_permuteLowerAdj(frontmtx, msglvl, msgFile) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(mark) ;
IVfree(incounts) ;
IVfree(inbuffer) ;
IVfree(outbuffer) ;

return ; }

/*--------------------------------------------------------------------*/
