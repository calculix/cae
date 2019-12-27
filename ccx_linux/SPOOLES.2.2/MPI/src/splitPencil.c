/*  setupPencil.c  */

#include "../spoolesMPI.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- to split a distributed Pencil object

   pencil     -- pointer to the local Pencil object
   mapIV      -- pointer to the map from vertices to processes
   firsttag   -- first tag value, two will be used, tag and tag+1
   stats[4]    -- statistics vector
      stats[0] -- # of messages sent
      stats[1] -- # of messages received
      stats[2] -- # of bytes sent
      stats[3] -- # of bytes received
   msglvl     -- local message level
   msgFile    -- local message file
   comm       -- MPI communication structure
 
   created  -- 98may20, cca
   --------------------------------------------------------------
*/
void
Pencil_MPI_split (
   Pencil     *pencil,
   IV         *mapIV,
   int        stats[],
   int        msglvl,
   FILE       *msgFile,
   int        firsttag,
   MPI_Comm   comm
) {
InpMtx   *inpmtxA, *inpmtxB, *keepmtx ;
int      tag ;
/*
   -----------------------
   check the range of tags
   -----------------------
*/
if ( firsttag < 0 || firsttag + 1 > maxTagMPI(comm) ) {
   fprintf(stderr, "\n fatal error in Pencil_MPI_split()"
           "\n range of tags is [%d,%d], tag_bound = %d",
           firsttag, firsttag + 1, maxTagMPI(comm)) ;
   exit(-1) ;
}
/*
   ------------------------------------
   split the DInpMtx object into pieces
   ------------------------------------
*/
tag = firsttag ;
if ( (inpmtxA = pencil->inpmtxA) != NULL ) {
   keepmtx = InpMtx_MPI_split(inpmtxA, mapIV, stats,
                              msglvl, msgFile, tag, comm) ;
   InpMtx_free(inpmtxA) ;
   pencil->inpmtxA = keepmtx ;
}
tag += 1 ;
if ( (inpmtxB = pencil->inpmtxB) != NULL ) {
   keepmtx = InpMtx_MPI_split(inpmtxB, mapIV, stats,
                              msglvl, msgFile, tag, comm) ;
   InpMtx_free(inpmtxB) ;
   pencil->inpmtxB = keepmtx ;
}
tag += 1 ;

return ; }

/*--------------------------------------------------------------------*/
