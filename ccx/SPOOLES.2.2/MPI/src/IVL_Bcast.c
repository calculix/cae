/*  IVL_Bcast.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   purpose -- to broadcast an IVL object from
              one processor to all the others

   created -- 98sep10, cca
   ------------------------------------------
*/
IVL *
IVL_MPI_Bcast (
   IVL        *ivl,
   int        root,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
int   ilist, myid, nlist ;
/*
   -------------
   find identity
   -------------
*/
MPI_Comm_rank(comm, &myid) ;
if ( myid == root ) {
/*
   -----------------------------------------------
   this process owns the front tree,
   broadcast the number of lists in the front tree
   -----------------------------------------------
*/
   nlist = ivl->nlist ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n broadcasting %d ", nlist) ;
      fflush(msgFile) ;
   }
   MPI_Bcast(&nlist, 1, MPI_INT, root, comm) ;
/*
   ------------------------
   broadcast the list sizes
   ------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n broadcasting sizes[]") ;
      fflush(msgFile) ;
   }
   MPI_Bcast(ivl->sizes, nlist, MPI_INT, root, comm) ;
/*
   -------------------------------------------
   loop over the lists and broadcast each list
   -------------------------------------------
*/
   for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n broadcasting list %d", ilist) ;
         fflush(msgFile) ;
      }
      MPI_Bcast(ivl->p_vec[ilist], ivl->sizes[ilist], MPI_INT, 
                root, comm) ;
   }
} else {
   int   *sizes ;
/*
   --------------
   clear the data
   --------------
*/
   if ( ivl == NULL ) {
      ivl = IVL_new() ;
   } else {
      IVL_clearData(ivl) ;
   }
/*
   ---------------------------------------------
   receive the number of lists in the front tree
   ---------------------------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n receiving nlist") ;
      fflush(msgFile) ;
   }
   MPI_Bcast(&nlist, 1, MPI_INT, root, comm) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, " %d", nlist) ;
      fflush(msgFile) ;
   }
/*
   ---------------------
   receive the list sizes
   ---------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n receiving sizes") ;
      fflush(msgFile) ;
   }
   sizes = IVinit(nlist, 0) ;
   MPI_Bcast(sizes, nlist, MPI_INT, root, comm) ;
   if ( msglvl > 2 ) {
      IVfprintf(msgFile, nlist, sizes) ;
      fflush(msgFile) ;
   }
/*
   -------------------------
   initialize the IVL object
   -------------------------
*/
   IVL_init3(ivl, IVL_CHUNKED, nlist, sizes) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n initialized IVL object") ;
      IVL_writeForHumanEye(ivl, msgFile) ;
      fflush(msgFile) ;
   }
   IVfree(sizes) ;
/*
   -----------------------------------------
   loop over the lists and receive each list
   -----------------------------------------
*/
   for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n receiving list %d, size %d, loc %p", 
                 ilist, ivl->sizes[ilist], ivl->p_vec[ilist]) ;
         fflush(msgFile) ;
      }
      MPI_Bcast(ivl->p_vec[ilist], ivl->sizes[ilist], MPI_INT, 
                root, comm) ;
      if ( msglvl > 2 ) {
         IVfprintf(msgFile, ivl->sizes[ilist], ivl->p_vec[ilist]) ;
         fflush(msgFile) ;
      }
   }
}
return(ivl) ; }

/*--------------------------------------------------------------------*/
