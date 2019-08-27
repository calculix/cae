/*  ETreeBcast.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to broadcast a front tree object 
              from one process to all the others

   created -- 98may21, cca
   ---------------------------------------------
*/
ETree *
ETree_MPI_Bcast (
   ETree      *etree,
   int        root,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
int   myid, nvtx, nfront, nint ;
int   *buffer ;
/*
   -------------
   find identity
   -------------
*/
MPI_Comm_rank(comm, &myid) ;
if ( myid == root ) {
/*
   --------------------------------------------
   this process owns the front tree, allocate a
   continuous buffer and load the data into it.
   --------------------------------------------
*/
   nfront = ETree_nfront(etree) ;
   nvtx   = ETree_nvtx(etree) ;
   nint   = 3 + 5*nfront + nvtx ;
   buffer = IVinit(nint, -1) ;
   buffer[0] = nfront ;
   buffer[1] = nvtx  ;
   buffer[2] = ETree_root(etree) ;
   IVcopy(nfront, buffer + 3,            ETree_par(etree)) ;
   IVcopy(nfront, buffer + 3 +   nfront, ETree_fch(etree)) ;
   IVcopy(nfront, buffer + 3 + 2*nfront, ETree_sib(etree)) ;
   IVcopy(nfront, buffer + 3 + 3*nfront, ETree_nodwghts(etree)) ;
   IVcopy(nfront, buffer + 3 + 4*nfront, ETree_bndwghts(etree)) ;
   IVcopy(nvtx,   buffer + 3 + 5*nfront, ETree_vtxToFront(etree)) ;
/*
   ------------------------------------
   send the size of the buffer and then 
   the buffer to the other processors
   ------------------------------------
*/
   MPI_Bcast(&nint,     1, MPI_INT, root, comm) ;
   MPI_Bcast(buffer, nint, MPI_INT, root, comm) ;
} else {
/*
   --------------------------------------------
   this process will receive the front tree.
   clear its data, receive the number of int's,
   then receive the buffer
   --------------------------------------------
*/
   if ( etree != NULL ) {
      ETree_free(etree) ;
   }
   MPI_Bcast(&nint, 1, MPI_INT, root, comm) ;
   buffer = IVinit(nint, -1) ;
   MPI_Bcast(buffer, nint, MPI_INT, root, comm) ;
/*
   ----------------------------------------
   create an ETree object and fill its data
   ----------------------------------------
*/
   etree  = ETree_new() ;
   nfront = buffer[0] ;
   nvtx   = buffer[1] ;
   ETree_init1(etree, nfront, nvtx) ;
   etree->tree->n    = nfront ;
   etree->tree->root = buffer[2] ;
   IVcopy(nfront, ETree_par(etree),        buffer + 3) ;
   IVcopy(nfront, ETree_fch(etree),        buffer + 3 +   nfront) ;
   IVcopy(nfront, ETree_sib(etree),        buffer + 3 + 2*nfront) ;
   IVcopy(nfront, ETree_nodwghts(etree),   buffer + 3 + 3*nfront) ;
   IVcopy(nfront, ETree_bndwghts(etree),   buffer + 3 + 4*nfront) ;
   IVcopy(nvtx,   ETree_vtxToFront(etree), buffer + 3 + 5*nfront) ;
}
/*
   ---------------
   free the buffer
   ---------------
*/
IVfree(buffer) ;

return(etree) ; }

/*--------------------------------------------------------------------*/
