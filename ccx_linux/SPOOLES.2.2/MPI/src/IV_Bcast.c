/*  IV_Bcast.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to broadcast an IV object from 
              one processor to all the others.

   created -- 98sep25, cca
   -------------------------------------------
*/
IV *
IV_MPI_Bcast (
   IV         *obj, 
   int        root,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( obj == NULL ) {
   fprintf(stderr, "\n fatal error in IV_MPI_Bcast()"
           "\n obj is NULL\n") ;
   exit(-1) ;
}
MPI_Comm_rank(comm, &myid) ;
if ( myid == root ) {
/*
   --------------------------------
   broadcast the size of the vector
   --------------------------------
*/
   MPI_Bcast((void *) &obj->size, 1, MPI_INT, root, comm) ;
   if ( obj->size > 0 ) {
/*
      ---------------------
      broadcast the entries
      ---------------------
*/
      MPI_Bcast((void *) obj->vec, obj->size, MPI_INT, root, comm) ;
   }
} else {
   int   size ;
/*
   --------------
   clear the data
   --------------
*/
   if ( obj == NULL ) {
      obj = IV_new() ;
   } else {
      IV_clearData(obj) ;
   }
/*
   ------------------------------
   receive the size of the vector
   ------------------------------
*/
   MPI_Bcast((void *) &size, 1, MPI_INT, root, comm) ;
   IV_setSize(obj, size) ;
   if ( size > 0 ) {
/*
      -------------------
      receive the entries
      -------------------
*/
      MPI_Bcast((void *) obj->vec, obj->size, MPI_INT, root, comm) ;
   }
}
return(obj) ; }

/*--------------------------------------------------------------------*/
