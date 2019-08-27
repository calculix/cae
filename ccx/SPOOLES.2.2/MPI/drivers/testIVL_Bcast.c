/*  testIVL_Bcast.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------------------------
   this program tests the IVL_MPI_Bcast() method

   (1) process root generates a random IVL object
       and computes its checksum
   (2) process root broadcasts the IVL object to the other processors
   (3) each process computes the checksum of its IVL object
   (4) the checksums are compared on root

   created -- 98sep10, cca
   ------------------------------------------------------------------
*/
{
char         *buffer ;
double       chksum, t1, t2 ;
double       *sums ;
Drand        drand ;
int          ilist, iproc, length, loc, maxlistsize, msglvl, myid, 
             nlist, nproc, root, seed, size ;
int          *list, *vec ;
IVL          *ivl ;
FILE         *msgFile ;
/*
   ---------------------------------------------------------------
   find out the identity of this process and the number of process
   ---------------------------------------------------------------
*/
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
fprintf(stdout, "\n process %d of %d, argc = %d", myid, nproc, argc) ;
fflush(stdout) ;
if ( argc != 7 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile nlist maxlistsize seed "
           "\n    msglvl      -- message level"
           "\n    msgFile     -- message file"
           "\n    nlist       -- number of lists in the IVL object"
           "\n    maxlistsize -- maximum size of a list"
           "\n    root        -- root processor for broadcast"
           "\n    seed        -- random number seed"
           "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else {
   length = strlen(argv[2]) + 1 + 4 ;
   buffer = CVinit(length, '\0') ;
   sprintf(buffer, "%s.%d", argv[2], myid) ;
   if ( (msgFile = fopen(buffer, "w")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], argv[2]) ;
      return(-1) ;
   }
   CVfree(buffer) ;
}
nlist       = atoi(argv[3]) ;
maxlistsize = atoi(argv[4]) ;
root        = atoi(argv[5]) ;
seed        = atoi(argv[6]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n nlist       -- %d" 
        "\n maxlistsize -- %d" 
        "\n root        -- %d" 
        "\n seed        -- %d" 
        "\n",
        argv[0], msglvl, argv[2], nlist, maxlistsize, root, seed) ;
fflush(msgFile) ;
/*
   --------------------------------------------
   set up the IVL object and fill owned entries
   --------------------------------------------
*/
MARKTIME(t1) ;
ivl = IVL_new() ;
if ( myid == root ) {
   IVL_init1(ivl, IVL_CHUNKED, nlist) ;
   vec = IVinit(maxlistsize, -1) ;
   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, maxlistsize) ;
   for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
      size = (int) Drand_value(&drand) ;
      Drand_fillIvector(&drand, size, vec) ;
      IVL_setList(ivl, ilist, size, vec) ;
   }
   IVfree(vec) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the IVL object", t2 - t1) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(ivl, msgFile) ;
} else {
   IVL_writeStats(ivl, msgFile) ;
}
fflush(msgFile) ;
if ( myid == root ) {
/*
   --------------------------------------
   compute the checksum of the ivl object
   --------------------------------------
*/
   for ( ilist = 0, chksum = 0.0 ; ilist < nlist ; ilist++ ) {
         IVL_listAndSize(ivl, ilist, &size, &list) ;
         chksum += 1 + ilist + size + IVsum(size, list) ;
   }
   fprintf(msgFile, "\n\n local chksum = %12.4e", chksum) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   broadcast the IVL object
   ------------------------
*/
MARKTIME(t1) ;
ivl = IVL_MPI_Bcast(ivl, root, msglvl, msgFile, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : broadcast the IVL object", t2 - t1) ;
/*
   --------------------------------------
   compute the checksum of the ivl object
   --------------------------------------
*/
for ( ilist = 0, chksum = 0.0 ; ilist < nlist ; ilist++ ) {
      IVL_listAndSize(ivl, ilist, &size, &list) ;
      chksum += 1 + ilist + size + IVsum(size, list) ;
}
fprintf(msgFile, "\n\n local chksum = %12.4e", chksum) ;
fflush(msgFile) ;
/*
   ---------------------------------------
   gather the checksums from the processes
   ---------------------------------------
*/
sums = DVinit(nproc, 0.0) ;
MPI_Gather((void *) &chksum, 1, MPI_DOUBLE, 
           (void *) sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n\n sums") ;
   DVfprintf(msgFile, nproc, sums) ;
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      sums[iproc] -= chksum ;
   }
   fprintf(msgFile, "\n\n errors") ;
   DVfprintf(msgFile, nproc, sums) ;
   fprintf(msgFile, "\n\n maxerror = %12.4e", DVmax(nproc, sums, &loc));
}
/*
   ----------------
   free the objects
   ----------------
*/
DVfree(sums) ;
IVL_free(ivl) ;
/*
   ------------------------
   exit the MPI environment
   ------------------------
*/
MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
