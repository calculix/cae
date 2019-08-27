/*  testIVL_alltoall.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   this program tests the IVL_MPI_alltoall() method.

   (1) each process creates a "send" IVL object 
       and fills the lists with random numbers
   (3) the processes all-to-all gather the lists 
       to create a "recv" IVL object

   created -- 98jul26, cca
   -------------------------------------------------
*/
{
char         *buffer ;
double       chksum, globalsum1, globalsum2, t1, t2 ;
Drand        drand ;
int          ilist, length, myid, msglvl, n, nlist, 
             nproc, rc, seed, size, tag ;
int          *list, *vec ;
int          stats[4], tstats[4] ;
IVL          *recvIVL, *sendIVL ;
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
if ( argc != 5 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile n seed "
           "\n    msglvl  -- message level"
           "\n    msgFile -- message file"
           "\n    n       -- maximum size of any list"
           "\n    seed    -- random number seed"
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
n     = atoi(argv[3]) ;
seed  = atoi(argv[4]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl  -- %d" 
        "\n msgFile -- %s" 
        "\n n       -- %d" 
        "\n seed    -- %d" 
        "\n",
        argv[0], msglvl, argv[2], n, seed) ;
fflush(msgFile) ;
/*
   ----------------------------
   set up the "send" IVL object
   ----------------------------
*/
MARKTIME(t1) ;
sendIVL = IVL_new() ;
nlist = nproc ;
IVL_init1(sendIVL, IVL_CHUNKED, nlist) ;
vec = IVinit(n, -1) ;
Drand_setDefaultFields(&drand) ;
Drand_setSeed(&drand, seed + myid) ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   Drand_setUniform(&drand, 0, n) ;
   size = (int) Drand_value(&drand) ;
   Drand_setUniform(&drand, 0, n*n) ;
   Drand_fillIvector(&drand, size, vec) ;
   IVqsortUp(size, vec) ;
   IVL_setList(sendIVL, ilist, size, vec) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the IVL object",
        t2 - t1) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(sendIVL, msgFile) ;
} else {
   IVL_writeStats(sendIVL, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------------
   compute the local checksum of the ivl object
   --------------------------------------------
*/
for ( ilist = 0, chksum = 0.0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(sendIVL, ilist, &size, &list) ;
   chksum += 1 + ilist + size + IVsum(size, list) ;
}
fprintf(msgFile, "\n\n local partial chksum = %12.4e", chksum) ;
fflush(msgFile) ;
/*
   -----------------------------
   get the first global checksum
   -----------------------------
*/
rc = MPI_Allreduce((void *) &chksum, (void *) &globalsum1, 
                   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
/*
   --------------------------------
   execute the all-to-all operation
   --------------------------------
*/
tag = 47 ;
IVzero(4, stats) ;
recvIVL = IVL_MPI_alltoall(sendIVL, NULL, stats, 
                           msglvl, msgFile, tag, MPI_COMM_WORLD) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n return from IVL_MPI_alltoall()") ;
   fprintf(msgFile, 
           "\n local send stats : %10d messages with %10d bytes"
           "\n local recv stats : %10d messages with %10d bytes",
           stats[0], stats[2], stats[1], stats[3]) ;
   fflush(msgFile) ;
}
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
           "\n total send stats : %10d messages with %10d bytes"
           "\n total recv stats : %10d messages with %10d bytes",
           tstats[0], tstats[2], tstats[1], tstats[3]) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n recvIVL") ;
   IVL_writeForHumanEye(recvIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   compute the second global checksum
   ----------------------------------
*/
for ( ilist = 0, chksum = 0.0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(recvIVL, ilist, &size, &list) ;
   chksum += 1 + ilist + size + IVsum(size, list) ;
}
rc = MPI_Allreduce((void *) &chksum, (void *) &globalsum2, 
                   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
fprintf(msgFile, 
        "\n globalsum1 = %12.4e, globalsum2 = %12.4e, error = %12.4e",
        globalsum1, globalsum2, fabs(globalsum1 - globalsum2)) ;
fflush(msgFile) ;
/*
   ----------------
   free the objects
   ----------------
*/
IVL_free(sendIVL) ;
IVL_free(recvIVL) ;
IVfree(vec) ;
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
