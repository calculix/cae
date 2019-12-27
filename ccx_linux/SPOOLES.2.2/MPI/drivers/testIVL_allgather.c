/*  testIVLallgather.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   this program tests the IVL_MPI_allgather() method

   (1) each process generates the same owners[n] map
   (2) each process creates an IVL object 
       and fills its owned lists with random numbers
   (3) the processes gather-all's the lists of ivl

   created -- 98apr03, cca
   -------------------------------------------------
*/
{
char         *buffer ;
double       chksum, globalsum, t1, t2 ;
Drand        drand ;
int          ilist, length, myid, msglvl, nlist, 
             nproc, rc, seed, size, tag ;
int          *list, *owners, *vec ;
int          stats[4], tstats[4] ;
IV           *ownersIV ;
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
if ( argc != 5 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile n seed "
           "\n    msglvl  -- message level"
           "\n    msgFile -- message file"
           "\n    nlist   -- number of lists in the IVL object"
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
nlist = atoi(argv[3]) ;
seed  = atoi(argv[4]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl  -- %d" 
        "\n msgFile -- %s" 
        "\n nlist   -- %d" 
        "\n seed    -- %d" 
        "\n",
        argv[0], msglvl, argv[2], nlist, seed) ;
fflush(msgFile) ;
/*
   ----------------------------
   generate the ownersIV object
   ----------------------------
*/
MARKTIME(t1) ;
ownersIV = IV_new() ;
IV_init(ownersIV, nlist, NULL) ;
owners = IV_entries(ownersIV) ;
Drand_setDefaultFields(&drand) ;
Drand_setSeed(&drand, seed) ;
Drand_setUniform(&drand, 0, nproc) ;
Drand_fillIvector(&drand, nlist, owners) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the ownersIV object",
        t2 - t1) ;
fflush(msgFile) ;
fprintf(msgFile, "\n\n ownersIV generated") ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(ownersIV, msgFile) ;
} else {
   IV_writeStats(ownersIV, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------------
   set up the IVL object and fill owned entries
   --------------------------------------------
*/
MARKTIME(t1) ;
ivl = IVL_new() ;
IVL_init1(ivl, IVL_CHUNKED, nlist) ;
vec = IVinit(nlist, -1) ;
Drand_setSeed(&drand, seed + myid) ;
Drand_setUniform(&drand, 0, nlist) ;
for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
   if ( owners[ilist] == myid ) {
      size = (int) Drand_value(&drand) ;
      Drand_fillIvector(&drand, size, vec) ;
      IVL_setList(ivl, ilist, size, vec) ;
   }
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the IVL object",
        t2 - t1) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(ivl, msgFile) ;
} else {
   IVL_writeStats(ivl, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------------
   compute the local checksum of the ivl object
   --------------------------------------------
*/
for ( ilist = 0, chksum = 0.0 ; ilist < nlist ; ilist++ ) {
   if ( owners[ilist] == myid ) {
      IVL_listAndSize(ivl, ilist, &size, &list) ;
      chksum += 1 + ilist + size + IVsum(size, list) ;
   }
}
fprintf(msgFile, "\n\n local partial chksum = %12.4e", chksum) ;
fflush(msgFile) ;
/*
   -----------------------
   get the global checksum
   -----------------------
*/
rc = MPI_Allreduce((void *) &chksum, (void *) &globalsum, 
                   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
/*
   --------------------------------
   execute the all-gather operation
   --------------------------------
*/
tag = 47 ;
IVzero(4, stats) ;
IVL_MPI_allgather(ivl, ownersIV, 
                  stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n return from IVL_MPI_allgather()") ;
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
   fprintf(msgFile, "\n\n ivl") ;
   IVL_writeForHumanEye(ivl, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------
   compute the checksum of the entire object
   -----------------------------------------
*/
for ( ilist = 0, chksum = 0.0 ; ilist < nlist ; ilist++ ) {
   IVL_listAndSize(ivl, ilist, &size, &list) ;
   chksum += 1 + ilist + size + IVsum(size, list) ;
}
fprintf(msgFile, 
        "\n globalsum = %12.4e, chksum = %12.4e, error = %12.4e",
        globalsum, chksum, fabs(globalsum - chksum)) ;
fflush(msgFile) ;
/*
   ----------------
   free the objects
   ----------------
*/
IV_free(ownersIV) ;
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
