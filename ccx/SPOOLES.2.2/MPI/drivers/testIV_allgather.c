/*  testIVallgather.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------
   this program tests the IV_MPI_allgather() method

   (1) each process generates the same owners[n] map
   (2) each process creates a vec[n] vector
       and fills its owned entries with random numbers
   (3) the processes gather-all the entries of vec[]

   created -- 98may20, cca
   -------------------------------------------------------
*/
{
char         *buffer ;
double       chksum, globalsum, t1, t2 ;
Drand        drand ;
int          ii, length, myid, msglvl, n, nproc, rc, seed, tag ;
int          *owners, *vec ;
int          stats[4] ;
IV           *objIV, *ownersIV ;
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
           "\n    n       -- number of entries in the vector"
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
n    = atoi(argv[3]) ;
seed = atoi(argv[4]) ;
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
   generate the ownersIV object
   ----------------------------
*/
MARKTIME(t1) ;
ownersIV = IV_new() ;
IV_init(ownersIV, n, NULL) ;
owners = IV_entries(ownersIV) ;
Drand_setDefaultFields(&drand) ;
Drand_setSeed(&drand, seed) ;
Drand_setUniform(&drand, 0, nproc) ;
Drand_fillIvector(&drand, n, owners) ;
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
   ----------------------------------------------
   set up the vec[] vector and fill owned entries
   ----------------------------------------------
*/
MARKTIME(t1) ;
objIV = IV_new() ;
IV_init(objIV, n, NULL) ;
vec = IV_entries(objIV) ;
IVfill(n, vec, -1) ;
Drand_setSeed(&drand, seed + myid) ;
Drand_setUniform(&drand, 0, n) ;
for ( ii = 0 ; ii < n ; ii++ ) {
   if ( owners[ii] == myid ) {
      vec[ii] = (int) Drand_value(&drand) ;
   }
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the objIV object",
        t2 - t1) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(objIV, msgFile) ;
} else {
   IV_writeStats(objIV, msgFile) ;
}
fflush(msgFile) ;
/*
   ----------------------------------------
   compute the local checksum of the vector
   ----------------------------------------
*/
for ( ii = 0, chksum = 0.0 ; ii < n ; ii++ ) {
   if ( owners[ii] == myid ) {
      chksum += (ii + 1)*(vec[ii] + 1) ;
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
IV_MPI_allgather(objIV, ownersIV, stats, 
                 msglvl, msgFile, tag, MPI_COMM_WORLD) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n return from IV_MPI_allgather()") ;
   fprintf(msgFile, "\n\n stats") ;
   IVfprintf(msgFile, 4, stats) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n objIV") ;
   IV_writeForHumanEye(objIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------
   compute the checksum of the entire vector
   -----------------------------------------
*/
for ( ii = 0, chksum = 0.0 ; ii < n ; ii++ ) {
   chksum += (ii + 1)*(vec[ii] + 1) ;
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
IV_free(objIV) ;
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
