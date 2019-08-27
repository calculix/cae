/*  testScatterDenseMtx.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------------------------
   this program tests the simple split and merge of a DenseMtx object

   (1) process root generates a nrow x ncol DenseMtx object
   (2) using a random map, scatter the rows into local matrices
   (3) gather the local matrices into a second global matrix
       and compare

   created -- 98sep26, cca
   ------------------------------------------------------------------
*/
{
char        *buffer ;
DenseMtx    *X, *Y, *Z ;
double      t1, t2 ;
double      checksums[3], Xchecksums[3], Ychecksums[3], Zchecksums[3] ;
Drand       drand ;
int         inc1, inc2, length, myid, msglvl, ncol, 
            nproc, nrow, root, seed, tag, type, v ;
int         *map ;
int         stats[4], tstats[4] ;
IV          *mapIV ;
FILE        *msgFile ;
/*
   ---------------------------------------------------------------
   find out the identity of this process and the number of process
   ---------------------------------------------------------------
*/
/*
fprintf(stdout, "\n stdout, MPI_COMM_WORLD = %p", MPI_COMM_WORLD) ;
fflush(stdout) ;
*/
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
/*
fprintf(stdout, "\n process %d of %d, argc = %d", myid, nproc, argc) ;
fflush(stdout) ;
*/
if ( argc != 10 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile type nrow ncol inc1 inc2 seed root"
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    type    -- type of entries"
"\n      1 -- real"
"\n      2 -- complex"
"\n    nrow -- number of rows"
"\n    ncol -- number of columns"
"\n    inc1 -- row increment"
"\n    inc2 -- column increment"
"\n    seed -- random number seed"
"\n    root -- root processor"
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
      return(0) ;
   }
   CVfree(buffer) ;
}
type = atoi(argv[3]) ;
nrow = atoi(argv[4]) ;
ncol = atoi(argv[5]) ;
inc1 = atoi(argv[6]) ;
inc2 = atoi(argv[7]) ;
seed = atoi(argv[8]) ;
root = atoi(argv[9]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl     -- %d" 
        "\n msgFile    -- %s" 
        "\n type       -- %d" 
        "\n nrow       -- %d" 
        "\n ncol       -- %d" 
        "\n inc1       -- %d" 
        "\n inc2       -- %d" 
        "\n seed       -- %d" 
        "\n root       -- %d" 
        "\n",
        argv[0], msglvl, argv[2], type, nrow, ncol, 
        inc1, inc2, seed, root) ;
fflush(msgFile) ;
/*
   -----------------------------
   generate the DenseMtx object
   -----------------------------
*/
if ( myid == root ) {
   X = DenseMtx_new() ;
   fprintf(msgFile, "\n X = %p", X) ;
   fflush(msgFile) ;
   MARKTIME(t1) ;
   DenseMtx_init(X, type, 1, -1, nrow, ncol, inc1, inc2) ;
   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0.0, 1.0) ;
   switch ( type ) {
   case SPOOLES_REAL :
      Drand_fillDvector(&drand, nrow*ncol, DenseMtx_entries(X)) ;
      break ;
   case SPOOLES_COMPLEX :
      Drand_fillDvector(&drand, 2*nrow*ncol, DenseMtx_entries(X)) ;
      break ;
   }
/*
   ------------------------------------------
   compute the serial checksums of the matrix
   ------------------------------------------
*/
   DenseMtx_checksums(X, Xchecksums) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : initialize the DenseMtx object",
           t2 - t1) ;
   fflush(msgFile) ;
   fprintf(msgFile, "\n\n DenseMtx generated") ;
   if ( msglvl > 2 ) {
      DenseMtx_writeForHumanEye(X, msgFile) ;
   } else {
      DenseMtx_writeStats(X, msgFile) ;
   }
   fflush(msgFile) ;
/*
   ---------------------
   generate a random map
   ---------------------
*/
   MARKTIME(t1) ;
   mapIV = IV_new() ;
   IV_init(mapIV, nrow, NULL) ;
   map = IV_entries(mapIV) ;
   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0.0, (double) nrow) ;
   for ( v = 0 ; v < nrow ; v++ ) {
      map[v] = ((int) Drand_value(&drand)) % nproc ;
      if ( map[v] == root ) {
         map[v] = nproc - 1 ;
      }
   }
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : generate random map", t2 - t1) ;
   fflush(msgFile) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n map") ;
      IV_writeForHumanEye(mapIV, msgFile) ;
      fflush(msgFile) ;
   }
} else {
   X = NULL ;
   mapIV = NULL ;
}
/*
   ------------------------------------------
   now scatter the matrix with the random map
   ------------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
tag = 1 ;
Y = DenseMtx_MPI_splitFromGlobalByRows(X, NULL, mapIV, root, stats, 
                              msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : scatter matrix ", t2 - t1) ;
fprintf(msgFile, 
        "\n send stats : %d messages with %d bytes"
        "\n recv stats : %d messages with %d bytes",
        stats[0], stats[2], stats[1], stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, root, MPI_COMM_WORLD) ;
if ( myid == root ) {
   fprintf(msgFile, 
           "\n total send stats : %d messages with %d bytes"
           "\n total recv stats : %d messages with %d bytes",
           tstats[0], tstats[2], tstats[1], tstats[3]) ;
   fflush(msgFile) ;
}
if ( Y != NULL ) {
   fprintf(msgFile, "\n\n local matrix") ;
   if ( msglvl > 2 ) {
      DenseMtx_writeForHumanEye(Y, msgFile) ;
   } else {
      DenseMtx_writeStats(Y, msgFile) ;
   }
} else {
   fprintf(msgFile, "\n\n local matrix is NULL") ;
}
fflush(msgFile) ;
/*
   -----------------------------------------------
   compute the checksums of the distributed matrix
   -----------------------------------------------
*/
if ( Y != NULL ) {
   DenseMtx_checksums(Y, checksums) ;
} else {
   checksums[0] = checksums[1] = checksums[2] = 0.0 ;
}
MPI_Reduce((void *) checksums, (void *) Ychecksums, 3, MPI_DOUBLE,
          MPI_SUM, root, MPI_COMM_WORLD) ;
if ( myid == root ) {
   fprintf(msgFile, 
           "\n\n checksums for original matrix    : %12.4e %12.4e"
           "\n checksums for distributed matrix : %12.4e %12.4e"
           "\n error in checksum                : %12.4e %12.4e",
         Xchecksums[0], Xchecksums[2], Ychecksums[0], Ychecksums[2],
         Ychecksums[0] - Xchecksums[0], Ychecksums[2] - Xchecksums[2]) ;
}
/*
   -----------------------------------------------------
   gather the local matrices into a second global matrix
   -----------------------------------------------------
*/
MARKTIME(t1) ;
stats[0] = stats[1] = stats[2] = stats[3] = 0 ;
if ( myid == root ) {
   Z = DenseMtx_new() ;
   fprintf(msgFile, "\n Z = %p", Z) ;
} else {
   Z = NULL ;
}
Z = DenseMtx_MPI_mergeToGlobalByRows(Z, Y, root, stats, msglvl,
                                     msgFile, tag, MPI_COMM_WORLD) ;
fprintf(msgFile, "\n Z = %p", Z) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : merge into global matrix ", t2 - t1) ;
fprintf(msgFile, "\n send stats : %d messages with %d bytes"
                 "\n recv stats : %d messages with %d bytes",
                 stats[0], stats[2], stats[1], stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, root, MPI_COMM_WORLD) ;
if ( myid == root ) {
   fprintf(msgFile, "\n total send stats : %d messages with %d bytes"
                    "\n total recv stats : %d messages with %d bytes",
                    tstats[0], tstats[2], tstats[1], tstats[3]) ;
   DenseMtx_sort(Z) ;
   fprintf(msgFile, "\n\n global Z matrix") ;
   if ( msglvl > 2 ) {
      DenseMtx_writeForHumanEye(Z, msgFile) ;
   } else {
      DenseMtx_writeStats(Z, msgFile) ;
   }
   fflush(msgFile) ;
   DenseMtx_checksums(Z, Zchecksums) ;
   fprintf(msgFile, 
           "\n\n checksums for original matrix    : %12.4e %12.4e"
           "\n checksums for distributed matrix : %12.4e %12.4e"
           "\n checksums for global matrix      : %12.4e %12.4e"
           "\n X x Y error in checksum          : %12.4e %12.4e"
           "\n X x Z error in checksum          : %12.4e %12.4e",
         Xchecksums[0], Xchecksums[2], 
         Ychecksums[0], Ychecksums[2],
         Zchecksums[0], Zchecksums[2],
         Xchecksums[0] - Ychecksums[0], Xchecksums[2] - Ychecksums[2],
         Xchecksums[0] - Zchecksums[0], Xchecksums[2] - Zchecksums[2]) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
if ( myid == root ) {
   IV_free(mapIV) ;
}
if ( X != NULL ) {
   DenseMtx_free(X) ;
}
if ( Y != NULL ) {
   DenseMtx_free(Y) ;
}
if ( Z != NULL ) {
   DenseMtx_free(Z) ;
}
MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
