/*  testSplitDenseMtx.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------------
   this program tests splitting a DenseMtx object

   (1) process 0 generates a nrow x ncol DenseMtx object
   (2) using a wrap map, distribute the object
   (3) using a random map, distribute the object again

   created -- 98may16, cca
   -----------------------------------------------------
*/
{
char        *buffer ;
DenseMtx    *mtx, *mtxkeep ;
double      t1, t2 ;
double      pchecksums[3], schecksums[3], tchecksums[3] ;
Drand       drand ;
int         inc1, inc2, length, myid, 
            msglvl, ncol, nproc, nrow, seed, tag, type, v ;
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
if ( argc != 9 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile type nrow ncol inc1 inc2 seed "
"\n    msglvl  -- message level"
"\n    msgFile -- message file"
"\n    type    -- type of entries"
"\n      1 -- real"
"\n      2 -- complex"
"\n    nrow    -- number of rows"
"\n    ncol    -- number of columns"
"\n    inc1    -- row increment"
"\n    inc2    -- column increment"
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
        "\n",
        argv[0], msglvl, argv[2], type, nrow, ncol, 
        inc1, inc2, seed) ;
fflush(msgFile) ;
/*
   -----------------------------
   generate the DenseMtx object
   -----------------------------
*/
mtx = DenseMtx_new() ;
fprintf(msgFile, "\n mtx = %p", mtx) ;
fflush(msgFile) ;
MARKTIME(t1) ;
if ( myid == 0 ) {
   DenseMtx_init(mtx, type, 1, -1, nrow, ncol, inc1, inc2) ;
   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0.0, 1.0) ;
   Drand_fillDvector(&drand, nrow*ncol, DenseMtx_entries(mtx)) ;
/*
   ------------------------------------------
   compute the serial checksums of the matrix
   ------------------------------------------
*/
   DenseMtx_checksums(mtx, schecksums) ;
} else {
   if ( inc1 == 1 ) {
      DenseMtx_init(mtx, type, 1, -1, 0, ncol, inc1, 0) ;
   } else {
      DenseMtx_init(mtx, type, 1, -1, 0, ncol, inc1, inc2) ;
   }
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize the DenseMtx object",
        t2 - t1) ;
fflush(msgFile) ;
fprintf(msgFile, "\n\n DenseMtx generated") ;
if ( msglvl > 2 ) {
   DenseMtx_writeForHumanEye(mtx, msgFile) ;
} else {
   DenseMtx_writeStats(mtx, msgFile) ;
}
fflush(msgFile) ;
/*
   ------------------------------------------
   set the initial owners IV to be a wrap map
   ------------------------------------------
*/
MARKTIME(t1) ;
mapIV = IV_new() ;
IV_init(mapIV, nrow, NULL) ;
map = IV_entries(mapIV) ;
for ( v = 0 ; v < nrow ; v++ ) {
   map[v] = v % nproc ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : set wrap map", t2 - t1) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n map") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------
   split the DenseMtx object into pieces
   --------------------------------------
*/
tag = 1 ;
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
mtxkeep = DenseMtx_MPI_splitByRows(mtx, mapIV, stats, msglvl, 
                                   msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : split matrix", t2 - t1) ;
fprintf(msgFile, 
        "\n send stats : %d messages with %d bytes"
        "\n recv stats : %d messages with %d bytes",
        stats[0], stats[2], stats[1], stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
           "\n total send stats : %d messages with %d bytes"
           "\n total recv stats : %d messages with %d bytes",
           tstats[0], tstats[2], tstats[1], tstats[3]) ;
   fflush(msgFile) ;
}
DenseMtx_free(mtx) ;
mtx = mtxkeep ;
tag += 2 ;
fprintf(msgFile, 
        "\n\n after splitting the DenseMtx object with the wrap map") ;
if ( msglvl > 2 ) {
   DenseMtx_writeForHumanEye(mtx, msgFile) ;
} else {
   DenseMtx_writeStats(mtx, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------
   generate a random map
   ---------------------
*/
MARKTIME(t1) ;
Drand_setDefaultFields(&drand) ;
Drand_setSeed(&drand, seed) ;
Drand_setUniform(&drand, 0.0, (double) nrow) ;
for ( v = 0 ; v < nrow ; v++ ) {
   map[v] = ((int) Drand_value(&drand)) % nproc ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : generate random map", t2 - t1) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n map") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   now split the matrix with the random map
   ----------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
MARKTIME(t1) ;
mtxkeep = DenseMtx_MPI_splitByRows(mtx, mapIV, stats, msglvl, 
                                   msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : split matrix ", t2 - t1) ;
fprintf(msgFile, 
        "\n send stats : %d messages with %d bytes"
        "\n recv stats : %d messages with %d bytes",
        stats[0], stats[2], stats[1], stats[3]) ;
fflush(msgFile) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
           "\n total send stats : %d messages with %d bytes"
           "\n total recv stats : %d messages with %d bytes",
           tstats[0], tstats[2], tstats[1], tstats[3]) ;
   fflush(msgFile) ;
}
DenseMtx_free(mtx) ;
mtx = mtxkeep ;
fprintf(msgFile, 
      "\n\n after splitting the DenseMtx object with the owners map") ;
if ( msglvl > 2 ) {
   DenseMtx_writeForHumanEye(mtx, msgFile) ;
} else {
   DenseMtx_writeStats(mtx, msgFile) ;
}
fflush(msgFile) ;
/*
   -----------------------------------------------
   compute the checksums of the distributed matrix
   -----------------------------------------------
*/
DenseMtx_checksums(mtx, pchecksums) ;
MPI_Reduce((void *) pchecksums, (void *) tchecksums, 3, MPI_DOUBLE,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
           "\n\n checksums for original matrix    : %12.4e %12.4e"
           "\n checksums for distributed matrix : %12.4e %12.4e"
           "\n error in checksum                : %12.4e %12.4e",
         schecksums[0], schecksums[2], tchecksums[0], tchecksums[2],
         tchecksums[0] - schecksums[0], tchecksums[2] - schecksums[2]) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
IV_free(mapIV) ;
DenseMtx_free(mtx) ;

MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
