/*  testSplitInpMtx.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   (1) process 0 reads in a InpMtx object. 
   (2) using a wrap map, distribute the object
   (3) using a random map, distribute the object again

   created -- 98may16, cca
   ---------------------------------------------------
*/
{
char         *buffer ;
InpMtx       *mtxA, *mtxB ;
double       t1, t2 ;
double       schecksums[3], pchecksums[3], tchecksums[3] ;
Drand        *drand ;
int          coordType, inputMode, iproc, length, myid, 
             msglvl, nent, neqns, nproc, rc, seed, tag, v ;
int          ibuffer[3], stats[4], tstats[4] ;
int          *map ;
IV           *mapIV ;
FILE         *msgFile ;
MPI_Status   status ;
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
if ( argc != 8 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile neqns nent seed "
      "\n         coordType inputMode "
      "\n    msglvl    -- message level"
      "\n    msgFile   -- message file"
      "\n    neqns     -- number of equations"
      "\n    nent      -- number of equations"
      "\n    seed      -- random number seed"
      "\n    coordType -- coordinate type"
      "\n       1 -- store by rows"
      "\n       2 -- store by columns"
      "\n       3 -- store by chevrons"
      "\n    inputMode -- input mode"
      "\n       0 -- indices only"
      "\n       1 -- indices and real entries"
      "\n       2 -- indices and complex entries"
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
neqns     = atoi(argv[3]) ;
nent      = atoi(argv[4]) ;
seed      = atoi(argv[5]) ;
coordType = atoi(argv[6]) ;
inputMode = atoi(argv[7]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl    -- %d" 
        "\n msgFile   -- %s" 
        "\n neqns     -- %d" 
        "\n nent      -- %d" 
        "\n seed      -- %d" 
        "\n coordType -- %d" 
        "\n inputMode -- %d" 
        "\n",
        argv[0], msglvl, argv[2], neqns, nent, seed, 
        coordType, inputMode) ;
fflush(msgFile) ;
/*
   ----------------------------------
   create the random number generator
   ----------------------------------
*/
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
if ( myid == 0 ) {
/*
   -------------------------------------------
   processor 0, generate a random input matrix
   -------------------------------------------
*/
   MARKTIME(t1) ;
   mtxA = InpMtx_new() ;
   InpMtx_init(mtxA, INPMTX_BY_ROWS, inputMode, nent, 0) ;
   Drand_setUniform(drand, 0, neqns) ;
   Drand_fillIvector(drand, nent, InpMtx_ivec1(mtxA)) ;
   Drand_fillIvector(drand, nent, InpMtx_ivec2(mtxA)) ;
   if ( inputMode == SPOOLES_REAL ) {
      Drand_setUniform(drand, -1.0, 1.0) ;
      Drand_fillDvector(drand, nent, InpMtx_dvec(mtxA)) ;
   } else if ( inputMode == SPOOLES_COMPLEX ) {
      Drand_setUniform(drand, -1.0, 1.0) ;
      Drand_fillDvector(drand, 2*nent, InpMtx_dvec(mtxA)) ;
   }
   mtxA->nent = nent ;
   InpMtx_sortAndCompress(mtxA) ;
   InpMtx_changeCoordType(mtxA, coordType) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : generate random matrix", t2 - t1) ;
/*
   -----------------------------------------------
   change the coordinate type and the storage mode
   -----------------------------------------------
*/
   InpMtx_changeCoordType(mtxA, coordType) ;
   InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
   mtxA->inputMode = inputMode ;
   fprintf(msgFile, "\n\n random input matrix") ;
   if ( msglvl > 2 ) {
      InpMtx_writeForHumanEye(mtxA, msgFile) ;
   } else {
      InpMtx_writeStats(mtxA, msgFile) ;
   }
   fflush(msgFile) ;
/*
   ------------------------------------------
   compute the serial checksums of the matrix
   ------------------------------------------
*/
   InpMtx_checksums(mtxA, schecksums) ;
/*
   -----------------------------------
   send the first message identifying 
   the coordinate type, and input mode
   -----------------------------------
*/
   tag = 1 ;
   ibuffer[0] = InpMtx_coordType(mtxA) ;
   ibuffer[1] = InpMtx_inputMode(mtxA) ;
   MARKTIME(t1) ;
   for ( iproc = 1 ; iproc < nproc ; iproc++ ) {
      MPI_Send((void *) ibuffer, 2, MPI_INT, 
               iproc, tag, MPI_COMM_WORLD);
   }
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : first message sent", t2 - t1) ;
} else {
/*
   --------------------------------------------
   receive the first message identifying the
   coordinate type, storage mode and input mode
   --------------------------------------------
*/
   tag = 1 ;
   MARKTIME(t1) ;
   MPI_Recv((void *) ibuffer, 2, MPI_INT, 
            0, tag, MPI_COMM_WORLD, &status) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : first message received", t2 - t1) ;
   mtxA = InpMtx_new() ;
   InpMtx_init(mtxA, ibuffer[0], ibuffer[1], 0, 0) ;
   InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
   fprintf(msgFile, "\n\n after receiving initial information") ;
   if ( msglvl > 2 ) {
      InpMtx_writeForHumanEye(mtxA, msgFile) ;
   } else {
      InpMtx_writeStats(mtxA, msgFile) ;
   }
}
/*
   ------------------------------------------
   set the initial owners IV to be a wrap map
   ------------------------------------------
*/
mapIV = IV_new() ;
IV_init(mapIV, neqns, NULL) ;
map = IV_entries(mapIV) ;
for ( v = 0 ; v < neqns ; v++ ) {
   map[v] = v % nproc ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n map") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   split the InpMtx object into pieces
   ------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
tag++ ;
MARKTIME(t1) ;
mtxB = InpMtx_MPI_split(mtxA, mapIV, stats,
                        msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : matrix split", t2 - t1) ;
fprintf(msgFile, 
        "\n send stats : %d messages with %d bytes"
        "\n recv stats : %d messages with %d bytes",
        stats[0], stats[2], stats[1], stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
           "\n total send stats : %d messages with %d bytes"
           "\n total recv stats : %d messages with %d bytes",
           tstats[0], tstats[2], tstats[1], tstats[3]) ;
   fflush(msgFile) ;
}
InpMtx_free(mtxA) ;
mtxA = mtxB ;
tag += 2 ;
if ( mtxA != NULL ) {
   InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
   fprintf(msgFile, 
          "\n\n after splitting the InpMtx object with the wrap map") ;
   if ( msglvl > 2 ) {
      InpMtx_writeForHumanEye(mtxA, msgFile) ;
   } else {
      InpMtx_writeStats(mtxA, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   ---------------------
   generate a random map
   ---------------------
*/
MARKTIME(t1) ;
Drand_setSeed(drand, seed + 1) ;
Drand_setUniform(drand, 0.0, (double) neqns) ;
for ( v = 0 ; v < neqns ; v++ ) {
   map[v] = ((int) Drand_value(drand)) % nproc ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : random map set", t2 - t1) ;
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
mtxB = InpMtx_MPI_split(mtxA, mapIV, stats,
                        msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : matrix split", t2 - t1) ;
fprintf(msgFile, 
        "\n send stats : %d messages with %d bytes"
        "\n recv stats : %d messages with %d bytes",
        stats[0], stats[2], stats[1], stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, 
           "\n total send stats : %d messages with %d bytes"
           "\n total recv stats : %d messages with %d bytes",
           tstats[0], tstats[2], tstats[1], tstats[3]) ;
   fflush(msgFile) ;
}
InpMtx_free(mtxA) ;
mtxA = mtxB ;
fprintf(msgFile, 
        "\n\n after splitting the InpMtx object with the owners map") ;
if ( msglvl > 2 ) {
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
} else {
   InpMtx_writeStats(mtxA, msgFile) ;
}
fflush(msgFile) ;
/*
   -----------------------------------------------
   compute the checksums of the distributed matrix
   -----------------------------------------------
*/
InpMtx_checksums(mtxA, pchecksums) ;
MPI_Reduce((void *) pchecksums, (void *) tchecksums, 3, MPI_DOUBLE,
           MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile,
          "\n\n checksums for original matrix    : %12.4e %12.4e %12.4e"
          "\n checksums for distributed matrix : %12.4e %12.4e %12.4e"
          "\n error in checksum                : %12.4e %12.4e %12.4e",
          schecksums[0], schecksums[1], schecksums[2], 
          tchecksums[0], tchecksums[1], tchecksums[2],
          tchecksums[0] - schecksums[0], 
          tchecksums[1] - schecksums[1], 
          tchecksums[2] - schecksums[2]) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
IV_free(mapIV) ;
InpMtx_free(mtxA) ;
Drand_free(drand) ;

MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
