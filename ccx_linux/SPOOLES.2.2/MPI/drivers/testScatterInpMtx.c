/*  testScatterInpMtx.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   (1) process root generates an InpMtx object
   (2) using a random map, it distributes the object 

   created -- 98sep27, cca
   ---------------------------------------------------
*/
{
char         *buffer ;
InpMtx       *mtxA, *mtxB ;
double       t1, t2 ;
double       schecksums[3], pchecksums[3], tchecksums[3] ;
Drand        *drand ;
int          coordType, inputMode, iproc, length, myid, 
             msglvl, nent, neqns, nproc, rc, root, seed, tag, v ;
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
if ( argc != 9 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile neqns nent seed "
      "\n         coordType inputMode root"
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
      "\n    root      -- root processor that creates the global matrix"
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
root      = atoi(argv[8]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl    -- %d" 
        "\n msgFile   -- %s" 
        "\n neqns     -- %d" 
        "\n nent      -- %d" 
        "\n seed      -- %d" 
        "\n coordType -- %d" 
        "\n inputMode -- %d" 
        "\n root      -- %d" 
        "\n",
        argv[0], msglvl, argv[2], neqns, nent, seed, 
        coordType, inputMode, root) ;
fflush(msgFile) ;
/*
   ----------------------------------
   create the random number generator
   ----------------------------------
*/
if ( myid == root ) {
/*
   ----------------------------------------------
   processor root, generate a random input matrix
   ----------------------------------------------
*/
   MARKTIME(t1) ;
   mtxA = InpMtx_new() ;
   InpMtx_init(mtxA, INPMTX_BY_ROWS, inputMode, nent, 0) ;
   drand = Drand_new() ;
   Drand_setSeed(drand, seed) ;
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
   fprintf(msgFile, "\n schecksums %12.4e %12.4e %12.4e",
           schecksums[0], schecksums[1], schecksums[2]) ;
   fflush(msgFile) ;
/*
   ---------------------
   generate a random map
   ---------------------
*/
   MARKTIME(t1) ;
   Drand_setSeed(drand, seed + 1) ;
   Drand_setUniform(drand, 0.0, (double) neqns) ;
   mapIV = IV_new() ;
   IV_init(mapIV, neqns, NULL) ;
   map = IV_entries(mapIV) ;
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
} else {
   mtxA  = NULL ;
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
mtxB = InpMtx_MPI_splitFromGlobal(mtxA, NULL, mapIV, root, stats, 
                                 msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : matrix split", t2 - t1) ;
fprintf(msgFile, 
        "\n send stats : %d messages with %d bytes"
        "\n recv stats : %d messages with %d bytes",
        stats[0], stats[2], stats[1], stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, root, MPI_COMM_WORLD) ;
if ( myid == root ) {
   fprintf(msgFile, 
           "\n total send stats : %d messages with %d bytes"
           "\n total recv stats : %d messages with %d bytes",
           tstats[0], tstats[2], tstats[1], tstats[3]) ;
   fflush(msgFile) ;
}
if ( mtxA != NULL ) {
   InpMtx_free(mtxA) ;
}
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
           MPI_SUM, root, MPI_COMM_WORLD) ;
if ( myid == root ) {
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
if ( myid == root ) {
   IV_free(mapIV) ;
   Drand_free(drand) ;
}
InpMtx_free(mtxA) ;

MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
