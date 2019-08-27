/*  testSymbFac.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------

   created -- 98may16, cca
   ---------------------------------------------------
*/
{
char         *buffer, *inETreeFileName, *inGraphFileName ;
double       t1, t2 ;
Drand        *drand ;
ETree        *frontETree ;
InpMtx       *inpmtx ;
int          firsttag, ient, ii, jent, J, length, myid, msglvl, nent, 
             nerror, nfront, nproc, nvtx, rc, seed, sizeJ1, sizeJ2, 
             v, vsize ;
int          stats[4], tstats[4] ;
int          *frontOwners, *indicesJ1, *indicesJ2, *ivec1, *ivec2, 
             *vadj, *vtxToFront ;
IV           *frontOwnersIV, *oldToNewIV ;
IVL          *symbfac1, *symbfac2 ;
FILE         *msgFile ;
Graph        *graph ;
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
if ( argc != 6 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile inETreeFile seed"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inGraphFile -- name of Graph file"
      "\n    inETreeFile -- name of ETree file"
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
inGraphFileName = argv[3] ;
inETreeFileName = argv[4] ;
seed            = atoi(argv[5]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inGraphFile -- %s" 
        "\n inETreeFile -- %s" 
        "\n seed        -- %d" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, 
        inETreeFileName, seed) ;
fflush(msgFile) ;
/*
   ----------------------------------
   create the random number generator
   ----------------------------------
*/
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
graph = Graph_new() ;
MARKTIME(t1) ;
rc = Graph_readFromFile(graph, inGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, graph, inGraphFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading Graph object from file %s",
        inGraphFileName) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
} else {
   Graph_writeStats(graph, msgFile) ;
}
fflush(msgFile) ;
nvtx = graph->nvtx ;
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
frontETree = ETree_new() ;
MARKTIME(t1) ;
rc = ETree_readFromFile(frontETree, inETreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in frontETree from file %s",
        t2 - t1, inETreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
           rc, frontETree, inETreeFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading ETree object from file %s",
        inETreeFileName) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(frontETree, msgFile) ;
} else {
   ETree_writeStats(frontETree, msgFile) ;
}
fflush(msgFile) ;
nfront     = ETree_nfront(frontETree) ;
vtxToFront = ETree_vtxToFront(frontETree) ;
/*
   -----------------------------------------
   create an InpMtx object to hold the graph
   -----------------------------------------
*/
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, 0, 0) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   Graph_adjAndSize(graph, v, &vsize, &vadj) ;
   InpMtx_inputRow(inpmtx, v, vsize, vadj) ;
}
InpMtx_sortAndCompress(inpmtx) ;
fprintf(msgFile, "\n\n InpMtx made from graph") ;
if ( msglvl > 2 ) {
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
} else {
   InpMtx_writeStats(inpmtx, msgFile) ;
}
/*
   ----------------------------------------------
   get the permutation vector from the front tree
   ----------------------------------------------
*/
oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
fprintf(msgFile, "\n\n old-to-new permutation vector") ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(oldToNewIV, msgFile) ;
} else {
   IV_writeStats(oldToNewIV, msgFile) ;
}
/*
   --------------------------------------
   permute the vertices in the front tree
   --------------------------------------
*/
ETree_permuteVertices(frontETree, oldToNewIV) ;
fprintf(msgFile, "\n\n front tree after permutation") ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(frontETree, msgFile) ;
} else {
   ETree_writeStats(frontETree, msgFile) ;
}
/*
   -------------------------------------------------
   permute the rows and columns in the InpMtx object
   and convert to chevron format
   -------------------------------------------------
*/
InpMtx_permute(inpmtx, IV_entries(oldToNewIV), IV_entries(oldToNewIV)) ;
InpMtx_changeCoordType(inpmtx, INPMTX_BY_CHEVRONS) ;
InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
fprintf(msgFile, "\n\n InpMtx after permutation") ;
if ( msglvl > 2 ) {
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
} else {
   InpMtx_writeStats(inpmtx, msgFile) ;
}
/*
   -------------------------------------------------------------
   compute the symbolic factorization using the entire structure
   -------------------------------------------------------------
*/
symbfac1 = SymbFac_initFromInpMtx(frontETree, inpmtx) ;
fprintf(msgFile, "\n\n symbolic factorization from global data") ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(symbfac1, msgFile) ;
} else {
   IVL_writeStats(symbfac1, msgFile) ;
}
/*
   ---------------------
   generate a random map
   ---------------------
*/
MARKTIME(t1) ;
Drand_setSeed(drand, seed + 1) ;
Drand_setUniform(drand, 0.0, (double) nproc) ;
frontOwnersIV = IV_new() ;
IV_init(frontOwnersIV, nfront, NULL) ;
frontOwners = IV_entries(frontOwnersIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   frontOwners[J] = (int) Drand_value(drand) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : random map set", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n frontOwners") ;
   IV_writeForHumanEye(frontOwnersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   throw away unowned chevrons
   ---------------------------
*/
InpMtx_changeStorageMode(inpmtx, INPMTX_RAW_DATA) ;
nent  = InpMtx_nent(inpmtx) ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
for ( ient = jent = 0 ; ient < nent ; ient++ ) {
   v = ivec1[ient] ;
   J = vtxToFront[v] ;
   if ( frontOwners[J] == myid ) {
      ivec1[jent] = ivec1[ient] ;
      ivec2[jent] = ivec2[ient] ;
      jent++ ;
   }
}
InpMtx_setNent(inpmtx, jent) ;
InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
fprintf(msgFile, "\n\n InpMtx after dropping unowned entries") ;
if ( msglvl > 2 ) {
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
} else {
   InpMtx_writeStats(inpmtx, msgFile) ;
}
/*
   ------------------------------------------------
   compute the symbolic factorization cooperatively
   ------------------------------------------------
*/
firsttag = 47 ;
symbfac2 = SymbFac_MPI_initFromInpMtx(frontETree, frontOwnersIV, inpmtx,
                     stats, msglvl, msgFile, firsttag, MPI_COMM_WORLD) ;
fprintf(msgFile, "\n\n symbolic factorization from cooperation") ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(symbfac2, msgFile) ;
} else {
   IVL_writeStats(symbfac2, msgFile) ;
}
/*
   ----------------------------------------
   check that the supported part of the 
   second symbolic factorization is correct
   ----------------------------------------
*/
for ( J = 0, nerror = 0 ; J < nfront ; J++ ) {
   IVL_listAndSize(symbfac2, J, &sizeJ2, &indicesJ2) ;
   if ( sizeJ2 > 0 ) {
      IVL_listAndSize(symbfac1, J, &sizeJ1, &indicesJ1) ;
      if ( sizeJ1 != sizeJ2 ) {
         fprintf(msgFile, "\n error, J = %d, sizeJ1 = %d, sizeJ2 = %d",
                 J, sizeJ1, sizeJ2) ;
         nerror++ ;
      } else {
         for ( ii = 0 ; ii < sizeJ1 ; ii++ ) {
            if ( indicesJ1[ii] != indicesJ2[ii] ) {
               fprintf(msgFile, "\n error, J = %d", J) ;
               fprintf(msgFile, "\n indicesJ1") ;
               IVfprintf(msgFile, sizeJ1, indicesJ1) ;
               fprintf(msgFile, "\n indicesJ2") ;
               IVfprintf(msgFile, sizeJ1, indicesJ2) ;
               nerror++ ;
               break ;
            }
         }
      }
   }
}
fprintf(msgFile, "\n\n %d errors in symbolic factorization", nerror) ;
/*
   ----------------
   free the objects
   ----------------
*/
IV_free(frontOwnersIV) ;
IV_free(oldToNewIV) ;
InpMtx_free(inpmtx) ;
ETree_free(frontETree) ;
Drand_free(drand) ;

MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
