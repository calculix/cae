/*  testMMM2.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

#define MMM_WITH_A   0
#define MMM_WITH_AT  1
#define MMM_WITH_AH  2

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------
   (1) process 0 creates Y, A and X
   (2) using a random map, it distributes X
   (3) using a random map, it distributes Y
   (4) using a random map, it distributes A

   created -- 98aug21, cca
   ---------------------------------------------------
*/
{
char         *buffer ;
DenseMtx     *X, *Xloc, *Y, *Ykeep, *Yloc, *Z ;
double       imag, real, t1, t2 ;
double       alpha[2] ;
Drand        *drand ;
InpMtx       *A, *Aloc ;
int          coordType, inputMode, length, myid, 
             msglvl, ncolA, ncolX, nentA, nproc, nrowA, 
             nrowX, nrowY, opflag, seed, symflag, tag ;
int          stats[4], tstats[4] ;
IV           *AmapIV, *mapIV, *XmapIV, *XownedIV, *YmapIV, *YownedIV ;
FILE         *msgFile ;
MatMulInfo   *info ;
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
if ( argc != 14 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile nrowA ncolA nentA ncolX"
      "\n         coordType inputMode symflag opflag seed real imag"
      "\n    msglvl    -- message level"
      "\n    msgFile   -- message file"
      "\n    nrowA     -- number of rows in A"
      "\n    ncolA     -- number of columns in A"
      "\n    nentA     -- number of entries in A"
      "\n    ncolX     -- number of columns in X"
      "\n    coordType -- coordinate type"
      "\n       1 -- store by rows"
      "\n       2 -- store by columns"
      "\n       3 -- store by chevrons"
      "\n    inputMode -- input mode"
      "\n       1 -- indices and real entries"
      "\n       2 -- indices and complex entries"
      "\n    symflag -- symmetry flag"
      "\n       0 -- symmetric"
      "\n       1 -- hermitian"
      "\n       2 -- nonsymmetric"
      "\n    opflag -- operations flag"
      "\n       0 -- multiply with A"
      "\n       1 -- multiply with A^T"
      "\n       2 -- multiply with A^H"
      "\n    seed -- random number seed"
      "\n    real -- real part of scalar alpha"
      "\n    imag -- real part of scalar alpha"
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
nrowA     = atoi(argv[3]) ;
ncolA     = atoi(argv[4]) ;
nentA     = atoi(argv[5]) ;
ncolX     = atoi(argv[6]) ;
coordType = atoi(argv[7]) ;
inputMode = atoi(argv[8]) ;
symflag   = atoi(argv[9]) ;
opflag    = atoi(argv[10]) ;
seed      = atoi(argv[11]) ;
real      = atof(argv[12]) ;
imag      = atof(argv[13]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl    -- %d" 
        "\n msgFile   -- %s" 
        "\n nrowA     -- %d" 
        "\n ncolA     -- %d" 
        "\n nentA     -- %d" 
        "\n ncolX     -- %d" 
        "\n coordType -- %d" 
        "\n inputMode -- %d" 
        "\n symflag   -- %d" 
        "\n opflag    -- %d" 
        "\n seed      -- %d" 
        "\n real      -- %e" 
        "\n imag      -- %e" 
        "\n",
        argv[0], msglvl, argv[2], nrowA, ncolA, nentA, ncolX,
        coordType, inputMode, symflag, opflag, seed, real, imag) ;
fflush(msgFile) ;
/*
   --------------------
   check the input data
   --------------------
*/
if ( nrowA <= 0 || ncolA <= 0 || nentA <= 0 || ncolX <= 0 ) {
   fprintf(stderr, "\n fatal error, bad dimensions") ;
   exit(-1) ;
}
switch ( coordType ) {
case INPMTX_BY_ROWS :
case INPMTX_BY_COLUMNS :
case INPMTX_BY_CHEVRONS :
   break ;
default :
   fprintf(stderr, "\n bad coordType") ;
   exit(-1) ;
   break ;
}
switch ( inputMode ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n bad type") ;
   exit(-1) ;
   break ;
}
switch ( symflag ) {
case SPOOLES_SYMMETRIC :
   if ( nrowA != ncolA ) {
      fprintf(stderr, 
              "\n fatal error, symmetric but nrowA = %d, ncolA = %d", 
              nrowA, ncolA) ;
      exit(-1) ;
   }
   nrowX = nrowA ;
   nrowY = nrowA ;
   break ;
case SPOOLES_HERMITIAN :
   if ( nrowA != ncolA ) {
      fprintf(stderr, 
              "\n fatal error, hermitian but nrowA = %d, ncolA = %d", 
              nrowA, ncolA) ;
      exit(-1) ;
   }
   if ( inputMode != SPOOLES_COMPLEX ) {
      fprintf(stderr, "\n fatal error, hermitian but not complex") ;
      exit(-1) ;
   }
   nrowX = nrowA ;
   nrowY = nrowA ;
   break ;
case SPOOLES_NONSYMMETRIC :
   if ( opflag == MMM_WITH_AT || opflag == MMM_WITH_AH ) {
      nrowX = nrowA ;
      nrowY = ncolA ;
   } else {
      nrowX = ncolA ;
      nrowY = nrowA ;
   }
   break ;
default :
   fprintf(stderr, "\n bad symflag") ;
   exit(-1) ;
   break ;
}
switch ( opflag ) {
case MMM_WITH_A :
case MMM_WITH_AT :
case MMM_WITH_AH :
   break ;
default :
   fprintf(stderr, "\n bad opflag") ;
   exit(-1) ;
   break ;
}
alpha[0] = real ;
alpha[1] = imag ;
/*
   ----------------------------------
   create the random number generator
   ----------------------------------
*/
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
if ( myid == 0 ) {
/*
   ---------------------------------------
   processor 0, generate a random matrices
   ---------------------------------------
*/
   MARKTIME(t1) ;
   A = InpMtx_new() ;
   InpMtx_init(A, INPMTX_BY_ROWS, inputMode, nentA, 0) ;
   Drand_setUniform(drand, 0, nrowA) ;
   Drand_fillIvector(drand, nentA, InpMtx_ivec1(A)) ;
   Drand_setUniform(drand, 0, ncolA) ;
   Drand_fillIvector(drand, nentA, InpMtx_ivec2(A)) ;
   Drand_setUniform(drand, -1.0, 1.0) ;
   if ( inputMode == SPOOLES_REAL ) {
      Drand_fillDvector(drand, nentA, InpMtx_dvec(A)) ;
   } else if ( inputMode == SPOOLES_COMPLEX ) {
      Drand_fillDvector(drand, 2*nentA, InpMtx_dvec(A)) ;
   }
   A->nent = nentA ;
   InpMtx_sortAndCompress(A) ;
   InpMtx_changeCoordType(A, coordType) ;
   InpMtx_changeStorageMode(A, INPMTX_BY_VECTORS) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %9.5f : generate random matrix A", t2 - t1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n global matrix A") ;
      InpMtx_writeForHumanEye(A, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -----------------
   create a matrix X
   -----------------
*/
   MARKTIME(t1) ;
   X = DenseMtx_new() ;
   DenseMtx_init(X, inputMode, 1, -1, nrowX, ncolX, 1, nrowX) ;
   if ( inputMode == SPOOLES_REAL ) {
      Drand_fillDvector(drand, nrowX*ncolX, DenseMtx_entries(X)) ;
   } else if ( inputMode == SPOOLES_COMPLEX ) {
      Drand_fillDvector(drand, 2*nrowX*ncolX, DenseMtx_entries(X)) ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n global matrix X") ;
      DenseMtx_writeForHumanEye(X, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -----------------
   create a matrix Y
   -----------------
*/
   MARKTIME(t1) ;
   Y = DenseMtx_new() ;
   DenseMtx_init(Y, inputMode, 1, -1, nrowY, ncolX, 1, nrowY) ;
   if ( inputMode == SPOOLES_REAL ) {
      Drand_fillDvector(drand, nrowY*ncolX, DenseMtx_entries(Y)) ;
   } else if ( inputMode == SPOOLES_COMPLEX ) {
      Drand_fillDvector(drand, 2*nrowY*ncolX, DenseMtx_entries(Y)) ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n global matrix Y") ;
      DenseMtx_writeForHumanEye(Y, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -------------------------------
   create Ykeep and copy Y into it
   -------------------------------
*/
   MARKTIME(t1) ;
   Ykeep = DenseMtx_new() ;
   DenseMtx_init(Ykeep, inputMode, 1, -1, nrowY, ncolX, 1, nrowY) ;
   if ( inputMode == SPOOLES_REAL ) {
      DVcopy(nrowY*ncolX, 
             DenseMtx_entries(Ykeep), DenseMtx_entries(Y)) ;
   } else if ( inputMode == SPOOLES_COMPLEX ) {
      DVcopy(2*nrowY*ncolX, 
             DenseMtx_entries(Ykeep), DenseMtx_entries(Y)) ;
   }
/*
   -----------------------------------------
   compute the serial matrix-matrix multiply
   -----------------------------------------
*/
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
      InpMtx_sym_mmm(A, Ykeep, alpha, X) ;
      break ;
   case SPOOLES_HERMITIAN :
      InpMtx_herm_mmm(A, Ykeep, alpha, X) ;
      break ;
   case SPOOLES_NONSYMMETRIC :
      switch ( opflag ) {
      case MMM_WITH_A :
         InpMtx_nonsym_mmm(A, Ykeep, alpha, X) ;
         break ;
      case MMM_WITH_AT :
         InpMtx_nonsym_mmm_T(A, Ykeep, alpha, X) ;
         break ;
      case MMM_WITH_AH :
         InpMtx_nonsym_mmm_H(A, Ykeep, alpha, X) ;
         break ;
      }
      break ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n Ykeep") ;
      DenseMtx_writeForHumanEye(Ykeep, msgFile) ;
      fflush(msgFile) ;
   }
} else {
/*
   -----------------------------
   initialize the empty matrices
   -----------------------------
*/
   A = InpMtx_new() ;
   InpMtx_init(A, coordType, inputMode, 0, 0) ;
   X = DenseMtx_new() ;
   DenseMtx_init(X, inputMode, 1, -1, 0, ncolX, 1, 0) ;
   Y = DenseMtx_new() ;
   DenseMtx_init(Y, inputMode, 1, -1, 0, ncolX, 1, 0) ;
   Ykeep = NULL ;
}
/*
   --------------------------------------
   set the A owners IV to be a random map
   --------------------------------------
*/
Drand_setSeed(drand, ++seed) ;
Drand_setUniform(drand, 0, nproc) ;
AmapIV = IV_new() ;
switch ( coordType ) {
case INPMTX_BY_ROWS :
   IV_init(AmapIV, nrowA, NULL) ;
   Drand_fillIvector(drand, nrowA, IV_entries(AmapIV)) ;
   break ;
case INPMTX_BY_COLUMNS :
   IV_init(AmapIV, ncolA, NULL) ;
   Drand_fillIvector(drand, ncolA, IV_entries(AmapIV)) ;
   break ;
case INPMTX_BY_CHEVRONS :
   IV_init(AmapIV, nrowA, NULL) ;
   Drand_fillIvector(drand, nrowA, IV_entries(AmapIV)) ;
   break ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n A's map") ;
   IV_writeForHumanEye(AmapIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   split the InpMtx object into pieces
   ------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
tag = 1 ;
MARKTIME(t1) ;
InpMtx_changeStorageMode(A, INPMTX_BY_VECTORS) ;
Aloc = InpMtx_MPI_split(A, AmapIV, 
                        stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : A matrix split", t2 - t1) ;
fprintf(msgFile, "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
InpMtx_free(A) ;
A = NULL ;
tag += 2 ;
fprintf(msgFile, "\n\n local InpMtx object ") ;
if ( msglvl > 1 ) {
   InpMtx_writeForHumanEye(Aloc, msgFile) ;
} else {
   InpMtx_writeStats(Aloc, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------
   set the Y owners IV to be a random map
   --------------------------------------
*/
Drand_setSeed(drand, ++seed) ;
Drand_setUniform(drand, 0, nproc) ;
YmapIV = IV_new() ;
IV_init(YmapIV, nrowY, NULL) ;
Drand_fillIvector(drand, nrowY, IV_entries(YmapIV)) ;
YownedIV = IV_targetEntries(YmapIV, myid) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n Y's map") ;
   IV_writeForHumanEye(YmapIV, msgFile) ;
   fprintf(msgFile, "\n\n owned rows of Y") ;
   IV_writeForHumanEye(YownedIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   split the Y DenseMtx object into pieces
   ---------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
tag++ ;
MARKTIME(t1) ;
Yloc = DenseMtx_MPI_splitByRows(Y, YmapIV, stats, msglvl, 
                                msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : Y matrix split", t2 - t1) ;
fprintf(msgFile, "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
DenseMtx_free(Y) ;
Y = NULL ;
tag += 2 ;
fprintf(msgFile, "\n\n local Y DenseMtx object ") ;
if ( msglvl > 1 ) {
   DenseMtx_writeForHumanEye(Yloc, msgFile) ;
} else {
   DenseMtx_writeStats(Yloc, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------
   set the X owners IV to be a random map
   --------------------------------------
*/
Drand_setSeed(drand, ++seed) ;
Drand_setUniform(drand, 0, nproc) ;
XmapIV = IV_new() ;
IV_init(XmapIV, nrowX, NULL) ;
Drand_fillIvector(drand, nrowX, IV_entries(XmapIV)) ;
XownedIV = IV_targetEntries(XmapIV, myid) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n X's map") ;
   IV_writeForHumanEye(XmapIV, msgFile) ;
   fprintf(msgFile, "\n\n owned rows of X") ;
   IV_writeForHumanEye(XownedIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   split the X DenseMtx object into pieces
   ---------------------------------------
*/
stats[0]  = stats[1]  = stats[2]  = stats[3]  = 0 ;
tstats[0] = tstats[1] = tstats[2] = tstats[3] = 0 ;
tag++ ;
MARKTIME(t1) ;
Xloc = DenseMtx_MPI_splitByRows(X, XmapIV, stats, msglvl, 
                                msgFile, tag, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : X matrix split", t2 - t1) ;
fprintf(msgFile, "\n local comm  : %6d %12d   %6d %12d",
        stats[0],  stats[2],  stats[1],  stats[3]) ;
MPI_Reduce((void *) stats, (void *) tstats, 4, MPI_INT,
          MPI_SUM, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n global comm : %6d %12d   %6d %12d",
           tstats[0],  tstats[2],  tstats[1],  tstats[3]) ;
   fflush(msgFile) ;
}
DenseMtx_free(X) ;
X = NULL ;
tag += 2 ;
fprintf(msgFile, 
      "\n\n after splitting the X DenseMtx object with the wrap map") ;
if ( msglvl > 1 ) {
   DenseMtx_writeForHumanEye(Xloc, msgFile) ;
} else {
   DenseMtx_writeStats(Xloc, msgFile) ;
}
fflush(msgFile) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   initialize the MatMulInfo data structure
   ----------------------------------------
*/
IVfill(4, stats, 0) ;
tag  = 1 ;
info = MatMul_MPI_setup(Aloc, symflag, opflag, XmapIV, YmapIV, 
                        stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
/*
   -------------------------------------------
   map the input matrix into local coordinates
   -------------------------------------------
*/
MatMul_setLocalIndices(info, Aloc) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n local A") ;
   InpMtx_writeForHumanEye(Aloc, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   compute the matrix-matrix multiply
   ----------------------------------
*/
MatMul_MPI_mmm(info, Yloc, alpha, Aloc, Xloc, 
               stats, msglvl, msgFile, tag, MPI_COMM_WORLD) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after mmm, local matrix Y") ;
   DenseMtx_writeForHumanEye(Yloc, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------------------
   gather all of Yloc onto processor 0 and compare with Ykeep
   ----------------------------------------------------------
*/
mapIV = IV_new() ;
IV_init(mapIV, nrowY, NULL) ;
IV_fill(mapIV, 0) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n mapIV to gather Y") ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
}
Z = DenseMtx_MPI_splitByRows(Yloc, mapIV, stats, msglvl, 
                             msgFile, tag, MPI_COMM_WORLD) ;
DenseMtx_free(Yloc) ;
Yloc = Z ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n Yloc after split") ;
   DenseMtx_writeForHumanEye(Yloc, msgFile) ;
   fflush(msgFile) ;
}
if ( myid == 0 ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n global Yloc") ;
      DenseMtx_writeForHumanEye(Yloc, msgFile) ;
      fflush(msgFile) ;
   }
   DenseMtx_sub(Ykeep, Yloc) ;
   fprintf(msgFile, "\n\n error = %12.4e", DenseMtx_maxabs(Ykeep)) ;
} 
/*
   ------------------------------
   clean up the MatMulInfo object
   ------------------------------
*/
MatMul_cleanup(info) ;
info = NULL ;
/*
   ----------------
   free the objects
   ----------------
*/
IV_free(AmapIV) ;
IV_free(XmapIV) ;
IV_free(XownedIV) ;
IV_free(YmapIV) ;
IV_free(YownedIV) ;
InpMtx_free(Aloc) ;
DenseMtx_free(Xloc) ;
DenseMtx_free(Yloc) ;
if ( Ykeep != NULL ) {
   DenseMtx_free(Ykeep) ;
} 
IV_free(mapIV) ;
Drand_free(drand) ;

MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
