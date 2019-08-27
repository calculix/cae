/*  testR2D.c  */

#include "../InpMtx.h"
#include "../../EGraph.h"
#include "../../Coords.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

void
generateMatrix (
   int      eadj[],
   Coords   *coords,
   double   mtxent[],
   int      msglvl,
   FILE     *msgFile
) ;
void
loadMatrices (
   EGraph   *egraph,
   Coords   *coords,
   InpMtx  *inpmtx,
   int      msglvl,
   FILE     *msgFile
) ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   test the InpMtx object

   (1) read in an EGraph file that contains the vertex lists
       for triangular element
   (2) read in a Coords file that contains the vertex coordinates
   (3) for each element
       (4) create an elemental matrix and insert into InpMtx
       (5) assemble the matrix

   created -- 96jul05, cca
   --------------------------------------------------------------
*/
int
main ( int argc, char *argv[] ) 
{
char       *CoordsFileName, *EGraphFileName, *outFileName ;
Coords     *coords ;
DenseMtx   *X, *Y, *Y2 ;
double     alpha[2] = {1.0, 0.0} ;
double     *x, *y ;
double     *mtxent ;
double     t1, t2 ;
Drand      *drand ;
EGraph     *egraph ;
FILE       *msgFile ;
int        coordType, esize, ielem, msglvl, nelem, nvtx, rc, seed ;
int        *colOldToNew, *eadj, *rowOldToNew ;
InpMtx     *inpmtx ;
IV         *colOldToNewIV, *rowOldToNewIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( argc != 8 ) {
   fprintf(stdout, 
       "\n\n usage : %s msglvl msgFile EGraphFile CoordsFile "
       "\n         coordType seed outFile"
       "\n    msglvl     -- message level"
       "\n    msgFile    -- message file"
       "\n    EGraphFile -- file that contains the element graph"
       "\n                  must be *.egraphf or *.egraphb"
       "\n    CoordsFile -- file that contains the coordinates object"
       "\n                  must be *.coordsf or *.coordsb"
       "\n    coordType  -- coordinate type"
       "\n       1 --> by rows"
       "\n       2 --> by columns"
       "\n       3 --> by chevrons"
       "\n    seed       -- random number seed"
       "\n    outFile    -- file to contain the InpMtx object"
       "\n                  must be *.inpmtxb or *.inpmtxf"
       "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
EGraphFileName = argv[3] ;
CoordsFileName = argv[4] ;
coordType      = atoi(argv[5]) ;
seed           = atoi(argv[6]) ;
outFileName    = argv[7] ;
/*
   --------------
   echo the input
   --------------
*/
fprintf(msgFile, "\n input to %s"
        "\n msglvl      = %d"
        "\n msgFile     = %s"
        "\n EGraphFile  = %s"
        "\n CoordsFile  = %s"
        "\n coordType   = %d"
        "\n seed        = %d"
        "\n",
        argv[0], msglvl, argv[2], EGraphFileName, CoordsFileName,
        coordType, seed) ;
fflush(msgFile) ;
/*
   -------------------------
   read in the EGraph object
   -------------------------
*/
egraph = EGraph_new() ;
MARKTIME(t1) ;
rc = EGraph_readFromFile(egraph, EGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in egraph from file %s",
        t2 - t1, EGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from EGraph_readFromFile(%p,%s)",
           rc, egraph, EGraphFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading EGraph object from file %s",
           EGraphFileName) ;
   EGraph_writeForHumanEye(egraph, msgFile) ;
   fflush(msgFile) ;
}
nvtx  = egraph->nvtx  ;
nelem = egraph->nelem ;
/*
   -------------------------
   read in the Coords object
   -------------------------
*/
coords = Coords_new() ;
MARKTIME(t1) ;
rc = Coords_readFromFile(coords, CoordsFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in egraph from file %s",
        t2 - t1, CoordsFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from Coords_readFromFile(%p,%s)",
           rc, coords, CoordsFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading Coords object from file %s",
           CoordsFileName) ;
   Coords_writeForHumanEye(coords, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   create the random number object
   -------------------------------
*/
drand = Drand_new() ;
Drand_init(drand) ;
Drand_setUniform(drand, -1.0, 1.0) ;
Drand_setSeed(drand, seed) ;
/*
   -------------------------
   create the InpMtx object
   -------------------------
*/
MARKTIME(t1) ;
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, coordType, INPMTX_REAL_ENTRIES, 
            egraph->adjIVL->tsize, 0) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : initialize the InpMtx object", 
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n InpMtx after initialization") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
}
/*
   --------------------------
   fill X with random numbers
   --------------------------
*/
nvtx = egraph->nvtx ;
X = DenseMtx_new() ;
DenseMtx_init(X, SPOOLES_REAL, 0, 0, nvtx, 1, 1, nvtx) ;
DenseMtx_fillRandomEntries(X, drand) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n solution vector") ;
   DenseMtx_writeForHumanEye(X, msgFile) ;
}
x = DenseMtx_entries(X) ;
/*
   -----------------
   fill Y with zeros
   -----------------
*/
Y = DenseMtx_new() ;
DenseMtx_init(Y, SPOOLES_REAL, 0, 0, nvtx, 1, 1, nvtx) ;
DenseMtx_zero(Y) ;
y = DenseMtx_entries(Y) ;
/*
   ---------------------------------
   load the matrices and compute y[]
   ---------------------------------
*/
MARKTIME(t1) ;
mtxent = DVinit(9, 0.0) ;
for ( ielem = 0 ; ielem < nelem ; ielem++ ) {
   IVL_listAndSize(egraph->adjIVL, ielem, &esize, &eadj) ;
   if ( msglvl > 3 && msgFile != NULL ) {
      fprintf(msgFile, "\n\n elemental matrix %d", ielem) ;
   }
   generateMatrix(eadj, coords, mtxent, msglvl, msgFile) ;
   InpMtx_inputRealMatrix(inpmtx, 3, 3, 1, 3, eadj, eadj, mtxent) ;
   if ( msglvl > 3 && msgFile != NULL ) {
      InpMtx_writeForHumanEye(inpmtx, stdout) ;
   }
   y[eadj[0]] += mtxent[0] * x[eadj[0]] 
              +  mtxent[3] * x[eadj[1]] + mtxent[6] * x[eadj[2]] ;
   y[eadj[1]] += mtxent[1] * x[eadj[0]] 
              +  mtxent[4] * x[eadj[1]] + mtxent[7] * x[eadj[2]] ;
   y[eadj[2]] += mtxent[2] * x[eadj[0]] 
              +  mtxent[5] * x[eadj[1]] + mtxent[8] * x[eadj[2]] ;
}
/*
   -----------------------------
   change to packed storage mode
   -----------------------------
*/
InpMtx_changeStorageMode(inpmtx, INPMTX_SORTED) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n InpMtx after matrices loaded") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n right hand side vector via elemental mvm") ;
   DVfprintf(msgFile, nvtx, y) ;
}
/*
   ---------------------------------------------------
   the matrix as constructed is singular,
   uncomment this code to make node 0 a dirichlet node
   and so make the matrix nonsingular
   ---------------------------------------------------
*/
/*
{
double   *dvec ;
int      ii, nent ;
int      *ivec1, *ivec2 ;

nent = inpmtx->nent ;
ivec1 = InpMtx_ivec1(inpmtx) ;
ivec2 = InpMtx_ivec2(inpmtx) ;
dvec  = InpMtx_dvec(inpmtx) ;
for ( ii = 0 ; ii < nent ; ii++ ) {
   if ( ivec1[ii] == 0 ) {
      if ( ivec2[ii] == 0 ) {
         dvec[ii] = 1.0 ;
      } else {
         dvec[ii] = 0.0 ;
      }
   }
}
}
*/
/*
   -------------------------------------------------
   optionally write out the InpMtx object to a file
   -------------------------------------------------
*/
if ( strcmp(outFileName, "none") != 0 ) {
   InpMtx_writeToFile(inpmtx, outFileName) ;
}
/*
   ---------------------
   change to vector mode
   ---------------------
*/
InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n InpMtx in vector mode") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
}
/*
   ---------------------
   compute y[] = A * x[]
   ---------------------
*/
Y2 = DenseMtx_new() ;
DenseMtx_init(Y2, SPOOLES_REAL, 0, 0, nvtx, 1, 1, nvtx) ;
DenseMtx_zero(Y2) ;
MARKTIME(t1) ;
InpMtx_nonsym_mmm(inpmtx, Y2, alpha, X) ;
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n CPU %8.3f : compute matrix-vector multiply, %.3f mflops",
        t2 - t1, 2*inpmtx->nent*1.e-6/(t2 - t1)) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n right hand side vector via assembled mvm") ;
   DenseMtx_writeForHumanEye(Y2, msgFile) ;
}
/*
   ----------------
   compare Y and Y2
   ----------------
*/
DenseMtx_sub(Y2, Y) ;
fprintf(msgFile, "\n ||error||_max = %12.4e", DenseMtx_maxabs(Y2)) ;
/*
   --------------------------------------
   get random row and column permutations
   --------------------------------------
*/
rowOldToNewIV = IV_new() ;
IV_init(rowOldToNewIV, nvtx, NULL) ;
IV_ramp(rowOldToNewIV, 0, 1) ;
IV_shuffle(rowOldToNewIV, seed + 1) ;
rowOldToNew = IV_entries(rowOldToNewIV) ;
colOldToNewIV = IV_new() ;
IV_init(colOldToNewIV, nvtx, NULL) ;
IV_ramp(colOldToNewIV, 0, 1) ;
IV_shuffle(colOldToNewIV, seed + 2) ;
colOldToNew = IV_entries(colOldToNewIV) ;
/*
   ----------------------------------
   permute the InpMtx object, X and Y
   ----------------------------------
*/
MARKTIME(t1) ;
InpMtx_permute(inpmtx, rowOldToNew, colOldToNew) ;
DenseMtx_permuteRows(X, colOldToNewIV) ;
DenseMtx_permuteRows(Y, rowOldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute matrix", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n permuted matrix") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   fprintf(msgFile, "\n permuted solution") ;
   DenseMtx_writeForHumanEye(X, msgFile) ;
   fprintf(msgFile, "\n permuted right hand side") ;
   DenseMtx_writeForHumanEye(Y, msgFile) ;
}
/*
   -------------------------------------------------------
   compute the matrix-vector product with permuted vectors
   -------------------------------------------------------
*/
DenseMtx_zero(Y2) ;
MARKTIME(t1) ;
InpMtx_nonsym_mmm(inpmtx, Y2, alpha, X) ;
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n CPU %8.3f : compute matrix-vector multiply, %.3f mflops",
        t2 - t1, 2*inpmtx->nent*1.e-6/(t2 - t1)) ;
DenseMtx_sub(Y2, Y) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n error") ;
   DenseMtx_writeForHumanEye(Y2, msgFile) ;
}
fprintf(msgFile, "\n ||error||_max = %12.4e", DenseMtx_maxabs(Y2)) ;
/*
   ----------------
   free the objects
   ----------------
*/
EGraph_free(egraph) ;
Coords_free(coords) ;
InpMtx_free(inpmtx) ;
Drand_free(drand) ;
DenseMtx_free(X) ;
DenseMtx_free(Y) ;
DenseMtx_free(Y2) ;
DVfree(mtxent) ;
IV_free(rowOldToNewIV) ;
IV_free(colOldToNewIV) ;

fprintf(msgFile, "\n") ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------
   generate a matrix
   -----------------
*/
void
generateMatrix (
   int      eadj[],
   Coords   *coords,
   double   mtxent[],
   int      msglvl,
   FILE     *msgFile
) {
double   area, dx10, dx20, dx21, dy10, dy20, dy21 ;
double   x0, x1, x2, y0, y1, y2 ;
int      i, v0, v1, v2 ;

v0 = eadj[0] ;
v1 = eadj[1] ;
v2 = eadj[2] ;
x0 = Coords_value(coords, 1, v0) ;
x1 = Coords_value(coords, 1, v1) ;
x2 = Coords_value(coords, 1, v2) ;
y0 = Coords_value(coords, 2, v0) ;
y1 = Coords_value(coords, 2, v1) ;
y2 = Coords_value(coords, 2, v2) ;
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n vertex %d at (%6.3f, %6.3f)", v0, x0, y0) ;
   fprintf(msgFile, "\n vertex %d at (%6.3f, %6.3f)", v1, x1, y1) ;
   fprintf(msgFile, "\n vertex %d at (%6.3f, %6.3f)", v2, x2, y2) ;
   fflush(msgFile) ;
}
dx10 = x1 - x0 ;
dx20 = x2 - x0 ;
dx21 = x2 - x1 ;
dy10 = y1 - y0 ;
dy20 = y2 - y0 ;
dy21 = y2 - y1 ;
area = 0.5*(dx10*dy20 - dx20*dy10) ;
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n area = %20.12e", area) ;
   fflush(msgFile) ;
}
if ( area > 1.e-12 ) {
   mtxent[0] =   dx21*dx21 + dy21*dy21 ;
   mtxent[1] = - dx21*dx20 - dy21*dy20 ;
   mtxent[2] =   dx21*dx10 + dy21*dy10 ;
   mtxent[3] =   mtxent[1] ;
   mtxent[4] =   dx20*dx20 + dy20*dy20 ;
   mtxent[5] = - dx20*dx10 - dy20*dy10 ;
   mtxent[6] =   mtxent[2] ;
   mtxent[7] =   mtxent[5] ;
   mtxent[8] =   dx10*dx10 + dy10*dy10 ;
   for ( i = 0 ; i < 9 ; i++ ) {
      mtxent[i] /= (4.*area) ;
   }
   if ( msglvl > 3 && msgFile != NULL ) {
      fprintf(msgFile, "\n matrix "
              "\n [ %20.12e %20.12e %20.12e ] "
              "\n [ %20.12e %20.12e %20.12e ] "
              "\n [ %20.12e %20.12e %20.12e ] ",
              mtxent[0], mtxent[3], mtxent[6],
              mtxent[1], mtxent[4], mtxent[7],
              mtxent[2], mtxent[5], mtxent[8]) ;
      fprintf(msgFile, "\n rowsums = [ %13.5e %13.5e %13.5e ]",
              mtxent[0] + mtxent[3] + mtxent[6],
              mtxent[1] + mtxent[4] + mtxent[7],
              mtxent[2] + mtxent[5] + mtxent[8]) ;
      fprintf(msgFile, "\n loading matrix into bag") ;
      fflush(msgFile) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   load the matrices into the bag
   ------------------------------
*/
void
loadMatrices (
   EGraph   *egraph,
   Coords   *coords,
   InpMtx  *inpmtx,
   int      msglvl,
   FILE     *msgFile
) {
double   area, dx10, dx20, dx21, dy10, dy20, dy21 ;
double   mtxent[9] ;
double   x0, x1, x2, y0, y1, y2 ;
int      esize, i, ielem, nelem, nvtx, v0, v1, v2 ;
int      *eadj ;
IVL      *adj ;

nvtx  = egraph->nvtx   ;
nelem = egraph->nelem  ;
adj   = egraph->adjIVL ;
/*
for ( ielem = 0 ; ielem < 10 ; ielem++ ) {
*/
for ( ielem = 0 ; ielem < nelem ; ielem++ ) {
   IVL_listAndSize(adj, ielem, &esize, &eadj) ;
   v0 = eadj[0] ;
   v1 = eadj[1] ;
   v2 = eadj[2] ;
   x0 = Coords_value(coords, 1, v0) ;
   x1 = Coords_value(coords, 1, v1) ;
   x2 = Coords_value(coords, 1, v2) ;
   y0 = Coords_value(coords, 2, v0) ;
   y1 = Coords_value(coords, 2, v1) ;
   y2 = Coords_value(coords, 2, v2) ;
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n vertex %d at (%6.3f, %6.3f)", v0, x0, y0) ;
      fprintf(msgFile, "\n vertex %d at (%6.3f, %6.3f)", v1, x1, y1) ;
      fprintf(msgFile, "\n vertex %d at (%6.3f, %6.3f)", v2, x2, y2) ;
      fflush(msgFile) ;
   }
   dx10 = x1 - x0 ;
   dx20 = x2 - x0 ;
   dx21 = x2 - x1 ;
   dy10 = y1 - y0 ;
   dy20 = y2 - y0 ;
   dy21 = y2 - y1 ;
   area = 0.5*(dx10*dy20 - dx20*dy10) ;
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n area = %20.12e", area) ;
      fflush(msgFile) ;
   }
   if ( area > 1.e-12 ) {
      mtxent[0] =   dx21*dx21 + dy21*dy21 ;
      mtxent[1] = - dx21*dx20 - dy21*dy20 ;
      mtxent[2] =   dx21*dx10 + dy21*dy10 ;
      mtxent[3] =   mtxent[1] ;
      mtxent[4] =   dx20*dx20 + dy20*dy20 ;
      mtxent[5] = - dx20*dx10 - dy20*dy10 ;
      mtxent[6] =   mtxent[2] ;
      mtxent[7] =   mtxent[5] ;
      mtxent[8] =   dx10*dx10 + dy10*dy10 ;
      for ( i = 0 ; i < 9 ; i++ ) {
         mtxent[i] /= (4.*area) ;
      }
      if ( msglvl > 3 && msgFile != NULL ) {
         fprintf(msgFile, "\n matrix %d"
                 "\n [ %20.12e %20.12e %20.12e ] "
                 "\n [ %20.12e %20.12e %20.12e ] "
                 "\n [ %20.12e %20.12e %20.12e ] ",
                 ielem,
                 mtxent[0], mtxent[3], mtxent[6],
                 mtxent[1], mtxent[4], mtxent[7],
                 mtxent[2], mtxent[5], mtxent[8]) ;
         fprintf(msgFile, "\n rowsums = [ %13.5e %13.5e %13.5e ]",
                 mtxent[0] + mtxent[3] + mtxent[6],
                 mtxent[1] + mtxent[4] + mtxent[7],
                 mtxent[2] + mtxent[5] + mtxent[8]) ;
         fprintf(msgFile, "\n loading matrix into bag") ;
         fflush(msgFile) ;
      }
      InpMtx_inputRealMatrix(inpmtx, 3, 3, 1, 3, eadj, eadj, mtxent) ;
      if ( msglvl > 3 && msgFile != NULL ) {
         InpMtx_writeForHumanEye(inpmtx, stdout) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
