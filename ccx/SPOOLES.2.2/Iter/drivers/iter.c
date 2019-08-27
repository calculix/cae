/*  iter.c  */
#include "../Iter.h"
#include "../../Graph.h"
#include <string.h>
#define METHODS 10

int 
InpMtx_readFromAIJ2file(InpMtx *mtxA, char *srcFileName);
int
InpMtx_createGraph(InpMtx *mtxA, Graph    *graph);

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------------
   test the factor method for a grid matrix
   (0) read in matrix from source file 
   (1) conver data matrix to InpMtx object if necessary
   (2) create Graph and ETree object if necessary
   (3) read in/create an ETree object
   (4) create a solution matrix object
   (5) multiply the solution with the matrix
       to get a right hand side matrix object
   (6) factor the matrix 
   (7) solve the system

   created   -- 98dec30, jwu
   -----------------------------------------------------
*/
{
char            etreeFileName[80], mtxFileName[80], *cpt, rhsFileName[80],
                srcFileName[80], ctemp[81], msgFileName[80], slnFileName[80] ;
Chv             *chv, *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxB, *mtxQ, *mtxX, *mtxZ ;
double          one[2] = { 1.0, 0.0 } ;
FrontMtx        *frontmtx ;
InpMtx          *mtxA ;
SubMtxManager   *mtxmanager ;
double          cputotal, droptol, conv_tol, factorops ;
double          cpus[9] ;
Drand           drand ;
double          nops, tau, t1, t2   ;
ETree           *frontETree   ;
Graph           *graph ;
FILE            *msgFile, *inFile ;
int             error, loc, msglvl, neqns, nzf, iformat, 
                pivotingflag, rc, seed, sparsityflag, symmetryflag, 
                method[METHODS], type, nrhs, etreeflag ;
int             stats[6] ;
int             nnzA, Ik, itermax, zversion, iterout ;
IV              *newToOldIV, *oldToNewIV ;
IVL             *symbfacIVL ;
int             i, j, k, m, n, imethod, maxdomainsize, maxzeros, maxsize;
int             nouter,ninner ;

if ( argc != 2 ) {
   fprintf(stdout, 
"\n\n usage : %s inFile"
"\n    inFile       -- input filename"
"\n", argv[0]) ;
   return(-1) ;
}

/* read input file */
inFile = fopen(argv[1], "r");
if (inFile == (FILE *)NULL) {
  fprintf(stderr, "\n fatal error in %s: unable to open file %s\n",
           argv[0], argv[1]) ;
  return(-1) ;
}

for (i=0; i<METHODS; i++) method[i]=-1; 
imethod=0;
k=0;
while (1) {
  fgets(ctemp, 80, inFile);
  if (ctemp[0] != '*') {
    /*printf("l=%2d:%s\n", strlen(ctemp),ctemp);*/
    if (strlen(ctemp)==80) {
      fprintf(stderr, "\n fatal error in %s: input line contains more than "
	      "80 characters.\n",argv[0]);
      exit(-1);
    }
    if (k==0) {
      sscanf(ctemp, "%d",  &iformat);
      if (iformat < 0 || iformat > 2) {
	fprintf(stderr, "\n fatal error in %s: "
		"invalid source matrix format\n",argv[0]) ;
	return(-1) ;
      }
    }
    else if (k==1)
      sscanf(ctemp, "%s", srcFileName);
    else if (k==2)
      sscanf(ctemp, "%s", mtxFileName);
    else if (k==3) {
      sscanf(ctemp, "%d",  &etreeflag);
      if (etreeflag < 0 || etreeflag > 4) {
	fprintf(stderr, "\n fatal error in %s: "
                        "invalid etree file status\n",argv[0]) ;
	return(-1) ;
      }
    }
    else if (k==4)
      sscanf(ctemp, "%s", etreeFileName);
    else if (k==5)
      sscanf(ctemp, "%s", rhsFileName);
    else if (k==6)
      sscanf(ctemp, "%s", slnFileName);
    else if (k==7){
      sscanf(ctemp, "%s", msgFileName);
      if ( strcmp(msgFileName, "stdout") == 0 ) {
	msgFile = stdout ;
      }
      else if ( (msgFile = fopen(msgFileName, "a")) == NULL ) {
	fprintf(stderr, "\n fatal error in %s"
		"\n unable to open file %s\n", argv[0], ctemp) ;
	return(-1) ;
      }
    }
    else if (k==8)
      sscanf(ctemp, "%d %d %d %d %d %d", 
	     &msglvl, &seed, &nrhs, &Ik, &itermax, &iterout);
    else if (k==9)
      sscanf(ctemp, "%d %d %d", &symmetryflag, &sparsityflag, &pivotingflag);
    else if (k==10)
      sscanf(ctemp, "%lf %lf %lf", &tau, &droptol, &conv_tol);
    else if (k==11) {
      /*
      for (j=0; j<strlen(ctemp); j++) {
	printf("j=%2d:%s",j,ctemp+j);
	if (ctemp[j] == ' ' && ctemp[j+1] != ' ') {
	  sscanf(ctemp+j, "%d", method+imethod);
          printf("method[%d]=%d\n",imethod,method[imethod]);
	  if (method[imethod] < 0) break;
	  imethod++;
	}
      }
      */
      imethod = sscanf(ctemp,"%d %d %d %d %d %d %d %d %d %d",
		       method, method+1, method+2, method+3, method+4,
		       method+5, method+6, method+7, method+8, method+9);
      /*printf("imethod=%d\n",imethod);*/
      for (j=0; j<imethod; j++) {
	/*printf("method[%d]=%d\n",j,method[j]);*/
	if (method[j]<0) {
	  imethod=j;
          break;
	}
      }
      if (imethod == 0) {
	fprintf(msgFile,"No method assigned in input file\n");
	return(-1);
      }
    }
    k++;
  }
  if (k==12) break;
}

fclose(inFile);

/* reset nrhs to 1 */
if (nrhs > 1) {
  fprintf(msgFile,"*** Multiple right-hand-side vectors is not allowed yet.\n");
  fprintf(msgFile,"*** nrhs is reset to 1.\n");
  nrhs =1;
}

fprintf(msgFile, 
        "\n %s "
        "\n srcFileName   -- %s"
        "\n mtxFileName   -- %s"
        "\n etreeFileName -- %s"
        "\n rhsFileName   -- %s"
        "\n msglvl        -- %d" 
        "\n seed          -- %d" 
        "\n symmetryflag  -- %d" 
        "\n sparsityflag  -- %d" 
        "\n pivotingflag  -- %d" 
        "\n tau           -- %e" 
        "\n droptol       -- %e" 
        "\n conv_tol      -- %e"
        "\n method        -- ",
        argv[0], srcFileName, mtxFileName, etreeFileName, rhsFileName,
	msglvl, seed, symmetryflag, sparsityflag, pivotingflag, 
        tau, droptol, conv_tol) ;
 
for (k=0; k<imethod; k++) 
  fprintf(msgFile, "%d ", method[k]);
fprintf(msgFile, "\n ", method[k]);

fflush(msgFile) ;

/*
   --------------------------------------
   initialize the random number generator
   --------------------------------------
*/
Drand_setDefaultFields(&drand) ;
Drand_init(&drand) ;
Drand_setSeed(&drand, seed) ;
/*Drand_setUniform(&drand, 0.0, 1.0) ;*/
Drand_setNormal(&drand, 0.0, 1.0) ;
/*
   ----------------------------------------------
   read in or convert source to the InpMtx object
   ----------------------------------------------
*/
rc = 1;

if ( strcmp(srcFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(-1) ;
}
mtxA = InpMtx_new() ;

MARKTIME(t1) ;
if (iformat == 0)  { /* InpMtx source format */
  rc = InpMtx_readFromFile(mtxA, srcFileName) ;
  strcpy(mtxFileName, srcFileName);
  if ( rc != 1 ) 
    fprintf(msgFile, "\n return value %d from InpMtx_readFromFile(%p,%s)",
	    rc, mtxA, srcFileName) ;
}
else if (iformat == 1) {  /* HBF source format */
  rc = InpMtx_readFromHBfile(mtxA, srcFileName) ;
  if ( rc != 1 ) 
    fprintf(msgFile, "\n return value %d from InpMtx_readFromHBfile(%p,%s)",
	    rc, mtxA, srcFileName) ;
}
else { /* AIJ2 source format */
  rc = InpMtx_readFromAIJ2file(mtxA, srcFileName) ;
  if ( rc != 1 ) 
    fprintf(msgFile, "\n return value %d from InpMtx_readFromAIJ2file(%p,%s)",
	    rc, mtxA, srcFileName) ;
}
MARKTIME(t2) ;
if (iformat>0 && strcmp(mtxFileName, "none") != 0 ) {
  rc = InpMtx_writeToFile(mtxA, mtxFileName) ;
  if ( rc != 1 )
    fprintf(msgFile, "\n return value %d from InpMtx_writeToFile(%p,%s)",
	    rc, mtxA, mtxFileName) ;
}

fprintf(msgFile, "\n CPU %8.3f : read in (+ convert to) mtxA from file %s",
	t2 - t1, mtxFileName) ;
if (rc != 1) {
  goto end_read;
}
type = mtxA->inputMode ;
neqns = 1 + IVmax(mtxA->nent, InpMtx_ivec1(mtxA), &loc) ;
if ( INPMTX_IS_BY_ROWS(mtxA) ) {
  fprintf(msgFile, "\n matrix coordinate type is rows") ;
} else if ( INPMTX_IS_BY_COLUMNS(mtxA) ) {
  fprintf(msgFile, "\n matrix coordinate type is columns") ;
} else if ( INPMTX_IS_BY_CHEVRONS(mtxA) ) {
  fprintf(msgFile, "\n matrix coordinate type is chevrons") ;
} else {
  fprintf(msgFile, "\n\n, error, bad coordinate type") ;
  rc=-1;
  goto end_read;
}
if ( INPMTX_IS_RAW_DATA(mtxA) ) {
  fprintf(msgFile, "\n matrix storage mode is raw data\n") ;
} else if ( INPMTX_IS_SORTED(mtxA) ) {
  fprintf(msgFile, "\n matrix storage mode is sorted\n") ;
} else if ( INPMTX_IS_BY_VECTORS(mtxA) ) {
  fprintf(msgFile, "\n matrix storage mode is by vectors\n") ;
} else {
  fprintf(msgFile, "\n\n, error, bad storage mode") ;
  rc=-1;
  goto end_read;
}

if ( msglvl > 1 ) {
  fprintf(msgFile, "\n\n after reading InpMtx object from file %s",
	  mtxFileName) ;
  if ( msglvl == 2 ) {
    InpMtx_writeStats(mtxA, msgFile) ;
  } else {
    InpMtx_writeForHumanEye(mtxA, msgFile) ;
  }
  fflush(msgFile) ;
}
/*
  Get the nonzeros in matrix A and print it
  */
nnzA  = InpMtx_nent( mtxA );
fprintf(msgFile, "\n\n Input matrix size  %d NNZ  %d",
	neqns, nnzA) ;

/*
   --------------------------------------------------------
   generate the linear system
   1. generate solution matrix and fill with random numbers
   2. generate rhs matrix and fill with zeros
   3. compute matrix-matrix multiply
   --------------------------------------------------------
*/
MARKTIME(t1) ;
mtxX = DenseMtx_new() ;
DenseMtx_init(mtxX, type, 0, -1, neqns, nrhs, 1, neqns) ;
mtxB = DenseMtx_new() ; 

if (strcmp(rhsFileName, "none")) {
  rc = DenseMtx_readFromFile(mtxB, rhsFileName) ;
  if ( rc != 1 )
    fprintf(msgFile, "\n return value %d from DenseMtx_readFromFile(%p,%s)",
	    rc, mtxB, rhsFileName) ;
  DenseMtx_zero(mtxX) ;
}
else {
  DenseMtx_init(mtxB, type, 1, -1, neqns, nrhs, 1, neqns) ;
  DenseMtx_fillRandomEntries(mtxX, &drand) ;
  DenseMtx_zero(mtxB) ;
  switch ( symmetryflag ) {
  case SPOOLES_SYMMETRIC : 
    InpMtx_sym_mmm(mtxA, mtxB, one, mtxX) ;
    break ;
  case SPOOLES_HERMITIAN :
    InpMtx_herm_mmm(mtxA, mtxB, one, mtxX) ;
    break ;
  case SPOOLES_NONSYMMETRIC :
    InpMtx_nonsym_mmm(mtxA, mtxB, one, mtxX) ;
    break ;
  default :
    break ;
  }
}
  
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : set up the solution and rhs ",
        t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n original mtxX") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fprintf(msgFile, "\n\n original mtxB") ;
   DenseMtx_writeForHumanEye(mtxB, msgFile) ;
   fflush(msgFile) ;
}
if (rc != 1) {
  InpMtx_free(mtxA);
  DenseMtx_free(mtxX);
  DenseMtx_free(mtxB);
  goto end_init;
}

/*
   ------------------------
   read in/create the ETree object
   ------------------------
*/

MARKTIME(t1) ;
if (etreeflag == 0) { /* read in ETree from file */
  if ( strcmp(etreeFileName, "none") == 0 ) 
    fprintf(msgFile, "\n no file to read from") ;
  frontETree = ETree_new() ;
  rc = ETree_readFromFile(frontETree, etreeFileName) ;
  if (rc!=1) 
    fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
	    rc, frontETree, etreeFileName) ;
}
else {
  graph = Graph_new() ;
  rc = InpMtx_createGraph(mtxA, graph);
  if (rc!=1) {
    fprintf(msgFile, "\n return value %d from InpMtx_createGraph(%p,%p)",
	    rc, mtxA, graph) ;
    Graph_free(graph);
    goto end_tree;
  }
  if (etreeflag == 1) { /* Via BestOfNDandMS */
    maxdomainsize = 500; maxzeros      = 1000; maxsize       = 64    ;
    frontETree = orderViaBestOfNDandMS(graph, maxdomainsize, maxzeros,
				       maxsize, seed, msglvl, msgFile) ;
  }
  else if (etreeflag == 2) { /* Via MMD */
    frontETree = orderViaMMD(graph, seed, msglvl, msgFile) ;        
  }
  else if (etreeflag == 3) { /* Via MS */
    maxdomainsize = 500;
    frontETree = orderViaMS(graph, maxdomainsize, seed, msglvl, msgFile) ;
  }
  else if (etreeflag == 4) { /* Via ND */
    maxdomainsize = 500;
    frontETree = orderViaND(graph, maxdomainsize, seed, msglvl, msgFile) ;
  }
  Graph_free(graph);

  /*    optionally write out the ETree object    */
  if ( strcmp(etreeFileName, "none") != 0 ) {
    fprintf(msgFile, "\n\n writing out ETree to file %s", 
	    etreeFileName) ;
    ETree_writeToFile(frontETree, etreeFileName) ;
  }
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in/create frontETree from file %s",
	t2 - t1, etreeFileName) ;
if ( rc != 1 ) {
  ETree_free(frontETree);
  goto end_tree;
}

ETree_leftJustify(frontETree) ;
if ( msglvl > 1 ) {
  fprintf(msgFile, "\n\n after reading ETree object from file %s",
	  etreeFileName) ;
  if ( msglvl == 2 ) {
    ETree_writeStats(frontETree, msgFile) ;
  } else {
    ETree_writeForHumanEye(frontETree, msgFile) ;
  }
}
fflush(msgFile) ;
/*
   --------------------------------------------------
   get the permutations, permute the matrix and the 
   front tree, and compute the symbolic factorization
   --------------------------------------------------
*/
MARKTIME(t1) ;
oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : get permutations", t2 - t1) ;
MARKTIME(t1) ;
ETree_permuteVertices(frontETree, oldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute front tree", t2 - t1) ;
MARKTIME(t1) ;
InpMtx_permute(mtxA, IV_entries(oldToNewIV), IV_entries(oldToNewIV)) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute mtxA", t2 - t1) ;
if (  symmetryflag == SPOOLES_SYMMETRIC
   || symmetryflag == SPOOLES_HERMITIAN ) {
   MARKTIME(t1) ;
   InpMtx_mapToUpperTriangle(mtxA) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : map to upper triangle", t2 - t1) ;
}
if ( ! INPMTX_IS_BY_CHEVRONS(mtxA) ) {
   MARKTIME(t1) ;
   InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : change coordinate type", t2 - t1) ;
}
if ( INPMTX_IS_RAW_DATA(mtxA) ) {
   MARKTIME(t1) ;
   InpMtx_changeStorageMode(mtxA, INPMTX_SORTED) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : sort entries ", t2 - t1) ;
}
if ( INPMTX_IS_SORTED(mtxA) ) {
   MARKTIME(t1) ;
   InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : convert to vectors ", t2 - t1) ;
}
MARKTIME(t1) ;
symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : symbolic factorization", t2 - t1) ;
MARKTIME(t1) ;
DenseMtx_permuteRows(mtxB, oldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : permute rhs", t2 - t1) ;

/*
   ------------------------------
   initialize the FrontMtx object
   ------------------------------
*/
MARKTIME(t1) ;
frontmtx   = FrontMtx_new() ;
mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(mtxmanager, NO_LOCK, 0) ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL,
              type, symmetryflag, sparsityflag, pivotingflag,
              NO_LOCK, 0, NULL, mtxmanager, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : initialize the front matrix",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile,
           "\n nendD  = %d, nentL = %d, nentU = %d",
           frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;
   SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n front matrix initialized") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------
   factor the matrix
   -----------------
*/
nzf       = ETree_nFactorEntries(frontETree, symmetryflag) ;
factorops = ETree_nFactorOps(frontETree, type, symmetryflag) ;
fprintf(msgFile, 
        "\n %d factor entries, %.0f factor ops, %8.3f ratio",
        nzf, factorops, factorops/nzf) ;
IVzero(6, stats) ;
DVzero(9, cpus) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, NO_LOCK, 1) ;
MARKTIME(t1) ;
rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, droptol, 
                                chvmanager, &error, cpus, 
                                stats, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : factor matrix, %8.3f mflops",
        t2 - t1, 1.e-6*factorops/(t2-t1)) ;
if ( rootchv != NULL ) {
   fprintf(msgFile, "\n\n factorization did not complete") ;
   for ( chv = rootchv ; chv != NULL ; chv = chv->next ) {
      fprintf(stdout, "\n chv %d, nD = %d, nL = %d, nU = %d",
              chv->id, chv->nD, chv->nL, chv->nU) ;
   }
}
if ( error >= 0 ) {
   fprintf(msgFile, "\n\n error encountered at front %d\n", error) ;
   rc=error ;
   goto end_front;
}
fprintf(msgFile,
        "\n %8d pivots, %8d pivot tests, %8d delayed rows and columns",
        stats[0], stats[1], stats[2]) ;
if ( frontmtx->rowadjIVL != NULL ) {
   fprintf(msgFile,
           "\n %d entries in rowadjIVL", frontmtx->rowadjIVL->tsize) ;
}
if ( frontmtx->coladjIVL != NULL ) {
   fprintf(msgFile,
           ", %d entries in coladjIVL", frontmtx->coladjIVL->tsize) ;
}
if ( frontmtx->upperblockIVL != NULL ) {
   fprintf(msgFile,
           "\n %d fronts, %d entries in upperblockIVL", 
           frontmtx->nfront, frontmtx->upperblockIVL->tsize) ;
}
if ( frontmtx->lowerblockIVL != NULL ) {
   fprintf(msgFile,
           ", %d entries in lowerblockIVL", 
           frontmtx->lowerblockIVL->tsize) ;
}
fprintf(msgFile,
        "\n %d entries in D, %d entries in L, %d entries in U",
        stats[3], stats[4], stats[5]) ;
fprintf(msgFile, "\n %d locks", frontmtx->nlocks) ;
if (  FRONTMTX_IS_SYMMETRIC(frontmtx)
   || FRONTMTX_IS_HERMITIAN(frontmtx) ) {
   int   nneg, npos, nzero ;

   FrontMtx_inertia(frontmtx, &nneg, &nzero, &npos) ;
   fprintf(msgFile, 
           "\n %d negative, %d zero and %d positive eigenvalues",
           nneg, nzero, npos) ;
   fflush(msgFile) ;
}
cputotal = cpus[8] ;
if ( cputotal > 0.0 ) {
   fprintf(msgFile,
   "\n    initialize fronts       %8.3f %6.2f"
   "\n    load original entries   %8.3f %6.2f"
   "\n    update fronts           %8.3f %6.2f"
   "\n    assemble postponed data %8.3f %6.2f"
   "\n    factor fronts           %8.3f %6.2f"
   "\n    extract postponed data  %8.3f %6.2f"
   "\n    store factor entries    %8.3f %6.2f"
   "\n    miscellaneous           %8.3f %6.2f"
   "\n    total time              %8.3f",
   cpus[0], 100.*cpus[0]/cputotal,
   cpus[1], 100.*cpus[1]/cputotal,
   cpus[2], 100.*cpus[2]/cputotal,
   cpus[3], 100.*cpus[3]/cputotal,
   cpus[4], 100.*cpus[4]/cputotal,
   cpus[5], 100.*cpus[5]/cputotal,
   cpus[6], 100.*cpus[6]/cputotal,
   cpus[7], 100.*cpus[7]/cputotal, cputotal) ;
}
if ( msglvl > 1 ) {
  SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
  ChvManager_writeForHumanEye(chvmanager, msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front factor matrix") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
}

/*
   ------------------------------
   post-process the factor matrix
   ------------------------------
*/
MARKTIME(t1) ;
FrontMtx_postProcess(frontmtx, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %8.3f : post-process the matrix", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front factor matrix after post-processing") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
}
fprintf(msgFile, "\n\n after post-processing") ;
if ( msglvl > 1 ) SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
/*
   ----------------
   solve the system
   ----------------
*/
neqns = mtxB->nrow ;
mtxZ  = DenseMtx_new() ;
DenseMtx_init(mtxZ, type, 0, 0, neqns, nrhs, 1, neqns) ;
zversion=INPMTX_IS_COMPLEX_ENTRIES(mtxA);

for (k=0; k<imethod; k++) {
  DenseMtx_zero(mtxZ) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n rhs") ;
    DenseMtx_writeForHumanEye(mtxB, msgFile) ;
    fflush(stdout) ;
  }
  fprintf(msgFile, "\n\n itemax  %d", itermax) ;
  DVzero(6, cpus) ;
  MARKTIME(t1) ;
  switch ( method[k] ) {
  case BiCGStabR :
    if (zversion)
      rc=zbicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		    itermax, conv_tol, msglvl, msgFile);
    else
      rc=bicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		   itermax, conv_tol, msglvl, msgFile);

    break;
  case BiCGStabL :
    if (zversion)
    rc=zbicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		  itermax, conv_tol, msglvl, msgFile);
    else
      rc=bicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		   itermax, conv_tol, msglvl, msgFile);
    break;
  case TFQMRR :
    if (zversion)
      rc=ztfqmrr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		 itermax, conv_tol, msglvl, msgFile);
    else
      rc=tfqmrr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		itermax, conv_tol, msglvl, msgFile);
    break;
  case TFQMRL :
    if (zversion)
      rc=ztfqmrl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		 itermax, conv_tol, msglvl, msgFile);
    else
      rc=tfqmrl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
		itermax, conv_tol, msglvl, msgFile);
    break;
  case PCGR :
    if (zversion)
      rc=zpcgr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
	       itermax, conv_tol, msglvl, msgFile);
    else
      rc=pcgr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
	      itermax, conv_tol, msglvl, msgFile);
    break;
  case PCGL :
    if (zversion)
      rc=zpcgl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
	       itermax, conv_tol, msglvl, msgFile);
    else
      rc=pcgl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
	      itermax, conv_tol, msglvl, msgFile);
    break;
  case MLBiCGStabR :
    mtxQ = DenseMtx_new() ;
    DenseMtx_init(mtxQ, type, 0, -1, neqns, Ik, 1, neqns) ;
    Drand_setUniform(&drand, 0.0, 1.0) ;
    DenseMtx_fillRandomEntries(mtxQ, &drand) ;
    if (zversion)
      rc=zmlbicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
		      mtxB, itermax, conv_tol, msglvl, msgFile);
    else
      rc=mlbicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
		     mtxB, itermax, conv_tol, msglvl, msgFile);
    DenseMtx_free(mtxQ) ;
    break;
  case MLBiCGStabL :
    mtxQ = DenseMtx_new() ;
    DenseMtx_init(mtxQ, type, 0, -1, neqns, Ik, 1, neqns) ;
    Drand_setUniform(&drand, 0.0, 1.0) ;
    DenseMtx_fillRandomEntries(mtxQ, &drand) ;
    if (zversion)
      rc=zmlbicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
		      mtxB, itermax, conv_tol, msglvl, msgFile);
    else
      rc=mlbicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
		     mtxB, itermax, conv_tol, msglvl, msgFile);
    DenseMtx_free(mtxQ) ;
    break;
  case BGMRESR:    
    if (zversion)
      fprintf(msgFile, "\n\n *** BGMRESR complex version is not available "
	      "at this moment.   ") ;
    else
      rc=bgmresr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ,
                 mtxB, iterout, itermax, &nouter, &ninner, conv_tol,
                 msglvl, msgFile);
    break;
  case BGMRESL:    
    if (zversion)
      fprintf(msgFile, "\n\n *** BGMRESR complex version is not available "
	      "at this moment.   ") ;
    else
      rc=bgmresl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ,
                 mtxB, iterout, itermax, &nouter, &ninner, conv_tol,
                 msglvl, msgFile);
    break;
  default:
    fprintf(msgFile, "\n\n *** Invalid method number   ") ;
  }
  
  MARKTIME(t2) ;
  fprintf(msgFile, "\n\n CPU %8.3f : solve the system", t2 - t1) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n computed solution") ;
    DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
    fflush(stdout) ;
  }
  
/*
  -------------------------------------------------------------
  permute the computed solution back into the original ordering
  -------------------------------------------------------------
*/
  MARKTIME(t1) ;
  DenseMtx_permuteRows(mtxZ, newToOldIV) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : permute solution", t2 - t1) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n permuted solution") ;
    DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
    fflush(stdout) ;
  }
/*
  -------------
  save solution
  -------------
*/
  if (  strcmp(slnFileName, "none") != 0 ) {
    DenseMtx_writeToFile(mtxZ, slnFileName) ;
  }
/*
  -----------------
  compute the error
  -----------------
*/
  if (!strcmp(rhsFileName, "none")) {    
    DenseMtx_sub(mtxZ, mtxX) ;
    if (method[k] <8) {
      mtxQ = DenseMtx_new() ;
      DenseMtx_init(mtxQ, type, 0, -1, neqns, 1, 1, neqns) ;
      rc=DenseMtx_initAsSubmatrix (mtxQ, mtxZ, 0, neqns-1, 0, 0);
      fprintf(msgFile, "\n\n maxabs error = %12.4e", DenseMtx_maxabs(mtxQ)) ;
      DenseMtx_free(mtxQ) ;
    }
    else
      fprintf(msgFile, "\n\n maxabs error = %12.4e", DenseMtx_maxabs(mtxZ)) ;

    if ( msglvl > 1 ) {
      fprintf(msgFile, "\n\n error") ;
      DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
      fflush(stdout) ;
    }
    if ( msglvl > 1 ) 
      SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
  }
  fprintf(msgFile, "\n---------  End of Method %d -------\n",method[k]) ;
      
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DenseMtx_free(mtxZ) ;

end_front:
ChvManager_free(chvmanager) ;
SubMtxManager_free(mtxmanager) ;
FrontMtx_free(frontmtx) ;
IVL_free(symbfacIVL) ;
IV_free(oldToNewIV) ;
IV_free(newToOldIV) ;

end_tree:
ETree_free(frontETree) ;

end_init:
DenseMtx_free(mtxB) ;
DenseMtx_free(mtxX) ;

end_read:
InpMtx_free(mtxA) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   InpMtx_readFromAIJ2file reads in an AIJ2 matrix from file and
   converts to InpMtx object. 
 
   created -- 98dec01, jwu
   -------------------------------------------------------------------
*/
int 
InpMtx_readFromAIJ2file(InpMtx *mtxA, char *srcFileName)
{

FILE  *srcFile;
int   rc, nrow, ncol, irow, icol, jrow, jcol, ient, nentInRow;
double ent;
/*
   ----------------------------
   open the input file and read
   #rows #columns 
   ----------------------------
*/
if ( (srcFile = fopen(srcFileName, "r")) == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_readFromAIJ2file"
           "\n unable to open file %s\n",
           srcFileName) ;
   return(-1) ;
}
rc = fscanf(srcFile, "%d %d", &nrow, &ncol) ;
if ( rc != 2 ) {
   fprintf(stderr, "\n fatal error in InpMtx_readFromAIJ2file"
           "\n %d of 2 fields read on first line of file %s",
           rc, srcFileName) ;
   return(-1) ;
}

/*
   --------------------------------------------------
   initialize the object
   set coordType = INPMTX_BY_ROWS --> row coordinates
   --------------------------------------------------
*/
InpMtx_init(mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, 10*nrow, 0) ;
/*
   -------------------------------------------------
   read in the entries and load them into the object
   -------------------------------------------------
*/
ient=0;
for ( irow = 0 ; irow < nrow ; irow++ ) {
   rc = fscanf(srcFile, "%d %d", &jrow, &nentInRow) ;
   if ( rc != 2 ) {
      fprintf(stderr, "\n fatal error in InpMtx_readFromAIJ2file"
              "\n %d of 2 fields read on entry %d of file %s",
              rc, ient, srcFileName) ;
      return(-1) ;
   }
   for ( icol = 0 ; icol < nentInRow ; icol++ ) {
      rc = fscanf(srcFile, "%d %lf", &jcol, &ent) ;
      if ( rc != 2 ) {
         fprintf(stderr, "\n fatal error in InpMtx_readFromAIJ2file"
                 "\n %d of 2 fields read on entry %d of file %s",
                 rc, ient, srcFileName) ;
         return(-1) ;
      }
      InpMtx_inputRealEntry(mtxA, jrow-1, jcol-1, ent) ;
      ient++;
   }
}
/*
   -----------------------------
   sort and compress the entries
   -----------------------------
*/
InpMtx_changeStorageMode(mtxA, 3) ;

fclose(srcFile);
return(1) ; }

/*--------------------------------------------------------------------*/

/*  createGraph.c  */
/*--------------------------------------------------------------------*/
int
InpMtx_createGraph (InpMtx *mtxA, Graph    *graph)
/*
   ----------------------------------------------------
   read in a InpMtx object and create the Graph object

   created -- 97feb14, cca
   ----------------------------------------------------
*/
{
int      count, msglvl, nvtx, rc ;
IVL      *adjIVL ;

InpMtx_changeStorageMode(mtxA, 3) ;
nvtx  = 1 + IV_max(&mtxA->ivec1IV) ;
count = 1 + IV_max(&mtxA->ivec2IV) ;
if ( nvtx < count ) {
   nvtx = count ;
}
/*
   ------------------------------------
   create the full adjacency IVL object
   ------------------------------------
*/
adjIVL = InpMtx_fullAdjacency(mtxA) ;
/*
   ---------------------
   fill the Graph object
   ---------------------
*/
Graph_init2(graph, 0, nvtx, 0, adjIVL->tsize, nvtx, adjIVL->tsize, 
            adjIVL, NULL, NULL) ;


return(1) ; }

/*--------------------------------------------------------------------*/


