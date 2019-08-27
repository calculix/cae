/*  testSolve.c  */

#include "../../timings.h"
#include "../ILUMtx.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------
   test the ILUMtx solve methods

   created -- 98sep04, cca
   -----------------------------
*/
{
double   ops, t1, t2 ;
Drand    *drand ;
DV       *rhsDV, *solDV ;
int      LstorageMode, msglvl, neqns, rc, seed, symmetryflag, type,
         UstorageMode ;
ILUMtx   *mtx ;
FILE     *matlabFile, *msgFile ;

if ( argc != 10 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile neqns type symmetryflag "
      "\n         LstorageMode UstorageMode seed matlabFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    neqns    -- number of equations"
      "\n    type     -- type of entries"
      "\n       1 -- real entries"
      "\n       2 -- complex entries"
      "\n    symmetryflag -- type of symmetry"
      "\n       0 -- symmetric"
      "\n       1 -- hermitian"
      "\n       2 -- nonsymmetric"
      "\n    LstorageMode -- type of storage for L"
      "\n       1 -- by rows"
      "\n       2 -- by columns"
      "\n    UstorageMode -- type of storage for U"
      "\n       1 -- by rows"
      "\n       2 -- by columns"
      "\n    seed -- random number seed"
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
neqns        = atoi(argv[3]) ;
type         = atoi(argv[4]) ;
symmetryflag = atof(argv[5]) ;
LstorageMode = atof(argv[6]) ;
UstorageMode = atof(argv[7]) ;
seed         = atoi(argv[8]) ;
if ( strcmp(argv[9], "stdout") == 0 ) {
   matlabFile = stdout ;
} else if ( strcmp(argv[2], argv[9]) == 0 ) {
   matlabFile = msgFile ;
} else if ( strcmp("none", argv[9]) != 0 ) {
   if ( (matlabFile = fopen(argv[9], "w")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], argv[9]) ;
      return(-1) ;
   }
} else {
   matlabFile = NULL ;
}
fprintf(msgFile, 
        "\n %% %s "
        "\n %% msglvl       -- %d" 
        "\n %% msgFile      -- %s" 
        "\n %% neqns        -- %d" 
        "\n %% type         -- %d" 
        "\n %% symmetryflag -- %d" 
        "\n %% LstorageMode -- %d" 
        "\n %% UstorageMode -- %d" 
        "\n %% seed         -- %d" 
        "\n %% matlabFile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], neqns, type, symmetryflag,
        LstorageMode, UstorageMode, seed, argv[9]) ;
fflush(msgFile) ;
/*
   ----------------------------
   initialize the ILUMtx object
   ----------------------------
*/
mtx = ILUMtx_new() ;
rc = ILUMtx_init(mtx, neqns, type, symmetryflag, 
                 LstorageMode, UstorageMode) ;
if ( rc != 1 ) { return(-1) ; }
/*
   --------------------------------------
   fill with random structure and entries
   --------------------------------------
*/
rc = ILUMtx_fillRandom(mtx, seed) ;
if ( rc != 1 ) { return(-1) ; }
if ( matlabFile != NULL ) {
/*
   ----------------------
   write to a matlab file
   ----------------------
*/
   rc = ILUMtx_writeForMatlab(mtx, "L", "D", "U", matlabFile) ;
   if ( rc != 1 ) { return(-1) ; }
}
/*
   -------------------------------
   generate a random rhs DV object
   -------------------------------
*/
drand = Drand_new() ;
Drand_setUniform(drand, 0, 1) ;
Drand_setSeed(drand, seed + 1) ;
rhsDV = DV_new() ;
if ( type == SPOOLES_REAL ) {
   DV_init(rhsDV, neqns, NULL) ;
   Drand_fillDvector(drand, neqns, DV_entries(rhsDV)) ;
} else {
   DV_init(rhsDV, 2*neqns, NULL) ;
   Drand_fillDvector(drand, 2*neqns, DV_entries(rhsDV)) ;
}
if ( matlabFile != NULL ) {
   fprintf(matlabFile, "\n B = zeros(%d,1) ;", neqns) ;
   if ( type == SPOOLES_REAL ) {
      DV_writeForMatlab(rhsDV, "B", matlabFile) ;
   } else {
      rhsDV->size = neqns ;
      ZV_writeForMatlab((ZV *) rhsDV, "B", matlabFile) ;
      rhsDV->size = 2*neqns ;
   }
}
/*
   -------------------------------------------
   get a solution vector by solving the system
   -------------------------------------------
*/
solDV = DV_new() ;
if ( type == SPOOLES_REAL ) {
   DV_init(solDV, neqns, NULL) ;
   Drand_fillDvector(drand, neqns, DV_entries(solDV)) ;
} else {
   DV_init(solDV, 2*neqns, NULL) ;
}
DV_zero(solDV) ;
ops = 0.0 ;
MARKTIME(t1) ;
ILUMtx_solveVector(mtx, solDV, rhsDV, rhsDV, &ops, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n\n CPU %8.3f : solve system, %.0f ops, %.2f mflops\n",
        t2 - t1, ops, 1.e-6*ops/(t2 - t1)) ;
if ( matlabFile != NULL ) {
   fprintf(matlabFile, "\n X = zeros(%d,1) ;", neqns) ;
   if ( type == SPOOLES_REAL ) {
      DV_writeForMatlab(solDV, "X", matlabFile) ;
   } else {
      solDV->size = neqns ;
      ZV_writeForMatlab((ZV *) solDV, "X", matlabFile) ;
      solDV->size = 2*neqns ;
   }
/*
   ---------------------------------------------
   write out the matlab expression for the error
   ---------------------------------------------
*/
   fprintf(matlabFile, "\n error = max(abs((U\\(D\\(L\\B))) - X))") ;
   fflush(matlabFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
ILUMtx_free(mtx) ;
Drand_free(drand) ;
DV_free(rhsDV) ;
DV_free(solDV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
