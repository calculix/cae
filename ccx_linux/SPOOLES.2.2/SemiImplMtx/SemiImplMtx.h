/*  SemiImplMtx.h */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   the SemiImplMtx object is used to perform "semi-implicit solves"
   of the following matrix factorization.

   A = P [ A11 A12 ] Q = P [ L11  0  ] [ D11  0  ] [ U11 U12 ] Q
         [ A21 A22 ]       [ L21 L22 ] [  0  D22 ] [  0  U22 ]

   where P and Q are chosen by pivoting for stability.

   the semi-implicit factorization drops L21 and U12, and trades
   operating with them with multiplies by A21 and A12 and extra
   solves involving A11.

   we want to solve [ A11 A12 ] [ X1 ] = [ B1 ]
                    [ A21 A22 ] [ X2 ] = [ B2 ]
   this can be done explicitly, as
      1. solve L11 Y1 = B1
      2. solve L22 Y2 = B2 - L21 Y1
      3. solve D11 Z1 = Y1
      4. solve D22 Z2 = Y2
      5. solve U22 X2 = Z2
      6. solve U11 X1 = Z1 - U12 X2
   or implicitly, as
      1. solve L11 D11 U11 T  = B1
      2. solve L22 D22 U22 X2 = B2 - A21 T
      3. solve L11 D11 U11 T  = B1 - A12 X2

   neqns -- # of equations
   type  -- type of entries
      1 (SPOOLES_REAL)    -- real entries
      2 (SPOOLES_COMPLEX) -- complex entries
   symmetryflag -- symmetry type
      0 (SPOOLES_SYMMETRIC)    -- symmetric matrix
      1 (SPOOLES_HERMITIAN)    -- hermitian matrix
      2 (SPOOLES_NONSYMMETRIC) -- nonsymmetric matrix
   ndomeqns    -- # of equations in the domains
   nschureqns  -- # of equations in the schur complement
   domainMtx   -- pointer to FrontMtx objec that holds 
                  the factorization of A11
   schurMtx    -- pointer to FrontMtx objec that holds 
                  the factorization of A22 - A21 * A11^{-1} * A12
   domRowsIV   -- vector of global ids of rows of A11
   schurRowsIV -- vector of global ids of rows of A22
   domColsIV   -- vector of global ids of columns of A11
   schurColsIV -- vector of global ids of columns of A22

   created -- 98sep14
   ----------------------------------------------------------------
*/
typedef struct _SemiImplMtx   SemiImplMtx ;
struct _SemiImplMtx {
   int        neqns        ;
   int        type         ;
   int        symmetryflag ;
   int        ndomeqns     ;
   int        nschureqns   ;
   FrontMtx   *domainMtx   ;
   FrontMtx   *schurMtx    ;
   InpMtx     *A21         ;
   InpMtx     *A12         ;
   IV         *domRowsIV   ;
   IV         *schurRowsIV ;
   IV         *domColsIV   ;
   IV         *schurColsIV ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98oct16, cca
   -----------------------
*/
SemiImplMtx *
SemiImplMtx_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   return code --
      1 -- normal return
     -1 -- mtx is NULL
 
   created -- 98oct16, cca
   -----------------------
*/
int
SemiImplMtx_setDefaultFields (
   SemiImplMtx   *mtx
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   return code --
      1 -- normal return
     -1 -- mtx is NULL
 
   created -- 98oct16, cca
   --------------------------------------------------
*/
int
SemiImplMtx_clearData (
   SemiImplMtx   *mtx
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   return code --
      1 -- normal return
     -1 -- mtx is NULL
 
   created -- 98oct16, cca
   ------------------------------------------
*/
int
SemiImplMtx_free (
   SemiImplMtx   *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- to initialize the semi-implicit matrix using as input a
              FrontMtx and a map from fronts to domains (map[J] != 0)
              or the schur complement (map[J] = 0)
 
   return value --
      1 -- normal return
     -1 -- semimtx is NULL
     -2 -- frontmtx is NULL
     -3 -- inpmtx is NULL
     -4 -- frontmapIV is NULL
     -5 -- frontmapIV is invalid
 
   created -- 98oct17, cca
   ------------------------------------------------------------------
*/
int
SemiImplMtx_initFromFrontMtx (
   SemiImplMtx   *semimtx,
   FrontMtx      *frontmtx,
   InpMtx        *inpmtx,
   IV            *frontmapIV,
   int           msglvl,
   FILE          *msgFile
) ;
/*
*/
int
FrontMtx_initFromSubmatrix (
   FrontMtx   *submtx,
   FrontMtx   *frontmtx,
   IV         *frontidsIV,
   IV         *rowsIV,
   IV         *colsIV,
   int        msglvl,
   FILE       *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solve.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   purpose -- to solve the linear system A X = B 
              using a semi-implicit factorization
 
   on return ---
     cpus[0] -- time to initialize working matrices
     cpus[1] -- time to load rhs 
     cpus[2] -- time for first solve with domains
     cpus[3] -- time to compute schur rhs
     cpus[4] -- time for schur solve
     cpus[5] -- time to compute domains' rhs
     cpus[6] -- time for second solve with domains
     cpus[7] -- time to store solution
     cpus[8] -- miscellaneous time 
     cpus[9] -- total time 
 
   return value --
      1 -- normal return
     -1 -- semimtx is NULL
     -2 -- mtxX is NULL
     -3 -- mtxB is NULL
     -4 -- mtxmanager is NULL
     -5 -- cpus is NULL
 
   created -- 98oct17, cca
   ------------------------------------------------
*/
int
SemiImplMtx_solve (
   SemiImplMtx     *semimtx,
   DenseMtx        *mtxX,
   DenseMtx        *mtxB,
   SubMtxManager   *mtxmanager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------
   fill a statistics array for a semi-implicit factorization
     stats[0]  -- # of equations
     stats[1]  -- # of equations in the (1,1) block
     stats[2]  -- # of equations in the (2,2) block
     stats[3]  -- # of entries in L11
     stats[4]  -- # of entries in D11
     stats[5]  -- # of entries in U11
     stats[6]  -- # of entries in L22
     stats[7]  -- # of entries in D22
     stats[8]  -- # of entries in U22
     stats[9]  -- # of entries in A12
     stats[10] -- # of entries in A21
     stats[11] -- total # of matrix entries
     stats[12] -- # of operations for a solve
 
   return value ---
      1 -- normal return
     -1 -- semimtx is NULL
     -2 -- stats is NULL
 
   created -- 98oct22, cca
   ---------------------------------------------------------
*/
int
SemiImplMtx_stats (
   SemiImplMtx   *semimtx,
   int           stats[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------
   purpose -- to write a SemiImplMtx to a file
              in a human readable format
 
   return values ---
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- type is invalid
     -3 -- symmetry flag is invalid
     -4 -- fp is NULL
 
   created -- 98oct16, cca
   -------------------------------------------
*/
int
SemiImplMtx_writeForHumanEye (
   SemiImplMtx   *mtx,
   FILE          *fp
) ;
/*--------------------------------------------------------------------*/
