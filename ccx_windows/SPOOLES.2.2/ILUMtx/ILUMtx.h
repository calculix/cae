/*  ILUMtx.h  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   the ILUMtx is used to compute, store and solve with a factorization 
     of the form A = (L+I)D(I+U), (U^T+I)D(I+U) or (U^H+I)D(I+U).
     the factorization is very simple, stored by vectors, and can
     be exact (up to roundoff) or approximate. it is a less
     complicated alternative to the FrontMtx object.

   neqns        - # of equations
   type         - SPOOLES_INDICES_ONLY, SPOOLES_REAL or SPOOLES_COMPLEX
   symmetryflag - SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN 
      or SPOOLES_COMPLEX
   UstorageMode -- storage mode for U
      SPOOLES_BY_ROWS or SPOOLES_BY_COLUMNS
   LstorageMode -- storage mode for L
      SPOOLES_BY_ROWS or SPOOLES_BY_COLUMNS
   sizesL -- vector of sizes of the L vectors
   sizesU -- vector of sizes of the U vectors
   p_indL -- vector of pointers to indices of the L vectors
   p_indU -- vector of pointers to indices of the U vectors
   entD   -- vector of entries in D
   p_entL -- vector of pointers to entries of the L vectors
   p_entU -- vector of pointers to entries of the U vectors

   created -- 98oct03, cca
   --------------------------------------------------------------------
*/
typedef struct _ILUMtx   ILUMtx ;
struct _ILUMtx {
  int      neqns        ;
  int      type         ;
  int      symmetryflag ;
  int      UstorageMode ;
  int      LstorageMode ;
  int      *sizesL      ;
  int      *sizesU      ;
  int      **p_indL     ;
  int      **p_indU     ;
  double   *entD        ;
  double   **p_entL     ;
  double   **p_entU     ;
} ;
/*--------------------------------------------------------------------*/
/*
   ------------
   handy macros
   ------------
*/
#define ILUMTX_IS_REAL(mtx)     ((mtx)->type == SPOOLES_REAL)
#define ILUMTX_IS_COMPLEX(mtx)  ((mtx)->type == SPOOLES_COMPLEX)

#define ILUMTX_IS_SYMMETRIC(mtx)    \
   ((mtx)->symmetryflag == SPOOLES_SYMMETRIC)
#define ILUMTX_IS_HERMITIAN(mtx)    \
   ((mtx)->symmetryflag == SPOOLES_HERMITIAN)
#define ILUMTX_IS_NONSYMMETRIC(mtx) \
   ((mtx)->symmetryflag == SPOOLES_NONSYMMETRIC)

#define ILUMTX_IS_L_BY_ROWS(mtx)     \
   ((mtx)->LstorageMode == SPOOLES_BY_ROWS)
#define ILUMTX_IS_L_BY_COLUMNS(mtx)     \
   ((mtx)->LstorageMode == SPOOLES_BY_COLUMNS)
#define ILUMTX_IS_U_BY_ROWS(mtx)     \
   ((mtx)->UstorageMode == SPOOLES_BY_ROWS)
#define ILUMTX_IS_U_BY_COLUMNS(mtx)     \
   ((mtx)->UstorageMode == SPOOLES_BY_COLUMNS)

/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98oct03, cca
   -----------------------
*/
ILUMtx *
ILUMtx_new (
   void
) ;
/*
   -----------------------
   set the default fields
 
   return code --
      1 -- normal return
     -1 -- mtx is NULL
 
   created -- 98oct03, cca
   -----------------------
*/
int
ILUMtx_setDefaultFields (
   ILUMtx   *mtx
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   return code --
      1 -- normal return
     -1 -- mtx is NULL
 
   created -- 98oct03, cca
   --------------------------------------------------
*/
int
ILUMtx_clearData (
   ILUMtx   *mtx
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   return code --
      1 -- normal return
     -1 -- mtx is NULL
 
   created -- 98oct03, cca
   ------------------------------------------
*/
int
ILUMtx_free (
   ILUMtx   *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------
   purpose -- initialize the ILUMtx object
 
   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtx
     -4  -- bad symmetryflag for mtx
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- matrix is symmetric or hermitian
            and storage modes are not compatible
 
   created -- 98oct03, cca
   ---------------------------------------------
*/
int
ILUMtx_init (
   ILUMtx   *mtx,
   int      neqns,
   int      type,
   int      symmetryflag,
   int      LstorageMode,
   int      UstorageMode
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in factor.c ----------------------------------------
------------------------------------------------------------------------
*/
/*

-------------------------------------------------------------------
  purpose -- to factor the linear system
      A = (L + I)D(I + U), A = (U^T + I)D(I + U) or
      A = (U^H + I)D(I + U). 
 
   if pops is not NULL, then on return *pops has been incremented
   by the number of operations performed in the solve.
 
   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtxLDU
     -4  -- bad symmetryflag for mtxLDU
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL
     -14 -- mtxA is NULL
     -15 -- types of mtxLDU and mtxA are not the same
     -16 -- mtxA is not in chevron mode
     -17 -- droptol < 0.0
     -18 -- msglvl > 0 and msgFile is NULL
     -19 -- singular pivot found
 
   created -- 98oct03, cca

-------------------------------------------------------------------
*/
int
ILUMtx_factor (
   ILUMtx   *mtxLDU,
   InpMtx   *mtxA,
   double   droptol,
   double   *pops,
   int      msglvl,
   FILE     *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solve.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------------
   purpose -- to solve the linear system
      (L + I)D(I + U) X = B, (U^T + I)D(I + U) X = B or
      (U^H + I)D(I + U) X = B. X and B are single vectors.
 
   workDV is a temporary vector. 
   note, if workDV is different than B, then B is unchanged on return.
   one can have X, B and workDV all point to the same object.
 
   if pops is not NULL, then on return *pops has been incremented
   by the number of operations performed in the solve.
 
   return values ---
      1  -- normal return
     -1  -- LDU is NULL
     -2  -- neqns <= 0
     -3  -- bad type for LDU
     -4  -- bad symmetryflag for LDU
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL
     -14 -- X is NULL
     -15 -- size of X is incorrect
     -16 -- entries of X are NULL
     -17 -- B is NULL
     -18 -- size of B is incorrect
     -19 -- entries of B are NULL
     -20 -- workDV is NULL
     -21 -- size of workDV != neqns
     -22 -- entries of workDV are NULL
     -23 -- msglvl > 0 and msgFile is NULL
 
   created -- 98oct03, cca
   -------------------------------------------------------------------
*/
int
ILUMtx_solveVector (
   ILUMtx   *LDU,
   DV       *X,
   DV       *B,
   DV       *workDV,
   double   *pops,
   int      msglvl,
   FILE     *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in misc.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   purpose -- to fill the indices and entries with random values
 
   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtx
     -4  -- bad symmetryflag for mtx
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL
 
   created -- 98oct03, cca
   -------------------------------------------------------------
*/
int
ILUMtx_fillRandom (
   ILUMtx   *mtx,
   int      seed
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------
   purpose -- to write an ILUMtx object to a file in matlab format
 
   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for LDU
     -4  -- bad symmetryflag for LDU
     -5  -- bad LstorageMode for LDU
     -6  -- bad UstorageMode for LDU
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL
     -14 -- Lname is NULL
     -15 -- Dname is NULL
     -16 -- Uname is NULL
     -17 -- fp is NULL
 
   created -- 98oct03, cca
   ---------------------------------------------------------------
*/
int
ILUMtx_writeForMatlab (
   ILUMtx   *mtx,
   char     *Lname,
   char     *Dname,
   char     *Uname,
   FILE     *fp
) ;
/*--------------------------------------------------------------------*/
