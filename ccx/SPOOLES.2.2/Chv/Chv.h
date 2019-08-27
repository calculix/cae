/*  Chv.h  */

#include "../SubMtx.h"
#include "../PatchAndGoInfo.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   basic chevron object, double complex

   id      -- object's id
   nD      -- number of rows and columns in the diagonal (1,1) block
   nL      -- number of rows in the lower (2,1) block
   nU      -- number of columns in the upper (1,2) block
   type    -- type of chevron
      SPOOLES_REAL    --> chevron has real entries
      SPOOLES_COMPLEX --> chevron has complex entries
   symflag -- symmetry flag
      SPOOLES_SYMMETRIC    --> chevron is symmetric
      SPOOLES_HERMITIAN    --> chevron is hermitian
      SPOOLES_NONSYMMETRIC --> chevron is nonsymmetric
   rowind  -- pointer to the row indices
   colind  -- pointer to the column indices
   entries -- pointer to the entries
   wrkDV   -- working storage DV object
   next    -- next Chv in a singly linked list
   --------------------------------------------------------------------
*/
typedef struct _Chv   Chv ;
struct _Chv {
   int      id       ;
   int      nD       ;
   int      nL       ;
   int      nU       ;
   int      type     ;
   int      symflag  ;
   int      *rowind  ;
   int      *colind  ;
   double   *entries ;
   DV       wrkDV    ;
   Chv      *next    ;
} ;
#define CHV_IS_REAL(chv)         ((chv)->type    == SPOOLES_REAL)
#define CHV_IS_COMPLEX(chv)      ((chv)->type    == SPOOLES_COMPLEX)
#define CHV_IS_SYMMETRIC(chv)    ((chv)->symflag == SPOOLES_SYMMETRIC)
#define CHV_IS_HERMITIAN(chv)    ((chv)->symflag == SPOOLES_HERMITIAN)
#define CHV_IS_NONSYMMETRIC(chv) ((chv)->symflag == SPOOLES_NONSYMMETRIC)
/*
   -------------------------------------------------
   example of storage layout for indices and entries 

   nonsymmetric case, nD = 6, nL = 4, nU = 5
   +---------------------------------------+
   |      10 11 12 13 14 15 16 17 18 19 20 |
   |   +-----------------------------------+
   |   | +---------------------------------+
   | 9 | | 9 10 11 12 13 14 15 16 17 18 19 |
   | 8 | | 8 28 29 30 31 32 33 34 35 36 37 |
   | 7 | | 7 27 45 46 47 48 49 50 51 52 53 |
   | 6 | | 6 26 44 60 61 62 63 64 65 66 67 |
   | 5 | | 5 25 43 59 73 74 75 76 77 78 79 |
   | 4 | | 4 24 42 58 72 84 85 86 87 88 89 |
   | 3 | | 3 23 41 57 71 83 +--------------+
   | 2 | | 2 22 40 56 70 82 |
   | 1 | | 1 21 39 55 69 81 |
   | 0 | | 0 20 38 54 68 80 |
   +---+ +------------------+

   symmetric case, nD = 6, nU = 5
         +---------------------------------+
         | 0  1  2  3  4  5  6  7  8  9 10 |
         +---------------------------------+
         +---------------------------------+
         | 0  1  2  3  4  5  6  7  8  9 10 |
         |   11 12 13 14 15 16 17 18 19 20 |
         |      21 22 23 44 25 26 27 28 29 |
         |         30 31 32 33 34 35 36 37 |
         |            38 39 40 41 42 43 44 |
         |               45 46 47 48 49 50 |
         +---------------------------------+
   -------------------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98apr30, cca
   -----------------------
*/
Chv *
Chv_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98apr30, cca
   -----------------------
*/
void
Chv_setDefaultFields (
   Chv   *chv
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98apr30, cca
   --------------------------------------------------
*/
void
Chv_clearData (
   Chv   *chv
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98apr30, cca
   ------------------------------------------
*/
void
Chv_free (
   Chv   *chv
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------
   return the number of bytes needed to store the chevron
 
   created -- 98apr30, cca
   ------------------------------------------------------
*/
int
Chv_nbytesNeeded (
   int   nD,
   int   nL,
   int   nU,
   int   type,
   int   symflag
) ;
/*
   ----------------------------------------------------------------
   return the number of bytes in the workspace owned by this object
 
   created -- 98apr30, cca
   ----------------------------------------------------------------
*/
int
Chv_nbytesInWorkspace (
   Chv   *chv
) ;
/*
   ----------------------------------------------------------------
   set the number of bytes in the workspace owned by this object
 
   created -- 98apr30, cca
   ----------------------------------------------------------------
*/
void
Chv_setNbytesInWorkspace (
   Chv   *chv,
   int    nbytes
) ;
/*
   ----------------------------
   purpose -- set the fields
 
   created -- 98apr30, cca
   ----------------------------
*/
void
Chv_setFields (
   Chv     *chv,
   int      id,
   int      nD,
   int      nL,
   int      nU,
   int      type,
   int      symflag
) ;
/*
   ----------------------------
   purpose -- basic initializer
 
   created -- 98apr30, cca
   ----------------------------
*/
void
Chv_init (
   Chv     *chv,
   int      id,
   int      nD,
   int      nL,
   int      nU,
   int      type,
   int      symflag
) ;
/*
   ------------------------------------
   purpose -- initializer with pointers
 
   created -- 98apr30, cca
   ------------------------------------
*/
void
Chv_initWithPointers (
   Chv     *chv,
   int      id,
   int      nD,
   int      nL,
   int      nU,
   int      type,
   int      symflag,
   int      *rowind,
   int      *colind,
   double   *entries
) ;
/*
   -------------------------------------------------------------
   purpose -- to initialize the object from its working storage,
              used when the object is an MPI message
 
   created -- 98apr30
   -------------------------------------------------------------
*/
void
Chv_initFromBuffer (
   Chv   *chv
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------
   return the id of the chevron
 
   created -- 98apr30, cca
   ----------------------------
*/
int
Chv_id (
   Chv   *chv
) ;
/*
   ------------------------------------------------------------
   return the type of the chevron
   return value = SPOOLES_REAL    --> chevron is real
   return value = SPOOLES_COMPLEX --> chevron is complex
 
   created -- 98apr30, cca
   ------------------------------------------------------------
*/
int
Chv_type (
   Chv   *chv
) ;
/*
   ---------------------------------------------------------------
   return the symmetry flag of the chevron
   return value = SPOOLES_SYMMETRIC --> chevron is symmetric
   return value = SPOOLES_HERMITIAN --> chevron is hermitian
   return value = SPOOLES_NONSYMMETRIC --> chevron is nonsymmetric
 
   created -- 98apr30, cca
   ---------------------------------------------------------------
*/
int
Chv_symmetryFlag (
   Chv   *chv
) ;
/*
   --------------------------------------------------
   fill *pnD with nD, *pnL with nL, and *pnU with nU.
 
   created -- 98apr30, cca
   --------------------------------------------------
*/
void
Chv_dimensions (
   Chv   *chv,
   int    *pnD,
   int    *pnL,
   int    *pnU
) ;
/*
   ----------------------------------------------
   fill *pnrow with nD + nL, *prowind with rowind
 
   created -- 98apr30, cca
   ----------------------------------------------
*/
void
Chv_rowIndices (
   Chv   *chv,
   int    *pnrow,
   int    **prowind
) ;
/*
   ----------------------------------------------
   fill *pncol with nD + nU, *pcolind with colind
 
   created -- 98apr30, cca
   ----------------------------------------------
*/
void
Chv_columnIndices (
   Chv   *chv,
   int    *pncol,
   int    **pcolind
) ;
/*
   ----------------------------
   return the number of entries
 
   created -- 98apr30, cca
   ----------------------------
*/
int
Chv_nent (
   Chv   *chv
) ;
/*
   --------------------------------------------
   fill *pentries with a pointer to the entries
 
   created -- 98apr30, cca
   --------------------------------------------
*/
double *
Chv_entries(
   Chv   *chv
) ;
/*
   -----------------------------------------
   return the location of the diagonal entry
   for the ichv'th chevron
 
   created -- 98apr30, cca
   -----------------------------------------
*/
double *
Chv_diagLocation(
   Chv   *chv,
   int    ichv
) ;
/*
   ----------------------------------------------
   return a pointer to the start of the workspace
 
   created -- 98apr30, cca
   ----------------------------------------------
*/
void *
Chv_workspace(
   Chv   *chv
) ;
/*
   ------------------------------------
   fill *pValue with entry (irow, jcol)
 
   created -- 98apr30, cca
   ------------------------------------
*/
void
Chv_realEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   *pValue
) ;
/*
   --------------------------------------------
   fill (*pReal,*pImag) with entry (irow, jcol)
 
   created -- 98apr30, cca
   --------------------------------------------
*/
void
Chv_complexEntry (
   Chv     *chv,
   int      irow,
   int      jcol,
   double   *pReal,
   double   *pImag
) ;
/*
   -----------------------------------------------------
   fill *ppValue with the location of entry (irow, jcol)
 
   created -- 98apr30, cca
   -----------------------------------------------------
*/
void
Chv_locationOfRealEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   **ppValue
) ;
/*
   ----------------------------------------------------------
   fill (*ppReal,*ppImag) with location of entry (irow, jcol)
 
   created -- 98apr30, cca
   ----------------------------------------------------------
*/
void
Chv_locationOfComplexEntry (
   Chv     *chv,
   int      irow,
   int      jcol,
   double   **ppReal,
   double   **ppImag
) ;
/*
   ------------------------------------
   set entry (irow, jcol) to value
 
   created -- 98apr30, cca
   ------------------------------------
*/
void
Chv_setRealEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   value
) ;
/*
   --------------------------------------------
   fill (*pReal,*pImag) with entry (irow, jcol)
 
   created -- 98apr30, cca
   --------------------------------------------
*/
void
Chv_setComplexEntry (
   Chv     *chv,
   int      irow,
   int      jcol,
   double   real,
   double   imag
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in update.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------------------
   purpose --  perform the hermitian factor update
     T_{\bnd{I} \cap J, \bnd{I} \cap J}
          -= U_{I, \bnd{I} \cap J}^H D_{I, I} U_{I, \bnd{I} \cap J}
   and
     T_{\bnd{I} \cap J, \bnd{I} \cap \bnd{J}}
         -= U_{I, \bnd{I} \cap J}^H D_{I, I} U_{I, \bnd{I} \cap \bnd{J}}
 
   created -- 98apr17, cca
   ---------------------------------------------------------------------
*/
void
Chv_updateH (
   Chv      *chvT,
   SubMtx   *mtxD,
   SubMtx   *mtxU,
   DV       *tempDV
) ;
/*
   ---------------------------------------------------------------------
   purpose --  perform the symmetric factor update
     T_{\bnd{I} \cap J, \bnd{I} \cap J}
          -= U_{I, \bnd{I} \cap J}^T D_{I, I} U_{I, \bnd{I} \cap J}
   and
     T_{\bnd{I} \cap J, \bnd{I} \cap \bnd{J}}
         -= U_{I, \bnd{I} \cap J}^T D_{I, I} U_{I, \bnd{I} \cap \bnd{J}}
 
   created -- 98apr17, cca
   ---------------------------------------------------------------------
*/
void
Chv_updateS (
   Chv      *chvT,
   SubMtx   *mtxD,
   SubMtx   *mtxU,
   DV       *tempDV
) ;
/*
   ---------------------------------------------------------------------
   purpose --  perform the nonsymmetric factor update
     T_{\bnd{I} \cap J, \bnd{I} \cap J}
          -= L_{\bnd{I} \cap J, I} D_{I, I} U_{I, \bnd{I} \cap J}
   and
     T_{\bnd{I} \cap J, \bnd{I} \cap \bnd{J}}
         -= L_{\bnd{I} \cap J, I} D_{I, I} U_{I, \bnd{I} \cap \bnd{J}}
   and
     T_{\bnd{I} \cap \bnd{J}, \bnd{I} \cap J}
         -= L_{\bnd{I} \cap \bnd{J}, I} D_{I, I} U_{I, \bnd{I} \cap J}
 
   created -- 98feb27, cca
   ---------------------------------------------------------------------
*/
void
Chv_updateN (
   Chv   *chvT,
   SubMtx   *mtxL,
   SubMtx   *mtxD,
   SubMtx   *mtxU,
   DV     *tempDV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in factor.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   purpose -- to factor the front without pivoting
 
   return value -- # of eliminated rows and columns
 
   created -- 98aug27, cca
   ------------------------------------------------
*/
int
Chv_factorWithNoPivoting (
   Chv              *chv,
   PatchAndGoInfo   *info
) ;
/*
   ------------------------------------------------------------------
   purpose -- factor the pivot chevron with pivoting
 
   ndelay -- number of delayed rows and columns
   pivotflag -- enable pivoting or not
      0 --> no pivoting
      1 --> enable pivoting
   pivotsizesIV -- IV object that holds the sizes of the pivots,
      used only when the front is symmetric or hermitian
      and pivoting is enabled
   workDV -- DV object used for working storage, resized as necessary
   tau    -- upper bound on the magnitude of the entries 
      in the factors, used only when pivoting is enabled
   pntest -- pointer to be incremented with the number of pivot tests
 
   return value -- # of eliminated rows and columns
 
   created -- 98aug27, cca
   ------------------------------------------------------------------
*/
int
Chv_factorWithPivoting (
   Chv     *chv,
   int      ndelay,
   int      pivotflag,
   IV       *pivotsizesIV,
   DV       *workDV,
   double   tau,
   int      *pntest
) ;
/*
   ---------------------------------------------------------
   perform a rank one update using the first row and column.
   this is used in the (L + I)D(I + U) factorization
 
   return code ---
      0 if the pivot was zero
      1 if the pivot was nonzero
 
   created -- 98jan23, cca
   ---------------------------------------------------------
*/
int
Chv_r1upd (
   Chv   *chv
) ;
/*
   ------------------------------------------------------------------
   perform a rank two update using the first two rows.
   used in the (U^T + I)D(I + U) and (U^H + I)D(I + U) factorizations
 
   return code ---
      0 if the pivot was zero
      1 if the pivot was nonzero
 
   created -- 98jan23, cca
   ------------------------------------------------------------------
*/
int
Chv_r2upd (
   Chv   *chv
) ;
/*
   ------------------------------------------------------------------
   purpose -- looking at just a single chevron inside the Chv object,
              find the absolute value of the diagonal element, and
              the maximum absolute values of the offdiagonal elements 
              in the chevron's row and column.
 
   created -- 98aug26, cca
   ------------------------------------------------------------------
*/
void
Chv_maxabsInChevron (
   Chv      *chv,
   int      ichv,
   double   *pdiagmaxabs,
   double   *prowmaxabs,
   double   *pcolmaxabs
) ;
/*
   -------------------------------------------------------
   purpose -- zero the offdiagonal entries of chevron ichv
 
   created -- 98aug26, cca
   -------------------------------------------------------
*/
void
Chv_zeroOffdiagonalOfChevron (
   Chv   *chv,
   int   ichv
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in swap.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   swap rows irow and jrow
 
   created -- 98apr30, cca
   -----------------------
*/
void
Chv_swapRows (
   Chv   *chv,
   int    irow,
   int    jrow
) ;
/*
   --------------------------
   swap columns icol and jcol
 
   created -- 98apr30, cca
   --------------------------
*/
void
Chv_swapColumns (
   Chv   *chv,
   int     icol,
   int     jcol
) ;
/*
   -------------------------------
   swap rows and columns ii and jj
 
   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_swapRowsAndColumns (
   Chv   *chv,
   int    ii,
   int    jj
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in search.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------
   find the first unmarked entry in
   the diagonal with largest magnitude
   if ( mark[jj] == tag ) then
      we can compare this entry
   endif
 
   created -- 98apr30, cca
   -----------------------------------
*/
int
Chv_maxabsInDiagonal11 (
   Chv     *chv,
   int      mark[],
   int      tag,
   double   *pmaxval
) ;
/*
   --------------------------------------------
   find the first unmarked entry in
   row irow with largest magnitude
   if ( colmark[jj] == tag ) then
      we can examined this entry
   endif
   only entries in the (1,1) block are examined
 
   created -- 98apr30, cca
   --------------------------------------------
*/
int
Chv_maxabsInRow11 (
   Chv     *chv,
   int      irow,
   int      colmark[],
   int      tag,
   double   *pmaxval
) ;
/*
   --------------------------------------------
   find the first unmarked entry in
   column jcol with largest magnitude
   if ( rowmark[ii] == tag ) then
      we can examined this entry
   endif
   only entries in the (1,1) block are examined
 
   created -- 98apr30, cca
   --------------------------------------------
*/
int
Chv_maxabsInColumn11 (
   Chv     *chv,
   int      jcol,
   int      rowmark[],
   int      tag,
   double   *pmaxval
) ;
/*
   --------------------------------------
   return the location of the first entry
   with largest magnitude in row irow.
   *pmaxval is filled with its magnitude.
 
   created -- 98apr30, cca
   --------------------------------------
*/
int
Chv_maxabsInRow (
   Chv     *chv,
   int      irow,
   double   *pmaxval
) ;
/*
   --------------------------------------
   return the location of the first entry
   with largest magnitude in column jcol.
   *pmaxval is filled with its magnitude.
 
   created -- 98apr30, cca
   --------------------------------------
*/
int
Chv_maxabsInColumn (
   Chv     *chv,
   int      jcol,
   double   *pmaxval
) ;
/*
   -------------------------------------------------------------
   return the magnitude of a quasimax entry from the unmarked
   rows and columns and fill *pirow and *pjcol with its location
 
   created -- 98apr30, cca
   -------------------------------------------------------------
*/
double
Chv_quasimax (
   Chv     *chv,
   int      rowmark[],
   int      colmark[],
   int      tag,
   int      *pirow,
   int      *pjcol
) ;
/*
   ---------------------------------------------------------------
   find a 1x1 or 2x2 pivot using the fast Bunch-Parlett algorithm.
   used only with symmetric chevrons.
 
   created -- 98apr30, cca
   ---------------------------------------------------------------
*/
void
Chv_fastBunchParlettPivot (
   Chv    *chv,
   int     mark[],
   int     tag,
   int     *pirow,
   int     *pjcol
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in findPivot.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   purpose -- find and test a pivot
 
   workDV -- object that contains work vectors
   tau    -- upper bound on magnitude of factor entries
   ndelay -- number of delayed rows and columns on input
   pirow  -- pointer to be filled with pivot row
   pjcol  -- pointer to be filled with pivot column
   pntest -- pointer to be incremented with the number of pivot tests
 
   return value -- size of pivot
     0 --> pivot not found
     1 --> 1x1 pivot in row *pirow and column *pjcol
     2 --> 2x2 pivot in rows and columns *pirow and *pjcol,
           symmetric front only
   
   created -- 98jan24, cca
   ------------------------------------------------------------------
*/
int
Chv_findPivot (
   Chv      *chv,
   DV       *workDV,
   double   tau,
   int      ndelay,
   int      *pirow,
   int      *pjcol,
   int      *pntest
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in fill.c ------------------------------------------
------------------------------------------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in assemble.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------------
   add a scaled multiple of a simple chevron to a Chv object.
   the indices are offsets.
   note: for this purpose, (assembling original entries into the
   matrix), the row and column indices of the chevron are identical.
   also, the indices of both the Chv object and the chvind[]
   vector are assumed to be in ascending order.
 
   created -- 98apr30, cca
   -----------------------------------------------------------------
*/
void
Chv_addChevron (
   Chv      *chv,
   double   alpha[],
   int      ichv,
   int      chvsize,
   int      chvind[],
   double   chvent[]
) ;
/*
   --------------------------------------------------------------
   assemble Chv object chvI into Chv object chvJ.
   note: the two objects must be of the same symmetry type,
         the row indices of chvI must nest into those of chvJ,
         the column indices of chvI must nest into those of chvJ.
 
   created -- 98apr30, cca
   --------------------------------------------------------------
*/
void
Chv_assembleChv (
   Chv   *chvJ,
   Chv   *chvI
) ;
/*
   ----------------------------------------------------------------
   purpose -- assemble the postponed data from the children
 
 
   newchv     -- Chv object to contain fully assembled front
   oldchv     -- Chv object that contains former front
   firstchild -- pointer to first child in the list of children
                 Chv objects to be merged into the new front
 
   return value -- # of delayed rows and columns added to the front

   created -- 98apr30, cca
   ----------------------------------------------------------------
*/
int
Chv_assemblePostponedData (
   Chv   *newchv,
   Chv   *oldchv,
   Chv   *firstchild
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in copy.c ------------------------------------------
------------------------------------------------------------------------
*/
#define CHV_STRICT_LOWER    1
#define CHV_DIAGONAL        2
#define CHV_STRICT_UPPER    3
#define CHV_STRICT_LOWER_11 4
#define CHV_LOWER_21        5
#define CHV_STRICT_UPPER_11 6
#define CHV_UPPER_12        7
#define CHV_BY_ROWS         0
#define CHV_BY_COLUMNS      1
/*
   -------------------------------------------------------------------
  purpose -- copy entries to a vector.
 
   length     -- length of dvec[]
   npivot     -- number of pivots, may be 0
   pivotsizes -- vector of pivot sizes, may be NULL
   dvec[]     -- vector to receive matrix entries
   copyflag   -- flag to denote what part of the entries to copy
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_DIAGONAL        --> copy diagonal entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
   storeflag  -- flag to denote how to store entries in dvec[]
      CHV_BY_ROWS    --> store by rows
      CHV_BY_COLUMNS --> store by columns
 
   return value -- number of entries copied
 
   created  -- 97jun05, cca
   modified -- 98feb27, cca
     cases 4-7 inserted
   -------------------------------------------------------------------
*/
int
Chv_copyEntriesToVector (
   Chv      *chv,
   int      npivot,
   int      pivotsizes[],
   int      length,
   double   *dvec,
   int      copyflag,
   int      storeflag
) ;
/*
   -------------------------------------------------------------------
   purpose -- copy large entries to a vector. the portion copied
              can be a union of the strict lower portion,
              the diagonal portion, and the strict upper
              portion. there is one restriction, if the strict
              lower and strict upper are to be copied, the
              diagonal will also be copied.
 
   npivot     -- number of pivots, may be 0
   pivotsizes -- vector of pivot sizes, may be NULL
   sizes[]    -- vector to receive row/column sizes
   ivec[]     -- vector to receive row/column indices
   dvec[]     -- vector to receive matrix entries
   copyflag   -- flag to denote what part of the entries to copy
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
   storeflag  -- flag to denote how to store entries in dvec[]
      CHV_BY_ROWS    --> store by rows
      CHV_BY_COLUMNS --> store by columns
   droptol    -- entry to be copied must be larger than this magnitude
 
   return value -- number of entries copied
 
   created  -- 97jun05, cca
   modified -- 97feb27, cca
      cases 4-7 inserted
   -------------------------------------------------------------------
*/
int
Chv_copyBigEntriesToVector (
   Chv     *chv,
   int      npivot,
   int      pivotsizes[],
   int      sizes[],
   int      ivec[],
   double   dvec[],
   int      copyflag,
   int      storeflag,
   double   droptol
) ;
/*
   -------------------------------------------------------------------
   purpose -- return the number of entries
              in a portion of the object
 
   countflag -- which entries to count
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
 
   created -- 98feb27, cca
   -------------------------------------------------------------------
*/
int
Chv_countEntries (
   Chv     *chv,
   int      npivot,
   int      pivotsizes[],
   int      countflag
) ;
/*
   -------------------------------------------------------------------
   purpose -- return the number of entries
   whose magnitude is larger than droptol.
 
   countflag -- which entries to count
      CHV_STRICT_LOWER    --> copy strict lower entries
      CHV_STRICT_UPPER    --> copy strict upper entries
      CHV_STRICT_LOWER_11 --> copy strict lower entries in (1,1) block
      CHV_LOWER_21        --> copy lower entries in (2,1) block
      CHV_STRICT_UPPER_11 --> copy strict upper entries in (1,1) block
      CHV_UPPER_12        --> copy upper entries in (1,2) block
 
   created  -- 97jun07, cca
   modified -- 98feb27, cca
     cases 4-7 inserted
   -------------------------------------------------------------------
*/
int
Chv_countBigEntries (
   Chv     *chv,
   int      npivot,
   int      pivotsizes[],
   int      countflag,
   double   droptol
) ;
/*
   ----------------------------------------------------------
   purpose -- copy the trailing chevron that starts at offset
 
   created -- 97may16, cca
   ----------------------------------------------------------
*/
void
Chv_copyTrailingPortion (
   Chv   *chvI,
   Chv   *chvJ,
   int    offset
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------
   shift the indices, entries and adjust the nD dimension.
   note: shift can be positive or negative
 
   created -- 98apr30, cca
   -------------------------------------------------------
*/
void
Chv_shift (
   Chv   *chv,
   int    shift
) ;
/*
   ----------------------------------------------------------
   return the maximum magnitude of the entries in the chevron
 
   created -- 98apr30, cca
   ----------------------------------------------------------
*/
double
Chv_maxabs (
   Chv   *chv
) ;
/*
   -------------------------------------------------------
   return the frobenius norm of the entries in the chevron
 
   created -- 98apr30, cca
   -------------------------------------------------------
*/
double
Chv_frobNorm (
   Chv   *chv
) ;
/*
   -----------------------
   subtract chvI from chvJ
 
   created -- 98apr30, cca
   -----------------------
*/
void
Chv_sub (
   Chv   *chvJ,
   Chv   *chvI
) ;
/*
   -------------------------------
   zero the entries in the chevron
 
   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_zero (
   Chv   *chv
) ;
/*
   -------------------------------
   fill A2 object with (1,1) block
 
   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_fill11block (
   Chv   *chv,
   A2    *mtx
) ;
/*
   -------------------------------
   fill A2 object with (1,2) block
 
   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_fill12block (
   Chv   *chv,
   A2    *mtx
) ;
/*
   -------------------------------
   fill A2 object with (2,1) block
 
   created -- 98apr30, cca
   -------------------------------
*/
void
Chv_fill21block (
   Chv   *chv,
   A2    *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form
 
   created -- 98apr30, cca
   ----------------------------------------
*/
void
Chv_writeForHumanEye (
   Chv    *chv,
   FILE   *fp
) ;
/*
   ------------------------------------------------
   purpose -- write out the entries in matlab style
 
   created -- 98apr30, cca
   ------------------------------------------------
*/
void
Chv_writeForMatlab (
   Chv    *chv,
   char   *chvname,
   FILE   *fp
) ;
/*--------------------------------------------------------------------*/
