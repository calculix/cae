/*  DenseMtx   */

#include "../cfiles.h"
#include "../Drand.h"
#include "../A2.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   this object handles a dense real and complex matrix.

   type  -- type of entries
      1 -- real
      2 -- complex
   rowid -- id for the rows in this matrix
   colid -- id for the columns in this matrix
   nrow  -- # of rows in the matrix
   ncol  -- # of columns in the matrix
   inc1  -- row increment for the entries
   inc2  -- column increment for the entries
   rowind  -- pointer to row indices
   colind  -- pointer to column indices
   entries -- pointer to matrix entries
   wrkDV   -- DV object to manage workspace
   next    -- pointer to next DenseMtx object in a singly linked list
   -------------------------------------------------------------------
*/
typedef struct _DenseMtx   DenseMtx ;
struct _DenseMtx {
   int        type     ;
   int        rowid    ;
   int        colid    ;
   int        nrow     ;
   int        ncol     ;
   int        inc1     ;
   int        inc2     ;
   int        *rowind  ;
   int        *colind  ;
   double     *entries ;
   DV         wrkDV    ;
   DenseMtx   *next    ;
} ;

#define DENSEMTX_IS_REAL(mtx)    ((mtx)->type == SPOOLES_REAL)
#define DENSEMTX_IS_COMPLEX(mtx) ((mtx)->type == SPOOLES_COMPLEX)
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98may02, cca
   -----------------------
*/
DenseMtx *
DenseMtx_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may02, cca
   -----------------------
*/
void
DenseMtx_setDefaultFields (
   DenseMtx   *mtx
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may02, cca
   --------------------------------------------------
*/
void
DenseMtx_clearData (
   DenseMtx   *mtx
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may02, cca
   ------------------------------------------
*/
void
DenseMtx_free (
   DenseMtx   *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------
   return the row id of the object
 
   created -- 98may02, cca
   -------------------------------
*/
int
DenseMtx_rowid (
   DenseMtx   *mtx
) ;
/*
   ----------------------------------
   return the column id of the object
 
   created -- 98may02, cca
   ----------------------------------
*/
int
DenseMtx_colid (
   DenseMtx   *mtx
) ;
/*
   ------------------------------------------
   fill *pnrow with nrow and *pncol with ncol
 
   created -- 98may02, cca
   ------------------------------------------
*/
void
DenseMtx_dimensions (
   DenseMtx   *mtx,
   int         *pnrow,
   int         *pncol
) ;
/*
   ------------------------
   return the row increment
 
   created -- 98may02, cca
   ------------------------
*/
int
DenseMtx_rowIncrement (
   DenseMtx   *mtx
) ;
/*
   ---------------------------
   return the column increment
 
   created -- 98may02, cca
   ---------------------------
*/
int
DenseMtx_columnIncrement (
   DenseMtx   *mtx
) ;
/*
   -------------------------------------------
   fill *pnrow with nrow, *prowind with rowind
 
   created -- 98may02, cca
   -------------------------------------------
*/
void
DenseMtx_rowIndices (
   DenseMtx   *mtx,
   int         *pnrow,
   int         **prowind
) ;
/*
   -------------------------------------------
   fill *pncol with ncol, *pcolind with colind
 
   created -- 98may02, cca
   -------------------------------------------
*/
void
DenseMtx_columnIndices (
   DenseMtx   *mtx,
   int         *pncol,
   int         **pcolind
) ;
/*
   --------------------------------------------
   fill *pentries with a pointer to the entries
 
   created -- 98may02, cca
   --------------------------------------------
*/
double *
DenseMtx_entries(
   DenseMtx   *mtx
) ;
/*
   ---------------------------------
   return a pointer to the workspace
 
   created -- 98may02, cca
   ---------------------------------
*/
void *
DenseMtx_workspace(
   DenseMtx   *mtx
) ;
/*
   ------------------------------------------
   fill *pValue with the entry in (irow,jcol)
 
   created -- 98jun05, cca
   ------------------------------------------
*/
void
DenseMtx_realEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     *pValue
) ;
/*
   ----------------------------------------------------
   fill *pReal and *pImag with the entry in (irow,jcol)
 
   created -- 98jun05, cca
   ----------------------------------------------------
*/
void
DenseMtx_complexEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     *pReal,
   double     *pImag
) ;
/*
   ------------------------------
   set entry (irow,jcol) to value
 
   created -- 98jun05, cca
   ------------------------------
*/
void
DenseMtx_setRealEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     value
) ;
/*
   ------------------------------------
   set entry (irow,jcol) to (real,imag)
 
   created -- 98jun05, cca
   ------------------------------------
*/
void
DenseMtx_setComplexEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     real,
   double     imag
) ;
/*
   ----------------------------------------------------------
   purpose -- fill *prowent with the base address of row irow

   return values --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- invalid type for mtx
     -3 -- irow is invalid
     -4 -- prowent is NULL

   created -- 98nov11, cca
   ----------------------------------------------------------
*/
int
DenseMtx_row (
   DenseMtx   *mtx,
   int        irow,
   double     **prowent
) ;
/*
   -------------------------------------------------------------
   purpose -- fill *pcolent with the base address of column jcol

   return values --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- invalid type for mtx
     -3 -- jcol is invalid
     -4 -- pcolent is NULL

   created -- 98nov11, cca
   -------------------------------------------------------------
*/
int
DenseMtx_column (
   DenseMtx   *mtx,
   int        jcol,
   double     **pcolent
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------
   return the number of bytes needed for an object of this size
 
   created -- 98may02, cca
   ------------------------------------------------------------
*/
int
DenseMtx_nbytesNeeded (
   int   type,
   int   nrow,
   int   ncol
) ;
/*
   ----------------------------------------------------------------
  return the number of bytes in the workspace owned by this object
 
   created -- 98may02, cca
   ----------------------------------------------------------------
*/
int
DenseMtx_nbytesInWorkspace (
   DenseMtx   *mtx
) ;
/*
   -------------------------------------------------------------
   set the number of bytes in the workspace owned by this object
 
   created -- 98may02, cca
   -------------------------------------------------------------
*/
void
DenseMtx_setNbytesInWorkspace (
   DenseMtx   *mtx,
   int        nbytes
) ;
/*
   ---------------------------------------
   purpose -- set the fields of the object
 
   created -- 98may02, cca
   ---------------------------------------
*/
void
DenseMtx_setFields (
   DenseMtx   *mtx,
   int        type,
   int        rowid,
   int        colid,
   int        nrow,
   int        ncol,
   int        inc1,
   int        inc2
) ;
/*
   ----------------------------
   purpose -- basic initializer
 
   created -- 98may02, cca
   ----------------------------
*/
void
DenseMtx_init (
   DenseMtx   *mtx,
   int        type,
   int        rowid,
   int        colid,
   int        nrow,
   int        ncol,
   int        inc1,
   int        inc2
) ;
/*
   ---------------------------------------------------------
   purpose -- initialize the object from its working storage
              used when the object is a MPI message
 
   created -- 98may02, cca
   ---------------------------------------------------------
*/
void
DenseMtx_initFromBuffer (
   DenseMtx   *mtx
) ;
/*
   ------------------------------------
   purpose -- initializer with pointers
 
   created -- 98may02, cca
   ------------------------------------
*/
void
DenseMtx_initWithPointers (
   DenseMtx   *mtx,
   int        type,
   int        rowid,
   int        colid,
   int        nrow,
   int        ncol,
   int        inc1,
   int        inc2,
   int        *rowind,
   int        *colind,
   double     *entries
) ;
/*
   -----------------------------------
   this method initializes a A2 object
   to point into the entries
 
   created -- 98may02, cca
   -----------------------------------
*/
void
DenseMtx_setA2 (
   DenseMtx   *mtx,
   A2         *a2
) ;
/*
   ----------------------------------------------------------------
  purpose -- initialize as a submatrix of another DenseMtx object.
         B = A(firstrow:lastrow, firstcol:lastcol)
      note, B only points into the storage of A.

   return values --
      1 -- normal return
     -1 -- B is NULL
     -2 -- A is NULL
     -3 -- A has invalid type
     -4 -- requested rows are invalid
     -5 -- requested columns are invalid

   created -- 98nov11, cca
   ----------------------------------------------------------------
*/
int
DenseMtx_initAsSubmatrix (
   DenseMtx   *B,
   DenseMtx   *A,
   int        firstrow,
   int        lastrow,
   int        firstcol,
   int        lastcol
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- to read in the object from a file
 
   input --
 
      fn -- filename, must be *.densemtxb or *.densemtxf
 
   return value -- 1 if success, 0 if failure
 
   created -- 98aug13
   -----------------------------------------------------------
*/
int
DenseMtx_readFromFile (
   DenseMtx   *mtx,
   char       *fn
) ;
/*
   --------------------------------------------------
   purpose -- to read an object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 98aug13
   --------------------------------------------------
*/
int
DenseMtx_readFromFormattedFile (
   DenseMtx   *mtx,
   FILE       *fp
) ;
/*
   -----------------------------------------------
   purpose -- to read an object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 98aug13
   -----------------------------------------------
*/
int
DenseMtx_readFromBinaryFile (
   DenseMtx   *mtx,
   FILE       *fp
) ;
/*
   ---------------------------------------
   purpose -- to write an object to a file
 
   input --
 
      fn -- filename
        *.densemtxb -- binary
        *.densemtxf -- formatted
        anything else -- for human eye
 
   error return --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fn is NULL
     -3 -- unable to open file
 
   created -- 98aug13
   ---------------------------------------
*/
int
DenseMtx_writeToFile ( 
   DenseMtx   *mtx, 
   char       *fn 
) ;
/*
   -------------------------------------------------
   purpose -- to write an object to a formatted file
 
   return value --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fn is NULL
 
   created -- 98aug13
   -------------------------------------------------
*/
int
DenseMtx_writeToFormattedFile (
   DenseMtx   *mtx,
   FILE       *fp
) ;
/*
   --------------------------------------------------------
   purpose -- to write an adjacency object to a binary file
 
   return value --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fn is NULL
 
   created -- 98aug13
   --------------------------------------------------------
*/
int
DenseMtx_writeToBinaryFile (
   DenseMtx   *mtx,
   FILE       *fp
) ;
/*
   -----------------------------------------------------
   purpose -- to write the object's statistics to a file
              in human readable form
 
   return value --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fp is NULL
 
   created -- 98may02, cca
   -----------------------------------------------------
*/
int
DenseMtx_writeStats (
   DenseMtx   *mtx,
   FILE       *fp
) ;
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form
 
   return value --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fp is NULL
 
   created -- 98may02, cca
   ----------------------------------------
*/
int
DenseMtx_writeForHumanEye (
   DenseMtx   *mtx,
   FILE       *fp
) ;
/*
   -----------------------------------------------
   purpose -- to write the object to a matlab file
 
   return value --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- mtx is NULL
     -3 -- fp is NULL
 
   created -- 98may02, cca
   -----------------------------------------------
*/
int
DenseMtx_writeForMatlab (
   DenseMtx   *mtx,
   char       *mtxname,
   FILE       *fp
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in permute.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------
   purpose -- to permute the rows of an object
 
   created -- 98may02, cca
   -------------------------------------------
*/
void
DenseMtx_permuteRows (
   DenseMtx   *mtx,
   IV         *oldToNewIV
) ;
/*
   ----------------------------------------------
   purpose -- to permute the columns of an object
 
   created -- 98may02, cca
   ----------------------------------------------
*/
void
DenseMtx_permuteColumns (
   DenseMtx   *mtx,
   IV         *oldToNewIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------------
   sort the rows so the row ids are in ascending order
   sort the columns so the column ids are in ascending order
 
   created -- 98may02, cca
   ---------------------------------------------------------
*/
void
DenseMtx_sort (
   DenseMtx   *mtx
) ;
/*
   -----------------------------------------------
   copy row irowA from mtxA into row irowB in mtxB
 
   created -- 98may02, cca
   -----------------------------------------------
*/
void
DenseMtx_copyRow (
   DenseMtx   *mtxB,
   int        irowB,
   DenseMtx   *mtxA,
   int        irowA
) ;
/*
   ---------------------------------------------------------------
   copy row irowA from mtxA into row irowB in mtxB
   and copy row index irowB from mtxB into row index irowA of mtxA
 
   created -- 98aug12, cca
   ---------------------------------------------------------------
*/
void
DenseMtx_copyRowAndIndex (
   DenseMtx   *mtxB,
   int        irowB,
   DenseMtx   *mtxA,
   int        irowA
) ;
/*
   ----------------------------------------------
   add row irowA from mtxA into row irowB in mtxB
 
   created -- 98aug12, cca
   ----------------------------------------------
*/
void
DenseMtx_addRow (
   DenseMtx   *mtxB,
   int        irowB,
   DenseMtx   *mtxA,
   int        irowA
) ;
/*
   -----------------------
   zero the entries
 
   created -- 98may16, cca
   -----------------------
*/
void
DenseMtx_zero (
   DenseMtx   *mtx
) ;
/*
   ------------------------
   fill with random entries
 
   created -- 98may16, cca
   ------------------------
*/
void
DenseMtx_fillRandomEntries (
   DenseMtx   *mtx,
   Drand      *drand
) ;
/*
   -----------------------------------
   compute three checksums
     sums[0] = sum of row indices
     sums[1] = sum of columns indices
     sums[2] = sum of entry magnitudes
 
   created -- 98may16, cca
   -----------------------------------
*/
void
DenseMtx_checksums (
   DenseMtx   *mtx,
   double     sums[]
) ;
/*
   -------------------------------------------
   return the maximum magnitude of the entries
 
   created -- 98may15, cca
   -------------------------------------------
*/
double
DenseMtx_maxabs (
   DenseMtx   *mtx
) ;
/*
   --------------------------------------------
   subtract one matrix from another, B := B - A
 
   created -- 98may25, cca
   --------------------------------------------
*/
void
DenseMtx_sub (
   DenseMtx   *mtxB,
   DenseMtx   *mtxA
) ;
/*
   ----------------------------------------------------
   purpose -- to copy a row of the matrix into a vector
 
   irow -- local row id
   vec  -- double vector to receive the row entries
 
   created -- 98jul31, cca
   ----------------------------------------------------
*/
void
DenseMtx_copyRowIntoVector (
   DenseMtx   *mtx,
   int        irow,
   double     *vec
) ;
/*
   ----------------------------------------------------
   purpose -- to copy a row of the matrix into a vector
 
   irow -- local row id
   vec  -- double vector to supply the row entries
 
   created -- 98jul31, cca
   ----------------------------------------------------
*/
void
DenseMtx_copyVectorIntoRow (
   DenseMtx   *mtx,
   int        irow,
   double     *vec 
) ;
/*
   ----------------------------------------------------
   purpose -- to add a row of the matrix into a vector
 
   irow -- local row id
   vec  -- double vector to supply the row entries
 
   created -- 98aug12, cca
   ----------------------------------------------------
*/
void
DenseMtx_addVectorIntoRow (
   DenseMtx   *mtx,
   int        irow,
   double     *vec
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in scale.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   purpose -- to scale a dense matrix object by a scalar
     A := alpha * A ;

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- A has invalid type
     -3 -- alpha is NULL

   created -- 98nov06, cca
   -----------------------------------------------------
*/
int
DenseMtx_scale (
   DenseMtx   *A,
   double     alpha[]
) ;
/*--------------------------------------------------------------------*/
