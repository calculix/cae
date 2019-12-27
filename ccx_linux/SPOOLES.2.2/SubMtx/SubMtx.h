/*  SubMtx.h  */

#include "../cfiles.h"
#include "../A2.h"
#include "../ZV.h"
#include "../DV.h"
#include "../SPOOLES.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   structure to hold a block of a block matrix.
   e.g., the block partition would be by supernodes.

   type -- type of matrix
      1 -- real
      2 -- complex
   mode -- storage mode of matrix
      0 -- dense stored by rows
      1 -- dense stored by columns
      2 -- sparse stored by rows
      3 -- sparse stored by columns
      4 -- sparse stored by triples
      5 -- sparse stored by dense subrows
      6 -- sparse stored by dense subcolumns
      7 -- diagonal
      8 -- block diagonal symmetric 
           (upper part of each block stored by rows)
      9 -- block diagonal hermitian (complex only)
           (upper part of each block stored by rows)
   rowid    -- block row to which this matrix belongs
   colid    -- block column to which this matrix belongs
   nrow     -- number of rows
   ncol     -- number of columns
   nent     -- number of entries
   entries  -- pointer to the entries
   wrkDV    -- DV object that manages workspace
   next     -- pointer to next SubMtx object in a singly linked list
   
   created -- 98jan29, cca
   -----------------------------------------------------------------
*/
typedef struct _SubMtx   SubMtx ;
struct _SubMtx {
   int      type     ;
   int      mode     ;
   int      rowid    ;
   int      colid    ;
   int      nrow     ;
   int      ncol     ;
   int      nent     ;
   double   *entries ;
   DV       wrkDV    ;
   SubMtx   *next    ;
} ;
/*
   ------------------------
   real/complex type macros
   ------------------------
*/
#define SUBMTX_IS_REAL(chv)    ((chv)->type == SPOOLES_REAL)
#define SUBMTX_IS_COMPLEX(chv) ((chv)->type == SPOOLES_COMPLEX)
/*
   -------------------
   storage mode macros
   -------------------
*/
#define SUBMTX_DENSE_ROWS           0
#define SUBMTX_DENSE_COLUMNS        1
#define SUBMTX_SPARSE_ROWS          2
#define SUBMTX_SPARSE_COLUMNS       3
#define SUBMTX_SPARSE_TRIPLES       4
#define SUBMTX_DENSE_SUBROWS        5
#define SUBMTX_DENSE_SUBCOLUMNS     6
#define SUBMTX_DIAGONAL             7
#define SUBMTX_BLOCK_DIAGONAL_SYM   8
#define SUBMTX_BLOCK_DIAGONAL_HERM  9

#define SUBMTX_IS_DENSE_ROWS(chv) \
   ((chv)->mode == SUBMTX_DENSE_ROWS)
#define SUBMTX_IS_DENSE_COLUMNS(chv) \
   ((chv)->mode == SUBMTX_DENSE_COLUMNS)
#define SUBMTX_IS_SPARSE_ROWS(chv) \
   ((chv)->mode == SUBMTX_SPARSE_ROWS)
#define SUBMTX_IS_SPARSE_COLUMNS(chv) \
   ((chv)->mode == SUBMTX_SPARSE_COLUMNS)
#define SUBMTX_IS_SPARSE_TRIPLES(chv) \
   ((chv)->mode == SUBMTX_SPARSE_TRIPLES)
#define SUBMTX_IS_DENSE_SUBROWS(chv) \
   ((chv)->mode == SUBMTX_DENSE_SUBROWS)
#define SUBMTX_IS_DENSE_SUBCOLUMNS(chv) \
   ((chv)->mode == SUBMTX_DENSE_SUBCOLUMNS)
#define SUBMTX_IS_DIAGONAL(chv) \
   ((chv)->mode == SUBMTX_DIAGONAL)
#define SUBMTX_IS_BLOCK_DIAGONAL_SYM(chv) \
   ((chv)->mode == SUBMTX_BLOCK_DIAGONAL_SYM)
#define SUBMTX_IS_BLOCK_DIAGONAL_HERM(chv) \
   ((chv)->mode == SUBMTX_BLOCK_DIAGONAL_HERM)
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98may01, cca
   -----------------------
*/
SubMtx *
SubMtx_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may01, cca
   -----------------------
*/
void
SubMtx_setDefaultFields (
   SubMtx   *mtx
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may01, cca
   --------------------------------------------------
*/
void
SubMtx_clearData (
   SubMtx   *mtx
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may01, cca
   ------------------------------------------
*/
void
SubMtx_free (
   SubMtx   *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------
   purpose -- fill *prowid with the row id
              and  *pcolid with the column id
 
   created -- 98may01, cca
   ------------------------------------------
*/
void
SubMtx_ids (
   SubMtx   *mtx,
   int      *prowid,
   int      *pcolid
) ;
/*
   -------------------------------------
   purpose -- set the row and column ids
 
   created -- 98may01, cca
   -------------------------------------
*/
void
SubMtx_setIds (
   SubMtx   *mtx,
   int      rowid,
   int      colid
) ;
/*
   ---------------------------------------------------
   purpose -- fill *pnrow with the # of rows
              and  *pncol with the # of columns
              and  *pnent with the # of matrix entries
 
   created -- 98may01, cca
   ---------------------------------------------------
*/
void
SubMtx_dimensions (
   SubMtx   *mtx,
   int      *pnrow,
   int      *pncol,
   int      *pnent
) ;
/*
   ------------------------------------------------------------
   purpose --
     fill *pnrow with the number of rows
     if prowind != NULL then
        fill *prowind with the base location of the row indices
     endif
 
   created -- 98may01, cca
   ------------------------------------------------------------
*/
void
SubMtx_rowIndices (
   SubMtx   *mtx,
   int      *pnrow,
   int      **prowind
) ;
/*
   ---------------------------------------------------------------
   purpose --
     fill *pncol with the number of column
     if pcolind != NULL then
        fill *prowind with the base location of the column indices
     endif
 
   created -- 98may01, cca
   ---------------------------------------------------------------
*/
void
SubMtx_columnIndices (
   SubMtx   *mtx,
   int      *pncol,
   int      **pcolind
) ;
/*
   ---------------------------------
   purpose -- for dense storage
     *pnrow    with mtx->nrow
     *pncol    with mtx->ncol
     *pinc1    with row increment
     *pinc2    with column increment
     *pentries with mtx->entries
 
   created -- 98may01, cca
   ---------------------------------
*/
void
SubMtx_denseInfo (
   SubMtx     *mtx,
   int        *pnrow,
   int        *pncol,
   int        *pinc1,
   int        *pinc2,
   double     **pentries
) ;
/*
   ---------------------------------------------
   purpose -- for sparse rows, fill
     *pnrow    with # of rows
     *pnent    with # of indices
     *psizes   with sizes[] of rows
     *pindices with indices[] for column indices
     *pentries with entries[] for matrix entries
 
   created -- 98may01, cca
   ---------------------------------------------
*/
void
SubMtx_sparseRowsInfo (
   SubMtx     *mtx,
   int        *pnrow,
   int        *pnent,
   int        **psizes,
   int        **pindices,
   double     **pentries
) ;
/*
   ----------------------------------------------
   purpose -- for sparse columns, fill
     *pncol    with # of columns
     *pnent    with # of matrix entries
     *psizes   with sizes[ncol], column sizes
     *pindices with indices[nent], matrix row ids
     *pentries with entries[nent], matrix entries
 
   created -- 98may01, cca
   ----------------------------------------------
*/
void
SubMtx_sparseColumnsInfo (
   SubMtx     *mtx,
   int        *pncol,
   int        *pnent,
   int        **psizes,
   int        **pindices,
   double     **pentries
) ;
/*
   ----------------------------------------------------
   purpose -- for sparse triples, fill
     *pnent    with # of matrix entries
     *prowids  with rowids[nent], row ids of entries
     *prowids  with colids[nent], column ids of entries
     *pentries with entries[nent], matrix entries
 
   created -- 98may01, cca
   ----------------------------------------------------
*/
void
SubMtx_sparseTriplesInfo (
   SubMtx     *mtx,
   int        *pnent,
   int        **prowids,
   int        **pcolids,
   double     **pentries
) ;
/*
   -----------------------------------------------------------
   purpose -- for dense subrows, fill
     *pnrow      with # of rows
     *pnent      with # of matrix entries
     *pfirstlocs with firstlocs[nrow], column of first nonzero
     *psizes     with sizes[nrow], number of nonzero columns
     *pentries   with entries[nent], matrix entries
 
   created -- 98may01, cca
   -----------------------------------------------------------
*/
void
SubMtx_denseSubrowsInfo (
   SubMtx     *mtx,
   int        *pnrow,
   int        *pnent,
   int        **pfirstlocs,
   int        **psizes,
   double     **pentries
) ;
/*
   -----------------------------------------------------------
   purpose -- for dense subcolumns, fill
     *pncol      with # of columns
     *pnent      with # of matrix entries
     *pfirstlocs with firstlocs[ncol], row of first nonzero
     *psizes     with sizes[ncol], number of nonzero rows
     *pentries   with entries[nent], matrix entries
 
   created -- 98may01, cca
   -----------------------------------------------------------
*/
void
SubMtx_denseSubcolumnsInfo (
   SubMtx   *mtx,
   int      *pncol,
   int      *pnent,
   int      **pfirstlocs,
   int      **psizes,
   double   **pentries
) ;
/*
   ------------------------------------------------
   purpose -- for a diagonal matrix, fill
     *pncol      with # of columns
     *pentries   with entries[nent], matrix entries
 
   created -- 98may01, cca
   ------------------------------------------------
*/
void
SubMtx_diagonalInfo (
   SubMtx   *mtx,
   int      *pncol,
   double   **pentries
) ;
/*
   ------------------------------------------------------
   purpose -- for a block diagonal symmetric matrix, fill
     *pncol        with # of columns
     *pnent        with # of entries
     *ppivotsizes  with pivotsizes[ncol]
     *pentries     with entries[nent], matrix entries
 
   created -- 98may01, cca
   ------------------------------------------------------
*/
void
SubMtx_blockDiagonalInfo (
   SubMtx   *mtx,
   int      *pncol,
   int      *pnent,
   int      **ppivotsizes,
   double   **pentries
) ;
/*
   -------------------------------------------------------
   purpose -- to find matrix entry (irow,jcol) if present.
 
   return value --
     if entry (irow,jcol) is not present then
        *pValue is 0.0
        return value is -1
     else entry (irow,jcol) is present then
        *pValue is the matrix entry
        return value is offset into entries array
     endif
 
   created -- 98may01, cca
   -------------------------------------------------------
*/
int
SubMtx_realEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   *pValue
) ;
/*
   -------------------------------------------------------
   purpose -- to find matrix entry (irow,jcol) if present.
 
   return value --
     if entry (irow,jcol) is not present then
        *pReal and *pImag are 0.0
        return value is -1
     else entry (irow,jcol) is present then
        (*pReal,*pImag) is the matrix entry
        return value is offset into entries array
     endif
 
   created -- 98may01, cca
   -------------------------------------------------------
*/
int
SubMtx_complexEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   *pReal,
   double   *pImag
) ;
/*
   -------------------------------------------------
   purpose -- to return a pointer to the location of
               matrix entry (irow,jcol) if present.
 
   if entry (irow,jcol) is not present then
      *ppValue is NULL
   else entry (irow,jcol) is present then
      *ppValue is the location of the matrix entry
   endif
 
   created -- 98may01, cca
   -------------------------------------------------
*/
void
SubMtx_locationOfRealEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   **ppValue
) ;
/*
   --------------------------------------------------------
   purpose -- to return a pointer to the location of
               matrix entry (irow,jcol) if present.
 
   if entry (irow,jcol) is not present then
      (*ppReal,*ppImag) is (NULL,NULL)
   else entry (irow,jcol) is present then
      (*ppReal,*ppImag) is the location of the matrix entry
   endif
 
   created -- 98may01, cca
   --------------------------------------------------------
*/
void
SubMtx_locationOfComplexEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   **ppReal,
   double   **ppImag
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
 
   created -- 98may01, cca
   ------------------------------------------------------------
*/
int
SubMtx_nbytesNeeded (
   int   type,
   int   mode,
   int   nrow,
   int   ncol,
   int   nent
) ;
/*
   ------------------------------------
   return the number of bytes in use in
   the workspace owned by this object
 
   created -- 98may01, cca
   ------------------------------------
*/
int
SubMtx_nbytesInUse (
   SubMtx   *mtx
) ;
/*
   ------------------------------------------------------
   return a pointer to the workspace owned by this object
 
   created -- 98may01, cca
   ------------------------------------------------------
*/
void *
SubMtx_workspace (
   SubMtx   *mtx
) ;
/*
   ----------------------------------------------------------------
  return the number of bytes in the workspace owned by this object
 
   created -- 98may01, cca
   ----------------------------------------------------------------
*/
int
SubMtx_nbytesInWorkspace (
   SubMtx   *mtx
) ;
/*
   -------------------------------------------------------------
   set the number of bytes in the workspace owned by this object
 
   created -- 98may01, cca
   -------------------------------------------------------------
*/
void
SubMtx_setNbytesInWorkspace (
   SubMtx   *mtx,
   int      nbytes
) ;
/*
   ---------------------------------------
   purpose -- set the fields of the object
 
   created -- 98may01, cca
   ---------------------------------------
*/
void
SubMtx_setFields (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent
) ;
/*
   ----------------------------
   purpose -- basic initializer
 
   created -- 98may01, cca
   ----------------------------
*/
void
SubMtx_init (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent
) ;
/*
   ---------------------------------------------------------
   purpose -- initialize the object from its working storage
              used when the object is a MPI message
 
   created -- 98may01, cca
   ---------------------------------------------------------
*/
void
SubMtx_initFromBuffer (
   SubMtx   *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in initRandom.c ------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------
   purpose -- to initialize a matrix object with 
      random entries and possibly random structure.
 
   created -- 98feb16, cca
   ------------------------------------------------
*/
void
SubMtx_initRandom (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent,
   int      seed
) ;
/*
   ------------------------------------------------
   purpose -- to initialize a matrix object with
      random entries in the upper triangle
 
   strict = 1 --> strict upper triangle
 
   created -- 98feb16, cca
   ------------------------------------------------
*/
void
SubMtx_initRandomUpperTriangle (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent,
   int      seed,
   int      strict
) ;
/*
   ------------------------------------------------
   purpose -- to initialize a matrix object with
      random entries in the lower triangle
 
   strict = 1 --> strict lower triangle
 
   created -- 98feb16, cca
   ------------------------------------------------
*/
void
SubMtx_initRandomLowerTriangle (
   SubMtx   *mtx,
   int      type,
   int      mode,
   int      rowid,
   int      colid,
   int      nrow,
   int      ncol,
   int      nent,
   int      seed,
   int      strict
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in scalevec.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------
   purpose -- compute [y0] := A * [x0]
 
   created -- 98apr17, cca
   -----------------------------------
*/
void
SubMtx_scale1vec (
   SubMtx   *mtxA,
   double   y0[],
   double   x0[]
) ;
/*
   -------------------------------------
   purpose -- compute [y0 y1] := [x0 x1]
 
   created -- 98apr17, cca
   -------------------------------------
*/
void
SubMtx_scale2vec (
   SubMtx     *mtxA,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[]
) ;
/*
   -----------------------------------------------
   purpose -- compute [y0 y1 y2] := A * [x0 x1 x2]
 
   created -- 98apr17, cca
   -----------------------------------------------
*/
void
SubMtx_scale3vec (
   SubMtx     *mtxA,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solve.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------------
   purpose -- solve A X = B, 
      where 
        (1) X overwrites B
        (2) A must be lower or upper triangular, 
            diagonal or block diagonal
        (3) if A is strict lower or upper triangular,
            then we solve (I + A) X = B
        (4) columns(A) = rows(X)
        (5) rows(A)    = rows(B)
        (6) B has type SUBMTX_DENSE_COLUMNS
        (7) if A is SUBMTX_DENSE_SUBROWS or SUBMTX_SPARSE_ROWS
            then A must be strict lower triangular
        (8) if A is SUBMTX_DENSE_SUBCOLUMNS or SUBMTX_SPARSE_COLUMNS
            then A must be strict upper triangular
        (9) A can be SUBMTX_DIAGONAL, SUBMTX_BLOCK_DIAGONAL_SYM
            or SUBMTX_BLOCK_DIAGONAL_HERM
 
   created -- 98may01, cca
   -----------------------------------------------------------------
*/
void
SubMtx_solve (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solveT.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------------
   purpose -- solve (A^T + I) X = B, 
      where 
        (1) X overwrites B
        (2) A must be strict lower or upper triangular
        (3) columns(A) = rows(X)
        (4) rows(A)    = rows(B)
        (5) B has type SUBMTX_DENSE_COLUMNS
        (6) if A is SUBMTX_DENSE_SUBROWS or SUBMTX_SPARSE_ROWS
            then A must be strict lower triangular
        (7) if A is SUBMTX_DENSE_SUBCOLUMNS or SUBMTX_SPARSE_COLUMNS
            then A must be strict upper triangular
 
   created -- 98may01, cca
   -----------------------------------------------------------------
*/
void
SubMtx_solveT (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solveH.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   purpose -- solve (A^H + I) X = B,
      where
        (1) X overwrites B
        (2) A must be strict lower or upper triangular
        (2) A, B and X are complex
        (4) columns(A) = rows(X)
        (5) rows(A)    = rows(B)
        (6) B has mode SUBMTX_DENSE_COLUMNS
        (7) if A is SUBMTX_DENSE_SUBROWS or SUBMTX_SPARSE_ROWS
            then A must be strict lower triangular
        (8) if A is SUBMTX_DENSE_SUBCOLUMNS or SUBMTX_SPARSE_COLUMNS
            then A must be strict upper triangular
 
   created -- 98may01, cca
   -------------------------------------------------------------
*/
void
SubMtx_solveH (
   SubMtx   *mtxA,
   SubMtx   *mtxB
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solveupd.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------
   purpose -- perform the matrix-matrix multiply 
      Y := Y - A * X used in the forward and backsolves
      where
        (1) rows(A) \subseteq rows(Y)
        (2) rows(A) are local w.r.t. rows(Y)
        (3) cols(A) \subseteq rows(X)
        (4) cols(A) are local w.r.t. rows(X)
        (5) cols(Y) = cols(X)
        (6) Y and X have mode SUBMTX_DENSE_COLUMNS
 
   created -- 98may02, cca
   ----------------------------------------------------
*/
void
SubMtx_solveupd (
   SubMtx   *mtxY,
   SubMtx   *mtxA,
   SubMtx   *mtxX
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solveupdT.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------
   purpose -- perform the matrix-matrix multiply 
      Y := Y - A^T * X used in the forward and backsolves
      where
        (1) rows(A) \subseteq rows(Y)
        (2) rows(A) are local w.r.t. rows(Y)
        (3) cols(A) \subseteq rows(X)
        (4) cols(A) are local w.r.t. rows(X)
        (5) cols(Y) = cols(X)
        (6) Y and X have type SUBMTX_DENSE_COLUMNS
 
   created -- 98may02, cca
   -----------------------------------------------------
*/
void
SubMtx_solveupdT (
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in solveupdH.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------
   purpose -- perform the matrix-matrix multiply 
      Y := Y - A^H * X used in the forward and backsolves
      where
        (1) rows(A) \subseteq rows(Y)
        (2) rows(A) are local w.r.t. rows(Y)
        (3) cols(A) \subseteq rows(X)
        (4) cols(A) are local w.r.t. rows(X)
        (5) cols(Y) = cols(X)
        (6) Y and X have type SUBMTX_DENSE_COLUMNS
 
   created -- 98may02, cca
   -----------------------------------------------------
*/
void
SubMtx_solveupdH (
   SubMtx     *mtxY,
   SubMtx     *mtxA,
   SubMtx     *mtxX
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in sort.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- sort the rows of the matrix into ascending order
 
   created -- 98mar02, cca
   -----------------------------------------------------------
*/
void
SubMtx_sortRowsUp (
   SubMtx   *mtx
) ;
/*
   --------------------------------------------------------------
   purpose -- sort the columns of the matrix into ascending order
 
   created -- 98mar02, cca
   --------------------------------------------------------------
*/
void
SubMtx_sortColumnsUp (
   SubMtx   *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------
   purpose -- to fill rowDV with the entries
              in row irow of the matrix
 
   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillRowDV (
   SubMtx   *mtx,
   int      irow,
   DV       *rowDV
) ;
/*
   -----------------------------------------
   purpose -- to fill rowZV with the entries
              in row irow of the matrix
 
   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillRowZV (
   SubMtx   *mtx,
   int      irow,
   ZV       *rowZV
) ;
/*
   -----------------------------------------
   purpose -- to fill colDV with the entries
              in column icol of the matrix
 
   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillColumnDV (
   SubMtx   *mtx,
   int      icol,
   DV       *colDV
) ;
/*
   -----------------------------------------
   purpose -- to fill colZV with the entries
              in column icol of the matrix
 
   created -- 98may01, cca
   -----------------------------------------
*/
void
SubMtx_fillColumnZV (
   SubMtx   *mtx,
   int      icol,
   ZV       *colZV
) ;
/*
   ------------------------------------------------------
   purpose -- return the magnitude of the largest element
 
   created -- 98may01, cca
   ------------------------------------------------------
*/
double
SubMtx_maxabs (
   SubMtx   *mtx
) ;
/*
   -----------------------------------------
   purpose -- zero the entries in the matrix
 
   created -- 98may04, cca
   -----------------------------------------
*/
void
SubMtx_zero (
   SubMtx   *mtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------
   purpose -- to read in an SubMtx object from a file
 
   input --
 
      fn -- filename, must be *.submtxb or *.submtxf
 
   return value -- 1 if success, 0 if failure
 
   created -- 98feb15, cca
   --------------------------------------------------
*/
int
SubMtx_readFromFile (
   SubMtx   *mtx,
   char     *fn
) ;
/*
   -------------------------------------------------------
   purpose -- to read an SubMtx object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 96jun23, cca
   -------------------------------------------------------
*/
int
SubMtx_readFromFormattedFile (
   SubMtx    *mtx,
   FILE   *fp
) ;
/*
   ------------------------------------------------------
   purpose -- to read an SubMtx object from a binary file
 
   return value -- 1 if success, 0 if failure
 
   created -- 96jun23, cca
   ------------------------------------------------------
*/
int
SubMtx_readFromBinaryFile (
   SubMtx   *mtx,
   FILE     *fp
) ;
/*
   ----------------------------------------------
   purpose -- to write an SubMtx object to a file
 
   input --
 
      fn -- filename
        *.submtxb -- binary
        *.submtxf -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 96jun23, cca
   ----------------------------------------------
*/
int
SubMtx_writeToFile (
   SubMtx   *mtx,
   char     *fn
) ;
/*
   --------------------------------------------------------
   purpose -- to write an SubMtx object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 96jun23, cca
   --------------------------------------------------------
*/
int
SubMtx_writeToFormattedFile (
   SubMtx   *mtx,
   FILE     *fp
) ;
/*
   -----------------------------------------------------
   purpose -- to write an SubMtx object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 96jun23, cca
   -----------------------------------------------------
*/
int
SubMtx_writeToBinaryFile (
   SubMtx   *mtx,
   FILE     *fp
) ;
/*
   -------------------------------------------------------
   purpose -- write a SubMtx object in human readable form
 
   created -- 98feb07, cca
   -------------------------------------------------------
*/
int
SubMtx_writeForHumanEye (
   SubMtx   *mtx,
   FILE     *fp
) ;
/*
   -------------------------------------------------------------------
   purpose -- write the header and scalar quantities for a SubMtx object
 
   created -- 98feb07, cca
   -------------------------------------------------------------------
*/
int
SubMtx_writeStats (
   SubMtx   *mtx,
   FILE   *fp
) ;
/*
   --------------------------------------------------
   purpose -- write the matrix entries out for matlab
 
   created -- 98feb07, cca
   --------------------------------------------------
*/
void
SubMtx_writeForMatlab (
   SubMtx   *mtx,
   char   *mtxname,
   FILE   *fp
) ;
/*--------------------------------------------------------------------*/
