/*  A2.h  */

#include "../cfiles.h"
#include "../SPOOLES.h"
#include "../ZV.h"
#include "../IV.h"
#include "../DV.h"

/*====================================================================*/
/*
   -----------------------------------------------------------
   the A2 structure holds a 2-d dense array of double complex

   type -- type of entries
      1 -- real
      2 -- complex
   n1 -- number of rows (first dimension)
   n2 -- number of columns (second dimension)
   inc1 -- 
      increment in storage along first dimension,
      when inc1 == 1, the storage is column major
   inc2 -- 
      increment in storage along second dimension,
      when inc2 == 1, the storage is row major
   nowned -- number of owned entries starting at entries
      when > 0, storage pointed to by entries
      has been allocated here and can be free'd.
      when = 0, storage pointed to by entries
      has not been allocated here and cannot be free'd.
   entries -- base address for the entries

   if real then
      entry(i,j) is found in entries[i*inc1 + j*inc2] 
   else
      entry(i,j) is found in entries[2*(i*inc1 + j*inc2)] 
      and entries[2*(i*inc1 + j*inc2)+1] 
   endif
   -----------------------------------------------------------
*/
typedef
struct _A2 {
   int      type     ;
   int      n1       ;
   int      n2       ;
   int      inc1     ;
   int      inc2     ;
   int      nowned   ;
   double   *entries ;
} A2 ;

#define A2_IS_REAL(chv)    ((chv)->type == SPOOLES_REAL)
#define A2_IS_COMPLEX(chv) ((chv)->type == SPOOLES_COMPLEX)

#define A2_STRICT_LOWER 1
#define A2_LOWER        2
#define A2_DIAGONAL     3
#define A2_UPPER        4
#define A2_STRICT_UPPER 5
#define A2_ALL_ENTRIES  6

#define A2_BY_ROWS     0
#define A2_BY_COLUMNS  1

/*====================================================================*/
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
A2 *
A2_new (
   void
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_setDefaultFields (
   A2   *mtx
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may01, cca
   --------------------------------------------------
*/
void
A2_clearData (
   A2   *mtx
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may01, cca
   ------------------------------------------
*/
void
A2_free (
   A2   *mtx
) ;
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------
   return the number of rows in the array
 
   created -- 98may01, cca
   --------------------------------------
*/
int
A2_nrow (
   A2   *mtx
) ;
/*
   -----------------------------------------
   return the number of columns in the array
 
   created -- 98may01, cca
   -----------------------------------------
*/
int
A2_ncol (
   A2   *mtx
) ;
/*
   --------------------------
   return the first increment
 
   created -- 98may01, cca
   --------------------------
*/
int
A2_inc1 (
   A2  *mtx
) ;
/*
   ---------------------------
   return the second increment
 
   created -- 98may01, cca
   ---------------------------
*/
int
A2_inc2 (
   A2  *mtx
) ;
/*
   -------------------------------
   return a pointer to the entries
 
   created -- 98may01, cca
   -------------------------------
*/
double *
A2_entries (
   A2      *mtx
) ;
/*
   --------------------------------------------
   return a pointer to the first entry in a row
 
   created -- 98may01, cca
   --------------------------------------------
*/
double *
A2_row (
   A2   *mtx,
   int   irow
) ;
/*
   -----------------------------------------------
   return a pointer to the first entry in a column
 
   created -- 98may01, cca
   -----------------------------------------------
*/
double *
A2_column (
   A2   *mtx,
   int   jcol
) ;
/*
   -------------------------------------------
   fill *pValue with the entry in (irow, jcol)
 
   created -- 98may01, cca
   -------------------------------------------
*/
void
A2_realEntry (
   A2       *mtx,
   int      irow,
   int      jcol,
   double   *pValue
) ;
/*
   ---------------------------------------------------
   fill (*pReal,*pImag) with the entry in (irow, jcol)
 
   created -- 98may01, cca
   ---------------------------------------------------
*/
void
A2_complexEntry (
   A2       *mtx,
   int      irow,
   int      jcol,
   double   *pReal,
   double   *pImag
) ;
/*
   -----------------------------------------
   set the entry in (irow, jcol) to be value
 
   created -- 98may01, cca
   -----------------------------------------
*/
void
A2_setRealEntry (
   A2       *mtx,
   int      irow,
   int      jcol,
   double   value
) ;
/*
   -----------------------------------------------
   set the entry in (irow, jcol) to be (real,imag)
 
   created -- 98may01, cca
   -----------------------------------------------
*/
void
A2_setComplexEntry (
   A2       *mtx,
   int      irow,
   int      jcol,
   double   real,
   double   imag
) ;
/*
   ---------------------------------------
   fill pointers to the matrix first entry
   in row irow and column jcol
 
   created -- 98may01, cca
   ---------------------------------------
*/
void
A2_pointerToRealEntry (
   A2       *mtx,
   int      irow,
   int      jcol,
   double   **ppValue
) ;
/*
   ---------------------------------------
   fill pointers to the matrix first entry
   in row irow and column jcol
 
   created -- 98may01, cca
   ---------------------------------------
*/
void
A2_pointerToComplexEntry (
   A2       *mtx,
   int      irow,
   int      jcol,
   double   **ppReal,
   double   **ppImag
) ;
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   initializer. sets the n1, n2, inc1 and inc2 fields.
   must have
      mtx != NULL
      type = SPOOLES_REAL or SPOOLES_COMPLEX
      n1, n2, inc1, inc2 > 0
      (inc1 = 1 and inc2 = nrow) or (inc1 = ncol and inc2 = 1)
 
   if entries is NULL then
      entries are allocated and zeroed.
   else
      mtx->nowned is set to 0
      mtx->entries is set to entries.
   endif
   n1, n2, inc1, inc2 set
 
   created -- 98apr15, cca
   ------------------------------------------------------------------
*/
void
A2_init (
   A2       *mtx,
   int      type,
   int      n1,
   int      n2,
   int      inc1,
   int      inc2,
   double   *entries
) ;
/*
   --------------------------------------------------
   submatrix initializer
 
   A(0:lastrow-firstrow,0:lastcol-firstcol)
              = B(firstrow:lastrow, firstcol:lastcol)
 
   created -- 98apr15, cca
   --------------------------------------------------
*/
void
A2_subA2 (
   A2   *mtxA,
   A2   *mtxB,
   int   firstrow,
   int   lastrow,
   int   firstcol,
   int   lastcol
) ;
/*
------------------------------------------------------------------------
----- methods found in norms.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------
   return the entry of maximum magnitude
 
   created -- 98apr15, cca
   -------------------------------------
*/
double
A2_maxabs (
   A2   *a
) ;
/*
   ---------------------------------------
   return the frobenius norm of the matrix
 
   created -- 98apr15, cca
   ---------------------------------------
*/
double
A2_frobNorm (
   A2   *mtx
) ;
/*
   ---------------------------------
   return the one-norm of the matrix
 
   created -- 98apr15, cca
   ---------------------------------
*/
double
A2_oneNorm (
   A2   *mtx
) ;
/*
   --------------------------------------
   return the infinity-norm of the matrix
 
   created -- 98apr15, cca
   --------------------------------------
*/
double
A2_infinityNorm (
   A2   *mtx
) ;
/*
   ----------------------------------
   return the one-norm of column jcol
 
   created -- 98apr15, cca
   ----------------------------------
*/
double
A2_oneNormOfColumn (
   A2   *mtx,
   int   jcol
) ;
/*
   ----------------------------------
   return the two-norm of column jcol
 
   created -- 98apr15, cca
   ----------------------------------
*/
double
A2_twoNormOfColumn (
   A2   *mtx,
   int   jcol
) ;
/*
   ---------------------------------------
   return the infinity-norm of column jcol
 
   created -- 98apr15, cca
   ---------------------------------------
*/
double
A2_infinityNormOfColumn (
   A2   *mtx,
   int   jcol
) ;
/*
   -------------------------------
   return the one-norm of row irow
 
   created -- 98apr15, cca
   -------------------------------
*/
double
A2_oneNormOfRow (
   A2   *mtx,
   int   irow
) ;
/*
   -------------------------------
   return the two-norm of row irow
 
   created -- 98apr15, cca
   -------------------------------
*/
double
A2_twoNormOfRow (
   A2   *mtx,
   int   irow
) ;
/*
   ------------------------------------
   return the infinity-norm of row irow
 
   created -- 98apr15, cca
   ------------------------------------
*/
double
A2_infinityNormOfRow (
   A2   *mtx,
   int   irow
) ;
/*
------------------------------------------------------------------------
----- methods found in sort.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------
   permute the rows of the matrix
   A(*,*) = A(index(*),*)
   this method calls A2_sortRowsUp
   but does not overwrite the index[] vector
 
   created -- 98apr15, cca
   -----------------------------------------
*/
void
A2_permuteRows (
   A2    *mtx,
   int   nrow,
   int   index[]
) ;
/*
   -----------------------------------------
   permute the columns of the matrix
   A(*,*) = A(*,index(*))
   this method calls A2_sortColumnsUp
   but does not overwrite the index[] vector
 
   created -- 98apr15, cca
   -----------------------------------------
*/
void
A2_permuteColumns (
   A2   *mtx,
   int   ncol,
   int   index[]
) ;
/*
   ----------------------------------------------
   sort the rows of the matrix in ascending order
   of the rowids[] vector. on return, rowids is
   in asending order. return value is the number
   of row swaps made.
 
   created -- 98apr15, cca
   ----------------------------------------------
*/
int
A2_sortRowsUp (
   A2    *mtx,
   int   nrow,
   int   rowids[]
) ;
/*
   -------------------------------------------------
   sort the columns of the matrix in ascending order
   of the colids[] vector. on return, colids is
   in asending order. return value is the number
   of column swaps made.
 
   created -- 98apr15, cca
   -------------------------------------------------
*/
int
A2_sortColumnsUp (
   A2   *mtx,
   int   ncol,
   int   colids[]
) ;
/*
------------------------------------------------------------------------
----- methods found in QRreduce.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------
   purpose -- compute A = QR, where Q is a product of householder
              vectors, (I - beta_j v_j v_j^T). on return, v_j is 
              found in the lower triangle of A, v_j(j) = 1.0.

   return value -- # of floating point operations

   created -- 98may25, cca
   --------------------------------------------------------------
*/
double
A2_QRreduce (
   A2       *mtxA,
   DV       *workDV,
   int      msglvl,
   FILE     *msgFile
) ;
/*
   -----------------------------------------------------------
   A contains the following data from the A = QR factorization

   A(1:ncolA,1:ncolA) = R
   A(j+1:nrowA,j) is v_j, the j-th householder vector,
       where v_j[j] = 1.0

   NOTE: A and Q must be column major

   created -- 98dec10, cca
   -----------------------------------------------------------
*/
void
A2_computeQ (
   A2     *Q,
   A2     *A,
   DV     *workDV,
   int    msglvl,
   FILE   *msgFile
) ;
/*
   -----------------------------------------------------------
   A contains the following data from the A = QR factorization

   A(1:ncolA,1:ncolA) = R
   A(j+1:nrowA,j) is v_j, the j-th householder vector, 
       where v_j[j] = 1.0

   we compute Y = Q^T X when A is real
          and Y = Q^H X when A is complex

   NOTE: A, Y and X must be column major.
   NOTE: Y and X can be the same object,
         in which case X is overwritten with Y

   created -- 98dec10, cca
   -----------------------------------------------------------
*/
void
A2_applyQT (
   A2     *Y,
   A2     *A,
   A2     *X,
   DV     *workDV,
   int    msglvl,
   FILE   *msgFile
) ;
/*
------------------------------------------------------------------------
----- methods found in copyEntriesToVector.c ---------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------------
  purpose -- copy entries to a vector. the portion copied
              can be a union of the strict lower portion,
              the diagonal portion, and the strict upper
              portion. there is one restriction, if the strict
              lower and strict upper are to be copied, the
              diagonal will also be copied.
 
   length -- length of dvec[]
   dvec[] -- vector to receive matrix entries
   copyflag  -- flag to denote what part of the entries to move
      A2_STRICT_LOWER --> move strict lower entries
      A2_LOWER        --> move lower entries (includes the diagonal)
      A2_DIAGONAL     --> move diagonal entries
      A2_UPPER        --> move upper entries (includes the diagonal)
      A2_STRICT_UPPER --> move strict upper entries
      A2_ALL_ENTRIES  --> move all entries
   storeflag -- flag to denote how to store entries in dvec[]
      A2_BY_ROWS    --> store by rows
      A2_BY_COLUMNS --> store by columns
 
   return value -- number of entries copied
 
   created  -- 97jun03, cca, dkw
   modified -- 98may25, cca
   ----------------------------------------------------------------
*/
int
A2_copyEntriesToVector (
   A2       *mtx,
   int      length,
   double   *dvec,
   int      copyflag,
   int      storeflag
) ;
/*
------------------------------------------------------------------------
----- methods found in makeStaircase.c ---------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------------
   purpose -- to permute the rows so the matrix is in staircase form
 
   created -- 98may25, cca
   -----------------------------------------------------------------
*/
void
A2_makeStaircase (
   A2   *mtxA
) ;
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object
 
   created -- 98may01, cca
   ----------------------------------------------
*/
int
A2_sizeOf (
   A2   *mtx
) ;
/*
   ---------------------------------------------------------------
   shift the base of the entries and adjust dimensions
 
   mtx(0:n1-rowoff-1,0:n2-coloff-1) = mtx(rowoff:n1-1,coloff:n2-1) 
 
   created -- 98may01, cca
   ---------------------------------------------------------------
*/
void
A2_shiftBase (
   A2   *mtx,
   int   rowoff,
   int   coloff
) ;
/*
   --------------------------------------------------------------
   returns 1 if the storage is row major, otherwise returns zero.
 
   created -- 98may01, cca
   --------------------------------------------------------------
*/
int
A2_rowMajor (
   A2   *mtx
) ;
/*
   -----------------------------------------------------------------
   returns 1 if the storage is column major, otherwise returns zero.
 
   created -- 98may01, cca
   -----------------------------------------------------------------
*/
int
A2_columnMajor (
   A2   *mtx
) ;
/*
   -----------------------
   transpose the matrix
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_transpose (
   A2   *mtx
) ;
/*
   ----------------------------
   extract row[*] = mtx(irow,*)
 
   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractRow (
   A2      *mtx,
   double   row[],
   int      irow
) ;
/*
   ----------------------------
   extract col[*] = mtx(*,jcol)
 
   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractColumn (
   A2      *mtx,
   double   col[],
   int      jcol
) ;
/*
   -----------------------
   set mtx(irow,*) = y[*]
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_setRow (
   A2      *mtx,
   double   row[],
   int      irow
) ;
/*
   -----------------------
   set mtx(*,jcol) = y[*]
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_setColumn (
   A2      *mtx,
   double   col[],
   int      jcol
) ;
/*
   ----------------------------
   extract row[*] = mtx(irow,*)
 
   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractRowDV (
   A2   *mtx,
   DV    *rowDV,
   int   irow
) ;
/*
   ----------------------------
   extract row[*] = mtx(irow,*)
 
   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractRowZV (
   A2   *mtx,
   ZV    *rowZV,
   int   irow
) ;
/*
   ----------------------------
   extract col[*] = mtx(*,jcol)
 
   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractColumnDV (
   A2   *mtx,
   DV    *colDV,
   int   jcol
) ;
/*
   ----------------------------
   extract col[*] = mtx(*,jcol)
 
   created -- 98may01, cca
   ----------------------------
*/
void
A2_extractColumnZV (
   A2   *mtx,
   ZV    *colZV,
   int   jcol
) ;
/*
   -----------------------
   set mtx(irow,*) = y[*]
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_setRowDV (
   A2      *mtx,
   DV       *rowDV,
   int      irow
) ;
/*
   -----------------------
   set mtx(irow,*) = y[*]
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_setRowZV (
   A2      *mtx,
   ZV       *rowZV,
   int      irow
) ;
/*
   -----------------------
   set mtx(*,jcol) = y[*]
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_setColumnDV (
   A2      *mtx,
   DV       *colDV,
   int      jcol
) ;
/*
   -----------------------
   set mtx(*,jcol) = y[*]
 
   created -- 98may01, cca
   -----------------------
*/
void
A2_setColumnZV (
   A2      *mtx,
   ZV       *colZV,
   int      jcol
) ;
/*
   -------------------------------------------------------------
   fill the matrix with uniform random numbers in [lower, upper]
 
   created -- 98may01, cca
   -------------------------------------------------------------
*/
void
A2_fillRandomUniform (
   A2       *a,
   double   lower,
   double   upper,
   int      seed
) ;
/*
   -----------------------------------------------
   fill the matrix with normal(0,1) random numbers
 
   created -- 98may01, cca
   -----------------------------------------------
*/
void
A2_fillRandomNormal (
   A2      *a,
   double   mean,
   double   variance,
   int      seed
) ;
/*
   ----------------------------------------
   fill the matrix with the identity matrix
 
   created -- 98may01, cca
   ----------------------------------------
*/
void
A2_fillWithIdentity (
   A2   *a
) ;
/*
   --------------------------
   fill the matrix with zeros
 
   created -- 98may01, cca
   --------------------------
*/
void
A2_zero (
   A2   *a
) ;
/*
   ----------------------------
   copy one matrix into another
      A := B
 
   created  -- 98may01, cca
   ----------------------------
*/
void
A2_copy (
   A2   *A,
   A2   *B
) ;
/*
   --------------------------------
   subtract one matrix from another
 
   A := A - B
 
   created -- 98may01, cca
   ----------------------------
*/
void
A2_sub (
   A2   *A,
   A2   *B
) ;
/*
   ---------------------------
   swap two rows of the matrix
 
   created -- 98may01, cca
   ---------------------------
*/
void
A2_swapRows (
   A2   *a,
   int   irow1,
   int   irow2
) ;
/*
   ------------------------------
   swap two columns of the matrix
 
   created -- 98may01, cca
   ------------------------------
*/
void
A2_swapColumns (
   A2   *a,
   int   jcol1,
   int   jcol2
) ;
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------
   purpose -- to read in the object from a file
 
   input --
 
      fn -- filename, must be *.a2b or *.a2f
 
   return value -- 1 if success, 0 if failure
 
   created -- 98may01, cca
   --------------------------------------------
*/
int
A2_readFromFile (
   A2     *mtx,
   char   *fn
) ;
/*
   --------------------------------------------------
   purpose -- to read an object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 98may01, cca
   --------------------------------------------------
*/
int
A2_readFromFormattedFile (
   A2    *mtx,
   FILE   *fp
) ;
/*
   -----------------------------------------------
   purpose -- to read an object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 98may01, cca
   -----------------------------------------------
*/
int
A2_readFromBinaryFile (
   A2    *mtx,
   FILE   *fp
) ;
/*
   ---------------------------------------
   purpose -- to write an object to a file
 
   input --
 
      fn -- filename
        *.a2b -- binary
        *.a2f -- formatted
        anything else -- for human eye
 
   created -- 98may01, cca
   ---------------------------------------
*/
void
A2_writeToFile (
   A2    *mtx,
   char   *fn
) ;
/*
   -------------------------------------------------
   purpose -- to write an object to a formatted file
 
   created -- 98may01, cca
   -------------------------------------------------
*/
void
A2_writeToFormattedFile (
   A2    *mtx,
   FILE   *fp
) ;
/*
   --------------------------------------------------------
   purpose -- to write an adjacency object to a binary file
 
   created -- 98may01, cca
   --------------------------------------------------------
*/
void
A2_writeToBinaryFile (
   A2    *mtx,
   FILE   *fp
) ;
/*
   ----------------------------------------------
   purpose -- to write the object for a human eye
 
   created -- 98may01, cca
   ----------------------------------------------
*/
void
A2_writeForHumanEye (
   A2    *mtx,
   FILE   *fp
) ;
/*
   --------------------------------------
   purpose -- to write out the statistics
 
   created -- 98may01, cca
   --------------------------------------
*/
void
A2_writeStats (
   A2    *mtx,
   FILE   *fp
) ;
/*
   -----------------------------------------------
   purpose -- to write the matrix in matlab format
 
   created -- 98may01, cca
   -----------------------------------------------
*/
void
A2_writeForMatlab (
   A2    *mtx,
   char   *mtxname,
   FILE   *fp
) ;
/*====================================================================*/
