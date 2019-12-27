/*  InpMtx.h  */

#include "../DenseMtx.h"
#include "../IV.h"
#include "../IVL.h"
#include "../DV.h"
#include "../ZV.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   coordType -- coordinate type
      0 -- no type specified
      1 -- row triples (i, j, a_{i,j})
      2 -- column triples (i, j, a_{j,i})
      3 -- chevron triples 
           (i, k, a_{i,i+k}) if k >= 0
           (i, k, a_{i-k,i}) if k < 0
           i is the chevron, k is the offset
      4 -- custom coordinate, e.g., one could store (I, k, a_{i,j})
           where I is the front where a_{i,j} will be assembled
           and k is the offset into the vector that holds the 
           entries in the front
   storageMode -- storage mode
      0 -- no mode specified
      1 -- filled with raw triples
      2 -- filled with sorted and distinct triples
      3 -- vectors by the first coordinate, ivec1[*] ignored
   inputMode -- input mode
      0 -- no input allowed
      1 -- indices only
      2 -- indices and real entries
      3 -- indices and complex entries
   maxnent -- present maximum number of entries
   nent    -- present number of entries
   resizeMultiple -- when resizing is done, 
      new maxnent = old maxnent * resizeMultiple
   ivec1IV    -- IV object that holds the first coordinates
   ivec2IV    -- IV object that holds the second coordinates
   dvecDV     -- DV object that holds the entries
   maxnvector -- present number of vectors
   nvector    -- present number of vectors
   vecidsIV   -- IV object that holds vector ids
   sizesIV    -- IV object that holds vector sizes
   offsetsIV  -- IV object that holds vector offsets

   created -- 98jan28, cca
   ------------------------------------------------------------------
*/
typedef struct _InpMtx   InpMtx ;
struct _InpMtx {
   int      coordType      ;
   int      storageMode    ;
   int      inputMode      ;
   int      maxnent        ;
   int      nent           ;
   double   resizeMultiple ;
   IV       ivec1IV        ;
   IV       ivec2IV        ;
   DV       dvecDV         ;
   int      maxnvector     ;
   int      nvector        ;
   IV       vecidsIV       ;
   IV       sizesIV        ;
   IV       offsetsIV      ;
} ;
#define INPMTX_NO_TYPE         0
#define INPMTX_BY_ROWS         1
#define INPMTX_BY_COLUMNS      2
#define INPMTX_BY_CHEVRONS     3
#define INPMTX_CUSTOM          4
#define INPMTX_IS_BY_ROWS(inpmtx)     \
   ((inpmtx)->coordType == INPMTX_BY_ROWS)
#define INPMTX_IS_BY_COLUMNS(inpmtx)  \
   ((inpmtx)->coordType == INPMTX_BY_COLUMNS)
#define INPMTX_IS_BY_CHEVRONS(inpmtx) \
   ((inpmtx)->coordType == INPMTX_BY_CHEVRONS)
#define INPMTX_IS_CUSTOM(inpmtx)      \
   ((inpmtx)->coordType == INPMTX_CUSTOM)

#define INPMTX_NO_MODE         0
#define INPMTX_RAW_DATA        1
#define INPMTX_SORTED          2
#define INPMTX_BY_VECTORS      3
#define INPMTX_IS_RAW_DATA(inpmtx)   \
   ((inpmtx)->storageMode == INPMTX_RAW_DATA)
#define INPMTX_IS_SORTED(inpmtx)     \
   ((inpmtx)->storageMode == INPMTX_SORTED)
#define INPMTX_IS_BY_VECTORS(inpmtx) \
   ((inpmtx)->storageMode == INPMTX_BY_VECTORS)

#define INPMTX_INDICES_ONLY    0
#define INPMTX_REAL_ENTRIES    1
#define INPMTX_COMPLEX_ENTRIES 2
#define INPMTX_IS_INDICES_ONLY(inpmtx) \
   ((inpmtx)->inputMode == INPMTX_INDICES_ONLY)
#define INPMTX_IS_REAL_ENTRIES(inpmtx) \
   ((inpmtx)->inputMode == SPOOLES_REAL)
#define INPMTX_IS_COMPLEX_ENTRIES(inpmtx) \
   ((inpmtx)->inputMode == SPOOLES_COMPLEX)
/*
#define INPMTX_IS_REAL_ENTRIES(inpmtx) \
   ((inpmtx)->inputMode == INPMTX_REAL_ENTRIES)
#define INPMTX_IS_COMPLEX_ENTRIES(inpmtx) \
   ((inpmtx)->inputMode == INPMTX_COMPLEX_ENTRIES)
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
 
   created -- 98jan28, cca
   -----------------------
*/
InpMtx *
InpMtx_new (
   void
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_setDefaultFields (
   InpMtx   *inpmtx
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98jan28, cca
   --------------------------------------------------
*/
void
InpMtx_clearData (
   InpMtx   *inpmtx
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98jan28, cca
   ------------------------------------------
*/
InpMtx *
InpMtx_free (
   InpMtx   *inpmtx
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in instance.c --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   returns coordinate type
 
   created -- 98jan28, cca
   -----------------------
*/
int
InpMtx_coordType (
   InpMtx  *inpmtx
) ;
/*
   -----------------------
   returns storage mode
 
   created -- 98jan28, cca
    ----------------------
*/
int
InpMtx_storageMode (
   InpMtx  *inpmtx
) ;
/*
   -----------------------
   returns input mode
 
   created -- 98jan28, cca
   -----------------------
*/
int
InpMtx_inputMode (
   InpMtx  *inpmtx
) ;
/*
   -----------------------
   returns inpmtx->mxnent
 
   created -- 98jan28, cca
   -----------------------
*/
int
InpMtx_maxnent (
   InpMtx  *inpmtx
) ;
/*
   -----------------------
   returns inpmtx->nent
 
   created -- 98jan28, cca
   -----------------------
*/
int
InpMtx_nent (
   InpMtx  *inpmtx
) ;
/*
   -------------------------
   returns inpmtx->mxnvector
 
   created -- 98jan28, cca
   -------------------------
*/
int
InpMtx_maxnvector (
   InpMtx  *inpmtx
) ;
/*
   -----------------------
   returns inpmtx->nvector
 
   created -- 98jan28, cca
   -----------------------
*/
int
InpMtx_nvector (
   InpMtx  *inpmtx
) ;
/*
   ------------------------------
   returns inpmtx->resizeMultiple
 
   created -- 98jan28, cca
   ------------------------------
*/
double
InpMtx_resizeMultiple (
   InpMtx  *inpmtx
) ;
/*
   ---------------------------------
   returns pointer to ivec1[] vector
 
   created -- 98jan28, cca
   ---------------------------------
*/
int *
InpMtx_ivec1 (
   InpMtx  *inpmtx
) ;
/*
   ---------------------------------
   returns pointer to ivec2[] vector
 
   created -- 98jan28, cca
   ---------------------------------
*/
int *
InpMtx_ivec2 (
   InpMtx  *inpmtx
) ;
/*
   --------------------------------
   returns pointer to dvec[] vector
 
   created -- 98jan28, cca
   --------------------------------
*/
double *
InpMtx_dvec (
   InpMtx  *inpmtx
) ;
/*
   ---------------------------------
   returns pointer to sizes[] vector
 
   created -- 98jan28, cca
   ---------------------------------
*/
int *
InpMtx_sizes (
   InpMtx  *inpmtx
) ;
/*
   ----------------------------------
   returns pointer to vecids[] vector
 
   created -- 98jan28, cca
   ----------------------------------
*/
int *
InpMtx_vecids (
   InpMtx  *inpmtx
) ;
/*
   -----------------------------------
   returns pointer to offsets[] vector
 
   created -- 98jan28, cca
   -----------------------------------
*/
int *
InpMtx_offsets (
   InpMtx  *inpmtx
) ;
/*
   ---------------------------------------
   retrieve requested vector
   set *pnent to # of entries
       *pindices to address of first index
 
   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_vector (
   InpMtx   *inpmtx,
   int       id,
   int       *pnent,
   int       **pindices
) ;
/*
   ---------------------------------------
   retrieve requested vector
   set *pnent to # of entries
       *pindices to address of first index
       *pentries to address of first entry
 
   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_realVector (
   InpMtx   *inpmtx,
   int       id,
   int       *pnent,
   int       **pindices,
   double    **pentries
) ;
/*
   ---------------------------------------
   retrieve requested vector
   set *pnent to # of entries
       *pindices to address of first index
       *pentries to address of first entry
 
   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_complexVector (
   InpMtx   *inpmtx,
   int       id,
   int       *pnent,
   int       **pindices,
   double    **pentries
) ;
/*
   --------------------------------------------------------------
   sets the maximum numnber of entries.  this methods resizes the
   ivec1[], ivece2[] and dvec[] vectors if newmaxnent != maxnent
 
   created -- 98jan28, cca
   --------------------------------------------------------------
*/
void
InpMtx_setMaxnent (
   InpMtx  *inpmtx,
   int      newmaxnent
) ;
/*
   ---------------------------------
   set the present number of entries
 
   created -- 98jan28, cca
   --------------------------------
*/
void
InpMtx_setNent (
   InpMtx   *inpmtx,
   int       newnent
) ;
/*
   --------------------------------------------------
   sets the maximum number of vectors.
   if newmaxnent != maxnent then this methods resizes
   the vecids[], sizes[] and offsets[] vectors
 
   created  -- 98jan28, cca
   --------------------------------------------------
*/
void
InpMtx_setMaxnvector (
   InpMtx  *inpmtx,
   int      newmaxnvector
) ;
/*
   ---------------------------------
   set the present number of vectors
 
   created -- 98jan28, cca
   ---------------------------------
*/
void
InpMtx_setNvector (
   InpMtx   *inpmtx,
   int       newnvector
) ;
/*
   ---------------------------
   sets inpmtx->resizeMultiple
 
   created -- 98jan28, cca
    ---------------------------
*/
void
InpMtx_setResizeMultiple (
   InpMtx   *inpmtx,
   double    resizeMultiple
) ;
/*
   --------------------------------------------
   sets coordType of InpMtx structure to
   allow user to define custom coordinate type.
   Note, new type must be > 3.
 
   created -- 98jan28, cca
    --------------------------------------------
*/
void
InpMtx_setCoordType (
   InpMtx  *inpmtx,
   int      type
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   initialize the object
 
   coordType -- coordinate type, input supported for types 1, 2 and 3
      1 -- row triples (i, j, a_{i,j})
      2 -- column triples (i, j, a_{j,i})
      3 -- chevron triples 
           (i, k, a_{i,i+k}) if k >= 0
           (i, k, a_{i-k,i}) if k < 0
           i is the chevron, k is the offset
      4 -- custom coordinate, e.g., one could store (I, k, a_{i,j})
          where I is the front where a_{i,j} will be assembled
           and k is the offset into the vector that holds the 
           entries in the front
   inputMode -- mode for input
      1 --> indices only
      2 --> indices and entries
   maxnent  -- upper bound on the number of entries,
      also equal to the workspace, so if the assembly includes
      overlapping data, give enough elbow room for efficiency.
   maxnvector -- upper bound on the number of vectors to be supported. 
      this may not be known ahead of time (e.g., the number of vectors 
      may be the number of fronts which is not known before the
      ordering is done and front tree constructed).
 
   created -- 98jan28, cca
   ------------------------------------------------------------------
*/
void
InpMtx_init (
  InpMtx   *inpmtx,
  int       coordType,
  int       inputMode,
  int       maxnent,
  int       maxnvector
) ;
/*
   --------------------------
   change the coordinate type
 
   created -- 98jan28, cca
   --------------------------
*/
void
InpMtx_changeCoordType (
   InpMtx   *inpmtx,
   int       newType
) ;
/*
   -----------------------
   change the storage mode
 
   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_changeStorageMode (
   InpMtx   *inpmtx,
   int       newMode
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in input.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------
   input a single entry in the matrix
 
   created -- 98jan28, cca
   ----------------------------------
*/
void
InpMtx_inputEntry (
   InpMtx   *inpmtx,
   int       row,
   int       col
) ;
/*
   ---------------------------------------
   input a single real entry in the matrix
 
   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_inputRealEntry (
   InpMtx   *inpmtx,
   int       row,
   int       col,
   double    value
) ;
/*
   ------------------------------------------
   input a single complex entry in the matrix
 
   created -- 98jan28, cca
   ------------------------------------------
*/
void
InpMtx_inputComplexEntry (
   InpMtx   *inpmtx,
   int       row,
   int       col,
   double    real,
   double    imag
) ;
/*
   ------------------------------
   input a real row in the matrix
 
   created -- 98jan28, cca
   ------------------------------
*/
void
InpMtx_inputRow (
   InpMtx   *inpmtx,
   int       row,
   int       rowsize,
   int       rowind[]
) ;
/*
   ------------------------------
   input a real row in the matrix
 
   created -- 98jan28, cca
   ------------------------------
*/
void
InpMtx_inputRealRow (
   InpMtx   *inpmtx,
   int       row,
   int       rowsize,
   int       rowind[],
   double    rowent[]
) ;
/*
   ---------------------------------
   input a complex row in the matrix
 
   created -- 98jan28, cca
   ---------------------------------
*/
void
InpMtx_inputComplexRow (
   InpMtx   *inpmtx,
   int       row,
   int       rowsize,
   int       rowind[],
   double    rowent[]
) ;
/*
   ----------------------------
   input a column in the matrix
 
   created -- 98jan28, cca
   ----------------------------
*/
void
InpMtx_inputColumn (
   InpMtx   *inpmtx,
   int       col,
   int       colsize,
   int       colind[]
) ;
/*
   ---------------------------------
   input a real column in the matrix
 
   created -- 98jan28, cca
   ---------------------------------
*/
void
InpMtx_inputRealColumn (
   InpMtx   *inpmtx,
   int       col,
   int       colsize,
   int       colind[],
   double    colent[]
) ;
/*
   ------------------------------------
   input a complex column in the matrix
 
   created -- 98jan28, cca
   ------------------------------------
*/
void
InpMtx_inputComplexColumn (
   InpMtx   *inpmtx,
   int       col,
   int       colsize,
   int       colind[],
   double    colent[]
) ;
/*
   -----------------------------
   input a chevron in the matrix
 
   created -- 98jan28, cca
   -----------------------------
*/
void
InpMtx_inputChevron (
   InpMtx   *inpmtx,
   int       chv,
   int       chvsize,
   int       chvind[]
) ;
/*
   -----------------------------
   input a chevron in the matrix
 
   created -- 98jan28, cca
   -----------------------------
*/
void
InpMtx_inputRealChevron (
   InpMtx   *inpmtx,
   int       chv,
   int       chvsize,
   int       chvind[],
   double    chvent[]
) ;
/*
   -----------------------------
   input a chevron in the matrix
 
   created -- 98jan28, cca
   -----------------------------
*/
void
InpMtx_inputComplexChevron (
   InpMtx   *inpmtx,
   int       chv,
   int       chvsize,
   int       chvind[],
   double    chvent[]
) ;
/*
   -----------------------
   input a matrix
 
   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_inputMatrix (
   InpMtx   *inpmtx,
   int       nrow,
   int       ncol,
   int       rowstride,
   int       colstride,
   int       rowind[],
   int       colind[]
) ;
/*
   -----------------------
   input a matrix
 
   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_inputRealMatrix (
   InpMtx   *inpmtx,
   int       nrow,
   int       ncol,
   int       rowstride,
   int       colstride,
   int       rowind[],
   int       colind[],
   double    mtxent[]
) ;
/*
   -----------------------
   input a matrix
 
   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_inputComplexMatrix (
   InpMtx   *inpmtx,
   int       nrow,
   int       ncol,
   int       rowstride,
   int       colstride,
   int       rowind[],
   int       colind[],
   double    mtxent[]
) ;
/*
   -------------------------------------------------------------
   input a number of (row,column, entry) triples into the matrix
 
   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
InpMtx_inputTriples (
   InpMtx   *inpmtx,
   int       ntriples,
   int       rowids[],
   int       colids[]
) ;
/*
   -------------------------------------------------------------
   input a number of (row,column, entry) triples into the matrix
 
   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
InpMtx_inputRealTriples (
   InpMtx   *inpmtx,
   int       ntriples,
   int       rowids[],
   int       colids[],
   double    entries[]
) ;
/*
   -------------------------------------------------------------
   input a number of (row,column, entry) triples into the matrix
 
   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
InpMtx_inputComplexTriples (
   InpMtx   *inpmtx,
   int       ntriples,
   int       rowids[],
   int       colids[],
   double    entries[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in permute.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   permute the entries

   created -- 96jul05, cca
   -----------------------
*/
void
InpMtx_permute (
   InpMtx   *inpmtx,
   int       rowOldToNew[],
   int       colOldToNew[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in extract.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------
   purpose -- to extract a submatrix, B = A(BrowsIV, BcolsIV)
 
   return values ---
      1 -- normal return
     -1 -- B is NULL
     -2 -- BcolsIV is NULL
     -3 -- BrowsIV is NULL
     -4 -- B is NULL
     -5 -- invalid input mode for A
     -6 -- invalid coordinate type for A
     -7 -- invalid symmetryflag
     -8 -- hermitian flag but not complex
     -9 -- msglvl > 0 and msgFile = NULL
 
   created -- 98oct15, cca
   ----------------------------------------------------------
*/
int
InpMtx_initFromSubmatrix (
   InpMtx   *B,
   InpMtx   *A,
   IV       *BrowsIV,
   IV       *BcolsIV,
   int      symmetryflag,
   int      msglvl,
   FILE     *msgFile
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------
   given the data is in raw triples,
   sort and compress the data
 
   created -- 98jan28, cca
   ---------------------------------
*/
void
InpMtx_sortAndCompress (
   InpMtx   *inpmtx
) ;
/*
   ----------------------------------------------------
   convert from sorted and compressed triples to vector
 
   created -- 98jan28, cca
   ----------------------------------------------------
*/
void
InpMtx_convertToVectors (
   InpMtx   *inpmtx
) ;
/*
   -------------------------
   drop off-diagonal entries
 
   created -- 98jan28, cca
   -------------------------
*/
void
InpMtx_dropOffdiagonalEntries (
   InpMtx   *inpmtx
) ;
/*
   ----------------------------------
   drop entries in the lower triangle
 
   created -- 98jan28, cca
   ----------------------------------
*/
void
InpMtx_dropLowerTriangle (
   InpMtx   *inpmtx
) ;
/*
   ----------------------------------
   drop entries in the upper triangle
 
   created -- 98jan28, cca
   ----------------------------------
*/
void
InpMtx_dropUpperTriangle (
   InpMtx   *inpmtx
) ;
/*
   -----------------------------------
   map entries into the lower triangle
 
   created -- 98jan28, cca
   -----------------------------------
*/
void
InpMtx_mapToLowerTriangle (
   InpMtx   *inpmtx
) ;
/*
   -----------------------------------
   map entries into the upper triangle
 
   created -- 98jan28, cca
   -----------------------------------
*/
void
InpMtx_mapToUpperTriangle (
   InpMtx   *inpmtx
) ;
/*
   -----------------------------------
   map entries into the upper triangle
   for a hermitian matrix
 
   created -- 98may15, cca
   -----------------------------------
*/
void
InpMtx_mapToUpperTriangleH (
   InpMtx   *inpmtx
) ;
/*
   -----------------------------------------------------------
   purpose -- compute the checksums of the indices and entries
 
      sums[0] = sum_{ii=0}^{nent} abs(ivec1[ii])
      sums[1] = sum_{ii=0}^{nent} abs(ivec2[ii])
      if real or complex entries then
         sums[2] = sum_{ii=0}^{nent} magnitudes of entries
      endif
 
   created -- 98may16, cca
   -----------------------------------------------------------
*/
void
InpMtx_checksums (
   InpMtx   *inpmtx,
   double   sums[]
) ;
/*
   ----------------------------------------------------------------
   purpose -- to create an InpMtx object filled with random entries
 
   input --
 
      mtx         -- matrix object, if NULL, it is created
      inputMode   -- input mode for the object,
                     indices only, real or complex entries
      coordType   -- coordinate type for the object,
                     by rows, by columns or by chevrons
      storageMode -- storage mode for the object,
                     raw data, sorted or by vectors
      nrow        -- # of rows
      ncol        -- # of columns
      symflag     -- symmetry flag for the matrix,
                     symmetric, hermitian or nonsymmetric
      nonzerodiag -- if 1, entries are placed on the diagonal
      nitem       -- # of items to be placed into the matrix
      seed        --  random number seed
 
   return value ---
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- bad input mode
     -3 -- bad coordinate type
     -4 -- bad storage mode
     -5 -- nrow or ncol <= 0
     -6 -- bad symmetry flag
     -7 -- hermitian matrix but not complex
     -8 -- symmetric or hermitian matrix but nrow != ncol
     -9 -- nitem < 0
   ----------------------------------------------------------------
*/
int
InpMtx_randomMatrix (
   InpMtx   *mtx,
   int      inputMode,
   int      coordType,
   int      storageMode,
   int      nrow,
   int      ncol,
   int      symflag,
   int      nonzerodiag,
   int      nitem,
   int      seed
) ;
/*
   ----------------------------------------------------
   determine the range of the matrix, 
   i.e., the minimum and maximum rows and columns
 
   if pmincol != NULL then *pmincol = minimum column id
   if pmaxcol != NULL then *pmaxcol = maximum column id
   if pminrow != NULL then *pminrow = minimum row id
   if pmaxrow != NULL then *pmaxrow = maximum row id
 
   return value ---
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- no entries in the matrix
     -3 -- invalid coordinate type
   
   created -- 98oct15, cca
   ----------------------------------------------------
*/
int
InpMtx_range (
   InpMtx   *mtx,
   int      *pmincol,
   int      *pmaxcol,
   int      *pminrow,
   int      *pmaxrow
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in fullAdj.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------------------------------
   purpose -- to return the full, symmetric adjacency IVL object
              for the graph of A + A^T
 
   created -- 98jan28, cca
   -------------------------------------------------------------
*/
IVL *
InpMtx_fullAdjacency (
   InpMtx   *inpmtx
) ;
/*
   -------------------------------------------------------------
   purpose -- to return the full, symmetric adjacency IVL object
              for the graph of (A + B) + (A + B)^T

   created -- 97nov05, cca
   -------------------------------------------------------------
*/
IVL *
InpMtx_fullAdjacency2 (
   InpMtx   *inpmtxA,
   InpMtx   *inpmtxB
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in adjForATA.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------
   return an IVL object that contains 
   the adjacency structure of A^TA.
 
   created -- 98jan28, cca
   ----------------------------------
*/
IVL *
InpMtx_adjForATA (
   InpMtx   *inpmtxA
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in profile.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------------------
   to fill xDV and yDV with a log10 profile of the magnitudes of
   the entries in the InpMtx object. tausmall and tau big provide
   cutoffs within which to examine the entries. pnsmall and pnbig 
   are address to hold the number of entries smaller than tausmall,

   and larger than taubig, respectively.
 
   created -- 97feb14, cca
   ------------------------------------------------------------------
*/
void
InpMtx_log10profile (
   InpMtx    *inpmtx,
   int        npts,
   DV         *xDV,
   DV         *yDV,
   double     tausmall,
   double     taubig,
   int        *pnzero,
   int        *pnsmall,
   int        *pnbig
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in gmmm.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
   where X and Y are DenseMtx objects, and X and Y
   must be column major.
 
   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- beta is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal
 
   created -- 98may02, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmmm (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A^T*X
 
   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- beta is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal
 
   created -- 98may02, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmmm_T (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A^H*X
 
   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- beta is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal
    -16 -- A, X and Y are real
 
   created -- 98may02, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmmm_H (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
              where A is symmetric
 
   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- beta is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal
 
   created -- 98nov06, cca
   --------------------------------------------------
*/
int
InpMtx_sym_gmmm (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
              where A is hermitian
 
   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- Y is NULL
     -6 -- type of Y is invalid
     -7 -- dimensions and strides of Y do not line up
     -8 -- entries of Y are NULL
     -9 -- beta is NULL
    -10 -- X is NULL
    -11 -- type of X is invalid
    -12 -- dimensions and strides of X do not line up
    -13 -- entries of X are NULL
    -14 -- types of A, X and Y are not identical
    -15 -- # of columns of X and Y are not equal
    -16 -- A, X and Y are real
 
   created -- 98nov06, cca
   --------------------------------------------------
*/
int
InpMtx_herm_gmmm (
   InpMtx     *A,
   double     beta[],
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in gmvm.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------
   purpose -- to compute y := beta*y + alpha*A*x
   where x and y are double vectors.

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL

   created -- 98nov14, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmvm (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) ;
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A^T*X

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL

   created -- 98may02, cca
   --------------------------------------------------
*/
int
InpMtx_nonsym_gmvm_T (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) ;
/*
   ---------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A^H*X

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL
    -10 -- A is real

   created -- 98may02, cca
   ---------------------------------------------------
*/
int
InpMtx_nonsym_gmvm_H (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) ;
/*
   ----------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
              where A is symmetric

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL

   created -- 98nov06, cca
   ----------------------------------------------------
*/
int
InpMtx_sym_gmvm (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) ;
/*
   --------------------------------------------------
   purpose -- to compute Y := beta*Y + alpha*A*X
              where A is hermitian

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- type of A is invalid
     -3 -- indices or entries of A are NULL
     -4 -- beta is NULL
     -5 -- ny <= 0
     -6 -- y is NULL
     -7 -- alpha is NULL
     -8 -- nx <= 0
     -9 -- x is NULL
    -10 -- A, X and Y are real

   created -- 98nov06, cca
   --------------------------------------------------
*/
int
InpMtx_herm_gmvm (
   InpMtx     *A,
   double     beta[],
   int        ny,
   double     y[],
   double     alpha[],
   int        nx,
   double     x[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in mmm.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_nonsym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^T*X
 
   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmm_T (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^H*X
 
   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmm_H (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
              where A is symmetric
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_sym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
              where A is hermitian
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_herm_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in mvmVector.c -------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_nonsym_mmmVector (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) ;
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^T*X
 
   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmmVector_T (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) ;
/*
   ------------------------------------------
   purpose -- to compute Y := Y + alpha*A^H*X
 
   created -- 98may28, cca
   ------------------------------------------
*/
void
InpMtx_nonsym_mmmVector_H (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
              where A is symmetric
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_sym_mmmVector (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) ;
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
              where A is hermitian
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
InpMtx_herm_mmmVector (
   InpMtx   *A,
   double   y[],
   double   alpha[],
   double   x[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in support.c ---------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   purpose -- 
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      nonsymmetric matrix.
 
      rowsupIV -- filled with row indices of y[] which will be updated.
      colsupIV -- filled with row indices of x[] which will be used.
 
   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportNonsym (
   InpMtx   *A,
   IV       *rowsupIV,
   IV       *colsupIV
) ;
/*
   --------------------------------------------------------------------
   purpose --
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      nonsymmetric matrix.
 
      rowsupIV -- filled with row indices of y[] which will be updated.
      colsupIV -- filled with row indices of x[] which will be used.
 
   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportNonsymT (
   InpMtx   *A,
   IV       *rowsupIV,
   IV       *colsupIV
) ;
/*
   --------------------------------------------------------------------
   purpose --
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      nonsymmetric matrix.
 
      rowsupIV -- filled with row indices of y[] which will be updated.
      colsupIV -- filled with row indices of x[] which will be used.
 
   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportNonsymH (
   InpMtx   *A,
   IV       *rowsupIV,
   IV       *colsupIV
) ;
/*
   --------------------------------------------------------------------
   purpose --
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      symmetric matrix.
 
      supIV -- filled with row indices of y[] which will be updated
               and row indices of x[] which will be used.
 
   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportSym (
   InpMtx   *A,
   IV       *supIV
) ;
/*
   --------------------------------------------------------------------
   purpose --
      this method is used to determine the support of this matrix
      for a matrix-vector multiply y[] = A * x[] when A is a
      Hermitian matrix.
 
      supIV -- filled with row indices of y[] which will be updated
               and row indices of x[] which will be used.
 
   created -- 98aug01, cca
   --------------------------------------------------------------------
*/
void
InpMtx_supportHerm (
   InpMtx   *A,
   IV       *supIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in map.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   --------------------------------------------------------------------
   purpose -- to map the coordinates of the entries of the matrix 
      into new coordinates. this method is used during the distributed
      matrix-vector multiply where a matrix local to a processor is
     mapped into a local coordinate system.
 
   (row,col) --> (rowmap[row],colmap[col])
 
   we check that row is in [0,nrow) and col is in [0,ncol),
   where nrow is the size of rowmapIV and ncol is the size of colmapIV.
 
   note, the storage mode is not changed. i.e., if the data is
   stored by vectors, it may be invalid after the indices have 
   been mapped. on the other hand, it may not, so it is the user's
   responsibility to reset the storage mode if necessary.
 
   created -- 98aug02, cca
   --------------------------------------------------------------------
*/
void
InpMtx_mapEntries (
   InpMtx   *inpmtx,
   IV       *rowmapIV,
   IV       *colmapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ---------------------------------------------------
   purpose -- to read in a InpMtx object from a file
 
   input --
 
      fn -- filename, must be *.inpmtxb or *.inpmtxf
 
   return value -- 1 if success, 0 if failure
 
   created -- 98jan28, cca
   --------------------------------------------------
*/
int
InpMtx_readFromFile ( 
   InpMtx   *inpmtx, 
   char      *fn 
) ;
/*
   --------------------------------------------------------
   purpose -- to read a InpMtx object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 98jan28, cca
   --------------------------------------------------------
*/
int
InpMtx_readFromFormattedFile (
   InpMtx   *inpmtx,
   FILE      *fp
) ;
/*
   ------------------------------------------------------
   purpose -- to read a InpMtx object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 98jan28, cca
   ------------------------------------------------------
*/
int
InpMtx_readFromBinaryFile (
   InpMtx   *inpmtx,
   FILE    *fp
) ;
/*
   ----------------------------------------------
   purpose -- to write a InpMtx object to a file
 
   input --
 
      fn -- filename
        *.inpmtxb -- binary
        *.inpmtxf -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan28, cca
   ----------------------------------------------
*/
int
InpMtx_writeToFile (
   InpMtx   *inpmtx,
   char    *fn
) ;
/*
   ------------------------------------------------------
   purpose -- to write a InpMtx object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan28, cca
   ------------------------------------------------------
*/
int
InpMtx_writeToFormattedFile (
   InpMtx   *inpmtx,
   FILE    *fp
) ;
/*
   ---------------------------------------------------
   purpose -- to write a InpMtx object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan28, cca
   ---------------------------------------------------
*/
int
InpMtx_writeToBinaryFile (
   InpMtx    *inpmtx,
   FILE   *fp
) ;
/*
   ----------------------------------------------------
   purpose -- to write a InpMtx object for a human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan28, cca
   ----------------------------------------------------
*/
int
InpMtx_writeForHumanEye (
   InpMtx    *inpmtx,
   FILE   *fp
) ;
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the InpMtx object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan28, cca
   -----------------------------------------------------------
*/
int
InpMtx_writeStats (
   InpMtx    *inpmtx,
   FILE   *fp
) ;
/*
   ----------------------------------------------------
   purpose -- to write a InpMtx object to a matlab file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98jan28, cca
   ----------------------------------------------------
*/
int
InpMtx_writeForMatlab ( 
   InpMtx   *inpmtx, 
   char     *mtxname,
   FILE     *fp 
) ;
/*
   ----------------------------------------------------------------
  purpose -- to read in a InpMtx object from a Harwell-Boeing file
 
   input --
 
      fn -- filename
 
   return value -- 1 if success, 0 if failure
 
   created -- 98sep11, cca
   ---------------------------------------------------------------
*/
int
InpMtx_readFromHBfile (
   InpMtx   *inpmtx,
   char     *fn
) ;
/*--------------------------------------------------------------------*/
