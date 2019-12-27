/*  DV.h  */

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to compute y[*] := y[*] + x[*]

   created -- 95sep22, cca
   -----------------------------------------
*/
void
DVadd ( 
   int      size, 
   double   y[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to compute y[*] := y[*] + alpha * x[*]

   created -- 95sep22, cca
   -------------------------------------------------
*/
void
DVaxpy ( 
   int      size, 
   double   y[], 
   double   alpha, 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to compute y[index[*]] := y[index[*]] + alpha * x[*]
 
   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
void
DVaxpyi ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   alpha, 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------
   purpose -- to copy y[*] := x[*]

   created -- 95sep22, cca
   -------------------------------
*/
void
DVcopy ( 
   int      size, 
   double   y[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- given the pair of arrays (x1[],y1[]), 
              create a pair of arrays (x2[],y2[]) whose
              entries are pairwise chosen from (x1[],y1[])
              and whose distribution is an approximation.

   return value -- the size of the (x2[],y2[]) arrays

   created -- 95sep22, cca
   -------------------------------------------------------
*/
int
DVcompress ( 
   int      size1, 
   double   x1[], 
   double   y1[],
   int      size2, 
   double   x2[], 
   double   y2[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- to copy sum{y[*] * x[*]}

   created -- 95sep22, cca
   -----------------------------------
*/
double
DVdot ( 
   int      size, 
   double   y[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to perform a indexed dot product
 
   sum = sum_k y[index[k]]*x[k]
 
   where y and x are real
 
   created -- 98may02, cca
   --------------------------------------------
*/
double
DVdoti ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[]
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to fill a double vector with a value

   created -- 95sep22, cca
   -----------------------------------------------
*/
void
DVfill ( 
   int      size, 
   double   y[], 
   double   dval 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to print out a double vector

   created -- 95sep22, cca
   -----------------------------------------
*/
void
DVfprintf ( 
   FILE     *fp, 
   int      size, 
   double   dvec[]
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to free storage taken by a double vector.
              note : dvec[] must have been returned by DVinit.

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
DVfree ( 
   double   dvec[] 
) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to read in a double vector
              return value -- # of entries read

   created -- 95sep22, cca
   --------------------------------------------
*/
int
DVfscanf ( 
   FILE     *fp, 
   int      size, 
   double   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- to gather y[*] = x[index[*]]

   created -- 95sep22, cca
   ---------------------------------------
*/
void
DVgather ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to gather add y[*] += x[index[*]] and zero x[index[*]]

   created -- 95sep22, cca
   -----------------------------------------------------------------
*/
void
DVgatherAddZero ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to gather y[*] = x[index[*]] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------
*/
void
DVgatherZero ( 
   int      size, 
   double   y[], 
   double   x[], 
   int      index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- allocate a double array with size entries
              and fill with value dval
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
double *
DVinit ( 
   int      size, 
   double   dval 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- allocate a double array with size entries
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
double *
DVinit2 ( 
   int   size 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------
   purpose -- to permute a vector
              y[index[*]] := y[*]

   created -- 95sep22, cca
   ------------------------------
*/
void
DVinvPerm ( 
   int      size, 
   double   y[], 
   int      index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to return the first entry of maximum value,
              *ploc contains its index 

   created -- 95sep22, cca
   ------------------------------------------------------
*/
double
DVmax ( 
   int      size, 
   double   y[], 
   int      *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of maximum absolute value,
              *ploc contains its index 

   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
double
DVmaxabs ( 
   int      size, 
   double   y[], 
   int      *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to return the first entry of minimum value,
              *ploc contains its index 

   created -- 95sep22, cca
   ------------------------------------------------------
*/
double
DVmin ( 
   int      size, 
   double   y[], 
   int      *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of minimum absolute value,
              *ploc contains its index 

   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
double
DVminabs ( 
   int      size, 
   double   y[], 
   int      *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------
   purpose -- to permute a vector
              y[*] := y[index[*]]

   created -- 95sep22, cca
   ------------------------------
*/
void
DVperm ( 
   int      size, 
   double   y[], 
   int      index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to fill a double vector with a ramp function

   created -- 95sep22, cca
   -------------------------------------------------------
*/
void
DVramp ( 
   int      size, 
   double   y[], 
   double   start, 
   double   inc 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to compute y[*] := y[*] - x[*]

   created -- 95sep22, cca
   -----------------------------------------
*/
void
DVsub ( 
   int      size, 
   double   y[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to scale a double vector by alpha

   created -- 95sep22, cca
   --------------------------------------------
*/
void
DVscale ( 
   int      size, 
   double   y[], 
   double   alpha 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to scatter y[index[*]] = x[*]

   created -- 95sep22, cca
   ----------------------------------------
*/
void
DVscatter ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to scatter add y[index[*]] += x[*]

   created -- 96aug17, cca
   ---------------------------------------------
*/
void
DVscatterAdd ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to scatter add y[index[*]] += x[*] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
DVscatterAddZero ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to scatter y[index[*]] = x[*] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------
*/
void
DVscatterZero ( 
   int      size, 
   double   y[], 
   int      index[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to return the sum of a double vector

   created -- 95sep22, cca
   -----------------------------------------------
*/
double
DVsum ( 
   int      size, 
   double   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return the sum of the absolute values 
              of the entries in a double vector

   created -- 95sep22, cca
   ---------------------------------------------------
*/
double
DVsumabs ( 
   int      size, 
   double   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- to swap y[*] and x[*]

   created -- 95sep22, cca
   --------------------------------
*/
void
DVswap ( 
   int      size, 
   double   y[], 
   double   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   purpose -- to zero a double vector

   created -- 95sep22, cca
   ----------------------------------
*/
void
DVzero ( 
   int      size, 
   double   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to permute an integer vector,
              procedure uses srand48() and drand48()

   input --

      size -- size of the vector
      y    -- vector to be permuted
      seed -- seed for the random number generator
              if seed <= 0, simple return

   created -- 95sep22, cca
   -------------------------------------------------
*/
void
DVshuffle ( 
   int      size, 
   double   y[], 
   int      seed 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to scale a double vector by a 2x2 matrix
 
   [ x ] := [ a b ] [ x ]
   [ y ]    [ c d ] [ y ]
 
   created -- 98jan23, cca
   ---------------------------------------------------
*/
void
DVscale2 (
   int      size,
   double   x[],
   double   y[],
   double   a,
   double   b,
   double   c,
   double   d
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to perform a axpy with two vectors
 
   z := z + a*x + b*y
 
   where y and x are double vectors
 
   created -- 98jan23, cca
   --------------------------------------------
*/
void
DVaxpy2 (
   int      size,
   double   z[],
   double   a,
   double   x[],
   double   b,
   double   y[]
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]
      sums[3] = row1[*] * col0[*]
      sums[4] = row1[*] * col1[*]
      sums[5] = row1[*] * col2[*]
      sums[6] = row2[*] * col0[*]
      sums[7] = row2[*] * col1[*]
      sums[8] = row2[*] * col2[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot33 (
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    col1[],
   double    col2[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]
      sums[3] = row1[*] * col0[*]
      sums[4] = row1[*] * col1[*]
      sums[5] = row1[*] * col2[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot23 (
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    col1[],
   double    col2[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot13 (
   int       n,
   double    row0[],
   double    col0[],
   double    col1[],
   double    col2[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row1[*] * col0[*]
      sums[3] = row1[*] * col1[*]
      sums[4] = row2[*] * col0[*]
      sums[5] = row2[*] * col1[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot32 (
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    col1[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row1[*] * col0[*]
      sums[3] = row1[*] * col1[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot22 (
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    col1[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot12 (
   int       n,
   double    row0[],
   double    col0[],
   double    col1[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[1] = row1[*] * col0[*]
      sums[2] = row2[*] * col0[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot31 (
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product
 
      sums[0] = row0[*] * col0[*]
      sums[2] = row1[*] * col0[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot21 (
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    sums[]
) ;
/*
   -----------------------------------------
   purpose -- compute a single dot product
 
      sums[0] = row0[*] * col0[*]
 
   created -- 98may02, cca
   -----------------------------------------
*/
void
DVdot11 (
   int       n,
   double    row0[],
   double    col0[],
   double    sums[]
) ;
