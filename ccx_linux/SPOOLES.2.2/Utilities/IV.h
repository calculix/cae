/*  IV.h  */

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- given the pair of arrays (x1[],y1[]), 
              create a pair of arrays (x2[],y2[]) whose
              entries are pairwise chosen from (x1[],y1[])
              and whose distribution is an approximation.

   return value -- the size of the (x2[],y2[]) arrays
   -------------------------------------------------------
*/
int
IVcompress ( 
   int   size1, 
   int   x1[], 
   int   y1[],
   int   size2,  
   int   x2[],  
   int   y2[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------
   purpose -- to copy y[*] := x[*]
   -------------------------------
*/
void
IVcopy ( 
   int   size, 
   int   y[], 
   int   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to fill a int vector with a value
   -----------------------------------------------
*/
void
IVfill ( 
   int   size, 
   int   y[], 
   int   ival 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   purpose -- to print out a int vector
   ------------------------------------
*/
void
IVfprintf ( 
   FILE   *fp, 
   int    size, 
   int    y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- to write out an integer vector with eighty column lines
 
   input --
 
      fp     -- file pointer, must be formatted and write access
      size   -- length of the vector
      y      -- integer vector
      column -- present column
      pierr  -- pointer to int to hold return value, 
                should be 1 if any print was successful,
                if fprintf() failed, then ierr = -1
  
   return value -- present column
 
   created -- 95sep22, cca
   mods    -- 95sep29, cca
      added ierr argument to handle error returns from fprintf()
   -------------------------------------------------------------------
*/
int
IVfp80 ( 
   FILE   *fp, 
   int    size, 
   int    y[], 
   int    column,
   int    *pierr
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- to free the storage taken by an integer vector.
              note : ivec must have been returned by IVinit
   ----------------------------------------------------------
*/
void
IVfree ( 
   int y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- to read in an int vector
   -----------------------------------
*/
int
IVfscanf ( 
   FILE   *fp, 
   int    size, 
   int    y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- to gather y[*] = x[index[*]]
   ---------------------------------------
*/
void
IVgather ( 
   int   size, 
   int   y[], 
   int   x[], 
   int   index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- allocate a int array with size entries
              and fill with value dval
   
   return value -- a pointer to the start of the array
 
   created : 95sep22, cca
   ---------------------------------------------------
*/
int *
IVinit ( 
   int   size, 
   int   ival 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- allocate a int array with size entries
   
   return value -- a pointer to the start of the array
 
   created : 95sep22, cca
   ---------------------------------------------------
*/
int *
IVinit2 ( 
   int   size 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- allocate an int array and fill with the inverse of
              y[]. note, y[] must be a permutation vector.

   created : 95sep22, cca
   -------------------------------------------------------------
*/
int *
IVinverse ( 
   int   size, 
   int   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------
   purpose -- to permute a vector
              y[index[*]] := y[*]
   ------------------------------
*/
void
IVinvPerm ( 
   int   size, 
   int   y[], 
   int   index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to return the first entry of maximum size,
              *ploc contains its index 
   -----------------------------------------------------
*/
int
IVmax ( 
   int   size, 
   int   y[], 
   int   *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of maximum absolute value,
              *ploc contains its index
 
   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
int
IVmaxabs (
   int   size,
   int   y[],
   int   *ploc
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to return the first entry of minimum value,
              *ploc contains its index
 
   created -- 95sep22, cca
   ------------------------------------------------------
*/
int
IVmin (
   int   size,
   int   y[],
   int   *ploc
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of minimum absolute value,
              *ploc contains its index
 
   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
int
IVminabs (
   int   size,
   int   y[],
   int   *ploc
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------
   purpose -- to permute a vector
              y[*] := y[index[*]]
   ------------------------------
*/
void
IVperm ( 
   int   size, 
   int   y[], 
   int   index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to fill a int vector with a ramp function
   ----------------------------------------------------
*/
void
IVramp ( 
   int   size, 
   int   y[], 
   int   start, 
   int   inc 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to scatter y[index[*]] = x[*]
   ----------------------------------------
*/
void
IVscatter ( 
   int   size, 
   int   y[], 
   int   index[], 
   int   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to return the sum of a int vector
   -----------------------------------------------
*/
int
IVsum ( 
   int   size, 
   int   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return the sum of the absolute values 
              of the entries in an int vector
   ---------------------------------------------------
*/
int
IVsumabs ( 
   int   size, 
   int   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- to swap y[*] and x[*]
   --------------------------------
*/
void
IVswap ( 
   int   size, 
   int   y[], 
   int   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   purpose -- to zero a int vector
   ----------------------------------
*/
void
IVzero ( 
   int   size, 
   int   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to permute an integer vector,
              procedure uses srand48() and drand48()

   input --

      size -- size of the vector
      ivec -- vector to be permuted
      seed -- seed for the random number generator
              if seed <= 0, simple return
   -------------------------------------------------
*/
void
IVshuffle ( 
   int   size, 
   int   y[], 
   int   seed 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   locate an instance of target in the vector ivec[size].
   we assume that ivec[] is sorted in ascending order
   so we can use binary search to locate target.
 
   return value --
      -1  -- if target not found in ivec[]
      loc -- where target = ivec[loc]
 
   created -- 96may27, cca
   ------------------------------------------------------
*/
int
IVlocateViaBinarySearch (
   int   size,
   int   ivec[],
   int   target
) ;
/*--------------------------------------------------------------------*/
