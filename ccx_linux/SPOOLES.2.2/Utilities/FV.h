/*  FV.h  */

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to compute y[*] := y[*] + x[*]

   created -- 95sep22, cca
   -----------------------------------------
*/
void
FVadd ( 
   int     size, 
   float   y[], 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to compute y[*] := y[*] + alpha * x[*]

   created -- 95sep22, cca
   -------------------------------------------------
*/
void
FVaxpy ( 
   int     size, 
   float   y[], 
   float   alpha, 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to compute y[index[*]] := y[index[*]] + alpha * x[*]
 
   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
void
FVaxpyi ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   alpha, 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------
   purpose -- to copy y[*] := x[*]

   created -- 95sep22, cca
   -------------------------------
*/
void
FVcopy ( 
   int     size, 
   float   y[], 
   float   x[] 
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
FVcompress ( 
   int     size1, 
   float   x1[], 
   float   y1[],
   int     size2, 
   float   x2[], 
   float   y2[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- to copy sum{y[*] * x[*]}

   created -- 95sep22, cca
   -----------------------------------
*/
float
FVdot ( 
   int     size, 
   float   y[], 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to fill a float vector with a value

   created -- 95sep22, cca
   -----------------------------------------------
*/
void
FVfill ( 
   int     size, 
   float   y[], 
   float   dval 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to print out a float vector

   created -- 95sep22, cca
   -----------------------------------------
*/
void
FVfprintf ( 
   FILE    *fp, 
   int     size, 
   float   y[]
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to free storage taken by a float vector.
              note : y[] must have been returned by FVinit.

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
FVfree ( 
   float   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to read in a float vector
              return value -- # of entries read

   created -- 95sep22, cca
   --------------------------------------------
*/
int
FVfscanf ( 
   FILE    *fp, 
   int     size, 
   float   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- to gather y[*] = x[index[*]]

   created -- 95sep22, cca
   ---------------------------------------
*/
void
FVgather ( 
   int     size, 
   float   y[], 
   float   x[], 
   int     index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to gather add y[*] += x[index[*]] and zero x[index[*]]

   created -- 95sep22, cca
   -----------------------------------------------------------------
*/
void
FVgatherAddZero ( 
   int     size, 
   float   y[], 
   float   x[], 
   int     index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to gather y[*] = x[index[*]] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------
*/
void
FVgatherZero ( 
   int     size, 
   float   y[], 
   float   x[], 
   int     index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- allocate a float array with size entries
              and fill with value dval
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
float *
FVinit ( 
   int     size, 
   float   dval 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- allocate a float array with size entries
   
   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------------
*/
float *
FVinit2 ( 
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
FVinvPerm ( 
   int     size, 
   float   y[], 
   int     index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to return the first entry of maximum value,
              *ploc contains its index 

   created -- 95sep22, cca
   ------------------------------------------------------
*/
float
FVmax ( 
   int     size, 
   float   y[], 
   int     *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of maximum absolute value,
              *ploc contains its index 

   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
float
FVmaxabs ( 
   int     size, 
   float   y[], 
   int     *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to return the first entry of minimum value,
              *ploc contains its index 

   created -- 95sep22, cca
   ------------------------------------------------------
*/
float
FVmin ( 
   int     size, 
   float   y[], 
   int     *ploc 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to return the first entry of minimum absolute value,
              *ploc contains its index 

   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
float
FVminabs ( 
   int     size, 
   float   y[], 
   int     *ploc 
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
FVperm ( 
   int     size, 
   float   y[], 
   int     index[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to fill a float vector with a ramp function

   created -- 95sep22, cca
   -------------------------------------------------------
*/
void
FVramp ( 
   int     size, 
   float   y[], 
   float   start, 
   float   inc 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to compute y[*] := y[*] - x[*]

   created -- 95sep22, cca
   -----------------------------------------
*/
void
FVsub ( 
   int     size, 
   float   y[], 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to scale a float vector by alpha

   created -- 95sep22, cca
   --------------------------------------------
*/
void
FVscale ( 
   int     size, 
   float   y[], 
   float   alpha 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to scatter y[index[*]] = x[*]

   created -- 95sep22, cca
   ----------------------------------------
*/
void
FVscatter ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to scatter add y[index[*]] += x[*] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------------
*/
void
FVscatterAddZero ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to scatter y[index[*]] = x[*] and zero x[*]

   created -- 95sep22, cca
   -----------------------------------------------------
*/
void
FVscatterZero ( 
   int     size, 
   float   y[], 
   int     index[], 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to return the sum of a float vector

   created -- 95sep22, cca
   -----------------------------------------------
*/
float
FVsum ( 
   int     size, 
   float   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to return the sum of the absolute values 
              of the entries in a float vector

   created -- 95sep22, cca
   ---------------------------------------------------
*/
float
FVsumabs ( 
   int     size, 
   float   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------
   purpose -- to swap y[*] and x[*]

   created -- 95sep22, cca
   --------------------------------
*/
void
FVswap ( 
   int     size, 
   float   y[], 
   float   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   purpose -- to zero a float vector

   created -- 95sep22, cca
   ----------------------------------
*/
void
FVzero ( 
   int     size, 
   float   y[] 
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
FVshuffle ( 
   int     size, 
   float   y[], 
   int     seed 
) ;
/*--------------------------------------------------------------------*/
