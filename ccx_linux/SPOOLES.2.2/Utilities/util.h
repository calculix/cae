/*  util.h  */

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   sort and count observations of two integers.

   input
      nobs -- # of observations
      x[]  -- vector of first components of observations
      y[]  -- vector of second components of observations

   output
      pndistinct -- pointer to hold # of distinct observations
      pxdistinct -- pointer to hold first components of the
                    distinct observations
      pydistinct -- pointer to hold second components of the
                    distinct observations
      pcounts    -- pointer to hold counts of the distinct observations

   created -- 95sep30, cca
   --------------------------------------------------------------------
*/
void
IV2sortAndCount (
   int   nobs,
   int   x[],
   int   y[],
   int   *pndistinct,
   int   **pxdistinct,
   int   **pydistinct,
   int   **pcounts
) ;
/*--------------------------------------------------------------------*/
