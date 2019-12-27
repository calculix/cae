/*  misc.h.c  */

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to print out a list of entries

   input --

      length  -- number of entries
      entries -- entries array
      fp      -- file pointer
   -----------------------------------------
*/
void
fp_dvec ( 
   int      length, 
   double   entries[], 
   FILE     *fp 
) ;
/*
   -----------------------------------------
   purpose -- to print out a list of entries

   input --

      id      -- id of the entries
      length  -- number of entries
      entries -- entries array
      fp      -- file pointer
   -----------------------------------------
*/
void
fp_entries ( 
   int      id, 
   int      length, 
   double   entries[], 
   FILE     *fp 
) ;
/*
   ------------------------------------------------------------------
   purpose -- to write out an integer vector with eighty column lines

   input --

      length -- length of the vector
      ivec   -- integer vector
      column -- present column
      fp     -- file pointer, must be formatted and write access
  
   return value -- present column
   ------------------------------------------------------------------
*/
int
fp_i80 ( 
   int    length, 
   int    ivec[], 
   int    column, 
   FILE   *fp 
) ;
/*
   ------------------------------------------------------------------
   purpose -- to print out a list of entries given by an indexed list

   input --

      id      -- id of the entries
      length  -- number of entries
      indices -- indexed list
      entries -- entries array
      fp      -- file pointer
   ------------------------------------------------------------------
*/
void
fp_i_entries ( 
   int      id, 
   int      length, 
   int      indices[], 
   double   entries[], 
   FILE     *fp 
) ;
/*
   ------------------------------------------------------------------
   purpose -- to print out a list of entries given by an indexed list
              print indices[index[*]]

   input --

      id      -- id of the entries
      length  -- number of entries
      indices -- indices 
      index   -- index array
      fp      -- file pointer
   ------------------------------------------------------------------
*/
void
fp_i_indices ( 
   int    id, 
   int    length, 
   int    indices[], 
   int    index[], 
   FILE   *fp 
) ;
/*
   -----------------------------------------
   purpose -- to print out a list of indices

   input --

      id      -- id of the indices
      length  -- number of indices
      indices -- index array
      fp      -- file pointer
   -----------------------------------------
*/
void
fp_indices ( 
   int    id, 
   int    length, 
   int    indices[], 
   FILE   *fp 
) ;
/*
   -----------------------------------------
   purpose -- to print out an integer vector

   input --

      length  -- number of indices
      indices -- index array
      fp      -- file pointer
   -----------------------------------------
*/
void
fp_ivec ( 
   int    length, 
   int    indices[], 
   FILE   *fp 
) ;
/*
   -----------------------------------------
   purpose -- to print out a long vector

   input --

      length    -- number of indices
      l_indices -- index array
      fp        -- file pointer
   -----------------------------------------
*/
void
fp_lvec ( 
   int    length, 
   long   l_indices[], 
   FILE   *fp 
) ;
/*
   ------------------------------------------
   purpose -- to print out a character vector

   input --

      length    -- number of indices
      l_indices -- index array
      fp        -- file pointer
   ------------------------------------------
*/
void
fp_cvec ( 
   int   length, 
   char   cvec[], 
   FILE   *fp 
) ;

/*--------------------------------------------------------------------*/

void   ND0 ( int n1, int n2, int n3, int new_to_old[], 
            int west, int east, int south, int north,
            int bottom, int top ) ;
void   fp2DGrid ( int n1, int n2, int ivec[], FILE *fp ) ;
void   fp3DGrid ( int n1, int n2, int n3, int ivec[], FILE *fp ) ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to sort int/double vectors in ascending order
              using a bubble sort algorithm
   --------------------------------------------------------
*/
void
IDbbsortUp ( 
   int      size, 
   int      ivec[], 
   double   dvec[] 
) ;
/*
   ---------------------------------------------------------
   purpose -- to sort int/double vectors in descending order
              using a bubble sort algorithm
   ---------------------------------------------------------
*/
void
IDbbsortDown ( 
   int      size, 
   int      ivec[], 
   double   dvec[] 
) ;
/*
   --------------------------------------------------------
   purpose -- to sort int/double vectors in ascending order
              using a bubble sort algorithm
   --------------------------------------------------------
*/
void
IDDbbsortUp ( 
   int      size, 
   int      ivec[], 
   double   dvec1[],
   double   dvec2[]
) ;
/*
   ---------------------------------------------------------
   purpose -- to sort int/double vectors in descending order
              using a bubble sort algorithm
   ---------------------------------------------------------
*/
void
IDDbbsortDown ( 
   int      size, 
   int      ivec[], 
   double   dvec1[],
   double   dvec2[]
) ;
/*--------------------------------------------------------------------*/
