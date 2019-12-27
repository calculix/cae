/*  CV.h  */

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   prototype definitions for character vectors.

   note, our use of character vectors is different from strings.
   most of their application is as marking vectors.
   -------------------------------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- allocate a char array with size entries 
              and fill with value c 

   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------
*/
char *
CVinit ( 
   int    size, 
   char   c 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- allocate a char array with size entries 

   return value -- a pointer to the start of the array

   created : 95sep22, cca
   ---------------------------------------------------
*/
char *
CVinit2 ( 
   int   size 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------
   purpose -- to copy y[*] := x[*]

   created -- 95sep22, cca
   -------------------------------
*/
void
CVcopy ( 
   int    size, 
   char   y[], 
   char   x[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------
   purpose -- to fill y[*] := c

   created -- 95sep22, cca
   ----------------------------
*/
void
CVfill ( 
   int    size, 
   char   y[], 
   char   c 
) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- to print out a char vector
      the vector starts on a new line 
      and takes 80 characters per line

   created -- 95sep22, cca
   -------------------------------------
*/
void
CVfprintf ( 
   FILE   *fp, 
   int    size, 
   char   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- to write out a char vector with eighty column lines

   input --

      fp     -- file pointer, must be formatted and write access
      size   -- length of the vector
      y[]    -- char vector
      column -- present column
      pierr  -- pointer to int to hold return value,
                should be 1 if any print was successful,
                if fprintf() failed, then ierr = -1
  
   return value -- present column

   created -- 95sep22, cca
   ------------------------------------------------------------------
*/
int
CVfp80 ( 
   FILE   *fp, 
   int    size, 
   char   y[], 
   int    column,
   int    *pierr
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to release storage for a character vector,
              note, should have been created using CVinit

   created -- 95sep22, cca
   ------------------------------------------------------
*/
void
CVfree ( 
   char y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- to read in a char vector
  
   return value -- number of characters read

   created -- 95sep22, cca
   -----------------------------------------
*/
int
CVfscanf ( 
   FILE   *fp, 
   int    size, 
   char   y[] 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 y2 ]^T [ x0 x1 x2]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 y2 ]^T [ x0 x1 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 y2 ]^T [ x0 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 ]^T [ x0 x1 x2]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU23 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 ]^T [ x0 x1 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU22 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 ]^T [ x0 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU21 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 ]^T [ x0 x1 x2]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU13 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 ]^T [ x0 x1 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU12 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 ]^T [ x0 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotU11 (
   int      n,
   double   y0[],
   double   x0[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 y2 ]^H [ x0 x1 x2]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 y2 ]^H [ x0 x1 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 y2 ]^H [ x0 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 ]^H [ x0 x1 x2]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC23 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 ]^H [ x0 x1 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC22 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 y1 ]^H [ x0 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC21 (
   int      n,
   double   y0[],
   double   y1[],
   double   x0[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 ]^H [ x0 x1 x2]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC13 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   x2[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 ]^H [ x0 x1 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC12 (
   int      n,
   double   y0[],
   double   x0[],
   double   x1[],
   double   sums[]
) ;
/*
   ----------------------------------------------
   purpose -- to compute the multiple dot product
 
     [ y0 ]^H [ x0 ]
 
   created -- 98apr17, cca
   ----------------------------------------------
*/
void
ZVdotC11 (
   int      n,
   double   y0[],
   double   x0[],
   double   sums[]
) ;
/*--------------------------------------------------------------------*/
