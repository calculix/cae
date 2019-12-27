/*  CV.c  */

#include "../Utilities.h"

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
) {
char   *y ;

if ( size <= 0 ) {
   y = NULL ;
} else {
   y = CVinit2(size) ;
   CVfill(size, y, c) ;
}
return(y) ; }

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
) {
char   *y ;
if ( size <= 0 ) {
   y = NULL ;
} else {
   ALLOCATE(y, char, size) ;
}
return(y) ; }

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
) {
if ( size <= 0 ) {
   return ;
} else if ( y == NULL || x == NULL ) {
   fprintf(stderr, 
           "\n fatal error in CVcopy, size = %d, y = %p, x = %p\n",
           size, y, x) ;
   exit(0) ;
} else { 
   int   i ;
   for ( i = 0 ; i < size ; i++ ) {
      y[i] = x[i] ; 
   }
}
return ; }

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
) {
if ( size <= 0 ) {
   return ;
} else if ( y == NULL ) {
   fprintf(stderr, 
           "\n fatal error in CVfill, size = %d, y = %p\n",
           size, y) ;
   exit(0) ;
} else { 
   int   i ;
   for ( i = 0 ; i < size ; i++ ) {
      y[i] = c ; 
   }
}
return ; }
   
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
) {
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, 
         "\n fatal error in CVfprintf, fp = %p, size = %d, y = %p\n",
         fp, size, y) ;
      exit(0) ;
   } else {
      int    i ;
      for ( i = 0 ; i < size ; i++ ) {
         if ( i % 80 == 0 ) fprintf(fp, "\n") ;
         fprintf(fp, "%c", y[i]) ; 
      }
   }
}
return ; }

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
) {
*pierr = 1 ;
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in CVfp80"
              "\n fp = %p, size = %d, y = %p, column = %d\n",
              fp, size, y, column) ;
      exit(0) ;
   } else {
      int    i ;
      for ( i = 0 ; i < size ; i++ ) {
         if ( ++column >= 80 ) {
            fprintf(fp, "\n") ;
            column = 0 ; 
         }
         if ( (*pierr = fprintf(fp, " %c", y[i])) < 0 ) {
/*
           --------------------
           error with fprintf()
           --------------------
*/
           break ;
         }
      }
   }
}
return(column) ; }

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
) {
if ( y != NULL ) {
   FREE(y) ;
}
return ; }

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
) {
int    i = 0 ;
if ( fp != NULL && size > 0 ) {
   if ( y == NULL ) {
      fprintf(stderr, "\n fatal error in CVfscanf"
              "\n fp = %p, size = %d, y = %p\n", fp, size, y) ;
      exit(0) ;
   } else {
      for ( i = 0 ; i < size ; i++ ) {
         if ( fscanf(fp, "%c", y + i) != 1 ) {
            break ; 
         } 
      }
   }
}
return(i) ; }

/*--------------------------------------------------------------------*/
