/*  mkGridGraph.c  */

#include "../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   generate a Graph object for a grid

   created -- 96mar15, cca
   -------------------------------------------------
*/
{
char     *outFileName ;
double   t1, t2 ;
FILE     *msgFile ;
Graph    *graph ;
int      count, msglvl, i, ij, ijk, imax, imin, inow, 
         j, jmax, jmin, jnow, k, kmax, kmin, know, 
         nvtx, n1, n2, n3, rc, stencil ;
int      *list ;
IVL      *adjIVL ;

if ( argc != 8 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile stencil n1 n2 n3 outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    stencil  -- type of stencil, 5, 9, 27 or 13"
      "\n    n1       -- # of grid points in first dimension"
      "\n    n2       -- # of grid points in second dimension"
      "\n    n3       -- # of grid points in third dimension"
      "\n    outFile  -- output file, must be *.graphf or *.graphb"
      "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
stencil = atoi(argv[3]) ;
switch ( stencil ) {
case  5 :
case  7 :
case  9 :
case 13 :
case 27 :
   break ;
default :
   fprintf(stderr, "\n fatal error in mkGridGraph"
           "\n stencil = %d, must be 5, 7, 9, 13 or 27\n", stencil) ;
   exit(-1) ;
}
n1 = atoi(argv[4]) ;
n2 = atoi(argv[5]) ;
n3 = atoi(argv[6]) ;
if ( n1 < 1 || n2 < 1 || n3 < 1 ) {
   fprintf(stderr, "\n fatal error in mkGridGraph"
           "\n n1 = %d, n2 = %d, n3 = %d, all must be positive\n",
           n1, n2, n3) ;
   exit(-1) ;
}
nvtx = n1 * n2 * n3 ;
outFileName = argv[7] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n stencil  -- %d" 
        "\n n1       -- %d" 
        "\n n2       -- %d" 
        "\n n3       -- %d" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], stencil, n1, n2, n3, outFileName) ;
fflush(msgFile) ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
graph = Graph_new() ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after setting default fields") ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   initialize the Graph object
   ---------------------------
*/
Graph_init1(graph, 0, nvtx, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
adjIVL = graph->adjIVL ;
list   = IVinit(nvtx, -1) ;
switch ( stencil ) {
case 5 :
   for ( j = 0 ; j < n2 ; j++ ) {
     for ( i = 0 ; i < n1 ; i++ ) {
         ij = i + j * n1 ;
         count = 0 ;
         if ( j >= 1 ) {
            list[count++] = ij - n1 ;
         }
         if ( i >= 1 ) {
            list[count++] = ij - 1 ;
         }
/*
         list[count++] = ij ;
*/
         if ( i <= n1 - 2 ) {
            list[count++] = ij + 1 ;
         }
         if ( j <= n2 - 2 ) {
            list[count++] = ij + n1 ;
         }
         IVqsortUp(count, list) ;
         IVL_setList(adjIVL, ij, count, list) ;
      }
   }
   break ;
case 7 :
   for ( k = 0 ; k < n3 ; k++ ) {
      for ( j = 0 ; j < n2 ; j++ ) {
        for ( i = 0 ; i < n1 ; i++ ) {
            ijk = i + j*n1 + k*n1*n2 ;
            count = 0 ;
            if ( k >= 1 ) {
               list[count++] = ijk - n1*n2 ;
            }
            if ( j >= 1 ) {
               list[count++] = ijk - n1 ;
            }
            if ( i >= 1 ) {
               list[count++] = ijk - 1 ;
            }
            list[count++] = ijk ;
            if ( i <= n1 - 2 ) {
               list[count++] = ijk + 1 ;
            }
            if ( j <= n2 - 2 ) {
               list[count++] = ijk + n1 ;
            }
            if ( k <= n3 - 2 ) {
               list[count++] = ijk + n1*n2 ;
            }
            IVqsortUp(count, list) ;
            IVL_setList(adjIVL, ijk, count, list) ;
         }
      }
   }
   break ;
case 9 :
   for ( j = 0 ; j < n2 ; j++ ) {
     for ( i = 0 ; i < n1 ; i++ ) {
         ij = i + j * n1 ;
         count = 0 ;
         if ( j >= 1 ) {
            if ( i >= 1 ) {
               list[count++] = ij - n1 - 1 ;
            }
            list[count++] = ij - n1 ;
            if ( i <= n1 - 2 ) {
               list[count++] = ij - n1 + 1 ;
            }
         }
         if ( i >= 1 ) {
            list[count++] = ij - 1 ;
         }
         list[count++] = ij ;
         if ( i <= n1 - 2 ) {
            list[count++] = ij + 1 ;
         }
         if ( j <= n2 - 2 ) {
            if ( i >= 1 ) {
               list[count++] = ij + n1 - 1 ;
            }
            list[count++] = ij + n1 ;
            if ( i <= n1 - 2 ) {
               list[count++] = ij + n1 + 1 ;
            }
         }
         IVqsortUp(count, list) ;
         IVL_setList(adjIVL, ij, count, list) ;
      }
   }
   break ;
case 13 :
   for ( j = 0 ; j < n2 ; j++ ) {
     for ( i = 0 ; i < n1 ; i++ ) {
         ij = i + j * n1 ;
         count = 0 ;
         if ( j >= 2 ) {
            list[count++] = ij - 2*n1 ;
         }
         if ( j >= 1 ) {
            if ( i >= 1 ) {
               list[count++] = ij - n1 - 1 ;
            }
            list[count++] = ij - n1 ;
            if ( i <= n1 - 2 ) {
               list[count++] = ij - n1 + 1 ;
            }
         }
         if ( i >= 2 ) {
            list[count++] = ij - 2 ;
         }
         if ( i >= 1 ) {
            list[count++] = ij - 1 ;
         }
         list[count++] = ij ;
         if ( i <= n1 - 2 ) {
            list[count++] = ij + 1 ;
         }
         if ( i <= n1 - 3 ) {
            list[count++] = ij + 2 ;
         }
         if ( j <= n2 - 2 ) {
            if ( i >= 1 ) {
               list[count++] = ij + n1 - 1 ;
            }
            list[count++] = ij + n1 ;
            if ( i <= n1 - 2 ) {
               list[count++] = ij + n1 + 1 ;
            }
         }
         if ( j <= n2 - 3 ) {
            list[count++] = ij + 2*n1 ;
         }
         IVqsortUp(count, list) ;
         IVL_setList(adjIVL, ij, count, list) ;
      }
   }
   break ;
case 27 :
   for ( k = 0 ; k < n3 ; k++ ) {
      kmin = (k > 0)    ? k-1 : k ;
      kmax = (k < n3-1) ? k+1 : k ;
      for ( j = 0 ; j < n2 ; j++ ) {
         jmin = (j > 0)    ? j-1 : j ;
         jmax = (j < n2-1) ? j+1 : j ;
         for ( i = 0 ; i < n1 ; i++ ) {
            ijk = i + j*n1 + k*n1*n2 ;
            imin = (i > 0)    ? i-1 : i ;
            imax = (i < n1-1) ? i+1 : i ;
            for ( know = kmin, count = 0 ; know <= kmax ; know++ ) {
               for ( jnow = jmin ; jnow <= jmax ; jnow++ ) {
                  for ( inow = imin ; inow <= imax ; inow++ ) {
                     list[count++] = inow + jnow*n1 + know*n1*n2 ;
                  }
               }
            }
            IVqsortUp(count, list) ;
            IVL_setList(adjIVL, ijk, count, list) ;
         }
      }
   }
   break ;
default :
   break ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n ") ;
   Graph_writeForHumanEye(graph, msgFile) ;
}
/*
   --------------------------
   write out the Graph object
   --------------------------
*/
if ( strcmp(outFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = Graph_writeToFile(graph, outFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write graph to file %s",
           t2 - t1, outFileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from Graph_writeToFile(%p,%s)",
              rc, graph, outFileName) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;
IVfree(list) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
