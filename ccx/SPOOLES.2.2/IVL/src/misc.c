/*  misc.c  */

#include "../IVL.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to make and return a 9-point adjacency structure

   input --

      n1    -- # of grid points in first direction
      n2    -- # of grid points in second direction
      ncomp -- # of components per grid point
   -----------------------------------------------------------
*/
IVL *
IVL_make9P ( 
   int   n1, 
   int   n2, 
   int   ncomp 
) {
IVL   *ivl ;
int   i, icomp, idof, ij, imax, imin, inow, j, jcomp, 
      jmax, jmin, jnow, k, naind, ndof, nvtx, size ;
int   *indices ;

if ( n1 <= 0 || n2 <= 0 || ncomp <= 0 ) {
   return(NULL) ;
}
/*
   ----------------------------------------
   initialize the adjacency graph structure
   ----------------------------------------
*/
nvtx  = n1 * n2 ;
ndof  = nvtx * ncomp ;
naind = ncomp*ncomp*((n1-2)*(n2-2)*9 + 2*(n1+n2-4)*6 + 4*4) ;
ivl   = IVL_new() ;
IVL_init2(ivl, IVL_CHUNKED, ndof, naind) ;
indices = IVinit(9*ncomp, -1) ;
/*
   ----------------------------
   fill the adjacency structure
   ----------------------------
*/
idof = 0 ;
k = 0 ;
for ( j = 0 ; j < n2 ; j++ ) {
  for ( i = 0 ; i < n1 ; i++ ) {
      ij = i + j * n1 ;
      imin = (i > 0)    ? i-1 : i ;
      imax = (i < n1-1) ? i+1 : i ;
      jmin = (j > 0)    ? j-1 : j ;
      jmax = (j < n2-1) ? j+1 : j ;
      for ( icomp = 0 ; icomp < ncomp ; icomp++, idof++ ) {
         size  = ncomp*(imax - imin + 1)*(jmax - jmin + 1) ;
         for ( jnow = jmin, k = 0 ; jnow <= jmax ; jnow++ ) {
            for ( inow = imin ; inow <= imax ; inow++ ) {
               for ( jcomp = 0 ; jcomp < ncomp ; jcomp++ ) {
                  indices[k++] = (inow + jnow*n1)*ncomp + jcomp ;
               }
            }
         }
         IVL_setList(ivl, idof, size, indices) ;
      }
   }
}
IVfree(indices) ;

return(ivl) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to make and return a 27-point adjacency structure

   input --

      n1    -- # of grid points in first direction
      n2    -- # of grid points in second direction
      n3    -- # of grid points in second direction
      ncomp -- # of components per grid point
   ------------------------------------------------------------
*/
IVL *
IVL_make27P ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   ncomp 
) {
IVL   *ivl ;
int   i, icomp, idof, imax, imin, inow, j, jcomp, jmax, jmin, jnow, 
      k, kmax, kmin, know, m, naind, ndof, nvtx, size ;
int   *indices ;

if ( n1 <= 0 || n2 <= 0 || n3 <= 0 || ncomp <= 0 ) {
   return(NULL) ;
}
/*
   ----------------------------------------
   initialize the adjacency graph structure
   ----------------------------------------
*/
nvtx    = n1 * n2 * n3 ;
ndof    = nvtx * ncomp ;
naind   = ncomp*ncomp*ncomp*(
            27*(n1-2)*(n2-2)*(n3-2)
          + 2*18*(n1-2)*(n2-2)
          + 2*18*(n1-2)*(n3-2)
          + 2*18*(n2-2)*(n3-2)
          + 4*12*(n1-2)
          + 4*12*(n2-2)
          + 4*12*(n3-2)
          + 8*8) ;
ivl = IVL_new() ;
IVL_init2(ivl, IVL_CHUNKED, ndof, naind) ;
indices = IVinit(27*ncomp, -1) ;
/*
   ----------------------------
   fill the adjacency structure
   ----------------------------
*/
idof = 0 ;
m = 0 ;
for ( k = 0 ; k < n3 ; k++ ) {
   kmin = (k > 0)    ? k-1 : k ;
   kmax = (k < n3-1) ? k+1 : k ;
   for ( j = 0 ; j < n2 ; j++ ) {
      jmin = (j > 0)    ? j-1 : j ;
      jmax = (j < n2-1) ? j+1 : j ;
      for ( i = 0 ; i < n1 ; i++ ) {
         imin = (i > 0)    ? i-1 : i ;
         imax = (i < n1-1) ? i+1 : i ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++, idof++ ) {
            size  = ncomp*(imax - imin + 1)
                                 *(jmax - jmin + 1)*(kmax - kmin + 1) ;
/*
fprintf(stdout, "\n k = %d, j = %d, i = %d, icomp = %d, size = %d", 
        k, j, i, icomp, size) ;
fflush(stdout) ;
fprintf(stdout, "\n dof %d : m %d :", idof, m) ;
*/
            for ( know = kmin, m = 0 ; know <= kmax ; know++ ) {
               for ( jnow = jmin ; jnow <= jmax ; jnow++ ) {
                  for ( inow = imin ; inow <= imax ; inow++ ) {
                     for ( jcomp = 0 ; jcomp < ncomp ; jcomp++ ) {
                        if ( m == naind ) {
                            fprintf(stderr, 
                            "\n error in IVL::IVLmake27P(%d,%d,%d,%d)"
                            "\n naind = %d, m = %d"
                            "\n (i,j,k) = (%d,%d,%d),"
                            " (inow,jnow,know) = (%d,%d,%d)",
                            n1, n2, n3, ncomp, naind, m, i, j, k, 
                            inow, jnow, know) ;
                            exit(-1) ;
                        }
                        indices[m++] = jcomp 
                                 + (inow + jnow*n1 + know*n1*n2)*ncomp ;
/*
                        fprintf(stdout, " %d", indices[m - 1]) ;
*/
                     }
                  }
               }
            }
            IVL_setList(ivl, idof, size, indices) ;
         }
      }
   }
}
IVfree(indices) ;

return(ivl) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to make and return a 13-point adjacency structure

   input --

      n1 -- # of grid points in first direction
      n2 -- # of grid points in second direction

   created -- 96feb01
   ------------------------------------------------------------
*/
IVL *
IVL_make13P ( 
   int   n1, 
   int   n2 
) {
IVL   *ivl ;
int   count, i, ij, j, nvtx ;
int   list[13] ;

if ( n1 <= 0 || n2 <= 0 ) {
   return(NULL) ;
}
/*
   ----------------------------------------
   initialize the adjacency graph structure
   ----------------------------------------
*/
nvtx  = n1 * n2 ;
ivl   = IVL_new() ;
IVL_init1(ivl, IVL_CHUNKED, nvtx) ;
/*
   ----------------------------
   fill the adjacency structure
   ----------------------------
*/
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
      IVL_setList(ivl, ij, count, list) ;
   }
}

return(ivl) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to make and return a 5-point adjacency structure

   input --

      n1 -- # of grid points in first direction
      n2 -- # of grid points in second direction

   created -- 96feb02
   -----------------------------------------------------------
*/
IVL *
IVL_make5P ( 
   int   n1, 
   int   n2 
) {
IVL   *ivl ;
int   count, i, ij, j, nvtx ;
int   list[5] ;

if ( n1 <= 0 || n2 <= 0 ) {
   return(NULL) ;
}
/*
   ----------------------------------------
   initialize the adjacency graph structure
   ----------------------------------------
*/
nvtx  = n1 * n2 ;
ivl   = IVL_new() ;
IVL_init1(ivl, IVL_CHUNKED, nvtx) ;
/*
   ----------------------------
   fill the adjacency structure
   ----------------------------
*/
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
      list[count++] = ij ;
      if ( i <= n1 - 2 ) {
         list[count++] = ij + 1 ;
      }
      if ( j <= n2 - 2 ) {
         list[count++] = ij + n1 ;
      }
/*
      if ( j >= 1 ) {
         list[count++] = ij - n1 ;
      } else {
         list[count++] = i + n1 * (n2 - 1) ;
      }
      if ( i >= 1 ) {
         list[count++] = ij - 1 ;
      } else {
         list[count++] = ij + n1 - 1 ;
      } 
      list[count++] = ij ;
      if ( i <= n1 - 2 ) {
         list[count++] = ij + 1 ;
      } else {
         list[count++] = ij - n1 + 1 ;
      }
      if ( j <= n2 - 2 ) {
         list[count++] = ij + n1 ;
      } else {
         list[count++] = i ;
      }
*/
      IVqsortUp(count, list) ;
      IVL_setList(ivl, ij, count, list) ;
   }
}

return(ivl) ; }

/*--------------------------------------------------------------------*/
