/*  test_sortAndCompress.c  */

#include "../Utilities.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------------------------
   purpose -- to test the sort routines

   created -- 98jan28, cca

   usage : test_sortAndCompress msglvl msgFile target n range mod seed

   msglvl -- message level
      msglvl == 1 --> just timings
      msglvl >  1 --> vector(s) before and after sort
   msgFile -- message file
   target -- domain of sort
      IV    -- int vector sorts
      IV2   -- (int,int) vector sorts
      IVDV  -- (int,double) vector sorts
      IV2DV -- (int,int,double) vector sorts
      IVZV  -- (int,complex) vector sorts
      IV2ZV -- (int,int,complex) vector sorts
   n     -- length of vector(s)
   range -- entries are in [1, range]
   mod   -- entries in vectors are uniform(0,n-1) % mod
   seed  -- seed for random number generator
   ------------------------------------------------------------------
*/
{
char     *target ;
double   t1, t2 ;
Drand    drand ;
FILE     *msgFile ;
int      i, ierr, ii, key1, key2, jj,
         msglvl, mod, n, nnew, range, rc, seed ;

if ( argc != 8 ) {
   fprintf(stdout, 
         "\n usage : %s msglvl msgFile target n range mod seed"
         "\n   msglvl -- message level"
         "\n      msglvl == 1 --> just timings"
         "\n      msglvl >  1 --> vector(s) before and after sort"
         "\n   msgFile -- message file"
         "\n   target -- domain of sort"
         "\n      IV    -- int vector sorts"
         "\n      IV2   -- (int,int) vector sorts"
         "\n      IVDV  -- (int,double) vector sorts"
         "\n      IV2DV -- (int,int,double) vector sorts"
         "\n      IVZV  -- (int,complex) vector sorts"
         "\n      IV2ZV -- (int,int,complex) vector sorts"
         "\n      DV    -- double vector sorts"
         "\n      DV2   -- (double,double) vector sorts"
         "\n      DVIV  -- (double,int) vector sorts"
         "\n   n     -- length of vector(s)"
         "\n   range -- entries in vectors are uniform(1,range) %% mod"
         "\n   mod   -- entries in vectors are uniform(1,range) %% mod"
         "\n   seed  -- seed for random number generator"
         "\n", argv[0]) ;
   exit(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp("stdout", argv[2]) == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n unable to open file %s\n", argv[2]) ;
   exit(-1) ;
}
target   = argv[3] ;
n        = atoi(argv[4]) ;
range    = atoi(argv[5]) ;
mod      = atoi(argv[6]) ;
seed     = atoi(argv[7]) ;
fprintf(msgFile, "\n\n ### target = %s, n = %d", target, n) ;
/*
   -------------------------
   switch over the sort type
   -------------------------
*/
Drand_setDefaultFields(&drand) ;
if ( strcmp(target, "IV") == 0 ) {
/*
   -------------------
   IVsortUpAndCompress
   -------------------
*/
   int   *ivec1 ;

   ivec1 = IVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 1, range+1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   for ( i = 0 ; i < n ; i++ ) {
      ivec1[i] = ivec1[i] % mod ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n IVsortAndCompress") ;
      fprintf(msgFile, "\n initial vector") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   nnew = IVsortUpAndCompress(n, ivec1) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %9.5f : time to sort and compress %d entries",
           t2 - t1, n) ;
   fprintf(msgFile, "\n %d entries in compressed vector", nnew) ;
   rc = IVisascending(n, ivec1) ;
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted vector") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "IV2") == 0 ) {
/*
   --------------------
   IV2sortUpAndCompress
   --------------------
*/
   int   *ivec1, *ivec2, *ivec3, *ivec4 ;

   Drand_setDefaultFields(&drand) ;
   ivec1 = IVinit2(n) ;
   ivec2 = IVinit2(n) ;
   ivec3 = IVinit2(n) ;
   ivec4 = IVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, n-1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_fillIvector(&drand, n, ivec2) ;
   for ( ii = 0 ; ii < n ; ii++ ) {
      ivec1[ii] = ivec1[ii] % mod ;
      ivec2[ii] = ivec2[ii] % mod ;
      ivec3[ii] = ivec1[ii] ;
      ivec4[ii] = ivec2[ii] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d initial pairs", n) ;
      for ( ii = 0 ; ii < n ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12d >", ivec1[ii], ivec2[ii]) ;
      }
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   nnew = IV2sortUpAndCompress(n, ivec1, ivec2) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : IV2sortUpAndCompress", t2 - t1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d sorted and compressed pairs", nnew) ;
      for ( ii = 0 ; ii < nnew ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12d >", ivec1[ii], ivec2[ii]) ;
      }
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < n ; ii++ ) {
      key1 = ivec3[ii] ;
      key2 = ivec4[ii] ;
      for ( jj = 0 ; jj < nnew ; jj++ ) {
         if ( key1 == ivec1[jj] && key2 == ivec2[jj] ) {
            if ( msglvl > 1 ) {
               fprintf(msgFile, "\n <%d,%d> found in entry %d",
                       key1, key2, jj) ;
            }
            break ;
         }
      }
      if ( jj == nnew ) {
         fprintf(msgFile, 
                 "\n error, <%d,%d> not found in compressed list",
                 key1, key2) ;
      }
   }
} else if ( strcmp(target, "IVDV") == 0 ) {
/*
   ---------------------
   IVDVsortUpAndCompress
   ---------------------
*/
   int    *ivec1 = IVinit2(n) ;
   int    *ivec2 = IVinit2(n) ;
   double *dvec1 = DVinit2(n) ;
   double *dvec2 = DVinit2(n) ;

   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, n-1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   DVfill(n, dvec1, 1.0) ;
   Drand_fillDvector(&drand, n, dvec1) ;
   IVcopy(n, ivec2, ivec1) ;
   DVcopy(n, dvec2, dvec1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d initial pairs", n) ;
      for ( ii = 0 ; ii < n ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12.4e >", ivec1[ii], dvec1[ii]) ;
      }
   fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   nnew = IVDVsortUpAndCompress(n, ivec1, dvec1) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %9.5f : IVDVsortUpAndCompress time", t2 - t1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d final pairs", nnew) ;
      for ( ii = 0 ; ii < nnew ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12.4e >", ivec1[ii], dvec1[ii]) ;
      }
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < n ; ii++ ) {
      key1 = ivec2[ii] ;
      for ( jj = 0 ; jj < nnew ; jj++ ) {
         if ( key1 == ivec1[jj] ) {
            dvec1[jj] -= dvec2[ii] ;
            break ;
         }
      }
      if ( jj == nnew ) {
         fprintf(msgFile, "\n pair < %d,%f> not found", 
                 ivec2[ii], dvec2[ii]) ;
      }
   }
   fprintf(msgFile, 
           "\n error norm = %12.4e", DVmaxabs(nnew, dvec1, &ii)) ;
   fprintf(msgFile, "\n dvec1[]") ;
   DVfprintf(msgFile, nnew, dvec1) ;
} else if ( strcmp(target, "IV2DV") == 0 ) {
/*
   ----------------------
   IV2DVsortUpAndCompress
   ----------------------
*/
   int      *ivec1 = IVinit2(n) ;
   int      *ivec2 = IVinit2(n) ;
   int      *ivec3 = IVinit2(n) ;
   int      *ivec4 = IVinit2(n) ;
   double   *dvec1 = DVinit2(n) ;
   double   *dvec2 = DVinit2(n) ;

   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, n-1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_fillIvector(&drand, n, ivec2) ;
   Drand_setUniform(&drand, 0, 1) ;
   Drand_fillDvector(&drand, n, dvec1) ;
   for ( ii = 0 ; ii < n ; ii++ ) {
      ivec1[ii] = ivec1[ii] % mod ;
      ivec2[ii] = ivec2[ii] % mod ;
   }
   IVcopy(n, ivec3, ivec1) ;
   IVcopy(n, ivec4, ivec2) ;
   DVcopy(n, dvec2, dvec1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d initial triples", n) ;
      for ( ii = 0 ; ii < n ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12d, %12.4e >", 
                 ivec1[ii], ivec2[ii], dvec1[ii]) ;
      }
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   nnew = IV2DVsortUpAndCompress(n, ivec1, ivec2, dvec1) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %9.5f : IV2DVsortUpAndCompress, n = %d, nnew = %d", 
           t2 - t1, n, nnew) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d sorted and compressed triples", nnew) ;
      for ( ii = 0 ; ii < nnew ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12d, %12.4e >", 
                 ivec1[ii], ivec2[ii], dvec1[ii]) ;
      }
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < n ; ii++ ) {
      key1 = ivec3[ii] ;
      key2 = ivec4[ii] ;
      for ( jj = 0 ; jj < nnew ; jj++ ) {
         if ( key1 == ivec1[jj] && key2 == ivec2[jj] ) {
            if ( msglvl > 1 ) {
               fprintf(msgFile, "\n <%d,%d> found in entry %d",
                       key1, key2, jj) ;
            }
            dvec1[jj] -= dvec2[ii] ;
            break ;
         }
      }
      if ( jj == nnew ) {
         fprintf(msgFile, 
                 "\n error, <%d,%d> not found in compressed list",
                 key1, key2) ;
      }
   }
   fprintf(msgFile, "\n max error = %12.4e\n",
           DVmaxabs(nnew, dvec1, &ii)) ;
   if ( msglvl > 1 ) {
      DVfprintf(msgFile, nnew, dvec1) ;
   }
} else if ( strcmp(target, "IV2ZV") == 0 ) {
/*
   ----------------------
   IV2ZVsortUpAndCompress
   ----------------------
*/
   int      *ivec1 = IVinit2(n) ;
   int      *ivec2 = IVinit2(n) ;
   int      *ivec3 = IVinit2(n) ;
   int      *ivec4 = IVinit2(n) ;
   double   *dvec1 = DVinit2(2*n) ;
   double   *dvec2 = DVinit2(2*n) ;

   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, n-1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_fillIvector(&drand, n, ivec2) ;
   Drand_setUniform(&drand, 0, 1) ;
   Drand_fillDvector(&drand, 2*n, dvec1) ;
   for ( ii = 0 ; ii < n ; ii++ ) {
      ivec1[ii] = ivec1[ii] % mod ;
      ivec2[ii] = ivec2[ii] % mod ;
   }
   IVcopy(n,   ivec3, ivec1) ;
   IVcopy(n,   ivec4, ivec2) ;
   DVcopy(2*n, dvec2, dvec1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d initial triples", n) ;
      for ( ii = 0 ; ii < n ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12d, %12.4e, %12.4e >", 
                 ivec1[ii], ivec2[ii], dvec1[2*ii], dvec1[2*ii+1]) ;
      }
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   nnew = IV2ZVsortUpAndCompress(n, ivec1, ivec2, dvec1) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %9.5f : IV2ZVsortUpAndCompress, n = %d, nnew = %d", 
           t2 - t1, n, nnew) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d sorted and compressed triples", nnew) ;
      for ( ii = 0 ; ii < nnew ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12d, %12.4e, %12.4e >", 
                 ivec1[ii], ivec2[ii], dvec1[2*ii], dvec1[2*ii+1]) ;
      }
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < n ; ii++ ) {
      key1 = ivec3[ii] ;
      key2 = ivec4[ii] ;
      for ( jj = 0 ; jj < nnew ; jj++ ) {
         if ( key1 == ivec1[jj] && key2 == ivec2[jj] ) {
            if ( msglvl > 1 ) {
               fprintf(msgFile, "\n <%d,%d> found in entry %d",
                       key1, key2, jj) ;
            }
            dvec1[2*jj]   -= dvec2[2*ii]   ;
            dvec1[2*jj+1] -= dvec2[2*ii+1] ;
            break ;
         }
      }
      if ( jj == nnew ) {
         fprintf(msgFile, 
                 "\n error, <%d,%d> not found in compressed list",
                 key1, key2) ;
      }
   }
   fprintf(msgFile, "\n max error = %12.4e\n",
           DVmaxabs(2*nnew, dvec1, &ii)) ;
   if ( msglvl > 1 ) {
      DVfprintf(msgFile, 2*nnew, dvec1) ;
   }
} else if ( strcmp(target, "IVZV") == 0 ) {
/*
   ---------------------
   IVZVsortUpAndCompress
   ---------------------
*/
   int    *ivec1 = IVinit2(n) ;
   int    *ivec2 = IVinit2(n) ;
   double *dvec1 = DVinit2(2*n) ;
   double *dvec2 = DVinit2(2*n) ;

   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, n-1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   DVfill(n, dvec1, 1.0) ;
   Drand_fillDvector(&drand, 2*n, dvec1) ;
   IVcopy(n,   ivec2, ivec1) ;
   DVcopy(2*n, dvec2, dvec1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d initial pairs", n) ;
      for ( ii = 0 ; ii < n ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12.4e, %12.4e >", 
                 ivec1[ii], dvec1[2*ii], dvec1[2*ii+1]) ;
      }
   fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   nnew = IVZVsortUpAndCompress(n, ivec1, dvec1) ;
   MARKTIME(t2) ;
   fprintf(msgFile, 
           "\n CPU %9.5f : IVZVsortUpAndCompress time", t2 - t1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d final pairs", nnew) ;
      for ( ii = 0 ; ii < nnew ; ii++ ) {
         fprintf(msgFile, "\n < %12d, %12.4e, %12.4e >", 
                 ivec1[ii], dvec1[2*ii], dvec1[2*ii+1]) ;
      }
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < n ; ii++ ) {
      key1 = ivec2[ii] ;
      for ( jj = 0 ; jj < nnew ; jj++ ) {
         if ( key1 == ivec1[jj] ) {
            dvec1[2*jj]   -= dvec2[2*ii] ;
            dvec1[2*jj+1] -= dvec2[2*ii+1] ;
            break ;
         }
      }
      if ( jj == nnew ) {
         fprintf(msgFile, "\n pair < %d,%f,%f> not found", 
                 ivec2[ii], dvec2[2*ii], dvec2[2*ii+1]) ;
      }
   }
   fprintf(msgFile, 
           "\n error norm = %12.4e", DVmaxabs(2*nnew, dvec1, &ii)) ;
   fprintf(msgFile, "\n dvec1[]") ;
   DVfprintf(msgFile, 2*nnew, dvec1) ;
}

exit(0) ; }

/*--------------------------------------------------------------------*/
