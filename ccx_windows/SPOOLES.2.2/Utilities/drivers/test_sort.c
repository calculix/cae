/*  test_sort.c  */

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

   usage : test_sort msglvl msgFile target sortType n range mod seed

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
      DV    -- double vector sorts
      DV2   -- (double,double) vector sorts
      DVIV  -- (double,int) vector sorts
   sortType -- type of sort
      IU -- ascending insert sort
      ID -- descending insert sort
      QU -- asscending quick sort
      QD -- descending quick sort
   n     -- length of vector(s)
   range -- entries are in [1, range]
   mod   -- entries in vectors are uniform(0,n-1) % mod
   seed  -- seed for random number generator
   ------------------------------------------------------------------
*/
{
char     *sortType, *target ;
double   t1, t2 ;
Drand    drand ;
FILE     *msgFile ;
int      i, ierr, msglvl, mod, n, range, rc, seed ;

if ( argc != 9 ) {
   fprintf(stdout, 
         "\n usage : %s msglvl msgFile target sortType n range mod seed"
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
         "\n   sortType -- type of sort"
         "\n      IU -- ascending insert sort"
         "\n      ID -- descending insert sort"
         "\n      QU -- asscending quick sort"
         "\n      QD -- descending quick sort"
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
sortType = argv[4] ;
n        = atoi(argv[5]) ;
range    = atoi(argv[6]) ;
mod      = atoi(argv[7]) ;
seed     = atoi(argv[8]) ;
fprintf(msgFile, "\n\n ### target = %s, sortType = %s, n = %d",
        target, sortType, n) ;
/*
   -------------------------
   switch over the sort type
   -------------------------
*/
Drand_setDefaultFields(&drand) ;
if ( strcmp(target, "IV") == 0 ) {
/*
   ------
   IVsort
   ------
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
      fprintf(msgFile, "\n IVsort") ;
      fprintf(msgFile, "\n initial vector") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      IVisortUp(n, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IVisortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      IVisortDown(n, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IVisortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      IVqsortUp(n, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IVqsortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      IVqsortDown(n, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IVqsortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   }
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
   -------
   IV2sort
   -------
*/
   double   sig1, sig2 ;
   int      *ivec1, *ivec2 ;

   ivec1 = IVinit2(n) ;
   ivec2 = IVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 1, range+1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_fillIvector(&drand, n, ivec2) ;
   for ( i = 0 ; i < n ; i++ ) {
      ivec1[i] = ivec1[i] % mod ;
      ivec2[i] = ivec2[i] % mod ;
   }
   for ( i = 0, sig1 = 0.0 ; i < n ; i++ ) {
      sig1 += ivec1[i]*ivec2[i] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n IV2sort, sig1 = %12.4e", sig1) ;
      fprintf(msgFile, "\n initial ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n initial ivec2") ;
      IVfp80(msgFile, n, ivec2, 80, &ierr) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      IV2isortUp(n, ivec1, ivec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IV2isortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      IV2isortDown(n, ivec1, ivec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IV2isortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      IV2qsortUp(n, ivec1, ivec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IV2qsortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      IV2qsortDown(n, ivec1, ivec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via IV2qsortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   }
   for ( i = 0, sig2 = 0.0 ; i < n ; i++ ) {
      sig2 += ivec1[i]*ivec2[i] ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   fprintf(msgFile, 
           "\n sig1 = %12.4e, sig2 = %12.4e, error/sig1 = %12.4e\n",
           sig1, sig2, (sig1 - sig2)/sig1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n companion ivec2") ;
      IVfp80(msgFile, n, ivec2, 80, &ierr) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "IVDV") == 0 ) {
/*
   --------
   IVDVsort
   --------
*/
   double   sig1, sig2 ;
   double   *dvec ;
   int      *ivec1 ;

   ivec1 = IVinit2(n) ;
   dvec  = DVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 1, range+1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_setUniform(&drand, 0.0, 1.0) ;
   Drand_fillDvector(&drand, n, dvec) ;
   for ( i = 0 ; i < n ; i++ ) {
      ivec1[i] = ivec1[i] % mod ;
   }
   for ( i = 0, sig1 = 0.0 ; i < n ; i++ ) {
      sig1 += ivec1[i]*dvec[i] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n IVDVsort, sig1 = %12.4e", sig1) ;
      fprintf(msgFile, "\n initial ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n initial dvec") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      IVDVisortUp(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVDVisortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      IVDVisortDown(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVDVisortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      IVDVqsortUp(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVDVqsortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      IVDVqsortDown(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVDVqsortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   }
   for ( i = 0, sig2 = 0.0 ; i < n ; i++ ) {
      sig2 += ivec1[i]*dvec[i] ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   fprintf(msgFile, 
           "\n sig1 = %12.4e, sig2 = %12.4e, error/sig1 = %12.4e\n",
           sig1, sig2, (sig1 - sig2)/sig1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n companion dvec") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "IV2DV") == 0 ) {
/*
   ---------
   IV2DVsort
   ---------
*/
   double   sig1, sig2 ;
   double   *dvec ;
   int      *ivec1, *ivec2 ;

   ivec1 = IVinit2(n) ;
   ivec2 = IVinit2(n) ;
   dvec  = DVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 1, range+1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_fillIvector(&drand, n, ivec2) ;
   Drand_setUniform(&drand, 0.0, 1.0) ;
   Drand_fillDvector(&drand, n, dvec) ;
   for ( i = 0 ; i < n ; i++ ) {
      ivec1[i] = ivec1[i] % mod ;
      ivec2[i] = ivec2[i] % mod ;
   }
   for ( i = 0, sig1 = 0.0 ; i < n ; i++ ) {
      sig1 += ivec1[i]*ivec2[i]*dvec[i] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n IV2DVsort, sig1 = %12.4e", sig1) ;
      fprintf(msgFile, "\n initial ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n initial ivec2") ;
      IVfp80(msgFile, n, ivec2, 80, &ierr) ;
      fprintf(msgFile, "\n initial dvec") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      IV2DVisortUp(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2DVisortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      IV2DVisortDown(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2DVisortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      IV2DVqsortUp(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2DVqsortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      IV2DVqsortDown(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2DVqsortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   }
   for ( i = 0, sig2 = 0.0 ; i < n ; i++ ) {
      sig2 += ivec1[i]*ivec2[i]*dvec[i] ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   fprintf(msgFile, 
           "\n sig1 = %12.4e, sig2 = %12.4e, error/sig1 = %12.4e\n",
           sig1, sig2, (sig1 - sig2)/sig1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n companion ivec2") ;
      IVfp80(msgFile, n, ivec2, 80, &ierr) ;
      fprintf(msgFile, "\n companion dvec") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "IVZV") == 0 ) {
/*
   --------
   IVZVsort
   --------
*/
   double   sig1, sig2 ;
   double   *dvec ;
   int      *ivec1 ;

   ivec1 = IVinit2(n) ;
   dvec  = DVinit2(2*n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 1, range+1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_setUniform(&drand, 0.0, 1.0) ;
   Drand_fillDvector(&drand, 2*n, dvec) ;
   for ( i = 0 ; i < n ; i++ ) {
      ivec1[i] = ivec1[i] % mod ;
   }
   for ( i = 0, sig1 = 0.0 ; i < n ; i++ ) {
      sig1 += ivec1[i]*dvec[2*i]*dvec[2*i+1] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n IVZVsort, sig1 = %12.4e", sig1) ;
      fprintf(msgFile, "\n initial ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n initial dvec") ;
      DVfprintf(msgFile, 2*n, dvec) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      IVZVisortUp(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVZVisortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      IVZVisortDown(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVZVisortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      IVZVqsortUp(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVZVqsortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      IVZVqsortDown(n, ivec1, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IVZVqsortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   }
   for ( i = 0, sig2 = 0.0 ; i < n ; i++ ) {
      sig2 += ivec1[i]*dvec[2*i]*dvec[2*i+1] ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   fprintf(msgFile, 
           "\n sig1 = %12.4e, sig2 = %12.4e, error/sig1 = %12.4e\n",
           sig1, sig2, (sig1 - sig2)/sig1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n companion dvec") ;
      DVfprintf(msgFile, 2*n, dvec) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "IV2ZV") == 0 ) {
/*
   ---------
   IV2ZVsort
   ---------
*/
   double   sig1, sig2 ;
   double   *dvec ;
   int      *ivec1, *ivec2 ;

   ivec1 = IVinit2(n) ;
   ivec2 = IVinit2(n) ;
   dvec  = DVinit2(2*n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 1, range+1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_fillIvector(&drand, n, ivec2) ;
   Drand_setUniform(&drand, 0.0, 1.0) ;
   Drand_fillDvector(&drand, 2*n, dvec) ;
   for ( i = 0 ; i < n ; i++ ) {
      ivec1[i] = ivec1[i] % mod ;
      ivec2[i] = ivec2[i] % mod ;
   }
   for ( i = 0, sig1 = 0.0 ; i < n ; i++ ) {
      sig1 += ivec1[i]*ivec2[i]*dvec[2*i]*dvec[2*i+1] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n IV2ZVsort, sig1 = %12.4e", sig1) ;
      fprintf(msgFile, "\n initial ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n initial ivec2") ;
      IVfp80(msgFile, n, ivec2, 80, &ierr) ;
      fprintf(msgFile, "\n initial dvec") ;
      DVfprintf(msgFile, 2*n, dvec) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      IV2ZVisortUp(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2ZVisortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      IV2ZVisortDown(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2ZVisortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      IV2ZVqsortUp(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2ZVqsortUp\n",
              t2 - t1, n) ;
      rc = IVisascending(n, ivec1) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      IV2ZVqsortDown(n, ivec1, ivec2, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via IV2ZVqsortDown\n",
              t2 - t1, n) ;
      rc = IVisdescending(n, ivec1) ;
   }
   for ( i = 0, sig2 = 0.0 ; i < n ; i++ ) {
      sig2 += ivec1[i]*ivec2[i]*dvec[2*i]*dvec[2*i+1] ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   fprintf(msgFile, 
           "\n sig1 = %12.4e, sig2 = %12.4e, error/sig1 = %12.4e\n",
           sig1, sig2, (sig1 - sig2)/sig1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n companion ivec2") ;
      IVfp80(msgFile, n, ivec2, 80, &ierr) ;
      fprintf(msgFile, "\n companion dvec") ;
      DVfprintf(msgFile, 2*n, dvec) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "DV") == 0 ) {
/*
   ------
   DVsort
   ------
*/
   double   *dvec ;

   dvec = DVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_fillDvector(&drand, n, dvec) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n DVsort") ;
      fprintf(msgFile, "\n initial vector") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      DVisortUp(n, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via DVisortUp\n",
              t2 - t1, n) ;
      rc = DVisascending(n, dvec) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      DVisortDown(n, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %9.5f : sort %d via DVisortDown\n",
              t2 - t1, n) ;
      rc = DVisdescending(n, dvec) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      DVqsortUp(n, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %9.5f : sort %d via DVqsortUp\n",
              t2 - t1, n) ;
      rc = DVisascending(n, dvec) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      DVqsortDown(n, dvec) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %9.5f : sort %d via DVqsortDown\n",
              t2 - t1, n) ;
      rc = DVisdescending(n, dvec) ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted vector") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "DVIV") == 0 ) {
/*
   -------
   DVIVsort
   -------
*/
   double   sig1, sig2 ;
   double   *dvec ;
   int      *ivec1 ;

   ivec1 = IVinit2(n) ;
   dvec  = DVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 1, range+1) ;
   Drand_fillIvector(&drand, n, ivec1) ;
   Drand_setUniform(&drand, 0, 1) ;
   Drand_fillDvector(&drand, n, dvec) ;
   for ( i = 0 ; i < n ; i++ ) {
      ivec1[i] = ivec1[i] % mod ;
   }
   for ( i = 0, sig1 = 0.0 ; i < n ; i++ ) {
      sig1 += ivec1[i]*dvec[i] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n DVIVsort, sig1 = %12.4e", sig1) ;
      fprintf(msgFile, "\n initial ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n initial dvec") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      DVIVisortUp(n, dvec, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d via DVIVisortUp\n",
              t2 - t1, n) ;
      rc = DVisascending(n, dvec) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      DVIVisortDown(n, dvec, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via DVIVisortDown\n",
              t2 - t1, n) ;
      rc = DVisdescending(n, dvec) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      DVIVqsortUp(n, dvec, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via DVIVqsortUp\n",
              t2 - t1, n) ;
      rc = DVisascending(n, dvec) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      DVIVqsortDown(n, dvec, ivec1) ;
      MARKTIME(t2) ;
      fprintf(msgFile, 
              "\n CPU %9.5f : sort %d integers via DVIVqsortDown\n",
              t2 - t1, n) ;
      rc = DVisdescending(n, dvec) ;
   }
   for ( i = 0, sig2 = 0.0 ; i < n ; i++ ) {
      sig2 += ivec1[i]*dvec[i] ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   fprintf(msgFile, 
           "\n sig1 = %12.4e, sig2 = %12.4e, error/sig1 = %12.4e\n",
           sig1, sig2, (sig1 - sig2)/sig1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted dvec") ;
      DVfprintf(msgFile, n, dvec) ;
      fprintf(msgFile, "\n companion ivec1") ;
      IVfp80(msgFile, n, ivec1, 80, &ierr) ;
      fprintf(msgFile, "\n ") ;
   }
} else if ( strcmp(target, "DV2") == 0 ) {
/*
   -------
   DV2sort
   -------
*/
   double   sig1, sig2 ;
   double   *dvec1, *dvec2 ;

   dvec1 = DVinit2(n) ;
   dvec2 = DVinit2(n) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, 1) ;
   Drand_fillDvector(&drand, n, dvec1) ;
   Drand_fillDvector(&drand, n, dvec2) ;
   for ( i = 0, sig1 = 0.0 ; i < n ; i++ ) {
      sig1 += dvec1[i]*dvec2[i] ;
   }
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n DV2sort, sig1 = %12.4e", sig1) ;
      fprintf(msgFile, "\n initial dvec1") ;
      DVfprintf(msgFile, n, dvec2) ;
      fprintf(msgFile, "\n initial dvec2") ;
      DVfprintf(msgFile, n, dvec2) ;
      fprintf(msgFile, "\n ") ;
      fflush(msgFile) ;
   }
   if ( strcmp(sortType, "IU") == 0 ) {
      MARKTIME(t1) ;
      DV2isortUp(n, dvec1, dvec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %9.5f : sort %d via DV2isortUp\n",
              t2 - t1, n) ;
      rc = DVisascending(n, dvec1) ;
   } else if ( strcmp(sortType, "ID") == 0 ) {
      MARKTIME(t1) ;
      DV2isortDown(n, dvec1, dvec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %9.5f : sort %d via DV2isortDown\n",
              t2 - t1, n) ;
      rc = DVisdescending(n, dvec1) ;
   } if ( strcmp(sortType, "QU") == 0 ) {
      MARKTIME(t1) ;
      DV2qsortUp(n, dvec1, dvec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %9.5f : sort %d via DV2qsortUp\n",
              t2 - t1, n) ;
      rc = DVisascending(n, dvec1) ;
   } else if ( strcmp(sortType, "QD") == 0 ) {
      MARKTIME(t1) ;
      DV2qsortDown(n, dvec1, dvec2) ;
      MARKTIME(t2) ;
      fprintf(msgFile, "\n CPU %9.5f : sort %d via DV2qsortDown\n",
              t2 - t1, n) ;
      rc = DVisdescending(n, dvec1) ;
   }
   for ( i = 0, sig2 = 0.0 ; i < n ; i++ ) {
      sig2 += dvec1[i]*dvec2[i] ;
   }
   if ( rc == 0 ) {
      fprintf(msgFile, "\n vector not sorted correctly") ;
   } else {
      fprintf(msgFile, "\n vector sorted correctly") ;
   }
   fprintf(msgFile, 
           "\n sig1 = %12.4e, sig2 = %12.4e, error/sig1 = %12.4e\n",
           sig1, sig2, (sig1 - sig2)/sig1) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n sorted dvec1") ;
      DVfprintf(msgFile, n, dvec1) ;
      fprintf(msgFile, "\n companion ivec2") ;
      DVfprintf(msgFile, n, dvec2) ;
      fprintf(msgFile, "\n ") ;
   }
}

exit(0) ; }

/*--------------------------------------------------------------------*/
