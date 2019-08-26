/*  mvm.c  */

#include "../spoolesMT.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
#define NONSYM   1
#define SYM      2
#define HERM     3
#define NONSYM_T 4
#define NONSYM_H 5
typedef struct _MTmvmObj   MTmvmObj ;
struct _MTmvmObj {
   InpMtx     *A ;
   DenseMtx   *Y ;
   double     alpha[2] ;
   DenseMtx   *X ;
} ;
static MTmvmObj * setup ( InpMtx *A, DenseMtx *Y, double alpha[],
                          DenseMtx *X, int nthread ) ;
static void InpMtx_MT_mmm ( int flag, InpMtx *A, DenseMtx *Y,
 double alpha[], DenseMtx *X, int nthread, int msglvl, FILE *msgFile ) ;
static void * worker_nonsym_mmm ( void *arg ) ;
static void * worker_sym_mmm ( void *arg ) ;
static void * worker_herm_mmm ( void *arg ) ;
static void * worker_nonsym_mmm_T ( void *arg ) ;
static void * worker_nonsym_mmm_H ( void *arg ) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_nonsym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  A == NULL || Y == NULL || alpha == NULL 
   || X == NULL || nthread < 1 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_MT_nonsym_mmm(%p,%p,%p,%p,%d)"
           "\n bad input\n", A, Y, alpha, X, nthread) ;
   exit(-1) ;
}
if ( nthread == 1 ) {
   InpMtx_nonsym_mmm(A, Y, alpha, X) ;
} else {
   InpMtx_MT_mmm(NONSYM, A, Y, alpha, X, nthread, msglvl, msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_sym_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  A == NULL || Y == NULL || alpha == NULL 
   || X == NULL || nthread < 1 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_MT_sym_mmm(%p,%p,%p,%p,%d)"
           "\n bad input\n", A, Y, alpha, X, nthread) ;
   exit(-1) ;
}
if ( nthread == 1 ) {
   InpMtx_sym_mmm(A, Y, alpha, X) ;
} else {
   InpMtx_MT_mmm(SYM, A, Y, alpha, X, nthread, msglvl, msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_herm_mmm (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  A == NULL || Y == NULL || alpha == NULL 
   || X == NULL || nthread < 1 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_MT_herm_mmm(%p,%p,%p,%p,%d)"
           "\n bad input\n", A, Y, alpha, X, nthread) ;
   exit(-1) ;
}
if ( nthread == 1 ) {
   InpMtx_herm_mmm(A, Y, alpha, X) ;
} else {
   InpMtx_MT_mmm(HERM, A, Y, alpha, X, nthread, msglvl, msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_nonsym_mmm_T (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  A == NULL || Y == NULL || alpha == NULL 
   || X == NULL || nthread < 1 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_MT_nonsym_mmm_T(%p,%p,%p,%p,%d)"
           "\n bad input\n", A, Y, alpha, X, nthread) ;
   exit(-1) ;
}
if ( nthread == 1 ) {
   InpMtx_nonsym_mmm_T(A, Y, alpha, X) ;
} else {
   InpMtx_MT_mmm(NONSYM_T, A, Y, alpha, X, nthread, msglvl, msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to compute Y := Y + alpha*A*X
 
   created -- 98jul09, cca
   ----------------------------------------
*/
void
InpMtx_MT_nonsym_mmm_H (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  A == NULL || Y == NULL || alpha == NULL 
   || X == NULL || nthread < 1 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_MT_nonsym_mmm_H(%p,%p,%p,%p,%d)"
           "\n bad input\n", A, Y, alpha, X, nthread) ;
   exit(-1) ;
}
if ( nthread == 1 ) {
   InpMtx_nonsym_mmm_H(A, Y, alpha, X) ;
} else {
   InpMtx_MT_mmm(NONSYM_H, A, Y, alpha, X, nthread, msglvl, msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
static void
InpMtx_MT_mmm (
   int        flag,
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread,
   int        msglvl,
   FILE       *msgFile
) {
double     t1, t2 ;
int        myid, nent, rc ;
MTmvmObj   *MTmvmObjs, *obj ;
/*
   -------------------------------
   set up the nthread data objects
   -------------------------------
*/
MARKTIME(t1) ;
MTmvmObjs = setup(A, Y, alpha, X, nthread) ;
MARKTIME(t2) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n %% CPU %8.3f : setup time", t2 - t1) ;
}
#if THREAD_TYPE == TT_POSIX 
{
pthread_t        *tids ;
pthread_attr_t   attr  ;
void             *status ;
/*
#####   NOTE: for SGI machines, this command must be present
#####         for the thread scheduling to be efficient.
#####         this is NOT a POSIX call, but SGI needs it anyway
pthread_setconcurrency(nthread) ;
*/
pthread_attr_init(&attr) ;
/*
pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM) ;
*/
pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS) ;
ALLOCATE(tids, pthread_t, nthread) ;
MARKTIME(t1) ;
for ( myid = 0, obj = MTmvmObjs ; myid < nthread ; myid++, obj++ ) {
   switch ( flag ) {
   case NONSYM :
      rc = pthread_create(&tids[myid], &attr, worker_nonsym_mmm, obj) ;
      break ;
   case SYM :
      rc = pthread_create(&tids[myid], &attr, worker_sym_mmm, obj) ;
      break ;
   case HERM :
      rc = pthread_create(&tids[myid], &attr, worker_herm_mmm, obj) ;
      break ;
   case NONSYM_T :
      rc = pthread_create(&tids[myid], &attr, worker_nonsym_mmm_T, obj);
      break ;
   case NONSYM_H :
      rc = pthread_create(&tids[myid], &attr, worker_nonsym_mmm_H, obj);
      break ;
   }
   if ( rc != 0 ) {
      fprintf(stderr, 
           "\n fatal error, myid = %d, rc = %d from pthread_create",
           myid, rc) ;
      exit(-1) ;
   } else if ( msglvl > 2 ) {
      fprintf(stderr, "\n %% thread %d created", tids[myid]) ;
   }
}
MARKTIME(t2) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n %% CPU %8.3f : thread creation time", t2 - t1) ;
}
MARKTIME(t1) ;
for ( myid = 0 ; myid < nthread ; myid++ ) {
   pthread_join(tids[myid], &status) ;
}
MARKTIME(t2) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n %% CPU %8.3f : thread join time", t2 - t1) ;
}
FREE(tids) ;
pthread_attr_destroy(&attr) ;
}
#endif
/*
   -------------------------------------
   accumulate the rhs hand side matrices
   -------------------------------------
*/
MARKTIME(t1) ;
nent = Y->nrow * Y->ncol ;
for ( myid = 1, obj = MTmvmObjs + 1 ; 
      myid < nthread ; 
      myid++, obj++ ) {
   if ( INPMTX_IS_REAL_ENTRIES(A) ) {
      DVadd(nent, DenseMtx_entries(Y), DenseMtx_entries(obj->Y)) ;
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
      DVadd(2*nent, DenseMtx_entries(Y), DenseMtx_entries(obj->Y)) ;
   }
}
MARKTIME(t2) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, 
           "\n %% CPU %8.3f : time to accumulate rhs", t2 - t1) ;
}
/*
   ---------------------------
   release the data structures
   ---------------------------
*/
MARKTIME(t1) ;
for ( myid = 0, obj = MTmvmObjs ; myid < nthread ; myid++, obj++ ) {
   InpMtx_free(obj->A) ;
   if ( myid > 0 ) {
      DenseMtx_free(obj->Y) ;
   }
}
FREE(MTmvmObjs) ;
MARKTIME(t2) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, 
           "\n %% CPU %8.3f : time to release and free data", t2 - t1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   set up the nthread MTmvmObj data structures
   -------------------------------------------
*/
static MTmvmObj *
setup (
   InpMtx     *A,
   DenseMtx   *Y,
   double     alpha[],
   DenseMtx   *X,
   int        nthread
) {
double     *dvec ;
int        ithread, nentA, nextra, nlocal, offset ;
int        *ivec1, *ivec2 ;
MTmvmObj   *MTmvmObjs, *obj ;
/*
   ---------------------------------
   allocate nthread MTmvmObj objects
   ---------------------------------
*/
ALLOCATE(MTmvmObjs, struct _MTmvmObj, nthread) ;
for ( ithread = 0, obj = MTmvmObjs ; 
      ithread < nthread ;
      ithread++, obj++ ) {
   obj->A = InpMtx_new() ;
   if ( ithread == 0 ) {
      obj->Y = Y ;
   } else {
      obj->Y = DenseMtx_new() ;
   }
   obj->alpha[0] = alpha[0] ;
   obj->alpha[1] = alpha[1] ;
   obj->X = X ;
}
/*
   ----------------------------------------
   set up and zero the replicated Y objects
   ----------------------------------------
*/
for ( ithread = 0, obj = MTmvmObjs ; 
      ithread < nthread ;
      ithread++, obj++ ) {
   if ( ithread > 0 ) {
      DenseMtx_init(obj->Y, Y->type, Y->rowid, Y->colid, 
                    Y->nrow, Y->ncol, Y->inc1, Y->inc2) ;
      DenseMtx_zero(obj->Y) ;
   }
}
/*
   -------------------------------------
   set up the partitioned InpMtx objects
   -------------------------------------
*/
nentA  = InpMtx_nent(A)  ;
nlocal = nentA / nthread ;
nextra = nentA % nthread ;
ivec1  = InpMtx_ivec1(A) ;
ivec2  = InpMtx_ivec2(A) ;
if ( INPMTX_IS_REAL_ENTRIES(A) || INPMTX_IS_COMPLEX_ENTRIES(A) ) {
   dvec = InpMtx_dvec(A) ;
} else {
   dvec = NULL ;
}
offset = 0 ;
for ( ithread = 0, obj = MTmvmObjs ; 
      ithread < nthread ;
      ithread++, obj++ ) {
   InpMtx_init(obj->A, A->coordType, A->inputMode, 0, 0) ;
   obj->A->storageMode = A->storageMode ;
   if ( ithread < nextra ) {
      obj->A->nent = nlocal + 1 ;
   } else {
      obj->A->nent = nlocal ;
   }
   IV_init(&(obj->A->ivec1IV), obj->A->nent, ivec1 + offset) ;
   IV_init(&(obj->A->ivec2IV), obj->A->nent, ivec2 + offset) ;
   if ( INPMTX_IS_REAL_ENTRIES(A) ) {
      DV_init(&(obj->A->dvecDV), obj->A->nent, dvec + offset) ;
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(A) ) {
      DV_init(&(obj->A->dvecDV), obj->A->nent, dvec + 2*offset) ;
   }
   offset += obj->A->nent ;
}
return(MTmvmObjs) ; }

/*--------------------------------------------------------------------*/
static void *
worker_nonsym_mmm (
   void   *arg
) {
MTmvmObj   *obj ;

obj = (MTmvmObj *) arg ;
InpMtx_nonsym_mmm(obj->A, obj->Y, obj->alpha, obj->X) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
static void *
worker_sym_mmm (
   void   *arg
) {
MTmvmObj   *obj ;

obj = (MTmvmObj *) arg ;
InpMtx_sym_mmm(obj->A, obj->Y, obj->alpha, obj->X) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
static void *
worker_herm_mmm (
   void   *arg
) {
MTmvmObj   *obj ;

obj = (MTmvmObj *) arg ;
InpMtx_herm_mmm(obj->A, obj->Y, obj->alpha, obj->X) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
static void *
worker_nonsym_mmm_T (
   void   *arg
) {
MTmvmObj   *obj ;

obj = (MTmvmObj *) arg ;
InpMtx_nonsym_mmm_T(obj->A, obj->Y, obj->alpha, obj->X) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
static void *
worker_nonsym_mmm_H (
   void   *arg
) {
MTmvmObj   *obj ;

obj = (MTmvmObj *) arg ;
InpMtx_nonsym_mmm_H(obj->A, obj->Y, obj->alpha, obj->X) ;

return(NULL) ; }

/*--------------------------------------------------------------------*/
