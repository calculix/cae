/*  util.c  */

#include "../Lock.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------
   lock the lock

   created -- 97aug22, cca
   -----------------------
*/
void
Lock_lock ( 
   Lock   *lock 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( lock == NULL ) {
   fprintf(stderr, "\n fatal error in Lock_lock(%p)"
           "\n bad input\n", lock) ;
   exit(-1) ;
}
/*
fprintf(stdout, "\n inside Lock_lock()") ;
fflush(stdout) ;
*/
#if THREAD_TYPE == TT_SOLARIS
#if MYDEBUG > 0
fprintf(stdout, "\n thread %d, lock %p : Lock_lock : solaris locking", 
        thr_self(), lock->mutex) ;
fflush(stdout) ;
#endif
mutex_lock(lock->mutex) ;
#if MYDEBUG > 0
fprintf(stdout, 
        "\n thread %d, lock %p : Lock_lock : solaris lock done", 
        thr_self(), lock->mutex) ;
fflush(stdout) ;
#endif
#endif

#if THREAD_TYPE == TT_POSIX
#if MYDEBUG > 0
fprintf(stdout, "\n Lock_lock : posix locking") ;
fflush(stdout) ;
#endif
pthread_mutex_lock(lock->mutex) ;
#if MYDEBUG > 0
fprintf(stdout, "\n Lock_lock : posix lock done") ;
fflush(stdout) ;
#endif
#endif
lock->nlocks++ ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   unlock the lock

   created -- 97aug22, cca
   -----------------------
*/
void
Lock_unlock ( 
   Lock   *lock 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( lock == NULL ) {
   fprintf(stderr, "\n fatal error in Lock_unlock(%p)"
           "\n bad input\n", lock) ;
   exit(-1) ;
}
lock->nunlocks++ ;
/*
fprintf(stdout, "\n inside Lock_unlock()") ;
fflush(stdout) ;
*/
#if THREAD_TYPE == TT_SOLARIS
#if MYDEBUG > 0
fprintf(stdout, 
        "\n thread %d, lock %p : Lock_unlock : solaris unlocking", 
        thr_self(), lock->mutex) ;
fflush(stdout) ;
#endif
mutex_unlock(lock->mutex) ;
#if MYDEBUG > 0
fprintf(stdout, 
        "\n thread %d, lock %p : Lock_unlock : solaris unlocking done", 
        thr_self(), lock->mutex) ;
fflush(stdout) ;
#endif
#endif

#if THREAD_TYPE == TT_POSIX
#if MYDEBUG > 0
fprintf(stdout, "\n Lock_unlock : posix unlocking") ;
fflush(stdout) ;
#endif
pthread_mutex_unlock(lock->mutex) ;
#if MYDEBUG > 0
fprintf(stdout, "\n Lock_unlock : posix unlocking done") ;
fflush(stdout) ;
#endif
#endif

return ; }

/*--------------------------------------------------------------------*/
