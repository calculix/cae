/*  init.c  */

#include "../Lock.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- basic initializer
 
   lockflag -- flag to specify lock status
      SOLARIS:
      lockflag = 0 --> mutex lock is not allocated or initialized
      lockflag = 1 --> mutex lock is allocated and it can synchronize
                       only threads in this process.
      lockflag = 2 --> mutex lock is allocated and it can synchronize
                       only threads in this and other processes.
      POSIX:
      lockflag = 0 --> mutex lock is not allocated or initialized
      lockflag = 1 --> mutex lock is allocated and it can synchronize
                       only threads in this process.
 
   created -- 97aug22, cca
   ------------------------------------------------------------------
*/
void
Lock_init (
   Lock   *lock,
   int    lockflag
) {
if ( lockflag > 0 ) {
/*
   -----------------
   allocate the lock
   -----------------
*/
#if THREAD_TYPE == TT_SOLARIS
   ALLOCATE(lock->mutex, mutex_t, 1) ;
   if ( lockflag == 1 ) {
      mutex_init(lock->mutex, USYNC_THREAD, NULL) ;
   } else if ( lockflag == 2 ) {
      mutex_init(lock->mutex, USYNC_PROCESS, NULL) ;
   }
#endif
#if THREAD_TYPE == TT_POSIX
   ALLOCATE(lock->mutex, pthread_mutex_t, 1) ;
   pthread_mutex_init(lock->mutex, NULL) ;
#endif
}
return ; }

/*--------------------------------------------------------------------*/
