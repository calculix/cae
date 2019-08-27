/*  Lock.h  */

#include "../cfiles.h"

#define TT_NONE    0
#define TT_SOLARIS 1
#define TT_POSIX   2
 
#define THREAD_TYPE TT_POSIX
 
#if THREAD_TYPE == TT_SOLARIS
#include <thread.h>
#include <synch.h>
#endif
#if THREAD_TYPE == TT_POSIX
#include <pthread.h>
#endif

#define NO_LOCK                  0
#define LOCK_IN_PROCESS          1
#define LOCK_OVER_ALL_PROCESSES  2

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   this structure contains a lock,
   presently solaris and posix thread packages are supported

   mutex    -- pointer to a lock
   nlocks   -- number of locks
   nunlocks -- number of unlocks

   created -- 97aug22, cca
   ---------------------------------------------------------
*/
typedef struct _Lock   Lock ;
struct _Lock {
#if THREAD_TYPE == TT_SOLARIS
   mutex_t   *mutex ;
#endif
#if THREAD_TYPE == TT_POSIX
   pthread_mutex_t   *mutex ;
#endif
#if THREAD_TYPE == TT_NONE
   void   *mutex ;
#endif
   int    nlocks   ;
   int    nunlocks ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method found in basics.c -----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor

   created -- 97aug22, cca
   -----------------------
*/
Lock *
Lock_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields

   created -- 97aug22, cca
   -----------------------
*/
void
Lock_setDefaultFields (
   Lock   *lock
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage

   created -- 97aug22, cca
   --------------------------------------------------
*/
void
Lock_clearData ( 
   Lock   *lock 
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data

   created -- 97aug22, cca
   ------------------------------------------
*/
void
Lock_free ( 
   Lock   *lock 
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method found in init.c -------------------------------------------
------------------------------------------------------------------------
*/
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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- method found in util.c -------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   lock the lock

   created -- 97aug22, cca
   -----------------------
*/
void
Lock_lock ( 
   Lock   *lock 
) ;
/*
   -----------------------
   unlock the lock

   created -- 97aug22, cca
   -----------------------
*/
void
Lock_unlock ( 
   Lock   *lock 
) ;
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
