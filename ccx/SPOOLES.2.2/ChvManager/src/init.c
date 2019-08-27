/*  init.c  */

#include "../ChvManager.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   simple initializer

   lockflag = 0 --> mutex lock is not allocated or initialized
   lockflag = 1 --> mutex lock is allocated and it can synchronize
                    omly threads in this process.
   lockflag = 2 --> mutex lock is allocated and it can synchronize
                    omly threads in this and other processes.

    mode = 0 --> free object and storage on release
    mode = 1 --> recycle object and storage on release
                        
   created -- 98may02, cca
   ---------------------------------------------------------------
*/
void
ChvManager_init (
   ChvManager   *manager,
   int          lockflag,
   int          mode
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( manager == NULL 
   || lockflag < 0 || lockflag > 2 
   || mode < 0 || mode > 1 ) {
   fprintf(stderr, "\n fatal error in ChvManager_init(%p,%d,%d)"
          "\n bad input\n", manager, lockflag, mode) ;
   exit(-1) ;
}
/*
   --------------------------------------------------
   clear any previous data and set the default fields
   --------------------------------------------------
*/
ChvManager_clearData(manager) ;
if ( lockflag > 0 ) {
/*
   ---------------------------
   initialize the mutex object
   ---------------------------
*/
   manager->lock = Lock_new() ;
   Lock_init(manager->lock, lockflag) ;
}
/*
   ------------
   set the mode
   ------------
*/
manager->mode = mode ;

return ; }

/*--------------------------------------------------------------------*/
