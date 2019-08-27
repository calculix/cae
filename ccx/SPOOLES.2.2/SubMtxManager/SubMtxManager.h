/*  SubMtxManager.h  */

#include "../SubMtx.h"
#include "../Lock.h"
 
/*--------------------------------------------------------------------*/
/*
*/
typedef struct _SubMtxManager  SubMtxManager ;
struct _SubMtxManager {
   SubMtx   *head           ;
   Lock     *lock           ;
   int      mode            ;
   int      nactive         ;
   int      nbytesactive    ;
   int      nbytesrequested ;
   int      nbytesalloc     ;
   int      nrequests       ;
   int      nreleases       ;
   int      nlocks          ;
   int      nunlocks        ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in basics.c ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 98may02, cca
   -----------------------
*/
SubMtxManager *
SubMtxManager_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may02, cca
   -----------------------
*/
void
SubMtxManager_setDefaultFields (
   SubMtxManager   *manager
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may02, cca
   --------------------------------------------------
*/
void
SubMtxManager_clearData (
   SubMtxManager   *manager
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may02, cca
   ------------------------------------------
*/
void
SubMtxManager_free (
   SubMtxManager   *manager
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in init.c ------------------------------------------
------------------------------------------------------------------------
*/
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
SubMtxManager_init (
   SubMtxManager   *manager,
   int             lockflag,
   int             mode
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   return a pointer to a SubMtx object that has 
   been initialized with the input parameters
 
   created -- 98may02, cca
   -----------------------------------------------
*/
SubMtx *
SubMtxManager_newObjectOfSizeNbytes (
   SubMtxManager   *manager,
   int             nbytesNeeded
) ;
/*
   ----------------------------
   release a SubMtx instance
 
   created -- 98may02, cca
   ----------------------------
*/
void
SubMtxManager_releaseObject (
   SubMtxManager   *manager,
   SubMtx          *mtx
) ;
/*
   -----------------------------------
   release a list of SubMtx objects
 
   created -- 98may02, cca
   -----------------------------------
*/
void
SubMtxManager_releaseListOfObjects (
   SubMtxManager   *manager,
   SubMtx          *head
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in IO.c --------------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form
 
   created -- 98may02, cca
   ----------------------------------------
*/
void
SubMtxManager_writeForHumanEye (
   SubMtxManager   *manager,
   FILE            *fp
) ;
/*--------------------------------------------------------------------*/
