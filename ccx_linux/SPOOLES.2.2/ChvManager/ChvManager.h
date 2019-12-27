/*  ChvManager.h  */

#include "../Chv.h"
#include "../Lock.h"
 
/*--------------------------------------------------------------------*/
/*
*/
typedef struct _ChvManager  ChvManager ;
struct _ChvManager {
   Chv       *head           ;
   Lock      *lock           ;
   int       mode            ;
   int       nactive         ;
   int       nbytesactive    ;
   int       nbytesrequested ;
   int       nbytesalloc     ;
   int       nrequests       ;
   int       nreleases       ;
   int       nlocks          ;
   int       nunlocks        ;
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
ChvManager *
ChvManager_new ( 
   void 
) ;
/*
   -----------------------
   set the default fields
 
   created -- 98may02, cca
   -----------------------
*/
void
ChvManager_setDefaultFields (
   ChvManager   *manager
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 98may02, cca
   --------------------------------------------------
*/
void
ChvManager_clearData (
   ChvManager   *manager
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 98may02, cca
   ------------------------------------------
*/
void
ChvManager_free (
   ChvManager   *manager
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
                    only threads in this process.
   lockflag = 2 --> mutex lock is allocated and it can synchronize
                    only threads in this and other processes.

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
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
----- methods found in util.c ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------
   return a pointer to a Chv object that has 
   been initialized with the input parameters
 
   created -- 98may02, cca
   ------------------------------------------
*/
Chv *
ChvManager_newObject (
   ChvManager   *manager,
   int           id,
   int           nD,
   int           nL,
   int           nU,
   int           symflag
) ;
/*
   ------------------------------------------
   return a pointer to a Chv object that has 
   been initialized with the input parameters
 
   created -- 98may02, cca
   ------------------------------------------
*/
Chv *
ChvManager_newObjectOfSizeNbytes (
   ChvManager   *manager,
   int           nbytesNeeded
) ;
/*
   -----------------------
   release a Chv instance
 
   created -- 98may02, cca
   -----------------------
*/
void
ChvManager_releaseObject (
   ChvManager   *manager,
   Chv          *chv
) ;
/*
   ------------------------------
   release a list of Chv objects
 
   created -- 98may02, cca
   ------------------------------
*/
void
ChvManager_releaseListOfObjects (
   ChvManager   *manager,
   Chv          *head
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
ChvManager_writeForHumanEye (
   ChvManager   *manager,
   FILE     *fp
) ;
/*--------------------------------------------------------------------*/
