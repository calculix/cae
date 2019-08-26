/*  Ideq.h  */

#include "../IV.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   Ideq -- dequeue with integer ids

   maxsize -- maxmimum size of the deq
   head    -- head of the list
   tail    -- tail of the list
   iv      -- IV object to manage dequeue

   created -- 96jun06, cca
   ------------------------------------------------------------------
*/
typedef struct _Ideq   Ideq ;
struct _Ideq {
   int    maxsize ;
   int    head    ;
   int    tail    ;
   IV     iv      ;
} ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
------ methods found in basics.c  --------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------
   create and return a new instance of the Ideq object

   created -- 96jun06, cca
   -----------------------------------------------------
*/
Ideq *
Ideq_new (
   void
) ;
/*
   -------------------------------------------
   set the default fields for an Ideq object

   created -- 96jun06, cca
   -------------------------------------------
*/
void
Ideq_setDefaultFields (
   Ideq   *deq
) ;
/*
   -----------------------
   clear the data fields

   created -- 96jun06, cca
   -----------------------
*/
void
Ideq_clearData (
   Ideq   *deq
) ;
/*
   -----------------------
   free the Ideq object

   created -- 96jun06, cca
   -----------------------
*/
void
Ideq_free (
   Ideq   *deq
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
------ methods found in resize.c  --------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------
   resize the deque
   if the new size is large enough then
      copy the old data
      return 1
   else
      error, return -1
   endif

   created -- 96jun06, cca
   ------------------------------------
*/
int
Ideq_resize (
   Ideq   *deq,
   int    newsize
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
------ methods found in util.c  ----------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   clear the dequeue,
  
   created -- 96jun06, cca
   -----------------------
*/
void
Ideq_clear (
   Ideq   *deq
) ;
/*
   ---------------------------------
   return the head of the dequeue,
   return -1 if the dequeue is empty
  
   created -- 96jun06, cca
   ---------------------------------
*/
int
Ideq_head (
   Ideq   *deq
) ;
/*
   ------------------------------------------
   return and remove the head of the dequeue,
   return -1 if the dequeue is empty
  
   created -- 96jun06, cca
   ------------------------------------------
*/
int
Ideq_removeFromHead (
   Ideq   *deq
) ;
/*
   ---------------------------------------
   insert value at head of dequeue
   return value
     1 --> value inserted
    -1 --> no room in dequeue, must resize
  
   created -- 96jun06, cca
   ---------------------------------------
*/
int
Ideq_insertAtHead (
   Ideq   *deq,
   int    val 
) ;
/*
   ---------------------------------
   return the tail of the dequeue,
   return -1 if the dequeue is empty
 
   created -- 96jun06, cca
   ---------------------------------
*/
int
Ideq_tail (
   Ideq   *deq
) ;
/*
   ------------------------------------------
   return and remove the tail of the dequeue,
   return -1 if the dequeue is empty
 
   created -- 96jun06, cca
   ------------------------------------------
*/
int
Ideq_removeFromTail (
   Ideq   *deq
) ;
/*
   ---------------------------------------
   insert value at tail of dequeue
   return value
     1 --> value inserted
    -1 --> no room in dequeue, must resize
 
   created -- 96jun06, cca
   ---------------------------------------
*/
int
Ideq_insertAtTail (
   Ideq   *deq,
   int    val
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
------ methods found in IO.c  ------------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------
   purpose -- write the contents of the dequeue to a file
              in a human readable format
 
   created -- 98feb11, cca
   ------------------------------------------------------
*/
void
Ideq_writeForHumanEye (
   Ideq   *dequeue,
   FILE   *fp
) ;
/*--------------------------------------------------------------------*/
