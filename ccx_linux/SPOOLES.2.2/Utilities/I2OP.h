/*  I2OP.h  */

/*--------------------------------------------------------------------*/

typedef struct _I2OP I2OP ;
struct _I2OP {
   int    value0  ;
   int    value1  ;
   void   *value2 ;
   I2OP   *next   ;
} ;

#define I2OP_NULL     0
#define I2OP_FORWARD  1
#define I2OP_BACKWARD 2

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   initializer. 
   create and return an array of n I2OP structures.
   the structures are linked in one of three ways.
   flag = 0 (I2OP_NULL)     --> ip->next = NULL
   flag = 1 (I2OP_FORWARD)  --> ip->next = successor in list
   flag = 2 (I2OP_BACKWARD) --> ip->next = predecessor in list
   
   created -- 98feb06, cca
   ---------------------------------------------------------
*/
I2OP *
I2OP_init ( 
   int   n, 
   int   flag 
) ;
/*
   ---------------------------------------------------------
   initializer. 
   create and return an array of n I2OP structures.
   the structures are linked in one of three ways.
   flag = 0 (I2OP_NULL)     --> ip->next = NULL
   flag = 1 (I2OP_FORWARD)  --> ip->next = successor in list
   flag = 2 (I2OP_BACKWARD) --> ip->next = predecessor in list
   
   created -- 98feb06, cca
   ---------------------------------------------------------
*/
void
I2OP_initStorage ( 
   int   n, 
   int   flag,
   I2OP   *base
) ;
/*
   -----------------------------------------------
   free the storage for an array of I2OP structures,
   must have been allocated by I2OP_init

   created -- 98feb06, cca
   -----------------------------------------------
*/
void
I2OP_free ( 
   I2OP   *ip
) ;
/*
   ----------------------------------
   purpose -- to print out a I2OP list
   
   created -- 98feb06, cca
   ----------------------------------
*/
void
I2OP_fprintf ( 
   FILE   *fp, 
   I2OP    *elem 
) ;
/*--------------------------------------------------------------------*/
