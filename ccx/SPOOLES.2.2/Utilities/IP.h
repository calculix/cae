/*  IP.h  */

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   the IP structure contains an Int field and an IP* field.
   it is the simplest singly linked list element, useful at times.
   ---------------------------------------------------------------
*/
typedef struct _IP IP ;
struct _IP {
   int   val   ;
   IP    *next ;
} ;
/*--------------------------------------------------------------------*/

#define IP_NULL     0
#define IP_FORWARD  1
#define IP_BACKWARD 2

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   purpose -- to print out a IP list

   created -- 95sep22, cca
   ---------------------------------
*/
void
IP_fprintf ( 
   FILE   *fp, 
   IP     *ip 
) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- to write out an integer list with eighty column lines

   input --

      fp     -- file pointer, must be formatted and write access
      ip     -- head of list
      column -- present column
  
   return value -- present column

   created -- 95sep22, cca
   ------------------------------------------------------------------
*/
int
IP_fp80 ( 
   FILE   *fp, 
   IP     *ip, 
   int    column 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   initializer. 
   create and return an array of n IP structures.
   the structures are linked in one of three ways.
   flag = 0 (IP_NULL)     --> ip->next = NULL
   flag = 1 (IP_FORWARD)  --> ip->next = successor in list
   flag = 2 (IP_BACKWARD) --> ip->next = predecessor in list

   created -- 95sep22, cca
   ---------------------------------------------------------
*/
IP *
IP_init ( 
   int   n, 
   int   flag 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   free the storage for an array of IP structures,
   must have been allocated by IP_init
 
   created -- 95sep22, cca
   -----------------------------------------------
*/
void
IP_free ( 
   IP   *ip
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   merge two lists in ascending order

   created -- 95sep22, cca
   ----------------------------------
*/
IP *
IP_mergeUp ( 
   IP   *ip1, 
   IP   *ip2 
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- to sort a singly linked list in
              ascending order using a radix sort

   created -- 95sep22, cca
   ---------------------------------------------
*/
IP *
IP_radixSortUp ( 
   IP   *ip 
) ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to sort a singly linked list in
              descending order using a radix sort

   created -- 95sep22, cca
   ----------------------------------------------
*/
IP *
IP_radixSortDown ( 
   IP   *ip 
) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   sort a list in ascending order using merge sort

   created -- 95sep22, cca
   -----------------------------------------------
*/
IP *
IP_mergeSortUp ( 
   IP   *ip0 
) ;
/*--------------------------------------------------------------------*/
