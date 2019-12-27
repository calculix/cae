/*  IP.c  */

#include "../Utilities.h"

#define MYDEBUG 0

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
) {
if ( fp != NULL && ip != NULL ) {
   int   i = 0 ;
   while ( ip != NULL ) {
      if ( i % 16 == 0 ) fprintf(fp, "\n ") ;
      fprintf(fp, " %4d", ip->val) ;
      ip = ip->next ;
      i++ ;
   }
}
return ; }

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
) {
if ( fp != NULL && ip != NULL ) {
   int    inum, nchar, pow ;
   while ( ip != NULL ) {
      inum = ip->val ;
      if ( inum < 0 ) {
         inum = -inum ; 
         nchar = 3 ; 
      } else {
         nchar = 2 ; 
      }
      for ( pow = 10 ; inum >= pow ; pow *= 10 ) {
         nchar++ ; 
      }
      if ( (column += nchar) >= 80 ) {
         fprintf(fp, "\n") ;
         column = nchar ; 
      }
      fprintf(fp, " %d", ip->val) ; 
      ip = ip->next ;
   }
}
return(column) ; }

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
) {
IP    *base = NULL ;
if ( n > 0 ) { 
   if ( flag != IP_NULL 
     && flag != IP_FORWARD 
     && flag != IP_BACKWARD ) {
      fprintf(stderr, "\n fatal error in IPinit, invalid data"
      "\n n = %d, flag = %d"
      "\n flag must be 0 (IP_NULL), 1 (IP_FORWARD) or 2(IP_BACKWARD)\n",
      n, flag) ;
      exit(-1) ;
   } else {
      int   i ;
      IP    *head, *ip, *tail ;
      ALLOCATE(base, struct _IP, n) ;
      switch ( flag ) {
      case IP_FORWARD :
         head = NULL ;
         for ( i = n - 1, ip = base + i ; i >= 0 ; i--, ip-- ) {
            ip->next = head ;
            ip->val  = 0 ;
            head     = ip ;
         }
         break ;
      case IP_BACKWARD :
         head = tail = base + n - 1 ;
         head->val = 0 ;
         for ( i = n - 2, ip = head + i ; i >= 0 ; i--, ip-- ) {
            tail->next = ip ;
            ip->val = 0 ;
         }
         tail->next = NULL ;
         break ;
      case IP_NULL :
         for ( i = 0, ip = base ; i < n ; i++, ip++ ) {
            ip->val  = 0 ;
            ip->next = NULL ;
         }
         break ;
      }
   }
}
return(base) ; }

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
) {
if ( ip != NULL ) {
   FREE(ip) ;
}
return ; }

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
) {
IP   *head, *tail ;
/*
   -------------------
   check for NULL list
   -------------------
*/
if ( ip1 == NULL ) {
   head = ip2 ;
} else if ( ip2 == NULL ) {
   head = ip1 ;
} else {
/*
   ------------------------------------------
   neither list is NULL, assign first element
   ------------------------------------------
*/
   if ( ip2->val < ip1->val ) {
      head = tail = ip2 ;
      ip2 = ip2->next ;
   } else {
      head = tail = ip1 ;
      ip1 = ip1->next ;
   }
/*
   --------------------------------------
   merge the lists until one is exhausted
   --------------------------------------
*/
   while ( ip1 != NULL && ip2 != NULL ) {
      if ( ip2->val < ip1->val ) {
         tail->next = ip2 ;
         tail = ip2 ;
         ip2 = ip2->next ;
      } else {
         tail->next = ip1 ;
         tail = ip1 ;
         ip1 = ip1->next ;
      }
   }
/*
   ----------------------------------
   add the remaining list to the tail
   ----------------------------------
*/
   if ( ip1 == NULL ) {
      tail->next = ip2 ;
   } else {
      tail->next = ip1 ;
   }
}

return(head) ; }

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
) {
int   b1, b2, d, dneg, dpos, i, j, negmin, posmax ;
IP    *head, *next, *tail ;
IP    *poshead, *neghead, *zerohead ;
IP    *postail, *negtail, *zerotail ;
#define BASE 10
IP    *heads[BASE], *tails[BASE] ;
/*
   --------------------------------------------------------
   split the list into negative, zero and positive sublists
   --------------------------------------------------------
*/
poshead = postail = neghead = negtail = zerohead = NULL ;
zerotail = NULL ;
posmax = negmin = 0 ;
while ( ip != NULL ) {
   next = ip->next ;
   if ( ip->val > 0 ) {
      ip->next = poshead, poshead = ip ;
      if ( posmax < ip->val ) {
         posmax = ip->val ;
      }
   } else if ( ip->val < 0 ) {
      ip->next = neghead, neghead = ip ;
      if ( negmin > ip->val ) {
         negmin = ip->val ;
      }
   } else {
      if ( zerohead == NULL ) {
         zerotail = ip ;
      }
      ip->next = zerohead, zerohead = ip ;
   }
   ip = next ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n positive list :") ;
IP_fp80(stdout, poshead, 16) ;
fprintf(stdout, "\n zero list :") ;
IP_fp80(stdout, zerohead, 16) ;
fprintf(stdout, "\n negative list :") ;
IP_fp80(stdout, neghead, 16) ;
fflush(stdout) ;
#endif
/*
   ---------------
   find the limits
   ---------------
*/
dpos = 0 ; 
while ( posmax > 0 ) {
   dpos++ ;
   posmax /= 10 ;
}
negmin = - negmin ;
dneg = 0 ; 
while ( negmin > 0 ) {
   dneg++ ;
   negmin /= 10 ;
}
if ( dpos > dneg ) {
   d = dpos ;
} else {
   d = dneg ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n dneg %d, dpos %d, d %d", dneg, dpos, d) ;
fflush(stdout) ;
#endif
/*
   ----------------------
   sort the positive list
   ----------------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n sorting the positive list") ;
#endif
for ( i = 0 ; i < BASE ; i++ ) {
   heads[i] = tails[i] = NULL ;
}
b1 = 1 ;
for ( i = 0 ; i < dpos ; i++ ) {
   b2 = BASE * b1 ;
   ip = poshead ; poshead = NULL ;
#if MYDEBUG > 0
   fprintf(stdout, "\n b1 %d, b2 %d", b1, b2) ;
#endif
   while ( ip != NULL ) {
      next = ip->next ;
      j = (ip->val % b2) / b1 ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    ip->val %d, j %d", ip->val, j) ;
#endif
      if ( heads[j] == NULL ) {
         heads[j] = ip ;
      } else {
         tails[j]->next = ip ;
      }
      tails[j] = ip ;
      ip = next ;
   }
   for ( j = 0 ; j < BASE ; j++ ) {
      if ( heads[j] != NULL ) {
         if ( poshead == NULL ) {
            poshead = heads[j] ;
         } else {
            postail->next = heads[j] ;
         }
         postail = tails[j] ;
         heads[j] = tails[j] = NULL ;
      }
   }
   postail->next = NULL ;
   b1 = b2 ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n positive list") ;
IP_fprintf(stdout, poshead) ;
#endif
/*
   ----------------------
   sort the negative list
   ----------------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n sorting the negative list") ;
#endif
b1 = 1 ;
for ( i = 0 ; i < dneg ; i++ ) {
   b2 = BASE * b1 ;
#if MYDEBUG > 0
   fprintf(stdout, "\n b1 %d, b2 %d", b1, b2) ;
#endif
   ip = neghead ; neghead = NULL ;
   while ( ip != NULL ) {
      next = ip->next ;
      j = ((-ip->val) % b2) / b1 ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    ip->val %d, j %d", ip->val, j) ;
#endif
      if ( heads[j] == NULL ) {
         heads[j] = ip ;
      } else {
         tails[j]->next = ip ;
      }
      tails[j] = ip ;
      ip = next ;
   }
   for ( j = 0 ; j < BASE ; j++ ) {
      if ( heads[j] != NULL ) {
         if ( neghead == NULL ) {
            neghead = heads[j] ;
         } else {
            negtail->next = heads[j] ;
         }
         negtail = tails[j] ;
         heads[j] = tails[j] = NULL ;
      }
   }
   negtail->next = NULL ;
   b1 = b2 ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n negative list") ;
IP_fprintf(stdout, neghead) ;
#endif
/*
   ---------------------------
   concatenate the three lists
   ---------------------------
*/
head = tail = ip = neghead ;
while ( ip != NULL ) {
   next = ip->next ;
   ip->next = head ;
   head = ip ;
   ip = next ;
}
if ( tail != NULL ) {
   tail->next = NULL ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n 1. head = %p, tail = %p", head, tail) ;
#endif
if ( zerohead != NULL ) {
   if ( tail != NULL ) {
      tail->next = zerohead ;
   } else {
      head = zerohead ;
   }
   tail = zerotail ;
}
if ( tail != NULL ) {
   tail->next = NULL ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n 2. head = %p, tail = %p", head, tail) ;
#endif
if ( poshead != NULL ) {
   if ( tail != NULL ) {
      tail->next = poshead ;
   } else {
      head = poshead ;
   }
   tail = postail ;
}
if ( tail != NULL ) {
   tail->next = NULL ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n 3. head = %p, tail = %p", head, tail) ;
IP_fprintf(stdout, head) ;
#endif

return(head) ; }
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
) {
IP   *ip1 = NULL ;
if ( ip != NULL ) {
   IP   *ip0 ;

   for ( ip0 = ip ; ip0 != NULL ; ip0 = ip0->next ) {
       ip0->val = -ip0->val ;
   }
   ip1 = IP_radixSortUp(ip) ;
   for ( ip0 = ip1 ; ip0 != NULL ; ip0 = ip0->next ) {
       ip0->val = -ip0->val ;
   }
}
return(ip1) ; }

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
) {
int   i, m, m1 ;
IP    *ip, *ip1, *ip2, *prev ;
/*
   ----------------------------------------
   count the number of elements in the list
   ----------------------------------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n inside IP_mergeSortUp :") ;
IP_fp80(stdout, ip0, 25) ;
fflush(stdout) ;
#endif
for ( ip = ip0, m = 0 ; ip != NULL ; ip = ip->next ) {
   m++ ;
}
if ( m <= 1 ) {
   return(ip0) ;
} else {
   m1 = m / 2 ;
#if MYDEBUG > 0
   fprintf(stdout, "\n m = %d, m1 = %d, m2 = %d", m, m1, m-m1) ;
   fflush(stdout) ;
#endif
   for ( i = 1, ip = ip0, prev = NULL ; i < m1 ; i++ ) {
      prev = ip ;
      ip = ip->next ;
   }
   ip2 = ip->next ;
   ip->next = NULL ;
   ip1 = ip0 ;
#if MYDEBUG > 0
   fprintf(stdout, "\n calling IP_mergeSortUp :") ;
   IP_fp80(stdout, ip1, 13) ;
   fflush(stdout) ;
#endif
   ip1 = IP_mergeSortUp(ip1) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n return from IP_mergeSortUp :") ;
   IP_fp80(stdout, ip1, 13) ;
   fprintf(stdout, "\n calling IPmergeSortUp :") ;
   IP_fp80(stdout, ip2, 13) ;
#endif
   ip2 = IP_mergeSortUp(ip2) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n return from IP_mergeSortUp :") ;
   IP_fp80(stdout, ip2, 13) ;
   fprintf(stdout, "\n calling IP_mergeUp :") ;
   fprintf(stdout, "\n first list") ;
   IP_fp80(stdout, ip1, 13) ;
   fprintf(stdout, "\n second list") ;
   IP_fp80(stdout, ip2, 13) ;
#endif
   ip  = IP_mergeUp(ip1, ip2) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n return from IP_mergeUp, sorted list : ") ;
   IP_fp80(stdout, ip, 40) ;
   fflush(stdout) ;
#endif
   return(ip) ; 
}
}
/*--------------------------------------------------------------------*/
