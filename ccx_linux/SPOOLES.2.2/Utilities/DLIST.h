#ifndef _DLIST_
#define _DLIST_
/*
   -------------------------------------------------------------------
   doubly linked list macros.
 
   all double linked list structures have the two prev and next fields
   struct delem {
      whatever other fields 
      struct delem   *prev, *next ; } ;
 
   all lists have a sentinel element which is the head of the list.
   when an element is not in a list, its link(s) point to itself
 
   example : allocate a set of doubly linked list elements
             and insert into a list
 
   struct delem   Head, *head, *elems ;
   head = &Head ; head->prev = head->next = head ;
   ALLOCATE(elems, struct delem, nelem, test) ;
   for ( ielem = 0, elem = elems ; ielem < nelem ; ielem++, elem++) {
      elem->prev = elem->next = elem ;
      DLIST_TAIL_INSERT(head, elem) ; }
 
   DLIST_DELETE      -- delete an element from a doubly linked list
   DLIST_TAIL_INSERT -- insert an element into a doubly linked list
                        at the tail of the list
   DLIST_HEAD_INSERT -- insert an element into a doubly linked list
                       at the head of the list
   DLIST_SWAP_HEADS  -- swaps the heads of two doubly linked lists
   DLIST_MERGE       -- merges two doubly linked lists
   -------------------------------------------------------------------
*/
#define DLIST_DELETE(elem) \
(elem)->prev->next = (elem)->next ; \
(elem)->next->prev = (elem)->prev ; \
(elem)->prev = (elem)->next = (elem) ;
 
#define DLIST_TAIL_INSERT(head, elem) \
((elem)->prev = (head)->prev)->next = (elem) ; \
((head)->prev = elem)->next = (head) ;
 
#define DLIST_HEAD_INSERT(head, elem) \
((elem)->next = (head)->next)->prev = (elem) ; \
((head)->next = elem)->prev = (head) ;
 
#define DLIST_SWAP_HEADS(head1, head2, type) \
if ( (head1)->next == head1 ) { \
   if ( (head2)->next != head2 ) { \
      ((head1)->next = (head2)->next)->prev  \
          = ((head1)->prev = (head2)->prev)->next = head1 ; \
      (head2)->prev = (head2)->next = head2 ; } } \
else if ( (head2)->next == head2 ) { \
   ((head2)->next = (head1)->next)->prev \
       = ((head2)->prev = (head1)->prev)->next = head2 ; \
   (head1)->prev = (head1)->next = head1 ; } \
else { \
   type *temp = head1->next ; \
   ((head1)->next = (head2)->next)->prev = head1 ; \
   ((head2)->next = temp)->prev = head2 ; \
   temp = head1->prev ; \
   ((head1)->prev = (head2)->prev)->next = head1 ; \
   ((head2)->prev = temp)->next = head2 ; }
 
#define DLIST_MERGE(head1, head2) \
if ( (head2)->next != head2 ) { \
   ((head2)->next->prev = (head1)->prev)->next = (head2)->next ; \
   ((head1)->prev = (head2)->prev)->next = head1 ; \
   (head2)->prev  = (head2)->next = head2 ; }

#endif
