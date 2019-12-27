/*  memory.h  */

/*====================================================================*/
/*
   -------------------------------------------------------
   flag to turn on the memory functions debugging commands
   -------------------------------------------------------
*/
#define MEMORY_DEBUG 0
/*
   ----------------------------------------------------------
   my memory allocation function

      ptr   -- variable to be given the address of the memory
      type  -- type of data, can be a struct
      count -- number of data elements
      proc  -- procedure name for error message

   e.g.,   

   int           *indices ;
   double        *entries ;
   struct elem   *elems
   ALLOCATE(indices, int, nindices) ;
   ALLOCATE(entries, int, nrow*ncol) ;
   ALLOCATE(elems, struct elem, nelem) ;

   created -- 95sep22, cca
   ----------------------------------------------------------
*/
#define ALLOCATE(ptr, type, count) \
if ( (count) > 0 ) { \
   if ( (ptr = (type *)malloc((unsigned long)((count)*sizeof(type)))) \
        == NULL ) {\
      fprintf(stderr, \
              "\n ALLOCATE failure : bytes %d, line %d, file %s", \
              (count)*sizeof(type), __LINE__, __FILE__) ; \
      exit(-1) ; } \
   else if ( MEMORY_DEBUG > 0 ) { \
      fprintf(stderr, \
              "\n ALLOCATE : address %p, bytes %d, line %d, file %s", \
              ptr, (count)*sizeof(type), __LINE__, __FILE__) ; } } \
else if ( (count) == 0 ) { \
   ptr = NULL ; } \
else { \
   fprintf(stderr, \
           "\n ALLOCATE error : bytes %d, line %d, file %s", \
           (count)*sizeof(type), __LINE__, __FILE__) ; \
   exit(-1) ; }
/*
   --------------------------------------------------------
   my function to free memory, all it does is check to see 
   if the pointer is NULL, calls system routine if not NULL
   --------------------------------------------------------
*/
#if MEMORY_DEBUG > 0
#define FREE(ptr) \
if ( (ptr) != NULL ) { \
   fprintf(stderr, "\n FREE, line %d, file %s : address %p", \
           __LINE__, __FILE__, ptr) ; \
   free((char *) (ptr)) ; \
   (ptr) = NULL ; }
#else
#define FREE(ptr) \
if ( (ptr) != NULL ) { \
   free((char *) (ptr)) ; \
   (ptr) = NULL ; }
#endif
/*====================================================================*/
