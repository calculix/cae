/*  IVL.h  */

#include "../IV.h"
#include "../cfiles.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   The IVL object stores and manipulates Lists of Int Vectors,
   so the acronym for Integer Vector Lists.

   The most common use of an IVL object is to represent a graph
   or the adjacency structure of a matrix.

   The IVL object supports the following functionality.
      1. initialization, # of lists must be known.
      2. clearing data and free'ing the object
      3. set the size of a list and/or its entries
      4. un-setting a list or releasing its data 
      5. returning the min or max entry in all the lists
      6. returning the maximum size of a list

   The number of lists is fixed on initialization, but the number
   of entries in any list or the total size of the lists need not
   be known at initialization. Storage for the lists is handled in
   one of three ways.
      1. chunked (IVL_CHUNKED = 1)
         A chunk of data is allocated by the object and lists point
         to data in the chunk. When the free space is not large
         enough to contain a new list, another chunk is allocated.
         The object keeps track of the chunks and free's all the
         storage when requested.
      2. solo (IVL_SOLO = 2)
         Each list is allocated separately using the IVinit() or
         IVinit2() function. When requested, each list is free'd
         using the IVfree() function.
      3. unknown (IVL_UNKNOWN = 3)
         Each list points to storage somewhere but the object does
         not keep track of any storage to be free'd. The application
         that gave rise to this storage mode has a subgraph "share"
         list storage with the larger graph that contains it. it made
         no sense to replicate the storage for a very large graph
         just to instantiate a subgraph. when the subgraph was free'd,
         it did not release any storage owned by the parent graph.

   created -- 95sep22, cca
   ---------------------------------------------------------------------
*/
/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   typedef for the IVL and Ichunk objects
   --------------------------------------
*/
typedef struct _IVL      IVL    ;
typedef struct _Ichunk   Ichunk ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   type -- type of integer vector list
      IVL_NOTYPE  -- no type
      IVL_CHUNKED -- list are in possibly many chunks of data
      IVL_SOLO    -- each list is created using IVinit()
                     and free'd using IVfree()
      IVL_UNKNOWN -- storage for lists unknown, 
                     it is the user's responsibility
   maxnlist -- maximum number of lists
   nlist    -- number of lists
   tsize    -- total size of the lists
   sizes    -- vector of list sizes, size nlist
   p_vec    -- vector of list pointers, size nlist
   incr     -- increment for chunks, number of entries allocated
               when a new chunk is necessary, used only when
               type = IVL_CHUNKED
   chunk    -- first Ichunk structure, 
               NULL unless type = IVL_CHUNKED
   -------------------------------------------------------------
*/
struct _IVL {
   int      type     ;
   int      maxnlist ;
   int      nlist    ;
   int      tsize    ;
   int      *sizes   ;
   int      **p_vec  ;
   int      incr     ;
   Ichunk   *chunk   ;
} ;
/*
   -------------------------------------------------------
   the Ichunk object is the data structure 
   that handles chunks of storage for mode IVL_CHUNKED

   size  -- number of entries in the chunk,
            also the dimension of the array base[]
   inuse -- the number of entries in use,
            size - inuse is the number of free entries
   base  -- base address for storage used for list entries
   next  -- pointer to the next Ichunk object
   -------------------------------------------------------
*/
struct _Ichunk {
   int      size  ;
   int      inuse ;
   int      *base ;
   Ichunk   *next ;
} ;
/*--------------------------------------------------------------------*/
/*
   ----------------------------
   definitions for storage type
   ----------------------------
*/
#define IVL_NOTYPE     (-1)
#define IVL_CHUNKED      1
#define IVL_SOLO         2
#define IVL_UNKNOWN      3
#define IVL_INCR      1024

/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in basics.c   -----------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------
   simplest constructor
 
   created -- 95sep22, cca
   -----------------------
*/
IVL *
IVL_new (
   void
) ;
/*
   -----------------------
   set the default fields
 
   created -- 95sep22, cca
   -----------------------
*/
void
IVL_setDefaultFields (
   IVL   *ivl
) ;
/*
   --------------------------------------------------
   clear the data fields, releasing allocated storage
 
   created -- 95sep22, cca
   --------------------------------------------------
*/
void
IVL_clearData (
   IVL   *ivl
) ;
/*
   ------------------------------------------
   destructor, free's the object and its data
 
   created -- 95sep22, cca
   ------------------------------------------
*/
IVL *
IVL_free (
   IVL   *ivl
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in instance.c   ---------------------------------
------------------------------------------------------------------------
*/
/*
   -------------------------------------
   return the storage type of the object
 
   created -- 96dec06, cca
   -------------------------------------
*/
int
IVL_type (
   IVL   *ivl
) ;
/*
   ----------------------------------
   return the maximum number of lists
 
   created -- 96dec06, cca
   ----------------------------------
*/
int
IVL_maxnlist (
   IVL   *ivl
) ;
/*
   --------------------------
   return the number of lists
 
   created -- 96dec06, cca
   --------------------------
*/
int
IVL_nlist (
   IVL   *ivl
) ;
/*
   ----------------------------------
   return the total size of the lists
 
   created -- 96dec06, cca
   ----------------------------------
*/
int
IVL_tsize (
   IVL   *ivl
) ;
/*
   ----------------------------
   return the storage increment
 
   created -- 96dec06, cca
   ----------------------------
*/
int
IVL_incr (
   IVL   *ivl
) ;
/*
   -------------------------
   set the storage increment
 
   created -- 96dec06, cca
   -------------------------
*/
void
IVL_setincr (
   IVL   *ivl,
   int   incr
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in init.c   -------------------------------------
------------------------------------------------------------------------
*/
/*
   ------------------------------------------------------
   initialize given the type and maximum number of lists.
   used for type IVL_CHUNKED, IVL_SOLO or IVL_UNKNOWN
 
   created -- 95sep22, cca
   ------------------------------------------------------
*/
void
IVL_init1 (
   IVL   *ivl,
   int   type,
   int   maxnlist
) ;
/*
   ---------------------------------------------------------------
   initialize given the type, number of lists and their total size
   only used when type is IVL_CHUNKED.
 
   created -- 95sep22, cca
   ---------------------------------------------------------------
*/
void
IVL_init2 (
   IVL   *ivl,
   int   type,
   int   maxnlist,
   int   tsize
) ;
/*
   --------------------------------------
   initialize from a vector of list sizes.
   used with IVL_SOLO or IVL_CHUNKED.
 
   created -- 95sep22, cca
   --------------------------------------
*/
void
IVL_init3 (
   IVL   *ivl,
   int   type,
   int   maxnlist,
   int   sizes[]
) ;
/*
   ------------------------------------------------
   this method resizes the maximum number of lists,
   replacing the old sizes[] and p_vec[] vectors
   as necessary. the nlist value is NOT reset.
 
   created -- 96dec05, cca
   ------------------------------------------------
*/
void
IVL_setMaxnlist (
   IVL   *ivl,
   int   newmaxnlist
) ;
/*
   ------------------------------------------------
   this method resizes the number of lists,
   replacing the old sizes[] and p_vec[] vectors
   as necessary.
 
   created -- 96dec05, cca
   ------------------------------------------------
*/
void
IVL_setNlist (
   IVL   *ivl,
   int   newnlist
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in subIVL.c   -----------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------------
  purpose -- initialize subIVL from ivl
     if keeplistIV is not NULL then
        keeplistIV contains the lists of ivl to be placed in subIVL
    else
        all lists are placed in subIVL
     endif
     if keepentriesIV is not NULL then
        keepentriesIV contains the entries in the lists to be kept
     else
        all entries in the kept lists are placed in subIVL
     endif
 
   return value ---
      1 -- normal return
     -1 -- subIVL is NULL
     -2 -- ivl is NULL
     -3 -- keeplistIV is invalid
 
   created -- 98oct16, cca
   ----------------------------------------------------------------
*/
int
IVL_initFromSubIVL (
   IVL   *subIVL,
   IVL   *ivl,
   IV    *keeplistIV,
   IV    *keepentriesIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in listmanip.c   --------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------------------------
   fills size of and pointer to the entries of a list
 
   set *psize = size of list ilist
   set *pivec = base address of list ilist
 
   use as follows :
 
   IVL_listAndSize(ivl, ilist, &isize, &ivec) ;
   for ( i = 0 ; i < isize ; i++ ) ;
      do something with ivec[i] ;
   }
 
   created -- 95sep22, cca
   ----------------------------------------------------------------
*/
void
IVL_listAndSize (
   IVL   *ivl,
   int   ilist,
   int   *psize,
   int   **pivec
) ;
/*
   ----------------------------------------------------
   returns a pointer to the first element in list ilist
 
   to be used as a general iterator, e.g.,
 
   for ( pi = IVL_firstInList(ivl, ilist) ;
         pi != NULL ;
         pi = IVL_nextInList(ivl, ilist, pi) ) ;
      do something ;
   }
 
   created -- 95sep27, cca
   ----------------------------------------------------
*/
int *
IVL_firstInList (
   IVL   *ivl,
   int   ilist
) ;
/*
   ----------------------------------------------------
   returns a pointer to the next element in list ilist
 
   to be used as a general iterator, e.g.,
 
   for ( pi = IVL_firstInList(ivl, ilist) ;
         pi != NULL ;
         pi = IVL_nextInList(ivl, ilist, pi) ) ;
      do something ;
   }
 
   created -- 95sep27, cca
   ----------------------------------------------------
*/
int *
IVL_nextInList (
   IVL   *ivl,
   int   ilist,
   int   *pi
) ;
/*
   -----------------------------------------------------------------
   purpose -- to set or reset a list.
 
   ilist -- list id to set or reset
   isize -- size of the list
      if the present size of list ilist is smaller than isize,
         the old list is free'd (if ivl->type = IVL_SOLO)
         or lost (if ivl->type = IVL_CHUNKED)
         or un-set (if ivl->type = IVL_UNKNOWN)
         and new storage is allocated (for IVL_SOLO and IVL_CHUNKED)
   ivec  -- list vector
      if ivl->type is IVL_UNKNOWN then
         if ivec != NULL then
            we set the ivl list pointer to be ivec
         endif
      else if ivec != NULL
         we copy ivec[] into ivl's storage for the list
      endif
 
   created   -- 95sep27, cca
   last mods -- 95oct06, cca
      type = IVL_UNKNOWN, p_vec[ilist] set to ivec
      only when ivec is not NULL
      bug fixed, ivl->sizes[ilist] = isize ;
   -----------------------------------------------------------------
*/
void
IVL_setList (
   IVL   *ivl,
   int   ilist,
   int   isize,
   int   ivec[]
) ;
/*
   ----------------------------------------------------------------
   set a pointer to a list but don't allocate storage for the list.
   this method was needed when we form a subgraph with a boundary.
   lists for the interior vertices point into the parent graph,
   but lists for the boundary vertices must be allocated and owned.
   used only for type = IVL_CHUNKED. at some point in the future we
   should rethink the storage semantics for the IVL object.
 
   created -- 95nov11, cca
   ----------------------------------------------------------------
*/
void
IVL_setPointerToList (
   IVL   *ivl,
   int   ilist,
   int   isize,
   int   ivec[]
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in misc.c   -------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------------------
   purpose -- to make and return a 9-point adjacency structure
 
   input --
 
      n1    -- # of grid points in first direction
      n2    -- # of grid points in second direction
      ncomp -- # of components per grid point
   -----------------------------------------------------------
*/
IVL *
IVL_make9P (
   int   n1,
   int   n2,
   int   ncomp
) ;
/*
   ------------------------------------------------------------
   purpose -- to make and return a 27-point adjacency structure
 
   input --
 
      n1    -- # of grid points in first direction
      n2    -- # of grid points in second direction
      n3    -- # of grid points in second direction
      ncomp -- # of components per grid point
   ------------------------------------------------------------
*/
IVL *
IVL_make27P (
   int   n1,
   int   n2,
   int   n3,
   int   ncomp
) ;
/*
   ------------------------------------------------------------
   purpose -- to make and return a 13-point adjacency structure
 
   input --
 
      n1 -- # of grid points in first direction
      n2 -- # of grid points in second direction
 
   created -- 96feb01
   ------------------------------------------------------------
*/
IVL *
IVL_make13P (
   int   n1,
   int   n2
) ;
/*
   -----------------------------------------------------------
   purpose -- to make and return a 5-point adjacency structure
 
   input --
 
      n1 -- # of grid points in first direction
      n2 -- # of grid points in second direction
 
   created -- 96feb02
   -----------------------------------------------------------
*/
IVL *
IVL_make5P (
   int   n1,
   int   n2
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in util.c   -------------------------------------
------------------------------------------------------------------------
*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object
 
   created -- 95sep22, cca
   ----------------------------------------------
*/
int
IVL_sizeOf (
   IVL   *ivl
) ;
/*
   ---------------------------------------------------
   purpose -- to return the minimum entry in the lists
 
   created -- 95sep22, cca
   ---------------------------------------------------
*/
int
IVL_min (
   IVL   *ivl
) ;
/*
   ---------------------------------------------------
   purpose -- to return the maximum entry in the lists
 
   created -- 95sep22, cca
   ---------------------------------------------------
*/
int
IVL_max (
   IVL   *ivl
) ;
/*
   ----------------------------
   return the maximum list size
 
   created -- 95sep22, cca
   ----------------------------
*/
int
IVL_maxListSize (
   IVL   *ivl
) ;
/*
   -------------------------------
   return the sum of all the lists
 
   created -- 95sep29, cca
   -------------------------------
*/
int
IVL_sum (
   IVL   *ivl
) ;
/*
   -------------------------------------------
   sort the adjacency lists in ascending order
 
   created -- 95sep22, cca
   -------------------------------------------
*/
void
IVL_sortUp (
   IVL   *ivl
) ;
/*
   -----------------------------------------------------
   create an equivalence map
   if ( map[j] == map[k] ) then
      the lists for j and k are identical
   endif
   NOTE : each empty list is mapped to a different value
 
   return value -- pointer to map[] vector
 
   created -- 95mar15, cca
   -----------------------------------------------------
*/
int *
IVL_equivMap1 (
   IVL   *ivl
) ;
/*
   -----------------------------------------------------
   create an equivalence map
   if ( map[j] == map[k] ) then
      the lists for j and k are identical
   endif
   NOTE : each empty list is mapped to a different value
 
   return value -- pointer to map IV object
 
   created -- 96mar15, cca
   -----------------------------------------------------
*/
IV *
IVL_equivMap2 (
   IVL   *ivl
) ;
/*
   ------------------------------------
   overwrite each list with new entries
 
   created -- 96oct03, cca
   ------------------------------------
*/
void
IVL_overwrite (
   IVL   *ivl,
   IV    *oldToNewIV
) ;
/*
   ---------------------------------------------------
   given an IVL object and a map from old list entries
   to new list entries, create a new IVL object whose
   new list entries contain no duplicates.
 
   return value -- pointer to new IVL object
 
   created -- 96nov07, cca
   ---------------------------------------------------
*/
IVL *
IVL_mapEntries (
   IVL   *ivl,
   IV    *mapIV
) ;
/*
   ----------------------------------------------------------------
   IVL object ivl1 absorbs the lists and data from IVL object ivl2.
   the lists in ivl2 are mapped into lists in ivl1 using the mapIV
   IV object.
 
   created -- 96dec06, cca
   ----------------------------------------------------------------
*/
void
IVL_absorbIVL (
   IVL   *ivl1,
   IVL   *ivl2,
   IV    *mapIV
) ;
/*
   -----------------------------------------------------------------
   purpose -- expand the lists in an IVL object.
 
   this method was created in support of a symbolic factorization.
   an IVL object is constructed using a compressed graph.
   it must be expanded to reflect the compressed graph.
   the number of lists does not change (there is one list per front)
   but the size of each list may change. so we create a new IVL
   object that contains entries for the uncompressed graph.
 
   created -- 97feb13, cca
   -----------------------------------------------------------------
*/
IVL *
IVL_expand (
   IVL   *ivl,
   IV    *eqmapIV
) ;
/*--------------------------------------------------------------------*/
/*
------------------------------------------------------------------------
---   procedures found in IO.c   ---------------------------------------
------------------------------------------------------------------------
*/
/*
   -----------------------------------------------
   purpose -- to read in an IVL object from a file
 
   input --
 
      fn -- filename, must be *.ivlb or *.ivlf
 
   return value -- 1 if success, 0 if failure
 
   created -- 95sep29, cca
   -----------------------------------------------
*/
int
IVL_readFromFile ( 
   IVL    *ivl, 
   char   *fn 
) ;
/*
   ------------------------------------------------------
   purpose -- to read an IVL object from a formatted file
 
   return value -- 1 if success, 0 if failure
 
   created -- 95sep29, cca
   ------------------------------------------------------
*/
int
IVL_readFromFormattedFile (
   IVL    *ivl,
   FILE   *fp
) ;
/*
   ---------------------------------------------------
   purpose -- to read an IVL object from a binary file
 
   return value -- 1 if success, 0  if failure
 
   created -- 95sep29, cca
   ---------------------------------------------------
*/
int
IVL_readFromBinaryFile (
   IVL    *ivl,
   FILE   *fp
) ;
/*
   -------------------------------------------
   purpose -- to write an IVL object to a file
 
   input --
 
      fn -- filename
        *.ivlb -- binary
        *.ivlf -- formatted
        anything else -- for human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95sep29, cca
   -------------------------------------------
*/
int
IVL_writeToFile (
   IVL    *ivl,
   char   *fn
) ;
/*
   -----------------------------------------------------
   purpose -- to write an IVL object to a formatted file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95sep29, cca
   -----------------------------------------------------
*/
int
IVL_writeToFormattedFile (
   IVL    *ivl,
   FILE   *fp
) ;
/*
   --------------------------------------------------
   purpose -- to write an IVL object to a binary file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95sep29, cca
   --------------------------------------------------
*/
int
IVL_writeToBinaryFile (
   IVL    *ivl,
   FILE   *fp
) ;
/*
   -------------------------------------------------
   purpose -- to write an IVL object for a human eye
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95sep29, cca
   -------------------------------------------------
*/
int
IVL_writeForHumanEye (
   IVL    *ivl,
   FILE   *fp
) ;
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the IVL object
 
   return value -- 1 if success, 0 otherwise
 
   created -- 95sep29, cca
   ---------------------------------------------------------
*/
int
IVL_writeStats (
   IVL    *ivl,
   FILE   *fp
) ;
/*--------------------------------------------------------------------*/
