/*  utilities.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   return the maximum tag value

   created -- 98jan08, cca
   ----------------------------
*/
int
maxTagMPI (
   MPI_Comm   comm
) {
int   *iptr, flag, rc, tag_bound ;

iptr = &tag_bound ;
rc = MPI_Attr_get(comm, MPI_TAG_UB, (void **) &iptr, &flag) ;
if ( flag == 0 ) {
   tag_bound = -1 ;
} else {
   tag_bound = *iptr ;
}
return(tag_bound) ; }
   
/*--------------------------------------------------------------------*/
