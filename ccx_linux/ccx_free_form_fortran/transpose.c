/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2019 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h> 
#include "CalculiX.h"
#include "mortar.h"
/**
 * \brief function to transpose sparse matrix
 * Author: Saskia Sitzmann
 * @param[in] au 	matrix A values
 * @param[in] jq	matrix A column pointer to irow
 * @param[in] irow	matrix A rows
 * @param[in] dim	dimension of matrix A (dim x dim)
 * @param[out] au_t 	matrix A^T values
 * @param[out] jq_t	matrix A^T column pointer to irow
 * @param[out] irow_t	matrix A^T rows
**/
void transpose(double *au,ITG  *jq, ITG  *irow,ITG  *dim,
	       double *au_t, ITG  *jq_t,ITG  *irow_t){
  
  ITG  i,j,k,nzs_t=jq[*dim]-1,*mast1=NULL;
  /*transpose*/	
  
  NNEW(mast1,ITG,nzs_t);
  
  for(j=0;j<nzs_t;j++){au_t[j]=au[j];}  
  /* mast1 contains the rows of the original matrix,
     irowbt the columns */  
  for(j=0;j<*dim;j++){
    for(i=jq[j]-1;i<jq[j+1]-1;i++){
      mast1[i]=irow[i];
      irow_t[i]=j+1;
    }
  }  
  /*for(j=0;j<nzs_t;j++){
    printf("\t transpose: column %d row %d v %e \n",mast1[j],irow_t[j],au_t[j]);
    }*/
  
  matrixsort(au_t,mast1,irow_t,jq_t,&nzs_t,dim); 
  free(mast1);
  return;
}
