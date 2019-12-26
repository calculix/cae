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
#include <string.h>
#include <time.h>
#include "CalculiX.h"
#include "mortar.h"
/** 
 * \brief sort unsorted sparse matrix au 
 * Author: Saskia Sitzmann
 *
 * @param [in,out] au		matrix values
 * @param [in,out] mast1	column numbers
 * @param [in,out] irow		row numbers
 * @param [out] jq 		column pointer to irow 
 * @param [in,out] nzs		number of non-zero values in au
 * @param [in,out] ndim		dimention of matrix au ndim x ndim
**/
void matrixsort(double *au,ITG *mast1,ITG *irow,ITG *jq, 
		ITG *nzs,ITG *ndim){
  
  ITG  i,j,jj,k,kk,l,m,ll,kflag,numb;
  
  /* Sort mast1, irow and au; 
     Outcome: the values in field au are sorted, column by
     column; no sorting is done within the columns */
  
  kflag=2;
  FORTRAN(isortiid,(mast1,irow,au,nzs,&kflag));
  
  /*  fill in jqbd
      jqbd(i): first element in field aubd belonging to column i  */
  
  j=0;
  for(i=0;i<*ndim;i++){
    if(j==*nzs){
      for(k=i;k<*ndim;k++) 
	jq[k]=*nzs+1;
      break;
    }
    
    if(mast1[j]!=i+1){
      jq[i]=j+1;
      continue;
    }
    
    jq[i]=j+1;
    
    while(1){
      j++;
      if(j==*nzs) break;
      if(mast1[j]!=i+1) break;
    }
  }
  
  jq[*ndim]=*nzs+1;
  
  /* Sorting of the rows*/
  
  for (i=0;i<*ndim;i++){
    if(jq[i+1]-jq[i]>0){
      numb=jq[i+1]-jq[i]; 
      FORTRAN(isortid,(&irow[jq[i]-1],&au[jq[i]-1],&numb,&kflag));
    }
  }
  
  return;
}
