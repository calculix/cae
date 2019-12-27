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
 * \brief sparse matrix multiplication a_r=a_1^T*a_2
 * 
 * @param [in] au_1 		values of matrix 1
 * @param [in] irow_1		rows of matrix 1
 * @param [in] jq_1		column pointer to iorw_1
 * @param [in]  n_1		number of rows matrix 1  
 * @param [in]  m_1		number of columns matrix 1
  * @param [in] au_2 		values of matrix 2
 * @param [in] irow_2		rows of matrix 2
 * @param [in] jq_2		column pointer to iorw_2
 * @param [in]  n_2		number of rows matrix 2  
 * @param [in]  m_2		number of columns matrix 2
 * @param [out] au_rp		values of result matrix
 * @param [out] irow_rp		rows of result matrix
 * @param [out] jq_r		column pointer to iorw_r
 * @param [out] nzs		number of non-zero entries in au_r
**/
void multi_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1, ITG m_1,
		double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2, ITG m_2,
		double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs){ 
  
  /*Result fields*/
  ITG *irow=NULL,ifree=1,numb,icol,i,j,k,l,m,carre=0,kflag=2,istart,icounter;
  ITG flag=0;
  double *autmp=NULL,value;
  clock_t debut;
  clock_t fin;
  /*Perform a_1T*a_2*/
  //printf("multi_rect: n1 %" ITGFORMAT " m1 %" ITGFORMAT " n2 %" ITGFORMAT " m2 %" ITGFORMAT " \n", n_1,m_1,n_2,m_2);
  debut=clock();
  if (n_1!=n_2) {
    printf("Error in mutli_rec : Matrix sizes are not compatible\n");
    return;
  } 
  
  //        nzs=n_1*m_2;
  irow=*irow_rp;
  autmp=*au_rp; 
  
  if (n_1==m_2) carre=1;
  
  jq_r[0]=1;
  for(j=0;j<m_2;j++){
    m=j+1;
    if(jq_2[j+1]-jq_2[j]>0){
      for (i=0;i<m_1;i++){
	l=i+1;
	flag=0;
	value=0.0;
	multi_scal(au_1,irow_1,jq_1,au_2,irow_2,jq_2,i,j,&value,&flag);
	if (flag!=0) insertas_ws(&irow,&l,&m,&ifree,nzs,&value,&autmp);
      }
    }
    jq_r[m]=ifree;
    //printf("multi_rect: column %" ITGFORMAT " jq %" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT " \n",m,jq_r[m-1],jq_r[m],m_2);
  }
  
  /* Sort the column and compute jq*/ 
  //printf("multi_rect: column %" ITGFORMAT " jq %" ITGFORMAT " %" ITGFORMAT " \n",m_2,ifree,jq_r[m_2]);
  *nzs=ifree-1;
  RENEW(autmp,double,*nzs);
  RENEW(irow,ITG,*nzs);
  *irow_rp=irow;*au_rp=autmp;   
  fin= clock();
  //printf("multi_rect : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
  return;
}
