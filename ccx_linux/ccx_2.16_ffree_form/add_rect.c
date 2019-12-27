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
 *Calculate C=A+B for sparse matrices A,B
 * 
 * @param [in] au_1		sparse matrix A
 * @param [in] irow_1           row numbers for matrix A
 * @param [in] jq_1             colum numbers for matrix A
 * @param [in] n_1              # rows  for matrix A 
 * @param [in] m_1              # columns  for matrix A
 * @param [in] au_2		sparse matrix B
 * @param [in] irow_2           row numbers for matrix B
 * @param [in] jq_2             colum numbers for matrix B
 * @param [in] n_2              # rows  for matrix B
 * @param [in] m_2              # columns  for matrix B
 * @param [out] au_rp		sparse matrix C
 * @param [out] irow_rp          row numbers for matrix C
 * @param [out] jq_r            colum numbers for matrix C
 * @param [out] nzs             # nonzero entries in C 
 *
*/
void add_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1, ITG m_1,
	      double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2, ITG m_2,
	      double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs){
  
  /*Result fields*/
  ITG *irow=NULL,ifree=1,numb,icol,i,j,k,l,m,carre=0,kflag=2,istart,icounter,
    pt1,pt2,row1,row2;
  double *autmp=NULL,value;
  clock_t debut;
  clock_t fin;
  
  debut=clock();
  if((m_1!=m_2)||(n_1!=n_2)){
    printf("Error in mutli_rec : Matrix sizes are not compatible\n");
    return;
  } 
  
  //        nzs=jq_1[m_1]+jq_2[m_2]-2;
  //printf("nzs add_rect = %d\n",nzs);
  //        irow=NNEW(int,nzs);
  //        au=NNEW(double,nzs);
  irow=*irow_rp;
  autmp=*au_rp; 
  
  jq_r[0]=1;
  
  for(j=0;j<m_2;j++){
    pt1=jq_1[j]-1;
    pt2=jq_2[j]-1;
    m=j+1; 
    while((pt1<jq_1[m]-1)||(pt2<jq_2[m]-1)){
      if ((pt1<jq_1[m]-1)&&(pt2<jq_2[m]-1)){ //normal case
	row1=irow_1[pt1];
	row2=irow_2[pt2];
	if (row1==row2){
	  value=au_1[pt1]+au_2[pt2];
	  //printf("\taddrect: r %d  c %d v %e \n",row1,m,value);
	  insertas_ws(&irow,&row1,&m,&ifree,nzs,&value,&autmp);
	  pt1++;
	  pt2++;
	}else{
	  if (row1<row2) {
	    //printf("\taddrect: r %d  c %d v %e \n",row1,m,au_1[pt1]);
	    insertas_ws(&irow,&row1,&m,&ifree,nzs,&au_1[pt1],&autmp);
	    pt1++;
	  }else{
				  //printf("\taddrect: r %d  c %d v %e \n",row2,m,au_2[pt2]);
	    insertas_ws(&irow,&row2,&m,&ifree,nzs,&au_2[pt2],&autmp);
	    pt2++;
	  }
	}
      }
      else{ 
	if (pt1<jq_1[m]-1){ //column 2 finished
	  //printf("\taddrect: r %d  c %d v %e \n",irow_1[pt1],m,au_1[pt1]);
	  insertas_ws(&irow,&irow_1[pt1],&m,&ifree,nzs,&au_1[pt1],&autmp);
	  pt1++;				
	}
	else{ //column 1 finished
	  //printf("\taddrect: r %d  c %d v %e \n",irow_2[pt2],m,au_2[pt2]);
	  insertas_ws(&irow,&irow_2[pt2],&m,&ifree,nzs,&au_2[pt2],&autmp);
	  pt2++;						
	}
      }
    }	   
    jq_r[m]=ifree;
    //printf("\t\t ifree %d\n",ifree);
  }
  
  *nzs=ifree-1;
  RENEW(irow,ITG,*nzs);
  RENEW(autmp,double,*nzs);
  *irow_rp=irow;*au_rp=autmp;
  fin= clock();
  //printf("add_rect : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
  return;
}
